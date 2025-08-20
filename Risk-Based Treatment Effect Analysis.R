# ================================================================
# Risk-Based Subgroup Analysis: Assumptions, Fixes, and Implementation
# ================================================================

library(readr)
# ==== LOAD DATA ==== #
df <- read.csv("C:/Users/SANJANA/Downloads/preprocessed_dataset (2).csv")

# ---- Packages ----
pkgs <- c(
  "dplyr","tidyr","purrr","stringr","ggplot2","forcats","broom",
  "car","bestNormalize","brglm2","logistf","glmnet","pROC","DescTools",
  "rsample","yardstick","rms"
)
new_pkgs <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if(length(new_pkgs)) install.packages(new_pkgs, dependencies = TRUE)
invisible(lapply(pkgs, library, character.only = TRUE))

if(!require(moments)) install.packages("moments")
library(moments)

# ================================================================
# ==== Harmonize/Assert key columns ==============================
# ================================================================
# Outcome column
outcome_col   <- "conversion"
if(!(outcome_col %in% names(df))) stop("`conversion` column not found.")

# Treatment column: prefer existing binary 'treatment'; if absent, infer from 'offer'
if("treatment" %in% names(df)) {
  df <- df %>% mutate(treatment = as.integer(treatment %in% c(1, "1", TRUE)))
} else if("offer" %in% names(df)) {
  # Any promo (Discount/BOGO) = 1, No Offer = 0
  df <- df %>%
    mutate(
      treatment = ifelse(str_to_lower(as.character(offer)) %in%
                           c("discount","buy one get one","bogo","buy-one-get-one","buy_one_get_one"),
                         1L, 0L)
    )
} else {
  stop("Neither `treatment` nor `offer` found to define treatment.")
}

# Ensure outcome is binary 0/1
df <- df %>% mutate(
  conversion = as.integer(conversion %in% c(1, "1", TRUE))
)

# Quick sanity check (per your preprocessing report — no missing expected) :contentReference[oaicite:3]{index=3}
miss_summary <- sapply(df, function(x) sum(is.na(x)))
if(any(miss_summary > 0)) {
  message("Warning: Missing values detected. Imputing simple strategies for predictors (median/mode).")
}

# ================================================================
# ==== Identify predictors for the baseline risk model ===========
# ================================================================
# We follow your covariate structure; exclude outcome & treatment.  (PDF method relies on baseline risk from controls) :contentReference[oaicite:4]{index=4}
candidate_x <- setdiff(names(df), c(outcome_col, "treatment","offer"))

# Keep only columns that are numeric, integer, or factor/logical/binary dummies
# Factors will be one-hot via model.matrix later if needed.
# We'll detect binary numeric dummies too.
is_binary01 <- function(v) {
  u <- unique(na.omit(v))
  is.numeric(v) && length(setdiff(u, c(0,1))) == 0
}

# Convert character categoricals to factor
df <- df %>%
  mutate(across(all_of(candidate_x), ~ if(is.character(.x)) factor(.x) else .x))

# Continuous numeric vars (non-binary numeric)
num_vars <- candidate_x[sapply(df[candidate_x], is.numeric)]
cont_vars <- num_vars[sapply(df[num_vars], function(v) !is_binary01(v))]

# Binary numeric predictors
bin_vars <- num_vars[sapply(df[num_vars], is_binary01)]

# Factor predictors
fac_vars <- candidate_x[sapply(df[candidate_x], is.factor)]

# ================================================================
# ==== Train/Test split is optional; we fit on all controls ======
# ================================================================
controls <- df %>% filter(treatment == 0)

# Impute simple strategies for controls if needed (should be none per preprocessing) :contentReference[oaicite:5]{index=5}
imp_num <- function(v) { v[is.na(v)] <- median(v, na.rm = TRUE); v }
imp_fac <- function(v) { v[is.na(v)] <- fct_explicit_na(v, na_level = "Missing"); v }

controls <- controls %>%
  mutate(across(all_of(num_vars), imp_num),
         across(all_of(fac_vars), imp_fac))

# ================================================================
# ==== Assumption checks for the baseline risk model (controls) ==
# ================================================================

# Build a model matrix (no outcome, no treatment) for VIF, etc.
# We'll start with a simple main-effects logistic model.
# We'll one-hot-encode factors safely.
mm_formula <- as.formula(paste(outcome_col, "~", paste(candidate_x, collapse = " + ")))
mm <- model.matrix(~ .,
                   data = controls %>% select(all_of(candidate_x)))
mm <- as.data.frame(mm[, -1, drop = FALSE])  # drop intercept column

# Align outcome
y_ctrl <- controls[[outcome_col]]

# ---- Multicollinearity (VIF) ----
# Fit a preliminary GLM to compute VIF on numeric frame
# car::vif expects an lm/glm with non-singular design
glm_pre <- glm(y_ctrl ~ ., data = mm, family = binomial())
vifs <- tryCatch(car::vif(glm_pre), error = function(e) NA)
print(vifs)

# Flag high VIF (>5 as a soft rule)
high_vif <- names(vifs)[which(vifs > 5)]
if(length(high_vif)) {
  message("High VIF detected for: ", paste(high_vif, collapse = ", "),
          ". Will consider dropping or ridge-penalizing.")
}

# ---- Linearity of the logit for continuous predictors ----
# Box-Tidwell test: requires positive values; use shifted if non-positive
bt_results <- list()
if(length(cont_vars)) {
  controls_bt <- controls
  # shift any non-positive variables slightly
  shifts <- sapply(cont_vars, function(v) {
    mn <- min(controls_bt[[v]], na.rm = TRUE)
    if(mn <= 0) abs(mn) + 1e-3 else 0
  })
  for(v in cont_vars) {
    x <- controls_bt[[v]] + shifts[[v]]
    # Box-Tidwell uses interaction with log(x)
    df_bt <- data.frame(y = controls_bt[[outcome_col]], x = x)
    # exclude zeros post-shift (shouldn't have any)
    fit_bt <- glm(y ~ x + I(x*log(x)), data = df_bt, family = binomial())
    p_nonlin <- summary(fit_bt)$coefficients["I(x * log(x))","Pr(>|z|)"]
    bt_results[[v]] <- p_nonlin
  }
  print(bt_results)
}

nonlinear_vars <- names(bt_results)[which(unlist(bt_results) < 0.05)]
if(length(nonlinear_vars)) {
  message("Nonlinearity in logit detected for: ", paste(nonlinear_vars, collapse = ", "),
          ". Will apply Yeo-Johnson normalization + standardization.")
}

# ---- Skewness check for continuous predictors ----

skews <- sapply(cont_vars, function(v) moments::skewness(controls[[v]], na.rm = TRUE))
print(skews)
skewed_vars <- names(skews)[which(abs(skews) > 1)]

# ================================================================
# ==== Transformations (if needed) ===============================
# ================================================================
to_transform <- union(nonlinear_vars, skewed_vars)

transforms <- list()
controls_tf <- controls

if(length(to_transform)) {
  for(v in to_transform) {
    # Yeo-Johnson handles zeros/negatives
    bn <- bestNormalize(controls_tf[[v]], quiet = TRUE)
    controls_tf[[v]] <- as.numeric(predict(bn))
    transforms[[v]] <- bn
  }
}
# Center & scale continuous predictors (stabilize / aid convergence)
scale_vars <- union(cont_vars, to_transform)
scale_params <- list()
if(length(scale_vars)) {
  for(v in intersect(scale_vars, names(controls_tf))) {
    m <- mean(controls_tf[[v]], na.rm = TRUE); s <- sd(controls_tf[[v]], na.rm = TRUE)
    if(is.finite(s) && s > 0) {
      controls_tf[[v]] <- (controls_tf[[v]] - m)/s
      scale_params[[v]] <- c(mean = m, sd = s)
    }
  }
}

# ================================================================
# ==== Separation check & choose modeling strategy ===============
# ================================================================
# Build final model.frame for controls after transforms
x_ctrl <- controls_tf %>% select(all_of(candidate_x))
y_ctrl <- controls_tf[[outcome_col]]

# detect separation
sep <- tryCatch(brglm2::detect_separation(
  formula = as.formula(paste(outcome_col, "~ .")),
  data = data.frame(controls_tf %>% select(all_of(c(outcome_col, candidate_x)))),
  family = binomial()
), error = function(e) NULL)

use_firth <- FALSE
if(!is.null(sep) && (isTRUE(sep$separation) || isTRUE(sep$quasi_separation))) {
  message("Separation detected. Using Firth-penalized logistic regression.")
  use_firth <- TRUE
} else if(length(high_vif)) {
  message("High VIF present; considering ridge-penalized logistic regression.")
}

# ================================================================
# ==== Fit baseline risk model on CONTROLS only ==================
# ================================================================
fit_ctrl <- NULL
mm_ctrl <- model.matrix(~ ., data = x_ctrl)[, -1, drop = FALSE]

if(use_firth) {
  d_ctrl <- data.frame(y = y_ctrl, x_ctrl)
  fit_ctrl <- logistf::logistf(
    formula = as.formula(paste("y ~", paste(names(x_ctrl), collapse = " + "))),
    data = d_ctrl
  )
} else if(length(high_vif)) {
  # Ridge-penalized logistic (alpha=0)
  cvfit <- cv.glmnet(x = as.matrix(mm_ctrl), y = y_ctrl, family = "binomial", alpha = 0)
  fit_ctrl <- list(type = "glmnet", cv = cvfit, xnames = colnames(mm_ctrl))
} else {
  d_ctrl <- data.frame(y = y_ctrl, x_ctrl)
  fit_ctrl <- glm(
    formula = as.formula(paste("y ~", paste(names(x_ctrl), collapse = " + "))),
    data = d_ctrl, family = binomial()
  )
}

# ---- Risk model performance (controls) ----
# AUC, Brier, optional calibration
pred_ctrl <- if(is.list(fit_ctrl) && !is.null(fit_ctrl$type) && fit_ctrl$type == "glmnet") {
  p <- predict(fit_ctrl$cv, newx = as.matrix(mm_ctrl), type = "response", s = "lambda.min")
  as.numeric(p)
} else {
  as.numeric(predict(fit_ctrl, type = "response"))
}
auc_ctrl <- pROC::auc(y_ctrl, pred_ctrl)
brier_ctrl <- mean((pred_ctrl - y_ctrl)^2)

cat(sprintf("Baseline risk model (controls) AUC = %.3f; Brier = %.4f\n", auc_ctrl, brier_ctrl))

# Optional: calibration (suppressed by default)
# rms::val.prob(pred_ctrl, y_ctrl)  # uncomment to see calibration plot in R

# ================================================================
# ==== Apply the same transforms to ALL data & score risk ========
# ================================================================
apply_transforms <- function(df_in) {
  df_out <- df_in
  # Yeo-Johnson transforms
  if(length(transforms)) {
    for(v in names(transforms)) {
      if(v %in% names(df_out)) df_out[[v]] <- as.numeric(predict(transforms[[v]], newdata = df_out[[v]]))
    }
  }
  # Center/scale
  if(length(scale_params)) {
    for(v in names(scale_params)) {
      if(v %in% names(df_out)) {
        m <- scale_params[[v]]["mean"]; s <- scale_params[[v]]["sd"]
        df_out[[v]] <- (df_out[[v]] - m)/s
      }
    }
  }
  df_out
}

df_tf <- apply_transforms(df)
X_all <- df_tf %>% select(all_of(candidate_x))
MM_all <- model.matrix(~ ., data = X_all)[, -1, drop = FALSE]

risk_hat <- if(is.list(fit_ctrl) && !is.null(fit_ctrl$type) && fit_ctrl$type == "glmnet") {
  as.numeric(predict(fit_ctrl$cv, newx = as.matrix(MM_all), type = "response", s = "lambda.min"))
} else {
  d_all <- data.frame(X_all)
  as.numeric(predict(fit_ctrl, newdata = d_all, type = "response"))
}

df_scored <- df_tf %>%
  mutate(baseline_risk = risk_hat)

# ================================================================
# ==== Create risk-based subgroups (quartiles by default) ========
# ================================================================
K <- 4  # change to 5 or deciles if desired
df_scored <- df_scored %>%
  mutate(risk_group = ntile(baseline_risk, K)) %>%
  mutate(risk_group = factor(risk_group, levels = 1:K,
                             labels = paste0("Q", 1:K)))

# ================================================================
# ==== Estimate treatment effects within each risk group ==========
# ================================================================
estimate_effects <- function(dat) {
  # Difference in conversion rates (treated - control), RR, OR
  d_t <- dat %>% filter(treatment == 1)
  d_c <- dat %>% filter(treatment == 0)
  n_t <- nrow(d_t); n_c <- nrow(d_c)
  y_t <- sum(d_t[[outcome_col]]); y_c <- sum(d_c[[outcome_col]])
  
  # Risk (proportion)
  p_t <- if(n_t > 0) y_t / n_t else NA_real_
  p_c <- if(n_c > 0) y_c / n_c else NA_real_
  diff <- p_t - p_c
  
  # CI for difference in proportions (approx via prop.test for each & delta)
  ci_t <- if(n_t > 0) prop.test(y_t, n_t)$conf.int else c(NA, NA)
  ci_c <- if(n_c > 0) prop.test(y_c, n_c)$conf.int else c(NA, NA)
  # Simple normal approx for diff
  se_diff <- sqrt( (p_t*(1-p_t))/max(1,n_t) + (p_c*(1-p_c))/max(1,n_c) )
  ci_diff <- if(is.finite(se_diff)) diff + c(-1,1)*1.96*se_diff else c(NA, NA)
  
  # Risk ratio & OR via small eps
  eps <- 1e-6
  rr <- (p_t + eps)/(p_c + eps)
  or <- ((y_t + 0.5)/(n_t - y_t + 0.5))/((y_c + 0.5)/(n_c - y_c + 0.5))
  
  # Logistic regression within group for treatment effect (OR w/ p-value)
  fit_g <- tryCatch(
    glm(as.formula(paste(outcome_col,"~ treatment")), data = dat, family = binomial()),
    error = function(e) NULL
  )
  or_glm <- pval <- NA_real_
  if(!is.null(fit_g)) {
    co <- summary(fit_g)$coefficients
    if("treatment" %in% rownames(co)) {
      or_glm <- exp(co["treatment","Estimate"])
      pval   <- co["treatment","Pr(>|z|)"]
    }
  }
  
  tibble(
    n_treated = n_t, conv_treated = y_t, rate_treated = p_t,
    n_control = n_c, conv_control = y_c, rate_control = p_c,
    diff_in_props = diff, ci_diff_low = ci_diff[1], ci_diff_high = ci_diff[2],
    risk_ratio = rr, odds_ratio_cc = or, odds_ratio_glm = or_glm,
    p_value_logit = pval
  )
}

effects_by_group <- df_scored %>%
  group_by(risk_group) %>%
  group_modify(~ estimate_effects(.x)) %>%
  ungroup()

print(effects_by_group)

# ================================================================
# ==== Optional: Visualization of treatment effect by risk group ==
# ================================================================
# Difference in proportions with 95% CI
ggplot(effects_by_group, aes(x = risk_group, y = diff_in_props, group = 1)) +
  geom_line() +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = ci_diff_low, ymax = ci_diff_high), width = 0.1) +
  labs(title = "Treatment Effect (Treated - Control) by Baseline Risk Quartile",
       x = "Risk Group (Quartiles of Baseline Risk)",
       y = "Difference in Conversion Rate") +
  theme_minimal()

# ================================================================
# ==== Outputs (tidy) ============================================
# ================================================================
assumption_report <- list(
  missing_summary = miss_summary,
  vif = vifs,
  box_tidwell_p = bt_results,
  skewness = skews,
  nonlinear_vars = nonlinear_vars,
  high_vif = high_vif,
  separation_detected = use_firth,
  auc_controls = as.numeric(auc_ctrl),
  brier_controls = brier_ctrl,
  transforms_applied = names(transforms),
  scaling_applied = names(scale_params)
)

print(assumption_report)

# A compact table you can export:
final_te_table <- effects_by_group %>%
  select(risk_group, n_treated, n_control,
         rate_treated, rate_control, diff_in_props,
         ci_diff_low, ci_diff_high, risk_ratio, odds_ratio_glm, p_value_logit)

print(final_te_table)


write.csv(final_te_table, "risk_based_subgroup_effects.csv", row.names = FALSE)

