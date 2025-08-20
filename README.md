
# Risk-Based Treatment Effect Analysis (Conversion Uplift)

This repository contains an end-to-end R workflow for estimating heterogeneous treatment effects using baseline risk stratification. Specifically, we build a baseline risk model using untreated customers, segment the population into quartiles of baseline conversion risk, and estimate how treatment effectiveness differs among these subgroups.

---

## 🎯 Objective

To evaluate whether a marketing treatment (e.g., discount, BOGO) produces different conversion lift depending on each customer's baseline risk of converting without treatment.

---

## ✅ Key Findings (Results)

| Risk Group | n_treated | n_control | rate_treated | rate_control | Difference (Treated - Control) |
|------------|-----------|-----------|---------------|----------------|-------------------------------|
| Q1 (lowest risk) | 9,628 | 4,722 | 0.110 | 0.074 | **+0.036** |
| Q2 | 9,597 | 4,811 | 0.122 | 0.061 | +0.061 |
| Q3 | 9,621 | 4,761 | 0.120 | 0.061 | +0.059 |
| Q4 (highest risk) | 9,585 | 4,764 | 0.248 | 0.149 | **+0.099** |

- The **treatment effect is positive in all quartiles**.
- The effect is **smallest in Q1 (+3.6%)** and **largest in Q4 (+9.9%)**.
- The plot shows an **increasing trend** in treatment effect from Q1 → Q4, even though there is uncertainty (wide error bars).

### Interpretation

Customers who already had a higher baseline probability of conversion benefited much more from receiving the treatment. This suggests targeting higher-risk or medium-high-risk customers leads to the greatest incremental gain.

---

## 📊 Plot Explanation

![Treatment Effect Plot] <img width="695" height="439" alt="Rplot" src="https://github.com/user-attachments/assets/fc11bc9d-a7a4-4530-9f80-45843bf1e24d" />


- **X-axis:** Risk Group (Quartiles of baseline risk)
- **Y-axis:** Difference in conversion rate between treated vs control
- **Points:** Estimated average uplift
- **Error bars:** 95% confidence intervals
- The line slopes gradually upward — meaning treatment becomes more effective as baseline risk increases.

---

## 🔁 Full Workflow Explained

1. **Load Data**
   - Preprocessed CSV is loaded into R.
   - Outcome (`conversion`) and `treatment` are coerced to numeric 0/1.

2. **Select Predictors**
   - Exclude conversion, treatment, and offer; everything else is a potential predictor.

3. **Check Assumptions (Logistic Regression)**
   - Multicollinearity: all VIF < 5 ✓
   - Linearity of logit (Box–Tidwell): recency & history are linear ✓
   - Skewness: `history` is skewed → apply Yeo-Johnson transform
   - Separation: none detected ✓

4. **Build Baseline Risk Model**
   ```r
   ControlsOnly <- df[df$treatment == 0, ]
   fit_ctrl <- glm(conversion ~ ., data = ControlsOnly, family = binomial())
## 🔍 What is Baseline Risk?

**Baseline risk** is the probability that a customer would convert **even without** receiving any treatment or promotion. It measures each customer's natural likelihood to buy on their own.

- High baseline risk = likely to convert anyway
- Low baseline risk = unlikely to convert without help

This is important because two customers can respond differently to the same treatment depending on their baseline risk.

---

## 🛠 How the Baseline Risk Model Was Built

We estimated baseline risk using a logistic regression model trained only on control (untreated) customers:

1. Filtered dataset to include rows where `treatment = 0`.
2. Ran a logistic regression:
   ```r
   fit_ctrl <- glm(conversion ~ ., data = controls_only, family = binomial())
# Risk-Based Subgroup Analysis

This repository implements **Risk-Based Subgroup Analysis** for analyzing heterogeneous treatment effects.

## Project Steps
1. Load and preprocess data
2. Fit baseline risk model on control group
3. Score baseline risk for all individuals
4. Create risk-based subgroups (quartiles by default)
5. Estimate treatment effects within each subgroup
6. Visualize treatment effect by risk group
7. Export results

## Requirements
- R >= 4.2
- Packages: dplyr, tidyr, purrr, stringr, ggplot2, forcats, broom,
  car, bestNormalize, brglm2, logistf, glmnet, pROC, DescTools,
  rsample, yardstick, rms, moments

