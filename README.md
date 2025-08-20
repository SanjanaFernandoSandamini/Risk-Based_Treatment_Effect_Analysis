
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
