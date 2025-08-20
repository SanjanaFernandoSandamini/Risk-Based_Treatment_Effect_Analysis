# Dataset Description and Preprocessing (README)

## 📊 Data Source Description

This dataset is from a **fictional marketing promotion campaign** created for educational and research purposes. Although not sourced from a real company, it simulates realistic marketing data and customer behavior patterns.

It includes:
- Multichannel customer interactions (Phone, Web, and Multichannel)
- Demographic + behavioral data
- Geographic coverage across **Suburban**, **Urban**, and **Rural** regions

---

## 🧾 Customer Behavior Variables
- **Recency**: Months since last purchase (customer engagement)
- **History**: Total past purchase value (lifetime value)
- **Used_discount**: Whether customer used discount offers before
- **Used_bogo**: Whether customer used buy-one-get-one offers

---

## 🌍 Demographic & Geographic Information
- **Zip_code**: Suburban / Urban / Rural
- **Is_referral**: Whether customer came via referral program

---

## 📣 Campaign-Specific Variables
- **Channel**: Preferred interaction method (Phone / Web / Multichannel)
- **Offer**: Type of promotion (Discount / BOGO / No Offer)
- **Conversion**: Target variable — 1 if customer made a purchase, 0 otherwise

---

## 🧹 Data Preprocessing & Cleaning

The original data had **64,000 records** with 9 columns.

### ✅ Duplicate Records
- Found **6,603 duplicate rows** (~10%)
- All duplicates removed
- Final dataset: **57,397 unique observations**

### ✅ Missing Values
- **No missing values** detected in any column
- No imputation required

---

## 🧪 Treatment Variable Creation

To study treatment effects, `offer` was converted into a binary variable:

| Group             | Count   | Percentage |
|------------------|---------|------------|
| Control (No Offer) | 19,072  | 33.2%      |
| Treatment (Any Offer) | 38,325 | 66.8%      |

---

## 🧾 Categorical Encoding (One-Hot)

### Zip Code:
- `zip_code_Suburban`
- `zip_code_Urban`
  - *Rural is the reference*

### Channel:
- `channel_Phone`
- `channel_Web`
  - *Multichannel is the reference*

After encoding, total variables increased from 9 to **12**.

---

## ✅ Statistical Assumption Checks

### Multicollinearity (VIF)
- Most predictors had **VIF < 5**
- Constant term had high VIF (36.36) — not problematic
- `used_discount` and `used_bogo`: VIF ≈ 3.13 (acceptable)

### Sample Size
- Total conversions = **8,961**
- Predictors = 10
- Events Per Variable (EPV) = **896.1**
  - Well above the required minimum (EPV > 10)

---

## 🔍 Key Variables for Heterogeneous Treatment Effects

- **Outcome**: Conversion (binary)
- **Treatment**: Offer vs No Offer (binary)
- **Effect Modifiers**:
  - Geography
  - Recency & Purchase History
  - Past discount/BOGO usage
  - Referral status
  - Preferred channel

---

## ✅ Summary

The cleaned, encoded dataset is **well-prepared** for:
✔ Regression with interaction terms  
✔ Risk-based subgroup analysis  
✔ Heterogeneous treatment effect estimation

With strong sample size, no missing data, low multicollinearity, and a clear treatment definition, this dataset is ideal for modeling how marketing offers work differently across customer segments.

---

