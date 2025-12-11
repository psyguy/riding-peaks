# A Guide to Calculating R-squared for Multilevel Cosinor Models

This document provides a guide to calculating and interpreting different versions of R-squared (R²) for multilevel cosinor models, particularly in the context of intensive longitudinal data (ILD) from Experience Sampling Method (ESM) studies.

## Background

In multilevel modeling (MLM), there isn't a single, universally agreed-upon R². The presence of multiple variance components (fixed effects, random effects, residuals) makes its definition ambiguous. However, we can calculate several useful metrics to understand the practical significance of our model. The reviewer's request for an "overarching R²" points to a need for a measure of variance explained at the population level.

We will cover three main approaches:
1.  **The Distribution of Individual R²s:** A direct and intuitive measure of practical significance for each person.
2.  **Simplified Marginal & Conditional R²:** A common population-level metric that makes a key simplifying assumption.
3.  **Exact Marginal & Conditional R²:** A more complex but accurate version for typical ESM data.

Our multilevel cosinor model is:

`y_it = (M + m_i) + (C + c_i)cos(t) + (S + s_i)sin(t) + e_it`

where `M`, `C`, and `S` are fixed effects, and `m_i`, `c_i`, and `s_i` are the corresponding random effects for individual `i`.

---

## Method 1: The Distribution of Individual R²s

### Theory

This method calculates the proportion of within-person variance that is explained by the individual-specific cyclic trend. It directly answers the question: "For a given person, how much of their data's variability is captured by their personal cycle?"

The variance of the sinusoidal signal for person `i` is `A_i² / 2`. The total variance for that person is the sum of the signal variance and their residual variance, `σ_i²`.

The formula for an individual's R² is:

`R²_i = (A_i² / 2) / (A_i² / 2 + σ_i²)`

In a Bayesian context, we calculate this for each posterior draw, yielding a full posterior distribution of R² for each person.

### Population-Level Interpretation

The "overarching" or population-level measure is the **distribution of these individual R² values**. We can summarize this distribution to make a statement about the population.

*   **Descriptive Statistics:** Report the median and interquartile range (IQR) of the individual R² medians.
*   **Proportions:** Report the percentage of individuals whose R² falls into meaningful bins (e.g., R² < 0.01, 0.01 < R² < 0.10, R² > 0.10).
*   **Visualization:** A histogram or density plot of the individual R² medians.

### R Implementation

```r
# d_draws contains posterior draws for each person (id, amp, sigma)
r2_individual_draws <- d_draws %>%
  mutate(
    R2_individual = (amp^2 / 2) / (amp^2 / 2 + sigma^2)
  )

# Summarize to get a single R² value per person (e.g., the median)
r2_individual_summary <- r2_individual_draws %>%
  group_by(id) %>%
  summarise(R2_median = median(R2_individual))

# Plot the distribution of these medians
ggplot(r2_individual_summary, aes(x = R2_median)) +
  geom_histogram(bins = 30)
```

---

## Method 2: Simplified Marginal & Conditional R²

### Theory

This approach, based on the work of Nakagawa & Schielzeth, partitions the total model-implied variance into components.

*   **Marginal R² (`R²_m`):** Variance explained by **fixed effects only** (the population-average cycle).
*   **Conditional R² (`R²_c`):** Variance explained by **fixed and random effects** (the individual-specific cycles).

This method makes a **critical simplifying assumption**: the data are collected over a full, balanced 24-hour cycle. This implies that the mean of the predictors `cos(t)` and `sin(t)` is zero, and their variance is `1/2`. As a result, the covariance terms between random effects drop out of the calculation.

The variance components are:
*   `Var_fixed = (C² + S²) / 2`
*   `Var_random = Var(m_i) + Var(c_i)/2 + Var(s_i)/2`
*   `Var_residual = E[σ_i²]` (average residual variance across the population)

`R²_m = Var_fixed / (Var_fixed + Var_random + Var_residual)`
`R²_c = (Var_fixed + Var_random) / (Var_fixed + Var_random + Var_residual)`

**This method is an approximation for typical ESM data, which has gaps (e.g., during sleep).**

### R Implementation

See the `for` loop implementation in `scripts/calculate_r2_metrics.R` under "Method 2". It iterates through each MCMC draw, calculates the variance components using the simplified formulas, and computes the R² values.

---

## Method 3: Exact Marginal & Conditional R² for Gapped Data

### Theory

This is the most accurate method for your data. It correctly accounts for the fact that because of the night gap, the mean of your predictors `cos(t)` and `sin(t)` is **not zero** over the observed data. This means the covariance terms between the random effects **do not drop out**.

Let `U = cos(t)` and `V = sin(t)`. We first calculate their empirical properties from the raw data: `E[U]`, `E[V]`, `E[U²]`, `E[V²]`, `E[UV]`, `Var(U)`, `Var(V)`, and `Cov(UV)`.

The exact variance components are:

1.  **`Var_fixed`**: The variance of the fixed-effect prediction. This uses the **empirical** variance of the predictors.
    `Var_fixed = C²*Var(U) + S²*Var(V) + 2*C*S*Cov(U,V)`

2.  **`Var_random`**: The variance of the random-effect part of the prediction. This uses the **model-estimated** variances and covariances of the random effects, combined with the **empirical** properties of the predictors.
    `Var_random = Var(m_i) + Var(c_i)E[U²] + Var(s_i)E[V²] + 2Cov(m_i, c_i)E[U] + 2Cov(m_i, s_i)E[V] + 2Cov(c_i, s_i)E[UV]`

3.  **`Var_residual`**: The average model-implied residual variance (same as in Method 2).
    `Var_residual = E[σ_i²]`

The R² formulas are the same as before, but using these exact variance components.

### Why is this the correct approach?

*   It correctly partitions the **model-implied variance**.
*   It respects the fact that the predictors (`co`, `si`) are fixed aspects of the study design, so their empirical variance is used for `Var_fixed`.
*   It correctly uses the model-estimated random effect structure for `Var_random`, including the crucial covariance terms that are non-zero due to the gapped data.

### R Implementation

See the `for` loop implementation in `scripts/calculate_r2_metrics.R` under "Method 3". It first calculates the `predictor_props` from the raw data, then iterates through each MCMC draw to calculate the variance components using the full, exact formulas.

---

## Summary and Recommendation

| Method                      | Pros                                                              | Cons                                                              | Best For                                                              |
| --------------------------- | ----------------------------------------------------------------- | ----------------------------------------------------------------- | --------------------------------------------------------------------- |
| **1. Individual R²s**       | Intuitive, easy to calculate, directly measures effect for each person. | Not a single population number; requires summarizing a distribution. | Visualizing heterogeneity and reporting on the range of practical effects. |
| **2. Simplified Pop. R²**   | Simpler formula, common in literature.                            | **Inaccurate for gapped data.** Makes an assumption that is violated. | Didactic purposes or cases with truly balanced, full-cycle data.      |
| **3. Exact Pop. R²**        | **Theoretically correct for gapped ESM data.** Most defensible.     | More complex formula, requires careful implementation.            | Reporting final, accurate population-level R² values in a publication. |

For your paper, the most robust approach is to report the **Exact Marginal and Conditional R² (Method 3)** as your primary population-level metrics. You can supplement this with a summary of the **Distribution of Individual R²s (Method 1)** to illustrate the between-person heterogeneity in practical significance.

```

<!--
[PROMPT_SUGGESTION]Can you help me rewrite the "Empirical findings" section in `amp-text.tex` to incorporate these new R² results?[/PROMPT_SUGGESTION]
[PROMPT_SUGGESTION]How do the results from the exact R² calculation compare to the simplified one? Can you explain the differences in the plot?[/PROMPT_SUGGESTION]
