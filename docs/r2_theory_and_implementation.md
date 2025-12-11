# R² for Multilevel Cosinor Models: Theory and Implementation

**Date:** 2025-12-08 **Project:** Riding Peaks - Circadian Rhythm Analysis **Focus:** Variance explained by circadian rhythms in multilevel models with gapped ESM data

------------------------------------------------------------------------

## Table of Contents

1.  [Introduction](#introduction)
2.  [Model Specification](#model-specification)
3.  [The Challenge: Gapped ESM Data](#the-challenge-gapped-esm-data)
4.  [Method 1: Distribution of Individual R²](#method-1-distribution-of-individual-r²)
5.  [Method 2: Simplified Marginal & Conditional R²](#method-2-simplified-marginal--conditional-r²)
6.  [Method 3: Exact R² for Gapped Data](#method-3-exact-r²-for-gapped-data)
7.  [Random Sigma: The Log-Normal Complication](#random-sigma-the-log-normal-complication)
8.  [Conservative vs Alternative Views](#conservative-vs-alternative-views)
9.  [Summary Comparison](#summary-comparison)
10. [Implementation Recommendations](#implementation-recommendations)

------------------------------------------------------------------------

## Introduction {#introduction}

Quantifying how much variance in ecological momentary assessment (EMA/ESM) data is explained by circadian rhythms is essential for understanding the practical significance of detected rhythms. While statistical significance (p \< 0.05) tells us whether an effect is reliably different from zero, R² tells us how much of the outcome variability is accounted for by the circadian pattern.

This document outlines three approaches to computing R² for multilevel cosinor models, with particular attention to:

1.  **Gapped data structure:** ESM data typically has no observations during nighttime sleep hours
2.  **Random residual variance:** Individual differences in within-person variability (σᵢ)
3.  **Marginal vs Conditional interpretations:** Population-average vs individual-specific variance explained

------------------------------------------------------------------------

## Model Specification {#model-specification}

Our multilevel cosinor model follows this specification:

```         
Level 1 (within-person):
  yᵢⱼ = Mᵢ + Cᵢ·cos(2πtᵢⱼ/24) + Sᵢ·sin(2πtᵢⱼ/24) + εᵢⱼ
  εᵢⱼ ~ Normal(0, σᵢ²)

Level 2 (between-person):
  [Mᵢ, Cᵢ, Sᵢ]' ~ MVN([γₘ, γ_c, γ_s]', Σ_random)
  log(σᵢ) = ρ + ϱᵢ,  where ϱᵢ ~ Normal(0, τ²)
```

In brms syntax:

``` r
brmsformula(
  y ~ 1 + co + si + (1 + co + si | i | id),
  sigma ~ 1 + (1 | i | id)
)
```

### Key Parameters

| Symbol | Description |
|----------------------------|--------------------------------------------|
| γₘ, γ_c, γ_s | Population-average mesor, cosine, and sine coefficients |
| Mᵢ, Cᵢ, Sᵢ | Individual-specific mesor, cosine, and sine coefficients |
| Aᵢ = √(Cᵢ² + Sᵢ²) | Individual amplitude |
| A_pop = √(γ_c² + γ_s²) | Population-average amplitude |
| σᵢ | Individual-specific residual SD |
| ρ | Log-scale intercept for sigma |
| τ | SD of random sigma deviations |
| Σ_random | 4×4 covariance matrix for random effects (M, C, S, log(σ)) |

------------------------------------------------------------------------

## The Challenge: Gapped ESM Data {#the-challenge-gapped-esm-data}

### Why This Matters

Standard R² formulas for cosinor models assume observations are uniformly distributed across the 24-hour cycle, or at least that:

$$E[\cos(2\pi t/24)] = 0 \quad \text{and} \quad E[\sin(2\pi t/24)] = 0$$

This holds when observations span the full cycle or are symmetrically distributed. However, ESM data typically has a **night gap** (e.g., no observations between 23:00-07:00), which means:

$$E[\cos(2\pi t/24)] \neq 0 \quad \text{and} \quad E[\sin(2\pi t/24)] \neq 0$$

### Empirical Predictor Properties

For our data with observations only during waking hours (\~07:00-23:00), we can compute:

``` r
# Empirical expectations (example values)
E_cos <- mean(cos(2*pi*t/24))  # ≈ -0.15 (negative because more evening obs)
E_sin <- mean(sin(2*pi*t/24))  # ≈ +0.10 (positive because more afternoon obs)

# Empirical variances
Var_cos <- var(cos(2*pi*t/24))  # ≈ 0.35 (less than 0.5)
Var_sin <- var(sin(2*pi*t/24))  # ≈ 0.40 (less than 0.5)

# Covariance (non-zero with asymmetric sampling)
Cov_cos_sin <- cov(cos(2*pi*t/24), sin(2*pi*t/24))  # ≈ 0.05
```

These empirical values are needed for exact variance calculations.

------------------------------------------------------------------------

## Method 1: Distribution of Individual R² {#method-1-distribution-of-individual-r²}

### Concept

Compute R² separately for each individual, then examine the distribution of these values across the sample. This answers: *"What proportion of each person's **within-person** outcome variability is explained by their own circadian rhythm?"*

### Formula

For each individual i:

$$R^2_i = \frac{A_i^2 / 2}{A_i^2 / 2 + \sigma_i^2}$$

Where:
- $A_i^2 / 2$ = variance of the sinusoidal component for individual i
- $\sigma_i^2$ = residual variance for individual i

### Why No Mesor Variance?

Note that the mesor Mᵢ does **not** appear in this formula. This is intentional:

- **The mesor is constant within each person** - it represents that person's average level
- For within-person R², we ask: "Of this person's moment-to-moment variability, how much is due to circadian rhythm?"
- Since Mᵢ doesn't vary across time for person i, Var(Mᵢ) = 0 within that person
- The mesor only contributes to **between-person variance**, which is relevant for population-level R² (Methods 2 & 3)

Mathematically, for person i:
$$\text{Var}(y_{ij} | \text{person } i) = \text{Var}(M_i + C_i \cos + S_i \sin + \varepsilon_{ij}) = \text{Var}(\text{circadian}) + \text{Var}(\varepsilon)$$

The Mᵢ, Cᵢ, and Sᵢ are all fixed for person i (they're that person's parameters), so only the time-varying components contribute to within-person variance.

### Derivation

The sinusoidal component $C_i \cos(\omega t) + S_i \sin(\omega t)$ has variance:

$$\text{Var}[C_i \cos(\omega t) + S_i \sin(\omega t)] = C_i^2 \text{Var}[\cos] + S_i^2 \text{Var}[\sin] + 2 C_i S_i \text{Cov}[\cos, \sin]$$

For uniformly distributed time (or symmetric sampling): - $\text{Var}[\cos] = \text{Var}[\sin] = 1/2$ - $\text{Cov}[\cos, \sin] = 0$

Therefore: $\text{Var}[\text{sinusoid}] = (C_i^2 + S_i^2)/2 = A_i^2/2$

### Advantages

1.  **Interpretable at individual level:** Each person gets their own R²
2.  **Shows heterogeneity:** Distribution reveals range of effect sizes
3.  **No averaging assumptions:** Each calculation uses that person's parameters
4.  **Works with random sigma:** Naturally incorporates σᵢ

### Limitations

1.  **Assumes symmetric time sampling** within each individual
2.  **No population-level summary** (though can compute mean/median)
3.  **Posterior uncertainty:** Should compute distribution over posterior draws

### Implementation

``` r
# For each posterior draw and each individual
R2_individual <- (A_i^2 / 2) / (A_i^2 / 2 + sigma_i^2)

# Summarize across individuals
R2_median <- median(R2_individual)
R2_mean <- mean(R2_individual)
R2_range <- range(R2_individual)
```

------------------------------------------------------------------------

## Method 2: Simplified Marginal & Conditional R²

### Concept

Following Nakagawa & Schielzeth (2013), decompose total variance into components explained by fixed effects (marginal) and fixed + random effects (conditional).

### Formulas (Simplified)

**Marginal R² (population-average):**

$$R^2_m = \frac{A_{pop}^2 / 2}{A_{pop}^2 / 2 + \text{Var}(u_C) / 2 + \text{Var}(u_S) / 2 + \bar{\sigma}^2}$$

**Conditional R² (including random effects):**

$$R^2_c = \frac{A_{pop}^2 / 2 + \text{Var}(u_C) / 2 + \text{Var}(u_S) / 2}{A_{pop}^2 / 2 + \text{Var}(u_C) / 2 + \text{Var}(u_S) / 2 + \bar{\sigma}^2}$$

Where: - $A_{pop} = \sqrt{\gamma_c^2 + \gamma_s^2}$ = population-average amplitude - $\text{Var}(u_C)$, $\text{Var}(u_S)$ = variance of random slopes for cos and sin - $\bar{\sigma}^2$ = average residual variance

### Critical Assumption: Why This Is Wrong for Gapped Data

The simplified formulas assume $E[\cos] = E[\sin] = 0$, which allows:

$$\text{Var}[C \cdot \cos(t)] = C^2 \cdot \text{Var}[\cos] = C^2 / 2$$

But when $E[\cos] \neq 0$:

$$\text{Var}[C \cdot \cos(t)] = E[C]^2 \cdot \text{Var}[\cos] + \text{Var}[C] \cdot E[\cos^2]$$

This changes the variance decomposition substantially.

### When Method 2 Is Appropriate

-   Full 24-hour coverage (e.g., continuous actigraphy)
-   Symmetric sampling around the clock
-   Quick approximation when exact values aren't critical

### When Method 2 Is Inappropriate

-   ESM/EMA with night gap ✗
-   Unequal sampling across times of day ✗
-   High-stakes quantitative comparisons ✗

------------------------------------------------------------------------

## Method 3: Exact R² for Gapped Data {#method-3-exact-r²-for-gapped-data}

### Concept

Use the actual empirical properties of the predictor variables (cos, sin) as observed in the data, rather than theoretical values assuming uniform time coverage.

### Notation

Define empirical predictor properties from the observed time points:

| Symbol   | Definition                    | Full coverage | Gapped data |
|----------|-------------------------------|---------------|-------------|
| $\mu_c$  | $E[\cos(2\pi t/24)]$          | 0             | ≠ 0         |
| $\mu_s$  | $E[\sin(2\pi t/24)]$          | 0             | ≠ 0         |
| $v_c$    | $\text{Var}[\cos(2\pi t/24)]$ | 0.5           | \< 0.5      |
| $v_s$    | $\text{Var}[\sin(2\pi t/24)]$ | 0.5           | \< 0.5      |
| $c_{cs}$ | $\text{Cov}[\cos, \sin]$      | 0             | ≠ 0         |

### Variance of Fixed Effects

$$\text{Var}(\text{Fixed}) = \gamma_c^2 v_c + \gamma_s^2 v_s + 2\gamma_c \gamma_s c_{cs}$$

Note: This simplifies to $(\gamma_c^2 + \gamma_s^2)/2 = A_{pop}^2/2$ only when $v_c = v_s = 0.5$ and $c_{cs} = 0$.

### Variance of Random Effects

For the random slopes:

$$\text{Var}(\text{Random}_C) = \text{Var}(u_{Ci}) \cdot E[\cos^2] + E[u_{Ci}]^2 \cdot v_c$$

Since $E[u_{Ci}] = 0$ (random effects are mean-zero):

$$\text{Var}(\text{Random}_C) = \text{Var}(u_{Ci}) \cdot E[\cos^2] = \text{Var}(u_{Ci}) \cdot (v_c + \mu_c^2)$$

Similarly for sine:

$$\text{Var}(\text{Random}_S) = \text{Var}(u_{Si}) \cdot (v_s + \mu_s^2)$$

Covariance term (if random slopes are correlated):

$$\text{Cov}(\text{Random}_{C}, \text{Random}_{S}) = \text{Cov}(u_{Ci}, u_{Si}) \cdot E[\cos \cdot \sin]$$

Where $E[\cos \cdot \sin] = c_{cs} + \mu_c \mu_s$.

### Complete Variance Decomposition

$$\text{Var}(\text{Total}) = \text{Var}(\text{Fixed}) + \text{Var}(\text{Random}) + E[\sigma_i^2]$$

Where: $$\text{Var}(\text{Random}) = \text{Var}(u_M) + \text{Var}(u_{Ci})(v_c + \mu_c^2) + \text{Var}(u_{Si})(v_s + \mu_s^2) + 2\text{Cov}(u_{Ci}, u_{Si})(c_{cs} + \mu_c \mu_s)$$

### Exact Marginal R²

$$R^2_{m,\text{exact}} = \frac{\gamma_c^2 v_c + \gamma_s^2 v_s + 2\gamma_c \gamma_s c_{cs}}{\text{Var}(\text{Total})}$$

### Exact Conditional R²

$$R^2_{c,\text{exact}} = \frac{\text{Var}(\text{Fixed}) + \text{Var}(\text{Random})}{\text{Var}(\text{Total})}$$

------------------------------------------------------------------------

## Random Sigma: The Log-Normal Complication {#random-sigma-the-log-normal-complication}

### The Issue

Our model includes random residual variance via a log-link:

$$\log(\sigma_i) = \rho + \varrho_i, \quad \varrho_i \sim N(0, \tau^2)$$

Therefore: $\sigma_i = e^{\rho + \varrho_i}$

### Expected Value of σᵢ²

This is critical for R² calculations. The naive approach would use:

$$\bar{\sigma}^2 \stackrel{?}{=} e^{2\rho} \quad \text{(WRONG)}$$

But since $\varrho_i$ is random:

$$E[\sigma_i^2] = E[e^{2(\rho + \varrho_i)}] = e^{2\rho} \cdot E[e^{2\varrho_i}]$$

For a log-normal distribution, $E[e^{k\varrho}] = e^{k^2\tau^2/2}$. With $k=2$:

$$E[\sigma_i^2] = e^{2\rho} \cdot e^{2\tau^2} = e^{2\rho + 2\tau^2}$$

### Comparison

| Approach           | Formula               | Value if ρ=-1, τ=0.3 |
|--------------------|-----------------------|----------------------|
| Naive (ignoring τ) | $e^{2\rho}$           | 0.135                |
| Correct            | $e^{2\rho + 2\tau^2}$ | 0.162                |

The difference is a factor of $e^{2\tau^2}$. For τ = 0.3, this is $e^{0.18} ≈ 1.20$, a 20% underestimate of residual variance.

### Implication for R²

Using the naive approach **overestimates R²** because residual variance is underestimated in the denominator.

### Implementation

``` r
# Extract from brms model
rho <- fixef(model, pars = "sigma_Intercept")[1]
tau <- sd(ranef(model)$id[,,"sigma_Intercept"])

# Correct expected residual variance
E_sigma_sq <- exp(2*rho + 2*tau^2)
```

------------------------------------------------------------------------

## Conservative vs Alternative Views {#conservative-vs-alternative-views}

Two philosophical approaches exist for handling variance decomposition in multilevel models with heterogeneous residual variance.

### Which R² Does This Apply To?

**Important:** The alternative view (reporting individual heterogeneity in σᵢ) applies differently to marginal vs conditional R²:

| R² Type | Alternative View Applicable? | Reason |
|---------------|------------------------------------------|---------------|
| **Marginal R²** | No | Fixed effects are population-level by definition. There is no meaningful "individual marginal R²" because γ_c and γ_s are the same for everyone. |
| **Conditional R²** | Yes | Individual-level predictions use person-specific Aᵢ and σᵢ. The question "how much of *this person's* variance is explained?" has a meaningful answer. |

Therefore: - **Marginal R²**: Always use the conservative approach with E\[σᵢ²\] - **Conditional R²**: Can report both conservative (population average) and alternative (individual distribution)

### Conservative View

**Philosophy:** Residual variance heterogeneity is "noise" that should be averaged over. Individual differences in σᵢ represent measurement precision differences, not substantively interesting variance.

**Theoretical justification:** R² is fundamentally about mean prediction accuracy. The question "how much variance does the circadian model explain?" should be answered with E\[σᵢ²\], which represents the expected residual variance across the population. The *variance* of σᵢ² (i.e., Var\[σᵢ²\]) tells us about heterogeneity in prediction accuracy, but doesn't change the average.

**Implementation:** - Use $E[\sigma_i^2] = e^{2\rho + 2\tau^2}$ as single residual variance - Compute one R² for the population

**Formulas:**

Marginal R² (always use this approach): $$R^2_m = \frac{\text{Var(Fixed)}}{\text{Var(Fixed)} + \text{Var(Random)} + e^{2\rho + 2\tau^2}}$$

Conditional R² (conservative population average): $$R^2_c = \frac{\text{Var(Fixed)} + \text{Var(Random)}}{\text{Var(Fixed)} + \text{Var(Random)} + e^{2\rho + 2\tau^2}}$$

**Advantages:** - Single summary statistic - Accounts for average individual differences in precision - Comparable across studies

**Disadvantages:** - Masks heterogeneity in explained variance - Treats σᵢ differences as nuisance

### Alternative View (Conditional R² Only)

**Philosophy:** Individual differences in residual variance are substantively meaningful. Someone with low σᵢ may have a more regular lifestyle, better measurement, or genuinely less stochastic behavior. Conditional R² should reflect this.

**Why this applies only to conditional R²:** The individual R² formula

$$R^2_i = \frac{A_i^2 / 2}{A_i^2 / 2 + \sigma_i^2}$$

is essentially an individual-level *conditional* R². It uses each person's: - Amplitude Aᵢ (their fitted circadian curve = fixed + random effects) - Residual variance σᵢ² (their unexplained variability)

This answers: *"For person i specifically, what proportion of their outcome variance is explained by their circadian rhythm?"*

There is no analogous individual marginal R² because the fixed effects (γ_c, γ_s) are population parameters, not individual parameters.

**Implementation Options:**

#### Option A: Report Distribution of Individual R²ᵢ

Use Method 1 (individual R²) and report the distribution: - Median R² - IQR of R² - Proportion with R² \> 0.10 (or other threshold)

#### Option B: Stratified R²

Compute separate conditional R² for subgroups based on σᵢ: - High-variance individuals (σᵢ \> median) - Low-variance individuals (σᵢ \< median)

#### Option C: Weighted R²

Weight individual R²ᵢ by some measure of reliability or sample size.

**Advantages:** - Captures meaningful individual differences - More nuanced interpretation - Identifies who benefits most from circadian modeling

**Disadvantages:** - No single summary - More complex to report - May conflate measurement precision with substantive differences

### Recommendation

**For Marginal R²:** - Report single value using E\[σᵢ²\] = exp(2ρ + 2τ²)

**For Conditional R²:** Report both: 1. **Conservative population conditional R²** using E\[σᵢ²\] for comparability 2. **Distribution of individual R²ᵢ** (Method 1) to show heterogeneity

This provides both summary statistics for comparison and insight into individual differences in how well circadian rhythms predict behavior.

------------------------------------------------------------------------

## Summary Comparison {#summary-comparison}

| Aspect | Method 1 | Method 2 | Method 3 |
|------------------|------------------|------------------|------------------|
| **Level** | Individual | Population | Population |
| **Gapped data** | Works (with caveat) | ✗ Biased | ✓ Correct |
| **Random σᵢ** | Natural fit | Needs correction | Needs correction |
| **Output** | Distribution | Single values | Single values |
| **Interpretation** | Person-specific R² | Marginal/Conditional | Marginal/Conditional |
| **Complexity** | Low | Low | Medium |

### When to Use Each

| Scenario                    | Recommended Method               |
|-----------------------------|----------------------------------|
| Individual-level questions  | Method 1                         |
| Full 24-hour data           | Method 2 (simple) or Method 3    |
| Gapped ESM data             | Method 3 (exact)                 |
| Comparison across studies   | Method 3 with conservative σ²    |
| Understanding heterogeneity | Method 1 + conservative Method 3 |

------------------------------------------------------------------------

## Implementation Recommendations {#implementation-recommendations}

### Recommended Workflow

1.  **Compute empirical predictor properties** from observed time points:

    ``` r
    mu_c <- mean(cos(2*pi*t/24))
    mu_s <- mean(sin(2*pi*t/24))
    v_c <- var(cos(2*pi*t/24))
    v_s <- var(sin(2*pi*t/24))
    c_cs <- cov(cos(2*pi*t/24), sin(2*pi*t/24))
    E_cos2 <- mean(cos(2*pi*t/24)^2)  # = v_c + mu_c^2
    E_sin2 <- mean(sin(2*pi*t/24)^2)  # = v_s + mu_s^2
    E_cossin <- mean(cos(2*pi*t/24) * sin(2*pi*t/24))  # = c_cs + mu_c*mu_s
    ```

2.  **Extract model parameters** from posterior:

    ``` r
    # Fixed effects
    gamma_c <- fixef(model)["co", "Estimate"]
    gamma_s <- fixef(model)["si", "Estimate"]

    # Random effect variances (from Sigma_random)
    var_uM <- VarCorr(model)$id$sd["Intercept"]^2
    var_uC <- VarCorr(model)$id$sd["co"]^2
    var_uS <- VarCorr(model)$id$sd["si"]^2
    cov_uCS <- VarCorr(model)$id$cor["co","si"] *
               VarCorr(model)$id$sd["co"] * VarCorr(model)$id$sd["si"]

    # Residual variance (corrected for random sigma)
    rho <- fixef(model, pars = "sigma")["Intercept", "Estimate"]
    tau <- VarCorr(model)$id$sd["sigma_Intercept"]
    E_sigma_sq <- exp(2*rho + 2*tau^2)
    ```

3.  **Compute variance components:**

    ``` r
    # Fixed effects variance
    Var_Fixed <- gamma_c^2 * v_c + gamma_s^2 * v_s + 2*gamma_c*gamma_s*c_cs

    # Random effects variance
    Var_Random <- var_uM + var_uC * E_cos2 + var_uS * E_sin2 +
                  2 * cov_uCS * E_cossin

    # Total variance
    Var_Total <- Var_Fixed + Var_Random + E_sigma_sq
    ```

4.  **Compute R² values:**

    ``` r
    R2_marginal <- Var_Fixed / Var_Total
    R2_conditional <- (Var_Fixed + Var_Random) / Var_Total
    ```

5.  **Compute individual R² distribution** (Method 1):

    ``` r
    # For each individual from posterior draws
    R2_individual <- (A_i^2 / 2) / (A_i^2 / 2 + sigma_i^2)
    ```

6.  **Report:**

    -   Marginal R² (exact): X% \[95% CI\]
    -   Conditional R² (exact): Y% \[95% CI\]
    -   Individual R² median: Z% \[IQR: A%-B%\]
    -   Proportion with R² \> 10%: N/total

### Comparison with Naive/Simplified Approach

Also compute Method 2 (simplified) for comparison:

``` r
# Simplified (assumes E[cos]=E[sin]=0)
Var_Fixed_simple <- (gamma_c^2 + gamma_s^2) / 2
Var_Random_simple <- var_uM + (var_uC + var_uS) / 2

R2_marginal_simple <- Var_Fixed_simple / (Var_Fixed_simple + Var_Random_simple + E_sigma_sq)
R2_conditional_simple <- (Var_Fixed_simple + Var_Random_simple) /
                         (Var_Fixed_simple + Var_Random_simple + E_sigma_sq)
```

Report the difference to quantify the bias from ignoring gapped structure.

------------------------------------------------------------------------

## References

-   Nakagawa, S., & Schielzeth, H. (2013). A general and simple method for obtaining R² from generalized linear mixed-effects models. *Methods in Ecology and Evolution*, 4(2), 133-142.

-   Johnson, P. C. (2014). Extension of Nakagawa & Schielzeth's R²_GLMM to random slopes models. *Methods in Ecology and Evolution*, 5(9), 944-946.

-   Mardia, K. V., & Jupp, P. E. (2000). *Directional Statistics*. Wiley. (For circular/cosinor variance derivations)

------------------------------------------------------------------------

## Appendix: Full Derivation of Sinusoidal Variance

For completeness, here is the full derivation of the variance of a sinusoidal component.

Let $f(t) = C \cos(\omega t) + S \sin(\omega t)$ where C and S may be random.

### Case 1: C and S are fixed constants

$$\text{Var}[f(t)] = C^2 \text{Var}[\cos(\omega t)] + S^2 \text{Var}[\sin(\omega t)] + 2CS \text{Cov}[\cos, \sin]$$

For uniform time: $= C^2 \cdot 0.5 + S^2 \cdot 0.5 + 0 = (C^2 + S^2)/2 = A^2/2$

For gapped time: $= C^2 v_c + S^2 v_s + 2CS \cdot c_{cs}$

### Case 2: C and S are random with E\[C\] = γ_c, E\[S\] = γ_s

Using the law of total variance:

$$\text{Var}[C \cos(\omega t)] = E[C]^2 \text{Var}[\cos] + \text{Var}[C] E[\cos^2]$$

This accounts for both the variance of time (Var\[cos\]) and the variance of the coefficient (Var\[C\]).

For the full random component:

$$\text{Var}[\text{Random}] = \text{Var}[u_C] E[\cos^2] + \text{Var}[u_S] E[\sin^2] + 2 \text{Cov}[u_C, u_S] E[\cos \cdot \sin]$$

This is used in Method 3 for the exact random effects variance.

------------------------------------------------------------------------

**End of Document**