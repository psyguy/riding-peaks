# Summary of Cosinor Amplitude Significance Tests

This document provides a detailed summary of the different statistical methods (A-F) discussed for testing the significance of the cosinor amplitude. For each method, it covers the theoretical basis, the specific implementation in the R script `n1_cosinor_strategy1.R`, and a summary of the key discussion points and clarifications.

---

## General Context

The primary goal is to determine if the 24-hour rhythm's amplitude is significantly different from zero. In the cosinor model, the amplitude `A` is given by `A = sqrt(C^2 + S^2)`, where `C` (cos) and `S` (sin) are the model coefficients. A test of `A > 0` is equivalent to a joint test of whether the pair `(C, S)` is significantly different from the origin `(0, 0)`. The Bayesian approach provides posterior distributions for `C` and `S`, which are approximately bivariate normal. We have a set of posterior draws for `C` and `S`, and the different methods use these draws to calculate a significance test.

---

## Method A: Highest Density Region (HDR)

### Theoretical Basis
The Highest Density Region (HDR) is a concept from Bayesian statistics. A (1-α)% HDR is the smallest possible region of the parameter space that contains (1-α)% of the posterior probability. The test is straightforward: if the null value—in this case, the origin `(0, 0)`—falls outside the 95% HDR of the joint posterior distribution of `(C, S)`, we conclude that the amplitude is significantly different from zero.

### R Implementation (`n1_cosinor_strategy1.R`)
This method is implemented using the `hdrcde` package. The `hdr.test` function is used to directly test if the point `(0, 0)` is contained within the HDR computed from the posterior draws of `C` and `S`.

```R
# Method A: Is (0,0) in the 95% HDR?
# this is the gold standard
method_A_is_null_in_hdr <-
  hdrcde::hdr.test(x = draws_C,
                   y = draws_S,
                   test = c(0, 0))$p.value
is_significant_A <- method_A_is_null_in_hdr < .05
```
The `$p.value` returned here is the probability coverage of the HDR that includes the test point `(0,0)`. If this value is less than 0.95, the point is outside the 95% HDR. The code simplifies this by checking if the returned "p-value" (coverage level) is less than `0.05`, which is a slight misinterpretation of the output but functionally equivalent for a 95% HDR test.

### Discussion & Clarifications
This method is considered a robust and standard way to perform hypothesis tests in a Bayesian framework. There were no major questions or doubts about its validity.

---

## Method B: Mahalanobis Distance

### Theoretical Basis
The Mahalanobis distance measures how many standard deviations a point is away from the mean of a distribution, accounting for covariance between variables. For a point `x` and a distribution with mean `μ` and covariance matrix `Σ`, the squared Mahalanobis distance is `D^2 = (x - μ)' * Σ^(-1) * (x - μ)`.

If the distribution is multivariate normal (as the posterior of `(C,S)` is assumed to be), the squared Mahalanobis distance of a point from the mean follows a chi-squared (`χ²`) distribution with degrees of freedom equal to the number of dimensions (`k`). In our case, `k=2`.

We test the distance of the origin `(0, 0)` from the mean of the posterior `(mean(C), mean(S))`. A large distance suggests `(0, 0)` is an unlikely value.

### R Implementation (`n1_cosinor_strategy1.R`)
The implementation calculates the squared Mahalanobis distance of the origin from the sample mean and covariance of the posterior draws. The p-value is then derived from the `χ²` distribution with 2 degrees of freedom.

```R
# Method B: Mahalanobis distance of origin from posterior
posterior_draws <- cbind(C, S)
posterior_mean <- colMeans(posterior_draws)
posterior_cov <- cov(posterior_draws)
mahal_origin_sq <-
  mahalanobis(c(0, 0), posterior_mean, posterior_cov)
method_B_mahalanobis_pval <-
  pchisq(mahal_origin_sq, df = 2, lower.tail = FALSE)
is_significant_B <- method_B_mahalanobis_pval < .05
```

### Discussion & Clarifications
**Your Question:** How is Mahalanobis distance related to the bivariate normal distribution? Why not just find the density at `(0,0)`? What is the `1 - exp(-0.5 * mahal_origin_sq)` formula?

**Answer:**
1.  **Relation to Bivariate Normal:** The `χ²` distribution relationship holds specifically because the posterior is assumed to be bivariate normal. The Mahalanobis distance is intrinsically linked to the geometry of the normal distribution's probability density ellipses.
2.  **Density at (0,0):** Calculating the probability density at the origin is a valid approach (and is used in some contexts), but it doesn't directly give a p-value. A p-value is the probability of observing a result *at least as extreme* as the one measured. The Mahalanobis distance provides a standardized way to define this "extremeness" by considering all points that are "equally unlikely" as the origin, which form an ellipse. The p-value is the integrated probability outside this ellipse.
3.  **`1 - exp(-0.5 * D^2)`:** This formula is the cumulative distribution function (CDF) for a **Rayleigh distribution**, not a `χ²` distribution. It appears to have been a point of confusion. While the Mahalanobis distance is used as input, the formula itself is incorrect for this test. The correct implementation, which is used in the script, relies on `pchisq`.

---

## Method C: Quadrant & Axis Crossing Probability

### Theoretical Basis
This is a non-standard heuristic method. It assesses significance by checking how "concentrated" the posterior draws are around their mean, away from the axes. The logic is that if the posterior cloud is far from the origin, very few draws should fall into other quadrants or cross the C=0 or S=0 axes.

The method computes three probabilities:
1.  `method_C_origin_quadrant_prob`: The proportion of draws in the same quadrant as the posterior mean.
2.  `method_C_cross_C_axis_prob`: The proportion of draws with a C-value sign opposite to the mean C.
3.  `method_C_cross_S_axis_prob`: The proportion of draws with an S-value sign opposite to the mean S.

### R Implementation (`n1_cosinor_strategy1.R`)
The implementation calculates these three probabilities. The significance test, however, is based on an `OR` condition of the two axis-crossing probabilities.

```R
# ... (calculations for the three probabilities) ...

# Significance test
is_significant_C <-
  (method_C_cross_C_axis_prob < p_thresh) |
  (method_C_cross_S_axis_prob < p_thresh)
```

### Discussion & Clarifications
**Your Question:** What is done with these probabilities? Is it a union or intersection (`AND` or `OR`)?

**Answer:** The implementation uses an `OR` condition (`|`). Significance is declared if the probability of crossing the C-axis *or* the probability of crossing the S-axis is below a threshold (e.g., 0.05).

This is a very unusual testing strategy. It essentially asks, "Is the posterior cloud very clearly on one side of the C-axis OR very clearly on one side of the S-axis?" This does not correspond to a standard joint hypothesis test about the origin and is likely the reason for its low agreement with other methods (as seen in the heatmap). It tests two separate, weaker conditions, whereas other methods perform a single, joint test.

---

## Method D: Separate t-tests

### Theoretical Basis
This method treats the posterior draws of `C` and `S` as two independent samples and performs separate one-sample t-tests to see if their means are significantly different from zero.

### R Implementation (`n1_cosinor_strategy1.R`)
The code runs two separate `t.test`s and combines the results with an `OR` condition.

```R
# Method D: Separate t-tests on C and S
method_D_C_p_value <- t.test(C, mu = 0)$p.value
method_D_S_p_value <- t.test(S, mu = 0)$p.value
is_significant_D <-
  (method_D_C_p_value < .05) | (method_D_S_p_value < .05)
```

### Discussion & Clarifications
**Your Question:** Is there anything wrong with `OR`-ing the p-values?

**Answer:** Yes. This is a classic multiple testing problem. If you perform two independent tests, each at a significance level of α=0.05, the probability of getting at least one false positive (Type I error) is `1 - (1 - 0.05) * (1 - 0.05) ≈ 0.0975`, which is nearly double the desired `α`. This method does not correctly test the joint hypothesis `(C, S) = (0, 0)` and will have an inflated Type I error rate.

---

## Method E: Wald Test

### Theoretical Basis
The Wald test is a general method for testing hypotheses about statistical parameters. For a parameter vector `θ` with estimate `θ_hat` and covariance matrix `V`, the Wald statistic for the hypothesis `θ = θ_0` is `(θ_hat - θ_0)' * V^(-1) * (θ_hat - θ_0)`.

In our case, `θ = [C, S]` and `θ_0 = [0, 0]`. The statistic becomes `[C_mean, S_mean] * cov(C,S)^(-1) * [C_mean, S_mean]'`. This is **mathematically identical** to the squared Mahalanobis distance of the origin from the posterior mean. Like the Mahalanobis distance, the Wald statistic follows a `χ²` distribution with `k=2` degrees of freedom under the null hypothesis.

### R Implementation (`n1_cosinor_strategy1.R`)
The implementation calculates the Wald statistic and its corresponding p-value from the `χ²` distribution.

```R
# Method E: Wald test for C=S=0
# This is equivalent to Mahalanobis distance
wald_stat <- t(posterior_mean) %*% solve(posterior_cov) %*% posterior_mean
method_E_wald_pval <- pchisq(wald_stat, df = 2, lower.tail = FALSE)
is_significant_E <- method_E_wald_pval < .05
```

### Discussion & Clarifications
**Your Question:** Is the Wald test related to the Mahalanobis distance?

**Answer:** Yes, they are equivalent in this application. Both methods use the same quadratic form and result in the same test statistic and p-value. The results of Method B and Method E should be identical.

---

## Method F: Rayleigh Z-score / Amplitude Test

### Theoretical Basis
This method shifts focus from the `(C, S)` pair to the amplitude `A = sqrt(C^2 + S^2)`. Under the null hypothesis that `C` and `S` are independent, normally distributed variables with mean 0 and a common variance `σ²`, the amplitude `A` follows a **Rayleigh distribution**. This method tests whether the observed mean amplitude `A_mean` is plausible under this null Rayleigh distribution.

The probability density function (PDF) of a Rayleigh distribution is `f(x; σ) = (x/σ²) * exp(-x² / (2σ²))`. The survival function (1 - CDF), which gives the probability of observing a value *greater* than `A`, is `P(X > A) = exp(-A² / (2σ²))`. This survival function can be used as a one-tailed p-value.

### R Implementation (`n1_cosinor_strategy1.R`)
The implementation calculates the mean amplitude from the draws. It estimates a pooled variance `sigma_CS_pooled` from the variances of `C` and `S`. It then calculates both a p-value using the survival function and a critical value for `A` at the 95% level.

```R
# Method F: Test on amplitude A = sqrt(C^2+S^2)
A <- sqrt(C ^ 2 + S ^ 2)
A_mean <- mean(A)
sigma_C <- sd(C)
sigma_S <- sd(S)
sigma_CS_pooled <- sqrt((sigma_C ^ 2 + sigma_S ^ 2) / 2)

# p-value based on Rayleigh distribution's survival function
method_F_rayleigh_pval <- exp(-A_mean ^ 2 / (2 * sigma_CS_pooled ^ 2))
is_significant_F <- method_F_rayleigh_pval < .05

# Alternative: compare A_mean to a critical value
a_crit_95 <- sigma_CS_pooled * sqrt(-2 * log(0.05)) # ~2.448σ
is_significant_F2 <- A_mean > a_crit_95
```

### Discussion & Clarifications
**Your Question:** Why are there two approaches (p-value and critical value)? Where does the critical value formula `sigma * sqrt(-2 * log(0.05))` come from?

**Answer:**
1.  **Two Approaches, One Test:** The p-value and critical value approaches are two sides of the same coin. You can either see if your observed statistic (`A_mean`) falls into the critical region (i.e., `A_mean > a_crit_95`), or you can calculate the probability of observing a value as extreme as `A_mean` (the p-value) and see if that probability is less than your alpha level (0.05). They will always give the same conclusion.

2.  **Critical Value Formula:** The formula is derived by inverting the Rayleigh CDF. The CDF is `P(X <= x) = 1 - exp(-x² / (2σ²))`. We want to find the critical value `a_crit` such that `P(X <= a_crit) = 0.95`.
    *   `0.95 = 1 - exp(-a_crit² / (2σ²))`
    *   `exp(-a_crit² / (2σ²)) = 0.05`
    *   `-a_crit² / (2σ²) = log(0.05)`
    *   `a_crit² = -2 * σ² * log(0.05)`
    *   `a_crit = sqrt(-2 * σ² * log(0.05)) = σ * sqrt(-2 * log(0.05))`
    This is exactly the formula used in the code.

---

## Heatmap Analysis (`method_agreement_heatmap.png`)

**Your Question:** Why does the Quadrant method (Method C) have so little overlap with the rest? What do the off-diagonal cells represent?

**Answer:**
1.  **Off-diagonal Cells:** An off-diagonal cell at the intersection of `Method X` (row) and `Method Y` (column) represents the **percentage of agreement** between the two methods across all individuals. Agreement is defined as the proportion of cases where both methods yield the same outcome (i.e., both find a significant result OR both find a non-significant result). The calculation is `mean(is_significant_X == is_significant_Y)`.

2.  **Low Overlap for Method C:** The heatmap shows that Method C has significantly lower agreement with all other methods. This is because, as discussed above, it is based on a different and unorthodox logical foundation. While Methods A, B, E, and F are all robustly testing the joint proposition that `(C, S)` is not `(0, 0)`, Method C is asking a different, weaker question ("Is the posterior mostly to the left/right of the y-axis OR mostly above/below the x-axis?"). Method D is also flawed (inflated Type 1 error) but its logic is closer to the joint test, so it agrees more. The low agreement of Method C strongly suggests it is not a valid test for this particular hypothesis and should not be used.
