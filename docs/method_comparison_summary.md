# Summary of n-of-1 Cosinor Analysis Comparison Methods

This document provides a detailed explanation of the different statistical methods used to test the significance of the amplitude in n-of-1 cosinor analysis. It summarizes the theoretical basis for each method, details its implementation in the R scripts, and addresses specific questions raised during our analysis of the results.

## User's Summary for Response Letter

As a starting point, here is the summary drafted for the response letter:

> For method B, we used Mahalanobis distance between the origin and the implied bivariate normal posterior distribution, which gives identical results as method E. Method C was more unorthodox than the HDR method, as to our understanding, it doesn’t correspond to an account involving probability or density at the origin itself. For method D we expect an inflated Type-1 error due to multiple testing. And method F was not a suitable test, as the distribution of the square root of the sum of squares of two normally-distributed variables is Rayleigh-distributed (not Gaussian) and does not have a zero mean.

## The Core Problem: Is the Amplitude Zero?

In cosinor analysis, the rhythm is modeled as:
`y(t) = M + A * cos(2πt/τ - φ)`

This can be re-parameterized into a linear model using sine and cosine components:
`y(t) = M + β * cos(2πt/τ) + γ * sin(2πt/τ)`

Where:
- Amplitude `A = sqrt(β² + γ²)`
- Acrophase `φ = atan2(γ, β)`

The null hypothesis `H₀: A = 0` is equivalent to the joint null hypothesis `H₀: β = 0 AND γ = 0`. Our goal is to test this hypothesis. The posterior draws of `β` and `γ` from the Bayesian model form a bivariate (typically normal) distribution. The challenge is to determine if the origin (0,0) is a credible value for the pair `(β, γ)`.

---

## Method A: Highest Density Region (HDR)

- **Theory:** The Highest Density Region contains the `x%` most probable values of a distribution. A `(1-α) * 100%` HDR is the smallest region containing `(1-α) * 100%` of the probability mass. To test the null hypothesis, we construct, for example, a 95% HDR of the joint posterior distribution of `(β, γ)` and check if the origin (0,0) falls inside it. If the origin is outside the HDR, we reject the null hypothesis.
- **Implementation:** This was the primary method used in the multilevel model (`brms`) and serves as the "gold standard" for comparison. It is computationally intensive for n-of-1 models, which is why other methods were explored.

---

## Method B: Mahalanobis Distance

- **Theory:** The Mahalanobis distance measures the distance between a point and the center of a distribution, scaled by the covariance of that distribution. For a point `(β, γ)` and a bivariate normal distribution centered at `(μ_β, μ_γ)` with covariance matrix `Σ`, the squared Mahalanobis distance is `D² = [(β, γ) - (μ_β, μ_γ)]ᵀ Σ⁻¹ [(β, γ) - (μ_β, μ_γ)]`.
- When testing if the posterior distribution is consistent with the origin (0,0), we calculate the distance from the posterior's mean `(β̄, γ̄)` to the origin. This value, `mahal_origin_sq`, follows a Chi-squared distribution with 2 degrees of freedom (`χ²(2)`).
- The p-value is the probability of observing a distance as large or larger than the one calculated, i.e., `P(χ²(2) > mahal_origin_sq)`. This can be calculated in R as `pchisq(mahal_origin_sq, df = 2, lower.tail = FALSE)`.
- **Relationship to Density:** This p-value is equivalent to `exp(-0.5 * mahal_origin_sq)`. This formula comes from the probability density function of the `χ²(2)` distribution. Calculating the probability at the tail is mathematically equivalent to evaluating the density at the origin, just expressed differently. So, this method is directly related to finding the density at (0,0).
- **Implementation (`n1_cosinor_strategy1.R`):**
    ```R
    # Calculate squared Mahalanobis distance from origin to posterior mean
    mahal_origin_sq <- mahalanobis(
      x = c(0, 0),
      center = c(beta_mean, gamma_mean),
      cov = cov_matrix
    )
    # Convert to p-value
    method_B_mahal_pval <- 1 - pchisq(mahal_origin_sq, df = 2)
    ```
    The calculation `1 - pchisq(...)` is equivalent to `pchisq(..., lower.tail=FALSE)`.

---

## Method C: Quadrant & Axis Crossing Probabilities

- **Theory:** This is a geometric heuristic based on the location of the posterior samples of `(β, γ)`. It does not directly calculate density at the origin.
    1.  **`method_C_origin_quadrant_prob`**: Calculates the proportion of posterior draws in the same quadrant as the mean of the draws. A significant amplitude would mean the posterior is concentrated in one quadrant, far from the origin. The test checks if this proportion is high (e.g., > 97.5%).
    2.  **`method_C_cross_S_axis_prob` & `method_C_cross_C_axis_prob`**: These calculate the probability that the confidence ellipse of the posterior crosses the sine (`β=0`) or cosine (`γ=0`) axes. If the posterior is far from the origin, it is less likely to cross the axes. The test checks if this probability is low.
- **Implementation (`n1_cosinor_strategy1.R`):**
    - The significance of these three metrics is determined.
    - The final verdict for Method C is `sig_quadrant OR sig_cross_S OR sig_cross_C`. This is a lenient criterion, as significance in *any* of these three sub-tests leads to rejecting the null hypothesis. This union of tests likely inflates the Type I error rate.

---

## Method D: Separate Univariate T-tests

- **Theory:** This method abandons the joint nature of the hypothesis `H₀: β = 0 AND γ = 0` and instead performs two separate one-sample t-tests: one for `β` and one for `γ`.
- The null hypothesis is rejected if *either* `β` is significantly different from zero *or* `γ` is significantly different from zero (`p_beta < 0.05 OR p_gamma < 0.05`).
- **Issues:** This approach suffers from the problem of multiple comparisons. By conducting two tests, the probability of a Type I error (falsely rejecting the null) is higher than the nominal `α` level of 0.05. It should be `1 - (1 - α)² ≈ 0.0975`.
- **Implementation (`n1_cosinor_strategy1.R`):**
    ```R
    method_D_S_pval <- t.test(beta_draws)$p.value
    method_D_C_pval <- t.test(gamma_draws)$p.value
    ...
    # Significance is determined if either p-value is below the threshold
    ```

---

## Method E: Wald Test

- **Theory:** The Wald test is a standard econometric and statistical test to assess the joint significance of multiple regression coefficients. In this context, it tests `H₀: β = 0 AND γ = 0`.
- The Wald statistic for this joint hypothesis is mathematically identical to the squared Mahalanobis distance between the coefficient estimates and the null hypothesis values (0,0).
- **Conclusion:** Method E (Wald test) and Method B (Mahalanobis distance) are theoretically and practically the same test. They produce identical statistics and p-values.

---

## Method F: Rayleigh "Z-score"

- **Theory:** This method incorrectly assumes a normal distribution for the amplitude `A`. The amplitude `A = sqrt(β² + γ²)`, where `β` and `γ` are normally distributed, actually follows a **Rayleigh distribution** (or a Rice distribution if the true amplitude is non-zero). A Rayleigh distribution is defined only for positive values and is skewed, not symmetric and zero-centered like a normal distribution. Therefore, applying a Z-test is inappropriate.
- **Implementation (`n1_cosinor_strategy1.R`):** The code contains two seemingly different calculations, but they are equivalent ways of performing the same (flawed) test against a Rayleigh distribution.
    1.  **P-value calculation:** `p = exp(-x² / (2σ²))`. This is the survival function (1 - CDF) of a Rayleigh distribution, giving the probability of observing a value greater than `x`.
        ```R
        # A_mean is the mean of the posterior amplitude draws
        # sigma_CS_pooled is the pooled standard deviation of β and γ draws
        method_F_rayleigh_pval <- exp(-A_mean^2 / (2 * sigma_CS_pooled^2))
        ```
    2.  **Critical value comparison:** The code also calculates a critical amplitude `a_crit_95` and checks if `A_mean` exceeds it.
        ```R
        a_crit_95 <- sigma_CS_pooled * sqrt(-2 * log(0.05))
        ```
        This critical value is derived by inverting the p-value formula. If we set the p-value `p` to our alpha level (0.05), we can solve for `x`, which becomes our critical value `a_crit_95`:
        `0.05 = exp(-a_crit_95² / (2σ²))`
        `log(0.05) = -a_crit_95² / (2σ²)`
        `-2 * log(0.05) = a_crit_95² / σ²`
        `σ² * -2 * log(0.05) = a_crit_95²`
        `a_crit_95 = σ * sqrt(-2 * log(0.05))`
        So, checking if `A_mean > a_crit_95` is identical to checking if `p < 0.05`.

---

## Heatmap Agreement Calculation

- **Question:** How is the agreement in the off-diagonal cells of `method_agreement_heatmap.png` calculated?
- **Analysis (`scripts/compare_multilevel_n1.R`):** The `agreement()` function in the script calculates the proportion of individuals for whom two methods yield the same outcome (i.e., both 'significant' or both 'non-significant').
    ```R
    # From scripts/compare_multilevel_n1.R
    agreement <- function(df, var1, var2) {
      mean(df[[var1]] == df[[var2]], na.rm = TRUE)
    }
    ```
- The diagonal of the heatmap shows 100% agreement because it compares a method to itself.
- The off-diagonal cells show the percentage of agreement between two different methods. For example, the cell at the intersection of "Method B" row and "Method D" column would show the proportion of all individuals where the significance decision of Method B matched that of Method D. It is **not** conditioned on significance status; it is the overall raw agreement across all cases.
- **Why does the Quadrant method (C) have low overlap?** Method C uses a liberal `OR` condition across three different heuristics. This makes it more likely to declare a result significant compared to other methods, thus reducing its agreement with them. It catches different "types" of significance that the other, more standard tests might miss, but likely at the cost of a higher false-positive rate.
