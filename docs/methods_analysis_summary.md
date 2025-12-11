# Analysis of Significance Testing Methods for Cosinor Models

This document summarizes the implementation and interpretation of six alternative methods (A-F), suggested during peer review, for testing the significance of circadian cycles in time-series data. It is based on an analysis of the project's R scripts and an interactive Q&A session.

## 1. Detailed Description of Methods (A-F)

The goal of each method is to test the null hypothesis that the circadian amplitude is zero (H₀: A = 0), which is equivalent to the joint hypothesis that the cosinor coefficients C and S are both zero (H₀: C=0 and S=0).

### Method A: Highest Density Region (HDR)
*   **Concept:** This non-parametric method checks if the point `(0,0)` in the `(C, S)` plane lies within the 95% Highest Density Region of the joint posterior distribution. The HDR is the smallest possible region that contains 95% of the posterior probability mass.
*   **Criterion:** If the origin `(0,0)` falls *outside* this 95% region, the result is significant. This is the primary method used in the original paper.

### Method B: Bivariate Normal Density
*   **Concept:** This method approximates the posterior distribution of `(C, S)` with a bivariate normal distribution. It then assesses how likely the origin `(0,0)` is under this idealized distribution.
*   **Criterion:** Significance is determined by calculating the probability mass of the distribution that is "more central" than the origin. If the origin lies in the outer 5% tail of the distribution (i.e., `prob_higher_dens > 0.95`), the result is significant.

### Method C: Quadrant Method
*   **Concept:** A simple, non-parametric test that checks how concentrated the posterior distribution is in a single quadrant. It calculates the proportion of posterior draws that fall in the quadrant diagonally opposite to the median `(C, S)` estimate.
*   **Criterion:** If the probability of draws falling in this "opposite" quadrant is very low (e.g., `< 0.05`), it implies the posterior is well-defined away from the origin, and the result is significant.

### Method D: Univariate Credible Intervals (CI)
*   **Concept:** This method breaks the joint hypothesis into two separate tests on the marginal distributions of C and S.
*   **Criterion:** The result is considered significant if the 95% credible interval for `C` excludes zero **OR** if the 95% credible interval for `S` excludes zero.

### Method E: Joint Wald Test
*   **Concept:** A standard, parametric joint hypothesis test. It calculates a single statistic based on how far the mean estimate `(mean(C), mean(S))` is from the origin `(0,0)`, accounting for the covariance between C and S.
*   **Criterion:** The resulting Wald statistic follows a chi-squared (χ²) distribution, which is used to derive a standard p-value. If `p < 0.05`, the result is significant.

### Method F: Simple Amplitude Z-score
*   **Concept:** This was the reviewer's suggestion to treat the amplitude `A` as a simple statistic and calculate a Z-score.
*   **Criterion:** `Z = mean(A) / sd(A) > 1.96`.
*   **Flaw:** As noted in the R scripts and our discussion, this test is theoretically flawed because the amplitude `A` under the null hypothesis follows a **Rayleigh distribution**, not a Gaussian one, and its mean is not zero.

### Corrected Method F: Rayleigh Test (RayleighZ)
*   **Concept:** A corrected version of the amplitude test that properly uses the Rayleigh distribution as the null distribution for amplitude.
*   **Criterion:** It compares the observed `mean(A)` to the true 95% critical value derived from the corresponding Rayleigh distribution. A significant result occurs if the observed mean amplitude is greater than this critical value.

## 2. R Implementation Analysis

The core logic for all tests is implemented in the `compute_alternative_tests` function within `scripts/n1_cosinor_strategy1.R` and is then reused in `scripts/compare_multilevel_n1.R`.

*   **Method A (HDR):** Implemented via the `compute_alpha_hdr` function using the `ks` package. The test is `amp_hdr > 95`.
*   **Method B (BivNorm):** Implemented by calculating `mahal_origin_sq` (the squared Mahalanobis distance) and converting it to a cumulative probability, `method_B_prob_higher_dens`, which is then checked to be `> 0.95`.
*   **Method C (Quadrant):** Implemented by calculating `method_C_origin_quadrant_prob` and checking if it is `< 0.05`.
*   **Method D (UnivarCI):** Implemented with a logical OR: `method_D_C_excl_zero_95 | method_D_S_excl_zero_95`.
*   **Method E (Wald):** Implemented by calculating `method_E_wald_stat` (which is identical to `mahal_origin_sq`) and deriving `method_E_pval` from a chi-squared distribution.
*   **RayleighZ Test:** Implemented by first estimating the Rayleigh scale parameter `sigma_CS_pooled`, then calculating the critical value `a_crit_95`, and finally checking if `A_mean > a_crit_95`.

## 3. Key Questions and Clarifications

This section summarizes the key questions and doubts raised during our discussion.

#### Q: How are Method B (BivNorm), Method E (Wald), the Mahalanobis distance, and the 95% Confidence Ellipse related?
**A:** They are all different names for or components of the same statistical test.
1.  The **95% Confidence Ellipse** is the boundary containing 95% of the probability.
2.  The **Mahalanobis distance** is the metric used to draw this ellipse. All points on the ellipse have the same Mahalanobis distance from the center.
3.  The **Wald test statistic (Method E)** is mathematically identical to the squared Mahalanobis distance of the origin from the posterior mean.
4.  **Method B** uses the same Mahalanobis distance to calculate a cumulative probability, which is just a different way of expressing the same result as the Wald test's p-value.

#### Q: For Method B, is the Mahalanobis distance used, or the density calculation?
**A:** Both. The Mahalanobis distance is the essential first step used to calculate both the density and, more importantly, the final cumulative probability (`prob_higher_dens`) used for the test.

#### Q: Where does the formula `1 - exp(-0.5 * mahal_origin_sq)` come from?
**A:** This is the textbook formula for the **Cumulative Distribution Function (CDF) of a chi-squared distribution with 2 degrees of freedom**. The squared Mahalanobis distance of a 2D variable follows this distribution, so this formula correctly converts the distance into a cumulative probability (from 0 to 1).

#### Q: Are diagnostic values like `dens_ratio` and `cross_axis_prob` used in the tests?
**A:** No. These values are calculated and saved for potential exploratory analysis but are not used in any of the final significance tests summarized in the heatmaps.

#### Q: Are the different methods (A-F) combined using an OR or AND condition?
**A:** Neither. They are treated as separate, alternative tests. The goal of the analysis is to compare their results, not to combine them into a single verdict. The only exception is *within* Method D, which is an OR condition of two sub-tests.

#### Q: Is there anything wrong with OR-ing the tests in Method D?
**A:** It's not "wrong" in the context of this comparison, but it has the well-known statistical implication of **inflating the Type I error rate** (increasing false positives) due to multiple comparisons. This makes Method D a more "liberal" test, a fact that is confirmed in the heatmap.

#### Q: How is `mean(A)` related to `S` and `C`?
**A:** There is no simple formula. `mean(A)` is the *mean of the posterior amplitudes*, while `sqrt(mean(C)² + mean(S)²) ` is the *amplitude of the mean posterior*. These are not the same. A large difference between them (specifically when `mean(A)` is large but the other is near zero) indicates high uncertainty in the cycle's phase.

#### Q: How are the off-diagonal cells in the heatmap calculated?
**A:** They represent the percentage of agreement across **all** individuals. An agreement is counted when two methods both say "significant" OR when they both say "not significant".

## 4. Interpretation of the Heatmap

The `method_agreement_heatmap.png` provides several key insights:

1.  **Methods Form "Families":** The tests cluster into three groups:
    *   **Core Block (HDR, BivNorm, Wald, RayleighZ):** These methods are robust and highly consistent with each other (>90% agreement).
    *   **Liberal Outlier (UnivarCI):** Method D is the most lenient test, finding the most significant results.
    *   **Conservative Outlier (Quadrant):** Method C is the most stringent test, finding the fewest significant results. This is because it is sensitive to the overall shape of the posterior (e.g., "banana" shapes) and not just the density at the origin.

2.  **Multilevel Modeling is More Conservative:** The multilevel model yields significantly fewer significant results across all methods compared to the N=1 models. This demonstrates the effect of Bayesian shrinkage, where uncertain individual estimates are pulled toward the group average, demanding stronger evidence for significance.

## 5. Summary

The project's R scripts correctly implement the six suggested testing procedures. The analysis reveals important differences in their statistical properties, with a "core" group of four methods providing the most consistent results. The comparison between N=1 and multilevel models clearly demonstrates the conservative nature of shrinkage in hierarchical modeling.
