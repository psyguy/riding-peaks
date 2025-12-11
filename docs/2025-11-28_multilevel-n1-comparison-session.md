# Session Notes: Multilevel vs N=1 Cosinor Model Comparison

**Date:** 2025-11-28
**Project:** Riding Peaks - Circadian Rhythm Analysis
**Focus:** Comparing amplitude estimates and significance testing methods between multilevel and N=1 cosinor models

---

## Background

This session continued work on comparing two approaches to analyzing individual circadian rhythms:

1. **Multilevel Model (Partial Pooling)**: Uses hierarchical Bayesian modeling where individual parameters are drawn from population distributions, allowing information sharing between individuals.

2. **N=1 Model (No Pooling)**: Fits separate models for each individual with no information sharing, treating each person's data completely independently.

The goal was to understand how these modeling approaches differ in their amplitude estimates and significance testing outcomes.

---

## Session Overview

### Files Modified
- `scripts/compare_multilevel_n1.R` - Main comparison script created and enhanced

### Related Files (Referenced)
- `scripts/n1_cosinor_strategy1.R` - N=1 model fitting (separate models per person)
- `scripts/n1_cosinor_strategy2.R` - N=1 model fitting (single model with factor structure)
- `fits/ests_level1_pa_n1_strategy1.rds` - N=1 estimates with alternative tests
- `fits/ests_level1_pa_random_var.rds` - Multilevel estimates
- `fits/draws_pa_random_var.rds` - Multilevel posterior draws
- `fits/draws_pa_n1_strategy1.rds` - N=1 posterior draws

---

## Key Decisions and Actions

### 1. Visualization Changes

#### Decision: Use Line Segments Instead of Rectangles for CIs
**Rationale:** Line segments provide clearer visualization of confidence interval crossings between models while maintaining readability.

**Implementation:**
```r
# Horizontal CI segments (Multilevel CI on x-axis)
geom_segment(
  aes(x = amp_ci95_lower_ml, xend = amp_ci95_upper_ml,
      y = amp_median_n1, yend = amp_median_n1),
  color = "steelblue", alpha = 0.4, linewidth = 0.5
)
# Vertical CI segments (N=1 CI on y-axis)
geom_segment(
  aes(x = amp_median_ml, xend = amp_median_ml,
      y = amp_ci95_lower_n1, yend = amp_ci95_upper_n1),
  color = "steelblue", alpha = 0.4, linewidth = 0.5
)
```

**Location:** `compare_multilevel_n1.R:217-227`

#### Decision: Remove Shrinkage Plot from Combined Output
**Rationale:** Focus on the most informative comparisons (amplitude, CI width, HDR, R²) rather than including all possible plots.

**Implementation:**
```r
# Changed from: (p1 | p2) / (p3 | p4) / (p5 | plot_spacer())
combined <- (p1 | p2) / (p3 | p5)
```

**Impact:** Reduced combined plot height from 14 to 10 inches for better proportions.

**Location:** `compare_multilevel_n1.R:317`

---

### 2. Significance Testing Methods

#### Context: Six Methods Compared
The analysis compares six methods for testing H₀: Amplitude = 0 at the individual level:

1. **Method A: HDR (Highest Density Region)**
   - Tests whether origin [0,0] falls outside 95% HDR in (C,S) space
   - Criterion: `amp_hdr > 95`

2. **Method B: Bivariate Normal**
   - Fits bivariate normal to (C,S) draws, computes Mahalanobis distance
   - Criterion: `method_B_pval < 0.05`

3. **Method C: Quadrant Method**
   - Probability that draws fall in same quadrant as median
   - Criterion: `method_C_same_quadrant_prob > 0.95`

4. **Method D: Univariate CI (C or S)**
   - Tests whether C or S 95% CI excludes zero
   - Criterion: `method_D_C_excl_zero_95 | method_D_S_excl_zero_95`

5. **Method E: Wald Test**
   - Joint test H₀: C = 0 AND S = 0
   - Criterion: `method_E_pval < 0.05`

6. **Method F: Amplitude Z-score**
   - Z = mean(A) / SD(A)
   - Criterion: `method_F_zscore > 1.96`

#### Decision: Add Method D to All Comparisons
**Issue Identified:** Method D (univariate CI test) was computed in the individual model scripts but not included in between-model or within-model comparisons.

**Actions Taken:**

1. **Between-Model Agreement Table** (`compare_multilevel_n1.R:354-367`)
   - Added `univar_n1` and `univar_ml` to `sig_compare` data frame
   - Added "univar" and "Univar CI (C|S)" to method names

2. **Within-Model Comparison Matrices** (`compare_multilevel_n1.R:472, 483`)
   - Added `UnivarCI` to both `sig_n1` and `sig_ml` data frames
   - Updated `method_names_short` vector

**Verification:** Method D was already being computed in `compute_alternative_tests()` function in both N=1 and multilevel model scripts, so no changes to those files were needed.

---

### 3. Enhanced Agreement Tables

#### Decision: Provide Both Count and Percentage Tables
**Rationale:** Different audiences prefer different presentations; researchers may want raw counts for reporting, percentages for interpretation.

**Implementation:**

**Counts Table:**
```r
agreement_table_counts <- tibble(
  Method = method_names,
  `N=1 sig` = ...,      # Count finding significance in N=1
  `ML sig` = ...,       # Count finding significance in Multilevel
  `Both sig` = ...,     # Count finding significance in both
  `Neither sig` = ...,  # Count finding significance in neither
  `N=1 only` = ...,     # Count finding significance only in N=1
  `ML only` = ...,      # Count finding significance only in Multilevel
  Agreement = ...       # Total count of agreements
)
```

**Percentages Table:**
```r
agreement_table_pct <- tibble(
  Method = method_names,
  `N=1 sig %` = ...,
  `ML sig %` = ...,
  # ... same structure as counts but in percentages
)
```

**Disagreement Summary:**
```r
disagree_summary <- tibble(
  Method = method_names,
  `Disagree (n)` = `N=1 only` + `ML only`,
  `N=1 finds sig` = "count (percent)",
  `ML finds sig` = "count (percent)"
)
```

**Purpose:** The disagreement summary helps identify which model is more conservative - when the two models disagree, which one tends to find significance?

**Location:** `compare_multilevel_n1.R:372-436`

---

### 4. Within-Model Method Comparisons

#### Decision: Add Agreement Matrices for Each Model
**Rationale:** Understanding how different significance methods agree *within* each model helps evaluate method consistency and identify which methods are more/less conservative.

**Implementation:**

**Agreement Matrix Function:**
```r
compute_agreement_matrix <- function(sig_data, method_cols) {
  # Returns matrix where:
  # - Diagonal: count of individuals deemed "significant"
  # - Off-diagonal [i,j]: count where methods i and j agree
}
```

**Generated Outputs:**
- N=1 Model Agreement Matrix (counts and percentages)
- Multilevel Model Agreement Matrix (counts and percentages)

**Interpretation:**
- High diagonal values = method finds many significant effects
- High off-diagonal values = methods agree frequently
- Low off-diagonal values = methods disagree (one more conservative)

**Location:** `compare_multilevel_n1.R:448-509`

---

### 5. Pairwise Distribution Plots

#### Decision: Create Within-Model Pairplots Using GGally
**Rationale:** Visualize relationships between key amplitude metrics (median, CI width, HDR %, R²) to understand how uncertainty relates to effect size and variance explained.

**Implementation:**
```r
library(GGally)

pairplot_data_n1 <- compare_wide %>%
  select(
    `Amp median` = amp_median_n1,
    `CI width` = amp_ci_width_n1,
    `HDR %` = amp_hdr_n1,
    `R2` = R2_median_n1
  )

pp_n1 <- ggpairs(
  pairplot_data_n1,
  title = "N=1 Model: Pairwise Distributions",
  lower = list(continuous = wrap("points", alpha = 0.5, size = 1)),
  diag = list(continuous = wrap("densityDiag", alpha = 0.5)),
  upper = list(continuous = wrap("cor", size = 3))
)
```

**Outputs Generated:**
- `figs/pairplot_n1_model.pdf` - N=1 model pairplot
- `figs/pairplot_ml_model.pdf` - Multilevel model pairplot
- `figs/pairplot_both_models.pdf` - Side-by-side comparison

**Expected Insights:**
- Relationship between amplitude and CI width (larger effects = wider CIs?)
- Relationship between HDR % and R² (do they measure similar concepts?)
- Distribution shapes for each metric

**Location:** `compare_multilevel_n1.R:518-576`

---

### 6. Method Agreement Heatmaps

#### Decision: Create Visual Heatmap Representation
**Rationale:** Matrices are informative but visual heatmaps make patterns more immediately apparent.

**Implementation:**
```r
p_heatmap <- ggplot(heatmap_data, aes(x = Method1, y = Method2, fill = Agreement_pct)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Agreement_pct, 0)), size = 3) +
  scale_fill_gradient2(
    low = "white", mid = "steelblue", high = "darkblue",
    midpoint = 50, limits = c(0, 100)
  ) +
  facet_wrap(~ Model)
```

**Output:** `figs/method_agreement_heatmap.pdf`

**Interpretation Guide:**
- Diagonal: % of individuals deemed significant by that method
- Off-diagonal: % agreement between two methods
- Color intensity: darker = higher agreement/significance rate

**Location:** `compare_multilevel_n1.R:585-630`

---

## Technical Details

### Cosinor Model Specification
Both models fit the same functional form at the individual level:

```
y_ij = M_i + C_i * cos(2πt/24) + S_i * sin(2πt/24) + ε_ij
```

Where:
- `M_i` = individual mesor (mean level)
- `C_i`, `S_i` = cosine and sine coefficients
- `A_i = √(C_i² + S_i²)` = amplitude
- `φ_i = atan2(S_i, C_i)` = acrophase (peak time)

**Multilevel Model:** Individual parameters drawn from population distributions
```
[M_i, C_i, S_i]' ~ MVN(μ_pop, Σ_pop)
```

**N=1 Model:** Individual parameters estimated independently (no population distribution)

### R² Calculation
Proportion of variance explained by circadian rhythm:

```r
R2 = (A²/2) / (A²/2 + σ²)
```

Where:
- `A²/2` = variance explained by sinusoidal component
- `σ²` = residual variance

---

## Outputs Generated

### Figures Created
1. `figs/multilevel_vs_n1_comparison.pdf` (12×10 inches)
   - Combined 4-panel plot: amplitude, CI width, HDR, R²

2. `figs/multilevel_vs_n1_comparison.png` (12×10 inches, 150 dpi)
   - PNG version of combined plot

3. `figs/amp_median_ml_vs_n1.pdf` (7×7 inches)
   - Standalone amplitude comparison plot

4. `figs/pairplot_n1_model.pdf` (8×8 inches)
   - N=1 model pairwise distributions

5. `figs/pairplot_ml_model.pdf` (8×8 inches)
   - Multilevel model pairwise distributions

6. `figs/pairplot_both_models.pdf` (16×8 inches)
   - Side-by-side pairplots

7. `figs/method_agreement_heatmap.pdf` (10×5 inches)
   - Within-model method agreement heatmap

### Data Saved
- `fits/comparison_multilevel_n1_pa.rds`
  - Contains `compare_wide` data frame with all estimates from both models

---

## Key Insights to Look For

When running the comparison script, examine:

### 1. Shrinkage Effects
- Do multilevel estimates tend to be smaller than N=1 estimates?
- Are individuals with small N=1 amplitudes shrunk upward?
- Are individuals with large N=1 amplitudes shrunk downward?

### 2. Uncertainty Quantification
- Are N=1 confidence intervals wider than multilevel CIs?
- Does multilevel modeling "borrow strength" to reduce uncertainty?

### 3. Significance Agreement
- Do the models agree on which individuals have significant rhythms?
- When they disagree, which model is more conservative?
- Is disagreement related to sample size, effect size, or data quality?

### 4. Method Consistency
- Do all six methods agree on significance decisions?
- Which methods are most conservative? Most liberal?
- Does method agreement differ between models?

### 5. Practical Significance (R²)
- What proportion of variance is explained by circadian rhythms?
- Is statistical significance (p < 0.05) associated with practical significance (R² > 0.10)?
- Do multilevel and N=1 models give similar R² estimates?

---

## Next Steps / Future Directions

### Potential Analyses
1. **Examine specific disagreement cases**
   - Which individuals show largest differences between models?
   - What data characteristics lead to disagreement?

2. **Simulate ground truth scenarios**
   - Generate data with known amplitudes
   - Compare which model recovers truth better under various conditions

3. **Level-2 (population-level) analysis**
   - What does multilevel model tell us about population means?
   - Are there subgroups with different circadian patterns?

4. **Sensitivity analyses**
   - How do results change with different priors?
   - Impact of different significance thresholds (90%, 99%)?

5. **Method validation**
   - Compare with published benchmark datasets
   - Cross-validation to assess predictive performance

### Documentation Improvements
- Add interpretation guide for each significance method
- Document expected type I error rates under different scenarios
- Create decision tree for which method to use when

---

## References

### Key Concepts
- **Partial pooling:** Compromise between complete pooling (all individuals identical) and no pooling (all individuals independent)
- **HDR:** Highest Density Region - smallest region containing X% of posterior probability
- **Mahalanobis distance:** Accounts for covariance between C and S
- **Acrophase:** Time of peak in circadian rhythm (in radians or hours)

### Related Documentation
- `scripts/n1_cosinor_strategy1.R` - Detailed comments on N=1 fitting approach
- `scripts/Tihomir review.pdf` - Reviewer suggestions that motivated methods A-F

---

## Reproducibility Notes

### R Package Versions Required
- `brms` - Bayesian regression models
- `tidyverse` - Data manipulation and visualization
- `patchwork` - Combining plots
- `GGally` - Pairwise plot matrices
- `circular` - Circular statistics for phase
- `ks` - Kernel density estimation for HDR
- `here` - Path management

### Computational Requirements
- N=1 strategy 1: ~1-2 minutes per individual (sequential fitting)
- N=1 strategy 2: ~30-60 minutes (single large model)
- Multilevel model: ~30-60 minutes
- Comparison script: ~5-10 minutes

### File Dependencies
The comparison script requires these files to exist:
```
fits/ests_level1_pa_n1_strategy1.rds
fits/ests_level1_pa_random_var.rds
fits/draws_pa_random_var.rds
fits/draws_pa_n1_strategy1.rds
```

---

## Session Git History

**Commit:** `ebc15aa`
**Message:** "Add multilevel vs N=1 comparison script with enhanced visualizations"
**Files Changed:** `scripts/compare_multilevel_n1.R` (646 insertions, new file)

**Key Features Committed:**
- Amplitude comparison plot with CI line segments
- Between-model agreement tables (counts and percentages)
- All 6 significance testing methods including Method D
- Within-model method comparison matrices
- Within-model pairwise distribution plots
- Method agreement heatmaps

---

## Session Timeline

1. **Modification Request:** Change CI visualization from rectangles to line segments
2. **Modification Request:** Remove shrinkage plot from combined output
3. **Enhancement Request:** Add within-model method comparison matrices
4. **Enhancement Request:** Add within-model pairplots using GGally
5. **Bug Fix:** Add Method D (univariate CI) to all comparisons
6. **Enhancement Request:** Add percentage tables alongside count tables
7. **Enhancement Request:** Add disagreement summary showing which model finds significance
8. **Commit:** Save all changes to git repository
9. **Documentation:** Create this comprehensive session notes document

---

**End of Session Notes**
