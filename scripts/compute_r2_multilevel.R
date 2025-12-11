# =============================================================================
# R² Calculations for Multilevel Cosinor Models
# =============================================================================
#
# This script implements the R² calculations described in:
#   docs/r2_theory_and_implementation.md
#
# Three methods are computed:
#   1. Distribution of Individual R² (Method 1)
#   2. Simplified Marginal & Conditional R² (Method 2) - for comparison only
#   3. Exact R² for Gapped Data (Method 3) - recommended for ESM data
#
# The script uses saved RDS objects from the multilevel model fitting.
#
# =============================================================================

library(tidyverse)
library(brms)
library(here)

# -----------------------------------------------------------------------------
# 1. Load saved model outputs
# -----------------------------------------------------------------------------

cat("Loading saved model outputs...\n")

# Choose item (pa or fitbit)
item_ <- "pa"

# Load posterior draws (contains individual-level parameters per MCMC draw)
d_draws <- readRDS(here("fits", paste0("draws_", item_, "_random_var.rds")))

# Load level-1 estimates (contains summaries like amp_median, R2_median)
ests_level1 <- readRDS(here("fits", paste0("ests_level1_", item_, "_random_var.rds")))

# Load the original data to compute empirical predictor properties
d <- readRDS(here("data", "d_leuven_joined.rds")) %>%
  filter(item == item_)

cat("Loaded data for item:", item_, "\n")
cat("  - d_draws:", nrow(d_draws), "rows (draws x individuals)\n")
cat("  - ests_level1:", nrow(ests_level1), "individuals\n")
cat("  - d:", nrow(d), "observations\n\n")

# -----------------------------------------------------------------------------
# 2. Compute empirical predictor properties (for Method 3)
# -----------------------------------------------------------------------------
# These account for the gapped ESM structure (no night observations)

cat("Computing empirical predictor properties...\n")

# Empirical expectations
mu_c <- mean(d$co)  # E[cos(2πt/24)]
mu_s <- mean(d$si)  # E[sin(2πt/24)]

# Empirical variances
v_c <- var(d$co)    # Var[cos(2πt/24)]
v_s <- var(d$si)    # Var[sin(2πt/24)]

# Covariance
c_cs <- cov(d$co, d$si)  # Cov[cos, sin]

# Second moments (needed for random effects variance)
E_cos2 <- mean(d$co^2)      # E[cos²] = v_c + mu_c²
E_sin2 <- mean(d$si^2)      # E[sin²] = v_s + mu_s²
E_cossin <- mean(d$co * d$si)  # E[cos·sin] = c_cs + mu_c·mu_s

cat("\nEmpirical predictor properties:\n")
cat("  E[cos] (mu_c):", round(mu_c, 4), "(theoretical for full cycle: 0)\n")
cat("  E[sin] (mu_s):", round(mu_s, 4), "(theoretical for full cycle: 0)\n")
cat("  Var[cos] (v_c):", round(v_c, 4), "(theoretical for full cycle: 0.5)\n")
cat("  Var[sin] (v_s):", round(v_s, 4), "(theoretical for full cycle: 0.5)\n")
cat("  Cov[cos,sin] (c_cs):", round(c_cs, 4), "(theoretical for full cycle: 0)\n")
cat("  E[cos²]:", round(E_cos2, 4), "(theoretical for full cycle: 0.5)\n")
cat("  E[sin²]:", round(E_sin2, 4), "(theoretical for full cycle: 0.5)\n")
cat("  E[cos·sin]:", round(E_cossin, 4), "(theoretical for full cycle: 0)\n\n")

# -----------------------------------------------------------------------------
# 3. Extract model parameters from posterior draws
# -----------------------------------------------------------------------------

cat("Extracting model parameters from posterior draws...\n")

# Number of individuals
n_ids <- n_distinct(d_draws$id)

# Get fixed effects (population-level) per iteration
# These are stored as b_co and b_si in d_draws
fixed_effects <- d_draws %>%
  select(iteration, b_co, b_si, b_logsigma) %>%
  distinct()

# Population-level amplitude
fixed_effects <- fixed_effects %>%
  mutate(
    gamma_c = b_co,
    gamma_s = b_si,
    A_pop = sqrt(gamma_c^2 + gamma_s^2),
    rho = b_logsigma  # log-scale intercept for sigma
  )

# Get random effect variances per iteration
# Random effects are: co - b_co, si - b_si, logsigma - b_logsigma
random_effect_vars <- d_draws %>%
  group_by(iteration) %>%
  summarise(
    # Variance of random intercept (mesor)
    var_uM = var(mesor - first(b_Intercept)),
    # Variance of random slopes
    var_uC = var(co - first(b_co)),
    var_uS = var(si - first(b_si)),
    # Covariance of random slopes
    cov_uCS = cov(co - first(b_co), si - first(b_si)),
    # Variance of random sigma (on log scale)
    var_logsigma = var(logsigma - first(b_logsigma)),
    # tau = SD of random sigma deviations
    tau = sd(logsigma - first(b_logsigma)),
    .groups = "drop"
  )

# Merge fixed and random effect summaries
model_params <- fixed_effects %>%
  left_join(random_effect_vars, by = "iteration")

cat("  Number of iterations:", nrow(model_params), "\n")
cat("  Population amplitude (median):", round(median(model_params$A_pop), 4), "\n")
cat("  tau (SD of random log-sigma, median):", round(median(model_params$tau), 4), "\n\n")

# -----------------------------------------------------------------------------
# 4. Method 1: Distribution of Individual R²
# -----------------------------------------------------------------------------
# R²_i = (A_i² / 2) / (A_i² / 2 + σ_i²)
# This is computed for each individual using their specific amplitude and sigma

cat("Computing Method 1: Distribution of Individual R²...\n")

# Compute R² for each individual and iteration
individual_r2 <- d_draws %>%
  mutate(
    # Individual amplitude (already computed in d_draws as 'amp')
    A_i = amp,
    # Individual residual variance
    sigma_i = sigma,
    sigma_i_sq = sigma^2,
    # Individual R² (assuming symmetric time sampling within person)
    R2_i = (A_i^2 / 2) / (A_i^2 / 2 + sigma_i_sq)
  )

# Summarize across iterations for each individual
individual_r2_summary <- individual_r2 %>%
  group_by(id) %>%
  summarise(
    R2_mean = mean(R2_i),
    R2_median = median(R2_i),
    R2_ci_lower = quantile(R2_i, 0.025),
    R2_ci_upper = quantile(R2_i, 0.975),
    A_median = median(A_i),
    sigma_median = median(sigma_i),
    .groups = "drop"
  )

# Distribution across individuals (using posterior medians)
cat("\nMethod 1 Results - Distribution of Individual R²:\n")
cat("  Median R² across individuals:", round(median(individual_r2_summary$R2_median), 4), "\n")
cat("  Mean R² across individuals:", round(mean(individual_r2_summary$R2_median), 4), "\n")
cat("  IQR of R²:", round(quantile(individual_r2_summary$R2_median, 0.25), 4), "-",
    round(quantile(individual_r2_summary$R2_median, 0.75), 4), "\n")
cat("  Range of R²:", round(min(individual_r2_summary$R2_median), 4), "-",
    round(max(individual_r2_summary$R2_median), 4), "\n")
cat("  Proportion with R² > 0.01:", round(mean(individual_r2_summary$R2_median > 0.01), 3), "\n")
cat("  Proportion with R² > 0.05:", round(mean(individual_r2_summary$R2_median > 0.05), 3), "\n")
cat("  Proportion with R² > 0.10:", round(mean(individual_r2_summary$R2_median > 0.10), 3), "\n\n")

# -----------------------------------------------------------------------------
# 5. Compute E[σ_i²] - The Correct Residual Variance
# -----------------------------------------------------------------------------
# For log-normal sigma: E[σ²] = exp(2ρ + 2τ²), NOT exp(2ρ)
# This is critical for accurate R² calculations

cat("Computing E[σ_i²] with log-normal correction...\n")

# Method A: From model parameters (theoretical)
model_params <- model_params %>%
  mutate(
    # Naive (incorrect) approach
    E_sigma_sq_naive = exp(2 * rho),
    # Correct approach accounting for random sigma
    E_sigma_sq_correct = exp(2 * rho + 2 * tau^2),
    # Correction factor
    correction_factor = exp(2 * tau^2)
  )

# Method B: Empirical (from actual posterior draws)
E_sigma_sq_empirical <- d_draws %>%
  group_by(iteration) %>%
  summarise(E_sigma_sq = mean(sigma^2), .groups = "drop")

model_params <- model_params %>%
  left_join(E_sigma_sq_empirical, by = "iteration")

cat("\nResidual Variance (E[σ_i²]):\n")
cat("  Naive (exp(2ρ)), median:", round(median(model_params$E_sigma_sq_naive), 4), "\n")
cat("  Correct (exp(2ρ + 2τ²)), median:", round(median(model_params$E_sigma_sq_correct), 4), "\n")
cat("  Empirical (mean(σ_i²)), median:", round(median(model_params$E_sigma_sq), 4), "\n")
cat("  Correction factor (exp(2τ²)), median:", round(median(model_params$correction_factor), 4), "\n")
cat("  Underestimate if using naive:",
    round((1 - median(model_params$E_sigma_sq_naive) / median(model_params$E_sigma_sq_correct)) * 100, 1), "%\n\n")

# -----------------------------------------------------------------------------
# 6. Method 2: Simplified Marginal & Conditional R² (for comparison)
# -----------------------------------------------------------------------------
# CAUTION: This assumes E[cos] = E[sin] = 0, which is WRONG for gapped ESM data
# Included for comparison only

cat("Computing Method 2: Simplified R² (for comparison - NOT recommended for gapped data)...\n")

method2_r2 <- model_params %>%
  mutate(
    # Simplified variance components (assuming E[cos]=E[sin]=0)
    Var_Fixed_simple = (gamma_c^2 + gamma_s^2) / 2,  # = A_pop² / 2
    Var_Random_simple = var_uM + (var_uC + var_uS) / 2,
    Var_Total_simple = Var_Fixed_simple + Var_Random_simple + E_sigma_sq_correct,

    # Simplified R² values
    R2_marginal_simple = Var_Fixed_simple / Var_Total_simple,
    R2_conditional_simple = (Var_Fixed_simple + Var_Random_simple) / Var_Total_simple
  )

cat("\nMethod 2 Results - Simplified R² (assumes E[cos]=E[sin]=0):\n")
cat("  Marginal R² (median):", round(median(method2_r2$R2_marginal_simple), 4), "\n")
cat("  Marginal R² (95% CI):", round(quantile(method2_r2$R2_marginal_simple, 0.025), 4), "-",
    round(quantile(method2_r2$R2_marginal_simple, 0.975), 4), "\n")
cat("  Conditional R² (median):", round(median(method2_r2$R2_conditional_simple), 4), "\n")
cat("  Conditional R² (95% CI):", round(quantile(method2_r2$R2_conditional_simple, 0.025), 4), "-",
    round(quantile(method2_r2$R2_conditional_simple, 0.975), 4), "\n\n")

# -----------------------------------------------------------------------------
# 7. Method 3: Exact R² for Gapped Data (RECOMMENDED)
# -----------------------------------------------------------------------------
# Uses empirical predictor properties to account for night gap

cat("Computing Method 3: Exact R² for Gapped Data (RECOMMENDED)...\n")

method3_r2 <- model_params %>%
  mutate(
    # Exact variance of fixed effects
    # Var(Fixed) = γ_c² v_c + γ_s² v_s + 2 γ_c γ_s c_cs
    Var_Fixed_exact = gamma_c^2 * v_c + gamma_s^2 * v_s + 2 * gamma_c * gamma_s * c_cs,

    # Exact variance of random effects
    # Var(Random) = Var(u_M) + Var(u_C) E[cos²] + Var(u_S) E[sin²] + 2 Cov(u_C, u_S) E[cos·sin]
    Var_Random_exact = var_uM +
                       var_uC * E_cos2 +
                       var_uS * E_sin2 +
                       2 * cov_uCS * E_cossin,

    # Total variance
    Var_Total_exact = Var_Fixed_exact + Var_Random_exact + E_sigma_sq_correct,

    # Exact R² values
    R2_marginal_exact = Var_Fixed_exact / Var_Total_exact,
    R2_conditional_exact = (Var_Fixed_exact + Var_Random_exact) / Var_Total_exact
  )

cat("\nMethod 3 Results - Exact R² for Gapped Data:\n")
cat("  Marginal R² (median):", round(median(method3_r2$R2_marginal_exact), 4), "\n")
cat("  Marginal R² (95% CI):", round(quantile(method3_r2$R2_marginal_exact, 0.025), 4), "-",
    round(quantile(method3_r2$R2_marginal_exact, 0.975), 4), "\n")
cat("  Conditional R² (median):", round(median(method3_r2$R2_conditional_exact), 4), "\n")
cat("  Conditional R² (95% CI):", round(quantile(method3_r2$R2_conditional_exact, 0.025), 4), "-",
    round(quantile(method3_r2$R2_conditional_exact, 0.975), 4), "\n\n")

# -----------------------------------------------------------------------------
# 8. Compare Methods 2 and 3 (Bias from ignoring gapped structure)
# -----------------------------------------------------------------------------

cat("Comparison: Method 2 vs Method 3 (bias from ignoring gapped structure)...\n")

comparison <- tibble(
  Measure = c("Var(Fixed)", "Var(Random)", "R² Marginal", "R² Conditional"),
  Method2_median = c(
    median(method2_r2$Var_Fixed_simple),
    median(method2_r2$Var_Random_simple),
    median(method2_r2$R2_marginal_simple),
    median(method2_r2$R2_conditional_simple)
  ),
  Method3_median = c(
    median(method3_r2$Var_Fixed_exact),
    median(method3_r2$Var_Random_exact),
    median(method3_r2$R2_marginal_exact),
    median(method3_r2$R2_conditional_exact)
  )
) %>%
  mutate(
    Difference = Method2_median - Method3_median,
    Pct_Difference = round((Method2_median - Method3_median) / Method3_median * 100, 1)
  )

print(comparison)
cat("\n")

# -----------------------------------------------------------------------------
# 9. Compile final R² summary
# -----------------------------------------------------------------------------

cat("=== FINAL R² SUMMARY ===\n\n")

# Create summary table
r2_summary <- tibble(
  Method = c(
    "Method 1: Individual R² (median across persons)",
    "Method 1: Individual R² (mean across persons)",
    "Method 2: Marginal R² (simplified)",
    "Method 2: Conditional R² (simplified)",
    "Method 3: Marginal R² (exact, RECOMMENDED)",
    "Method 3: Conditional R² (exact, RECOMMENDED)"
  ),
  Estimate = c(
    median(individual_r2_summary$R2_median),
    mean(individual_r2_summary$R2_median),
    median(method2_r2$R2_marginal_simple),
    median(method2_r2$R2_conditional_simple),
    median(method3_r2$R2_marginal_exact),
    median(method3_r2$R2_conditional_exact)
  ),
  CI_lower = c(
    quantile(individual_r2_summary$R2_median, 0.025),
    NA,  # No CI for mean of medians
    quantile(method2_r2$R2_marginal_simple, 0.025),
    quantile(method2_r2$R2_conditional_simple, 0.025),
    quantile(method3_r2$R2_marginal_exact, 0.025),
    quantile(method3_r2$R2_conditional_exact, 0.025)
  ),
  CI_upper = c(
    quantile(individual_r2_summary$R2_median, 0.975),
    NA,
    quantile(method2_r2$R2_marginal_simple, 0.975),
    quantile(method2_r2$R2_conditional_simple, 0.975),
    quantile(method3_r2$R2_marginal_exact, 0.975),
    quantile(method3_r2$R2_conditional_exact, 0.975)
  )
) %>%
  mutate(
    Estimate = round(Estimate, 4),
    CI_lower = round(CI_lower, 4),
    CI_upper = round(CI_upper, 4),
    CI = paste0("[", CI_lower, ", ", CI_upper, "]")
  )

print(r2_summary %>% select(Method, Estimate, CI))

# -----------------------------------------------------------------------------
# 10. Save results
# -----------------------------------------------------------------------------

cat("\nSaving results...\n")

results <- list(
  # Empirical predictor properties
  predictor_properties = tibble(
    mu_c = mu_c, mu_s = mu_s,
    v_c = v_c, v_s = v_s,
    c_cs = c_cs,
    E_cos2 = E_cos2, E_sin2 = E_sin2, E_cossin = E_cossin
  ),

  # Individual R² (Method 1)
  individual_r2 = individual_r2_summary,

  # Population R² per iteration (Methods 2 & 3)
  population_r2_draws = method3_r2 %>%
    left_join(
      method2_r2 %>% select(iteration, R2_marginal_simple, R2_conditional_simple),
      by = "iteration"
    ),

  # Summary table
  summary = r2_summary,

  # Comparison of methods
  method_comparison = comparison
)

saveRDS(results, here("fits", paste0("r2_results_", item_, ".rds")))

cat("Results saved to: fits/r2_results_", item_, ".rds\n")

# -----------------------------------------------------------------------------
# 11. Visualization
# -----------------------------------------------------------------------------

cat("\nCreating visualizations...\n")

library(patchwork)

# Plot 1: Distribution of Individual R²
p1 <- ggplot(individual_r2_summary, aes(x = R2_median)) +
  geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7, color = "white") +
  geom_vline(xintercept = median(individual_r2_summary$R2_median),
             color = "red", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = c(0.01, 0.10),
             color = "gray40", linetype = "dotted", linewidth = 0.8) +
  annotate("text", x = 0.01, y = Inf, label = "1%", vjust = 2, hjust = -0.2, size = 3) +
  annotate("text", x = 0.10, y = Inf, label = "10%", vjust = 2, hjust = -0.2, size = 3) +
  labs(
    title = "Method 1: Distribution of Individual R²",
    subtitle = paste0("Median = ", round(median(individual_r2_summary$R2_median) * 100, 1), "%"),
    x = expression(R[i]^2),
    y = "Count"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

# Plot 2: Individual R² vs Amplitude
p2 <- ggplot(individual_r2_summary, aes(x = A_median, y = R2_median)) +
  geom_point(alpha = 0.5, color = "steelblue") +
  geom_smooth(method = "loess", color = "red", se = FALSE) +
  labs(
    title = "Individual R² vs Amplitude",
    x = "Amplitude (median)",
    y = expression(R[i]^2)
  ) +
  theme_minimal()

# Plot 3: Individual R² vs Residual SD
p3 <- ggplot(individual_r2_summary, aes(x = sigma_median, y = R2_median)) +
  geom_point(alpha = 0.5, color = "steelblue") +
  geom_smooth(method = "loess", color = "red", se = FALSE) +
  labs(
    title = "Individual R² vs Residual SD",
    x = expression(sigma[i] ~ "(median)"),
    y = expression(R[i]^2)
  ) +
  theme_minimal()

# Plot 4: Posterior distribution of population R²
r2_long <- method3_r2 %>%
  select(iteration, R2_marginal_exact, R2_conditional_exact) %>%
  pivot_longer(
    cols = c(R2_marginal_exact, R2_conditional_exact),
    names_to = "Type",
    values_to = "R2"
  ) %>%
  mutate(Type = recode(Type,
                       "R2_marginal_exact" = "Marginal R²",
                       "R2_conditional_exact" = "Conditional R²"))

p4 <- ggplot(r2_long, aes(x = R2, fill = Type)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Marginal R²" = "coral", "Conditional R²" = "steelblue")) +
  labs(
    title = "Method 3: Population R² (Exact for Gapped Data)",
    x = expression(R^2),
    y = "Density",
    fill = NULL
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

# Combine plots
combined_plot <- (p1 | p4) / (p2 | p3) +
  plot_annotation(
    title = paste0("R² Analysis for ", toupper(item_), " Data"),
    theme = theme(plot.title = element_text(size = 14, face = "bold"))
  )

ggsave(
  here("figs", paste0("r2_analysis_", item_, ".pdf")),
  combined_plot,
  width = 12, height = 10
)

cat("Figure saved to: figs/r2_analysis_", item_, ".pdf\n")

# -----------------------------------------------------------------------------
# 12. Additional: Compare with naive sigma (without correction)
# -----------------------------------------------------------------------------

cat("\n=== Sensitivity Analysis: Effect of Log-Normal Correction ===\n")

naive_r2 <- model_params %>%
  mutate(
    Var_Fixed_exact = gamma_c^2 * v_c + gamma_s^2 * v_s + 2 * gamma_c * gamma_s * c_cs,
    Var_Random_exact = var_uM + var_uC * E_cos2 + var_uS * E_sin2 + 2 * cov_uCS * E_cossin,
    # Using NAIVE sigma (no correction for random sigma)
    Var_Total_naive = Var_Fixed_exact + Var_Random_exact + E_sigma_sq_naive,
    R2_marginal_naive = Var_Fixed_exact / Var_Total_naive,
    R2_conditional_naive = (Var_Fixed_exact + Var_Random_exact) / Var_Total_naive
  )

cat("\nIf we ignore the log-normal correction (use exp(2ρ) instead of exp(2ρ + 2τ²)):\n")
cat("  Marginal R² would be:", round(median(naive_r2$R2_marginal_naive), 4),
    "instead of", round(median(method3_r2$R2_marginal_exact), 4), "\n")
cat("  Conditional R² would be:", round(median(naive_r2$R2_conditional_naive), 4),
    "instead of", round(median(method3_r2$R2_conditional_exact), 4), "\n")
cat("  Overestimation in marginal R²:",
    round((median(naive_r2$R2_marginal_naive) - median(method3_r2$R2_marginal_exact)) /
            median(method3_r2$R2_marginal_exact) * 100, 1), "%\n")
cat("  Overestimation in conditional R²:",
    round((median(naive_r2$R2_conditional_naive) - median(method3_r2$R2_conditional_exact)) /
            median(method3_r2$R2_conditional_exact) * 100, 1), "%\n")

cat("\nDone!\n")
beepr::beep(5)
