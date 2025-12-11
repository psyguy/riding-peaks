# =============================================================================
# Compare Multilevel vs N=1 Cosinor Estimates
# =============================================================================
#
# This script compares individual-level amplitude estimates from:
# - Multilevel model (partial pooling / shrinkage)
# - N=1 separate models (no pooling)
#
# Tasks:
# 1. Compute alternative significance tests for multilevel model
# 2. Compare amp_mean, amp_median, amp_ci_width, amp_hdr between methods
# 3. Pairwise distributions of estimates
# 4. Joint distribution plot with CIs
# 5. Shrinkage visualization
# =============================================================================

library(tidyverse)
library(here)
library(ks)
library(patchwork)

# -----------------------------------------------------------------------------
# 1. Load data
# -----------------------------------------------------------------------------

cat("Loading data...\n")

# N=1 estimates (already has alternative tests)
n1_ests <- readRDS(here("fits", "ests_level1_pa_n1_strategy1.rds"))
cat("  N=1 estimates:", nrow(n1_ests), "individuals,", ncol(n1_ests), "columns\n")

# N=1 base estimates (has phase information)
n1_ests_base <- readRDS(here("fits", "ests_level1_pa_n1_strategy1_base.rds"))
# Add phase columns from base file if not present
if (!"phi_circ_ci_width" %in% names(n1_ests)) {
  n1_ests <- n1_ests %>%
    left_join(
      n1_ests_base %>% select(id, starts_with("phi")),
      by = "id"
    )
  cat("  Added phase columns from base file\n")
}

# Multilevel estimates (base only, needs alternative tests)
ml_ests_base <- readRDS(here("fits", "ests_level1_pa_random_var.rds"))
cat("  Multilevel base estimates:", nrow(ml_ests_base), "individuals,", ncol(ml_ests_base), "columns\n")

# Multilevel draws (for computing alternative tests and CI bounds)
ml_draws <- readRDS(here("fits", "draws_pa_random_var.rds"))
cat("  Multilevel draws:", nrow(ml_draws), "rows\n")

# -----------------------------------------------------------------------------
# 2. Compute alternative tests for multilevel model
# -----------------------------------------------------------------------------

cat("\nComputing alternative tests for multilevel model...\n")

compute_alternative_tests <- function(draws_i, y_i = NULL) { # This function is now synchronized with n1_cosinor_strategy1.R
  # Extract C and S samples for this person
  C <- draws_i$co
  S <- draws_i$si
  A <- draws_i$amp
  M <- draws_i$mesor
  sigma <- draws_i$sigma
  n_draws <- length(C)

  mean_CS <- c(mean(C), mean(S))
  cov_CS <- cov(cbind(C, S))

  # Method B: Bivariate normal DENSITY method
  method_B_results <- tryCatch({ # This block is correct and complete
    cov_inv <- solve(cov_CS)
    cov_det <- det(cov_CS)
    norm_const <- 1 / (2 * pi * sqrt(cov_det))
    mahal_origin_sq <- as.numeric(t(mean_CS) %*% cov_inv %*% mean_CS)
    dens_origin <- norm_const * exp(-0.5 * mahal_origin_sq)
    dens_mode <- norm_const
    dens_ratio <- dens_origin / dens_mode
    prob_higher_dens <- 1 - exp(-0.5 * mahal_origin_sq)
    list( # Return all calculated values
      dens_origin = dens_origin,
      dens_mode = dens_mode,
      dens_ratio = dens_ratio,
      prob_higher_dens = prob_higher_dens,
      mahal_sq = mahal_origin_sq
    )
  }, error = function(e) {
    list(dens_origin = NA, dens_mode = NA, dens_ratio = NA,
         prob_higher_dens = NA, mahal_sq = NA)
  })

  # Method C: Quadrant method (CORRECTED)
  median_C <- median(C)
  median_S <- median(S)
  method_C_same_quadrant_prob <- mean(
    (sign(C) == sign(median_C)) & (sign(S) == sign(median_S))
  )

  # This logic correctly finds the quadrant opposite the median
  if (median_C >= 0 & median_S >= 0) {
    # Median in Q1 (+,+), origin in Q3 (-,-)
    method_C_origin_quadrant_prob <- mean(C < 0 & S < 0)
  } else if (median_C < 0 & median_S >= 0) {
    # Median in Q2 (-,+), origin in Q4 (+,-)
    method_C_origin_quadrant_prob <- mean(C > 0 & S < 0)
  } else if (median_C < 0 & median_S < 0) {
    # Median in Q3 (-,-), origin in Q1 (+,+)
    method_C_origin_quadrant_prob <- mean(C > 0 & S > 0)
  } else {
    # Median in Q4 (+,-), origin in Q2 (-,+)
    method_C_origin_quadrant_prob <- mean(C < 0 & S > 0)
  }

  # These are useful diagnostics but not the primary test
  method_C_cross_C_axis_prob <- mean(sign(C) != sign(median_C))
  method_C_cross_S_axis_prob <- mean(sign(S) != sign(median_S))

  # Method D: Univariate CIs
  C_ci_90 <- quantile(C, c(0.05, 0.95))
  C_ci_95 <- quantile(C, c(0.025, 0.975))
  C_ci_99 <- quantile(C, c(0.005, 0.995))
  S_ci_90 <- quantile(S, c(0.05, 0.95))
  S_ci_95 <- quantile(S, c(0.025, 0.975))
  S_ci_99 <- quantile(S, c(0.005, 0.995))
  C_excl_zero_95 <- (C_ci_95[1] > 0) | (C_ci_95[2] < 0)
  S_excl_zero_95 <- (S_ci_95[1] > 0) | (S_ci_95[2] < 0)
  C_excl_zero_90 <- (C_ci_90[1] > 0) | (C_ci_90[2] < 0)
  C_excl_zero_99 <- (C_ci_99[1] > 0) | (C_ci_99[2] < 0)
  S_excl_zero_90 <- (S_ci_90[1] > 0) | (S_ci_90[2] < 0)
  S_excl_zero_99 <- (S_ci_99[1] > 0) | (S_ci_99[2] < 0)

  # Method E: Joint Wald test
  wald_results <- tryCatch({
    cov_inv <- solve(cov_CS)
    wald_stat <- as.numeric(t(mean_CS) %*% cov_inv %*% mean_CS)
    list(
      stat = wald_stat,
      pval = 1 - pchisq(wald_stat, df = 2)
    )
  }, error = function(e) list(stat = NA, pval = NA))
  method_E_wald_stat <- wald_results$stat # This was missing
  method_E_pval <- wald_results$pval     # This was missing

  # Method F: Amplitude Z-score
  A_mean <- mean(A)
  A_se <- sd(A)
  method_F_zscore <- A_mean / A_se
  sigma_C <- sd(C) # This block for the corrected Z-score was missing
  sigma_S <- sd(S)
  sigma_CS_pooled <- sqrt((sigma_C^2 + sigma_S^2) / 2)
  A_expected_null <- sigma_CS_pooled * sqrt(pi / 2)
  A_sd_null <- sigma_CS_pooled * sqrt(2 - pi / 2)
  method_F_rayleigh_z <- (A_mean - A_expected_null) / A_sd_null
  A_ci_95 <- quantile(A, c(0.025, 0.975))
  method_F_ci_lower_over_mean <- A_ci_95[1] / A_mean # This was missing
  method_F_cv <- A_se / A_mean # This was missing
  A_ci_90 <- quantile(A, c(0.05, 0.95)) # This was missing
  A_ci_99 <- quantile(A, c(0.005, 0.995)) # This was missing

  # R² for practical significance (WITHIN-PERSON variance explained)
  # Note: Mesor (M) is NOT included because it's constant within each person.
  # This R² answers: "Of this person's moment-to-moment variability, how much
  # is explained by their circadian rhythm?" See docs/r2_theory_and_implementation.md
  R2_draws <- (A^2 / 2) / (A^2 / 2 + sigma^2)
  R2_mean <- mean(R2_draws)
  R2_median <- median(R2_draws)
  R2_ci <- quantile(R2_draws, c(0.025, 0.975)) # This was missing

  tibble(
    # Method B (Corrected to return all values)
    method_B_dens_origin = method_B_results$dens_origin,
    method_B_dens_mode = method_B_results$dens_mode,
    method_B_dens_ratio = method_B_results$dens_ratio,
    method_B_prob_higher_dens = method_B_results$prob_higher_dens,
    method_B_mahal_sq = method_B_results$mahal_sq,
    # Method C (Corrected to return all values)
    method_C_origin_quadrant_prob = method_C_origin_quadrant_prob,
    method_C_same_quadrant_prob = method_C_same_quadrant_prob,
    method_C_cross_C_axis_prob = method_C_cross_C_axis_prob,
    method_C_cross_S_axis_prob = method_C_cross_S_axis_prob,
    # Method D (Corrected to return all values)
    method_D_C_ci95_lower = C_ci_95[1],
    method_D_C_ci95_upper = C_ci_95[2],
    method_D_S_ci95_lower = S_ci_95[1],
    method_D_S_ci95_upper = S_ci_95[2],
    method_D_C_excl_zero_95 = C_excl_zero_95,
    method_D_S_excl_zero_95 = S_excl_zero_95,
    method_D_C_excl_zero_90 = C_excl_zero_90,
    method_D_C_excl_zero_99 = C_excl_zero_99,
    method_D_S_excl_zero_90 = S_excl_zero_90,
    method_D_S_excl_zero_99 = S_excl_zero_99,
    # Method E (Corrected to return all values)
    method_E_wald_stat = method_E_wald_stat,
    method_E_pval = method_E_pval,
    # Method F (Corrected to return all values)
    method_F_amp_mean = A_mean,
    method_F_amp_se = A_se,
    method_F_zscore = method_F_zscore,
    method_F_rayleigh_z = method_F_rayleigh_z,
    method_F_cv = method_F_cv,
    method_F_ci_lower_ratio = method_F_ci_lower_over_mean,
    # Amplitude CIs (Corrected to return all values)
    amp_ci90_lower = A_ci_90[1],
    amp_ci90_upper = A_ci_90[2],
    amp_ci95_lower = A_ci_95[1],
    amp_ci95_upper = A_ci_95[2],
    amp_ci99_lower = A_ci_99[1],
    amp_ci99_upper = A_ci_99[2],
    # R-squared (Corrected to return all values)
    R2_mean = R2_mean,
    R2_median = R2_median,
    R2_ci_lower = R2_ci[1],
    R2_ci_upper = R2_ci[2]
  )
}

ml_alt_tests <- ml_draws %>%
  group_by(id) %>%
  group_modify(~ compute_alternative_tests(.x)) %>%
  ungroup()

cat("  Computed for", nrow(ml_alt_tests), "individuals\n")

# Merge with base estimates
ml_ests <- ml_ests_base %>%
  left_join(ml_alt_tests, by = "id")

cat("  Multilevel full estimates:", ncol(ml_ests), "columns\n")

# -----------------------------------------------------------------------------
# 3. Create comparison dataset
# -----------------------------------------------------------------------------

cat("\nCreating comparison dataset...\n")
# Select key columns and add method indicator
n1_compare <- n1_ests %>%
  # Select all relevant columns to ensure they are carried over
  select(id, item, starts_with("amp_"), starts_with("method_"), starts_with("R2_"),
         phi_circ_ci_width) %>%
  mutate(method = "N=1")

ml_compare <- ml_ests %>%
  # Select all relevant columns to ensure they are carried over
  select(id, item, starts_with("amp_"), starts_with("method_"), starts_with("R2_"),
         phi_circ_ci_width) %>%
  mutate(method = "Multilevel")

# Wide format for pairwise comparisons
compare_wide <- n1_compare %>%
  select(-method) %>%
  rename_with(~ paste0(.x, "_n1"), -c(id, item)) %>%
  left_join(
    ml_compare %>%
      select(-method) %>%
      rename_with(~ paste0(.x, "_ml"), -c(id, item)),
    by = c("id", "item") # Join by item as well for safety
  )

cat("  Comparison dataset:", nrow(compare_wide), "individuals\n")

# -----------------------------------------------------------------------------
# 4. Summary statistics
# -----------------------------------------------------------------------------

cat("\n=== Comparison Summary ===\n\n")

# Amplitude median comparison
cat("Amplitude median:\n")
cat("  N=1:        mean =", round(mean(compare_wide$amp_median_n1), 3),
    ", SD =", round(sd(compare_wide$amp_median_n1), 3), "\n")
cat("  Multilevel: mean =", round(mean(compare_wide$amp_median_ml), 3),
    ", SD =", round(sd(compare_wide$amp_median_ml), 3), "\n")
cat("  Correlation:", round(cor(compare_wide$amp_median_n1, compare_wide$amp_median_ml), 3), "\n\n")

# CI width comparison
cat("95% CI width:\n")
cat("  N=1:        mean =", round(mean(compare_wide$amp_ci_width_n1), 3),
    ", SD =", round(sd(compare_wide$amp_ci_width_n1), 3), "\n")
cat("  Multilevel: mean =", round(mean(compare_wide$amp_ci_width_ml), 3),
    ", SD =", round(sd(compare_wide$amp_ci_width_ml), 3), "\n")
cat("  N=1 CIs wider on average:",
    round(mean(compare_wide$amp_ci_width_n1 > compare_wide$amp_ci_width_ml) * 100, 1), "%\n\n")

# HDR comparison
cat("HDR > 95% (significant amplitude):\n")
cat("  N=1:       ", sum(compare_wide$amp_hdr_n1 > 95), "/", nrow(compare_wide),
    "(", round(mean(compare_wide$amp_hdr_n1 > 95) * 100, 1), "%)\n")
cat("  Multilevel:", sum(compare_wide$amp_hdr_ml > 95), "/", nrow(compare_wide),
    "(", round(mean(compare_wide$amp_hdr_ml > 95) * 100, 1), "%)\n")
cat("  Agreement: ", sum((compare_wide$amp_hdr_n1 > 95) == (compare_wide$amp_hdr_ml > 95)),
    "/", nrow(compare_wide), "\n\n")

# -----------------------------------------------------------------------------
# 5. Plots
# -----------------------------------------------------------------------------

cat("Creating plots...\n")

# Load viridis for consistent coloring
library(viridis)

# Theme for plots
theme_comparison <- theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 12)
  )

# Add significance indicator for border coloring
# Use average HDR across both models for coloring, and flag significant cases
compare_wide <- compare_wide %>%
  mutate(
    # Average HDR for fill color
    avg_hdr = (amp_hdr_n1 + amp_hdr_ml) / 2,
    # Significance flags
    sig_n1 = amp_hdr_n1 > 95,
    sig_ml = amp_hdr_ml > 95,
    # Significance categories for border color
    sig_category = case_when(
      sig_n1 & sig_ml ~ "both",      # Black border
      sig_n1 & !sig_ml ~ "n1_only",  # Red border
      !sig_n1 & sig_ml ~ "ml_only",  # Blue border
      TRUE ~ "neither"               # Gray border
    )
  )

# --- Plot 1: Amplitude median comparison with CIs (line segments) ---
# CI segments colored by HDR
p1 <- ggplot(compare_wide, aes(x = amp_median_ml, y = amp_median_n1)) +
  # Horizontal CI segments (Multilevel CI on x-axis) - colored by ML HDR
  geom_segment(
    aes(x = amp_ci95_lower_ml, xend = amp_ci95_upper_ml,
        y = amp_median_n1, yend = amp_median_n1,
        color = amp_hdr_ml),
    alpha = 0.5, linewidth = 0.5
  ) +
  # Vertical CI segments (N=1 CI on y-axis) - colored by N=1 HDR
  geom_segment(
    aes(x = amp_median_ml, xend = amp_median_ml,
        y = amp_ci95_lower_n1, yend = amp_ci95_upper_n1,
        color = amp_hdr_n1),
    alpha = 0.5, linewidth = 0.5
  ) +
  # Point estimates - neither significant (gray border)
  geom_point(
    data = compare_wide %>% filter(sig_category == "neither"),
    aes(fill = avg_hdr),
    shape = 21, size = 2.5, alpha = 0.8, stroke = 0.3, color = "gray50"
  ) +
  # Point estimates - N=1 only significant (red border)
  geom_point(
    data = compare_wide %>% filter(sig_category == "n1_only"),
    aes(fill = avg_hdr),
    shape = 21, size = 2.5, alpha = 0.9, stroke = 0.8, color = "red"
  ) +
  # Point estimates - ML only significant (blue border)
  geom_point(
    data = compare_wide %>% filter(sig_category == "ml_only"),
    aes(fill = avg_hdr),
    shape = 21, size = 2.5, alpha = 0.9, stroke = 0.8, color = "blue"
  ) +
  # Point estimates - both significant (black border)
  geom_point(
    data = compare_wide %>% filter(sig_category == "both"),
    aes(fill = avg_hdr),
    shape = 21, size = 2.5, alpha = 0.9, stroke = 0.8, color = "black"
  ) +
  # Identity line
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  # Color scales
  scale_color_viridis_c(name = "HDR %", begin = 0.3, end = 0.95, guide = "none") +
  scale_fill_viridis_c(name = "HDR %", begin = 0.3, end = 0.95) +
  # Labels
  labs(
    x = "Multilevel amplitude (median)",
    y = "N=1 amplitude (median)",
    title = "Amplitude estimates: Multilevel vs N=1",
    subtitle = "Border: black=both sig, red=N=1 only, blue=ML only (HDR > 95%)"
  ) +
  theme_comparison

# --- Plot 2: CI width comparison ---
p2 <- ggplot(compare_wide, aes(x = amp_ci_width_ml, y = amp_ci_width_n1)) +
  # Neither significant (gray border)
  geom_point(
    data = compare_wide %>% filter(sig_category == "neither"),
    aes(fill = avg_hdr),
    shape = 21, size = 2.5, alpha = 0.8, stroke = 0.3, color = "gray50"
  ) +
  # N=1 only significant (red border)
  geom_point(
    data = compare_wide %>% filter(sig_category == "n1_only"),
    aes(fill = avg_hdr),
    shape = 21, size = 2.5, alpha = 0.9, stroke = 0.8, color = "red"
  ) +
  # ML only significant (blue border)
  geom_point(
    data = compare_wide %>% filter(sig_category == "ml_only"),
    aes(fill = avg_hdr),
    shape = 21, size = 2.5, alpha = 0.9, stroke = 0.8, color = "blue"
  ) +
  # Both significant (black border)
  geom_point(
    data = compare_wide %>% filter(sig_category == "both"),
    aes(fill = avg_hdr),
    shape = 21, size = 2.5, alpha = 0.9, stroke = 0.8, color = "black"
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  scale_fill_viridis_c(name = "HDR %", begin = 0.3, end = 0.95) +
  labs(
    x = "Multilevel 95% CI width",
    y = "N=1 95% CI width",
    title = "Amplitude CI width comparison",
    subtitle = "Points above line = N=1 has wider CIs"
  ) +
  theme_comparison

# --- Plot: Phase CI width comparison ---
# Convert from radians to hours: multiply by 12/π (since 2π radians = 24 hours)
p_phi_ci_width <- ggplot(compare_wide,
                         aes(x = phi_circ_ci_width_ml * 12 / pi,
                             y = phi_circ_ci_width_n1 * 12 / pi)) +
  # Neither significant (gray border)
  geom_point(
    data = compare_wide %>% filter(sig_category == "neither"),
    aes(fill = avg_hdr),
    shape = 21, size = 2.5, alpha = 0.8, stroke = 0.3, color = "gray50"
  ) +
  # N=1 only significant (red border)
  geom_point(
    data = compare_wide %>% filter(sig_category == "n1_only"),
    aes(fill = avg_hdr),
    shape = 21, size = 2.5, alpha = 0.9, stroke = 0.8, color = "red"
  ) +
  # ML only significant (blue border)
  geom_point(
    data = compare_wide %>% filter(sig_category == "ml_only"),
    aes(fill = avg_hdr),
    shape = 21, size = 2.5, alpha = 0.9, stroke = 0.8, color = "blue"
  ) +
  # Both significant (black border)
  geom_point(
    data = compare_wide %>% filter(sig_category == "both"),
    aes(fill = avg_hdr),
    shape = 21, size = 2.5, alpha = 0.9, stroke = 0.8, color = "black"
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  scale_fill_viridis_c(name = "HDR %", begin = 0.3, end = 0.95) +
  labs(
    x = "Multilevel phase 95% CI width (hours)",
    y = "N=1 phase 95% CI width (hours)",
    title = "Phase CI width comparison",
    subtitle = "Points above line = N=1 has wider CIs"
  ) +
  theme_comparison

# --- Plot 3: HDR comparison ---
p3 <- ggplot(compare_wide, aes(x = amp_hdr_ml, y = amp_hdr_n1)) +
  # Neither significant (gray border)
  geom_point(
    data = compare_wide %>% filter(sig_category == "neither"),
    aes(fill = avg_hdr),
    shape = 21, size = 2.5, alpha = 0.8, stroke = 0.3, color = "gray50"
  ) +
  # N=1 only significant (red border)
  geom_point(
    data = compare_wide %>% filter(sig_category == "n1_only"),
    aes(fill = avg_hdr),
    shape = 21, size = 2.5, alpha = 0.9, stroke = 0.8, color = "red"
  ) +
  # ML only significant (blue border)
  geom_point(
    data = compare_wide %>% filter(sig_category == "ml_only"),
    aes(fill = avg_hdr),
    shape = 21, size = 2.5, alpha = 0.9, stroke = 0.8, color = "blue"
  ) +
  # Both significant (black border)
  geom_point(
    data = compare_wide %>% filter(sig_category == "both"),
    aes(fill = avg_hdr),
    shape = 21, size = 2.5, alpha = 0.9, stroke = 0.8, color = "black"
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  geom_hline(yintercept = 95, linetype = "dotted", color = "darkgreen") +
  geom_vline(xintercept = 95, linetype = "dotted", color = "darkgreen") +
  scale_fill_viridis_c(name = "HDR %", begin = 0.3, end = 0.95) +
  labs(
    x = "Multilevel HDR (%)",
    y = "N=1 HDR (%)",
    title = "HDR comparison",
    subtitle = "Green lines at 95% threshold"
  ) +
  coord_fixed(xlim = c(0, 100), ylim = c(0, 100)) +
  theme_comparison

# --- Plot 4: Shrinkage plot ---
# Calculate group mean amplitude from multilevel draws
ml_group_mean <- ml_draws %>%
  group_by(iteration) %>%
  summarise(amp = mean(amp)) %>%
  summarise(group_amp_median = median(amp)) %>%
  pull(group_amp_median)

shrinkage_data <- compare_wide %>%
  mutate(
    shrinkage = amp_median_n1 - amp_median_ml,
    distance_from_mean = amp_median_n1 - ml_group_mean
  )

p4 <- ggplot(shrinkage_data, aes(x = amp_median_n1, y = amp_median_ml)) +
  # Vertical segments showing shrinkage (from identity line to ML estimate)
  geom_segment(
    aes(xend = amp_median_n1, yend = amp_median_n1),
    color = "gray60", alpha = 0.5
  ) +
  # Neither significant (gray border)
  geom_point(
    data = shrinkage_data %>% filter(sig_category == "neither"),
    aes(fill = avg_hdr),
    shape = 21, size = 2.5, alpha = 0.8, stroke = 0.3, color = "gray50"
  ) +
  # N=1 only significant (red border)
  geom_point(
    data = shrinkage_data %>% filter(sig_category == "n1_only"),
    aes(fill = avg_hdr),
    shape = 21, size = 2.5, alpha = 0.9, stroke = 0.8, color = "red"
  ) +
  # ML only significant (blue border)
  geom_point(
    data = shrinkage_data %>% filter(sig_category == "ml_only"),
    aes(fill = avg_hdr),
    shape = 21, size = 2.5, alpha = 0.9, stroke = 0.8, color = "blue"
  ) +
  # Both significant (black border)
  geom_point(
    data = shrinkage_data %>% filter(sig_category == "both"),
    aes(fill = avg_hdr),
    shape = 21, size = 2.5, alpha = 0.9, stroke = 0.8, color = "black"
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  geom_hline(yintercept = ml_group_mean, linetype = "dotted", color = "purple", linewidth = 0.8) +
  scale_fill_viridis_c(name = "HDR %", begin = 0.3, end = 0.95) +
  labs(
    x = "N=1 amplitude (median)",
    y = "Multilevel amplitude (median)",
    title = "Shrinkage toward group mean",
    subtitle = paste0("Purple line = group mean (", round(ml_group_mean, 3),
                      "); vertical lines show shrinkage")
  ) +
  coord_fixed() +
  theme_comparison

# --- Plot 5: R² comparison ---
p5 <- ggplot(compare_wide, aes(x = R2_median_ml, y = R2_median_n1)) +
  # Neither significant (gray border)
  geom_point(
    data = compare_wide %>% filter(sig_category == "neither"),
    aes(fill = avg_hdr),
    shape = 21, size = 2.5, alpha = 0.8, stroke = 0.3, color = "gray50"
  ) +
  # N=1 only significant (red border)
  geom_point(
    data = compare_wide %>% filter(sig_category == "n1_only"),
    aes(fill = avg_hdr),
    shape = 21, size = 2.5, alpha = 0.9, stroke = 0.8, color = "red"
  ) +
  # ML only significant (blue border)
  geom_point(
    data = compare_wide %>% filter(sig_category == "ml_only"),
    aes(fill = avg_hdr),
    shape = 21, size = 2.5, alpha = 0.9, stroke = 0.8, color = "blue"
  ) +
  # Both significant (black border)
  geom_point(
    data = compare_wide %>% filter(sig_category == "both"),
    aes(fill = avg_hdr),
    shape = 21, size = 2.5, alpha = 0.9, stroke = 0.8, color = "black"
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  geom_hline(yintercept = 0.10, linetype = "dotted", color = "darkgreen") +
  geom_vline(xintercept = 0.10, linetype = "dotted", color = "darkgreen") +
  scale_fill_viridis_c(name = "HDR %", begin = 0.3, end = 0.95) +
  labs(
    x = "Multilevel R² (median)",
    y = "N=1 R² (median)",
    title = "R² (practical significance)",
    subtitle = "Green lines at 10% variance explained"
  ) +
  coord_fixed(xlim = c(0, max(c(compare_wide$R2_median_ml, compare_wide$R2_median_n1)) * 1.1),
              ylim = c(0, max(c(compare_wide$R2_median_ml, compare_wide$R2_median_n1)) * 1.1)) +
  theme_comparison

# Combine plots (excluding shrinkage plot p4)
combined <- (p1 | p2) / (p3 | p5) +
  plot_annotation(
    title = "Multilevel vs N=1 Cosinor: Individual Amplitude Estimates",
    subtitle = "Border: black=both sig, red=N=1 only, blue=ML only (HDR > 95%)",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11)
    )
  )

# Save
ggsave(here("figs", "multilevel_vs_n1_comparison.pdf"),
       combined, width = 12, height = 10)
ggsave(here("figs", "multilevel_vs_n1_comparison.png"),
       combined, width = 12, height = 10, dpi = 150)

cat("  Saved to figs/multilevel_vs_n1_comparison.pdf\n")

# Also save individual key plots
ggsave(here("figs", "amp_median_ml_vs_n1.pdf"), p1, width = 7, height = 7)
cat("  Saved to figs/amp_median_ml_vs_n1.pdf\n")

# Save shrinkage plot separately
ggsave(here("figs", "shrinkage_plot.pdf"), p4, width = 8, height = 8)
ggsave(here("figs", "shrinkage_plot.png"), p4, width = 8, height = 8, dpi = 150)
cat("  Saved to figs/shrinkage_plot.pdf\n")

# Save phase CI width plot separately
ggsave(here("figs", "phi_ci_width_ml_vs_n1.pdf"), p_phi_ci_width, width = 7, height = 7)
ggsave(here("figs", "phi_ci_width_ml_vs_n1.png"), p_phi_ci_width, width = 7, height = 7, dpi = 150)
cat("  Saved to figs/phi_ci_width_ml_vs_n1.pdf\n")

# Combined 6-panel plot: (p1 | p2) / (p4 | p_phi_ci_width) / (p3 | p5)
# Use design matrix with a spacer column (#) for horizontal breathing room
design <- "
A#B
C#D
E#F
"
combined_6 <- wrap_plots(
  A = p1, B = p2,
  C = p4, D = p_phi_ci_width,
  E = p3, F = p5,
  design = design,
  guides = "collect"
) +
  plot_layout(widths = c(1, 0.05, 1)) +
  plot_annotation(
    title = "Multilevel vs N=1 Cosinor: Comprehensive Comparison",
    subtitle = "Border: black=both sig, red=N=1 only, blue=ML only (HDR > 95%)",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11)
    )
  ) &
  theme(legend.position = "right")

ggsave(here("figs", "multilevel_vs_n1_comparison_6panel.pdf"),
       combined_6, width = 12, height = 15)
ggsave(here("figs", "multilevel_vs_n1_comparison_6panel.png"),
       combined_6, width = 12, height = 15, dpi = 150)
cat("  Saved to figs/multilevel_vs_n1_comparison_6panel.pdf\n")

# -----------------------------------------------------------------------------
# 6. Significance method agreement table
# -----------------------------------------------------------------------------
#
# NOTE on method interpretation (corrected per reviewer's intent):
#   A. HDR: HIGH value (>95) = significant (origin excluded from 95% HDR)
#   B. BivNorm: HIGH prob_higher_dens (>0.95) = significant (origin in low-density region)
#   C. Quadrant: LOW origin_quadrant_prob (<0.05) = significant (origin unlikely)
#   D. Univariate CI: CI excludes 0 = significant
#   E. Wald: LOW p-value (<0.05) = significant
#   F. Z-score: HIGH Z (>1.96) = significant (but see theoretical caveats)
# -----------------------------------------------------------------------------

cat("\n=== Significance Test Agreement (at alpha = 0.05) ===\n\n")

# Check which columns are available - N=1 may have different column names
available_cols <- names(compare_wide)
has_prob_higher_dens_n1 <- "method_B_prob_higher_dens_n1" %in% available_cols
has_prob_higher_dens_ml <- "method_B_prob_higher_dens_ml" %in% available_cols
has_rayleigh_n1 <- "method_F_rayleigh_z_n1" %in% available_cols
has_rayleigh_ml <- "method_F_rayleigh_z_ml" %in% available_cols

sig_compare <- compare_wide %>%
  transmute(
    id = id,
    # Method A: HDR > 95 (HIGH = significant)
    hdr_n1 = amp_hdr_n1 > 95,
    hdr_ml = amp_hdr_ml > 95,

    # Method B: Bivariate normal density - prob_higher_dens > 0.95 (HIGH = significant)
    # This is analogous to HDR level
    # Fall back to method_B_pval < 0.05 if prob_higher_dens not available
    bivnorm_n1 = if (has_prob_higher_dens_n1) method_B_prob_higher_dens_n1 > 0.95
                 else if ("method_B_pval_n1" %in% available_cols) method_B_pval_n1 < 0.05
                 else NA,
    bivnorm_ml = if (has_prob_higher_dens_ml) method_B_prob_higher_dens_ml > 0.95
                 else if ("method_B_pval_ml" %in% available_cols) method_B_pval_ml < 0.05
                 else NA,

    # Method C: Quadrant - origin_quadrant_prob < 0.05 (LOW = significant)
    quad_n1 = method_C_origin_quadrant_prob_n1 < 0.05,
    quad_ml = method_C_origin_quadrant_prob_ml < 0.05,

    # Method D: C or S 95% CI excludes zero
    univar_n1 = method_D_C_excl_zero_95_n1 | method_D_S_excl_zero_95_n1,
    univar_ml = method_D_C_excl_zero_95_ml | method_D_S_excl_zero_95_ml,

    # Method E: Wald p < 0.05 (LOW = significant)
    wald_n1 = method_E_pval_n1 < 0.05,
    wald_ml = method_E_pval_ml < 0.05,

    # Method F: Z > 1.96 (original, has theoretical issues)
    zscore_n1 = method_F_zscore_n1 > 1.96,
    zscore_ml = method_F_zscore_ml > 1.96,

    # Method F (Rayleigh-corrected): Better theoretical basis
    rayleigh_n1 = if (has_rayleigh_n1) method_F_rayleigh_z_n1 > 1.96 else NA,
    rayleigh_ml = if (has_rayleigh_ml) method_F_rayleigh_z_ml > 1.96 else NA
  )

# Agreement summary
methods <- c("hdr", "bivnorm", "quad", "univar", "wald", "zscore", "rayleigh")
method_names <- c("HDR > 95%", "BivNorm (dens)", "Quadrant < 5%",
                  "Univar CI (C|S)", "Wald test", "Z-score > 1.96", "Rayleigh-Z > 1.96")

n_total <- nrow(sig_compare)

# Counts table
agreement_table_counts <- tibble(
  Method = method_names,
  `N=1 sig` = sapply(methods, function(m) sum(sig_compare[[paste0(m, "_n1")]], na.rm = TRUE)),
  `ML sig` = sapply(methods, function(m) sum(sig_compare[[paste0(m, "_ml")]], na.rm = TRUE)),
  `Both sig` = sapply(methods, function(m) {
    sum(sig_compare[[paste0(m, "_n1")]] & sig_compare[[paste0(m, "_ml")]], na.rm = TRUE)
  }),
  `Neither sig` = sapply(methods, function(m) {
    sum(!sig_compare[[paste0(m, "_n1")]] & !sig_compare[[paste0(m, "_ml")]], na.rm = TRUE)
  }),
  `N=1 only` = sapply(methods, function(m) {
    sum(sig_compare[[paste0(m, "_n1")]] & !sig_compare[[paste0(m, "_ml")]], na.rm = TRUE)
  }),
  `ML only` = sapply(methods, function(m) {
    sum(!sig_compare[[paste0(m, "_n1")]] & sig_compare[[paste0(m, "_ml")]], na.rm = TRUE)
  }),
  Agreement = sapply(methods, function(m) {
    sum(sig_compare[[paste0(m, "_n1")]] == sig_compare[[paste0(m, "_ml")]], na.rm = TRUE)
  })
)

cat("--- Counts (N =", n_total, ") ---\n")
print(agreement_table_counts)

# Percentages table
agreement_table_pct <- tibble(
  Method = method_names,
  `N=1 sig %` = sapply(methods, function(m) {
    round(sum(sig_compare[[paste0(m, "_n1")]], na.rm = TRUE) / n_total * 100, 1)
  }),
  `ML sig %` = sapply(methods, function(m) {
    round(sum(sig_compare[[paste0(m, "_ml")]], na.rm = TRUE) / n_total * 100, 1)
  }),
  `Both %` = sapply(methods, function(m) {
    round(sum(sig_compare[[paste0(m, "_n1")]] & sig_compare[[paste0(m, "_ml")]], na.rm = TRUE) / n_total * 100, 1)
  }),
  `Neither %` = sapply(methods, function(m) {
    round(sum(!sig_compare[[paste0(m, "_n1")]] & !sig_compare[[paste0(m, "_ml")]], na.rm = TRUE) / n_total * 100, 1)
  }),
  `N=1 only %` = sapply(methods, function(m) {
    round(sum(sig_compare[[paste0(m, "_n1")]] & !sig_compare[[paste0(m, "_ml")]], na.rm = TRUE) / n_total * 100, 1)
  }),
  `ML only %` = sapply(methods, function(m) {
    round(sum(!sig_compare[[paste0(m, "_n1")]] & sig_compare[[paste0(m, "_ml")]], na.rm = TRUE) / n_total * 100, 1)
  }),
  `Agreement %` = sapply(methods, function(m) {
    round(sum(sig_compare[[paste0(m, "_n1")]] == sig_compare[[paste0(m, "_ml")]], na.rm = TRUE) / n_total * 100, 1)
  })
)

cat("\n--- Percentages ---\n")
print(agreement_table_pct)

# Compact summary: disagreement breakdown
cat("\n--- Disagreement Summary ---\n")
cat("(When methods disagree, which model finds significance?)\n\n")
disagree_summary <- tibble(
  Method = method_names,
  `Disagree (n)` = agreement_table_counts$`N=1 only` + agreement_table_counts$`ML only`,
  `N=1 finds sig` = paste0(agreement_table_counts$`N=1 only`, " (",
                           round(agreement_table_counts$`N=1 only` / (agreement_table_counts$`N=1 only` + agreement_table_counts$`ML only` + 0.001) * 100, 0), "%)"),
  `ML finds sig` = paste0(agreement_table_counts$`ML only`, " (",
                          round(agreement_table_counts$`ML only` / (agreement_table_counts$`N=1 only` + agreement_table_counts$`ML only` + 0.001) * 100, 0), "%)")
)
print(disagree_summary)

# -----------------------------------------------------------------------------
# 7. Within-model method comparison matrices
# -----------------------------------------------------------------------------
#
# For each model (N=1 and Multilevel), create a matrix showing how often
# different significance methods agree with each other.
# -----------------------------------------------------------------------------

cat("\n=== Within-Model Method Comparison Matrices ===\n")

# Helper function to compute pairwise agreement matrix
compute_agreement_matrix <- function(sig_data, method_cols) {
  n_methods <- length(method_cols)
  mat <- matrix(NA, nrow = n_methods, ncol = n_methods)
  rownames(mat) <- colnames(mat) <- method_cols

  for (i in 1:n_methods) {
    for (j in 1:n_methods) {
      if (i == j) {
        mat[i, j] <- sum(sig_data[[method_cols[i]]], na.rm = TRUE)  # Count of "significant"
      } else {
        mat[i, j] <- sum(sig_data[[method_cols[i]]] == sig_data[[method_cols[j]]], na.rm = TRUE)
      }
    }
  }
  return(mat)
}

# Create significance indicators for N=1 (at 95% level)
# Using CORRECTED method definitions
sig_n1 <- compare_wide %>%
  transmute(
    HDR = amp_hdr_n1 > 95,                                    # HIGH = sig
    BivNorm = if (has_prob_higher_dens_n1) method_B_prob_higher_dens_n1 > 0.95
              else if ("method_B_pval_n1" %in% available_cols) method_B_pval_n1 < 0.05
              else NA,
    Quadrant = method_C_origin_quadrant_prob_n1 < 0.05,       # LOW = sig (CORRECTED)
    UnivarCI = method_D_C_excl_zero_95_n1 | method_D_S_excl_zero_95_n1,
    Wald = method_E_pval_n1 < 0.05,                           # LOW = sig
    Zscore = method_F_zscore_n1 > 1.96,                       # HIGH = sig (has issues)
    RayleighZ = if (has_rayleigh_n1) method_F_rayleigh_z_n1 > 1.96 else NA
  )

# Create significance indicators for Multilevel (at 95% level)
sig_ml <- compare_wide %>%
  transmute(
    HDR = amp_hdr_ml > 95,
    BivNorm = if (has_prob_higher_dens_ml) method_B_prob_higher_dens_ml > 0.95
              else if ("method_B_pval_ml" %in% available_cols) method_B_pval_ml < 0.05
              else NA,
    Quadrant = method_C_origin_quadrant_prob_ml < 0.05,       # LOW = sig (CORRECTED)
    UnivarCI = method_D_C_excl_zero_95_ml | method_D_S_excl_zero_95_ml,
    Wald = method_E_pval_ml < 0.05,
    Zscore = method_F_zscore_ml > 1.96,
    RayleighZ = if (has_rayleigh_ml) method_F_rayleigh_z_ml > 1.96 else NA
  )

method_names_short <- c("HDR", "BivNorm", "Quadrant", "UnivarCI", "Wald", "Zscore", "RayleighZ")

# N=1 agreement matrix
n1_agree_mat <- compute_agreement_matrix(sig_n1, method_names_short)
cat("\nN=1 Model - Method Agreement Matrix:\n")
cat("(Diagonal: count of 'significant'; Off-diagonal: count of matching decisions)\n")
print(n1_agree_mat)

# Multilevel agreement matrix
ml_agree_mat <- compute_agreement_matrix(sig_ml, method_names_short)
cat("\nMultilevel Model - Method Agreement Matrix:\n")
print(ml_agree_mat)

# Convert to agreement percentages (off-diagonal only)
n1_agree_pct <- round(n1_agree_mat / nrow(compare_wide) * 100, 1)
ml_agree_pct <- round(ml_agree_mat / nrow(compare_wide) * 100, 1)

cat("\nN=1 Model - Agreement % (diagonal = % significant):\n")
print(n1_agree_pct)

cat("\nMultilevel Model - Agreement % (diagonal = % significant):\n")
print(ml_agree_pct)

# -----------------------------------------------------------------------------
# 8. Within-model pairwise plots (using GGally::ggpairs)
# -----------------------------------------------------------------------------

cat("\nCreating within-model pairwise distribution plots...\n")

# Check if GGally is available
if (requireNamespace("GGally", quietly = TRUE)) {
  library(GGally)

  # Prepare data for N=1 pairwise plot
  pairplot_data_n1 <- compare_wide %>%
    select(
      `Amp median` = amp_median_n1,
      `CI width` = amp_ci_width_n1,
      `HDR %` = amp_hdr_n1,
      `R2` = R2_median_n1
    )

  # N=1 pairplot
  pp_n1 <- ggpairs(
    pairplot_data_n1,
    title = "N=1 Model: Pairwise Distributions",
    lower = list(continuous = wrap("points", alpha = 0.5, size = 1)),
    diag = list(continuous = wrap("densityDiag", alpha = 0.5)),
    upper = list(continuous = wrap("cor", size = 3))
  ) +
    theme_minimal(base_size = 10)

  ggsave(here("figs", "pairplot_n1_model.pdf"), pp_n1, width = 8, height = 8)
  ggsave(here("figs", "pairplot_n1_model.png"), pp_n1, width = 8, height = 8, dpi = 150)
  cat("  Saved N=1 pairplot to figs/pairplot_n1_model.pdf\n")

  # Prepare data for Multilevel pairwise plot
  pairplot_data_ml <- compare_wide %>%
    select(
      `Amp median` = amp_median_ml,
      `CI width` = amp_ci_width_ml,
      `HDR %` = amp_hdr_ml,
      `R2` = R2_median_ml
    )

  # Multilevel pairplot
  pp_ml <- ggpairs(
    pairplot_data_ml,
    title = "Multilevel Model: Pairwise Distributions",
    lower = list(continuous = wrap("points", alpha = 0.5, size = 1)),
    diag = list(continuous = wrap("densityDiag", alpha = 0.5)),
    upper = list(continuous = wrap("cor", size = 3))
  ) +
    theme_minimal(base_size = 10)

  ggsave(here("figs", "pairplot_ml_model.pdf"), pp_ml, width = 8, height = 8)
  ggsave(here("figs", "pairplot_ml_model.png"), pp_ml, width = 8, height = 8, dpi = 150)
  cat("  Saved Multilevel pairplot to figs/pairplot_ml_model.pdf\n")

  # Combined side-by-side (using patchwork if both are ggplot-compatible)
  # GGally::ggpairs returns a ggmatrix, which works with patchwork
  combined_pairplots <- (pp_n1 | pp_ml)
  ggsave(here("figs", "pairplot_both_models.pdf"), combined_pairplots, width = 16, height = 8)
  cat("  Saved combined pairplot to figs/pairplot_both_models.pdf\n")

} else {
  cat("  GGally package not available. Skipping pairplots.\n")
  cat("  Install with: install.packages('GGally')\n")
}

# -----------------------------------------------------------------------------
# 9. Method agreement heatmaps (visualization of matrices)
# -----------------------------------------------------------------------------

cat("\nCreating method agreement heatmaps...\n")

# Reshape agreement matrices for ggplot
reshape_for_heatmap <- function(mat, model_name) {
  as_tibble(mat, rownames = "Method1") %>%
    pivot_longer(-Method1, names_to = "Method2", values_to = "Agreement") %>%
    mutate(
      Model = model_name,
      # Make it a percentage
      Agreement_pct = Agreement / nrow(compare_wide) * 100,
      # Keep diagonal as count
      is_diagonal = Method1 == Method2
    )
}

heatmap_data <- bind_rows(
  reshape_for_heatmap(n1_agree_mat, "N=1"),
  reshape_for_heatmap(ml_agree_mat, "Multilevel")
) %>%
  mutate(
    Method1 = factor(Method1, levels = method_names_short),
    Method2 = factor(Method2, levels = rev(method_names_short))
  )

# Create heatmap
p_heatmap <- ggplot(heatmap_data, aes(x = Method1, y = Method2, fill = Agreement_pct)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Agreement_pct, 1)), size = 3) +
  scale_fill_gradient2(
    low = "white", mid = "steelblue", high = "darkblue",
    midpoint = 50, limits = c(0, 100),
    name = "Agreement %"
  ) +
  facet_wrap(~ Model) +
  labs(
    title = "Within-Model Method Agreement",
    subtitle = "Diagonal: % significant; Off-diagonal: % agreement between methods",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) +
  coord_fixed()

ggsave(here("figs", "method_agreement_heatmap.pdf"), p_heatmap, width = 10, height = 5)
ggsave(here("figs", "method_agreement_heatmap.png"), p_heatmap, width = 10, height = 5, dpi = 150)
cat("  Saved method agreement heatmap to figs/method_agreement_heatmap.pdf\n")

# -----------------------------------------------------------------------------
# 10. Save comparison data
# -----------------------------------------------------------------------------

saveRDS(compare_wide, here("fits", "comparison_multilevel_n1_pa.rds"))
cat("\nSaved comparison data to fits/comparison_multilevel_n1_pa.rds\n")

cat("\nDone!\n")
beepr::beep(5)
