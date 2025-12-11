# =============================================================================
# Replicated N=1 Cosinor Model using Strategy 1 (Separate Models per Person)
# =============================================================================
#
# This script fits a separate cosinor model for each individual.
# Each person's model is completely independent - true N=1 analyses.
#
# Model for each person i:
#   y_ij ~ Normal(M_i + C_i * cos(2*pi*t/24) + S_i * sin(2*pi*t/24), sigma_i)
#
# Parallelization strategy:
#   - For testing: sequential loop, brm uses multiple cores for chains
#   - For production: parallel over individuals, brm runs chains sequentially
# =============================================================================

# -----------------------------------------------------------------------------
# 1. Load required packages
# -----------------------------------------------------------------------------

library(brms)
library(tidyverse)
library(here)
library(circular)

# For parallel execution over individuals (optional)
# library(furrr)
# library(future)

# -----------------------------------------------------------------------------
# 2. Load and prepare data
# -----------------------------------------------------------------------------

item_ <- "pa"
d <- readRDS(here("data",
                  "d_leuven_joined.rds")) %>%
  filter(item == item_)

# Get unique IDs
ids <- unique(d$id)
n_ids <- length(ids)
cat("Total number of individuals:", n_ids, "\n")

# -----------------------------------------------------------------------------
# 3. Settings
# -----------------------------------------------------------------------------

# For testing: use a subset of individuals
TEST_MODE <- FALSE
N_TEST <- 10  # Number of individuals for testing

if (TEST_MODE) {
  ids_to_fit <- ids[1:N_TEST]
  cat("TEST MODE: Fitting", N_TEST, "individuals\n")
} else {
  ids_to_fit <- ids
  cat("FULL MODE: Fitting all", n_ids, "individuals\n")
}

# Parallelization mode:
# - "chains": sequential over individuals, brm parallelizes chains (safer, good for testing)
# - "individuals": parallel over individuals, brm runs sequentially (faster for many individuals)
PARALLEL_MODE <- "chains"

# brm settings
ITER <- 3000
WARMUP <- 1500
CHAINS <- 4

# -----------------------------------------------------------------------------
# 4. Define the model fitting function
# -----------------------------------------------------------------------------

fit_single_person <- function(person_id, data, parallel_mode = "chains") {

  # Subset data for this person
  d_i <- data %>% filter(id == person_id)
  n_obs <- nrow(d_i)

  cat("Fitting person", person_id, "(", n_obs, "observations)...\n")

  # Set cores based on parallelization mode
  if (parallel_mode == "chains") {
    # Let brm parallelize chains
    n_cores <- CHAINS
  } else {
    # We're parallelizing over individuals, so brm runs sequentially
    n_cores <- 1
  }

  # Fit the model with person-specific sigma
  fit <- brm(
    brmsformula(
      y ~ 1 + co + si,
      sigma ~ 1  # Estimates person-specific sigma (just intercept = log(sigma_i))
    ),
    data = d_i,
    backend = "cmdstanr",
    chains = CHAINS,
    iter = ITER,
    warmup = WARMUP,
    cores = n_cores,
    silent = 2,        # Suppress most output
    refresh = 0        # Don't print iteration progress
  )

  return(fit)
}

# -----------------------------------------------------------------------------
# 5. Fit models
# -----------------------------------------------------------------------------

cat("\n=== Starting model fitting ===\n")
cat("Parallelization mode:", PARALLEL_MODE, "\n")
cat("Iterations:", ITER, "| Warmup:", WARMUP, "| Chains:", CHAINS, "\n\n")

(st <- Sys.time())

if (PARALLEL_MODE == "chains") {
  # -----------------------------
  # Sequential over individuals
  # -----------------------------
  fits <- list()

  for (i in seq_along(ids_to_fit)) {
    person_id <- ids_to_fit[i]
    cat(sprintf("[%d/%d] ", i, length(ids_to_fit)))
    fits[[as.character(person_id)]] <- fit_single_person(person_id, d, "chains")
  }

} else if (PARALLEL_MODE == "individuals") {
  # -----------------------------
  # Parallel over individuals
  # -----------------------------
  # Uncomment the following to enable:

  # library(furrr)
  # library(future)
  #
  # # Set up parallel backend (adjust workers based on your CPU)
  # plan(multisession, workers = 4)
  #
  # fits <- future_map(
  #   ids_to_fit,
  #   ~ fit_single_person(.x, d, "individuals"),
  #   .options = furrr_options(seed = TRUE),
  #   .progress = TRUE
  # )
  # names(fits) <- as.character(ids_to_fit)
  #
  # # Reset to sequential
  # plan(sequential)

  stop("Parallel mode not yet enabled. Uncomment the code above to use it.")
}

elapsed <- Sys.time() - st
cat("\n=== Model fitting complete ===\n")
cat("Total time:", round(as.numeric(elapsed, units = "mins"), 1), "minutes\n")
cat("Average per person:", round(as.numeric(elapsed, units = "secs") / length(ids_to_fit), 1), "seconds\n")

beepr::beep(1)

# Save the fits
# saveRDS(fits,
#         here("fits",
#              paste0("brms_", item_, "_n1_strategy1_fits.rds")))

# -----------------------------------------------------------------------------
# 6. Extract posterior samples from all fits
# -----------------------------------------------------------------------------

cat("\nExtracting posterior samples...\n")

# Function to extract draws from a single fit
extract_draws_single <- function(fit, person_id) {
  draws <- as_draws_df(fit) %>%
    mutate(
      iteration = row_number(),
      id = person_id,
      # Extract parameters
      mesor = b_Intercept,
      co = b_co,
      si = b_si,
      logsigma = b_sigma_Intercept,
      # Derived quantities
      sigma = exp(logsigma),
      sigma2 = sigma^2,
      amp = sqrt(co^2 + si^2),
      phi = atan2(si, co) %% (2 * pi),
      phi.begins6 = (phi - pi/2) %% (2 * pi),
      phi.begins12 = (phi - pi) %% (2 * pi)
    ) %>%
    select(iteration, id, mesor, co, si, logsigma, sigma, sigma2,
           amp, phi, phi.begins6, phi.begins12)

  return(draws)
}

# Extract from all fits
d_draws <- map2_dfr(
  fits,
  names(fits),
  ~ extract_draws_single(.x, as.integer(.y))
)

# Save draws
d_draws %>%
  mutate(item = item_, .before = 1) %>%
  saveRDS(here("fits",
               paste0("draws_", item_, "_n1_strategy1.rds")))

cat("Posterior samples extracted:", nrow(d_draws), "rows\n")

# -----------------------------------------------------------------------------
# 7. Compute Level-1 estimates
# -----------------------------------------------------------------------------

cat("Computing level-1 estimates...\n")

# Load ks for HDR computation
require(ks)

# Function to compute alpha-HDR level
compute_alpha_hdr <- function(d, x_0 = 0, y_0 = 0, n_grid = 200) {
  fhat <- kde(d, gridsize = c(n_grid, n_grid))

  dens_vec <- as.vector(fhat$estimate)
  dx <- diff(fhat$eval.points[[1]])[1]
  dy <- diff(fhat$eval.points[[2]])[1]

  dens_sorted <- sort(dens_vec, decreasing = TRUE)
  cumulative_mass <- cumsum(dens_sorted * dx * dy)
  cumulative_mass <- cumulative_mass / max(cumulative_mass)

  f_point <- predict(fhat, x = cbind(x_0, y_0))
  idx <- which(dens_sorted < f_point)[1]

  if (is.na(idx)) {
    alpha_hdr <- 1
  } else if (idx > 1) {
    alpha_hdr <- cumulative_mass[idx - 1]
  } else {
    alpha_hdr <- cumulative_mass[1]
  }

  return(round(alpha_hdr * 100, 1))
}

# Helper functions
compute_linear_stats <- function(data, variable) {
  tibble(
    mean = mean(data[[variable]], na.rm = TRUE),
    median = median(data[[variable]], na.rm = TRUE),
    ci_width = quantile(data[[variable]], 0.975, na.rm = TRUE) -
               quantile(data[[variable]], 0.025, na.rm = TRUE)
  ) %>%
    rename_with(~ paste0(variable, "_", .))
}

compute_circular_stats <- function(data, variable) {
  x <- circular(data[[variable]])
  ci_2.5 <- quantile.circular(x, 0.025)
  ci_97.5 <- quantile.circular(x, 0.975)

  tibble(
    circ_mean = (mean.circular(x) %% (2 * pi)) %>% as.numeric(),
    circ_median = (median.circular(x) %% (2 * pi)) %>% as.numeric(),
    circ_ci_2.5 = (ci_2.5 %% (2 * pi)) %>% as.numeric(),
    circ_ci_97.5 = (ci_97.5 %% (2 * pi)) %>% as.numeric(),
    circ_ci_width = ((ci_97.5 - ci_2.5) %% (2 * pi)) %>% as.numeric()
  ) %>%
    rename_with(~ paste0(variable, "_", .))
}

# Compute estimates
ests_level1 <- d_draws %>%
  group_by(id) %>%
  summarise(
    bind_cols(
      compute_linear_stats(cur_data(), "mesor"),
      compute_linear_stats(cur_data(), "co"),
      compute_linear_stats(cur_data(), "si"),
      compute_linear_stats(cur_data(), "amp"),
      compute_linear_stats(cur_data(), "logsigma"),
      compute_linear_stats(cur_data(), "sigma"),
      compute_linear_stats(cur_data(), "sigma2"),
      compute_linear_stats(cur_data(), "phi"),
      compute_linear_stats(cur_data(), "phi.begins6"),
      compute_linear_stats(cur_data(), "phi.begins12"),
      compute_circular_stats(cur_data(), "phi"),
      compute_circular_stats(cur_data(), "phi.begins6"),
      compute_circular_stats(cur_data(), "phi.begins12")
    ),
    amp_hdr = compute_alpha_hdr(cbind(cur_data()$co, cur_data()$si))
  ) %>%
  relocate(amp_hdr, .before = logsigma_mean)

# Save BASE estimates (before alternative tests are added)
# This prevents corruption from repeated left_joins
ests_level1 %>%
  mutate(item = item_, .before = 1) %>%
  saveRDS(here("fits",
               paste0("ests_level1_", item_, "_n1_strategy1_base.rds")))

# -----------------------------------------------------------------------------
# 8. Alternative methods for testing H0: A_i = 0 (Reviewer's suggestions)
# -----------------------------------------------------------------------------
#
# The reviewer (Tihomir) suggested comparing multiple methods for testing whether
# individual amplitudes are significantly different from zero. This section
# implements methods A-F from the review:
#
# A. HDR method (already computed): density at (0,0) in HPD area
#    - Criterion: HDR > 95% means (0,0) excluded from 95% HDR → significant
#
# B. Bivariate normal density method: evaluate density at (0,0) relative to
#    the fitted bivariate normal distribution
#    - Different from Wald test: compares density at origin to density at mode
#    - Criterion: low density ratio → significant
#
# C. Quadrant method (reviewer's "smallest quadrant"): P(draws in origin quadrant)
#    - The origin (0,0) lies in the quadrant opposite to the posterior median
#    - Criterion: LOW probability → significant (origin quadrant is unlikely)
#    - Note: This is the OPPOSITE logic from what was previously implemented!
#
# D. Univariate CI method: 95% CI for C excludes 0 OR 95% CI for S excludes 0
#    - Union test: reject if either marginal rejects
#    - Criterion: either CI excludes 0 → significant
#
# E. Joint Wald test: chi-square test for H0: C=0 AND S=0
#    - Uses Mahalanobis distance with posterior covariance
#    - Criterion: p < 0.05 → significant
#    - Note: Mathematically similar to B but with different theoretical basis
#
# F. Amplitude Z-score: mean(A) / SD(A)
#    - CAUTION: A = sqrt(C² + S²) follows Rice/Rayleigh distribution under null
#    - E[A] > 0 even when C=S=0, so this is NOT a proper Z-score
#    - We also compute corrected versions using Rayleigh distribution theory
#    - Criterion: Z > 1.96 (but see caveats about distribution)
#
# G. R-squared: variance explained by the cosinor component
#    - Practical significance measure, not a hypothesis test
#
# We compute raw values and apply multiple thresholds for comparison.
#
# NOTE: This function can also be applied to multilevel posterior draws!
# Just pass d_draws from the multilevel model instead of the N=1 model.
# -----------------------------------------------------------------------------

cat("\nComputing alternative significance tests (reviewer methods A-F) and R²...\n")

compute_alternative_tests <- function(draws_i, y_i = NULL) {

  # Extract C and S samples for this person
  C <- draws_i$co
  S <- draws_i$si
  A <- draws_i$amp
  M <- draws_i$mesor
  sigma <- draws_i$sigma
  n_draws <- length(C)

  # ---------------------------
  # Method A: HDR (computed separately, values stored for reference)
  # ---------------------------
  # amp_hdr is computed elsewhere using compute_alpha_hdr()
  # Interpretation: HDR > 95 means (0,0) is excluded from the 95% HDR

  # ---------------------------
  # Method B: Bivariate normal DENSITY method
  # ---------------------------
  # Fit bivariate normal to (C, S) draws, evaluate density at (0,0)
  # Compare to density at the mode (posterior mean)
  #
  # This is DIFFERENT from the Wald test (Method E):
  # - Method B: How likely is (0,0) as a point under the fitted distribution?
  # - Method E: How far is the distribution center from (0,0)?
  #
  # We compute:
  # 1. Density at origin: f(0,0)
  # 2. Density at mode: f(mean_C, mean_S)
  # 3. Density ratio: f(0,0) / f(mode) - low ratio means origin is unlikely
  # 4. Percentile: what % of draws have lower density than origin?

  mean_CS <- c(mean(C), mean(S))
  cov_CS <- cov(cbind(C, S))

  method_B_results <- tryCatch({
    # Bivariate normal density function
    # f(x) = (2π)^(-1) |Σ|^(-1/2) exp(-0.5 * (x-μ)' Σ^(-1) (x-μ))

    cov_inv <- solve(cov_CS)
    cov_det <- det(cov_CS)
    norm_const <- 1 / (2 * pi * sqrt(cov_det))

    # Mahalanobis distance squared from origin to mean
    mahal_origin_sq <- as.numeric(t(mean_CS) %*% cov_inv %*% mean_CS)

    # Density at origin (0,0)
    dens_origin <- norm_const * exp(-0.5 * mahal_origin_sq)

    # Density at mode (mean) - this is the maximum density
    dens_mode <- norm_const  # exp(-0.5 * 0) = 1

    # Density ratio: how much lower is origin density compared to mode?
    dens_ratio <- dens_origin / dens_mode  # = exp(-0.5 * mahal_origin_sq)

    # What percentile of the distribution has density <= density at origin?
    # For bivariate normal, points with density <= d form an ellipse
    # The probability inside is 1 - exp(-0.5 * mahal_sq) where mahal_sq gives density d
    # So probability of density <= dens_origin is exp(-0.5 * mahal_origin_sq)
    # And probability of density > dens_origin is 1 - exp(-0.5 * mahal_origin_sq)
    prob_higher_dens <- 1 - exp(-0.5 * mahal_origin_sq)

    list(
      dens_origin = dens_origin,
      dens_mode = dens_mode,
      dens_ratio = dens_ratio,
      prob_higher_dens = prob_higher_dens,  # Like HDR level for bivariate normal
      mahal_sq = mahal_origin_sq
    )
  }, error = function(e) {
    list(dens_origin = NA, dens_mode = NA, dens_ratio = NA,
         prob_higher_dens = NA, mahal_sq = NA)
  })

  # ---------------------------
  # Method C: Quadrant method (CORRECTED)
  # ---------------------------
  # As per reviewer's suggestion: compute P(draws in origin's quadrant)
  # The origin (0,0) is in the quadrant OPPOSITE to where the median lies
  #
  # If median is in quadrant I (+,+), origin is in quadrant III (-,-)
  # If median is in quadrant II (-,+), origin is in quadrant IV (+,-)
  # If median is in quadrant III (-,-), origin is in quadrant I (+,+)
  # If median is in quadrant IV (+,-), origin is in quadrant II (-,+)
  #
  # LOW probability = origin quadrant is unlikely = SIGNIFICANT
  # This is the "smallest quadrant" the reviewer mentioned
  #
  # We also compute the probability in the median's quadrant for comparison

  median_C <- median(C)
  median_S <- median(S)

  # Probability of draws in SAME quadrant as median (for comparison)
  method_C_same_quadrant_prob <- mean(
    (sign(C) == sign(median_C)) & (sign(S) == sign(median_S))
  )

  # Probability of draws in ORIGIN's quadrant (opposite to median)
  # This is the reviewer's suggested method - LOW = significant
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

  # Also compute: probability that draws cross zero on EITHER axis
  # (i.e., draws in adjacent quadrants to median)
  method_C_cross_C_axis_prob <- mean(sign(C) != sign(median_C))
  method_C_cross_S_axis_prob <- mean(sign(S) != sign(median_S))

  # ---------------------------
  # Method D: Univariate CI method
  # ---------------------------
  # Test H0: C=0 and H0: S=0 separately using credible intervals
  # Reject overall null if EITHER C or S CI excludes zero (union test)
  # Multiple alpha levels: 90%, 95%, 99%

  C_ci_90 <- quantile(C, c(0.05, 0.95))
  C_ci_95 <- quantile(C, c(0.025, 0.975))
  C_ci_99 <- quantile(C, c(0.005, 0.995))

  S_ci_90 <- quantile(S, c(0.05, 0.95))
  S_ci_95 <- quantile(S, c(0.025, 0.975))
  S_ci_99 <- quantile(S, c(0.005, 0.995))

  C_excl_zero_90 <- (C_ci_90[1] > 0) | (C_ci_90[2] < 0)
  C_excl_zero_95 <- (C_ci_95[1] > 0) | (C_ci_95[2] < 0)
  C_excl_zero_99 <- (C_ci_99[1] > 0) | (C_ci_99[2] < 0)

  S_excl_zero_90 <- (S_ci_90[1] > 0) | (S_ci_90[2] < 0)
  S_excl_zero_95 <- (S_ci_95[1] > 0) | (S_ci_95[2] < 0)
  S_excl_zero_99 <- (S_ci_99[1] > 0) | (S_ci_99[2] < 0)

  # ---------------------------
  # Method E: Joint Wald test = Likelihood Ratio Test (for bivariate normal)
  # ---------------------------
  # H0: C = 0 AND S = 0 (joint hypothesis)
  #
  # WALD TEST:
  #   W = (θ̂ - θ₀)' [Var(θ̂)]^(-1) (θ̂ - θ₀)
  #   For θ = (C, S) and θ₀ = (0, 0): W = Mahalanobis²
  #   Under H0, W ~ χ²(2)
  #
  # LIKELIHOOD RATIO TEST:
  #   LRT = -2 ln[L(θ₀)/L(θ̂)] = -2 ln[f(0,0)/f(μ)]
  #   For bivariate normal: LRT = -2 × (-0.5 × Mahal²) = Mahal²
  #
  # KEY INSIGHT: For bivariate normal, Wald = LRT = Mahalanobis²
  # Both reference χ²(2), giving identical p-values.
  #
  # Method B (prob_higher_dens) is related: p_B = 1 - exp(-0.5 × Mahal²)
  # which equals pchisq(Mahal², df=2), so p_E = 1 - p_B
  #
  # The p-values will be related: p_E = 1 - (1 - density_ratio_B)
  # But we keep both for conceptual clarity and to match reviewer's list

  wald_results <- tryCatch({
    cov_inv <- solve(cov_CS)
    wald_stat <- as.numeric(t(mean_CS) %*% cov_inv %*% mean_CS)
    list(
      stat = wald_stat,
      pval = 1 - pchisq(wald_stat, df = 2)
    )
  }, error = function(e) list(stat = NA, pval = NA))

  method_E_wald_stat <- wald_results$stat
  method_E_pval <- wald_results$pval

  # ---------------------------
  # Method F: Amplitude-based tests (WITH PROPER RAYLEIGH QUANTILES)
  # ---------------------------
  # The reviewer suggested: Z = A / SE(A) > 1.96
  #
  # THEORETICAL PROBLEM:
  # Under H0 (C=0, S=0), if C ~ N(0, σ²) and S ~ N(0, σ²) independently,
  # then A = sqrt(C² + S²) follows a RAYLEIGH distribution with:
  #   E[A] = σ · sqrt(π/2) ≈ 1.253σ  (NOT zero!)
  #   SD[A] = σ · sqrt(2 - π/2) ≈ 0.655σ
  #   Mode[A] = σ
  #   CDF: F(a) = 1 - exp(-a²/(2σ²))
  #
  # So A/SD(A) does NOT have mean 0 under the null, and 1.96 is not
  # the correct critical value.
  #
  # We compute several versions:
  # F1. Original Z-score as reviewer suggested (for comparison)
  # F2. Rayleigh p-value using actual CDF (CORRECT METHOD)
  # F3. Rayleigh critical value comparison
  # F4. Coefficient of variation

  A_mean <- mean(A)
  A_se <- sd(A)

  # F1: Original Z-score (has theoretical issues but included per reviewer)
  method_F_zscore <- A_mean / A_se

  # Estimate σ for Rayleigh null distribution from posterior spread of C and S
  sigma_C <- sd(C)
  sigma_S <- sd(S)
  sigma_CS_pooled <- sqrt((sigma_C^2 + sigma_S^2) / 2)

  # F2: Rayleigh p-value using actual CDF
  # Under H0, A ~ Rayleigh(σ), so P(A > a) = exp(-a²/(2σ²))
  # p-value = P(A ≥ observed | H0) = exp(-A_mean²/(2σ²))
  method_F_rayleigh_pval <- exp(-A_mean^2 / (2 * sigma_CS_pooled^2))

  # F3: Compare to Rayleigh critical values
  # Critical value at α: a_crit = σ × sqrt(-2 × ln(α))
  # For α = 0.05: a_crit = σ × sqrt(-2 × ln(0.05)) = σ × 2.448
  # For α = 0.01: a_crit = σ × sqrt(-2 × ln(0.01)) = σ × 3.035
  a_crit_90 <- sigma_CS_pooled * sqrt(-2 * log(0.10))  # 2.146σ
  a_crit_95 <- sigma_CS_pooled * sqrt(-2 * log(0.05))  # 2.448σ
  a_crit_99 <- sigma_CS_pooled * sqrt(-2 * log(0.01))  # 3.035σ

  method_F_rayleigh_exceeds_90 <- A_mean > a_crit_90
  method_F_rayleigh_exceeds_95 <- A_mean > a_crit_95
  method_F_rayleigh_exceeds_99 <- A_mean > a_crit_99

  # Also compute the old "z-score" style statistics for comparison
  A_expected_null <- sigma_CS_pooled * sqrt(pi / 2)
  A_sd_null <- sigma_CS_pooled * sqrt(2 - pi / 2)
  method_F_rayleigh_z <- (A_mean - A_expected_null) / A_sd_null

  # F4: Coefficient of variation (CV) - lower CV = more precise estimate
  method_F_cv <- A_se / A_mean

  # Amplitude credible intervals at multiple levels
  A_ci_90 <- quantile(A, c(0.05, 0.95))
  A_ci_95 <- quantile(A, c(0.025, 0.975))
  A_ci_99 <- quantile(A, c(0.005, 0.995))

  # Ratio of CI lower bound to mean (precision indicator)
  method_F_ci_lower_over_mean <- A_ci_95[1] / A_mean

  # ---------------------------
  # R-squared: Variance explained by cosinor
  # ---------------------------
  # R² = Var(signal) / Var(total) = (A²/2) / (A²/2 + σ²)
  #
  # This is PRACTICAL significance, not statistical significance.
  # Reviewer's guidelines:
  #   R² < 0.01: Don't bother with cycles
  #   R² > 0.10: Must include cycles
  #   0.01 < R² < 0.10: Check statistical significance

  R2_draws <- (A^2 / 2) / (A^2 / 2 + sigma^2)
  R2_mean <- mean(R2_draws)
  R2_median <- median(R2_draws)
  R2_ci <- quantile(R2_draws, c(0.025, 0.975))

  # ---------------------------
  # Return all results (raw values, not significance flags)
  # ---------------------------
  tibble(
    # Method B: Bivariate normal density
    method_B_dens_origin = method_B_results$dens_origin,
    method_B_dens_mode = method_B_results$dens_mode,
    method_B_dens_ratio = method_B_results$dens_ratio,
    method_B_prob_higher_dens = method_B_results$prob_higher_dens,  # Like HDR level
    method_B_mahal_sq = method_B_results$mahal_sq,

    # Method C: Quadrant (CORRECTED - origin quadrant prob, LOW = significant)
    method_C_origin_quadrant_prob = method_C_origin_quadrant_prob,
    method_C_same_quadrant_prob = method_C_same_quadrant_prob,  # For comparison
    method_C_cross_C_axis_prob = method_C_cross_C_axis_prob,
    method_C_cross_S_axis_prob = method_C_cross_S_axis_prob,

    # Method D: Univariate CIs
    method_D_C_ci95_lower = C_ci_95[1],
    method_D_C_ci95_upper = C_ci_95[2],
    method_D_S_ci95_lower = S_ci_95[1],
    method_D_S_ci95_upper = S_ci_95[2],
    method_D_C_excl_zero_90 = C_excl_zero_90,
    method_D_C_excl_zero_95 = C_excl_zero_95,
    method_D_C_excl_zero_99 = C_excl_zero_99,
    method_D_S_excl_zero_90 = S_excl_zero_90,
    method_D_S_excl_zero_95 = S_excl_zero_95,
    method_D_S_excl_zero_99 = S_excl_zero_99,

    # Method E: Wald test (joint test, different interpretation from B)
    method_E_wald_stat = method_E_wald_stat,
    method_E_pval = method_E_pval,

    # Method F: Amplitude-based tests
    method_F_amp_mean = A_mean,
    method_F_amp_se = A_se,
    method_F_zscore = method_F_zscore,              # Original (has issues)
    method_F_rayleigh_pval = method_F_rayleigh_pval,  # Proper Rayleigh p-value
    method_F_rayleigh_exceeds_90 = method_F_rayleigh_exceeds_90,
    method_F_rayleigh_exceeds_95 = method_F_rayleigh_exceeds_95,
    method_F_rayleigh_exceeds_99 = method_F_rayleigh_exceeds_99,
    method_F_rayleigh_z = method_F_rayleigh_z,      # Z-score style (for comparison)
    method_F_cv = method_F_cv,                       # Coefficient of variation
    method_F_ci_lower_ratio = method_F_ci_lower_over_mean,

    # Amplitude CIs
    amp_ci90_lower = A_ci_90[1],
    amp_ci90_upper = A_ci_90[2],
    amp_ci95_lower = A_ci_95[1],
    amp_ci95_upper = A_ci_95[2],
    amp_ci99_lower = A_ci_99[1],
    amp_ci99_upper = A_ci_99[2],

    # R-squared (variance explained by cycle)
    R2_mean = R2_mean,
    R2_median = R2_median,
    R2_ci_lower = R2_ci[1],
    R2_ci_upper = R2_ci[2]
  )
}

# Apply to all individuals
cat("  Computing for each individual...\n")
alternative_tests <- d_draws %>%
  group_by(id) %>%
  group_modify(~ compute_alternative_tests(.x)) %>%
  ungroup()

cat("  alternative_tests dimensions:", nrow(alternative_tests), "x", ncol(alternative_tests), "\n")

# IMPORTANT: Start with the BASE ests_level1 (without alternative tests)
# to prevent column duplication from repeated runs
ests_level1_base <- ests_level1 %>%
  select(id:amp_hdr)  # Keep only base columns (up to amp_hdr)

# Merge with alternative test results
ests_level1 <- ests_level1_base %>%
  left_join(alternative_tests, by = "id")

cat("  After join, ests_level1 has", ncol(ests_level1), "columns\n")

# Save the FULL estimates (base + alternative tests)
ests_level1 %>%
  mutate(item = item_, .before = 1) %>%
  saveRDS(here("fits",
               paste0("ests_level1_", item_, "_n1_strategy1.rds")))

# ---------------------------
# Summary of method agreement across multiple thresholds
# ---------------------------

cat("\n=== Comparison of significance testing methods ===\n")
cat("(Testing H0: A_i = 0, i.e., no cyclic trend)\n\n")

# Create significance flags at multiple thresholds
# NOTE on interpretation:
#   - Method A (HDR): HIGH value = significant (origin outside HDR)
#   - Method B (BivNorm density): HIGH prob_higher_dens = significant (like HDR)
#   - Method C (Quadrant): LOW origin_quadrant_prob = significant (origin unlikely)
#   - Method D (Univariate CI): CI excludes 0 = significant
#   - Method E (Wald): LOW p-value = significant
#   - Method F (Z-score): HIGH Z = significant (but see caveats)

method_comparison <- ests_level1 %>%
  transmute(
    id = id,

    # Method A: HDR at multiple thresholds
    # HIGH HDR level = significant (origin excluded from that % region)
    A_HDR_80 = amp_hdr > 80,
    A_HDR_85 = amp_hdr > 85,
    A_HDR_90 = amp_hdr > 90,
    A_HDR_95 = amp_hdr > 95,

    # Method B: Bivariate normal density
    # HIGH prob_higher_dens = significant (origin in low-density region)
    # This is analogous to HDR but assuming bivariate normality
    B_BivNorm_80 = method_B_prob_higher_dens > 0.80,
    B_BivNorm_90 = method_B_prob_higher_dens > 0.90,
    B_BivNorm_95 = method_B_prob_higher_dens > 0.95,

    # Method C: Quadrant method (CORRECTED)
    # LOW origin_quadrant_prob = significant (few draws in origin's quadrant)
    # Thresholds: < 0.10, < 0.05, < 0.025 (one-sided)
    C_Quad_10 = method_C_origin_quadrant_prob < 0.10,
    C_Quad_05 = method_C_origin_quadrant_prob < 0.05,
    C_Quad_025 = method_C_origin_quadrant_prob < 0.025,

    # Method D: Univariate CIs (C or S excludes 0)
    # Either CI excludes 0 = significant
    D_CorS_90 = method_D_C_excl_zero_90 | method_D_S_excl_zero_90,
    D_CorS_95 = method_D_C_excl_zero_95 | method_D_S_excl_zero_95,
    D_CorS_99 = method_D_C_excl_zero_99 | method_D_S_excl_zero_99,

    # Method D variant: Both C and S exclude 0 (more stringent)
    D_CandS_90 = method_D_C_excl_zero_90 & method_D_S_excl_zero_90,
    D_CandS_95 = method_D_C_excl_zero_95 & method_D_S_excl_zero_95,
    D_CandS_99 = method_D_C_excl_zero_99 & method_D_S_excl_zero_99,

    # Method E: Wald test at multiple alpha levels
    # LOW p-value = significant
    E_Wald_10 = method_E_pval < 0.10,
    E_Wald_05 = method_E_pval < 0.05,
    E_Wald_01 = method_E_pval < 0.01,

    # Method F: Amplitude-based tests
    # F1: Original Z-score (has theoretical issues - see documentation)
    F_Z_165 = method_F_zscore > 1.645,  # ~90% one-sided
    F_Z_196 = method_F_zscore > 1.96,   # ~95% one-sided
    F_Z_258 = method_F_zscore > 2.576,  # ~99% one-sided

    # F2: Rayleigh test with PROPER quantiles (A > σ√(-2ln(α)))
    # This is the CORRECT test using actual Rayleigh distribution
    F_Rayleigh_90 = method_F_rayleigh_exceeds_90,
    F_Rayleigh_95 = method_F_rayleigh_exceeds_95,
    F_Rayleigh_99 = method_F_rayleigh_exceeds_99,

    # F3: Rayleigh p-value thresholds (alternative to critical values)
    F_Rayleigh_pval_10 = method_F_rayleigh_pval < 0.10,
    F_Rayleigh_pval_05 = method_F_rayleigh_pval < 0.05,
    F_Rayleigh_pval_01 = method_F_rayleigh_pval < 0.01,

    # R² thresholds (practical significance, per reviewer's guidelines)
    R2_gt_01 = R2_median > 0.01,  # Minimum to consider
    R2_gt_05 = R2_median > 0.05,  # Moderate effect
    R2_gt_10 = R2_median > 0.10   # Must include cycles
  )

# Summary table
method_summary <- method_comparison %>%
  select(-id) %>%
  summarise(across(everything(), ~ sum(.x, na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "Method", values_to = "n_significant") %>%
  mutate(
    pct_significant = round(n_significant / nrow(method_comparison) * 100, 1),
    Method = factor(Method, levels = Method)
  )

cat("Number and percentage of individuals meeting each criterion:\n\n")
print(method_summary, n = 40)

# ---------------------------
# Agreement matrix for 95%-level methods
# ---------------------------

cat("\n\nPairwise agreement between methods at ~95% level (% agreeing):\n")
cat("Note: For Method C, we use < 0.05 threshold (LOW = significant)\n\n")

sig_cols_95 <- method_comparison %>%
  select(A_HDR_95, B_BivNorm_95, C_Quad_05, D_CorS_95, E_Wald_05, F_Rayleigh_95, F_Rayleigh_pval_05)

agreement_matrix <- matrix(NA, ncol(sig_cols_95), ncol(sig_cols_95))
colnames(agreement_matrix) <- names(sig_cols_95)
rownames(agreement_matrix) <- names(sig_cols_95)

for (i in 1:ncol(sig_cols_95)) {
  for (j in 1:ncol(sig_cols_95)) {
    agreement_matrix[i, j] <- mean(sig_cols_95[[i]] == sig_cols_95[[j]], na.rm = TRUE) * 100
  }
}

print(round(agreement_matrix, 1))

# ---------------------------
# Method-specific notes
# ---------------------------

cat("\n\n=== Method Notes ===\n")
cat("A (HDR): Non-parametric, uses kernel density estimation\n")
cat("B (BivNorm): Parametric (assumes bivariate normal), density-based = LRT\n")
cat("C (Quadrant): Non-parametric, counts draws in origin's quadrant\n")
cat("D (Univariate): Marginal tests, union (C OR S) or intersection (C AND S)\n")
cat("E (Wald): Joint test = LRT for bivariate normal; Mahalanobis² ~ χ²(2)\n")
cat("F (Rayleigh): CORRECT amplitude test using Rayleigh quantiles\n")
cat("   Critical value at α: a_crit = σ × sqrt(-2 × ln(α))\n")
cat("   p-value: P(A ≥ observed | H0) = exp(-A²/(2σ²))\n")
cat("R²: Practical significance (reviewer: <1% ignore, >10% must include)\n")

# ---------------------------
# R² summary (practical significance)
# ---------------------------

cat("\n\n=== R² Summary (Practical Significance) ===\n")
cat("Variance explained by cosinor cycle:\n\n")

cat("R² distribution:\n")
print(summary(ests_level1$R2_median))

cat("\nIndividuals by R² threshold:\n")
cat("  R² > 1%: ", sum(ests_level1$R2_median > 0.01),
    " (", round(mean(ests_level1$R2_median > 0.01) * 100, 1), "%)\n", sep = "")
cat("  R² > 5%: ", sum(ests_level1$R2_median > 0.05),
    " (", round(mean(ests_level1$R2_median > 0.05) * 100, 1), "%)\n", sep = "")
cat("  R² > 10%: ", sum(ests_level1$R2_median > 0.10),
    " (", round(mean(ests_level1$R2_median > 0.10) * 100, 1), "%)\n", sep = "")

# ---------------------------
# Cross-tabulation: Statistical vs Practical significance
# ---------------------------

cat("\n\n=== Statistical vs Practical Significance ===\n")
cat("(HDR > 95% vs R² > 5%)\n\n")

crosstab <- table(
  `HDR > 95%` = method_comparison$A_HDR_95,
  `R² > 5%` = method_comparison$R2_gt_05
)
print(crosstab)

cat("\nInterpretation:\n")
cat("- TRUE/TRUE: Statistically AND practically significant cycles\n")
cat("- TRUE/FALSE: Statistically significant but tiny effect (R² < 5%)\n")
cat("- FALSE/TRUE: Large effect but uncertain (wide posterior)\n")
cat("- FALSE/FALSE: No evidence of meaningful cycles\n")

# -----------------------------------------------------------------------------
# 9. Summary
# -----------------------------------------------------------------------------

cat("\n=== Summary ===\n")
cat("Individuals fitted:", length(ids_to_fit), "\n")
cat("Posterior samples per person:", CHAINS * (ITER - WARMUP), "\n")

cat("\nAmplitude (median):\n")
print(summary(ests_level1$amp_median))

cat("\nAmplitude HDR levels:\n")
cat("  > 95%:", mean(ests_level1$amp_hdr > 95) * 100, "%\n")
cat("  > 90%:", mean(ests_level1$amp_hdr > 90) * 100, "%\n")

cat("\nPhase CI width (hours):\n")
print(summary(ests_level1$phi_circ_ci_width * 12 / pi))

cat("\nDone!\n")
beepr::beep(5)
