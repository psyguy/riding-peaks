# =============================================================================
# Replicated N=1 Cosinor Model using Strategy 2 (No Pooling)
# =============================================================================
#
# This script fits a cosinor model where each individual gets their own
# parameters (M_i, C_i, S_i, sigma_i) estimated as fixed effects.
# No information is shared between individuals (no pooling/shrinkage).
#
# Model specification:
#   y_ij ~ Normal(M_i + C_i * cos(2*pi*t/24) + S_i * sin(2*pi*t/24), sigma_i)
#
# This is achieved in brms using factor(id) for all parameters, including
# the residual standard deviation (sigma).
# =============================================================================

# -----------------------------------------------------------------------------
# 1. Load required packages
# -----------------------------------------------------------------------------

library(brms)
library(tidyverse)
library(here)
library(circular)  # For circular statistics on phase

# -----------------------------------------------------------------------------
# 2. Load and prepare data
# -----------------------------------------------------------------------------

item_ <- "pa"
d <- readRDS(here("data",
                  "d_leuven_joined.rds")) %>%
  filter(item == item_)

# Quick check of the data structure
# Expects columns: id, y, t, co, si (and possibly others)
glimpse(d)

# Get unique IDs for later use
ids <- unique(d$id)
n_ids <- length(ids)
cat("Number of individuals:", n_ids, "\n")

# -----------------------------------------------------------------------------
# 3. Fit the no-pooling cosinor model (Strategy 2)
# -----------------------------------------------------------------------------
#
# Key differences from the multilevel model:
# - factor(id) creates separate fixed effects for each person
# - No random effects structure, so no population-level covariance estimated
# - sigma ~ 0 + factor(id) gives each person their own residual SD
# - Each person's estimates come purely from their own data
# -----------------------------------------------------------------------------

cat("Fitting no-pooling cosinor model...\n")
(st <- Sys.time())

m_n1 <- brm(
  brmsformula(
    # Mean structure: person-specific intercept and slopes
    # 0 + factor(id) removes global intercept, giving each person their own M_i
    # factor(id):co and factor(id):si give person-specific C_i and S_i
    y ~ 0 + factor(id) + factor(id):co + factor(id):si,
    # Residual SD: person-specific sigma_i (on log scale internally)
    sigma ~ 0 + factor(id)
  ),
  data = d,
  backend = "cmdstanr",
  chains = 4,
  iter = 3000,
  warmup = 1000,
  threads = threading(2),
  cores = 8
)

cat("Model fitting time:", Sys.time() - st, "\n")

# Optional: beep when done
beepr::beep(1)

# Save the fitted model
saveRDS(m_n1,
        here("fits",
             paste0("brms_", item_, "_n1_strategy2.rds")))

# -----------------------------------------------------------------------------
# 4. Extract posterior samples
# -----------------------------------------------------------------------------
#
# The parameter names in Strategy 2 differ from the multilevel model:
# - b_factoridX = M_i (intercept for person X)
# - b_factoridX:co = C_i (cosine coefficient for person X)
# - b_factoridX:si = S_i (sine coefficient for person X)
# - b_sigma_factoridX = log(sigma_i) for person X
# -----------------------------------------------------------------------------

# Load model if not in memory
# m_n1 <- readRDS(here("fits", paste0("brms_", item_, "_n1_strategy2.rds")))

cat("Extracting posterior samples...\n")

# Get all draws as a data frame
draws_raw <- as_draws_df(m_n1)

# Extract and reshape the draws
# This is more complex than the multilevel model because parameter names
# encode the person ID differently

d_draws <- draws_raw %>%
  mutate(iteration = row_number()) %>%
  # Select only the coefficient columns (start with "b_")
  select(iteration, starts_with("b_")) %>%
  # Pivot to long format for easier manipulation

  pivot_longer(
    cols = -iteration,
    names_to = "param",
    values_to = "value"
  ) %>%
  # Parse parameter names to extract id and parameter type
  mutate(
    # Check if it's a sigma parameter
    is_sigma = str_detect(param, "^b_sigma_"),
    # Extract the id number from the parameter name
    id = case_when(
      # Sigma parameters: b_sigma_factoridX
      is_sigma ~ str_extract(param, "(?<=b_sigma_factorid)\\d+"),
      # Interaction terms: b_factoridX:co or b_factoridX:si
      str_detect(param, ":") ~ str_extract(param, "(?<=b_factorid)\\d+(?=:)"),
      # Main effects (intercepts): b_factoridX
      TRUE ~ str_extract(param, "(?<=b_factorid)\\d+")
    ),
    # Determine parameter type
    param_type = case_when(
      is_sigma ~ "logsigma",
      str_detect(param, ":co$") ~ "co",
      str_detect(param, ":si$") ~ "si",
      TRUE ~ "mesor"
    )
  ) %>%
  # Remove the helper column

  select(-is_sigma, -param) %>%
  # Pivot back to wide format with one column per parameter type
  pivot_wider(
    names_from = param_type,
    values_from = value
  ) %>%
  # Convert id to integer for consistency
  mutate(id = as.integer(id)) %>%
  # Compute derived quantities
  mutate(
    # Residual SD and variance
    sigma = exp(logsigma),
    sigma2 = sigma^2,
    # Amplitude: distance from origin in (C, S) space
    amp = sqrt(co^2 + si^2),
    # Phase: angle in (C, S) space
    # atan2 handles all quadrants correctly
    phi = atan2(si, co) %% (2 * pi),
    # Alternative phase calculations with different cycle start times
    # phi.begins6: cycle starts at 6:00 instead of midnight
    phi.begins6 = (phi - pi/2) %% (2 * pi),
    # phi.begins12: cycle starts at 12:00 (noon)
    phi.begins12 = (phi - pi) %% (2 * pi)
  )

# Save the draws
d_draws %>%
  mutate(item = item_, .before = 1) %>%
  saveRDS(here("fits",
               paste0("draws_", item_, "_n1_strategy2.rds")))

cat("Posterior samples extracted and saved.\n")

# -----------------------------------------------------------------------------
# 5. Compute Level-1 (person-specific) estimates
# -----------------------------------------------------------------------------
#
# For each person, we compute:
# - Point estimates (mean, median)
# - 95% credible interval width
# - For phase (phi): both linear and circular statistics
# - HDR level for amplitude (where does [0,0] fall in the C,S posterior?)
# -----------------------------------------------------------------------------

# Load required packages for HDR computation
require(ks)

# Function to compute alpha-HDR level for amplitude
# Returns the HDR level at which the origin [0,0] falls
compute_alpha_hdr <- function(d,
                              x_0 = 0,
                              y_0 = 0,
                              n_grid = 200) {
  # 1. Kernel density estimation on (C, S) samples
  fhat <- kde(d, gridsize = c(n_grid, n_grid))

  xseq <- fhat$eval.points[[1]]
  yseq <- fhat$eval.points[[2]]
  zmat <- fhat$estimate

  dx <- diff(xseq)[1]
  dy <- diff(yseq)[1]
  cell_area <- dx * dy

  # 2. Compute cumulative probability for HDR thresholds
  dens_vec <- as.vector(zmat)
  dens_sorted <- sort(dens_vec, decreasing = TRUE)
  cumulative_mass <- cumsum(dens_sorted * cell_area)
  cumulative_mass <- cumulative_mass / max(cumulative_mass)

  # 3. Evaluate density at the origin
  f_point <- predict(fhat, x = cbind(x_0, y_0))

  # Find where the origin's density fits in the distribution
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

# Helper function for linear statistics
compute_linear_stats <- function(data, variable) {
  tibble(
    mean = mean(data[[variable]], na.rm = TRUE),
    median = median(data[[variable]], na.rm = TRUE),
    ci_width = quantile(data[[variable]], 0.975, na.rm = TRUE) -
               quantile(data[[variable]], 0.025, na.rm = TRUE)
  ) %>%
    rename_with(~ paste0(variable, "_", .))
}

# Helper function for circular statistics (for phase variables)
compute_circular_stats <- function(data, variable) {
  ci_2.5 <- quantile.circular(circular(data[[variable]]), 0.025)
  ci_97.5 <- quantile.circular(circular(data[[variable]]), 0.975)


  tibble(
    circ_mean = (mean.circular(circular(data[[variable]])) %% (2 * pi)) %>% as.numeric(),
    circ_median = (median.circular(circular(data[[variable]])) %% (2 * pi)) %>% as.numeric(),
    circ_ci_2.5 = (ci_2.5 %% (2 * pi)) %>% as.numeric(),
    circ_ci_97.5 = (ci_97.5 %% (2 * pi)) %>% as.numeric(),
    circ_ci_width = ((ci_97.5 - ci_2.5) %% (2 * pi)) %>% as.numeric()
  ) %>%
    rename_with(~ paste0(variable, "_", .))
}

cat("Computing level-1 estimates...\n")

# Compute person-specific summary statistics
ests_level1 <- d_draws %>%
  group_by(id) %>%
  summarise(
    bind_cols(
      # Linear parameters: mean, median, CI width
      compute_linear_stats(cur_data(), "mesor"),
      compute_linear_stats(cur_data(), "co"),
      compute_linear_stats(cur_data(), "si"),
      compute_linear_stats(cur_data(), "amp"),
      compute_linear_stats(cur_data(), "logsigma"),
      compute_linear_stats(cur_data(), "sigma"),
      compute_linear_stats(cur_data(), "sigma2"),
      # Phase with linear methods (can be misleading for circular data!)
      compute_linear_stats(cur_data(), "phi"),
      compute_linear_stats(cur_data(), "phi.begins6"),
      compute_linear_stats(cur_data(), "phi.begins12"),
      # Phase with circular methods (preferred for circular data)
      compute_circular_stats(cur_data(), "phi"),
      compute_circular_stats(cur_data(), "phi.begins6"),
      compute_circular_stats(cur_data(), "phi.begins12")
    ),
    # HDR level: what % HDR contains the origin in (C, S) space?
    # High values (close to 100) suggest amplitude is "significantly" > 0
    amp_hdr = compute_alpha_hdr(cbind(cur_data()$co, cur_data()$si))
  ) %>%
  relocate(amp_hdr, .before = logsigma_mean)

# Save level-1 estimates
ests_level1 %>%
  mutate(item = item_, .before = 1) %>%
  saveRDS(here("fits",
               paste0("ests_level1_", item_, "_n1_strategy2.rds")))

cat("Level-1 estimates computed and saved.\n")

# -----------------------------------------------------------------------------
# 6. Quick summary of results
# -----------------------------------------------------------------------------

cat("\n=== Summary of Level-1 Estimates ===\n")
cat("Number of individuals:", nrow(ests_level1), "\n\n")

cat("Amplitude statistics:\n")
print(summary(ests_level1$amp_median))

cat("\nAmplitude HDR levels (% with origin outside X% HDR):\n")
cat("  > 95%:", mean(ests_level1$amp_hdr > 95) * 100, "%\n")
cat("  > 90%:", mean(ests_level1$amp_hdr > 90) * 100, "%\n")
cat("  > 80%:", mean(ests_level1$amp_hdr > 80) * 100, "%\n")

cat("\nPhase (circular median) statistics (in hours, 0 = midnight):\n")
phi_hours <- ests_level1$phi_circ_median * 12 / pi
print(summary(phi_hours))

cat("\n95% CI width for phase (in hours):\n")
ci_width_hours <- ests_level1$phi_circ_ci_width * 12 / pi
print(summary(ci_width_hours))

cat("\nDone!\n")
