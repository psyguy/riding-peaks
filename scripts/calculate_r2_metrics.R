# =============================================================================
# Calculate and Compare R-squared Metrics for Multilevel Cosinor Model
# =============================================================================
#
# This script calculates and compares three different approaches to R-squared
# for a multilevel cosinor model, as discussed in response to the reviewer's
# comments.
#
# 1. Distribution of Individual R²s:
#    - R² calculated for each person, based on their specific parameters.
#    - The "population" R² is the distribution of these individual values.
#
# 2. Simplified Marginal & Conditional R²:
#    - Based on Nakagawa & Schielzeth's framework.
#    - Assumes a full 24-hour cycle, causing random effect covariance terms to be zero.
#    - This is a common approach but is an approximation for gapped data.
#
# =============================================================================
# Calculate and Compare R-squared Metrics for Multilevel Cosinor Model
# =============================================================================
#
# This script calculates and compares three different approaches to R-squared
# for a multilevel cosinor model, as discussed in response to the reviewer's
# comments. This version uses explicit `for` loops for clarity and robustness
# to avoid issues with vectorized operations that can lead to NaNs.
#
# 1. Distribution of Individual R²s:
#    - R² calculated for each person, based on their specific parameters.
#    - The "population" R² is the distribution of these individual values.
#
# 2. Simplified Marginal & Conditional R²:
#    - Based on Nakagawa & Schielzeth's framework.
#    - Assumes a full 24-hour cycle, causing random effect covariance terms to be zero.
#    - This is a common approach but is an approximation for gapped data.
#
# 3. Exact Marginal & Conditional R²:
#    - A more accurate version for gapped data (e.g., waking hours only).
#    - Accounts for non-zero means of cos(t) and sin(t) over the observed
#      data window and correctly includes the random effect covariance terms.
#
# =============================================================================

# -----------------------------------------------------------------------------
# 1. Setup: Load packages and data
# -----------------------------------------------------------------------------

library(brms)
library(tidyverse)
library(here)

cat("Loading data...\n")

# Define the dataset to analyze
item_ <- "pa"

# Load the raw data (needed for exact R² calculation)
d <- readRDS(here("data", "d_leuven_joined.rds")) %>%
  filter(item == item_)

# Load the brms model fit
m_fit <- readRDS(here("fits", paste0("brms_", item_, "_random_var.rds")))

# Load the posterior draws (individual-level)
d_draws <- readRDS(here("fits", paste0("draws_", item_, "_random_var.rds")))

# Load the posterior draws (population-level fixed and random effects)
pop_draws <- as_draws_df(m_fit)

cat("Data loaded successfully.\n")

# -----------------------------------------------------------------------------
# 2. Method 1: Distribution of Individual R²s
# -----------------------------------------------------------------------------
# This calculates R² = Var(signal) / (Var(signal) + Var(residual)) for each
# person at each posterior draw. The "population" measure is the distribution
# of these individual R² values.

cat("Calculating Method 1: Distribution of Individual R²s...\n")

r2_individual_draws <- d_draws %>%
  mutate(
    R2_individual = (amp^2 / 2) / (amp^2 / 2 + sigma^2)
  ) %>%
  select(id, iteration, R2_individual)

# We can summarize this by taking the median R² for each person
r2_individual_summary <- r2_individual_draws %>%
  group_by(id) %>%
  summarise(R2_median = median(R2_individual))

# -----------------------------------------------------------------------------
# 3. Method 2: Simplified Marginal & Conditional R²
# -----------------------------------------------------------------------------
# This method assumes a full 24-hour cycle, where E[cos(t)]=0, E[sin(t)]=0,
# Var(cos(t))=1/2, and Var(sin(t))=1/2. This causes the random effect
# covariance terms to drop out of the Var_random calculation.

cat("Calculating Method 2: Simplified Marginal & Conditional R²...\n")

# Pre-allocate a list to store results for each MCMC draw
results_simplified <- vector("list", nrow(pop_draws))

# Get the column names of the individual random sigma effects
sigma_re_cols <- grep("r_id__sigma_Intercept", names(pop_draws), value = TRUE)

for (i in 1:nrow(pop_draws)) {
  # Extract the current row (MCMC draw)
  row <- pop_draws[i, ]
  
  # --- Calculate average residual variance for this draw ---
  # Extract the random sigma effects for all individuals for this one draw
  sigma_re_values <- as.numeric(row[sigma_re_cols])
  # Calculate sigma_i for each person
  sigma_i_values <- exp(row$b_sigma_Intercept + sigma_re_values)
  # The average residual variance is the mean of the squared sigmas
  var_residual <- mean(sigma_i_values^2)
  
  # --- Calculate other variance components ---
  var_fixed <- (row$b_co^2 + row$b_si^2) / 2
  var_random <- row$sd_id__Intercept^2 + (row$sd_id__co^2 / 2) + (row$sd_id__si^2 / 2)
  var_total <- var_fixed + var_random + var_residual
  
  # --- Store results ---
  results_simplified[[i]] <- tibble(
    R2_marginal_simplified = var_fixed / var_total,
    R2_conditional_simplified = (var_fixed + var_random) / var_total
  )
}

# Combine the list of tibbles into a single data frame
r2_simplified_draws <- bind_rows(results_simplified)

# -----------------------------------------------------------------------------
# 4. Method 3: Exact Marginal & Conditional R² for Gapped Data
# -----------------------------------------------------------------------------
# This method correctly accounts for the properties of the observed predictors
# (cos(t) and sin(t)) over the actual data window, which has gaps.

cat("Calculating Method 3: Exact Marginal & Conditional R²...\n")

# First, calculate the empirical properties of the predictors from the raw data
predictor_props <- d %>%
  summarise(
    var_U = var(co),
    var_V = var(si),
    cov_UV = cov(co, si),
    E_U = mean(co),
    E_V = mean(si),
    E_U2 = mean(co^2),
    E_V2 = mean(si^2),
    E_UV = mean(co * si)
  )

# Pre-allocate a list for the results
results_exact <- vector("list", nrow(pop_draws))

for (i in 1:nrow(pop_draws)) {
  # Extract the current row (MCMC draw)
  row <- pop_draws[i, ]
  
  # --- Calculate average residual variance (same as before) ---
  sigma_re_values <- as.numeric(row[sigma_re_cols])
  sigma_i_values <- exp(row$b_sigma_Intercept + sigma_re_values)
  var_residual <- mean(sigma_i_values^2)
  
  # --- Variance of fixed effects prediction ---
  var_fixed <- row$b_co^2 * predictor_props$var_U + row$b_si^2 * predictor_props$var_V + 2 * row$b_co * row$b_si * predictor_props$cov_UV
  
  # --- Variance of random effects (EXACT formula) ---
  var_m <- row$sd_id__Intercept^2
  var_c <- row$sd_id__co^2
  var_s <- row$sd_id__si^2
  cov_mc <- row$cor_id__Intercept__co * row$sd_id__Intercept * row$sd_id__co
  cov_ms <- row$cor_id__Intercept__si * row$sd_id__Intercept * row$sd_id__si
  cov_cs <- row$cor_id__co__si * row$sd_id__co * row$sd_id__si
  
  var_random <- var_m +
    var_c * predictor_props$E_U2 +
    var_s * predictor_props$E_V2 +
    2 * cov_mc * predictor_props$E_U +
    2 * cov_ms * predictor_props$E_V +
    2 * cov_cs * predictor_props$E_UV
    
  var_total <- var_fixed + var_random + var_residual
  
  # --- Store results ---
  results_exact[[i]] <- tibble(
    R2_marginal_exact = var_fixed / var_total,
    R2_conditional_exact = (var_fixed + var_random) / var_total
  )
}

# Combine the list of tibbles into a single data frame
r2_exact_draws <- bind_rows(results_exact)

# -----------------------------------------------------------------------------
# 5. Combine and Plot for Comparison
# -----------------------------------------------------------------------------

cat("Combining results and generating plot...\n")

# Combine all R² posterior distributions into a long-format data frame
comparison_data <- bind_cols(
  r2_simplified_draws,
  r2_exact_draws
) %>%
  pivot_longer(everything(), names_to = "Metric", values_to = "R2")

# Create the comparison plot
p_r2_comparison <- ggplot(comparison_data, aes(x = R2, fill = Metric)) +
  geom_density(alpha = 0.6, na.rm = TRUE) +
  scale_fill_viridis_d(option = "D", name = "R² Metric") +
  labs(
    title = "Comparison of Population-Level R² Metrics",
    subtitle = "Posterior distributions for different calculation methods",
    x = "R² Value",
    y = "Density"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

print(p_r2_comparison)

# Save the plot
ggsave(here("figs", "r2_metric_comparison.pdf"), p_r2_comparison, width = 8, height = 6, device = cairo_pdf)
ggsave(here("figs", "r2_metric_comparison.png"), p_r2_comparison, width = 8, height = 6, dpi = 150)

cat("Plot saved to figs/r2_metric_comparison.pdf\n")
cat("Done!\n")
#    - A more accurate version for gapped data (e.g., waking hours only).
#    - Accounts for non-zero means of cos(t) and sin(t) over the observed
#      data window and correctly includes the random effect covariance terms.
#
# =============================================================================

# -----------------------------------------------------------------------------
# 1. Setup: Load packages and data
# -----------------------------------------------------------------------------

library(brms)
library(tidyverse)
library(here)

cat("Loading data...\n")

# Define the dataset to analyze
item_ <- "pa"

# Load the raw data (needed for exact R² calculation)
d <- readRDS(here("data", "d_leuven_joined.rds")) %>%
  filter(item == item_)

# Load the brms model fit
m_fit <- readRDS(here("fits", paste0("brms_", item_, "_random_var.rds")))

# Load the posterior draws (individual-level)
d_draws <- readRDS(here("fits", paste0("draws_", item_, "_random_var.rds")))

# Load the posterior draws (population-level fixed and random effects)
pop_draws <- as_draws_df(m_fit)

cat("Data loaded successfully.\n")

# -----------------------------------------------------------------------------
# 2. Method 1: Distribution of Individual R²s
# -----------------------------------------------------------------------------
# This calculates R² = Var(signal) / (Var(signal) + Var(residual)) for each
# person at each posterior draw. The "population" measure is the distribution
# of these individual R² values.

cat("Calculating Method 1: Distribution of Individual R²s...\n")

r2_individual_draws <- d_draws %>%
  mutate(
    R2_individual = (amp^2 / 2) / (amp^2 / 2 + sigma^2)
  ) %>%
  select(id, iteration, R2_individual)

# We can summarize this by taking the median R² for each person
r2_individual_summary <- r2_individual_draws %>%
  group_by(id) %>%
  summarise(R2_median = median(R2_individual))

# -----------------------------------------------------------------------------
# 3. Method 2: Simplified Marginal & Conditional R²
# -----------------------------------------------------------------------------
# This method assumes a full 24-hour cycle, where E[cos(t)]=0, E[sin(t)]=0,
# Var(cos(t))=1/2, and Var(sin(t))=1/2. This causes the random effect
# covariance terms to drop out of the Var_random calculation.

cat("Calculating Method 2: Simplified Marginal & Conditional R²...\n")

# Helper function to process one row (one MCMC draw)
calculate_r2_simplified_row <- function(row) {
  # Get all individual random sigma effects for this draw
  sigma_re_cols <- grep("r_id__sigma_Intercept", names(row), value = TRUE)
  sigma_re_values <- as.numeric(row[sigma_re_cols])
  
  # Calculate average residual variance for this draw
  sigma_i_values <- exp(row$b_sigma_Intercept + sigma_re_values)
  var_residual <- mean(sigma_i_values^2)
  
  # Calculate other variance components
  var_fixed <- (row$b_co^2 + row$b_si^2) / 2
  var_random <- row$sd_id__Intercept^2 + (row$sd_id__co^2 / 2) + (row$sd_id__si^2 / 2)
  var_total <- var_fixed + var_random + var_residual
  
  # Return as a tibble
  tibble(
    R2_marginal_simplified = var_fixed / var_total,
    R2_conditional_simplified = (var_fixed + var_random) / var_total
  )
}

# Apply the function to each row of the posterior draws
# `purrr::map_dfr` iterates over rows and binds the results into a data frame
r2_simplified_draws <- purrr::map_dfr(
  1:nrow(pop_draws),
  ~ calculate_r2_simplified_row(pop_draws[.x, ])
)

# Check for issues
if(any(is.na(r2_simplified_draws))) {
  warning("NaNs detected in simplified R2 calculation!")
  print(summary(r2_simplified_draws))
}

# -----------------------------------------------------------------------------
# 4. Method 3: Exact Marginal & Conditional R² for Gapped Data
# -----------------------------------------------------------------------------
# This method correctly accounts for the properties of the observed predictors
# (cos(t) and sin(t)) over the actual data window, which has gaps.

cat("Calculating Method 3: Exact Marginal & Conditional R²...\n")

# First, calculate the empirical properties of the predictors from the raw data
predictor_props <- d %>%
  summarise(
    var_U = var(co),
    var_V = var(si),
    cov_UV = cov(co, si),
    E_U = mean(co),
    E_V = mean(si),
    E_U2 = mean(co^2),
    E_V2 = mean(si^2),
    E_UV = mean(co * si)
  )

# Helper function for the exact calculation
calculate_r2_exact_row <- function(row, props) {
  # Get all individual random sigma effects for this draw
  sigma_re_cols <- grep("r_id__sigma_Intercept", names(row), value = TRUE)
  sigma_re_values <- as.numeric(row[sigma_re_cols])
  
  # Calculate average residual variance
  sigma_i_values <- exp(row$b_sigma_Intercept + sigma_re_values)
  var_residual <- mean(sigma_i_values^2)
  
  # Variance of fixed effects
  var_fixed <- row$b_co^2 * props$var_U + row$b_si^2 * props$var_V + 2 * row$b_co * row$b_si * props$cov_UV
  
  # Variance of random effects
  var_m <- row$sd_id__Intercept^2
  var_c <- row$sd_id__co^2
  var_s <- row$sd_id__si^2
  cov_mc <- row$cor_id__Intercept__co * row$sd_id__Intercept * row$sd_id__co
  cov_ms <- row$cor_id__Intercept__si * row$sd_id__Intercept * row$sd_id__si
  cov_cs <- row$cor_id__co__si * row$sd_id__co * row$sd_id__si
  
  var_random <- var_m +
    var_c * props$E_U2 +
    var_s * props$E_V2 +
    2 * cov_mc * props$E_U +
    2 * cov_ms * props$E_V +
    2 * cov_cs * props$E_UV
    
  var_total <- var_fixed + var_random + var_residual
  
  tibble(
    R2_marginal_exact = var_fixed / var_total,
    R2_conditional_exact = (var_fixed + var_random) / var_total
  )
}

# Apply the function to each row
r2_exact_draws <- purrr::map_dfr(
  1:nrow(pop_draws),
  ~ calculate_r2_exact_row(pop_draws[.x, ], predictor_props)
)

if(any(is.na(r2_exact_draws))) {
  warning("NaNs detected in exact R2 calculation!")
  print(summary(r2_exact_draws))
}

# -----------------------------------------------------------------------------
# 5. Combine and Plot for Comparison
# -----------------------------------------------------------------------------

cat("Combining results and generating plot...\n")

# Combine all R² posterior distributions into a long-format data frame
comparison_data <- bind_cols(
  r2_simplified_draws %>% select(starts_with("R2")),
  r2_exact_draws %>% select(starts_with("R2"))
) %>%
  pivot_longer(everything(), names_to = "Metric", values_to = "R2")

# Create the comparison plot
p_r2_comparison <- ggplot(comparison_data, aes(x = R2, fill = Metric)) +
  geom_density(alpha = 0.6) +
  scale_fill_viridis_d(option = "D", name = "R² Metric") +
  labs(
    title = "Comparison of Population-Level R² Metrics",
    subtitle = "Posterior distributions for different calculation methods",
    x = "R² Value",
    y = "Density"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

print(p_r2_comparison)

# Save the plot
ggsave(here("figs", "r2_metric_comparison.pdf"), p_r2_comparison, width = 8, height = 6, device = cairo_pdf)
ggsave(here("figs", "r2_metric_comparison.png"), p_r2_comparison, width = 8, height = 6, dpi = 150)

cat("Plot saved to figs/r2_metric_comparison.pdf\n")
cat("Done!\n")