# AR(1) Simulation for Tihomir's Review Comments
# Simulating N=200 individuals with T=500 time points each
#
# Data generating model:
# - Random means (mu_i)
# - Random AR parameter (phi_i): 80% near zero, 20% around 0.5
# - Random residual variance via rho_i = log(SD_i)
#
# Target correlations:
# - cor(phi_i, rho_i) = 0.20
# - cor(mu_i, rho_i) = -0.25

library(MASS)
library(tidyverse)
library(patchwork)

set.seed(20241211)

# -----------------------------------------------------------------------------
# Parameters
# -----------------------------------------------------------------------------
N <- 200        # Number of individuals
T_obs <- 500    # Time points per individual

# Group assignments: 80% low AR, 20% moderate AR
n_low <- 160    # Individuals 1-160: phi near 0
n_high <- 40    # Individuals 161-200: phi around 0.5

# -----------------------------------------------------------------------------
# Step 1: Generate phi_i (AR parameters) from mixture
# -----------------------------------------------------------------------------
# Low AR group: N(0, 0.05^2) - very close to zero
phi_low <- rnorm(n_low, mean = 0, sd = 0.05)

# High AR group: N(0.5, 0.10^2) - around 0.5
phi_high <- rnorm(n_high, mean = 0.5, sd = 0.10)

# Combine
phi <- c(phi_low, phi_high)

# -----------------------------------------------------------------------------
# Step 2: Generate rho_i and mu_i with EXACT target sample correlations
# -----------------------------------------------------------------------------
# Target correlations:
#   cor(phi, rho) = 0.20
#   cor(mu, rho) = -0.25
#   cor(mu, phi) = free (will be determined by the structure)
#
# Using post-hoc adjustment to achieve exact sample correlations

target_cor_phi_rho <- 0.20
target_cor_mu_rho <- -0.25

mu_rho <- 0       # Mean of log(SD), so SD ~ exp(0) = 1
sd_rho <- 0.3     # Variation in log(SD)
mu_mu <- 50       # Mean of the means
sd_mu <- 10       # SD of the means

# Generate independent standard normal variables
z1 <- rnorm(N)
z2 <- rnorm(N)

# Standardize phi (make it mean=0, sd=1 for correlation construction)
phi_std <- scale(phi)[, 1]

# Orthogonalize z1 with respect to phi_std to get truly independent component
z1_resid <- residuals(lm(z1 ~ phi_std))
z1_orth <- scale(z1_resid)[, 1]

# Construct rho with exact correlation to phi
# rho_std = r * phi_std + sqrt(1-r^2) * z_orth
rho_std <- target_cor_phi_rho * phi_std + sqrt(1 - target_cor_phi_rho^2) * z1_orth
rho <- mu_rho + sd_rho * rho_std

# Orthogonalize z2 with respect to rho_std
z2_resid <- residuals(lm(z2 ~ rho_std))
z2_orth <- scale(z2_resid)[, 1]

# Construct mu with exact correlation to rho
mu_std <- target_cor_mu_rho * scale(rho)[, 1] + sqrt(1 - target_cor_mu_rho^2) * z2_orth
mu <- mu_mu + sd_mu * mu_std

cat("Achieved cor(phi, rho):", round(cor(phi, rho), 3), "\n")
cat("Achieved cor(mu, rho):", round(cor(mu, rho), 3), "\n")
cat("Achieved cor(mu, phi):", round(cor(mu, phi), 3), "\n")

# -----------------------------------------------------------------------------
# Step 4: Create parameter table
# -----------------------------------------------------------------------------
params <- tibble(
  id = 1:N,
  group = c(rep("low_AR", n_low), rep("high_AR", n_high)),
  mu_i = mu,
  phi_i = phi,
  rho_i = rho,
  sd_i = exp(rho)  # Convert back to SD scale
)

# Print summary
cat("\n--- Parameter Summary ---\n")
print(summary(params[, c("mu_i", "phi_i", "rho_i", "sd_i")]))

cat("\n--- Correlation Matrix ---\n")
print(round(cor(params[, c("mu_i", "phi_i", "rho_i")]), 3))

cat("\n--- Summary by Group ---\n")
params %>%
  group_by(group) %>%
  summarise(
    n = n(),
    mean_phi = mean(phi_i),
    sd_phi = sd(phi_i),
    mean_rho = mean(rho_i),
    mean_sd = mean(sd_i),
    .groups = "drop"
  ) %>%
  print()

# -----------------------------------------------------------------------------
# Step 5: Visualize pairwise distributions
# -----------------------------------------------------------------------------

# Marginal distributions
p_phi <- ggplot(params, aes(x = phi_i, fill = group)) +
  geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
  scale_fill_manual(values = c("low_AR" = "#1b9e77", "high_AR" = "#d95f02")) +
  labs(x = expression(phi[i]), y = "Count", title = "AR(1) Parameter Distribution") +
  theme_minimal() +
  theme(legend.position = "bottom")

p_rho <- ggplot(params, aes(x = rho_i, fill = group)) +
  geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
  scale_fill_manual(values = c("low_AR" = "#1b9e77", "high_AR" = "#d95f02")) +
  labs(x = expression(rho[i] == log(SD[i])), y = "Count", title = "Log-SD Distribution") +
  theme_minimal() +
  theme(legend.position = "bottom")

p_mu <- ggplot(params, aes(x = mu_i, fill = group)) +
  geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
  scale_fill_manual(values = c("low_AR" = "#1b9e77", "high_AR" = "#d95f02")) +
  labs(x = expression(mu[i]), y = "Count", title = "Mean Distribution") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Pairwise scatterplots
p_phi_rho <- ggplot(params, aes(x = phi_i, y = rho_i, color = group)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
  scale_color_manual(values = c("low_AR" = "#1b9e77", "high_AR" = "#d95f02")) +
  labs(
    x = expression(phi[i]),
    y = expression(rho[i]),
    title = paste0("cor = ", round(cor(phi, rho), 2))
  ) +
  theme_minimal() +
  theme(legend.position = "none")

p_mu_rho <- ggplot(params, aes(x = mu_i, y = rho_i, color = group)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
  scale_color_manual(values = c("low_AR" = "#1b9e77", "high_AR" = "#d95f02")) +
  labs(
    x = expression(mu[i]),
    y = expression(rho[i]),
    title = paste0("cor = ", round(cor(mu, rho), 2))
  ) +
  theme_minimal() +
  theme(legend.position = "none")

p_mu_phi <- ggplot(params, aes(x = mu_i, y = phi_i, color = group)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
  scale_color_manual(values = c("low_AR" = "#1b9e77", "high_AR" = "#d95f02")) +
  labs(
    x = expression(mu[i]),
    y = expression(phi[i]),
    title = paste0("cor = ", round(cor(mu, phi), 2))
  ) +
  theme_minimal() +
  theme(legend.position = "none")

# Combine plots
p_combined <- (p_phi | p_rho | p_mu) / (p_phi_rho | p_mu_rho | p_mu_phi) +
  plot_annotation(
    title = "AR(1) Simulation: Parameter Distributions",
    subtitle = paste0("N = ", N, " individuals (", n_low, " low AR, ", n_high, " high AR)"),
    theme = theme(plot.title = element_text(size = 14, face = "bold"))
  )

print(p_combined)

# Save the parameter visualization
ggsave("figs/ar1_simulation_params.pdf", p_combined, width = 12, height = 8)
ggsave("figs/ar1_simulation_params.png", p_combined, width = 12, height = 8, dpi = 150)

# -----------------------------------------------------------------------------
# Step 6: Simulate AR(1) time series for each individual
# -----------------------------------------------------------------------------
cat("\n--- Simulating AR(1) time series ---\n")

simulate_ar1 <- function(T, mu, phi, sigma) {
  # AR(1) process: y_t = mu + phi * (y_{t-1} - mu) + epsilon_t
  # where epsilon_t ~ N(0, sigma^2)

  y <- numeric(T)
  # Initialize at stationary mean
  y[1] <- mu + rnorm(1, 0, sigma / sqrt(1 - phi^2 + 0.001))  # Small constant for stability


  for (t in 2:T) {
    y[t] <- mu + phi * (y[t-1] - mu) + rnorm(1, 0, sigma)
  }

  return(y)
}

# Generate data for all individuals
data_list <- vector("list", N)

for (i in 1:N) {
  # Ensure phi is in stationary range
  phi_i_bounded <- pmin(pmax(params$phi_i[i], -0.99), 0.99)

  y <- simulate_ar1(
    T = T_obs,
    mu = params$mu_i[i],
    phi = phi_i_bounded,
    sigma = params$sd_i[i]
  )

  data_list[[i]] <- tibble(
    id = i,
    group = params$group[i],
    time = 1:T_obs,
    y = y
  )
}

# Combine into single data frame
sim_data <- bind_rows(data_list)

cat("Simulated data dimensions:", nrow(sim_data), "rows x", ncol(sim_data), "columns\n")
cat("Total observations:", N, "individuals x", T_obs, "time points =", N * T_obs, "\n")

# -----------------------------------------------------------------------------
# Step 7: Verify simulation by estimating parameters
# -----------------------------------------------------------------------------
cat("\n--- Verifying simulation: Estimating AR(1) parameters ---\n")

estimate_ar1 <- function(y) {
  # Simple AR(1) estimation
  n <- length(y)
  y_mean <- mean(y)
  y_centered <- y - y_mean

  # Estimate phi using autocorrelation at lag 1
  phi_hat <- cor(y_centered[-n], y_centered[-1])

  # Estimate residual SD
  residuals <- y_centered[-1] - phi_hat * y_centered[-n]
  sigma_hat <- sd(residuals)

  return(c(mu_hat = y_mean, phi_hat = phi_hat, sigma_hat = sigma_hat))
}

# Estimate for each individual
estimates <- sim_data %>%
  group_by(id) %>%
  summarise(
    est = list(estimate_ar1(y)),
    .groups = "drop"
  ) %>%
  mutate(
    mu_hat = map_dbl(est, ~ .x["mu_hat"]),
    phi_hat = map_dbl(est, ~ .x["phi_hat"]),
    sigma_hat = map_dbl(est, ~ .x["sigma_hat"])
  ) %>%
  select(-est)

# Merge with true parameters
comparison <- params %>%
  left_join(estimates, by = "id")

# Compare true vs estimated
cat("\nCorrelation between true and estimated parameters:\n")
cat("  mu:    ", round(cor(comparison$mu_i, comparison$mu_hat), 3), "\n")
cat("  phi:   ", round(cor(comparison$phi_i, comparison$phi_hat), 3), "\n")
cat("  sigma: ", round(cor(comparison$sd_i, comparison$sigma_hat), 3), "\n")

# Visualize recovery
p_recovery_mu <- ggplot(comparison, aes(x = mu_i, y = mu_hat, color = group)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("low_AR" = "#1b9e77", "high_AR" = "#d95f02")) +
  labs(x = expression(True~mu[i]), y = expression(Estimated~hat(mu)[i]), title = "Mean Recovery") +
  theme_minimal() +
  theme(legend.position = "none")

p_recovery_phi <- ggplot(comparison, aes(x = phi_i, y = phi_hat, color = group)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("low_AR" = "#1b9e77", "high_AR" = "#d95f02")) +
  labs(x = expression(True~phi[i]), y = expression(Estimated~hat(phi)[i]), title = "AR(1) Recovery") +
  theme_minimal() +
  theme(legend.position = "none")

p_recovery_sigma <- ggplot(comparison, aes(x = sd_i, y = sigma_hat, color = group)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("low_AR" = "#1b9e77", "high_AR" = "#d95f02")) +
  labs(x = expression(True~sigma[i]), y = expression(Estimated~hat(sigma)[i]), title = "SD Recovery") +
  theme_minimal() +
  theme(legend.position = "bottom")

p_recovery <- p_recovery_mu | p_recovery_phi | p_recovery_sigma
p_recovery <- p_recovery + plot_annotation(
  title = "Parameter Recovery Check",
  subtitle = "True vs. Estimated Parameters"
)

print(p_recovery)
ggsave("figs/ar1_simulation_recovery.pdf", p_recovery, width = 12, height = 4)
ggsave("figs/ar1_simulation_recovery.png", p_recovery, width = 12, height = 4, dpi = 150)

# -----------------------------------------------------------------------------
# Step 8: Plot example time series
# -----------------------------------------------------------------------------
# Select a few example individuals from each group
example_ids <- c(sample(1:n_low, 3), sample((n_low+1):N, 3))

p_ts <- sim_data %>%
  filter(id %in% example_ids) %>%
  left_join(params %>% select(id, group, phi_i), by = c("id", "group")) %>%
  mutate(label = paste0("ID ", id, " (phi = ", round(phi_i, 2), ")")) %>%
  ggplot(aes(x = time, y = y, color = group)) +
  geom_line(alpha = 0.7) +
  facet_wrap(~ label, scales = "free_y", ncol = 2) +
  scale_color_manual(values = c("low_AR" = "#1b9e77", "high_AR" = "#d95f02")) +
  labs(
    x = "Time",
    y = "y",
    title = "Example AR(1) Time Series",
    subtitle = "3 individuals from each group"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p_ts)
ggsave("figs/ar1_simulation_timeseries.pdf", p_ts, width = 10, height = 8)
ggsave("figs/ar1_simulation_timeseries.png", p_ts, width = 10, height = 8, dpi = 150)

# -----------------------------------------------------------------------------
# Step 9: Save outputs
# -----------------------------------------------------------------------------
# Save parameter table
write_csv(params, "fits/ar1_simulation_params.csv")

# Save simulated data (warning: large file)
write_csv(sim_data, "fits/ar1_simulation_data.csv")

cat("\n--- Outputs saved ---\n")
cat("Parameter table: fits/ar1_simulation_params.csv\n
")
cat("Simulated data: fits/ar1_simulation_data.csv\n")
cat("Figures: figs/ar1_simulation_*.pdf/png\n")

cat("\n--- Simulation complete! ---\n")

# =============================================================================
# PART 2: Fit Multilevel AR(1) Model with brms
# =============================================================================

library(brms)

# -----------------------------------------------------------------------------
# Step 10: Prepare data for brms AR(1) model
# -----------------------------------------------------------------------------
# For AR(1) modeling, we need to create lagged observations
# The model is: y_t = mu_i + phi_i * (y_{t-1} - mu_i) + epsilon_it
# Which can be rewritten as: y_t = mu_i * (1 - phi_i) + phi_i * y_{t-1} + epsilon_it

cat("\n--- Preparing data for brms ---\n")

# Create lagged variable within each individual
sim_data_lagged <- sim_data %>%
  arrange(id, time) %>%
  group_by(id) %>%
  mutate(
    y_lag = lag(y, 1)
  ) %>%
  ungroup() %>%
  filter(!is.na(y_lag))  # Remove first observation per person (no lag available)

cat("Data prepared:", nrow(sim_data_lagged), "observations\n")
cat("(Removed", N, "first observations due to lagging)\n")

# -----------------------------------------------------------------------------
# Step 11: Fit brms model with random mean, AR, and residual variance
# -----------------------------------------------------------------------------
# Model specification:
#   y ~ 1 + y_lag + (1 + y_lag | id)  for the mean structure
#   sigma ~ 1 + (1 | id)              for the residual variance
#
# This parameterizes:
#   y_t = (Intercept_i) + (phi_i) * y_{t-1} + epsilon_it
#   log(sigma_it) = rho_i
#
# Note: The intercept here is mu_i * (1 - phi_i), not mu_i directly
# We'll recover mu_i in post-processing

cat("\n--- Fitting brms model ---\n")
cat("This may take a while...\n")

(st <- Sys.time())

# Fit the model
m_ar1 <- brm(
 brmsformula(
    y ~ 1 + y_lag + (1 + y_lag | i | id),
    sigma ~ 1 + (1 | i | id)
  ),
  data = sim_data_lagged,
  backend = "cmdstanr",
  chains = 4,
  iter = 2000,
  warmup = 1000,
  threads = threading(2),
  cores = 8,
  seed = 2025-12-11,
  file = "fits/brms_ar1_simulation"
)

elapsed <- Sys.time() - st
cat("\nModel fitting completed in", round(as.numeric(elapsed, units = "mins"), 1), "minutes\n")

# Save the model
saveRDS(m_ar1, "fits/brms_ar1_simulation.rds")

# -----------------------------------------------------------------------------
# Step 12: Extract and compare estimated parameters
# -----------------------------------------------------------------------------
cat("\n--- Extracting model estimates ---\n")

# Get summary of fixed effects
cat("\nFixed Effects:\n")
print(fixef(m_ar1))

# Get summary of random effects SDs and correlations
cat("\nRandom Effects (Group-Level):\n")
print(VarCorr(m_ar1))

# Extract individual-specific estimates (posterior means of random effects)
ranef_summary <- ranef(m_ar1, summary = TRUE)

# Extract random effects for the mean equation
ranef_mean <- ranef_summary$id[, , "Intercept"]
ranef_phi <- ranef_summary$id[, , "y_lag"]
ranef_sigma <- ranef_summary$id[, , "sigma_Intercept"]

# Get fixed effects
fixef_vals <- fixef(m_ar1)
b_intercept <- fixef_vals["Intercept", "Estimate"]
b_phi <- fixef_vals["y_lag", "Estimate"]
b_logsigma <- fixef_vals["sigma_Intercept", "Estimate"]

# Compute individual-specific parameter estimates
# Note: The model parameterizes y = a + phi*y_lag, where a = mu*(1-phi)
# So mu = a / (1 - phi) for each individual

estimates_brms <- tibble(
  id = 1:N,
  # Raw estimates from model
  intercept_hat = b_intercept + ranef_mean[, "Estimate"],
  phi_hat_brms = b_phi + ranef_phi[, "Estimate"],
  logsigma_hat = b_logsigma + ranef_sigma[, "Estimate"],
  sigma_hat_brms = exp(logsigma_hat),
  # Recover mu: since y = mu*(1-phi) + phi*y_lag, we have intercept = mu*(1-phi)
  mu_hat_brms = intercept_hat / (1 - phi_hat_brms)
)

# Merge with true parameters
comparison_brms <- params %>%
  left_join(estimates_brms, by = "id")

# -----------------------------------------------------------------------------
# Step 13: Compare true vs estimated parameters
# -----------------------------------------------------------------------------
cat("\n--- Parameter Recovery (brms estimates) ---\n")
cat("\nCorrelation between true and estimated parameters:\n")
cat("  mu:    ", round(cor(comparison_brms$mu_i, comparison_brms$mu_hat_brms), 3), "\n")
cat("  phi:   ", round(cor(comparison_brms$phi_i, comparison_brms$phi_hat_brms), 3), "\n")
cat("  sigma: ", round(cor(comparison_brms$sd_i, comparison_brms$sigma_hat_brms), 3), "\n")

# Compare with simple OLS estimates
cat("\nComparison with simple OLS estimates:\n")
comparison_both <- comparison_brms %>%
  left_join(comparison %>% select(id, mu_hat, phi_hat, sigma_hat), by = "id")

cat("  mu (OLS):    ", round(cor(comparison_both$mu_i, comparison_both$mu_hat), 3), "\n")
cat("  mu (brms):   ", round(cor(comparison_both$mu_i, comparison_both$mu_hat_brms), 3), "\n")
cat("  phi (OLS):   ", round(cor(comparison_both$phi_i, comparison_both$phi_hat), 3), "\n")
cat("  phi (brms):  ", round(cor(comparison_both$phi_i, comparison_both$phi_hat_brms), 3), "\n")
cat("  sigma (OLS): ", round(cor(comparison_both$sd_i, comparison_both$sigma_hat), 3), "\n")
cat("  sigma (brms):", round(cor(comparison_both$sd_i, comparison_both$sigma_hat_brms), 3), "\n")

# -----------------------------------------------------------------------------
# Step 14: Visualize brms parameter recovery
# -----------------------------------------------------------------------------
cat("\n--- Creating brms recovery plots ---\n")

p_brms_mu <- ggplot(comparison_brms, aes(x = mu_i, y = mu_hat_brms, color = group)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("low_AR" = "#1b9e77", "high_AR" = "#d95f02")) +
  labs(x = expression(True~mu[i]), y = expression(brms~hat(mu)[i]),
       title = "Mean Recovery (brms)") +
  theme_minimal() +
  theme(legend.position = "none")

p_brms_phi <- ggplot(comparison_brms, aes(x = phi_i, y = phi_hat_brms, color = group)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("low_AR" = "#1b9e77", "high_AR" = "#d95f02")) +
  labs(x = expression(True~phi[i]), y = expression(brms~hat(phi)[i]),
       title = "AR(1) Recovery (brms)") +
  theme_minimal() +
  theme(legend.position = "none")

p_brms_sigma <- ggplot(comparison_brms, aes(x = sd_i, y = sigma_hat_brms, color = group)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("low_AR" = "#1b9e77", "high_AR" = "#d95f02")) +
  labs(x = expression(True~sigma[i]), y = expression(brms~hat(sigma)[i]),
       title = "SD Recovery (brms)") +
  theme_minimal() +
  theme(legend.position = "bottom")

p_brms_recovery <- p_brms_mu | p_brms_phi | p_brms_sigma
p_brms_recovery <- p_brms_recovery + plot_annotation(
  title = "Parameter Recovery Check (brms multilevel model)",
  subtitle = "True vs. Estimated Parameters"
)

print(p_brms_recovery)
ggsave("figs/ar1_simulation_recovery_brms.pdf", p_brms_recovery, width = 12, height = 4)
ggsave("figs/ar1_simulation_recovery_brms.png", p_brms_recovery, width = 12, height = 4, dpi = 150)

# -----------------------------------------------------------------------------
# Step 15: Check recovery of correlations between parameters
# -----------------------------------------------------------------------------
cat("\n--- Correlation Recovery ---\n")

cat("\nTrue correlations:\n")
cat("  cor(phi, rho):  ", round(cor(params$phi_i, params$rho_i), 3), "\n")
cat("  cor(mu, rho):   ", round(cor(params$mu_i, params$rho_i), 3), "\n")
cat("  cor(mu, phi):   ", round(cor(params$mu_i, params$phi_i), 3), "\n")

cat("\nEstimated correlations (brms):\n")
cat("  cor(phi, logsigma): ", round(cor(comparison_brms$phi_hat_brms, comparison_brms$logsigma_hat), 3), "\n")
cat("  cor(mu, logsigma):  ", round(cor(comparison_brms$mu_hat_brms, comparison_brms$logsigma_hat), 3), "\n")
cat("  cor(mu, phi):       ", round(cor(comparison_brms$mu_hat_brms, comparison_brms$phi_hat_brms), 3), "\n")

# Save comparison table
write_csv(comparison_brms, "fits/ar1_simulation_comparison_brms.csv")

cat("\n--- All analyses complete! ---\n")
