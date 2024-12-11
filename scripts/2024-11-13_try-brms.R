library(brms)
library(lme4)
library(psych)
library(Rfast)
library(tidyverse)
library(data.table)
library(circular)
library(parallel)
library(future.apply)

d <- d_l %>%
  filter(item == "hap") %>%
  mutate(t = hour_in_day,
         co = cos(2*pi/24 * t),
         si = sin(2*pi/24 * t)) %>%
  select(id, t, y, co, si)

(t_start <- Sys.time())

# Fixed residual variance
m_hap_fixed_var <- brm(
  y ~ 1 + co + si | id,
  data = d,
  backend = "cmdstanr",
  chains = 4,
  iter = 3000,
  threads = threading(2),
  cores = 8)

# Random residual variance
m_hap_random_var <- brm(
  y ~ 1 + co + si | id,
  data = d,
  backend = "cmdstanr",
  chains = 4,
  iter = 3000,
  threads = threading(2),
  cores = 8)
Sys.time() - t_start

# saveRDS(m_hap_meh, "brms_hap.rds")

m_hap <- read_rds("brms_hap.rds")

# x <- as_draws_df(m_hap)

# Taking out the individual-specific sampled parameter estimates
d_draws <- as_draws_df(mm) %>%
  select(contains("r_id[")) %>%
  mutate(iteration = 1:n()) %>%
  pivot_longer(
    cols = starts_with("r_id"),
    names_to = "id_x_variable",
    values_to = "value"
  ) %>%
  # Extract id (number between [ and ,)
  mutate(
    id = as.numeric(str_extract(id_x_variable, "(?<=\\[)[0-9]+(?=,)")),
    # Extract parameter (characters between , and ])
    parameter = str_extract(id_x_variable, "(?<=,)[^\\]]+(?=\\])"),
    # Clean up parameter by removing whitespace
    parameter = str_trim(parameter)
  ) %>%
  # Remove original variable column and reorder
  select(iteration, id, parameter, value) %>%
  pivot_wider(names_from = "parameter",
              values_from = "value") %>%
  mutate(phi = atan2(si, co) %% (2*pi),
         phi_atan = atan(si/co) %% (2*pi),
         psi = phi*24/(2*pi),
         psi_atan = phi_atan*24/(2*pi),
         amp = sqrt(si^2 + co^2))

names(d_draws) <- tolower(names(d_draws))

d_level1_level2 <- d_draws %>%
  group_by(iteration) %>%
  mutate(phi.plus6 = (phi + pi/2) %% (2*pi),
         phi.plus12 = (phi + pi) %% (2*pi)) %>%
  mutate(
    # Ordinary correlations
    cor_si_intercept = cor(si, intercept),
    cor_co_intercept = cor(co, intercept),
    cor_si_co = cor(si, co),
    # Linear correlations: intercept with phi
    cor_intercept_phi_lin = cor(phi, intercept),
    cor_intercept_phi_lin_pval = cor.test(phi, intercept)$p.value,
    cor_intercept_phi.plus6_lin = cor(phi.plus6, intercept),
    cor_intercept_phi.plus6_lin_pval = cor.test(phi.plus6, intercept)$p.value,
    cor_intercept_phi.plus12_lin = cor(phi.plus12, intercept),
    cor_intercept_phi.plus12_lin_pval = cor.test(phi.plus12, intercept)$p.value,
    # Linear correlations: amp with phi
    cor_amp_phi_lin = cor(phi, amp),
    cor_amp_phi_lin_pval = cor.test(phi, amp)$p.value,
    cor_amp_phi.plus6_lin = cor(phi.plus6, amp),
    cor_amp_phi.plus6_lin_pval = cor.test(phi.plus6, amp)$p.value,
    cor_amp_phi.plus12_lin = cor(phi.plus12, amp),
    cor_amp_phi.plus12_lin_pval = cor.test(phi.plus12, amp)$p.value,
    # Circular correlations: intercept with phi
    cor_intercept_phi_circ = circlin.cor(phi, intercept)[1],
    cor_intercept_phi_circ_pval = circlin.cor(phi, intercept)[2],
    cor_intercept_phi.plus6_circ = circlin.cor(phi.plus6, intercept)[1],
    cor_intercept_phi.plus6_circ_pval = circlin.cor(phi.plus6, intercept)[2],
    cor_intercept_phi.plus12_circ = circlin.cor(phi.plus12, intercept)[1],
    cor_intercept_phi.plus12_circ_pval = circlin.cor(phi.plus12, intercept)[2],
    # Circular correlations: amp with phi
    cor_amp_phi_circ = circlin.cor(phi, amp)[1],
    cor_amp_phi_circ_pval = circlin.cor(phi, amp)[2],
    cor_amp_phi.plus6_circ = circlin.cor(phi.plus6, amp)[1],
    cor_amp_phi.plus6_circ_pval = circlin.cor(phi.plus6, amp)[2],
    cor_amp_phi.plus12_circ = circlin.cor(phi.plus12, amp)[1],
    cor_amp_phi.plus12_circ_pval = circlin.cor(phi.plus12, amp)[2]
  ) %>%
  pivot_longer(
    cols = starts_with("cor_") &
      -contains("cor_si_intercept") &
      -contains("cor_co_intercept") & -contains("cor_si_co"),
    names_to = "what",
    values_to = "value"
  ) %>%
  mutate(
    linear_criterion = case_when(
      str_detect(what, "cor_amp") ~ "amp",
      str_detect(what, "cor_intercept") ~ "intercept",
      TRUE ~ NA_character_
    ),
    angle = case_when(
      str_detect(what, "phi.plus12") ~ "phi.plus12",
      str_detect(what, "phi.plus6") ~ "phi.plus6",
      str_detect(what, "phi") ~ "phi",
      TRUE ~ NA_character_
    ),
    correlation_type = case_when(
      str_detect(what, "_lin") ~ "lin",
      str_detect(what, "_circ") ~ "circ",
      TRUE ~ NA_character_
    ),
    value_type = if_else(str_detect(what, "pval"), "p_value", "cor"),
    .before = "value"
  ) %>%
  select(-what)

# saveRDS(d_level1_level2, "hap_level1_level2.rds")

d_level2 <- d_level1_level2 %>%
  filter(id == 1) %>%
  group_by(iteration) %>%
  select(cor_si_intercept:value) %>%
  pivot_wider(names_from = value_type,
              values_from = value)

# saveRDS(d_level2, "hap_level2.rds")

d_level1_level2 <- read_rds("hap_level1_level2.rds")
setDT(d_level1_level2)

library(parallel)

# Define function for each group
compute_stats <- function(data) {
  data[, .(
    ## Linears
    intercept_mean = mean(intercept),
    intercept_median = median(intercept),
    intercept_ci_width = quantile(intercept, 0.975) - quantile(intercept, 0.025),
    si_mean = mean(si),
    si_median = median(si),
    si_ci_width = quantile(si, 0.975) - quantile(si, 0.025),
    co_mean = mean(co),
    co_median = median(co),
    co_ci_width = quantile(co, 0.975) - quantile(co, 0.025),
    amp_mean = mean(amp),
    amp_median = median(amp),
    amp_ci_width = quantile(amp, 0.975) - quantile(amp, 0.025),
    ## Phi without offset
    phi_lin_mean = mean(phi),
    phi_lin_median = median(phi),
    phi_lin_ci_width = quantile(phi, 0.975) - quantile(phi, 0.025),
    phi_circ_mean = (mean.circular(phi) %% (2 * pi)) %>% as.numeric(),
    phi_circ_median = (median.circular(phi) %% (2 * pi)) %>% as.numeric(),
    phi_circ_ci_width = (
      (quantile.circular(phi, 0.975) - quantile.circular(phi, 0.025)) %% (2 * pi)
    ) %>% as.numeric(),
    ## Phi with quarter cycle (6h) offset
    phi.plus6_lin_mean = mean(phi.plus6),
    phi.plus6_lin_median = median(phi.plus6),
    phi.plus6_lin_ci_width = quantile(phi.plus6, 0.975) - quantile(phi.plus6, 0.025),
    phi.plus6_circ_mean = (mean.circular(phi.plus6) %% (2 * pi)) %>% as.numeric(),
    phi.plus6_circ_median = (median.circular(phi.plus6) %% (2 * pi)) %>% as.numeric(),
    phi.plus6_circ_ci_width = (
      (quantile.circular(phi.plus6, 0.975) - quantile.circular(phi.plus6, 0.025)) %% (2 * pi)
    ) %>% as.numeric(),
    ## Phi with half cycle (12h) offset
    phi.plus12_lin_mean = mean(phi.plus12),
    phi.plus12_lin_median = median(phi.plus12),
    phi.plus12_lin_ci_width = quantile(phi.plus12, 0.975) - quantile(phi.plus12, 0.025),
    phi.plus12_circ_mean = (mean.circular(phi.plus12) %% (2 * pi)) %>% as.numeric(),
    phi.plus12_circ_median = (median.circular(phi.plus12) %% (2 * pi)) %>% as.numeric(),
    phi.plus12_circ_ci_width = (
      (quantile.circular(phi.plus12, 0.975) - quantile.circular(phi.plus12, 0.025)) %% (2 * pi)
    ) %>% as.numeric()
  )]
}

(t_start <- Sys.time())
# Split data by 'id' for parallel processing
split_data <- d_level1_level2 %>%
  split(by = "id")
(t_start <- Sys.time())

# plan(multisession, workers = 8)

(t_start <- Sys.time())
# Apply the function in parallel using 11 cores
d_level1_first <-
  rbindlist(mclapply(split_data,
                     compute_stats,
                     mc.cores = 11)) %>%
  cbind(id = 1:202, .)
Sys.time() - t_start
Sys.time()



result %>%
  cbind(id = 1:202) %>%
  saveRDS("hap_level1_imoprtant.rds")


