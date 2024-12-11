library(brms)
library(lme4)
library(psych)
library(Rfast)
library(tidyverse)
library(rlang)
library(data.table)
library(circular)
library(parallel)
library(future.apply)

d_l <- readRDS("d_leuven_long.rds")

d <- d_l %>%
  filter(item == "hap") %>%
  mutate(t = hour_in_day,
         co = cos(2*pi/24 * t),
         si = sin(2*pi/24 * t)) %>%
  select(id, t, y, co, si)

Sys.time()
m_hap_random_var <- brm(
  brmsformula(
    y ~ 1 + co + si | id,
    sigma ~ 1 | id),
  data = d,
  backend = "cmdstanr",
  chains = 4,
  iter = 3000,
  threads = threading(2),
  cores = 8)

# saveRDS(m_hap_random_var,
#         here::here(
#                    "brms_hap_random_var.rds"))

# Selecting the second model from now on
m_hap <- m_hap_random_var

# Taking out the individual-specific sampled parameter estimates
d_draws <- as_draws_df(m_hap) %>%
  select(contains("r_id[") | contains("r_id__sigma[")) %>%
  mutate(iteration = 1:n()) %>%
  pivot_longer(
    cols = starts_with("r_id"),
    names_to = "id_x_variable",
    values_to = "value"
  ) %>%
  # Extract id (number between [ and ,)
  mutate(
    id_x_variable =
      str_replace(id_x_variable,
                  "r_id__sigma\\[(\\d+),Intercept\\]",
                  "r_id__sigma[\\1,sigmaIntercept]"),
    id = str_extract(id_x_variable,
                     "(?<=\\[)[0-9]+(?=,)") %>%
      as.numeric(),
    # Extract parameter (characters between , and ])
    parameter = str_extract(id_x_variable,
                            "(?<=,)[^\\]]+(?=\\])"),
    # Clean up parameter by removing whitespace
    parameter = str_trim(parameter)
  ) %>%
  # Remove original variable column and reorder
  select(iteration, id, parameter, value) %>%
  pivot_wider(names_from = "parameter",
              values_from = "value") %>%
  mutate(
    ## Not sure if sigmaIntercept is on log scale or not
    ## Commenting it out now
    log_sigma = sigmaIntercept,
    sigma2 = exp(sigmaIntercept)^2,
    # sigma2 = sigmaIntercept^2,
    # log_sigma2 = log(sigma2),
    mesor = Intercept,
    amp = sqrt(si^2 + co^2),
    phi = atan2(si, co) %% (2*pi),
    phi_atan = atan(si/co) %% (2*pi),
    psi = phi*24/(2*pi),
    psi_atan = phi_atan*24/(2*pi)
  ) %>%
  mutate(phi.begins6 = (phi - pi/2) %% (2*pi),
         phi.begins12 = (phi - pi) %% (2*pi))


library(dplyr)
library(purrr)
library(circular)

# Helper function for linear stats
compute_linear_stats <- function(data, variable) {
  tibble(
    mean = mean(data[[variable]], na.rm = TRUE),
    median = median(data[[variable]], na.rm = TRUE),
    ci_width = quantile(data[[variable]], 0.975, na.rm = TRUE) - quantile(data[[variable]], 0.025, na.rm = TRUE)
  ) %>%
    rename_with(~ paste0(variable, "_", .))
}

# Helper function for circular stats
compute_circular_stats <- function(data, variable) {
  tibble(
    circ_mean = (mean.circular(data[[variable]]) %% (2 * pi)) %>% as.numeric(),
    circ_median = (median.circular(data[[variable]]) %% (2 * pi)) %>% as.numeric(),
    circ_ci_width = (
      (quantile.circular(data[[variable]], 0.975) - quantile.circular(data[[variable]], 0.025)) %% (2 * pi)
    ) %>% as.numeric()
  ) %>%
    rename_with(~ paste0(variable, "_", .))
}

# Apply calculations for all variables
ests_level1 <- d_draws %>%
  group_by(id) %>%
  summarise(
    bind_cols(
      # Linear calculations
      compute_linear_stats(cur_data(), "mesor"),
      compute_linear_stats(cur_data(), "si"),
      compute_linear_stats(cur_data(), "co"),
      compute_linear_stats(cur_data(), "amp"),
      compute_linear_stats(cur_data(), "sigma2"),
      compute_linear_stats(cur_data(), "log_sigma"),
      compute_linear_stats(cur_data(), "phi"),
      compute_linear_stats(cur_data(), "phi.begins6"),
      compute_linear_stats(cur_data(), "phi.begins12"),

      # Circular calculations
      compute_circular_stats(cur_data(), "phi"),
      compute_circular_stats(cur_data(), "phi.begins6"),
      compute_circular_stats(cur_data(), "phi.begins12")
    )
  )

saveRDS(ests_level1, "ests_level1_hap_random_var.rds")

# Create correlation function with linear and circular methods
compute_correlations <- function(x, y) {
  tibble(
    # Linear correlations
    lin_cor = cor(x, y),
    lin_pval = cor.test(x, y)$p.value,
    # Circular correlations (if circlin.cor function exists)
    circ_cor = tryCatch(
      circlin.cor(x, y)[1],
      error = function(e)
        NA_real_
    ),
    circ_pval = tryCatch(
      circlin.cor(x, y)[2],
      error = function(e)
        NA_real_
    )
  )
}

# Compute correlations for all combinations
cors_lin_angle <- d_draws %>%
  group_by(iteration) %>%
  summarise(
    crossing(
      linear_variable = c('mesor', 'log_sigma', 'sigma2', 'amp'),
      circular_variable = c('phi', 'phi.begins6', 'phi.begins12')
    ) %>%
      mutate(corr_data = map2(
        .x = syms(circular_variable),
        .y = syms(linear_variable),
        .f = ~ compute_correlations(eval_tidy(.x, data = cur_data()),
                                    eval_tidy(.y, data = cur_data()))
      )) %>%
      unnest(corr_data) %>%
      rename(par2 = circular_variable,
             par1 = linear_variable)
  ) %>%
  unnest(cols = last_col()) %>%
  pivot_longer(lin_cor:circ_pval,
               names_to = "cortypeXmeasure",
               values_to = "value")

cors_lin_lin <- d_draws %>%
  group_by(iteration) %>%
  summarise(
    pair_data = list(
      tibble(
        par1 = c("mesor", "mesor", "mesor", "amp", "amp"),
        par2 = c("amp", "log_sigma", "sigma2", "log_sigma", "sigma2")
      )
    )
  ) %>%
  unnest(pair_data) %>%
  mutate(
    x = map(par1, ~ pull(filter(d_draws, iteration == cur_group_id()), all_of(.))),
    y = map(par2, ~ pull(filter(d_draws, iteration == cur_group_id()), all_of(.))),
    lin_cor = map2_dbl(x, y, ~ cor(.x, .y, use = "complete.obs")),
    lin_pval = map2_dbl(x, y, ~ cor.test(.x, .y)$p.value)
  ) %>%
  select(-x, -y) %>%
  unnest(cols = last_col()) %>%
  pivot_longer(lin_cor:lin_pval,
               names_to = "cortypeXmeasure",
               values_to = "value")


ests_level2 <- cors_lin_angle %>%
  rbind(cors_lin_lin) %>%
  rename(cortypeXmeasure = what) %>%
  group_by(cortypeXmeasure) %>%
  mutate(
    correlation_type = str_split(cortypeXmeasure, "_")[[1]][1],
    measure = str_split(cortypeXmeasure, "_")[[1]][2],
    .before = value
  ) %>%
  ungroup() %>%
  mutate(what = paste(par1, correlation_type, par2, sep = "_"),
         .after = correlation_type) %>%
  select(-cortypeXmeasure)

saveRDS(ests_level2, "ests_level2_hap_random_var.rds")
