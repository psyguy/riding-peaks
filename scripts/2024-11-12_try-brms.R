library(brms)
library(lme4)
library(psych)
library(Rfast)
library(tidyverse)

d <- d_l %>%
  filter(item == "hap") %>%
  mutate(t = hour_in_day,
         co = cos(2*pi/24 * t),
         si = sin(2*pi/24 * t)) %>%
  select(id, t, y, co, si)

(t_start <- Sys.time())
m_hap <- brm(
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

x <- as_draws_df(m_hap)

# Taking out the individual-specific sampled parameter estimates
d_draws <- as_draws_df(m_hap) %>%
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

d_cors <- d_draws %>%
  group_by(iteration) %>%
  mutate(phi_plus6 = (phi + pi/2) %% (2*pi),
         phi_plus12 = (phi + pi) %% (2*pi)) %>%
  mutate(
    # Ordinary correlations
    cor_si_intercept = cor(si, intercept),
    cor_co_intercept = cor(co, intercept),
    cor_si_co = cor(si, co),
    # Linear correlations: intercept with phi
    cor_intercept_phi_lin = cor(phi, intercept),
    cor_intercept_phi_lin_pval = cor.test(phi, intercept)$p.value,
    cor_intercept_phi_plus6_lin = cor(phi_plus6, intercept),
    cor_intercept_phi_plus6_lin_pval = cor.test(phi_plus6, intercept)$p.value,
    cor_intercept_phi_plus12_lin = cor(phi_plus12, intercept),
    cor_intercept_phi_plus12_lin_pval = cor.test(phi_plus12, intercept)$p.value,
    # Linear correlations: amp with phi
    cor_amp_phi_lin = cor(phi, amp),
    cor_amp_phi_lin_pval = cor.test(phi, amp)$p.value,
    cor_amp_phi_plus6_lin = cor(phi_plus6, amp),
    cor_amp_phi_plus6_lin_pval = cor.test(phi_plus6, amp)$p.value,
    cor_amp_phi_plus12_lin = cor(phi_plus12, amp),
    cor_amp_phi_plus12_lin_pval = cor.test(phi_plus12, amp)$p.value,
    # Circular correlations: intercept with phi
    cor_intercept_phi_circ = circlin.cor(phi, intercept)[1],
    cor_intercept_phi_circ_pval = circlin.cor(phi, intercept)[2],
    cor_intercept_phi_plus6_circ = circlin.cor(phi_plus6, intercept)[1],
    cor_intercept_phi_plus6_circ_pval = circlin.cor(phi_plus6, intercept)[2],
    cor_intercept_phi_plus12_circ = circlin.cor(phi_plus12, intercept)[1],
    cor_intercept_phi_plus12_circ_pval = circlin.cor(phi_plus12, intercept)[2],
    # Circular correlations: amp with phi
    cor_amp_phi_circ = circlin.cor(phi, amp)[1],
    cor_amp_phi_circ_pval = circlin.cor(phi, amp)[2],
    cor_amp_phi_plus6_circ = circlin.cor(phi_plus6, amp)[1],
    cor_amp_phi_plus6_circ_pval = circlin.cor(phi_plus6, amp)[2],
    cor_amp_phi_plus12_circ = circlin.cor(phi_plus12, amp)[1],
    cor_amp_phi_plus12_circ_pval = circlin.cor(phi_plus12, amp)[2]
  )

dd <- d_cors %>% filter(id < 50, iteration < 100)


ddl <- dd %>%
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
      str_detect(what, "phi_plus12") ~ "phi_plus12",
      str_detect(what, "phi_plus6") ~ "phi_plus6",
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
  select(-what) %>%
  # pivot_wider(names_from = value_type,
  #             values_from = value) %>%
  group_by(id) %>%
  mutate(
    phi_lin_mean = mean(phi),
    phi_lin_median = median(phi),
    phi_lin_ci_width = quantile(phi, 0.975) - quantile(phi, 0.025),
    phi_circ_mean = mean.circular(phi) %% (2 * pi) %>% as.numeric(),
    phi_circ_median = median.circular(phi) %% (2 * pi) %>% as.numeric(),
    phi_circ_ci_width = (
      quantile.circular(phi, 0.975) - quantile.circular(phi, 0.025)
    ) %% (2 * pi) %>% as.numeric()
  ) %>%
  group_by(iteration) %>%
  mutate(cor_summary_mean = mean(value),
         cor_summary_median = median(value))

# saveRDS(ddl-meh, "processed_summary_hap-meh.rds")
Sys.time()

#%>%
  mutate(
    value_type = ifelse(is.na(value_type), "correlation_value", "p_value")
  ) %>%
  pivot_wider(
    names_from = value_type,
    values_from = value
  )
