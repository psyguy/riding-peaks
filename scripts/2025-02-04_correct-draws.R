item_ <- "fitbit"

library(dplyr)
library(tidyr)

dr_new <- dr %>%
  mutate(iteration = row_number()) %>%
  pivot_longer(
    cols = starts_with("r_id"),
    names_pattern = "r_id(_{0,2}sigma)?\\[(\\d+),(.*)\\]",
    names_to = c("sigma_part", "id", "term"),
    values_to = "value"
  ) %>%
  mutate(term = if_else(sigma_part == "__sigma", paste0("sigma_", term), term)) %>%
  pivot_wider(
    id_cols = c(iteration, id, b_Intercept, b_co, b_si, b_sigma_Intercept),
    names_from = term,
    values_from = value
  ) %>%
  mutate(
    mesor        = Intercept + b_Intercept,
    co               = co + b_co,
    si               = si + b_si,
    logsigma  = sigma_Intercept + b_sigma_Intercept
  ) %>%
  select(iteration, id, Intercept, co, si, logsigma,
         b_Intercept, b_co, b_si, b_sigma_Intercept) %>%
  rename(b_logsigma = b_sigma_Intercept) %>%
  mutate(
    sigma = exp(logsigma),
    sigma2 = exp(logsigma)^2,
    amp = sqrt(si^2 + co^2),
    phi = atan2(si, co) %% (2*pi),
    phi_atan = atan(si/co) %% (2*pi)
  ) %>%
  mutate(phi.begins6 = (phi - pi/2) %% (2*pi),
         phi.begins12 = (phi - pi) %% (2*pi))

