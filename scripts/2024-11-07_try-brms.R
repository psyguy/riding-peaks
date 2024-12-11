library(brms)
library(lme4)
library(psych)

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

# saveRDS(m_hap, "brms_hap.rds")

m %>% brms::as_draws()

x_draws <- x %>%
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
              values_from = "value")

names(x_draws) <- tolower(names(x_draws))

xx <- x_draws %>%
  mutate(psi = (24/(2*pi)*atan2(si, co)) %% 24,
         psi_atan = 24/(2*pi)*atan(si/co),
         A = sqrt(si^2 + co^2)) %>%
  group_by(iteration) %>%
  mutate(cor_si_intercept = cor(si, intercept),
         cor_co_intercept = cor(co, intercept),
         cor_si_co = cor(si, co),
         cor_psi_intercept = cor(psi, intercept),
         cor_psi_intercept_circ = psych::circadian.linear.cor(psi, intercept),
         cor_psi_plus6_intercept = cor((psi + 6)%%24, intercept),
         cor_psi_plus6_intercept_circ = psych::circadian.linear.cor((psi + 6)%%24, intercept)
         )

x_psi <- x_draws %>%
  mutate(phi = atan2(si, co) %% (2*pi),
         phi_atan = atan(si/co) %% (2*pi),
         psi = phi*24/(2*pi),
         psi_atan = phi_atan*24/(2*pi),
         A = sqrt(si^2 + co^2)) %>%
  # select(iteration, id, intercept, A, psi) %>%
  group_by(iteration) %>%
  mutate(cor_circ = Rfast::circlin.cor(phi, A)[1],
         cor_circ_pval = Rfast::circlin.cor(phi, A)[2],
         cor_circ_int = Rfast::circlin.cor(phi, intercept)[1],
         cor_circ_pval_int = Rfast::circlin.cor(phi, intercept)[2],
         # cor_plus1_circ = Rfast::circlin.cor((phi + 1)%%(2*pi), intercept),
         cor_lin = cor(phi, A),
         cor_plus1_lin = cor((phi + 2)%%(2*pi), A)
         )
  # mutate(
  #   cor_lin = cor(psi, intercept),
  #   cor_circ = circadian.linear.cor(psi, intercept),
  #   cor_plus6_lin = cor((psi + 6)%%24, intercept),
  #   cor_plus6_circ = circadian.linear.cor((psi + 6)%%24, intercept)
  # ) # %>%
  # group_by(id) %>%
  # mutate(
  #   mean_psi_lin = mean(psi),
  #   mean_psi_circ = mean.circular(psi, units = "hours"),
  #   median_psi_lin = median(psi),
  #   median_psi_circ = median.circular(psi, units = "hours")
  # )

x_psi_small <- x_psi %>% filter(id == 1)

plot(x_psi_small$cor_circ_int, x_psi_small$cor_circ_pval_int)

mean(x_psi_small$cor_circ_pval); median(x_psi_small$cor_circ_pval)

mean(x_psi_small$cor_circ_pval_int); median(x_psi_small$cor_circ_pval_int)


xx_phi <- x_psi$phi %>% head(1000) %>% as.circular()
xx_a <- x_psi$A %>% head(1000)
xx_int <- x_psi$intercept %>% head(1000)

circadian.linear.cor((xx_phi - 2)%%(2*pi), c(xx_int, xx_a), hours = FALSE)
cor((xx_phi)%%(2*pi), xx_int)

Rfast::circlin.cor((xx_phi - 12)%%(2*pi), xx_a)



mm <- lme4::lmer(y ~ (1 + co + si | id),
           data = d %>% filter(id <= 50))

ranefs <- ranef(mm)$id

ranefs$A <- sqrt(ranefs$co^2 + ranefs$si^2)
ranefs$psi <- (-10 + atan2(ranefs$si, ranefs$co) * 24/(2*pi)) %% 24

(correlation <- cor(ranefs$A, ranefs$psi))



# # GPT steps ---------------------------------------------------------------
#
# model <- brm(
#   formula = y ~ 1 + co + si | id,  # random intercept and random slopes for co and si by participant
#   data = d,
#   backend = "cmdstanr",
#   chains = 2,
#   cores = 2,
#   iter = 2000/2,  # number of iterations for MCMC
#   warmup = 1000/2
# )
#
# # Extract random slopes for 'co' and 'si' by participant
# random_co_draws <- posterior_samples(model, pars = "^r_id\\[.*,co\\]")
# random_si_draws <- posterior_samples(model, pars = "^r_id\\[.*,si\\]")
