library(tidyverse)
library(GGally)

d_level1 <- read_rds("hap_person-estimates.rds")
d_level1 %>%
  cbind(id = 1:202) %>%
  saveRDS("hap_level1.rds")
d_level1_level2 <- read_rds("d_cors_without-summaries.rds")

d_level2 %>%
  filter()
  saveRDS("hap")

d_level1 %>%
  cbind(id = 1:202) %>%
  select(
    #intercept_median, amp_median,
    id,
    phi_lin_mean,
    phi_lin_median,
    phi_lin_ci_width,
    phi_circ_mean,
    phi_circ_median,
    phi_circ_ci_width) %>%
  pivot_longer(cols = contains("_"),
               names_to = "measure",
               values_to = "value") %>%
  mutate(
    value = value /(2*pi/24),
    type = case_when(
      str_detect(measure, "_lin") ~ "lin",
      str_detect(measure, "_circ") ~ "circ",
      TRUE ~ NA_character_),
    what = case_when(
      str_detect(measure, "mean") ~ "mean",
      str_detect(measure, "median") ~ "median",
      TRUE ~ "Width of 95% CI")
  ) %>%
  # select(-measure) %>%
  # pivot_wider(names_from = what,
  #             values_from = value) %>%
  ggplot() +
  aes(x = value) +
  geom_histogram(aes(y = after_stat(count / 202)),
                 bins = 50) +
  scale_y_continuous(
    #trans='sqrt',
    labels = scales::percent) +
  xlab("Hour") +
  ylab("Percent") +
  theme_light() +
  facet_grid(
    #scales = "free",
    cols = vars(type),
    rows = vars(what)
    )

d_level2 %>%
  filter(id == 1) %>%
  filter(iteration <= 1000) %>%
  pivot_wider(names_from = value_type,
              values_from = value) %>%
  mutate(significance_alpha = case_when(
    p_value < 0.01 ~ 1,
    p_value >= 0.01 & p_value < 0.05 ~ 0.7,
    TRUE ~ 0.4
    )
  ) %>%
  # filter(value_type == "cor") %>%
  ggplot() +
  aes(x = cor,
      color = significance_alpha,
      fill = significance_alpha) +
  geom_histogram(#aes(y = after_stat(count / 202)),
                 bins = 100) +
  facet_grid()

