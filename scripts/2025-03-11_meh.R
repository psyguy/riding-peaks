Sys.time()
fixed_phi_at_12 <- d_draws %>%
  group_by(item, iteration) %>%
  summarise(
    bind_cols(
      # Linear calculations
      compute_linear_stats(cur_data(), "phi.begins12"),
      # Circular calculations
      compute_circular_stats(cur_data(), "phi.begins12")
    )) %>%
  select(!contains("mean")) %>%
  pivot_longer(!contains("width") & contains("phi"),
               names_to = "type_across_i",
               values_to = "phi_hat") %>%
  mutate(type_across_i = case_when(
    grepl("circ", type_across_i) ~ "circ_across_i",
    TRUE ~ "lin_across_i"
  ),
  phi_hat = (phi_hat + pi) %% (2*pi)) %>%
  select(!contains("width")) %>%
  group_by(item, type_across_i) %>%
  summarise(
    bind_cols(
      # Linear calculations
      compute_linear_stats(cur_data(), "phi_hat"),
      # Circular calculations
      compute_circular_stats(cur_data(), "phi_hat")
    )) %>%
  select(!contains("mean"))  %>%
  pivot_longer(
    cols = starts_with("phi_hat"),
    names_to = c("circ_flag", "metric"),
    names_pattern = "phi_hat(_circ)?_(.*)",
    values_drop_na = FALSE
  ) %>%
  mutate(
    metric_type = ifelse(circ_flag == "_circ", "circ", "lin")
  ) %>%
  select(-circ_flag) %>%
  pivot_wider(names_from = "metric",
              values_from = "value") %>%
  mutate(median = (median + pi) %>% `%%`(2*pi) %>% `*`(12/pi) %>% round(2),
         ci_width = ci_width*12/pi)
beep(1)
