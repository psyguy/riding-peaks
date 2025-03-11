el1.5 %>%
  pivot_longer(
    cols = matches("^phi"),
    names_to = "phi_type",
    values_to = "phi"
  ) %>%
  mutate(
    method = ifelse(str_detect(phi_type, "_circ"), "circular", "linear"),
    starting_hour = case_when(
      str_detect(phi_type, "begins6") ~ 6,
      str_detect(phi_type, "begins12") ~ 12,
      TRUE ~ 0
    ),
    starting_hour_factor = case_when(
      str_detect(phi_type, "begins6") ~ "Starting at 6:00",
      str_detect(phi_type, "begins12") ~ "Starting at 12:00",
      TRUE ~ "Starting at 00:00"
    ),
    measure = case_when(
      str_detect(phi_type, "_mean$") ~ "Mean",
      str_detect(phi_type, "_median$") ~ "Median",
      str_detect(phi_type, "_ci_width$") ~ "95%CI width",
      TRUE ~ NA_character_
    ),
    # Streamlined psi calculation:
    psi = ifelse(
      measure == "95%CI width",
      phi * 12 / pi,
      ((phi * 12 / pi) + starting_hour) %% 24
    )
  ) %>%
  filter(!is.na(measure)) %>%
  select(
    item,
    # id,
    method,
    starting_hour,
    starting_hour_factor,
    measure,
    phi,
    psi
  )
