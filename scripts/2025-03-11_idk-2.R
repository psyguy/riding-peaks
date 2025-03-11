el2 %>%
  filter(grepl("phi", par2)) %>%
  # filter(!grepl("log", par1)) %>%
  pivot_wider(names_from = measure,
              values_from = value) %>%
  mutate(cor_abs = abs(cor)) %>%
  pivot_longer(cor:cor_abs,
               names_to = "measure",
               values_to = "value") %>%
  mutate(
    variable = case_when(
      grepl("6", par2) ~ "Starting at 6:00",
      grepl("12", par2) ~ "Starting at 12:00",
      TRUE ~ "Starting at 0:00"
    ) %>%
      factor(levels = c("Starting at 0:00",
                        "Starting at 6:00",
                        "Starting at 12:00")),,
    correlation_type =
      factor(
        correlation_type,
        levels = c("pearson", "spearman", "kendall", "jwm", "mard"),
        labels = c(
          "Pearson",
          "Spearman",
          "Kendall",
          "Johnson–Wehrly–Mardia",
          "Mardia"
        )
      )
  ) %>%
  group_by(item, correlation_type, par1, measure, variable) %>%
  filter(
    !grepl("Kendall", correlation_type, ignore.case = TRUE),
    measure == "cor",
    par1 %in% c("amp", "mesor", "logsigma")
  ) %>%
  mutate(mm = median (value),
         value = case_when(
           correlation_type == "Johnson–Wehrly–Mardia" ~ sqrt(value),
           TRUE ~ value
         ),
         item = case_when(item == "pa" ~ "Dataset 1",
                          TRUE ~ "Dataset 2"),
         par1 = case_when(par1 == "amp" ~ "Amplitude",
                          par1 == "mesor" ~ "MESOR",
                          par1 == "logsigma" ~ "log(SD)",
                          TRUE ~ par1) %>%
           factor(levels = c("Amplitude", "MESOR", "log(SD)"))) -> ee
