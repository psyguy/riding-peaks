library(dplyr)
library(tidyr)
library(stringr)

el1.5_long <- el1.5e_level1.5 %>%
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

el1_long %>%
  mutate(measure = factor(measure,
                          levels = c("Mean",
                                     "Median",
                                     "95%CI width")),
         item = case_when(
           item == "pa" ~ "Dataset 1",
           TRUE ~ "Dataset 2"
         )
         ) %>%
  ggplot() +
  aes(x = psi,
      fill = method) +
  geom_histogram(
    aes(y = ..ncount..),
    bins = 24,
    position = "identity",
    breaks = seq(0, 24, by = 0.5),
    alpha = 0.6
  ) +
  facet_nested(starting_hour_factor ~ item + measure,
               scales = "fixed") +
  xlab("Hours") +
  scale_x_continuous(breaks = (0:4) * 6) +
  scale_fill_manual(values = c("cornflowerblue", "brown1")) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.spacing.x = unit(1.5, "lines"),
    ggh4x.facet.nestline = element_line(linewidth = 0.5),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
  legend.text = element_text(size = 8),
  strip.text = element_text(size = 6),
  axis.title.x = element_text(size = 8),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank()
  )
