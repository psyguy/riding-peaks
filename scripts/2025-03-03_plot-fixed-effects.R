el1.5 %>%
  # select(!contains(c("co_",
  #                    "si_",
  #                    "ci_",
  #                    "hdr",
  #                    "mean"))) %>%
  select(item,
         iteration,
         contains(c("mesor",
                    "amp",
                    "logsigma",
                    "phi")),
         -amp_hdr) %>%
  select(!contains("ci_")) %>%
  pivot_longer(
    cols = matches("^phi") | matches("mean") | matches("median"),
    names_to = "par",
    values_to = "value"
  ) %>%
  mutate(
    based_on = case_when(
      str_detect(par, "mean") ~ "Mean",
      str_detect(par, "median") ~ "Median",
      TRUE ~ par
    ),
    method = ifelse(str_detect(par, "_circ"),
                    "circular",
                    "linear"),
    starting_hour = case_when(
      str_detect(par, "begins6") ~ 6,
      str_detect(par, "begins12") ~ 12,
      TRUE ~ 0
    ),
    starting_hour_factor = case_when(
      str_detect(par, "begins6") ~ "Starting at 6:00",
      str_detect(par, "begins12") ~ "Starting at 12:00",
      TRUE ~ "Starting at 00:00"
    ),
    par = case_when(
      str_detect(par, "mesor") ~ "MESOR",
      str_detect(par, "amp") ~ "Amplitude",
      str_detect(par, "logsigma") ~ "log(SD)",
      str_detect(par, "phi") ~ "phi",
      TRUE ~ par
    )
  ) %>%
  filter() -> ee

ee %>%
  filter(par != "phi") %>%
ggplot(aes(x = value)) +
  geom_histogram(aes(y = ..ncount..),
                 alpha = 0.6,
                 bins = 300,
                 fill = "chartreuse3") +
  theme_minimal() +
  labs(x = "", y = "") +
  scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) +
  facet_nested(based_on ~ item + par, scales = "free",
               solo_line = TRUE) +
  # Add line under dataset name in the facet title
  annotate("segment",
           x = Inf,
           xend = -Inf,
           y = Inf,
           yend = Inf,
           color = "black",
           lwd = 0.75) +
  theme(panel.spacing = unit(3, "lines"),
        strip.text = element_text(size = 12),
        legend.position = "none",
        axis.text = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank(),
        plot.margin = margin(0,0,0,0))


# Plot four circles -------------------------------------------------------




