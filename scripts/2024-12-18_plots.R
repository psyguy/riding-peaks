# library(ggh4x)

el1 <-
  ests_level1 <-
  readRDS(here::here("fits/ests_level1_hap_random_var.rds"))
el2 <-
  ests_level2 <-
  readRDS(here::here("fits/ests_level2_hap_random_var.rds"))

## Level-1 plots

names(el1)

el1 <- ests_level1 %>%
  select(!contains("log")) %>%
  group_by(id) %>%
  pivot_longer(-id,
               names_to = "what",
               values_to = "value") %>%
  mutate(what = gsub("circ_", "circ.", what) %>%
           gsub("ci_", "ci.", .)) %>%
  group_by(id) %>%
  separate(what,
           c("par", "stat"),
           sep = "_") %>%
  group_by(par) %>%
  pivot_wider(names_from = stat,
              values_from = value)

el1 %>%
  filter(grepl("phi", par)) %>%
  ggplot() +
  aes(x = mean,
      y = circ.mean) +
  geom_point() +
  facet_grid(rows = ~par)


# Level 1 -----------------------------------------------------------------

# Filter and reshape the data
ests_level1 %>%
  select(contains("phi")) %>%
  #matches("phi_|phi\\.begins6_|phi\\.begins12_")) %>%  # Select relevant columns
  pivot_longer(
    cols = everything(),
    names_to = "variable",
    values_to = "value"
  ) %>%
  group_by(variable) %>%
  mutate(
    # Identify whether the variable is _mean or _circ_mean
    lin_or_circ = case_when(
      grepl("_circ", variable) ~ "circ",
      TRUE ~ "lin"
    ) %>%
      factor(levels = c("lin", "circ"),
             labels = c("Treated as linear variable",
                        "Treated as circular variable")),
    stat = case_when(
      grepl("mean", variable) ~ "mean",
      grepl("median", variable) ~ "median",
      grepl("width", variable) ~ "95%CI width",
      TRUE ~ NA
    ) %>%
      factor(levels = c("mean", "median", "95%CI width"),
             labels = c("Person mean",
                        "Person median",
                        "Person 95%CI width")),
    # variable = stringr::word(variable, sep = "_")[1],
    variable = case_when(
      grepl("6", variable) ~ "Starting at 06:00",
      grepl("12", variable) ~ "Starting at 12:00",
      TRUE ~ "Starting at 00:00"
    ),
    value = case_when(
      grepl("6", variable) & stat != "Person 95%CI width" ~ (value + pi/2) %% (2*pi),
      grepl("12", variable) & stat != "Person 95%CI width" ~ (value + pi) %% (2*pi),
      TRUE ~ value
    ) * 24/(2*pi)
  ) %>%
  # Create the ggplot
  ggplot(aes(x = value)) +
  geom_histogram(aes(fill = lin_or_circ) , bins = 48,
                 alpha = 0.7,
                 position="identity") +
  facet_nested(variable ~ stat) +
  labs(title = "Estiamed person-specific peak shifts",
       x = "Hour",
       y = "Frequency") +
  scale_x_continuous(breaks = (0:6) * 4) +
  scale_fill_manual(values = c("brown1", "cornflowerblue")) +
  ggthemes::theme_tufte() +
  # theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    strip.text = element_text(size = 16) # Increase panel label size
  ) +
  labs(fill = NULL)


# Level 2 -----------------------------------------------------------------


ests_level2 %>%
  filter(grepl("phi", par2)) %>%
  filter(!grepl("log", par1)) %>%
  pivot_wider(names_from = measure,
              values_from = value) %>%
  mutate(cor_abs = abs(cor)) %>%
  pivot_longer(cor:cor_abs,
               names_to = "measure",
               values_to = "value") %>%
  mutate(
    variable = case_when(
      grepl("6", par2) ~ "Starting at 06:00",
      grepl("12", par2) ~ "Starting at 12:00",
      TRUE ~ "Starting at 00:00"
    ),
    par1 = case_when(
      par1 == "mesor" ~ "Correlated with MASOR",
      par1 == "amp" ~ "Correlated with Amplitude",
      par1 == "sigma2" ~ "Correlated with random variance",
    ),
    correlation_type =
      factor(
        correlation_type,
        levels = c("pearson", "spearman", "kendall", "jwm", "mard"),
        labels = c(
          "Pearson (lin-lin)",
          "Spearman (lin-lin rank)",
          "Kendall (lin-lin rank)",
          "Johnson–Wehrly–Mardia (circ-lin)",
          "Mardia (circ-lin rank)"
        )
      )
  ) %>%
  # filter(correlation_type %in% c("Pearson (lin-lin)", "Mardia (circ-lin rank)")) %>%
  group_by(correlation_type, par1, measure,variable) %>%
  mutate(mm = median (value)) %>%
  ggplot(aes(x = value, fill = correlation_type)) +
  geom_histogram(bins = 50,
                 alpha = 0.7,
                 position="identity") +
  facet_nested(variable ~ par1 + measure,
               scales = "fixed") +
  geom_vline(aes(xintercept = mm,
                 color = correlation_type)) +
  labs(title = "Histograms circular-linear correlations",
       x = "Value",
       y = "Frequency") +
  scale_fill_manual(values = c("brown1",
                               "brown3",
                               "brown4",
                               "lightskyblue",
                               "cornflowerblue")) +
  scale_color_manual(values = c("brown1",
                                "brown3",
                                "brown4",
                                "lightskyblue",
                                "cornflowerblue"),
                     guide = "none") +
  ggthemes::theme_tufte() +
  # theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    strip.text = element_text(size = 16) # Increase panel label size
  ) +
  labs(fill = NULL)

