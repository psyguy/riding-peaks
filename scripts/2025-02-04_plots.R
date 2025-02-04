library(ggrepel)

## Reading and combining level-1 and level-2 plots

el1 <- e_level1 <-
  rbind(
    readRDS("fits/ests_level1_pa_random_var.rds"),
    # readRDS("fits/ests_level1_na_random_var.rds"),
    # readRDS("fits/ests_level1_hap_random_var.rds"),
    # readRDS("fits/ests_level1_sad_random_var.rds"),
    # readRDS("fits/ests_level1_fitbit.cov_random_var.rds),
    readRDS("fits/ests_level1_fitbit_random_var.rds")
    )

el2 <- e_level2 <-
  rbind(readRDS("fits/ests_level2_pa_random_var.rds"),
        # readRDS("fits/ests_level2_na_random_var.rds"),
        # readRDS("fits/ests_level2_hap_random_var.rds"),
        # readRDS("fits/ests_level2_sad_random_var.rds"),
        # readRDS("fits/ests_level2_fitbit.cov_random_var.rds"),
        readRDS("fits/ests_level2_fitbit_random_var.rds")
        )

# Selecting only pa, fitbit, fitbit.cov

el1 <- e_level1 %>%
  filter(item %in% c("pa", "fitbit", "fitbit.cov"))

el2 <- e_level2 %>%
  filter(item %in% c("pa", "fitbit", "fitbit.cov"))

# Level 1 -----------------------------------------------------------------

# Filter and reshape the data
ee1 <- el1 %>%
  select(contains("phi"), id, item, amp_hdr) %>%
  mutate(amp_hdr = amp_hdr/100) %>%
  #matches("phi_|phi\\.begins6_|phi\\.begins12_")) %>%  # Select relevant columns
  pivot_longer(
    cols = -item,
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
      grepl("hdr", variable) ~ "HDR%",
      TRUE ~ NA
    ) %>%
      factor(levels = c("mean", "median", "95%CI width", "HDR%"),
             labels = c("Person mean",
                        "Person median",
                        "Person 95%CIW",
                        "Person HDR%")),
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
  na.omit()


  # Create the ggplot
ee1 %>%
  # filter(item != "fitbit") %>%
  ggplot(aes(x = value)) +
  geom_histogram(aes(fill = lin_or_circ) , bins = 48,
                 alpha = 0.7,
                 position="identity") +
  facet_nested(variable ~ item + stat,
               scales = "free") +
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


# Create the ggplot for HDR% vs psi, for each individual
el1 %>%
  # filter(item != "fitbit") %>%
  ggplot(aes(y = amp_hdr,
             x = phi_circ_mean*24/(2*pi))) +
  geom_point() +
  facet_grid(cols = vars(item)) +
  # scale_x_continuous(breaks = (0:10) * 10) +
  geom_hline(yintercept = 95,
             color = "cornflowerblue") +
  geom_hline(yintercept = 90,
             color = "brown1") +
  geom_hline(yintercept = 85,
             color = "brown2") +
  geom_hline(yintercept = 80,
             color = "brown3") +
  scale_x_continuous(breaks = (0:6) * 4) +
  scale_fill_manual(values = c("brown1", "cornflowerblue")) +
  # ggthemes::theme_tufte() +
  theme_minimal() +
  labs(fill = NULL)


# Level 2 -----------------------------------------------------------------


ee2 <-
  el2 %>%
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
      par1 == "mesor" ~ "cor(.,MESOR)",
      par1 == "amp" ~ "cor(.,Amplitude)",
      par1 == "sigma2" ~ "cor(.,random variance)",
      par1 == "bl_pa" ~ "cor(.,baseline PA)",
      par1 == "bl_na" ~ "cor(.,baseline NA)"
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
  group_by(item, correlation_type, par1, measure,variable) %>%
  mutate(mm = median (value))


ee2 %>%
  filter(correlation_type %in% c("Pearson (lin-lin)",
                                 "Johnson–Wehrly–Mardia (circ-lin)",
                                 "Mardia (circ-lin rank)")) %>%
  # filter(item == "fitbit") %>%
  filter(#measure != "cor_abs",
         measure == "pval") %>%
  filter(par1 %in% c("cor(.,MESOR)","cor(.,Amplitude)")) %>%
  ggplot(aes(x = value, fill = correlation_type)) +
  geom_vline(xintercept = 0.05, linetype = "dashed") +
  geom_histogram(bins = 200,
                 alpha = 0.7,
                 position="identity") +
  # facet_nested(item + variable ~ par1 + measure,
  facet_nested(par1 + variable ~ item + measure,
               scales = "free") +
  geom_vline(aes(xintercept = mm,
                 color = correlation_type)) +
  labs(title = "Histograms circular-linear correlations",
       x = "Value",
       y = "Frequency") +
  scale_fill_manual(values = c("brown1",
                               # "brown3",
                               # "brown4",
                               "lightskyblue",
                               "cornflowerblue")) +
  scale_color_manual(values = c("brown1",
                                # "brown3",
                                # "brown4",
                                "lightskyblue",
                                "cornflowerblue"),
                     guide = "none") +
  ggthemes::theme_tufte() +
  # theme_minimal() +
  scale_y_continuous(transform = "sqrt") +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    strip.text = element_text(size = 16) # Increase panel label size
  ) +
  labs(fill = NULL)

# Scatter plots of cor and p-value


ee2 %>%
  filter(correlation_type %in% c("Pearson (lin-lin)",
                                 "Mardia (circ-lin rank)")) %>%
  # filter(item == "fitbit",
  #        iteration < 1000) %>%
  filter(measure != "cor_abs") %>%
  select(-mm) %>%
  pivot_wider(names_from = measure,
              values_from = value) %>%
  ggplot(aes(x = pval, y = cor, color = correlation_type)) +
  facet_nested(par1 + variable ~ item,# + measure,
               scales = "free") +
  geom_point() +
  geom_vline(aes(xintercept = c(0.05))) +
  labs(title = "Histograms circular-linear correlations",
       x = "p-value",
       y = "Correlation") +
  geom_label_repel(
    data = . %>%
      group_by(item, variable, correlation_type, par1) %>%
      summarise(percent_significant = paste0(
        round(100*sum(pval <= 0.05)/6000),
        "% significant")),
    aes(x = 0.5, y = 0.5,
        label = percent_significant,
        color = correlation_type),
    inherit.aes = FALSE,
    box.padding = 0.5,   # Add space between text and the plot area
    point.padding = 0.5, # Add space between text and data points
    max.overlaps = 20,   # Limit number of overlaps to avoid crowded areas
    direction = "y",  # Allow text to move in both x and y directions
    nudge_y = 0.05,     # Nudge labels slightly for better spacing
    nudge_x = 0.05,
    min.segment.length = 100
  ) +
  scale_fill_manual(values = c("brown1",
                               # "brown3",
                               # "brown4",
                               # "lightskyblue",
                               "cornflowerblue")) +
  scale_color_manual(values = c("brown1",
                                # "brown3",
                                # "brown4",
                                # "lightskyblue",
                                "cornflowerblue")) +
  ggthemes::theme_tufte() +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    # axis.text.y = element_blank(),
    # axis.ticks.y = element_blank(),
    # axis.title.y = element_blank(),
    strip.text = element_text(size = 16) # Increase panel label size
  ) +
  labs(fill = NULL)

