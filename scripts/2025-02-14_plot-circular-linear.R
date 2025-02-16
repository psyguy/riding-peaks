library(ggforce)

phi <- draws %>%
  filter(item == "pa", id == 26) %>%
  pull(phi)

f_plot_cirular_linear <- function(phi) {

# one bin for every minutes
bins <- 24 * 60 / 1

arc <- (2 * pi) / bins
breaks <- seq(0, 2 * pi, length.out = (bins + 1))
bins.count <- hist.default(phi,
                           breaks = breaks,
                           plot = FALSE,
                           right = TRUE)$counts
mids <- seq(arc / 2, 2 * pi - pi / bins, length = bins)

inc <- .2 / max(bins.count)

d <- data.frame(r = rep(NA, 1),
                x = rep(NA, 1),
                y = rep(NA, 1))
d_tmp <- d

for (i in 1:bins) {
  if (bins.count[i] != 0) {
    for (j in 0:(bins.count[i] - 1)) {
      d_tmp$r <- 1 + j * inc
      d_tmp$x <- d_tmp$r * cos(mids[i])
      d_tmp$y <- d_tmp$r * sin(mids[i])
      d <- rbind(d, d_tmp)
    }
  }
}

# Consolidated Data Construction
d_summary <- data.frame(
  type = c("linear", "circular"),
  mean = c(phi %>% mean(), phi %>% circular %>% mean() %>% `%%`(2 * pi)),
  median = c(phi %>% median(), phi %>% circular %>% median() %>% `%%`(2 *
                                                                        pi)),
  ci_2.5 = c(phi %>% quantile(0.025), phi %>% circular %>% quantile(0.025)),
  ci_97.5 = c(phi %>% quantile(0.975), phi %>% circular %>% quantile(0.975)),
  r = c(0.7, 0.9)
) %>%
  mutate(ci_2.5 = if_else(ci_2.5 > ci_97.5, ci_2.5 - 2 * pi, ci_2.5))

d_estimates <- d_summary %>%
  pivot_longer(
    cols = c(mean, median),
    names_to = "measure",
    values_to = "value"
  ) %>%
  mutate(x = cos(value) * 1.4, y = sin(value) * 1.4)

p_circle <- ggplot() +
  # Background points
  geom_point(data = d,
             aes(x = x, y = y),
             size = 0.5,
             color = "gray50") +

  # Cross lines
  geom_vline(xintercept = 0,
             linewidth = 0.3,
             color = "gray60") +
  geom_hline(yintercept = 0,
             linewidth = 0.3,
             color = "gray60") +

  # Outer circle
  geom_circle(
    aes(x0 = 0, y0 = 0, r = 1),
    data = data.frame(),
    inherit.aes = FALSE,
    linetype = "dashed",
    color = "gray70"
  ) +

  # Confidence Bands as Slices
  geom_polygon(
    data = d_summary %>%
      group_by(type) %>%
      summarise(x = c(0, cos(
        seq(ci_2.5, ci_97.5, length.out = 100)
      ) * r), y = c(0, sin(
        seq(ci_2.5, ci_97.5, length.out = 100)
      ) * r)),
    aes(x = x, y = y, fill = type),
    alpha = 0.5
  ) +

  # Rays for mean/median
  geom_segment(
    aes(
      x = 0,
      y = 0,
      xend = cos(value) * r,
      yend = sin(value) * r,
      color = type,
      linetype = measure
    ),
    data = d_summary %>%
      pivot_longer(
        cols = c(mean, median),
        names_to = "measure",
        values_to = "value"
      ),
    linewidth = 1.2,
    alpha = 0.7
  ) +

  # Asterisk Markers
  geom_point(
    aes(
      x = cos(value) * r,
      y = sin(value) * r,
      shape = measure,
      color = type
    ),
    size = 5,
    stroke = 1.2,
    data = d_summary %>%
      pivot_longer(
        cols = c(mean, median),
        names_to = "measure",
        values_to = "value"
      )
  ) +

  # Color, Shape, and Line Scales
  scale_fill_manual(values = c("circular" = "#377eb8",
                               "linear" = "#e41a1c")) +
  scale_color_manual(values = c("circular" = "#377eb8",
                                "linear" = "#e41a1c")) +
  scale_shape_manual(values = c("mean" = 1,
                                "median" = 8)) +
  scale_linetype_manual(values = c("mean" = "solid",
                                   "median" = "dashed")) +

  # Axis and Layout
  coord_fixed(xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5)) +

  # Theme
  theme_minimal() +
  theme(
    legend.position = "none",
    # legend.title = element_blank(),
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank()
  )


p_hist <-
  data.frame(x = phi) %>%
  ggplot() +
  aes(x = x) +
  geom_histogram(aes(y = ..ncount..),
                 bins = 24 * 10) +
  geom_segment(
    aes(
      x = value,
      xend = value,
      y = 0,
      yend = r * 0.9,
      color = type,
      linetype = measure
    ),
    lwd = 1.5,
    alpha = 0.6,
    data = d_estimates,
    inherit.aes = FALSE
  ) +
  geom_point(
    aes(
      x = value,
      y = r * 0.9,
      shape = measure,
      color = type
    ),
    size = 5,
    stroke = 1.2,
    data = d_estimates,
    inherit.aes = FALSE
  ) +
  geom_rect(aes(xmin = ci_2.5,
                xmax = ci_97.5,
                ymin = 0,
                ymax = r * 0.9,
                fill = type),
            data = d_estimates %>%
              mutate(
                ci_2.5= ci_2.5 %% (2*pi),
                ci_97.5 = ci_97.5 %% (2*pi),
                new_ci_2.5 = case_when(
                  measure == "mean" & ci_2.5 > ci_97.5 ~ 0,
                  measure == "median" & ci_2.5 > ci_97.5 ~ ci_2.5,
                  TRUE ~ ci_2.5
                ),
                new_ci_97.5 = case_when(
                  measure == "mean" & ci_2.5 > ci_97.5 ~ ci_97.5,
                  measure == "median" & ci_2.5 > ci_97.5 ~ 2 * pi,
                  TRUE ~ ci_97.5
                )
              ) %>%
              select(-ci_2.5, -ci_97.5) %>%
              rename(ci_2.5 = new_ci_2.5, ci_97.5 = new_ci_97.5) %>%
              select(type, ci_2.5, ci_97.5, r) %>%
              distinct(),
            alpha = 0.4,
            inherit.aes = FALSE
                  ) +
  theme_minimal() +
  theme(#legend.position = "none",
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank()) +
  xlab("Position in cycle") +
  ylab(NULL) +
  # Color, Shape, and Line Scales
  scale_fill_manual(values = c("circular" = "#377eb8",
                               "linear" = "#e41a1c")) +
  scale_color_manual(values = c("circular" = "#377eb8",
                                "linear" = "#e41a1c")) +
  scale_shape_manual(values = c("mean" = 1,
                                "median" = 8)) +
  scale_linetype_manual(values = c("mean" = "solid",
                                   "median" = "dashed")) +
    expand_limits(x = c(0, 2 * pi))

(p_circle / p_hist) +
  plot_layout(heights = c(1, 0.5))

}
