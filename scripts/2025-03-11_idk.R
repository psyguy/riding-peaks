dd <- draws %>%
  filter(iteration == 12, item == "PA") %>%
  mutate(z = amp,
         x = phi)

r_p <- dd$phi %>% rank()
r_z <- dd$amp %>% rank()
phi <- (2*pi*r_p/nrow(dd)) %>% rep(times = r_z)
title <- "meh"
# Set up bins: one per minute in 24 hours (1440 bins)
bins <- 24 * 60
arc <- 2 * pi / bins
breaks <- seq(0, 2 * pi, length.out = bins + 1)
bins.count <- hist.default(phi, breaks = breaks, plot = FALSE, right = TRUE)$counts
mids <- seq(arc / 2, 2 * pi - arc / 2, length.out = bins)
inc <- 0.2 / max(bins.count)
# Create background points (circular histogram) using a vectorised approach
d <- do.call(rbind, lapply(seq_len(bins), function(i) {
  count <- bins.count[i]
  if (count == 0) return(NULL)
  j <- seq(0, count - 1)
  r <- 1 + j * inc
  data.frame(r = r, x = r * cos(mids[i]), y = r * sin(mids[i]))
}))
# Summarise estimates for both linear and circular methods
d_summary <- data.frame(
  type    = c("linear", "circular"),
  mean    = c(mean(phi), (mean(circular(phi)) %% (2 * pi))),
  median  = c(median(phi), (median(circular(phi)) %% (2 * pi))),
  ci_2.5  = c(quantile(phi, 0.025), quantile(circular(phi), 0.025)),
  ci_97.5 = c(quantile(phi, 0.975), quantile(circular(phi), 0.975)),
  r       = c(0.85, 1)
) %>%
  mutate(ci_2.5 = if_else(ci_2.5 > ci_97.5, ci_2.5 - 2 * pi, ci_2.5))
d_estimates <- d_summary %>%
  pivot_longer(cols = c(mean, median), names_to = "measure", values_to = "value") %>%
  mutate(x = cos(value) * 1.4, y = sin(value) * 1.4)
# Circular plot


  ggplot() +
  # Cross lines
  geom_vline(xintercept = 0, linewidth = 0.3, color = "gray60") +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "gray60") +
  # Outer circle
  geom_circle(aes(x0 = 0, y0 = 0, r = 1), color = "gray70", inherit.aes = FALSE) +
  # Background points
  geom_point(data = d, aes(x = x, y = y), size = 0.2, color = "azure4") +
  coord_fixed(1) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.text  = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank()
  )put change the layout a bit:
