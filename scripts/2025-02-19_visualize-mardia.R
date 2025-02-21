# Mardia circular-linear rank correlation
cor_mardia <- function(theta, x) {
  # Remove NAs
  valid <- complete.cases(theta, x)
  theta <- theta[valid]
  x <- x[valid]
  n <- length(theta)
  if (n == 0) stop("No valid data after removing NAs.")

  # Rank and compute theta star
  r_theta <- rank(theta, ties.method = "average")
  r_x <- rank(x, ties.method = "average")
  r_theta_star <- r_theta * 2 * pi / n

  # Precompute sine and cosine
  cos_theta_star <- cos(r_theta_star)
  sin_theta_star <- sin(r_theta_star)

  # Compute T_c and T_s
  T_c <- sum(r_x * cos_theta_star)
  T_s <- sum(r_x * sin_theta_star)

  # Calculate coefficient 'a'
  a <- ifelse(
    n %% 2 == 0,
    1 / (1 + 5 / tan(pi / n)^2 + 4 / tan(pi / n)^4),
    2 * sin(pi / n)^4 / (1 + cos(pi / n))^3
  )

  # Final D value (rank correlation, values between 0 and 1)
  D <- a * (T_c^2 + T_s^2)

  # Calculate U_n and p-value
  U_n <- 24 * (T_c^2 + T_s^2) / (n^2 * (n + 1))
  p.value <- 1 - pchisq(U_n, df = 2)

  # Return results as a list
  return(list(estimate = D, p.value = p.value, U_n = U_n))
}

df <- draws %>%
  filter(item == "pa",
         iteration == 1) %>%
  rename(x = amp,
         theta = phi) %>%
  select(id, theta, x) %>%
  mutate(r_theta = rank(theta),
         r_x = rank(x))

n <- nrow(df)

r_theta <- df$r_theta
r_x <- df$r_theta
r_theta_star <- r_theta * 2 * pi / n

phi <- rep(r_theta_star, r_x)

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


ggplot() +
  # Cross lines
  geom_vline(xintercept = 0, linewidth = 0.3, color = "gray60") +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "gray60") +
  # Outer circle
  geom_circle(aes(x0 = 0, y0 = 0, r = 1), color = "gray70", inherit.aes = FALSE) +
  # Background points
  geom_point(data = d, aes(x = x, y = y), size = 0.2, color = "azure4") +
  # Rays for mean/median
  geom_segment(
    data = d_summary %>% pivot_longer(cols = c(mean, median), names_to = "measure", values_to = "value"),
    aes(x = 0, y = 0, xend = cos(value) * r, yend = sin(value) * r,
        color = type, linetype = measure),
    linewidth = 1.2, alpha = 0.7
  ) +
  # Asterisk markers
  geom_point(
    data = d_summary %>% pivot_longer(cols = c(mean, median), names_to = "measure", values_to = "value"),
    aes(x = cos(value) * r, y = sin(value) * r, shape = measure, color = type),
    size = 3, stroke = 1
  ) +
  coord_fixed(xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1)) +
  theme_minimal()
