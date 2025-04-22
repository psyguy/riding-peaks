library(shiny)
library(dplyr)
library(ggplot2)
library(ggforce)       # for geom_circle
library(circular)      # for circular::mean.circular

point_size <- 1
main_title <- "Raw data"

f_plot_correlation_demonstration <- function(ed,
                                             main_title = "",
                                             # subtitle = "",
                                             point_size = 1){

# -------------- Main Plot --------------

# Linear: z ~ phi
fit_lin <- lm(z ~ phi, data = ed)
phi_grid <- seq(0, 2 * pi, length.out = 200)
df_lin   <- data.frame(phi = phi_grid)
df_lin$z <- predict(fit_lin, newdata = df_lin)

# Cosinor: z ~ cos(phi) + sin(phi)
fit_cos <- lm(z ~ cos(phi) + sin(phi), data = ed)
df_cos  <- data.frame(phi = phi_grid)
df_cos$z <- predict(fit_cos, newdata = df_cos)
# MESOR line
M <- coef(fit_cos)[1]

corr <- suppressWarnings(compute_correlations(ed$phi, ed$z))
lbl <- paste(
  sprintf("Pearson:     %.2f (p=%.2f)", corr$pearson_cor, corr$pearson_pval),
  sprintf("Spearman:  %.2f (p=%.2f)",
    corr$spearman_cor,
    corr$spearman_pval
  ),
  sprintf("JWM:           %.2f (p=%.2f)", corr$jwm_cor, corr$jwm_pval),
  sprintf("Mardia:        %.2f (p=%.2f)", corr$mard_cor, corr$mard_pval),
  sep = "\n"
)


p_main <- ggplot(ed, aes(x = phi, y = z)) +
  geom_point(color = "azure4", size = point_size) +
  labs(x = TeX("Circular variable \\phi"),
       y = "Linear variable Z",
       title = main_title) +
  geom_line(
    data = df_lin,
    aes(phi, z),
    color = "#9C179EFF",
    size = 1
  ) +
  geom_line(
    data = df_cos,
    aes(phi, z),
    color = "#54C568FF",
    size = 1
  ) +
  geom_hline(yintercept = M,
             color = "#A5DB36FF",
             linetype = "dashed") +
  # ylim(-26, +44) +
  annotate(
    "label",
    x = -.1,
    y = 1.1*max(ed$z),
    label = lbl,
    fill = "white",
    alpha = 0.3,
    hjust = 0,
    vjust = 1,
    size = 4
  ) +
  scale_x_continuous(
    breaks = c(0, pi/2, pi, 3*pi/2, 2*pi),
    labels = c("0", expression(pi/2), expression(pi), expression(3*pi/2), expression(2*pi))
  ) +
  theme_few() +
  theme(plot.title = element_text(size = 16))
# -------------- Polar Plot --------------

df_cart <- ed %>% mutate(x_eff = z * cos(phi), y_eff = z * sin(phi))

lim_polar <- ed %>% abs() %>% max()

# Linear fit curve
phi_grid <- seq(0, 2 * pi, length.out = 200)
z_lin    <- predict(fit_lin, newdata = data.frame(phi = phi_grid))
df_lin   <- data.frame(x_eff = z_lin * cos(phi_grid),
                       y_eff = z_lin * sin(phi_grid))

# Cosinor fit curve
z_cos   <- predict(fit_cos, newdata = data.frame(phi = phi_grid))
df_cos  <- data.frame(x_eff = z_cos * cos(phi_grid),
                      y_eff = z_cos * sin(phi_grid))

p_polar <- ggplot(df_cart, aes(x = x_eff, y = y_eff)) +
  # xlim(-lim_max(), lim_max()) +
  # ylim(-lim_max(), lim_max()) +
  geom_hline(yintercept = 0,
             # linetype = "dashed",
             color = "gray") +
  geom_vline(xintercept = 0,
             # linetype = "dashed",
             color = "gray") +
  geom_point(size = point_size, color = "azure4") +
  geom_path(
    data = df_lin,
    aes(x_eff, y_eff),
    color = "#9C179EFF",
    size = 1
  ) +
  geom_path(
    data = df_cos,
    aes(x_eff, y_eff),
    color = "#54C568FF",
    size = 1
  ) +
  # annotate(
  #   "label",
  #   x = -0.95 * 10, #lim_max(),
  #   y = 0.95 * 10, #lim_max(),
  #   label = lbl,
  #   fill = "white",
  #   alpha = 0.5,
  #   hjust = 0,
  #   vjust = 1,
  #   size = 3
  # ) +
  # labs(title = "Polar Plot") +
coord_fixed(1,
            xlim = c(-lim_polar, lim_polar),
            ylim = c(-lim_polar, lim_polar)) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid      = element_blank(),
    axis.text       = element_blank(),
    axis.ticks      = element_blank(),
    axis.title      = element_blank()
  )

# -------------- Mardia Plot --------------

# Build histogram from rank-transform
r_phi  <- rank(ed$phi)
r_z  <- rank(ed$z)
phi_vec <- rep(2 * pi * r_phi / nrow(ed), times = r_z)

bins  <- min(length(r_phi), 24 * 60)
arc   <- 2 * pi / bins
brks  <- seq(0, 2 * pi, length.out = bins + 1)
h     <- hist.default(phi_vec,
                      breaks = brks,
                      plot = FALSE,
                      right = TRUE)
mids  <- seq(arc / 2, 2 * pi - arc / 2, length.out = bins)

bins_count <- h$counts
inc <- 0.5*point_size / max(bins_count, na.rm = TRUE)  # spacing outward

d <- do.call(rbind, lapply(seq_len(bins), function(i) {
  count <- bins_count[i]
  if (count == 0)
    return(NULL)
  j <- seq(0, count - 1)
  r <- 1 + j * inc
  data.frame(r = r,
             x = r * cos(mids[i]),
             y = r * sin(mids[i]))
}))

lim_mardia <- d %>% abs() %>% max()

T_c <- r_z * cos(r_phi * 2 * pi / nrow(ed))
T_s <- r_z * sin(r_phi * 2 * pi / nrow(ed))

resultant_l <- cor_mardia(ed$phi, ed$z)$estimate
resultant_angle <- atan2(sum(T_s), sum(T_c))

p_mardia <- ggplot() +
  geom_vline(xintercept = 0,
             linewidth = 0.3,
             color = "gray60") +
  geom_hline(yintercept = 0,
             linewidth = 0.3,
             color = "gray60") +
  geom_circle(aes(x0 = 0, y0 = 0, r = 1),
              color = "gray70",
              inherit.aes = FALSE) +
  geom_point(
    data = d,
    aes(x = x, y = y),
    size = 0.5*point_size,
    color = "azure4"
  ) +
  geom_segment(aes(
    x = 0,
    y = 0,
    xend = resultant_l * cos(resultant_angle),
    yend = resultant_l * sin(resultant_angle)
  ),
  arrow = arrow(type = "closed",
                        length = unit(point_size*10, "points")),
  color = "cornflowerblue"
  ) +
  coord_fixed(1,
              xlim = c(-lim_mardia, lim_mardia),
              ylim = c(-lim_mardia, lim_mardia)) +
  geom_circle(aes(x0 = 0, y0 = 0, r = resultant_l),
              linetype = "dashed",
              color = "cornflowerblue", inherit.aes = FALSE) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid      = element_blank(),
    axis.text       = element_blank(),
    axis.ticks      = element_blank(),
    axis.title      = element_blank()
  )


# p_main + p_polar + p_mardia +
#   plot_layout(widths = c(2.5, 1.5, 1.5))

p_main + p_mardia +
  plot_layout(widths = c(2, 1.5))

}



set.seed(2025-03-07)
N <- 200
phi_sim <- runif(N, 0, 2*pi)
rn <- rnorm(N)
z_sim <- 10 + 8*cos(phi_sim) +
  (1.5*cos(2*phi_sim))^3 +
  cos(3*phi_sim) +
  5*rn
# phi_sim <- (phi_sim + 1.7*pi) %% (2*pi)


dd <- data.frame(phi = (phi_sim + 1.7*pi) %% (2*pi),
                 phi_sim = phi_sim,
                 z = 80*cos(phi_sim + pi/3*cos(phi_sim)) +
                   0*(1.5*cos(2*phi_sim))^3 +
                   0*cos(3*phi_sim) +
                   .5*rn,
                 rn = rn)

set.seed(2025-03-07)
N <- 300
rn <- rnorm(N)
phi_sim <- runif(N, 0, 2*pi)
data.frame(phi = (phi_sim + 1.04*pi) %% (2*pi),
           z = 10*cos(phi_sim + pi/3*cos(phi_sim)) +
             1*rn) %>%
  f_plot_correlation_demonstration()


p_1 <- f_plot_correlation_demonstration(
  data.frame(phi = (phi_sim + 0.54*pi) %% (2*pi),
             z = cos(phi_sim + pi/3*cos(phi_sim)) +  0.1*rn),
  "Raw simulated data")

p_2 <- f_plot_correlation_demonstration(
  data.frame(phi = (phi_sim + 0.04*pi) %% (2*pi),
             z = cos(phi_sim + pi/3*cos(phi_sim)) +  0.1*rn),
  expression(phi ~ "shifted by" ~ -pi/2))

p_3 <- f_plot_correlation_demonstration(
  data.frame(phi = (phi_sim + 1.04*pi) %% (2*pi),
             z = cos(phi_sim + pi/3*cos(phi_sim)) +  0.1*rn),
  expression(phi ~ "shifted by" ~ +pi/2))

p_4 <- f_plot_correlation_demonstration(
  data.frame(phi = (phi_sim + 0.54*pi) %% (2*pi),
             z = cos(phi_sim + pi/3*cos(phi_sim)) +  1*rn),
  "More noise added Z")

p_5 <- f_plot_correlation_demonstration(
  data.frame(phi = (phi_sim + 1.04*pi) %% (2*pi),
             z = cos(phi_sim + pi/3*cos(phi_sim)) +  1*rn),
  expression("More noise added to Z," ~ phi ~ "shifted by" ~ +pi/2))


(p_1 / plot_spacer() /
    p_2 / plot_spacer() /
    p_3 / plot_spacer() /
    p_4 / plot_spacer() /
    p_5 ) +
  plot_layout(heights = c(1, .05, 1, .05, 1, .05, 1, .05, 1))

