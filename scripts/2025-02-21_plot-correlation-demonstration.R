# Data is from the ozone dataset, but shifted a bit

dw <- data.frame(
  wind_degree = c(
    327, 91, 88, 305, 344, 270, 67, 21, 281, 8,
    204, 86, 333, 18, 57, 6, 11, 27, 84),
  ozone = c(
    28.0, 85.2, 80.5, 4.7, 45.9, 12.7, 72.5, 56.6, 31.5, 112.0,
    20.0, 72.5, 16.0, 45.9, 32.6, 56.6, 52.6, 91.8, 55.2)
) %>%
  mutate(phi = (wind_degree * pi/180 + 1.5) %% (2*pi),
         z = ozone) %>%
  mutate(r_phi = rank(phi),
         r_z = rank(z)) %>%
  mutate(ang = 2*pi*r_phi/19,
         x = cos(ang),
         y = sin(ang))

fit_lin <- lm(z ~ phi, data = dw)
fit_cos <- lm(z ~ cos(phi) + sin(phi), data = dw)

pp <- seq(0, 2*pi, length.out = 200)
curves <- data.frame(
  phi = pp,
  z_lin = predict.lm(fit_lin, newdata = data.frame(phi = pp)),
  z_cos = predict.lm(fit_cos, newdata = data.frame(phi = pp))
) #
# %>%
#   pivot_longer(z_lin:z_cos,
#                names_to = "model",
#                values_to = "z") %>%
#   mutate(color = case_when(
#     model == "z_lin" ~ "blue",
#     TRUE ~ "red"
#    ))


dw %>%
  ggplot() +
  aes(x = x, y = y,
      color = r_phi,
      fill = r_phi,
      size = r_z) +
  # Outer circle
  geom_circle(aes(x0 = 0, y0 = 0, r = 1),
              color = "gray70",
              inherit.aes = FALSE) +
  geom_line(aes(x = phi, y = z_lin),
            data = curves,
            color = "blue",
            show.legend = FALSE,
            linetype = "dashed",
            inherit.aes = FALSE) +
  geom_line(aes(x = phi, y = z_cos),
            data = curves,
            color = "red",
            show.legend = FALSE,
            linetype = "dashed",
            inherit.aes = FALSE) +
  geom_point(color = "black",
             stroke = 0.75,
             shape = 21) +
  # lollipops
  geom_point(aes(x = phi, y = z, color = r_phi),
             size = 2) +
  geom_segment(aes(x = phi, xend = phi,
                   y = 0,
                   yend = z, color = r_phi),
               size = 0.5) +
  # polar plot
  geom_point(aes(x = z*cos(phi), y = z*sin(phi), color = r_phi),
             size = 2) +
  geom_segment(aes(x = 0, xend = z*cos(phi),
                   y = 0, yend = z*sin(phi), color = r_phi),
               size = .5) +
  # coord_fixed(xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1)) +
  scale_fill_viridis(end = 0.9) +
  scale_color_viridis(end = 0.9) +
  theme_clean() +
  # theme(
  # #  legend.position = "none",
  #   panel.grid = element_blank(),
  #   axis.text  = element_blank(),
  #   axis.ticks = element_blank(),
  #   axis.title = element_blank()
  # ) +
  ggtitle("Meh")
