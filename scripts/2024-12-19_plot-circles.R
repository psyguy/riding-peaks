ee <-
  cbind(item = "hap",
        readRDS("fits/ests_level1_hap_random_var.rds")) %>%
  rbind(cbind(item = "sad",
              readRDS("fits/ests_level1_sad_random_var.rds"))) %>%
  rbind(cbind(item = "pa",
              readRDS("fits/ests_level1_pa_random_var.rds"))) %>%
  rbind(cbind(item = "na",
              readRDS("fits/ests_level1_na_random_var.rds"))) %>%
  rbind(cbind(item = "fitbit",
              readRDS("fits/ests_level1_fitbit_random_var.rds"))) %>%
  ungroup() %>%
  select(id,
         item,
         amp_median,
         mesor_median,
         sigma2_median,
         phi_median,
         phi_circ_median) %>%
  pivot_longer(amp_median:sigma2_median,
               names_to = "lin_variable",
               values_to = "lin_value") %>%
  pivot_longer(phi_median:phi_circ_median,
               names_to = "phi_version",
               values_to = "phi_value") %>%
  # mutate(lin_value = case_when(grepl("mesor", lin_variable) ~ lin_value /
  #                                10,
  #                              TRUE ~ lin_value)) %>%
  group_by(lin_variable, phi_version, item) %>%
  mutate(
    lin_value = lin_value/max(abs(lin_value)),
    x = lin_value * cos(phi_value),
    y = lin_value * sin(phi_value),
    r = median(lin_value),
    r_mean = mean(lin_value),
    x_mean = mean(x),
    y_mean = mean(y)
  )



cor_results <- ee %>%
  group_by(lin_variable, phi_version, item) %>%
  summarize(correlations = list(compute_correlations(phi_value, lin_value)),
            .groups = "drop") %>%
  unnest_wider(correlations) %>%
  mutate(
    cor_text = paste(
      "Pearson:",
      round(pearson_cor, 2),
      get_significance(pearson_pval),
      "\n",
      "Spearman:",
      round(spearman_cor, 2),
      get_significance(spearman_pval),
      "\n",
      "Kendall:",
      round(kendall_cor, 2),
      get_significance(kendall_pval),
      "\n",
      "JWM:",
      round(jwm_cor, 2),
      get_significance(jwm_pval),
      "\n",
      "Mardia:",
      round(mard_cor, 2),
      get_significance(mard_pval)
    )
  )

ee %>%
  left_join(cor_results,
            by = c("lin_variable",
                   "phi_version",
                   "item")) %>%
  ggplot() +
  aes(x = x, y = y) +
  geom_point(size = 0.5,
             alpha = 0.5) +
  ggthemes::theme_tufte() +
  theme(aspect.ratio = 1) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  ggforce::geom_circle(aes(x0 = 0, y0 = 0, r = r),
                       alpha = 0.4,
                       color = "red",
                       linewidth = 0.2) +
  ggforce::geom_circle(aes(x0 = 0, y0 = 0, r = r_mean),
                       alpha = 0.4,
                       color = "blue",
                       linewidth = 0.2) +
  facet_nested(lin_variable + phi_version ~ item,
               scales = "fixed") +
  geom_text(aes(x = 1, y = 0.8, label = cor_text),
            size = 2,
            hjust = 1) +
  xlim(-1.1, 1.1) +
  ylim(-1.1, 1.1)

# Save the above as 25cmX30cm (10inchX12inch) PDF


#######################

cor_text <- function(theta, x){
cor_values <- compute_correlations(theta, x)
# Format the text
paste(
  " Pearson:", round(cor_values$pearson_cor, 2), get_significance(cor_values$pearson_pval), "\n",
  "Spearman:", round(cor_values$spearman_cor, 2), get_significance(cor_values$spearman_pval), "\n",
  "Kendall:", round(cor_values$kendall_cor, 2), get_significance(cor_values$kendall_pval), "\n",
  "JWM:", round(cor_values$jwm_cor, 2), get_significance(cor_values$jwm_pval), "\n",
  "Mardia:", round(cor_values$mard_cor, 2), get_significance(cor_values$mard_pval)
)
# %>%
#   cat()
}



set.seed(2024-12-19)
n <- 1000
theta <- runif(n, 0, 2*pi)
co <- cos(theta)
er <- rnorm(n)

f_plotcor <- function(constant = 1,
                      coefficient = 1,
                      error_var = 0.3,
                      theta_min = 0,
                      theta_max = 2*pi,
                      meh = 1,
                      n = 1000,
                      return_df = FALSE){

  theta <- runif(n, theta_min, theta_max)
  co <- cos(theta)
  er <- rnorm(n)

df <- data.frame(theta = theta,
                 x = constant +
                   max(coefficient*cos(meh*theta), 0.3) +
                   sqrt(error_var)*er)
# %>%
#   filter(theta > theta_min,
#          theta < theta_max)
if(return_df == TRUE) return(df)

xmax <- max(abs(df$x))
p_polar <- df %>%
  ggplot() +
  aes(x = x * cos(theta), y = x * sin(theta)) +
  geom_point(size = 0.5,
             alpha = 0.5) +
  ggthemes::theme_tufte() +
  theme(aspect.ratio = 1) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  ggforce::geom_circle(
    aes(x0 = 0, y0 = 0, r = median(x)),
    alpha = 0.4,
    color = "red",
    linewidth = 0.5
  ) +
  ggforce::geom_circle(
    aes(x0 = 0, y0 = 0, r = mean(x)),
    alpha = 0.4,
    color = "blue",
    linewidth = 0.5
  ) +
  annotate(
    "text",
    x = -0.7 * xmax,
    y = 0.7 * xmax,
    label = cor_text(df$theta,
                     df$x)
  ) +
  annotate(
    "text",
    x = 0.7 * xmax,
    y = 0.7 * xmax,
    label = cor_text((df$theta + pi / 2) %% (2 * pi),
                     df$x),
    color = "orange"
  ) +
  xlim(-xmax, xmax) +
  ylim(-xmax, xmax) +
  ggtitle(paste0("x = ",
                 constant,
                 " + ",
                 coefficient,
                 " * cos(theta) + er, er ~ N(0,",
                 error_var,
                 ")"),
          paste0((theta_min/pi) %>% round(2),
                 " < theta <",
                 (theta_max/pi) %>% round(2),
                 "*pi"))

p_scatter1 <- df %>%
  ggplot(aes(x = theta, y = x)) +
  geom_point(size = 0.5, alpha = 0.5) +
  ggthemes::theme_tufte() +
  xlim(0, 2*pi)

p_scatter2 <- df %>%
  ggplot(aes(x = (theta + pi / 2) %% (2 * pi), y = x)) +
  geom_point(size = 0.5, alpha = 0.5, color = "orange") +
  ggthemes::theme_tufte() +
  xlim(0, 2*pi)

p_polar + (p_scatter1/p_scatter2)

}
