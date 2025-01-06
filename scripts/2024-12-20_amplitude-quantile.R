f_ampquant <- function(x, y, x_test = 0, y_test = 0, gridsize = 151){

  data <- cbind(x, y)
  # Perform kernel density estimation
  kde_result <- ks::kde(data, gridsize = gridsize)
  # Extract grid points and density matrix
  gx <- kde_result$eval.points[[1]]
  gy <- kde_result$eval.points[[2]]
  m <- kde_result$estimate
  # Normalise the density matrix to represent probabilities
  m_normalised <- m / sum(m)
  # Flatten the grid and density matrix into a data frame
  grid_df <- expand.grid(gx = gx, gy = gy) %>%
    mutate(density = as.vector(m_normalised))
  # Compute cumulative probabilities for contours
  grid_df <- grid_df %>%
    arrange(desc(density)) %>%
    mutate(cumulative_prob = cumsum(density) * 100) # Cumulative probabilities in percentages
  # Define percentile levels
  contour_levels <- c(50, 75, 80, 85, 90, 95, 99)
  # Find the quantile for the point (x, y)
  x_index <- which.min(abs(gx - x_test))
  y_index <- which.min(abs(gy - y_test))
  density_at_point <- m_normalised[y_index, x_index]
  sorted_densities <- sort(as.vector(m_normalised), decreasing = TRUE)
  cumulative_prob <- cumsum(sorted_densities) * 100
  density_rank <- which(sorted_densities == density_at_point)
  point_quantile <- cumulative_prob[density_rank]

  return(point_quantile %>% round(1))
}


print("pa")
(tt <- Sys.time())

d_draws_fitbit <- as_draws_df(readRDS("fits/brms_fitbit_random_var.rds")) %>%
  select(contains("r_id[") | contains("r_id__sigma[")) %>%
  mutate(iteration = 1:n()) %>%
  pivot_longer(
    cols = starts_with("r_id"),
    names_to = "id_x_variable",
    values_to = "value"
  ) %>%
  # Extract id (number between [ and ,)
  mutate(
    id_x_variable =
      str_replace(id_x_variable,
                  "r_id__sigma\\[(\\d+),Intercept\\]",
                  "r_id__sigma[\\1,sigmaIntercept]"),
    id = str_extract(id_x_variable,
                     "(?<=\\[)[0-9]+(?=,)") %>%
      as.numeric(),
    # Extract parameter (characters between , and ])
    parameter = str_extract(id_x_variable,
                            "(?<=,)[^\\]]+(?=\\])"),
    # Clean up parameter by removing whitespace
    parameter = str_trim(parameter)
  ) %>%
  # Remove original variable column and reorder
  select(iteration, id, parameter, value) %>%
  pivot_wider(names_from = "parameter",
              values_from = "value") %>%
  mutate(
    ## Not sure if sigmaIntercept is on log scale or not
    ## Commenting it out now
    log_sigma = sigmaIntercept,
    sigma2 = exp(sigmaIntercept)^2,
    mesor = Intercept,
    amp = sqrt(si^2 + co^2),
    phi = atan2(si, co) %% (2*pi),
    phi_atan = atan(si/co) %% (2*pi),
    psi = phi*24/(2*pi),
    psi_atan = phi_atan*24/(2*pi)
  ) %>%
  mutate(phi.begins6 = (phi - pi/2) %% (2*pi),
         phi.begins12 = (phi - pi) %% (2*pi)) %>%
  # filter(id < 101, id>70) %>%
  group_by(id)



# %>%
#   summarize(amp_hdr = compute_alpha_hdr(cbind(co,si))) %>%
#   cbind(item = "na", .)

Sys.time() - tt
beep(10)
