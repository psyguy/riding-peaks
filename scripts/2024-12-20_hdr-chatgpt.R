# Load required packages
require(ks)
require(ggplot2)
require(metR)

compute_alpha_hdr <-
  function(d,
           x_0 = 0,
           y_0 = 0,
           n_grid = 200,
           plot = FALSE) {
    # 1. Kernel density estimation
    fhat <- kde(d, gridsize = c(n_grid, n_grid))

    xseq <- fhat$eval.points[[1]]
    yseq <- fhat$eval.points[[2]]
    zmat <- fhat$estimate

    dx <- diff(xseq)[1]
    dy <- diff(yseq)[1]
    cell_area <- dx * dy

    # 2. Compute cumulative probability distribution to determine HDR thresholds
    dens_vec <- as.vector(zmat)
    dens_sorted <- sort(dens_vec, decreasing = TRUE)
    cumulative_mass <- cumsum(dens_sorted * cell_area)
    cumulative_mass <-
      cumulative_mass / max(cumulative_mass)  # normalize

    # 3. Evaluate density at the given point (x_0, y_0)
    f_point <- predict(fhat, x = cbind(x_0, y_0))

    # Find where f_point fits in the density distribution
    # Identify the position of f_point in the sorted densities
    idx <- which(dens_sorted < f_point)[1]

    if (is.na(idx)) {
      # f_point is >= all densities, alpha_hdr ~ 1
      alpha_hdr <- 1
    } else {
      # f_point lies between dens_sorted[idx-1] and dens_sorted[idx]
      # if idx > 1, cumulative_mass[idx-1] gives us a close approximation
      if (idx > 1) {
        alpha_hdr <- cumulative_mass[idx - 1]
      } else {
        # idx = 1 means f_point < dens_sorted[1], so alpha is slightly less than cumulative_mass[1]
        alpha_hdr <- cumulative_mass[1]
      }
    }

    # If no plot is needed, just return alpha_hdr
    if (!plot) {
      return(round(alpha_hdr * 100, 1))
    }

    # If plot = TRUE, create the plot

    # 4. Prepare data for ggplot
    dens_df <- expand.grid(x = xseq, y = yseq)
    dens_df$z <- as.vector(zmat)

    # Desired HDR contour percentages and corresponding density thresholds
    cont_levels <- c(80, 85, 90, 95, 99)
    hdr_levels <- contourLevels(fhat, cont = cont_levels)

    # Create the ggplot
    p <- ggplot(dens_df, aes(x = x, y = y)) +
      geom_raster(aes(fill = z)) +
      scale_fill_viridis_c(option = "viridis",
                           direction = 1,
                           name = "Density") +
      coord_fixed() +
      labs(title = "HDR Contours with Heatmap",
           x = "X", y = "Y") +
      # Add HDR contour lines
      geom_contour(aes(z = z),
                   breaks = hdr_levels,
                   colour = "black",
                   size = 0.7) +
      # Point and its HDR label
      geom_point(aes(x = x_0, y = y_0), colour = "red", size = 3) +
      geom_text(
        aes(
          x = x_0,
          y = y_0,
          label = paste0(round(alpha_hdr * 100, 1), "%")
        ),
        colour = "red",
        hjust = -0.1,
        vjust = 0.5,
        size = 5
      ) +
      # Label the contour lines with their corresponding HDR percentages
      geom_text_contour(
        aes(z = z),
        breaks = hdr_levels,
        label_placer = label_placer_n(n = 1),
        size = 3,
        stroke = 0.2,
        label_format = function(breaks) {
          sapply(breaks, function(brk) {
            hdr_levels[which(hdr_levels == brk)]
          })
        }
      )

    print(p)

    # Return the alpha_hdr value
    return(round(alpha_hdr * 100, 1))
  }
