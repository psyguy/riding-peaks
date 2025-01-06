## In this I put codes for investigating correlations

library(dplyr)
library(metR) # Required for geom_text_contour
plot_density_with_contours <- function(data, x_test = 0, y_test = 0, gridsize = 151) {
  # Perform kernel density estimation
  kde_result <- ks::kde(data, gridsize = gridsize)

  # Extract grid points and density matrix
  gx <- kde_result$eval.points[[1]]
  gy <- kde_result$eval.points[[2]]
  m <- kde_result$estimate

  # Normalise the density matrix
  m_normalised <- m / sum(m)

  # Find the closest grid indices for (x_test, y_test)
  x_index <- which.min(abs(gx - x_test))
  y_index <- which.min(abs(gy - y_test))

  # Density at the specified point
  density_at_point <- m_normalised[y_index, x_index]

  # Sort all densities in descending order
  sorted_densities <- sort(as.vector(m_normalised), decreasing = TRUE)
  cumulative_prob <- cumsum(sorted_densities) * 100

  # Find the rank of the point's density approximately
  density_rank <- which.min(abs(sorted_densities - density_at_point))

  # Corresponding "quantile"
  point_quantile <- cumulative_prob[density_rank]
  # Plot the density and contours
  p <- ggplot(grid_df, aes(x = gx, y = gy)) +
    geom_raster(aes(fill = density)) +
    scale_fill_viridis_c(option = "viridis") +  # Use a colour palette for the density
    geom_contour(aes(z = cumulative_prob), breaks = contour_levels, colour = "white") +
    geom_text_contour(aes(z = cumulative_prob), breaks = contour_levels, colour = "black", size = 5) + # Larger text size
    labs(
      title = "Kernel Density Estimation with Percentile Contours",
      x = "X-axis",
      y = "Y-axis",
      fill = "Density"
    ) +
    theme_minimal()
  # Add coordinate system with native axes crossing at (0, 0)
  p <- p +
    geom_hline(yintercept = 0, colour = "black", size = 0.8) +
    geom_vline(xintercept = 0, colour = "black", size = 0.8) +
    coord_fixed(ratio = 1, expand = FALSE, xlim = range(gx), ylim = range(gy)) +
    theme(
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      panel.grid = element_blank()
    )
  # Add the point (x, y) with quantile label
  p <- p +
    annotate("point", x = x_test, y = y_test, color = "red", size = 3) +
    annotate("text", x = x_test, y = y_test,
             label = paste0("Quantile: ", round(point_quantile, 1), "%"),
             hjust = -0.2, colour = "red", size = 5)
  # Return the plot
  return(p)
}
# Example usage
set.seed(123)
data <- matrix(rnorm(200), ncol = 2)
plot <- plot_density_with_contours(data)
print(plot)
plot + theme_classic()
plot + theme_bw()
