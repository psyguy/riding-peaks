

# Generate the cylinder coordinates
cylinder_radius <- max(df$x)  # Adjust the radius as needed
# Create x, y, and z coordinates
x <- df$x * cos(df$theta)
y <- df$x * sin(df$theta)
z <- df$x
# Create the plot in 3D
open3d()
par3d(scale = c(1, 1, 1),
      aspect = c(1, 1, 1))  # Maintain aspect ratio across all axes
plot3d(x, y, z, type = "p", col = "blue", size = 3,
       xlab = "X", ylab = "Y", zlab = "Height")
# Add the cylinder mesh for better visualisation
theta <- seq(0, 2 * pi, length.out = 100)
z_cylinder <- seq(min(z), max(z), length.out = 100)
# Create a grid of points for the cylinder surface
grid_theta <- matrix(rep(theta, each = length(z_cylinder)), nrow = length(z_cylinder))
grid_z <- matrix(rep(z_cylinder, times = length(theta)), nrow = length(z_cylinder))
# Create the mesh for the cylinder
mesh_x <- cylinder_radius * cos(grid_theta)
mesh_y <- cylinder_radius * sin(grid_theta)
mesh_z <- grid_z
# Add the cylinder surface to the plot
surface3d(mesh_x, mesh_y, mesh_z, alpha = 0.5, color = "lightblue")
# Draw the filled cylinder using quads3d
quads3d(
  c(mesh_x[-1, -1], mesh_x[-nrow(mesh_x), -1], mesh_x[-nrow(mesh_x), -ncol(mesh_x)], mesh_x[-1, -ncol(mesh_x)]),
  c(mesh_y[-1, -1], mesh_y[-nrow(mesh_y), -1], mesh_y[-nrow(mesh_y), -ncol(mesh_y)], mesh_y[-1, -ncol(mesh_y)]),
  c(mesh_z[-1, -1], mesh_z[-nrow(mesh_z), -1], mesh_z[-nrow(mesh_z), -ncol(mesh_z)], mesh_z[-1, -ncol(mesh_z)]),
  col = "lightblue", alpha = 0.5
)
# Add the original points on top of the filled cylinder
points3d(x, y, z, col = "blue", size = 3)
# Ensure the viewport aspect ratio is consistent
decorate3d(aspect = c(1, 1, 1))
