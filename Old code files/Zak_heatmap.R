# Plotting matrices as raster images.

## EXAMPLE 1: using the outer function

# make the length of these different so we can see if it is the wrong way
x_values <- seq(1,10,length.out = 10)
y_values <- seq(1,5,length.out = 5)

z <- outer(x_values, y_values, "+")

fields::image.plot(x_values, y_values, z, xlab = "Easting (m)", ylab = "Northing (m)")
# Note: Each row is one x value - use transpose if you store in "sensible" way.

## EXAMPLE 2: Distance from (\pi, \pi)

distance_to_pi <- function(xy){
  sqrt(sum((xy - c(pi, pi))^2))
}

xy_pairs <- expand.grid(x_val = x_values, y_val =  y_values)

z_values <- apply(xy_pairs, MARGIN = 1, FUN = distance_to_pi)

z_mat <- matrix(
  data = z_values,
  nrow = length(x_values),
  ncol = length(y_values),
  byrow = FALSE)

z_mat

fields::image.plot(x_values, y_values, z_mat)

## EXAMPLE 3: Distance from (\pi, \pi) using `outer`

d_to_pi <- function(x, y){
  sqrt((x-pi)^2 + (y- pi)^2)
}

z_mat_3 <- outer(x_values, y_values, d_to_pi)
fields::image.plot(x_values, y_values, z_mat_3)
