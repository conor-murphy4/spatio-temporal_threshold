library(tidyverse)

# Points on x and y axis at which to evaluate intensity
axis_ticks <- seq(from = -10, to = 10, length.out = 501)

# Construct covariate grid
a <- 2
b <- 4
covariate <- tibble(
  x = rep(axis_ticks, times = length(axis_ticks)),
  y = rep(axis_ticks, each  = length(axis_ticks)),
  r = sqrt(x^2 + y^2),
  theta = atan2(y, x), # angle in radians - not actually needed
  g = a * abs(sin(b * r))
)

covariate_for_conor <- covariate %>% select(x,y,g)

covariate_for_conor %>%
  ggplot(mapping = aes(x = x, y = y)) +
  geom_raster(mapping = aes(fill = g))


# TODO: Simulate point patterns from the true intensity lambda(r,theta) = beta_0 + beta_1 g, where beta_0 = 4 and beta_1 = 20

true_intensity <- tibble(covariate_for_conor, lambda = 4 + 20*g)
true_intensity %>%
  ggplot(mapping = aes(x = x, y = y)) +
  geom_raster(mapping = aes(fill = lambda))

max_intensity <- max(true_intensity$lambda)

num_points <- rpois(1, max_intensity*(20^2))
HPP_points <- tibble(
  x = runif(num_points, -10, 10),
  y = runif(num_points, -10, 10)
)

true_intensity %>%
  ggplot(mapping = aes(x = x, y = y)) +
  geom_raster(mapping = aes(fill = lambda))+
  geom_point(data=HPP_points)

lambda_fun <- function(x,y){
  r <- sqrt(x^2 + y^2)
  return(4+20*(2*abs(sin(4*r))))
}

HPP_points <- tibble(HPP_points, p_thin=lambda_fun(HPP_points$x, HPP_points$y)/max_intensity)
u <- runif(num_points)
keep <- HPP_points$p_thin > u

NHPP_points <- tibble(HPP_points[keep,])

true_intensity %>%
  ggplot(mapping = aes(x = x, y = y)) +
  geom_raster(mapping = aes(fill = lambda))+
  geom_point(data=HPP_points)+
  geom_point(data=NHPP_points, aes(colour="red"))



# TODO: Use one simulated point pattern to estimate parameters beta_0 and beta_1.

# TODO: Visualise your point estimate of the fitted intensity function

# TODO: Check that over many simulated point patterns your estimation method recovers the true parameters

# TODO: Visualise the standard error of the fitted intensity function



