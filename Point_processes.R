
#-----------------------Poisson process simulations with covariate-dependent intensity ---------------------

#--------------------------------------1-D-------------------------------------------

#time-frame
t_min=pi
t_max=5*pi


#Covariate depending on time 
a <- 50
b <- 40

z_fun <- function(t){
  cov <- a + b*sin(t)
  return(cov)
}

#Intensity function
Beta <- c(0.5)
lambda_z <- function(t){
  beta_0 <- Beta[1]
  z <- z_fun(t)
  intensity <- beta_0*z
  return(intensity)  
}

#Finding max of lambda
neg_lambda <- function(t){
  return(-lambda_z(t))
}
opt_max <- optimize(neg_lambda, c(t_min, t_max))
lambda_max <- -opt_max$objective

#Simulate homo PP
num_points <- rpois(1, lambda_max*(t_max-t_min))
HPP_points <- runif(num_points, t_min, t_max)

thin <- function(t){
  return(lambda_z(t)/lambda_max)
}

#thinning probs
p_thin <- thin(HPP_points)
u <- runif(num_points)

#thinning
keep <- p_thin > u

NHPP_points <- HPP_points[keep]
n_NHPP <- length(NHPP_points)

#Plotting intensity with simulated points
t <- seq(t_min, t_max, length.out=1000)
plot(t, lambda_z(t), ylim = c(0, lambda_max))

points(NHPP_points, rep(0.5, n_NHPP), col="red", pch=19)


#-------------Simulating covariate on grid----------

#time-frame
t_min=3
t_max=16


#Covariate depending on time 
a <- 50
b <- 50

z_t <- function(t){
  z <- a + b*sin(t)
  return(z)
}
#Grid of time-points where covariate was recorded
t_grid <- seq(t_min, t_max, by=1)


#Intensity function
beta_0 <- 0.5
lambda_z <- function(t){
  intensity <- beta_0*z_t(t)
  return(intensity)  
}

#Finding max of lambda over grid
z_grid <- z_t(t_grid)
max_rec_lambda <- max(lambda_z(t_grid))

#Simulating HPP
num_points <- rpois(1, max_rec_lambda*(t_max-t_min))
HPP_points <- runif(num_points, t_min, t_max)
plot(t_grid, z_t(t_grid),col="red", pch=19 ,xlim=c(t_min,t_max),ylim=c(0, max(z_t(t_grid))))
points(t_grid, lambda_z(t_grid), pch=19)
interp_z <- approx(t_grid, z_grid, xout = HPP_points)
points(interp_z$x, interp_z$y, col="blue", pch=19, cex=0.2 )
interp_lambda <- approx(t_grid, lambda_z(t_grid), xout=HPP_points)
points(interp_lambda$x, interp_lambda$y, col="green", cex=0.2)

#Approximating lambda function
lambda_approx <- approxfun(t_grid, lambda_z(t_grid))
curve(lambda_approx, add=TRUE)

#Thinning
thin <- function(lambda){
  return(lambda/max_rec_lambda)
}

#Thinning
p_thin <- thin(interp_lambda$y)
u <- runif(num_points)

keep <- p_thin > u 

NHPP_points <- HPP_points[keep]
n_NHPP <- length(NHPP_points)

points(NHPP_points, rep(0, n_NHPP), col="grey", cex=2)


#------------------------ Estimation -------------------------------------------

Lambda_z <- function(lambda_fun){
  Lambda <- integrate(lambda_fun, t_min, t_max)$value
  return(Lambda)
}

NHPP_neg_LL <- function(par, points, z_t){
  beta_0 <- par[1]
  Lambda <- Lambda_z(function(t){beta_0*z_t(t)})
  term_1 <- log(length(points)) + Lambda
  term_2 <- -sum(log(beta_0*z_t(points)))
  return(term_1 + term_2)
}

(opt <- optimize(NHPP_neg_LL, interval = c(0,5), points=NHPP_points, z_t=z_t))


#------------------------ 2-D covarite dependent point pattern-------------------------------------

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


#--------------- Intensities from Zak's thesis ---------------------------------

#Simulated covariate

# Points on x and y axis at which to evaluate intensity
lower_limit<- -10
upper_limit<- 10
axis_ticks <- seq(from = lower_limit, to = upper_limit, length.out = 501)
area <- (u-l)^2
# Construct covariate grid
a <- 10
b <- 4
covariate <- tibble(
  x = rep(axis_ticks, times = length(axis_ticks)),
  y = rep(axis_ticks, each  = length(axis_ticks)),
  r = sqrt(x^2 + y^2),
  theta = atan2(y, x), # angle in radians - not actually needed
  g = a * abs(sin(b * r))
)

covariate_sim<- covariate %>% select(x,y,g)

covariate_sim %>%
  ggplot(mapping = aes(x = x, y = y)) +
  geom_raster(mapping = aes(fill = g))
#---------------------Model S1-----------------------

true_intensity_S1 <- tibble(covariate_sim, lambda = 0.33*g)
true_intensity_S1 %>%
  ggplot(mapping = aes(x = x, y = y)) +
  geom_raster(mapping = aes(fill = lambda))

max_intensity <- max(true_intensity_S1$lambda)

num_points <- rpois(1, max_intensity*area)
HPP_points <- tibble(
  x = runif(num_points, lower_limit, upper_limit),
  y = runif(num_points, lower_limit, upper_limit)
)

true_intensity_S1 %>%
  ggplot(mapping = aes(x = x, y = y)) +
  geom_raster(mapping = aes(fill = lambda))+
  geom_point(data=HPP_points)

lambda_fun <- function(x,y){
  r <- sqrt(x^2 + y^2)
  g <- a*abs(sin(b*r))
  return(0.33*g)
}

HPP_points <- tibble(HPP_points, p_thin=lambda_fun(HPP_points$x, HPP_points$y)/max_intensity)
U <- runif(num_points)
keep <- HPP_points$p_thin > U

NHPP_points <- tibble(HPP_points[keep,])

true_intensity_S1 %>%
  ggplot(mapping = aes(x = x, y = y)) +
  geom_raster(mapping = aes(fill = lambda))+
  geom_point(data=HPP_points)+
  geom_point(data=NHPP_points, aes(colour="red"))

#-------------------- Model S2----------------------

covariate_sim <- covariate_sim %>% mutate(g_shift = 0.3+g) #Shifting as dividing by small values was giving an infinite rate
covariate_sim <- covariate_sim %>% mutate(rate=(g_shift-lag(g_shift))/lag(g_shift))

true_intensity_S2 <- tibble(covariate_sim, lambda = (9.6e-2)*abs(rate) )
true_intensity_S2 <- true_intensity_S2[!is.na(true_intensity_S2$rate),]
true_intensity_S2 %>%
  ggplot(mapping = aes(x = x, y = y)) +
  geom_raster(mapping = aes(fill = lambda))

max_intensity_S2 <- max(true_intensity_S2$lambda)


num_points <- rpois(1, max_intensity_S2*area)
HPP_points <- tibble(
  x = runif(num_points, lower_limit, upper_limit),
  y = runif(num_points, lower_limit, upper_limit)
)

true_intensity_S2 %>%
  ggplot(mapping = aes(x = x, y = y)) +
  geom_raster(mapping = aes(fill = lambda))+
  geom_point(data=HPP_points)

#Need to think of different way to do this!
lambda_fun <- function(x,y){
  r <- sqrt(x^2 + y^2)
  g <- a*abs(sin(b*r))
  g_shift <- 0.3 + g
  rate <- (g_shift-lag(g_shift))/lag(g_shift)
  return((9.6e-2)*abs(rate))
}

HPP_points <- tibble(HPP_points, p_thin=lambda_fun(HPP_points$x, HPP_points$y)/max_intensity_S2)
U <- runif(num_points)
keep <- HPP_points$p_thin > U

NHPP_points <- tibble(HPP_points[keep,])

true_intensity_S2 %>%
  ggplot(mapping = aes(x = x, y = y)) +
  geom_raster(mapping = aes(fill = lambda))+
  geom_point(data=HPP_points)+
  geom_point(data=NHPP_points, aes(colour="red"))


#----------------------Model S3---------------------
shift <- 0.02
covariate_sim <- covariate_sim %>% mutate(g_shift = shift+g) #Shifting as dividing by small values was giving an infinite rate
covariate_sim <- covariate_sim %>% mutate(rate=(g_shift-lag(g_shift))/lag(g_shift))

true_intensity_S3 <- tibble(covariate_sim, lambda = (7e-4)*abs(rate)*(1 + 21.8*g) )
true_intensity_S3 <- true_intensity_S3[!is.na(true_intensity_S3$rate),]
true_intensity_S3 %>%
  ggplot(mapping = aes(x = x, y = y)) +
  geom_raster(mapping = aes(fill = lambda))

max_intensity_S3 <- max(true_intensity_S3$lambda)


num_points <- rpois(1, max_intensity_S3*area)
HPP_points <- tibble(
  x = runif(num_points, lower_limit, upper_limit),
  y = runif(num_points, lower_limit, upper_limit)
)

true_intensity_S3 %>%
  ggplot(mapping = aes(x = x, y = y)) +
  geom_raster(mapping = aes(fill = lambda))+
  geom_point(data=HPP_points)

