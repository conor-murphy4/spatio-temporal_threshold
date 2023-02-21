
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

