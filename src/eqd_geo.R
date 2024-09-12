# Load helper functions for working with the generalised Pareto distribution
source('src/helper_functions.R')

#-----------------------------------------------------------------------

#' Threshold selection method based on development of geophone network for extreme value modelling in the specific context of induced seismicity in Groningen in the Netherlands.
#'
#' 'eqd_geo' selects a spatio-temporal threshold (based on the evolution of the geophone network in Groningen) above which the data can be most closely modelled by a non-stationary Generalised Pareto distribution.
#'
#' @author Conor Murphy
#'
#' @param data A numeric vector of seismic magnitudes.
#' @param thresh A 2-d matrix of combinations of parameters for a proposed thresholds of the form a+bV where V denotes the third_nearest_distance function for the given set of magnitudes.
#' @param third_nearest_distance A numeric vector of the third nearest distances (to geophones) for the given set of magnitudes.
#' @param k  A positive integer denoting the number of bootstraps.
#' @param m A positive integer denoting the number of equally-spaced probabilities at which to evaluate quantiles.
#'
#' @returns A list containing the chosen threshold, the parameters of the fitted GPD, the number of observations above the chosen thresholds and the metric values 'd' corresponding to each proposed threshold.


eqd_geo <- function(data, thresh, third_nearest_distance, k = 100, m = 500, min_dist=0, max_dist = 33.782, underlying_thresh = 0){

  # Check inputs are valid
  if (!is.numeric(data)) stop("Data must be a vector")
  if (!is.matrix(thresh)) stop("Thresholds must be a matrix")
  if (!is.numeric(third_nearest_distance)) stop("Third nearest distances must be a vector")
  if (length(third_nearest_distance) != length(data)) stop("Length of third nearest distances must match length of data")
  if (k <= 0 | k %% 1 != 0) stop("Number of bootstrapped samples must be a positive integer")
  if (m <= 0 | m %% 1 != 0) stop("Number of equally spaced probabilities must be a positive integer")

  meandistances <- xis <- sigmas <- num_excess <- numeric(nrow(thresh))
  for (i in 1:nrow(thresh)) {
    print(i)
    u <- thresh[i,1] + thresh[i,2] * third_nearest_distance
    if(any(u < underlying_thresh)) stop(paste("Candidate thresholds must be above the underlying threshold of ", underlying_thresh))
    excess <- data[data > u] - u[data > u]
    third_nearest_excess <- third_nearest_distance[data > u]
    num_excess[i] <- length(excess)
    if (num_excess[i] > 10) {
      mle0 <- mean(excess)
      init.fit <- optim(GPD_LL_given_third_nearest, excess = excess, par = c(mle0,0.1), control = list(fnscale = -1), thresh_par=thresh[i,], third_nearest = third_nearest_excess, min_dist = min_dist, max_dist = max_dist)
      xis[i] <- init.fit$par[[2]]
      sigmas[i] <- init.fit$par[[1]]
      distances <- numeric(k)
      for (j in 1:k) {
        #possible way to bootstrap
        sampled_indices <- sample(1:num_excess[i], num_excess[i], replace = TRUE)
        excess_boot <- excess[sampled_indices]
        third_nearest_excess_boot <- third_nearest_excess[sampled_indices]
        init_mle <- mean(excess_boot)
        ifelse(xis[i] < 0, pars_init <-  c(init_mle, 0.1) ,pars_init <- c(sigmas[i], xis[i]) )
        gpd.fit <- optim(GPD_LL_given_third_nearest, excess = excess_boot, par = pars_init, control = list(fnscale = -1), thresh_par=thresh[i,], third_nearest = third_nearest_excess_boot, min_dist = min_dist, max_dist = max_dist)
        sigma_given_third_nearest <- gpd.fit$par[[1]] + gpd.fit$par[[2]] * (thresh[i,1] + thresh[i,2]*third_nearest_excess_boot)
        transformed_excess_boot <- transform_to_exp(y = excess_boot, sig = sigma_given_third_nearest, xi = gpd.fit$par[[2]])
        distances[j] <- mean(abs(quantile(excess_boot, probs = (1:m) / (m+1)) - qexp((1:m) / (m+1), rate = 1)))
      }
      meandistances[i] <- mean(distances)
    }
    else{
      meandistances[i] <- NA
    }
  }
  chosen_index <- which.min(meandistances)
  chosen_threshold_par <- thresh[chosen_index,]
  xi <- xis[chosen_index]
  sigma <- sigmas[chosen_index]
  len <- num_excess[chosen_index]
  result <- list(thresh_par = chosen_threshold_par, par = c(sigma,xi), num_excess = len, dists = meandistances)
  return(result)
}

eqd_geo_unconstrained <- function(data, thresh, third_nearest_distance, k = 100, m = 500, underlying_thresh = 0){
  
  # Check inputs are valid
  if (!is.numeric(data)) stop("Data must be a vector")
  if (!is.matrix(thresh)) stop("Thresholds must be a matrix")
  if (!is.numeric(third_nearest_distance)) stop("Third nearest distances must be a vector")
  if (length(third_nearest_distance) != length(data)) stop("Length of third nearest distances must match length of data")
  if (k <= 0 | k %% 1 != 0) stop("Number of bootstrapped samples must be a positive integer")
  if (m <= 0 | m %% 1 != 0) stop("Number of equally spaced probabilities must be a positive integer")
  
  meandistances <- xis <- sigmas <- num_excess <- numeric(nrow(thresh))
  for (i in 1:nrow(thresh)) {
    print(i)
    u <- thresh[i,1] + thresh[i,2] * third_nearest_distance
    if(any(u < underlying_thresh)) stop(paste("Candidate thresholds must be above the underlying threshold of ", underlying_thresh))
    excess <- data[data > u] - u[data > u]
    third_nearest_excess <- third_nearest_distance[data > u]
    num_excess[i] <- length(excess)
    if (num_excess[i] > 10) {
      mle0 <- mean(excess)
      init.fit <- optim(GPD_LL_given_third_nearest_unconstrained, excess = excess, par = c(mle0,0.1), control = list(fnscale = -1), thresh_par=thresh[i,], third_nearest = third_nearest_excess)
      xis[i] <- init.fit$par[[2]]
      sigmas[i] <- init.fit$par[[1]]
      distances <- numeric(k)
      for (j in 1:k) {
        #possible way to bootstrap
        sampled_indices <- sample(1:num_excess[i], num_excess[i], replace = TRUE)
        excess_boot <- excess[sampled_indices]
        third_nearest_excess_boot <- third_nearest_excess[sampled_indices]
        init_mle <- mean(excess_boot)
        ifelse(xis[i] < 0, pars_init <-  c(init_mle, 0.1) ,pars_init <- c(sigmas[i], xis[i]) )
        gpd.fit <- optim(GPD_LL_given_third_nearest_unconstrained, excess = excess_boot, par = pars_init, control = list(fnscale = -1), thresh_par=thresh[i,], third_nearest = third_nearest_excess_boot)
        sigma_given_third_nearest <- gpd.fit$par[[1]] + gpd.fit$par[[2]] * (thresh[i,1] + thresh[i,2]*third_nearest_excess_boot)
        transformed_excess_boot <- transform_to_exp(y = excess_boot, sig = sigma_given_third_nearest, xi = gpd.fit$par[[2]])
        distances[j] <- mean(abs(quantile(excess_boot, probs = (1:m) / (m+1)) - qexp((1:m) / (m+1), rate = 1)))
      }
      meandistances[i] <- mean(distances)
    }
    else{
      meandistances[i] <- NA
    }
  }
  chosen_index <- which.min(meandistances)
  chosen_threshold_par <- thresh[chosen_index,]
  xi <- xis[chosen_index]
  sigma <- sigmas[chosen_index]
  len <- num_excess[chosen_index]
  result <- list(thresh_par = chosen_threshold_par, par = c(sigma,xi), num_excess = len, dists = meandistances)
  return(result)
}

