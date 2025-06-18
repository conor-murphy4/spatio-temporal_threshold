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
#' @param thresh A 2-d matrix of combinations of parameters for a proposed thresholds of the form a+bV where V denotes the distance_to_geo function for the given set of magnitudes.
#' @param distance_to_geo A numeric vector of the third nearest distances (to geophones) for the given set of magnitudes.
#' @param k  A positive integer denoting the number of bootstraps.
#' @param m A positive integer denoting the number of equally-spaced probabilities at which to evaluate quantiles.
#'
#' @returns A list containing the chosen threshold, the parameters of the fitted GPD, the number of observations above the chosen thresholds and the metric values 'd' corresponding to each proposed threshold.


eqd_geo <- function(data, thresh, distance_to_geo, k = 100, m = 500, min_dist=0, max_dist = 33.782, underlying_thresh = 0){

  # Check inputs are valid
  if (!is.numeric(data)) stop("Data must be a vector")
  if (!is.matrix(thresh)) stop("Thresholds must be a matrix")
  if (!is.numeric(distance_to_geo)) stop("Third nearest distances must be a vector")
  if (length(distance_to_geo) != length(data)) stop("Length of distances_to_geo must match length of data")
  if (k <= 0 | k %% 1 != 0) stop("Number of bootstrapped samples must be a positive integer")
  if (m <= 0 | m %% 1 != 0) stop("Number of equally spaced probabilities must be a positive integer")

  meandistances <- xis <- sigmas <- num_excess <- numeric(nrow(thresh))
  for (i in 1:nrow(thresh)) {
    print(i)
    u <- thresh[i,1] + thresh[i,2] * distance_to_geo
    if(any(u < underlying_thresh)) stop(paste("Candidate thresholds must be above the underlying threshold of ", underlying_thresh))
    excess <- data[data > u] - u[data > u]
    V_excess <- distance_to_geo[data > u]
    num_excess[i] <- length(excess)
    if (num_excess[i] > 20) {
      mle0 <- mean(excess)
      init.fit <- optim(GPD_LL_given_V, excess = excess, par = c(mle0,0.1), control = list(fnscale = -1), thresh_par=thresh[i,], V = V_excess, min_dist = min_dist, max_dist = max_dist)
      xis[i] <- init.fit$par[[2]]
      sigmas[i] <- init.fit$par[[1]]
      distances <- numeric(k)
      j <- 1
      while(j  <= k) {
        #possible way to bootstrap
        sampled_indices <- sample(1:num_excess[i], num_excess[i], replace = TRUE)
        excess_boot <- excess[sampled_indices]
        V_excess_boot <- V_excess[sampled_indices]
        init_mle <- mean(excess_boot)
        ifelse(xis[i] < 0, pars_init <-  c(init_mle, 0.1) ,pars_init <- c(sigmas[i], xis[i]) )
        gpd.fit <- optim(GPD_LL_given_V, excess = excess_boot, par = pars_init, control = list(fnscale = -1), thresh_par=thresh[i,], V = V_excess_boot, min_dist = min_dist, max_dist = max_dist)
        sigma_given_V <- gpd.fit$par[[1]] + gpd.fit$par[[2]] * (thresh[i,1] + thresh[i,2]*V_excess_boot)
        if(any((1 + gpd.fit$par[[2]] * excess_boot/sigma_given_V) < 1e-323)) next
        else{
          transformed_excess_boot <- transform_to_exp(y = excess_boot, sig = sigma_given_V, xi = gpd.fit$par[[2]])
          sample_quants <- quantile(transformed_excess_boot, probs = (1:m) / (m+1))
          model_quants <- qexp((1:m) / (m+1), rate = 1)
          distances[j] <- mean(abs(sample_quants - model_quants))
          j <- j + 1
        }
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


eqd_geo_unconstrained <- function(data, thresh, distance_to_geo, k = 100, m = 500, underlying_thresh = 0){
  
  # Check inputs are valid
  if (!is.numeric(data)) stop("Data must be a vector")
  if (!is.matrix(thresh)) stop("Thresholds must be a matrix")
  if (!is.numeric(distance_to_geo)) stop("Third nearest distances must be a vector")
  if (length(distance_to_geo) != length(data)) stop("Length of third nearest distances must match length of data")
  if (k <= 0 | k %% 1 != 0) stop("Number of bootstrapped samples must be a positive integer")
  if (m <= 0 | m %% 1 != 0) stop("Number of equally spaced probabilities must be a positive integer")
  
  meandistances <- xis <- sigmas <- num_excess <- numeric(nrow(thresh))
  for (i in 1:nrow(thresh)) {
    print(i)
    u <- thresh[i,1] + thresh[i,2] * distance_to_geo
    if(any(u < underlying_thresh)) stop(paste("Candidate thresholds must be above the underlying threshold of ", underlying_thresh))
    excess <- data[data > u] - u[data > u]
    V_excess <- distance_to_geo[data > u]
    num_excess[i] <- length(excess)
    if (num_excess[i] > 20) {
      mle0 <- mean(excess)
      init.fit <- optim(GPD_LL_given_V_unconstrained, excess = excess, par = c(mle0,0.1), control = list(fnscale = -1), thresh_par=thresh[i,], V = V_excess)
      xis[i] <- init.fit$par[[2]]
      sigmas[i] <- init.fit$par[[1]]
      distances <- numeric(k)
      j <- 1
      while(j  <= k) {
        #possible way to bootstrap
        sampled_indices <- sample(1:num_excess[i], num_excess[i], replace = TRUE)
        excess_boot <- excess[sampled_indices]
        V_excess_boot <- V_excess[sampled_indices]
        init_mle <- mean(excess_boot)
        ifelse(xis[i] < 0, pars_init <-  c(init_mle, 0.1) ,pars_init <- c(sigmas[i], xis[i]) )
        gpd.fit <- optim(GPD_LL_given_V_unconstrained, excess = excess_boot, par = pars_init, control = list(fnscale = -1), thresh_par=thresh[i,], V = V_excess_boot)
        sigma_given_V <- gpd.fit$par[[1]] + gpd.fit$par[[2]] * (thresh[i,1] + thresh[i,2]*V_excess_boot)
        if(any((1 + gpd.fit$par[[2]] * excess_boot/sigma_given_V) < 1e-323)) next
        else{
          transformed_excess_boot <- transform_to_exp(y = excess_boot, sig = sigma_given_V, xi = gpd.fit$par[[2]])
          sample_quants <- quantile(transformed_excess_boot, probs = (1:m) / (m+1))
          model_quants <- qexp((1:m) / (m+1), rate = 1)
          distances[j] <- mean(abs(sample_quants - model_quants))
          j <- j + 1
        }
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

eqd_stepped <- function(data, thresh, change_index, k = 100, m = 500, underlying_thresh = 0){
  
  # Check inputs are valid
  if (!is.numeric(data)) stop("Data must be a vector")
  if (!is.matrix(thresh)) stop("Thresholds must be a matrix")
  if (k <= 0 | k %% 1 != 0) stop("Number of bootstrapped samples must be a positive integer")
  if (m <= 0 | m %% 1 != 0) stop("Number of equally spaced probabilities must be a positive integer")
  
  meandistances <- xis <- sigmas <- num_excess <- numeric(nrow(thresh))
  for (i in 1:nrow(thresh)) {
    print(i)
    u <- c(rep(thresh[i,1], change_index), rep(thresh[i,2], length(data) - change_index))
    if(any(u < underlying_thresh)) stop(paste("Candidate thresholds must be above the underlying threshold of ", underlying_thresh))
    excess <- data[data > u] - u[data > u]
    num_excess[i] <- length(excess)
    if (num_excess[i] > 10) {
      mle0 <- mean(excess)
      init.fit <- optim(GPD_LL_step, excess = excess, par = c(mle0,0.1), control = list(fnscale = -1), thresh=u)
      xis[i] <- init.fit$par[[2]]
      sigmas[i] <- init.fit$par[[1]]
      distances <- numeric(k)
      j <- 1
      while(j  <= k) {
        #possible way to bootstrap
        sampled_indices <- sample(1:num_excess[i], num_excess[i], replace = TRUE)
        excess_boot <- excess[sampled_indices]
        thresh_boot <- u[sampled_indices]
        init_mle <- mean(excess_boot)
        ifelse(xis[i] < 0, pars_init <-  c(init_mle, 0.1) ,pars_init <- c(sigmas[i], xis[i]) )
        gpd.fit <- optim(GPD_LL_step, excess = excess_boot, par = pars_init, control = list(fnscale = -1), thresh=thresh_boot)
        sigma_tilde <- gpd.fit$par[[1]] + gpd.fit$par[[2]] * (thresh_boot)
        if(any((1 + gpd.fit$par[[2]] * excess_boot/sigma_given_V) < 1e-323)) next
        else{
          transformed_excess_boot <- transform_to_exp(y = excess_boot, sig = sigma_tilde, xi = gpd.fit$par[[2]])
          sample_quants <- quantile(transformed_excess_boot, probs = (1:m) / (m+1))
          model_quants <- qexp((1:m) / (m+1), rate = 1)
          distances[j] <- mean(abs(sample_quants - model_quants))
          j <- j + 1
        }
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


eqd_geo_ics_max_dist <- function(data, thresh, distance_to_geo, ics, k = 100, m = 500, max_dist = 33.782, min_ics = 0, underlying_thresh = 0){
  # TODO Update code to allow for differing underlying threshold
  # TODO Assess chosen values of max_dist and min_ics
  
  # Check inputs are valid
  if (!is.numeric(data)) stop("Data must be a vector")
  if (!is.matrix(thresh)) stop("Thresholds must be a matrix")
  if (!is.numeric(distance_to_geo)) stop("distances_to_geo must be a vector")
  if (length(distance_to_geo) != length(data)) stop("Length of distances_to_geo must match length of data")
  if (!is.numeric(ics)) stop("ICS must be a vector")
  if (length(ics) != length(data)) stop("Length of ICS must match length of data")
  if (k <= 0 | k %% 1 != 0) stop("Number of bootstrapped samples must be a positive integer")
  if (m <= 0 | m %% 1 != 0) stop("Number of equally spaced probabilities must be a positive integer")
  
  meandistances <- xis <- num_excess <- numeric(nrow(thresh))
  sigma_par <- matrix(NA, nrow = nrow(thresh), ncol = 2)
  for (i in 1:nrow(thresh)) {
    print(i)
    u <- thresh[i,1] + thresh[i,2] * distance_to_geo
    if(any(u < underlying_thresh)) stop(paste("Candidate thresholds must be above the underlying threshold of ", underlying_thresh))
    excess <- data[data > u] - u[data > u]
    V_excess <- distance_to_geo[data > u]
    ICS_excess <- ics[data > u]
    num_excess[i] <- length(excess)
    if (num_excess[i] > 20) {
      mle0 <- mean(excess)
      init.fit <- optim(GPD_LL_given_V_ICS, excess = excess, par = c(mle0, 0, 0.1), control = list(fnscale = -1), thresh_par=thresh[i,], V = V_excess, ics = ICS_excess, max_dist = max_dist, min_ics = min_ics)
      xis[i] <- init.fit$par[3]
      sigma_par[i,] <- init.fit$par[1:2]
      distances <- numeric(k)
      j <- 1
      while(j  <= k) {
        #possible way to bootstrap
        sampled_indices <- sample(1:num_excess[i], num_excess[i], replace = TRUE)
        excess_boot <- excess[sampled_indices]
        V_excess_boot <- V_excess[sampled_indices]
        ICS_excess_boot <- ICS_excess[sampled_indices]
        init_mle <- mean(excess_boot)
        ifelse(xis[i] < 0, pars_init <-  c(init_mle, 0, 0.1) ,pars_init <- c(sigma_par[i,], xis[i]) )
        gpd.fit <- optim(GPD_LL_given_V_ICS, excess = excess_boot, par = pars_init, control = list(fnscale = -1), thresh_par=thresh[i,], V = V_excess_boot, 
                         ics = ICS_excess_boot, max_dist = max_dist, min_ics = min_ics)
        sigma_given_V <- gpd.fit$par[[1]] + gpd.fit$par[[2]]*ICS_excess_boot + gpd.fit$par[[3]] * (thresh[i,1] + thresh[i,2]*V_excess_boot)
        if(any((1 + gpd.fit$par[[3]] * excess_boot/sigma_given_V) < 1e-323)) next
        else{
          transformed_excess_boot <- transform_to_exp(y = excess_boot, sig = sigma_given_V, xi = gpd.fit$par[[3]])
          sample_quants <- quantile(transformed_excess_boot, probs = (1:m) / (m+1))
          model_quants <- qexp((1:m) / (m+1), rate = 1)
          distances[j] <- mean(abs(sample_quants - model_quants))
          j <- j + 1
        }
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
  sigmas <- sigma_par[chosen_index,]
  len <- num_excess[chosen_index]
  result <- list(thresh_par = chosen_threshold_par, par = c(sigmas,xi), num_excess = len, dists = meandistances)
  return(result)
}


eqd_geo_ics <- function(data, thresh, distance_to_geo, ics, k = 100, m = 500, underlying_thresh = 0){
  
  # Check inputs are valid
  if (!is.numeric(data)) stop("Data must be a vector")
  if (!is.matrix(thresh)) stop("Thresholds must be a matrix")
  if (!is.numeric(distance_to_geo)) stop("Third nearest distances must be a vector")
  if (length(distance_to_geo) != length(data)) stop("Length of third nearest distances must match length of data")
  if (k <= 0 | k %% 1 != 0) stop("Number of bootstrapped samples must be a positive integer")
  if (m <= 0 | m %% 1 != 0) stop("Number of equally spaced probabilities must be a positive integer")
  
  meandistances <- xis <- num_excess <- numeric(nrow(thresh))
  sigma_par <- matrix(NA, nrow = nrow(thresh), ncol = 2)
  ########################## extra code for Inf investigation ##################
  # inf_list_mean <- list()
  # inf_list_boot <- list()
  ########################## extra code for Inf investigation ##################
  for (i in 1:nrow(thresh)) {
    #print(i)
    u <- thresh[i,1] + thresh[i,2] * distance_to_geo
    if(any(u < underlying_thresh)) stop(paste("Candidate thresholds must be above the underlying threshold of ", underlying_thresh))
    excess <- data[data > u] - u[data > u]
    V_excess <- distance_to_geo[data > u]
    ICS_excess <- ics[data > u]
    num_excess[i] <- length(excess)
    if (num_excess[i] > 20) {
      mle0 <- mean(excess)
      init.fit <- optim(GPD_LL_given_V_ICS, excess = excess, par = c(mle0, 0, 0.1), control = list(fnscale = -1), 
                        thresh_par=thresh[i,], V = V_excess, ics=ICS_excess)
      xis[i] <- init.fit$par[3]
      sigma_par[i,] <- init.fit$par[1:2]
      distances <- numeric(k)
      j <- 1
      while(j  <= k) {
        #possible way to bootstrap
        sampled_indices <- sample(1:num_excess[i], num_excess[i], replace = TRUE)
        excess_boot <- excess[sampled_indices]
        V_excess_boot <- V_excess[sampled_indices]
        ICS_excess_boot <- ICS_excess[sampled_indices]
        init_mle <- mean(excess_boot)
        ifelse(xis[i] < 0, pars_init <-  c(init_mle, 0, 0.1) ,pars_init <- c(sigma_par[i,], xis[i]) )
        gpd.fit <- optim(GPD_LL_given_V_ICS, excess = excess_boot, par = pars_init, control = list(fnscale = -1), 
                         thresh_par=thresh[i,], V = V_excess_boot, ics = ICS_excess_boot)
        sigma_given_V <- gpd.fit$par[1] + gpd.fit$par[2]*ICS_excess_boot + gpd.fit$par[3] * (thresh[i,1] + thresh[i,2]*V_excess_boot)
        if(any((1 + gpd.fit$par[3] * excess_boot/sigma_given_V) < 1e-323)) next
        else{
          transformed_excess_boot <- transform_to_exp(y = excess_boot, sig = sigma_given_V, xi = gpd.fit$par[3])
          sample_quants <- quantile(transformed_excess_boot, probs = (1:m) / (m+1))
          model_quants <- qexp((1:m) / (m+1), rate = 1)
          distances[j] <- mean(abs(sample_quants - model_quants))
          ########################## extra code for Inf investigation ##########
          # if(distances[j] == Inf){
          #   inf_list_current <- list(i = i, j = j, excess_boot = excess_boot, transformed_excess = transformed_excess_boot, est = gpd.fit$par, sigma_var = sigma_given_V)
          #   inf_list_boot <- c(inf_list_boot, list(inf_list_current))
          # }
          ########################## extra code for Inf investigation ##########
          j <- j + 1
        }
      }
      meandistances[i] <- mean(distances)
      ########################## extra code for Inf investigation ##############
      # if(meandistances[i] == Inf){
      #   inf_list_mean <- c(inf_list_mean, list(list(i = i, excess = excess, par_ests = init.fit$par, dists = distances)))
      # }
      ########################## extra code for Inf investigation ##############
    }
    else{
      meandistances[i] <- NA
    }
  }
  chosen_index <- which.min(meandistances)
  chosen_threshold_par <- thresh[chosen_index,]
  xi <- xis[chosen_index]
  sigmas <- sigma_par[chosen_index,]
  len <- num_excess[chosen_index]
  result <- list(thresh_par = chosen_threshold_par, par = c(sigmas,xi), num_excess = len, dists = meandistances)
  ########################## extra code for Inf investigation ##################
  # result <- list(thresh_par = chosen_threshold_par, par = c(sigmas,xi), num_excess = len, dists = meandistances, inf_list_mean = inf_list_mean, inf_list_boot = inf_list_boot)
  ########################## extra code for Inf investigation ##################
  return(result)
}


