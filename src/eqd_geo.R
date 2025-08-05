# Load helper functions for working with the generalised Pareto distribution
source('src/helper_functions.R')

#-----------------------------------------------------------------------

#' Threshold selection method based on development of geophone network for extreme value modelling in the specific context of induced seismicity in Groningen in the Netherlands.
#'
#' 'eqd_geo_ics' selects a spatio-temporal threshold (based on the evolution of the geophone network in Groningen) above which the data can be most closely modelled 
#' by a covariate-dependent Generalised Pareto distribution.
#'
#' @author Conor Murphy
#'
#' @param data A numeric vector of earthquake magnitudes.
#' @param thresh A 2-d matrix of combinations of parameters (a,b) for proposed thresholds of the form a+bV where V denotes the distance_to_geo function for the given set of magnitudes.
#' @param distance_to_geo A numeric vector of the distances (to geophones) for the given set of earthquakes.
#' @param ics A numeric vector of the stress values for the given set of earthquakes. 
#' @param k  A positive integer denoting the number of bootstraps.
#' @param m A positive integer denoting the number of equally-spaced probabilities at which to evaluate quantiles.
#'
#' @returns A list containing the chosen threshold, the parameters of the fitted GPD, the number of observations above the chosen thresholds and the metric values 'd' corresponding to each proposed threshold.


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
  for (i in 1:nrow(thresh)) {
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


