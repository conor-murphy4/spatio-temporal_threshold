
# Functions for estimation of intensity of earthquakes above a threshold

#' Helper function to compute grid box length in one direction
#' 
#' @param direction Numeric vector with coordinates in one direction (Easting or Northing)
#' 
#' @return Numeric value with grid box length in that direction
grid_box <- function(direction){
  (max(unique(direction)) - min(unique(direction))) / length(unique(direction))
}

#' Function for estimation of intensity above 0 given geophone based threshold
#' 
#' @param par Numeric vector with parameter values (gamma_0, gamma_1)
#' @param data Data frame with earthquake catalog including observed locations/magnitudes/covariates
#' @param covariates Data frame with covariates (ICS, threshold, V1-V4, etc.) across whole space and time
#' @param threshold_obs Numeric vector with threshold values for observed exceedances
#' @param covariates_threshold Numeric vector with threshold values at covariates grid points
#' @param thresh_fit List with threshold selection object GPD and threshold parameters
#' 
#' @return Numeric value of Poisson process log-likelihood
Poisson_process_LL_icsmax <- function(par, data, covariates, threshold_obs, covariates_threshold, thresh_fit){

  # Extract parameters
  gamma_0 <- par[1]
  gamma_1 <- par[2]
  
  exceedances <- data[data$Magnitude > threshold_obs, ]
  threshold_ex <- threshold_obs[data$Magnitude > threshold_obs]
  
  # Removing locations where dsmaxdt = 0
  exceedances_ds <- exceedances[exceedances$dsmaxdt > 0, ]
  threshold_ex_ds <- threshold_ex[exceedances$dsmaxdt > 0]
  
  # Compute intensity on observed exceedances
  intensity_0_obs <- exceedances_ds$dsmaxdt * exp(gamma_0 + gamma_1 * exceedances_ds$ICS_max)
  sigma_var_0_obs <- thresh_fit$par[1] + thresh_fit$par[2] * exceedances_ds$ICS_max
  intensity_above_threshold_obs <- intensity_0_obs * (1 + thresh_fit$par[3] * (threshold_ex_ds) / sigma_var_0_obs)^(-1 / thresh_fit$par[3])
  
  # Compute log-likelihood term 1
  LL_1 <- sum(log(intensity_above_threshold_obs))
  
  # Compute integrated intensity across space and time
  intensity_0 <- covariates$dsmaxdt * exp(gamma_0 + gamma_1 * covariates$ICS_max)
  sigma_var_0 <- thresh_fit$par[1] + thresh_fit$par[2] * covariates$ICS_max
  intensity_above_threshold <- intensity_0 * (1 + thresh_fit$par[3] * (covariates_threshold) / sigma_var_0)^(-1 / thresh_fit$par[3])
  intensity_above_threshold[1 + thresh_fit$par[3] * (covariates_threshold) / sigma_var_0 < 0] <- 0

  grid_box_E <- grid_box(covariates$Easting)
  grid_box_N <- grid_box(covariates$Northing)
  grid_box_area <- grid_box_E / 1000 * grid_box_N / 1000
  integrated_intensity <- sum(intensity_above_threshold, na.rm = TRUE) * grid_box_area
  
  # Compute log-likelihood term 2
  LL_2 <- -integrated_intensity
  
  return(LL_1 + LL_2)
}

#' Function for intensity estimation above the conservative or other constant threshold
#' 
#' @param par Numeric vector with parameter values (gamma_0, gamma_1)
#' @param data Data frame with earthquake catalog including observed locations/magnitudes/covariates
#' @param covariates Data frame with covariates (ICS, threshold, V1-V4, etc.) across whole space and time
#' @param threshold Numeric value of constant threshold
#' 
#' @return Numeric value of Poisson process log-likelihood
Poisson_process_LL_const_thresh <- function(par, data, covariates, threshold=1.45){

  # Extract parameters
  gamma_0 <- par[1]
  gamma_1 <- par[2]
  
  # earthquakes exceeding threshild magnitude at locations where dsmaxdt > 0
  exceedances <- data[data$Magnitude > threshold, ]
  exceedances_ds <- exceedances[exceedances$dsmaxdt > 0, ]  
  
  # Compute intensity on observed exceedances
  intensity_thresh_obs <- exceedances_ds$dsmaxdt * exp(gamma_0 + gamma_1 * exceedances_ds$ICS_max)
  
  # Compute log-likelihood term 1
  LL_1 <- sum(log(intensity_thresh_obs))
  
  # Compute integrated intensity across space and time
  intensity_thresh <- covariates$dsmaxdt * exp(gamma_0 + gamma_1 * covariates$ICS_max)
  
  grid_box_E <- grid_box(covariates$Easting)
  grid_box_N <- grid_box(covariates$Northing)
  grid_box_area <- grid_box_E / 1000 * grid_box_N / 1000
  
  integrated_intensity <- sum(intensity_thresh) * grid_box_area
  
  # Compute log-likelihood term 2
  LL_2 <- -integrated_intensity
  
  return(LL_1 + LL_2)
}

#' Function to compute resulting intensity on covariates grid
#' 
#' @param opt_PP List with optimization output of Poisson process intensity estimation
#' @param covariates Data frame with covariates (ICS, threshold, V1-V4, etc.) across whole space and time
#' @param covariates_threshold Numeric vector with threshold values at covariates grid points
#' @param thresh_fit List with threshold selection object GPD and threshold parameters
#' 
#' @return Numeric vector with resulting intensity values at covariates grid points
resulting_intensity_icsmax <- function(opt_PP, covariates, covariates_threshold, thresh_fit){
  
  intensity_0 <- covariates$dsmaxdt * exp(opt_PP$par[1] + opt_PP$par[2] * covariates$ICS_max)
  
  sigma_var_0 <- thresh_fit$par[1] + thresh_fit$par[2] * covariates$ICS_max
  
  intensity_above_threshold <- intensity_0 * (1 + thresh_fit$par[3] * (covariates_threshold) / sigma_var_0)^(-1 / thresh_fit$par[3])
  
  intensity_above_threshold[1 + thresh_fit$par[3] * (covariates_threshold) / sigma_var_0 < 0] <- 0
  return(intensity_above_threshold)
}


