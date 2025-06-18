

Poisson_process_LL_icsmax <- function(par, data, covariates, threshold_obs, covariates_threshold, thresh_fit){
  # par: parameter vector
  # data: gron_eq_cat including observed locations/mags/covariates
  # covariates: ICS, threshold, V1-V4, etc. across space and time
  # threshold_obs: threshold for observed exceedances (done outside LL so 
  # that I can specify whether it is based on V1/2/3/4 before fitting LL)
  # covariates_threshold: threshold value at covariates grid points
  # thresh_fit: threshold selection object
  
  # Extract parameters
  gamma_0 <- par[1]
  gamma_1 <- par[2]
  
  exceedances <- data[data$Magnitude > threshold_obs, ]
  threshold_ex <- threshold_obs[data$Magnitude > threshold_obs]
  
  # Removing locations where dsmaxdt= 0
  exceedances_ds <- exceedances[exceedances$dsmaxdt > 0, ]
  threshold_ex_ds <- threshold_ex[exceedances$dsmaxdt > 0]
  
  # Compute intensity on observed exceedances
  intensity_0_obs <- exceedances_ds$dsmaxdt*exp(gamma_0 + gamma_1 * exceedances_ds$ICS_max)
  sigma_var_0_obs <- thresh_fit$par[1] + thresh_fit$par[2] * exceedances_ds$ICS_max
  intensity_above_threshold_obs <- intensity_0_obs * (1 + thresh_fit$par[3] * (threshold_ex_ds) / sigma_var_0_obs)^(-1 / thresh_fit$par[3])
  
  # Compute log-likelihood term 1
  LL_1 <- sum(log(intensity_above_threshold_obs))
  
  # Compute integrated intensity across space and time
  intensity_0 <- covariates$dsmaxdt*exp(gamma_0 + gamma_1 * covariates$ICS_max)
  sigma_var_0 <- thresh_fit$par[1] + thresh_fit$par[2] * covariates$ICS_max
  intensity_above_threshold <- intensity_0 * (1 + thresh_fit$par[3] * (covariates_threshold) / sigma_var_0)^(-1 / thresh_fit$par[3])
  intensity_above_threshold[1 + thresh_fit$par[3] * (covariates_threshold) / sigma_var_0 < 0] <- 0
  grid_box_E <- (max(unique(covariates$Easting))-min(unique(covariates$Easting)))/length(unique(covariates$Easting))
  grid_box_N <- (max(unique(covariates$Northing))-min(unique(covariates$Northing)))/length(unique(covariates$Northing))
  grid_box_area <- grid_box_E/1000* grid_box_N/1000
  integrated_intensity <- sum(intensity_above_threshold, na.rm=T)*grid_box_area
  
  
  # Compute log-likelihood term 2
  LL_2 <- -integrated_intensity
  
  # Compute log-likelihood
  LL <- LL_1 + LL_2
  
  return(LL)
}

# Function for intensity estimation above the conservative or other constant threshold

Poisson_process_LL_const_thresh <- function(par, data, covariates, threshold=1.45){
  # par: parameter vector
  # data: gron_eq_cat including observed locations/mags/covariates
  # covariates: ICS, threshold, V1-V4, etc. across space and time
  # threshold_obs: threshold for observed exceedances (done outside LL so 
  # that I can specify whether it is based on V1/2/3/4 before fitting LL)
  # thresh_fit: threshold selection object
  
  # Extract parameters
  gamma_0 <- par[1]
  gamma_1 <- par[2]
  
  exceedances <- data[data$Magnitude > threshold, ]
  exceedances_ds <- exceedances[exceedances$dsmaxdt > 0, ]  # Removing locations where dsmaxdt = 0
  
  # Compute intensity on observed exceedances
  intensity_thresh_obs <- exceedances_ds$dsmaxdt*exp(gamma_0 + gamma_1 * exceedances_ds$ICS_max)
  
  # Compute log-likelihood term 1
  LL_1 <- sum(log(intensity_thresh_obs))
  
  # Compute integrated intensity across space and time
  intensity_thresh <- covariates$dsmaxdt*exp(gamma_0 + gamma_1 * covariates$ICS_max)
  
  grid_box_E <- (max(unique(covariates$Easting))-min(unique(covariates$Easting)))/length(unique(covariates$Easting))
  grid_box_N <- (max(unique(covariates$Northing))-min(unique(covariates$Northing)))/length(unique(covariates$Northing))
  grid_box_area <- grid_box_E/1000* grid_box_N/1000
  
  integrated_intensity <- sum(intensity_thresh)*grid_box_area
  
  
  # Compute log-likelihood term 2
  LL_2 <- -integrated_intensity
  
  # Compute log-likelihood
  LL <- LL_1 + LL_2
  
  return(LL)
}

resulting_intensity_icsmax <- function(opt_PP, covariates, covariates_threshold, thresh_fit){
  intensity_0 <- covariates$dsmaxdt*exp(opt_PP$par[1] + opt_PP$par[2] * covariates$ICS_max)
  sigma_var_0 <- thresh_fit$par[1] + thresh_fit$par[2] * covariates$ICS_max
  intensity_above_threshold <- intensity_0 * (1 + thresh_fit$par[3] * (covariates_threshold) / sigma_var_0)^(-1 / thresh_fit$par[3])
  intensity_above_threshold[1 + thresh_fit$par[3] * (covariates_threshold) / sigma_var_0 < 0] <- 0
  return(intensity_above_threshold)
}


# Old intensity functions -------------------------------------------------


Poisson_process_LL1 <- function(par, data, covariates, threshold_obs, thresh_fit){
  # par: parameter vector
  # data: gron_eq_cat including observed locations/mags/covariates
  # covariates: ICS, threshold, V1-V4, etc. across space and time
  # threshold_obs: threshold for observed exceedances (done outside LL so 
  # that I can specify whether it is based on V1/2/3/4 before fitting LL)
  # thresh_fit: threshold selection object
  
  # Extract parameters
  gamma_0 <- par[1]
  gamma_1 <- par[2]
  
  exceedances <- data[data$Magnitude > threshold_obs, ]
  threshold_ex <- threshold_obs[data$Magnitude > threshold_obs]
  
  # Compute intensity on observed exceedances
  intensity_0_obs <- exp(gamma_0 + gamma_1 * exceedances$ICS)
  sigma_var_0_obs <- thresh_fit$par[1] + thresh_fit$par[2] * exceedances$ICS
  intensity_above_threshold_obs <- intensity_0_obs * (1 + thresh_fit$par[3] * (threshold_ex) / sigma_var_0_obs)^(-1 / thresh_fit$par[3])
  
  # Compute log-likelihood term 1
  LL_1 <- sum(log(intensity_above_threshold_obs))
  
  # Compute integrated intensity across space and time
  intensity_0 <- exp(gamma_0 + gamma_1 * covariates$ICS)
  sigma_var_0 <- thresh_fit$par[1] + thresh_fit$par[2] * covariates$ICS
  intensity_above_threshold <- intensity_0 * (1 + thresh_fit$par[3] * (covariates$threshold) / sigma_var_0)^(-1 / thresh_fit$par[3])
  intensity_above_threshold[1 + thresh_fit$par[3] * (covariates$threshold) / sigma_var_0 < 0] <- 0
  grid_box_E <- (max(unique(covariates$Easting))-min(unique(covariates$Easting)))/length(unique(covariates$Easting))
  grid_box_N <- (max(unique(covariates$Northing))-min(unique(covariates$Northing)))/length(unique(covariates$Northing))
  grid_box_area <- grid_box_E/1000* grid_box_N/1000
  integrated_intensity <- sum(intensity_above_threshold)*grid_box_area
  
  
  # Compute log-likelihood term 2
  LL_2 <- -integrated_intensity
  
  # Compute log-likelihood
  LL <- LL_1 + LL_2
  
  return(LL)
}

Poisson_process_LL2 <- function(par, data, covariates, threshold_obs, thresh_fit){
  # par: parameter vector
  # data: gron_eq_cat including observed locations/mags/covariates
  # covariates: ICS, threshold, V1-V4, etc. across space and time
  # threshold_obs: threshold for observed exceedances (done outside LL so 
  # that I can specify whether it is based on V1/2/3/4 before fitting LL)
  # thresh_fit: threshold selection object
  
  # Extract parameters
  gamma_0 <- par[1]
  gamma_1 <- par[2]
  gamma_2 <- par[3]
  
  exceedances <- data[data$Magnitude > threshold_obs, ]
  threshold_ex <- threshold_obs[data$Magnitude > threshold_obs]
  
  # Compute intensity on observed exceedances
  intensity_0_obs <- max(exceedances$dsdt,gamma_2)*exp(gamma_0 + gamma_1 * exceedances$ICS)
  sigma_var_0_obs <- thresh_fit$par[1] + thresh_fit$par[2] * exceedances$ICS
  intensity_above_threshold_obs <- intensity_0_obs * (1 + thresh_fit$par[3] * (threshold_ex) / sigma_var_0_obs)^(-1 / thresh_fit$par[3])
  
  # Compute log-likelihood term 1
  LL_1 <- sum(log(intensity_above_threshold_obs))
  
  # Compute integrated intensity across space and time
  intensity_0 <- max(covariates$dsdt,gamma_2)*exp(gamma_0 + gamma_1 * covariates$ICS)
  sigma_var_0 <- thresh_fit$par[1] + thresh_fit$par[2] * covariates$ICS
  intensity_above_threshold <- intensity_0 * (1 + thresh_fit$par[3] * (covariates$threshold) / sigma_var_0)^(-1 / thresh_fit$par[3])
  intensity_above_threshold[1 + thresh_fit$par[3] * (covariates$threshold) / sigma_var_0 < 0] <- 0
  grid_box_E <- (max(unique(covariates$Easting))-min(unique(covariates$Easting)))/length(unique(covariates$Easting))
  grid_box_N <- (max(unique(covariates$Northing))-min(unique(covariates$Northing)))/length(unique(covariates$Northing))
  grid_box_area <- grid_box_E/1000* grid_box_N/1000
  integrated_intensity <- sum(intensity_above_threshold)*grid_box_area
  
  
  # Compute log-likelihood term 2
  LL_2 <- -integrated_intensity
  
  # Compute log-likelihood
  LL <- LL_1 + LL_2
  
  return(LL)
}

resulting_intensity_LL1 <- function(opt_PP, covariates, thresh_fit){
  intensity_0 <- exp(opt_PP$par[1] + opt_PP$par[2] * covariates$ICS)
  sigma_var_0 <- thresh_fit$par[1] + thresh_fit$par[2] * covariates$ICS
  intensity_above_threshold <- intensity_0 * (1 + thresh_fit$par[3] * (covariates$threshold) / sigma_var_0)^(-1 / thresh_fit$par[3])
  intensity_above_threshold[1 + thresh_fit$par[3] * (covariates$threshold) / sigma_var_0 < 0] <- 0
  return(intensity_above_threshold)
}

resulting_intensity_LL2 <- function(opt_PP, covariates, thresh_fit){
  intensity_0 <- max(covariates$dsdt, opt_PP$par[3])*exp(opt_PP$par[1] + opt_PP$par[2] * covariates$ICS)
  sigma_var_0 <- thresh_fit$par[1] + thresh_fit$par[2] * covariates$ICS
  intensity_above_threshold <- intensity_0 * (1 + thresh_fit$par[3] * (covariates$threshold) / sigma_var_0)^(-1 / thresh_fit$par[3])
  intensity_above_threshold[1 + thresh_fit$par[3] * (covariates$threshold) / sigma_var_0 < 0] <- 0
  return(intensity_above_threshold)
}

resulting_intensity_const_thresh <- function(opt_PP, covariates){
  intensity_thresh <- max(covariates$dsdt, opt_PP$par[3])*exp(opt_PP$par[1] + opt_PP$par[2] * covariates$ICS)
  return(intensity_thresh)
}
