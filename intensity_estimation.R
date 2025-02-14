
Poisson_process_LL <- function(par, data, covariates, threshold_obs, thresh_fit){
  # par: parameter vector
  # data: gron_eq_cat including observed locations/mags/covariates
  # covariates: ICS, threshold, V1-V4, etc. across space and time
  # threshold_obs: threshold for observed exceedances (done outside LL so 
            # that I can specify whether it is based on V1/2/3/4 before fitting LL)
  # thresh_fit: threshold selection object
  
  # Extract parameters
  gamma_1 <- par[1]
  gamma_2 <- par[2]
  
  exceedances <- data[data$Magnitude > threshold_obs, ]
  threshold_ex <- threshold_obs[data$Magnitude > threshold_obs]
  
  # Compute intensity on observed exceedances
  intensity_0_obs <- exp(gamma_1 + gamma_2 * exceedances$ICS)
  sigma_var_0_obs <- thresh_fit$par[1] + thresh_fit$par[2] * exceedances$ICS
  intensity_above_threshold_obs <- intensity_0_obs * (1 + thresh_fit$par[3] * (threshold_ex) / sigma_var_0_obs)^(-1 / thresh_fit$par[3])
  
  # Compute log-likelihood term 1
  LL_1 <- sum(log(intensity_above_threshold_obs))
  
  # Compute integrated intensity across space and time
  intensity_0 <- exp(gamma_1 + gamma_2 * covariates$ICS)
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


(opt_PP <- optim(Poisson_process_LL, par = c(0.1, 0), data = gron_eq_cat, covariates = covariates, threshold_obs = threshold_obs_V1, thresh_fit = thresh_fit_V1_ics, control=list(fnscale=-1)))

# Resulting intensity
resulting_intensity <- function(opt_PP, covariates, thresh_fit){
  intensity_0 <- exp(opt_PP$par[1] + opt_PP$par[2] * covariates$ICS)
  sigma_var_0 <- thresh_fit$par[1] + thresh_fit$par[2] * covariates$ICS
  intensity_above_threshold <- intensity_0 * (1 + thresh_fit$par[3] * (covariates$threshold) / sigma_var_0)^(-1 / thresh_fit$par[3])
  intensity_above_threshold[1 + thresh_fit$par[3] * (covariates$threshold) / sigma_var_0 < 0] <- 0
  return(intensity_above_threshold)
}

covariates$intensity_above_threshold_V1 <- resulting_intensity(opt_PP, covariates, thresh_fit_V1_ics)

# write.csv(covariates, "Data/covariates/covariates_1995-2024.csv", row.names = F) 
