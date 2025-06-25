

# Alg3: Model selection for non-parametric bootstraps --------------------------


# # Functions
# 
# transform_to_exp <- function (y, sig, xi){
#   std_exp <- (1 / xi) * log( 1 + xi * (y/sig))  
#   return(std_exp)
# }
# 
# GPD_LL_given_V_ICS <- function(par, excess, thresh_par, V, ics){
# 
#   if(length(par)!=3) stop("par must be a vector of length 3")
#   if(length(thresh_par)!=2) stop("thresh must be a vector of length 2")
#   if (!is.numeric(excess)) stop("excess must be a vector")
#   if (!is.numeric(V)) stop("V must be vector")
#   if (!is.numeric(ics)) stop("ics must be vector")
#   if(length(excess) != length(V)) stop("excess and V must be the same length")
#   if(length(excess) != length(ics)) stop("excess and ics must be the same length")
# 
#   sigma_par <- par[1:2]
#   xi <- par[3]
# 
#   sigma_ics <- sigma_par[1] + sigma_par[2]*ics
# 
#   sigma_tilde <- sigma_ics + xi*(thresh_par[[1]] + thresh_par[[2]]*V)
# 
#   sigma_check <- c(sigma_ics, sigma_tilde)
# 
#   if(all(sigma_check > 0) & xi > -0.9){
#     if(abs(xi) < 1e-10){
#       return(-sum(log(sigma_tilde)) - sum(excess/sigma_tilde))
#     }
#     else{
#       if(all(1+(xi*excess)/sigma_tilde > 0)){
#         return(-sum(log(sigma_tilde)) - (1+1/xi)*(sum(log(1+(xi*excess)/sigma_tilde))))
#       }
#       else{
#         return(-1e6)
#       }
#     }
#   }
#   else{
#     return(-1e7)
#   }
# }
# 
# 
# eqd_geo_ics <- function(data, thresh, distance_to_geo, ics, k = 100, m = 500, underlying_thresh = 0){
#   
#   # Check inputs are valid
#   if (!is.numeric(data)) stop("Data must be a vector")
#   if (!is.matrix(thresh)) stop("Thresholds must be a matrix")
#   if (!is.numeric(distance_to_geo)) stop("Third nearest distances must be a vector")
#   if (length(distance_to_geo) != length(data)) stop("Length of third nearest distances must match length of data")
#   if (k <= 0 | k %% 1 != 0) stop("Number of bootstrapped samples must be a positive integer")
#   if (m <= 0 | m %% 1 != 0) stop("Number of equally spaced probabilities must be a positive integer")
#   
#   meandistances <- xis <- num_excess <- numeric(nrow(thresh))
#   sigma_par <- matrix(NA, nrow = nrow(thresh), ncol = 2)
#   ########################## extra code for Inf investigation ##################
#   # inf_list_mean <- list()
#   # inf_list_boot <- list()
#   ########################## extra code for Inf investigation ##################
#   for (i in 1:nrow(thresh)) {
#     #print(i)
#     u <- thresh[i,1] + thresh[i,2] * distance_to_geo
#     #if(any(u < underlying_thresh)) stop(paste("Candidate thresholds must be above the underlying threshold of ", underlying_thresh))
#     if(any(u < underlying_thresh)) {
#       meandistances[i] <- NA
#       next
#       }
#     excess <- data[data > u] - u[data > u]
#     V_excess <- distance_to_geo[data > u]
#     ICS_excess <- ics[data > u]
#     num_excess[i] <- length(excess)
#     if (num_excess[i] > 20) {
#       mle0 <- mean(excess)
#       init.fit <- optim(GPD_LL_given_V_ICS, excess = excess, par = c(mle0, 0, 0.1), control = list(fnscale = -1), 
#                         thresh_par=thresh[i,], V = V_excess, ics=ICS_excess)
#       xis[i] <- init.fit$par[3]
#       sigma_par[i,] <- init.fit$par[1:2]
#       distances <- numeric(k)
#       j <- 1
#       while(j  <= k) {
#         #possible way to bootstrap
#         sampled_indices <- sample(1:num_excess[i], num_excess[i], replace = TRUE)
#         excess_boot <- excess[sampled_indices]
#         V_excess_boot <- V_excess[sampled_indices]
#         ICS_excess_boot <- ICS_excess[sampled_indices]
#         init_mle <- mean(excess_boot)
#         ifelse(xis[i] < 0, pars_init <-  c(init_mle, 0, 0.1) ,pars_init <- c(sigma_par[i,], xis[i]) )
#         gpd.fit <- optim(GPD_LL_given_V_ICS, excess = excess_boot, par = pars_init, control = list(fnscale = -1), 
#                          thresh_par=thresh[i,], V = V_excess_boot, ics = ICS_excess_boot)
#         sigma_given_V <- gpd.fit$par[1] + gpd.fit$par[2]*ICS_excess_boot + gpd.fit$par[3] * (thresh[i,1] + thresh[i,2]*V_excess_boot)
#         if(any((1 + gpd.fit$par[3] * excess_boot/sigma_given_V) < 1e-323)) next
#         else{
#           transformed_excess_boot <- transform_to_exp(y = excess_boot, sig = sigma_given_V, xi = gpd.fit$par[3])
#           sample_quants <- quantile(transformed_excess_boot, probs = (1:m) / (m+1))
#           model_quants <- qexp((1:m) / (m+1), rate = 1)
#           distances[j] <- mean(abs(sample_quants - model_quants))
#           ########################## extra code for Inf investigation ##########
#           # if(distances[j] == Inf){
#           #   inf_list_current <- list(i = i, j = j, excess_boot = excess_boot, transformed_excess = transformed_excess_boot, est = gpd.fit$par, sigma_var = sigma_given_V)
#           #   inf_list_boot <- c(inf_list_boot, list(inf_list_current))
#           # }
#           ########################## extra code for Inf investigation ##########
#           j <- j + 1
#         }
#       }
#       meandistances[i] <- mean(distances)
#       ########################## extra code for Inf investigation ##############
#       # if(meandistances[i] == Inf){
#       #   inf_list_mean <- c(inf_list_mean, list(list(i = i, excess = excess, par_ests = init.fit$par, dists = distances)))
#       # }
#       ########################## extra code for Inf investigation ##############
#     }
#     else{
#       meandistances[i] <- NA
#     }
#   }
#   chosen_index <- which.min(meandistances)
#   chosen_threshold_par <- thresh[chosen_index,]
#   xi <- xis[chosen_index]
#   sigmas <- sigma_par[chosen_index,]
#   len <- num_excess[chosen_index]
#   result <- list(thresh_par = chosen_threshold_par, par = c(sigmas,xi), num_excess = len, dists = meandistances)
#   ########################## extra code for Inf investigation ##################
#   # result <- list(thresh_par = chosen_threshold_par, par = c(sigmas,xi), num_excess = len, dists = meandistances, inf_list_mean = inf_list_mean, inf_list_boot = inf_list_boot)
#   ########################## extra code for Inf investigation ##################
#   return(result)
# }
# 
# threshold_selection_varying_formulation <- function(data, threshold_matrix, num_boot=200, SEED=11111, output_all = FALSE) {
#   message(Sys.time(),"Starting selection")
#   nearest_dist_matrix <- matrix(c(data$V_1, data$V_2, data$V_3, data$V_4), ncol=4, byrow=F)
#   results <- vector("list", 12)
#   for (ii in c(1:4)) {
#     nearest_dist <- nearest_dist_matrix[,ii]
#     log_nearest_dist <- log(nearest_dist)
#     sqrt_nearest_dist <- sqrt(nearest_dist)
#     
#     # V
#     set.seed(SEED)
#     results[[ii]] <- eqd_geo_ics(data=data$Magnitude, thresh = threshold_matrix, distance_to_geo = nearest_dist, k=num_boot, ics = data$ICS_max)
#     
#     # log(V)
#     set.seed(SEED)
#     results[[ii+4]] <- eqd_geo_ics(data=data$Magnitude, thresh = threshold_matrix, distance_to_geo =  log_nearest_dist, k=num_boot, ics = data$ICS_max)
#     
#     # sqrt(V)
#     set.seed(SEED)
#     results[[ii+8]] <- eqd_geo_ics(data=data$Magnitude, thresh = threshold_matrix, distance_to_geo = sqrt_nearest_dist, k=num_boot, ics = data$ICS_max)
#   }
#   chosen_model_index <- which.min(lapply(results, function(x){
#     if (is.null(x) || is.null(x$dists) || all(is.na(x$dists))) return(Inf)
#     min(x$dists, na.rm = T)
#     } ))
#   chosen_model <- results[[chosen_model_index]]
#   if (output_all) {
#     return(list(chosen_form = chosen_model_index, model_results=results))
#   }
#   else {
#     return(list(chosen_form = chosen_model_index, model_results=chosen_model))
#   }
# }
# 
# non_parametric_samples <- readRDS("STORM_input/non_parametric_samples.rds")
# 
# 
# intercepts <- seq(0, 1.5, by=0.02)
# slopes <- seq(0,0.8, by=0.02)
# threshold_matrix <- as.matrix(expand.grid(intercepts, slopes))
# 
# library(parallel)
# 
# cl <- makeCluster(100)
# clusterExport(cl, varlist = c("threshold_selection_varying_formulation", 
# "threshold_matrix", "eqd_geo_ics", "GPD_LL_given_V_ICS", "transform_to_exp"))
# clusterSetRNGStream(cl, 11111)
# 
# bootstrap_results <- parLapply(cl, non_parametric_samples, function(sample) {
#   threshold_selection_varying_formulation(
#     data = sample,
#     threshold_matrix = threshold_matrix,
#     num_boot = 100,
#     SEED = 11111,
#     output_all = FALSE
#   )
# })
# stopCluster(cl)
# 
# saveRDS(bootstrap_results, file="STORM_output/model_uncertainty_results.rds")

# Alg2 and 3: Intensity estimation for bootstrapped samples --------------------------

# For Alg 2, use threshold_uncertainty results and fix chosen form
# For Alg 3, use model_uncertainty results and extract chosen form from each

#Functions
# Poisson_process_LL_icsmax <- function(par, data, covariates, threshold_obs, covariates_threshold, thresh_fit){
#   # par: parameter vector
#   # data: gron_eq_cat including observed locations/mags/covariates
#   # covariates: ICS, threshold, V1-V4, etc. across space and time
#   # threshold_obs: threshold for observed exceedances (done outside LL so
#   # that I can specify whether it is based on V1/2/3/4 before fitting LL)
#   # covariates_threshold: threshold value at covariates grid points
#   # thresh_fit: threshold selection object
# 
#   # Extract parameters
#   gamma_0 <- par[1]
#   gamma_1 <- par[2]
# 
#   exceedances <- data[data$Magnitude > threshold_obs, ]
#   threshold_ex <- threshold_obs[data$Magnitude > threshold_obs]
# 
#   # Removing locations where dsmaxdt= 0
#   exceedances_ds <- exceedances[exceedances$dsmaxdt > 0, ]
#   threshold_ex_ds <- threshold_ex[exceedances$dsmaxdt > 0]
# 
#   # Compute intensity on observed exceedances
#   intensity_0_obs <- exceedances_ds$dsmaxdt*exp(gamma_0 + gamma_1 * exceedances_ds$ICS_max)
#   sigma_var_0_obs <- thresh_fit$par[1] + thresh_fit$par[2] * exceedances_ds$ICS_max
#   intensity_above_threshold_obs <- intensity_0_obs * (1 + thresh_fit$par[3] * (threshold_ex_ds) / sigma_var_0_obs)^(-1 / thresh_fit$par[3])
# 
#   # Compute log-likelihood term 1
#   LL_1 <- sum(log(intensity_above_threshold_obs))
# 
#   # Compute integrated intensity across space and time
#   intensity_0 <- covariates$dsmaxdt*exp(gamma_0 + gamma_1 * covariates$ICS_max)
#   sigma_var_0 <- thresh_fit$par[1] + thresh_fit$par[2] * covariates$ICS_max
#   intensity_above_threshold <- intensity_0 * (1 + thresh_fit$par[3] * (covariates_threshold) / sigma_var_0)^(-1 / thresh_fit$par[3])
#   intensity_above_threshold[1 + thresh_fit$par[3] * (covariates_threshold) / sigma_var_0 < 0] <- 0
#   grid_box_E <- (max(unique(covariates$Easting))-min(unique(covariates$Easting)))/length(unique(covariates$Easting))
#   grid_box_N <- (max(unique(covariates$Northing))-min(unique(covariates$Northing)))/length(unique(covariates$Northing))
#   grid_box_area <- grid_box_E/1000* grid_box_N/1000
#   integrated_intensity <- sum(intensity_above_threshold, na.rm = TRUE)*grid_box_area
# 
# 
#   # Compute log-likelihood term 2
#   LL_2 <- -integrated_intensity
# 
#   # Compute log-likelihood
#   LL <- LL_1 + LL_2
# 
#   return(LL)
# }
# 
# resulting_intensity_icsmax <- function(opt_PP_par, covariates, covariates_threshold, thresh_fit){
#   intensity_0 <- covariates$dsmaxdt*exp(opt_PP_par[1] + opt_PP_par[2] * covariates$ICS_max)
#   sigma_var_0 <- thresh_fit$par[1] + thresh_fit$par[2] * covariates$ICS_max
#   intensity_above_threshold <- intensity_0 * (1 + thresh_fit$par[3] * (covariates_threshold) / sigma_var_0)^(-1 / thresh_fit$par[3])
#   intensity_above_threshold[1 + thresh_fit$par[3] * (covariates_threshold) / sigma_var_0 < 0] <- 0
#   return(intensity_above_threshold)
# }
# 
# 
# intensity_estimation <- function(data, covariates, thresh_fit, distance_obs, distance_covariate, grid_box_area=0.2445286) {
#   threshold_obs <- thresh_fit$thresh_par[1] + thresh_fit$thresh_par[2] * distance_obs
#   covariates_threshold <- thresh_fit$thresh_par[1] + thresh_fit$thresh_par[2] * distance_covariate
# 
#   PP_fit <- optim(Poisson_process_LL_icsmax, par = c(0.1, 0), data = data, covariates = covariates,
#                   threshold_obs = threshold_obs, covariates_threshold = covariates_threshold,
#                   thresh_fit = thresh_fit, control=list(fnscale=-1))
# 
#   return(PP_fit$par)
# }
# 
# non_parametric_samples <- readRDS("STORM_input/non_parametric_samples.rds")
# covariates <- read.csv("STORM_input/covariates_1995-2024.csv", header=T)
# threshold_uncertainty_results <- readRDS("STORM_input/threshold_values_uncertainty_results_Alg2.rds")
# 
# covariates_distances <- cbind(covariates$V1, covariates$V2, covariates$V3, covariates$V4,
#                               log(covariates$V1), log(covariates$V2), log(covariates$V3), log(covariates$V4),
#                               sqrt(covariates$V1), sqrt(covariates$V2), sqrt(covariates$V3), sqrt(covariates$V4))
# 
# library(parallel)
# 
# cl <- makeCluster(100)  # Leave one core free
# clusterExport(cl, varlist = c("Poisson_process_LL_icsmax", "resulting_intensity_icsmax", "intensity_estimation",
#                               "covariates", "covariates_distances", "threshold_uncertainty_results", "non_parametric_samples"))
# 
# bootstrap_intensity_results <- parLapply(cl, seq_along(non_parametric_samples), function(i) {
#   sample <- non_parametric_samples[[i]]
#   dist_matrix <- cbind(sample$V_1, sample$V_2, sample$V_3, sample$V_4, log(sample$V_1),
#                        log(sample$V_2), log(sample$V_3), log(sample$V_4),
#                        sqrt(sample$V_1), sqrt(sample$V_2), sqrt(sample$V_3), sqrt(sample$V_4))
#   # model_results <- model_uncertainty_results[[i]]
#   thresh_fit <- threshold_uncertainty_results[[i]]
#   # chosen_model_index <- model_results$chosen_form
#   chosen_model_index <- 2 # If using Alg 2, fix to chosen model (V2)
#   # thresh_fit <- model_results$model_results
#   intensity_estimation(
#     data = sample,
#     covariates = covariates,
#     thresh_fit = thresh_fit,
#     distance_obs = dist_matrix[, chosen_model_index],
#     distance_covariate = covariates_distances[, chosen_model_index]
#   )
# })
# stopCluster(cl)
# 
# bootstrap_intensity_fits <- do.call(rbind, bootstrap_intensity_results)
# 
# saveRDS(bootstrap_intensity_fits, file="STORM_output/bootstrap_intensity_fits_Alg2.rds")
# 



# Estimates for Alg 1 on its own and with Alg 2/3 results----------------------------------------

# For Alg 2, use threshold_uncertainty results and fix chosen form
# For Alg 3, use model_uncertainty results and extract chosen form from each


# Intensity estimation --------------------------

# Functions
# Poisson_process_LL_icsmax <- function(par, data, covariates, threshold_obs, covariates_threshold, thresh_fit){
#   # par: parameter vector
#   # data: gron_eq_cat including observed locations/mags/covariates
#   # covariates: ICS, threshold, V1-V4, etc. across space and time
#   # threshold_obs: threshold for observed exceedances (done outside LL so
#   # that I can specify whether it is based on V1/2/3/4 before fitting LL)
#   # covariates_threshold: threshold value at covariates grid points
#   # thresh_fit: threshold selection object
# 
#   # Extract parameters
#   gamma_0 <- par[1]
#   gamma_1 <- par[2]
# 
#   exceedances <- data[data$Magnitude > threshold_obs, ]
#   threshold_ex <- threshold_obs[data$Magnitude > threshold_obs]
# 
#   # Removing locations where dsmaxdt= 0
#   exceedances_ds <- exceedances[exceedances$dsmaxdt > 0, ]
#   threshold_ex_ds <- threshold_ex[exceedances$dsmaxdt > 0]
# 
#   # Compute intensity on observed exceedances
#   intensity_0_obs <- exceedances_ds$dsmaxdt*exp(gamma_0 + gamma_1 * exceedances_ds$ICS_max)
#   sigma_var_0_obs <- thresh_fit$par[1] + thresh_fit$par[2] * exceedances_ds$ICS_max
#   intensity_above_threshold_obs <- intensity_0_obs * (1 + thresh_fit$par[3] * (threshold_ex_ds) / sigma_var_0_obs)^(-1 / thresh_fit$par[3])
# 
#   # Compute log-likelihood term 1
#   LL_1 <- sum(log(intensity_above_threshold_obs))
# 
#   # Compute integrated intensity across space and time
#   intensity_0 <- covariates$dsmaxdt*exp(gamma_0 + gamma_1 * covariates$ICS_max)
#   sigma_var_0 <- thresh_fit$par[1] + thresh_fit$par[2] * covariates$ICS_max
#   intensity_above_threshold <- intensity_0 * (1 + thresh_fit$par[3] * (covariates_threshold) / sigma_var_0)^(-1 / thresh_fit$par[3])
#   intensity_above_threshold[1 + thresh_fit$par[3] * (covariates_threshold) / sigma_var_0 < 0] <- 0
#   grid_box_E <- (max(unique(covariates$Easting))-min(unique(covariates$Easting)))/length(unique(covariates$Easting))
#   grid_box_N <- (max(unique(covariates$Northing))-min(unique(covariates$Northing)))/length(unique(covariates$Northing))
#   grid_box_area <- grid_box_E/1000* grid_box_N/1000
#   integrated_intensity <- sum(intensity_above_threshold, na.rm = TRUE)*grid_box_area
# 
# 
#   # Compute log-likelihood term 2
#   LL_2 <- -integrated_intensity
# 
#   # Compute log-likelihood
#   LL <- LL_1 + LL_2
# 
#   return(LL)
# }
# 
# resulting_intensity_icsmax <- function(opt_PP_par, covariates, covariates_threshold, thresh_fit){
#   intensity_0 <- covariates$dsmaxdt*exp(opt_PP_par[1] + opt_PP_par[2] * covariates$ICS_max)
#   sigma_var_0 <- thresh_fit$par[1] + thresh_fit$par[2] * covariates$ICS_max
#   intensity_above_threshold <- intensity_0 * (1 + thresh_fit$par[3] * (covariates_threshold) / sigma_var_0)^(-1 / thresh_fit$par[3])
#   intensity_above_threshold[1 + thresh_fit$par[3] * (covariates_threshold) / sigma_var_0 < 0] <- 0
#   return(intensity_above_threshold)
# }
# 
# #
# intensity_estimation <- function(data, covariates, thresh_fit, distance_obs, distance_covariate, grid_box_area=0.2445286) {
#   threshold_obs <- thresh_fit$thresh_par[1] + thresh_fit$thresh_par[2] * distance_obs
#   covariates_threshold <- thresh_fit$thresh_par[1] + thresh_fit$thresh_par[2] * distance_covariate
# 
#   PP_fit <- optim(Poisson_process_LL_icsmax, par = c(0.1, 0), data = data, covariates = covariates,
#                   threshold_obs = threshold_obs, covariates_threshold = covariates_threshold,
#                   thresh_fit = thresh_fit, control=list(fnscale=-1))
# 
#   return(PP_fit$par)
# }
# 
# # GPD fits --------------------------
# 
# GPD_LL_given_V_ICS <- function(par, excess, thresh_par, V, ics){
# 
#   if(length(par)!=3) stop("par must be a vector of length 3")
#   if(length(thresh_par)!=2) stop("thresh must be a vector of length 2")
#   if (!is.numeric(excess)) stop("excess must be a vector")
#   if (!is.numeric(V)) stop("V must be vector")
#   if (!is.numeric(ics)) stop("ics must be vector")
#   if(length(excess) != length(V)) stop("excess and V must be the same length")
#   if(length(excess) != length(ics)) stop("excess and ics must be the same length")
# 
#   sigma_par <- par[1:2]
#   xi <- par[3]
# 
#   sigma_ics <- sigma_par[1] + sigma_par[2]*ics
# 
#   sigma_tilde <- sigma_ics + xi*(thresh_par[[1]] + thresh_par[[2]]*V)
# 
#   sigma_check <- c(sigma_ics, sigma_tilde)
# 
#   if(all(sigma_check > 0) & xi > -0.9){
#     if(abs(xi) < 1e-10){
#       return(-sum(log(sigma_tilde)) - sum(excess/sigma_tilde))
#     }
#     else{
#       if(all(1+(xi*excess)/sigma_tilde > 0)){
#         return(-sum(log(sigma_tilde)) - (1+1/xi)*(sum(log(1+(xi*excess)/sigma_tilde))))
#       }
#       else{
#         return(-1e6)
#       }
#     }
#   }
#   else{
#     return(-1e7)
#   }
# }
# 
# 
# parametric_bootstrap_fits <- function(boot_sample, thresh_fit, covariates_distances, chosen_form) {
# 
# 
#   sampled_covariates <- boot_sample$sampled_covariates
#   V_form <- boot_sample$V_form
#   thresh_vals <- boot_sample$thresh_vals
#   boot_excesses <- boot_sample$boot_excesses
# 
#   # Fit GPD model
#   xi <- thresh_fit$par[3]
#   init_mle <- mean(boot_excesses)
#   pars_init <- if (xi < 0) c(init_mle, 0, 0.1) else thresh_fit$par
#   gpd_fit <- optim(
#     GPD_LL_given_V_ICS,
#     par = pars_init,
#     excess = boot_excesses,
#     thresh_par = thresh_fit$thresh_par,
#     V = V_form,
#     ics = sampled_covariates$ICS_max,
#     control = list(fnscale = -1)
#   )
# 
#   # Add magnitudes and re-fit intensity model
#   sampled_covariates$Magnitude <- boot_excesses + thresh_vals
#   intensity_boot <- intensity_estimation(
#     sampled_covariates, covariates, thresh_fit,
#     distance_obs = V_form,
#     distance_covariate = covariates_distances[, chosen_form]
#   )
# 
#   boot_estimates <- list(
#     fit = gpd_fit$par,
#     intensity_fit = intensity_boot
#   )
# 
#   return(boot_estimates)
# }
# 
# # Generation of parmetric samples for combining Alg 1 and Alg 3
# 
# rgpd <- function(n, shape, scale = NULL, nu = NULL, mu = 0){
#   ## Input checks
#   # one and only one of {nu, scale} may be specified
#   if (is.null(scale) & is.null(nu)) {
#     stop('Define one of the parameters nu or scale.')
#   }
#   if (!is.null(scale) & !is.null(nu)) {
#     stop('Define only one of the parameters nu and scale.')
#   }
#   # Calculate scale from nu if required
#   if (!is.null(nu) & is.null(scale)) {
#     scale <- nu / (1 + shape)
#     if (any(scale <= 0)) {
#       stop('Implied scale parameter(s) must be positive.')
#     }
#   }
#   # Check that scale value(s) are positive
#   if (any(scale <= 0)) {
#     stop('Scale parameter(s) must be positive.')
#   }
#   # Ensure q, scale, shape and mu are of same length.
#   if ((length(scale) == 1) & (n > 1)) {
#     scale <- rep(scale, n)
#   }
#   if ((length(shape) == 1) & (n > 1)) {
#     shape <- rep(shape, n)
#   }
#   if ((length(mu) == 1) & (n > 1)) {
#     mu <- rep(mu, n)
#   }
# 
#   #simulate sample
#   sample <- mu + (scale/shape) * ((1 - stats::runif(n))^(-shape) - 1)
#   #correct sample values where xi = 0
#   #ex <- which(shape ==0)
#   if (any(abs(shape) < 1e-10)) {
#     ex <- which(abs(shape) < 1e-10)
#     sample[ex] <- mu[ex] +
#       stats::rexp(n = length(ex),rate = 1/scale[ex])
#   }
#   return(sample)
# }
# 
# 
# generate_parametric_bootstrap_samples <- function(probs, aggregated_intensity, thresh_fit, covariates, chosen_form, grid_box_area=0.2445286, n_samples=200, SEED=11111) {
# 
#   n_cells <- length(probs)
#   boot_samples <- vector("list", n_samples)
#   set.seed(SEED)
#   for (i in seq_len(n_samples)) {
#     # 1. Sample number of events
#     n_obs <- rpois(1, aggregated_intensity)
#     sampled_indices <- sample.int(n_cells, n_obs, replace = TRUE, prob = probs)
#     sampled_covariates <- covariates[sampled_indices, , drop = FALSE]
# 
#     # 2. Calculate distances and threshold values
#     V_form <- covariates_distances[sampled_indices, chosen_form]
#     thresh_vals <- thresh_fit$thresh_par[1] + thresh_fit$thresh_par[2] * V_form
#     sigma_vals <- thresh_fit$par[1] + thresh_fit$par[2]*sampled_covariates$ICS_max + thresh_fit$par[3]*thresh_vals
# 
#     # 3. Simulate GPD excesses
#     xi <- thresh_fit$par[3]
#     boot_excesses <- rgpd(n_obs, scale = sigma_vals, shape = xi)
# 
#     boot_samples[[i]] <- list(
#       sampled_covariates = sampled_covariates,
#       V_form = V_form,
#       thresh_vals = thresh_vals,
#       boot_excesses = boot_excesses
#     )
#   }
# 
#   return(boot_samples)
# }


# With Alg 2 and 3  -----------------------


# parametric_samples <- readRDS("STORM_input/parametric_samples_Alg1.rds")
# covariates <- read.csv("STORM_input/covariates_1995-2024.csv", header=T)
# thresh_fit_A2 <- readRDS("STORM_input/geo_thresh_fit_V2.rds")
# bootstrap_intensity_results <- readRDS("STORM_input/bootstrap_intensity_fits_Alg3.rds") #Alg3
# bootstrap_intensity_results <- readRDS("STORM_input/bootstrap_intensity_fits_Alg2.rds") #Alg2
# bootstrap_model_results <- readRDS("STORM_input/bootstrap_model_selection_results_Alg3.rds") #Alg3
# bootstrap_model_results <- readRDS("STORM_input/threshold_values_uncertainty_results_Alg2.rds") #Alg2

# covariates_distances <- cbind(covariates$V1, covariates$V2, covariates$V3, covariates$V4,
#                               log(covariates$V1), log(covariates$V2), log(covariates$V3), log(covariates$V4),
#                               sqrt(covariates$V1), sqrt(covariates$V2), sqrt(covariates$V3), sqrt(covariates$V4))
# grid_box_area <- 0.2445286  # Area of grid box in km^2
# library(parallel)
# 
# bootstrap_fits <- vector("list", length(bootstrap_model_results))
# 
# for (j in seq_along(bootstrap_model_results)) {
#   # boot_model <- bootstrap_model_results[[j]] #Alg 3
#   # chosen_form <- boot_model$chosen_form # Alg 3
#   chosen_form <- 2 # If using Alg 2, fix to chosen model (V2)
#   # thresh_fit <- boot_model$model_results # Alg 3
#   thresh_fit <- bootstrap_model_results[[j]] # Alg 2
#   intensity_par <- bootstrap_intensity_results[j,]
#   covariates_threshold <- thresh_fit$thresh_par[1] + thresh_fit$thresh_par[2] * covariates_distances[, chosen_form]
# 
#   estimated_intensity <- resulting_intensity_icsmax(intensity_par, covariates, covariates_threshold, thresh_fit)
#   aggregated_intensity <- sum(estimated_intensity, na.rm = TRUE)*grid_box_area
#   probs <- estimated_intensity / aggregated_intensity
#   n_cells <- length(probs)
# 
# 
#   # Generate parametric bootstrap samples
#   bootstrap_samples_Alg1 <- generate_parametric_bootstrap_samples(
#     probs = probs,
#     aggregated_intensity = aggregated_intensity,
#     thresh_fit = thresh_fit,
#     covariates = covariates,
#     chosen_form = chosen_form,
#     n_samples = 200,
#     SEED = 11111
#   )
#   cl <- makeCluster(100)  # Match this with your SBATCH -c value minus 1 ideally
# 
#   clusterExport(cl, varlist = c("parametric_bootstrap_fits","GPD_LL_given_V_ICS",
#                                 "intensity_estimation", "Poisson_process_LL_icsmax","covariates",
#                                 "covariates_distances", "thresh_fit",
#                                 "chosen_form" ), envir = environment())
#   bootstrap_fits_Alg1 <- parLapply(cl, bootstrap_samples_Alg1, function(sample) {
#     # Generate parametric bootstrap samples
#     parametric_bootstrap_fits(
#       boot_sample = sample,
#       thresh_fit = thresh_fit,
#       covariates_distances = covariates_distances,
#       chosen_form = chosen_form
#     )
#   })
#   stopCluster(cl)
# 
#   bootstrap_fits[[j]] <- bootstrap_fits_Alg1
# 
# }
# 
# 
# saveRDS(bootstrap_fits, file="STORM_output/bootstrap_fits_Alg1and2.rds")


# # Alg 1 on its own -----------------------------
# 
# covariates <- read.csv("STORM_input/covariates_1995-2024.csv", header=T)
# thresh_fit <- readRDS("STORM_input/geo_thresh_fit_sqrtV2.rds")
# 
# observed_sample <- read.csv("STORM_input/unrounded_after_1995_in_polygon_with_covariates.csv", header=T)
# 
# covariates_distances <- cbind(covariates$V1, covariates$V2, covariates$V3, covariates$V4,
#                               log(covariates$V1), log(covariates$V2), log(covariates$V3), log(covariates$V4),
#                               sqrt(covariates$V1), sqrt(covariates$V2), sqrt(covariates$V3), sqrt(covariates$V4))
# grid_box_area <- 0.2445286  # Area of grid box in km^2
# 
# library(parallel)
# 
# chosen_form <- 10 # If using Alg 2, fix to chosen model (V2)
# 
# intensity_par <- intensity_estimation(observed_sample, covariates, thresh_fit, log(observed_sample$V_1), covariates_distances[,chosen_form])
# covariates_threshold <- thresh_fit$thresh_par[1] + thresh_fit$thresh_par[2] * covariates_distances[, chosen_form]
# 
# estimated_intensity <- resulting_intensity_icsmax(intensity_par, covariates, covariates_threshold, thresh_fit)
# aggregated_intensity <- sum(estimated_intensity, na.rm = TRUE)*grid_box_area
# probs <- estimated_intensity / aggregated_intensity
# 
# 
# # Generate parametric bootstrap samples
# bootstrap_samples_Alg1 <- generate_parametric_bootstrap_samples(
#   probs = probs,
#   aggregated_intensity = aggregated_intensity,
#   thresh_fit = thresh_fit,
#   covariates = covariates,
#   chosen_form = chosen_form,
#   n_samples = 200,
#   SEED = 11111
# )
# cl <- makeCluster(100)  # Match this with your SBATCH -c value minus 1 ideally
# 
# clusterExport(cl, varlist = c("parametric_bootstrap_fits","GPD_LL_given_V_ICS",
#                               "intensity_estimation", "Poisson_process_LL_icsmax","covariates",
#                               "covariates_distances", "thresh_fit",
#                               "chosen_form" ), envir = environment())
# bootstrap_fits_Alg1 <- parLapply(cl, bootstrap_samples_Alg1, function(sample) {
#   # Generate parametric bootstrap samples
#   parametric_bootstrap_fits(
#     boot_sample = sample,
#     thresh_fit = thresh_fit,
#     covariates_distances = covariates_distances,
#     chosen_form = chosen_form
#   )
# })
# stopCluster(cl)
# 
# 
# 
# saveRDS(bootstrap_fits_Alg1, file="STORM_output/bootstrap_fits_Alg1_C2.rds")
# 
# 

# Future inference --------------------------------------------------------



# # Functions to evaluate endpoints -------------
# evaluate_endpoints <- function(future_covariates, intensity_par, gpd_par, grid_box_area = 0.2445286) {
#   
#   if(gpd_par[3] >= 0) {
#     endpoint_summaries <- c(NA, NA)
#   }
#   else{
#     endpoints <- -(gpd_par[1] + gpd_par[2] * future_covariates$ICS_max)/gpd_par[3]
#   
#     endpoint_max <- max(endpoints, na.rm = TRUE)
#     
#     # weighted mean
#     intensity_above0 <- future_covariates$dsmaxdt * exp(intensity_par[1] + intensity_par[2] * future_covariates$ICS_max) 
#     
#     agg_intensity <- sum(intensity_above0, na.rm = TRUE) 
#     
#     # Calculate the weighted mean endpoint
#     endpoint_wm <- sum(endpoints * intensity_above0/agg_intensity , na.rm = TRUE)
#     
#     endpoint_summaries <- c(endpoint_max, endpoint_wm)
#   }  
#   return(endpoint_summaries)
# }
# 
# 
# # Functions for evaluated v level  ------------
# 
# # v-level extreme event estimation
# v_level_extreme_event <- function(v, prob, future_covariates, intensity_par, gpd_par, grid_box_area = 0.2445286) {
#   
#   intensity_above0 <- future_covariates$dsmaxdt * exp(intensity_par[1] + intensity_par[2] * future_covariates$ICS_max)
#   sigma0 <- gpd_par[1] + gpd_par[2] * future_covariates$ICS_max
#   intensity_above_v <-  intensity_above0 * (1 + gpd_par[3] * (v) /sigma0 )^(-1 / gpd_par[3])
#   # Do we need below?
#   intensity_above_v[1 + gpd_par[3] * (v) / sigma0 < 0] <- 0
#   # Calculate aggregated intensity
#   agg_intensity_v <- sum(intensity_above_v, na.rm = TRUE) * grid_box_area
#   
#   function_to_solve <- (exp(-agg_intensity_v) -prob)^2
#   
#   return(function_to_solve)
# }
# 
# # v_level solver function
# v_level_solver <- function(prob, future_covariates,  intensity_par, gpd_par, upper_limit = 10, grid_box_area = 0.2445286) {
#   
#   if(gpd_par[3] >= 0) {
#     upper_limit <- 13
#   }
#   else{
#     # Calculate the endpoints
#     endpoints <- -(gpd_par[1] + gpd_par[2] * future_covariates$ICS_max)/gpd_par[3]
#     
#     endpoint_max <- max(endpoints, na.rm = TRUE)
#     v_solution <- optimize(v_level_extreme_event, interval = c(0, endpoint_max), prob = 0.9387, 
#                            future_covariates = future_covariates, intensity_par = intensity_par, gpd_par=gpd_par)$minimum
#   }
#   return(v_solution)
# }
# 
# 
# # Main script to run the uncertainty algorithm for STORM -----------------
# future_covariates <- read.csv("STORM_input/future_covariates_2024-2055.csv", header=T)
# bootstrap_estimates_list <- readRDS("STORM_input/bootstrap_fits_Alg1and3.rds") #
# #bootstrap_estimates_list <- readRDS("STORM_input/bootstrap_fits_Alg1and2.rds")
# 
# library(parallel)
# 
# summaries <- matrix(NA, nrow = length(bootstrap_estimates_list)*200, ncol = 3)
# colnames(summaries) <- c("v_level", "endpoint_max", "endpoint_wm")
# 
# for (j in seq_along(bootstrap_estimates_list)) {
#   bootstrap_estimates <- bootstrap_estimates_list[[j]]
# 
#   cl <- makeCluster(100)  # Match this with your SBATCH -c value minus 1 ideally
# 
#   clusterExport(cl, varlist = c("future_covariates", "evaluate_endpoints", "v_level_solver", "v_level_extreme_event"), envir = environment())
#   
#   v_levels_vec <- unlist(parLapply(cl, bootstrap_estimates, function(current_fits) {
#     v_level_solver(prob = 0.9387,
#                    future_covariates = future_covariates,
#                    intensity_par = current_fits$intensity_fit,
#                    gpd_par = current_fits$fit)
#   }))
#   
#   endpoints_list <- parLapply(cl, bootstrap_estimates, function(current_fits) {
#     evaluate_endpoints(future_covariates, current_fits$intensity_fit, current_fits$fit)
#   })
#   
#   endpoints_mat <- do.call(rbind, endpoints_list)
#   
#   stopCluster(cl)
# 
#   summaries[((j-1)*200+1):(j*200), 1] <- v_levels_vec
#   summaries[((j-1)*200+1):(j*200), 2] <- endpoints_mat[, 1]  # max endpoint
#   summaries[((j-1)*200+1):(j*200), 3] <- endpoints_mat[, 2]  # weighted mean endpoint
# 
# }
# 
# 
# saveRDS(summaries, file="STORM_output/future_inferences_Alg3.rds")







# Alg 1 on its own with conservative threshold    -----------------------------

# Functions
# Poisson_process_LL_icsmax <- function(par, data, covariates, threshold_obs, covariates_threshold, thresh_fit){
#   # par: parameter vector
#   # data: gron_eq_cat including observed locations/mags/covariates
#   # covariates: ICS, threshold, V1-V4, etc. across space and time
#   # threshold_obs: threshold for observed exceedances (done outside LL so
#   # that I can specify whether it is based on V1/2/3/4 before fitting LL)
#   # covariates_threshold: threshold value at covariates grid points
#   # thresh_fit: threshold selection object
# 
#   # Extract parameters
#   gamma_0 <- par[1]
#   gamma_1 <- par[2]
# 
#   exceedances <- data[data$Magnitude > threshold_obs, ]
#   threshold_ex <- threshold_obs[data$Magnitude > threshold_obs]
# 
#   # Removing locations where dsmaxdt= 0
#   exceedances_ds <- exceedances[exceedances$dsmaxdt > 0, ]
#   threshold_ex_ds <- threshold_ex[exceedances$dsmaxdt > 0]
# 
#   # Compute intensity on observed exceedances
#   intensity_0_obs <- exceedances_ds$dsmaxdt*exp(gamma_0 + gamma_1 * exceedances_ds$ICS_max)
#   sigma_var_0_obs <- thresh_fit$par[1] + thresh_fit$par[2] * exceedances_ds$ICS_max
#   intensity_above_threshold_obs <- intensity_0_obs * (1 + thresh_fit$par[3] * (threshold_ex_ds) / sigma_var_0_obs)^(-1 / thresh_fit$par[3])
# 
#   # Compute log-likelihood term 1
#   LL_1 <- sum(log(intensity_above_threshold_obs))
# 
#   # Compute integrated intensity across space and time
#   intensity_0 <- covariates$dsmaxdt*exp(gamma_0 + gamma_1 * covariates$ICS_max)
#   sigma_var_0 <- thresh_fit$par[1] + thresh_fit$par[2] * covariates$ICS_max
#   intensity_above_threshold <- intensity_0 * (1 + thresh_fit$par[3] * (covariates_threshold) / sigma_var_0)^(-1 / thresh_fit$par[3])
#   intensity_above_threshold[1 + thresh_fit$par[3] * (covariates_threshold) / sigma_var_0 < 0] <- 0
#   grid_box_E <- (max(unique(covariates$Easting))-min(unique(covariates$Easting)))/length(unique(covariates$Easting))
#   grid_box_N <- (max(unique(covariates$Northing))-min(unique(covariates$Northing)))/length(unique(covariates$Northing))
#   grid_box_area <- grid_box_E/1000* grid_box_N/1000
#   integrated_intensity <- sum(intensity_above_threshold, na.rm = TRUE)*grid_box_area
# 
# 
#   # Compute log-likelihood term 2
#   LL_2 <- -integrated_intensity
# 
#   # Compute log-likelihood
#   LL <- LL_1 + LL_2
# 
#   return(LL)
# }
# 
# resulting_intensity_icsmax <- function(opt_PP_par, covariates, covariates_threshold, thresh_fit){
#   intensity_0 <- covariates$dsmaxdt*exp(opt_PP_par[1] + opt_PP_par[2] * covariates$ICS_max)
#   sigma_var_0 <- thresh_fit$par[1] + thresh_fit$par[2] * covariates$ICS_max
#   intensity_above_threshold <- intensity_0 * (1 + thresh_fit$par[3] * (covariates_threshold) / sigma_var_0)^(-1 / thresh_fit$par[3])
#   intensity_above_threshold[1 + thresh_fit$par[3] * (covariates_threshold) / sigma_var_0 < 0] <- 0
#   return(intensity_above_threshold)
# }
# 
# #
# intensity_estimation <- function(data, covariates, thresh_fit, distance_obs, distance_covariate, grid_box_area=0.2445286) {
#   threshold_obs <- thresh_fit$thresh_par[1] + thresh_fit$thresh_par[2] * distance_obs
#   covariates_threshold <- thresh_fit$thresh_par[1] + thresh_fit$thresh_par[2] * distance_covariate
# 
#   PP_fit <- optim(Poisson_process_LL_icsmax, par = c(0.1, 0), data = data, covariates = covariates,
#                   threshold_obs = threshold_obs, covariates_threshold = covariates_threshold,
#                   thresh_fit = thresh_fit, control=list(fnscale=-1))
# 
#   return(PP_fit$par)
# }
# 
# # GPD fits --------------------------
# 
# GPD_LL_given_V_ICS <- function(par, excess, thresh_par, V, ics){
# 
#   if(length(par)!=3) stop("par must be a vector of length 3")
#   if(length(thresh_par)!=2) stop("thresh must be a vector of length 2")
#   if (!is.numeric(excess)) stop("excess must be a vector")
#   if (!is.numeric(V)) stop("V must be vector")
#   if (!is.numeric(ics)) stop("ics must be vector")
#   if(length(excess) != length(V)) stop("excess and V must be the same length")
#   if(length(excess) != length(ics)) stop("excess and ics must be the same length")
# 
#   sigma_par <- par[1:2]
#   xi <- par[3]
# 
#   sigma_ics <- sigma_par[1] + sigma_par[2]*ics
# 
#   sigma_tilde <- sigma_ics + xi*(thresh_par[[1]] + thresh_par[[2]]*V)
# 
#   sigma_check <- c(sigma_ics, sigma_tilde)
# 
#   if(all(sigma_check > 0) & xi > -0.9){
#     if(abs(xi) < 1e-10){
#       return(-sum(log(sigma_tilde)) - sum(excess/sigma_tilde))
#     }
#     else{
#       if(all(1+(xi*excess)/sigma_tilde > 0)){
#         return(-sum(log(sigma_tilde)) - (1+1/xi)*(sum(log(1+(xi*excess)/sigma_tilde))))
#       }
#       else{
#         return(-1e6)
#       }
#     }
#   }
#   else{
#     return(-1e7)
#   }
# }
# 
# 
# parametric_bootstrap_fits <- function(boot_sample, thresh_fit, covariates_distances, chosen_form) {
# 
# 
#   sampled_covariates <- boot_sample$sampled_covariates
#   V_form <- boot_sample$V_form
#   thresh_vals <- boot_sample$thresh_vals
#   boot_excesses <- boot_sample$boot_excesses
# 
#   # Fit GPD model
#   xi <- thresh_fit$par[3]
#   init_mle <- mean(boot_excesses)
#   pars_init <- if (xi < 0) c(init_mle, 0, 0.1) else thresh_fit$par
#   gpd_fit <- optim(
#     GPD_LL_given_V_ICS,
#     par = pars_init,
#     excess = boot_excesses,
#     thresh_par = thresh_fit$thresh_par,
#     V = V_form,
#     ics = sampled_covariates$ICS_max,
#     control = list(fnscale = -1)
#   )
# 
#   # Add magnitudes and re-fit intensity model
#   sampled_covariates$Magnitude <- boot_excesses + thresh_vals
#   intensity_boot <- intensity_estimation(
#     sampled_covariates, covariates, thresh_fit,
#     distance_obs = V_form,
#     distance_covariate = covariates_distances[, chosen_form]
#   )
# 
#   boot_estimates <- list(
#     fit = gpd_fit$par,
#     intensity_fit = intensity_boot
#   )
# 
#   return(boot_estimates)
# }
# 
# # Generation of parmetric samples for combining Alg 1 and Alg 3
# 
# rgpd <- function(n, shape, scale = NULL, nu = NULL, mu = 0){
#   ## Input checks
#   # one and only one of {nu, scale} may be specified
#   if (is.null(scale) & is.null(nu)) {
#     stop('Define one of the parameters nu or scale.')
#   }
#   if (!is.null(scale) & !is.null(nu)) {
#     stop('Define only one of the parameters nu and scale.')
#   }
#   # Calculate scale from nu if required
#   if (!is.null(nu) & is.null(scale)) {
#     scale <- nu / (1 + shape)
#     if (any(scale <= 0)) {
#       stop('Implied scale parameter(s) must be positive.')
#     }
#   }
#   # Check that scale value(s) are positive
#   if (any(scale <= 0)) {
#     stop('Scale parameter(s) must be positive.')
#   }
#   # Ensure q, scale, shape and mu are of same length.
#   if ((length(scale) == 1) & (n > 1)) {
#     scale <- rep(scale, n)
#   }
#   if ((length(shape) == 1) & (n > 1)) {
#     shape <- rep(shape, n)
#   }
#   if ((length(mu) == 1) & (n > 1)) {
#     mu <- rep(mu, n)
#   }
# 
#   #simulate sample
#   sample <- mu + (scale/shape) * ((1 - stats::runif(n))^(-shape) - 1)
#   #correct sample values where xi = 0
#   #ex <- which(shape ==0)
#   if (any(abs(shape) < 1e-10)) {
#     ex <- which(abs(shape) < 1e-10)
#     sample[ex] <- mu[ex] +
#       stats::rexp(n = length(ex),rate = 1/scale[ex])
#   }
#   return(sample)
# }
# 
# 
# generate_parametric_bootstrap_samples <- function(probs, aggregated_intensity, thresh_fit, covariates, chosen_form, grid_box_area=0.2445286, n_samples=200, SEED=11111) {
# 
#   n_cells <- length(probs)
#   boot_samples <- vector("list", n_samples)
#   set.seed(SEED)
#   for (i in seq_len(n_samples)) {
#     # 1. Sample number of events
#     n_obs <- rpois(1, aggregated_intensity)
#     sampled_indices <- sample.int(n_cells, n_obs, replace = TRUE, prob = probs)
#     sampled_covariates <- covariates[sampled_indices, , drop = FALSE]
# 
#     # 2. Calculate distances and threshold values
#     V_form <- covariates_distances[sampled_indices, chosen_form]
#     thresh_vals <- thresh_fit$thresh_par[1] + thresh_fit$thresh_par[2] * V_form
#     sigma_vals <- thresh_fit$par[1] + thresh_fit$par[2]*sampled_covariates$ICS_max + thresh_fit$par[3]*thresh_vals
# 
#     # 3. Simulate GPD excesses
#     xi <- thresh_fit$par[3]
#     boot_excesses <- rgpd(n_obs, scale = sigma_vals, shape = xi)
# 
#     boot_samples[[i]] <- list(
#       sampled_covariates = sampled_covariates,
#       V_form = V_form,
#       thresh_vals = thresh_vals,
#       boot_excesses = boot_excesses
#     )
#   }
# 
#   return(boot_samples)
# }
# 
# covariates <- read.csv("STORM_input/covariates_1995-2024.csv", header=T)
# threshold <- 1.45
# 
# observed_sample <- read.csv("STORM_input/unrounded_after_1995_in_polygon_with_covariates.csv", header=T)
# 
# excess_data <- observed_sample[observed_sample$Magnitude > threshold,]
# 
# fit_obs <- optim(GPD_LL_given_V_ICS, par = c(0.1, 0, 0.1), excess = excess_data$Magnitude - threshold,
#                  thresh_par = c(1.45, 0), V = excess_data$V_1, ics = excess_data$ICS_max,
#                  control = list(fnscale = -1))
# 
# thresh_fit <- list(thresh_par = c(1.45, 0), par = fit_obs$par)
#                    
#                    
# covariates_distances <- cbind(covariates$V1, covariates$V2, covariates$V3, covariates$V4,
#                               log(covariates$V1), log(covariates$V2), log(covariates$V3), log(covariates$V4),
#                               sqrt(covariates$V1), sqrt(covariates$V2), sqrt(covariates$V3), sqrt(covariates$V4))
# grid_box_area <- 0.2445286  # Area of grid box in km^2
# 
# library(parallel)
# 
# chosen_form <- 1 # If using Alg 2, fix to chosen model (V2)
# 
# intensity_par <- intensity_estimation(observed_sample, covariates, thresh_fit, log(observed_sample$V_1), covariates_distances[,chosen_form])
# covariates_threshold <- thresh_fit$thresh_par[1] + thresh_fit$thresh_par[2] * covariates_distances[, chosen_form]
# 
# estimated_intensity <- resulting_intensity_icsmax(intensity_par, covariates, covariates_threshold, thresh_fit)
# aggregated_intensity <- sum(estimated_intensity, na.rm = TRUE)*grid_box_area
# probs <- estimated_intensity / aggregated_intensity
# 
# 
# # Generate parametric bootstrap samples
# bootstrap_samples_Alg1 <- generate_parametric_bootstrap_samples(
#   probs = probs,
#   aggregated_intensity = aggregated_intensity,
#   thresh_fit = thresh_fit,
#   covariates = covariates,
#   chosen_form = chosen_form,
#   n_samples = 200,
#   SEED = 11111
# )
# cl <- makeCluster(100)  # Match this with your SBATCH -c value minus 1 ideally
# 
# clusterExport(cl, varlist = c("parametric_bootstrap_fits","GPD_LL_given_V_ICS",
#                               "intensity_estimation", "Poisson_process_LL_icsmax","covariates",
#                               "covariates_distances", "thresh_fit",
#                               "chosen_form" ), envir = environment())
# bootstrap_fits_Alg1 <- parLapply(cl, bootstrap_samples_Alg1, function(sample) {
#   # Generate parametric bootstrap samples
#   parametric_bootstrap_fits(
#     boot_sample = sample,
#     thresh_fit = thresh_fit,
#     covariates_distances = covariates_distances,
#     chosen_form = chosen_form
#   )
# })
# stopCluster(cl)
# 
# 
# 
# saveRDS(bootstrap_fits_Alg1, file="STORM_output/bootstrap_fits_Alg1_conservative.rds")
# 
# 
