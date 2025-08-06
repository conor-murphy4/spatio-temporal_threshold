source("src/eqd_geo.R")
source("intensity_estimation.R")

generate_non_par_bootstrap_samples <- function(data, n_samples) {
  n <- nrow(data)
  samples <- vector("list", n_samples)
  
  for (i in seq_len(n_samples)) {
    indices <- sample(1:n, replace = TRUE)
    samples[[i]] <- data[indices, ]
  }
  
  return(samples)
}

observed_sample <- gron_eq_cat[,c("Date","Easting", "Northing", "Depth", "Magnitude", "V_1", "V_2" , "V_3", "V_4", "ICS_max", "dsmaxdt")]
# Generate non-parametric bootstrap samples
set.seed(11111)
non_parametric_samples <- generate_non_par_bootstrap_samples(observed_sample, 200)
#saveRDS(non_parametric_samples, file="threshold_results/uncertainty/non_parametric_samples.rds")

# Function to apply threshold selection 
# V-based thresholds
intercepts <- seq(0, 1.5, by=0.02)
slopes <- seq(0,0.8, by=0.02)
threshold_matrix <- as.matrix(expand.grid(intercepts, slopes))

threshold_selection_varying_formulation <- function(data, threshold_matrix, num_boot=200, SEED=11111, output_all = FALSE) {
  nearest_dist_matrix <- matrix(c(data$V_1, data$V_2, data$V_3, data$V_4), ncol=4, byrow=F)
    results <- vector("list", 12)
    for (ii in c(1:4)) {
      nearest_dist <- nearest_dist_matrix[,ii]
      log_nearest_dist <- log(nearest_dist)
      sqrt_nearest_dist <- sqrt(nearest_dist)
      
      # V
      set.seed(SEED)
      results[[ii]] <- eqd_geo_ics(data=data$Magnitude, thresh = threshold_matrix, distance_to_geo = nearest_dist, k=num_boot, ics = data$ICS_max)
      
      # log(V)
      set.seed(SEED)
      results[[ii+4]] <- eqd_geo_ics(data=data$Magnitude, thresh = threshold_matrix, distance_to_geo =  log_nearest_dist, k=num_boot, ics = data$ICS_max)
      
      # sqrt(V)
      set.seed(SEED)
      results[[ii+8]] <- eqd_geo_ics(data=data$Magnitude, thresh = threshold_matrix, distance_to_geo = sqrt_nearest_dist, k=num_boot, ics = data$ICS_max)
    }
    chosen_model_index <- which.min(lapply(results, function(x){min(x$dists)} ))
    chosen_model <- results[[chosen_model_index]]
    if (output_all) {
      return(list(chosen_form = chosen_model_index, model_results=results))
    }
    else {
      return(list(chosen_form = chosen_model_index, model_results=chosen_model))
    }
}

threshold_selection_fixed_formulation <- function(data, threshold_matrix, num_boot=200, SEED=11111) {
    # A2
    set.seed(SEED)
    results <- eqd_geo_ics(data=data$Magnitude, thresh = threshold_matrix, distance_to_geo = data$V_2, k=num_boot, ics = data$ICS_max)
    return(results)
}

bootstrap_results <- vector("list", length(non_parametric_samples))
for (i in seq_along(non_parametric_samples)) {
  print(i)
  bootstrap_results[[i]] <- threshold_selection_fixed_formulation(data=non_parametric_samples[[i]], threshold_matrix=threshold_matrix, num_boot=100, SEED=11111)
}

saveRDS(bootstrap_results, file="threshold_results/uncertainty/threshold_values_uncertainty_results.rds")


# TODO - Fitting function for GPD fits 
# TODO - function for endpoint estimates
# TODO - function for quantile estimates

# Alg 1

# Generate parametric bootstrap samples

dist_matrix <- cbind(gron_eq_cat$V_1, gron_eq_cat$V_2, gron_eq_cat$V_3, gron_eq_cat$V_4, log(gron_eq_cat$V_1), 
                        log(gron_eq_cat$V_2), log(gron_eq_cat$V_3), log(gron_eq_cat$V_4), 
                        sqrt(gron_eq_cat$V_1), sqrt(gron_eq_cat$V_2), sqrt(gron_eq_cat$V_3), sqrt(gron_eq_cat$V_4))
covariates_distances <- cbind(covariates$V1, covariates$V2, covariates$V3, covariates$V4, 
                              log(covariates$V1), log(covariates$V2), log(covariates$V3), log(covariates$V4), 
                              sqrt(covariates$V1), sqrt(covariates$V2), sqrt(covariates$V3), sqrt(covariates$V4))
distance_obs <- dist_matrix[,chosen_model_index] 
distance_covariate <- covariates_distances[,chosen_model_index]

intensity_estimation <- function(data, covariates, thresh_fit, distance_obs, distance_covariate, grid_box_area=0.2445286) {
  threshold_obs <- thresh_fit$thresh_par[1] + thresh_fit$thresh_par[2] * distance_obs
  covariates_threshold <- thresh_fit$thresh_par[1] + thresh_fit$thresh_par[2] * distance_covariate
  
  PP_fit <- optim(Poisson_process_LL_icsmax, par = c(0.1, 0), data = data, covariates = covariates, 
                  threshold_obs = threshold_obs, covariates_threshold = covariates_threshold, 
                  thresh_fit = thresh_fit, control=list(fnscale=-1))
  
  estimated_intensity <- resulting_intensity_icsmax(PP_fit, covariates, covariates_threshold, thresh_fit)
  aggregated_intensity <- sum(estimated_intensity, na.rm = TRUE)*grid_box_area
  return(list(fit=PP_fit, est_intensity = estimated_intensity, agg_intensity = aggregated_intensity))
}


intensity_fit_obs <- intensity_estimation(gron_eq_cat, covariates, thresh_fit_A2, dist_matrix[,2], covariates_distances[,2])


generate_parametric_bootstrap_estimates <- function(intensity_fit, thresh_fit, covariates, chosen_form, n_samples) {
  
  boot_estimates <- vector("list", length = n_samples)
  
  for (i in seq_len(n_samples)) {
    n_obs <- rpois(1, intensity_fit$agg_intensity)
    probs <- intensity_fit$est_intensity / intensity_fit$agg_intensity
    sampled_indices <- sample(1:length(intensity_fit$est_intensity), size = n_obs, replace = TRUE, prob = probs)
    sampled_covariates <- covariates[sampled_indices, ]
    sampled_covariates_distances <- cbind(sampled_covariates$V_1, sampled_covariates$V_2, sampled_covariates$V_3, sampled_covariates$V_4, 
                                          log(sampled_covariates$V_1), log(sampled_covariates$V_2), log(sampled_covariates$V_3), log(sampled_covariates$V_4), 
                                          sqrt(sampled_covariates$V_1), sqrt(sampled_covariates$V_2), sqrt(sampled_covariates$V_3), sqrt(sampled_covariates$V_4))
    
    threshold_vals <- thresh_fit$thresh_par[1] + thresh_fit$thresh_par[2] * sampled_covariates_distances[,chosen_form]
    sigma_vals <- thresh_fit$par[1] + thresh_fit$par[2]*sampled_covariates$ICS_max + thresh_fit$par[3]*threshold_vals
    boot_excesses <- numeric(n_obs)
    for(j in seq_len(n_obs)) {
      boot_excesses[j] <- rgpd(1, scale = sigma_vals[j], shape = thresh_fit$par[3])
    }
    ifelse(thresh_fit$par[3] < 0, pars_init <-  c(mean(boot_excesses), 0, 0.1) ,pars_init <- c(thresh_fit$par) )
    gpd_fit <- optim(GPD_LL_given_V_ICS, par = pars_init, 
                     excess = boot_excesses, thresh_par=thresh_fit$thresh_par, V=sampled_covariates_distances[, chosen_form], 
                     ics= sampled_covariates$ICS_max, control=list(fnscale=-1))
    sampled_covariates$Magnitude <- boot_excesses + threshold_vals
    intensity_fit <- intensity_estimation(sampled_covariates, covariates, thresh_fit, 
                                        distance_obs = sampled_covariates_distances[, chosen_form], 
                                        distance_covariate = covariates_distances[, chosen_form])
    boot_estimates[[i]] <- list(fit = gpd_fit$par, intensity_fit = intensity_fit$PP_fit$par)
    
  }
  
  return(boot_estimates)
}


generate_parametric_bootstrap_estimates <- function(intensity_fit, thresh_fit, covariates, chosen_form, n_samples, SEED=11111) {
  
  probs <- intensity_fit$est_intensity / intensity_fit$agg_intensity
  n_cells <- length(probs)
  
  covariates_distances <- cbind(covariates$V1, covariates$V2, covariates$V3, covariates$V4,
                                log(covariates$V1), log(covariates$V2), log(covariates$V3), log(covariates$V4), 
                                sqrt(covariates$V1), sqrt(covariates$V2), sqrt(covariates$V3), sqrt(covariates$V4))
    
  
  boot_samples <- vector("list", n_samples)
  set.seed(SEED)
  for (i in seq_len(n_samples)) {
    # 1. Sample number of events
    n_obs <- rpois(1, intensity_fit$agg_intensity)
    sampled_indices <- sample.int(n_cells, n_obs, replace = TRUE, prob = probs)
    sampled_covariates <- covariates[sampled_indices, , drop = FALSE]
    
    # 2. Calculate distances and threshold values
    V_form <- covariates_distances[sampled_indices, chosen_form]
    thresh_vals <- thresh_fit$thresh_par[1] + thresh_fit$thresh_par[2] * V_form
    sigma_vals <- thresh_fit$par[1] + thresh_fit$par[2]*sampled_covariates$ICS_max + thresh_fit$par[3]*thresh_vals
    
    # 3. Simulate GPD excesses
    xi <- thresh_fit$par[3]
    boot_excesses <- rgpd(n_obs, scale = sigma_vals, shape = xi)
    
    boot_samples[[i]] <- list(
      sampled_covariates = sampled_covariates,
      V_form = V_form,
      thresh_vals = thresh_vals,
      boot_excesses = boot_excesses
    )
  }
  
  return(boot_samples)
}

parametric_bootstrap_samples <- generate_parametric_bootstrap_estimates(intensity_fit = intensity_fit_obs, thresh_fit = thresh_fit_A2, 
                                        covariates = covariates, chosen_form = 2, n_samples = 200)

saveRDS(parametric_bootstrap_samples, file="threshold_results/uncertainty/parametric_samples_Alg1.rds")

parametric_bootstrap_fits <- function(boot_samples, thresh_fit, covariates_distances, chosen_form) {
  

    sampled_covariates <- boot_samples$sampled_covariates
    V_form <- boot_samples$V_form
    thresh_vals <- boot_samples$thresh_vals
    boot_excesses <- boot_samples$boot_excesses
    
    # Fit GPD model
    xi <- thresh_fit$par[3]
    init_mle <- mean(boot_excesses)
    pars_init <- if (xi < 0) c(init_mle, 0, 0.1) else thresh_fit$par
    gpd_fit <- optim(
      GPD_LL_given_V_ICS, 
      par = pars_init, 
      excess = boot_excesses, 
      thresh_par = thresh_fit$thresh_par, 
      V = V_form, 
      ics = sampled_covariates$ICS_max, 
      control = list(fnscale = -1)
    )
    
    # Add magnitudes and re-fit intensity model
    sampled_covariates$Magnitude <- boot_excesses + thresh_vals
    intensity_boot <- intensity_estimation(
      sampled_covariates, covariates, thresh_fit,
      distance_obs = V_form,
      distance_covariate = covariates_distances[, chosen_form]
    )
    
    boot_estimates <- list(
      fit = gpd_fit$par, 
      intensity_fit = intensity_boot$fit$par
    )
  
  return(boot_estimates)
}

parametric_bootstrap_fits(boot_samples = parametric_bootstrap_samples[[1]], 
                          thresh_fit = thresh_fit_A2, 
                          covariates_distances = covariates_distances, 
                          chosen_form = 2)

# 4. Fit GPD model
init_mle <- mean(boot_excesses)
pars_init <- if (xi < 0) c(init_mle, 0, 0.1) else thresh_fit$par
gpd_fit <- optim(
  GPD_LL_given_V_ICS, 
  par = pars_init, 
  excess = boot_excesses, 
  thresh_par = thresh_fit$thresh_par, 
  V = V_form, 
  ics = sampled_covariates$ICS_max, 
  control = list(fnscale = -1)
)

# 5. Add magnitudes and re-fit intensity model
sampled_covariates$Magnitude <- boot_excesses + thresh_vals
intensity_boot <- intensity_estimation(
  sampled_covariates, covariates, thresh_fit,
  distance_obs = V_form,
  distance_covariate = covariates_distances[, chosen_form]
)

boot_estimates[[i]] <- list(
  fit = gpd_fit$par, 
  intensity_fit = intensity_boot$PP_fit$par
)


# Standard errors from bootstrap estimates for Alg 1
bootstrap_estimates_Alg1 <- readRDS("uncertainty/bootstrap_fits_Alg1_C2.rds")

standard_errors_Alg1 <- function(boot_estimates) {
  n_boot <- length(boot_estimates)
  par_matrix <- do.call(rbind, lapply(boot_estimates, function(x) x$fit))
  intensity_matrix <- do.call(rbind, lapply(boot_estimates, function(x) x$intensity_fit))
  
  se_params <- apply(par_matrix, 2, sd)
  se_intensity <- apply(intensity_matrix, 2, sd)
  
  return(list(se_params = se_params, se_intensity = se_intensity))
}

(se_results_Alg1 <- standard_errors_Alg1(bootstrap_estimates_Alg1))

CIs_Alg1 <- function(boot_estimates, alpha = 0.05) {
  n_boot <- length(boot_estimates)
  par_matrix <- do.call(rbind, lapply(boot_estimates, function(x) x$fit))
  intensity_matrix <- do.call(rbind, lapply(boot_estimates, function(x) x$intensity_fit))
  
  lower_bound_params <- apply(par_matrix, 2, quantile, probs = alpha / 2)
  upper_bound_params <- apply(par_matrix, 2, quantile, probs = 1 - alpha / 2)
  
  lower_bound_intensity <- apply(intensity_matrix, 2, quantile, probs = alpha / 2)
  upper_bound_intensity <- apply(intensity_matrix, 2, quantile, probs = 1 - alpha / 2)
  
  return(list(
    params_CI = data.frame(lower = lower_bound_params, upper = upper_bound_params),
    intensity_CI = data.frame(lower = lower_bound_intensity, upper = upper_bound_intensity)
  ))
}
(CIs_results_Alg1 <- CIs_Alg1(bootstrap_estimates_Alg1))



# Conservative threshold with Alg 1

conservative_threshold <- 1.45
excesses_above_threshold <- gron_eq_cat$Magnitude[gron_eq_cat$Magnitude > conservative_threshold] - conservative_threshold
mean_excess <- mean(excesses_above_threshold, na.rm = TRUE)

gpd_fit_observed_sample <- optim(
  GPD_LL, 
  par = c(mean_excess, 0.1), 
  z = excesses_above_threshold, 
  control = list(fnscale = -1)
)

intensity_fit_observed_sample <- optim(
  Poisson_process_LL_const_thresh, 
  par = c(0.1, 0), 
  data = gron_eq_cat, 
  covariates = covariates,
  control = list(fnscale = -1)
)

grid_box_area <- 0.2445286
num_samples <- 200
set.seed(11111)
for(ii in 1:num_samples){
  estimated_intensity <- resulting_intensity_const_thresh(intensity_fit_observed_sample, covariates)
  aggregated_intensity <- sum(estimated_intensity, na.rm = TRUE) * grid_box_area
  num_obs <- rpois(1, aggregated_intensity)
  ###### ????????
}

resulting_intensity_const_thresh <- function(opt_PP, covariates){
  intensity_thresh <- covariates$dsmaxdt*exp(opt_PP$par[1] + opt_PP$par[2] * covariates$ICS)
  return(intensity_thresh)
}



# Endpoint estimates with CIs for magnitude 3.6
gron_eq_cat <- read.csv("Data/Events/unrounded_after_1995_in_polygon_with_covariates.csv", header=T)

bootstrap_fits_Alg1 <- readRDS("uncertainty/bootstrap_fits_Alg1.rds")
bootstrap_fits_Alg1_B1 <- readRDS("uncertainty/bootstrap_fits_Alg1_B1.rds")
bootstrap_fits_Alg1_C2 <- readRDS("uncertainty/bootstrap_fits_Alg1_C2.rds")

par_ests_mat_A2 <- do.call(rbind, lapply(bootstrap_fits_Alg1, function(x) x$fit))
par_ests_mat_B1 <- do.call(rbind, lapply(bootstrap_fits_Alg1_B1, function(x) x$fit))
par_ests_mat_C2 <- do.call(rbind, lapply(bootstrap_fits_Alg1_C2, function(x) x$fit))

gron_largest_eq <- gron_eq_cat[gron_eq_cat$Magnitude > 3.5,]

endpoints_A2_at_largest <- -(par_ests_mat_A2[,1] + par_ests_mat_A2[,2] * gron_largest_eq$ICS_max) / par_ests_mat_A2[,3]
endpoints_B1_at_largest <- -(par_ests_mat_B1[,1] + par_ests_mat_B1[,2] * gron_largest_eq$ICS_max) / par_ests_mat_B1[,3]
endpoints_C2_at_largest <- -(par_ests_mat_C2[,1] + par_ests_mat_C2[,2] * gron_largest_eq$ICS_max) / par_ests_mat_C2[,3]

(CI_endpoints_A2 <- quantile(endpoints_A2_at_largest, probs = c(0.025, 0.975)))
(CI_endpoints_B1 <- quantile(endpoints_B1_at_largest, probs = c(0.025, 0.975)))
(CI_endpoints_C2 <- quantile(endpoints_C2_at_largest, probs = c(0.025, 0.975)))
