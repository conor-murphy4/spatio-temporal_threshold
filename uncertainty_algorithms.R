source("src/eqd_geo.R")

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

  