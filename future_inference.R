
source("src/distance_to_nearest_geo.R")
source("intensity_estimation.R")

library(dplyr)
library(pracma)
library(ggplot2)
library(ggspatial)
library(pracma)
library(cowplot)
library(dplyr)
library(purrr)
library(patchwork)

gron_eq_cat <- read.csv("Data/Events/unrounded_after_1995_in_polygon_with_covariates.csv", header=T)
geophones_deepest <- read.csv("Data/Geophones/Geophones_processed_03-07-2024_deepest_only.csv", header=T, row.names = 1)
covariates <- read.csv("Data/covariates/covariates_1995-2024.csv", header=T)
future_covariates <- read.csv("Data/covariates/future_covariates_2024-2055.csv", header=T)
gron_outline <- read.csv('Data/Geophones/Groningen_Field_outline.csv', header=T)
thresh_fit_A2 <- readRDS("threshold_results/icsmax/geo_thresh_fit_V2.rds")

# Future period of 30 years from 01-01-2025 to 01-01-2055
# future_covariates <- future_covariates %>%
#   filter(Date >= "2025-01-01" & Date < "2055-01-01")
# future_covariates$Year <- stringr::str_sub(future_covariates$Date, 1, 4) %>% as.numeric()

# future_covariates$V2 <- distance_to_nearest(future_covariates, geophones_deepest,2)
# future_covariates$V3 <- distance_to_nearest(future_covariates, geophones_deepest,3)
# Best threshold fit


# Poisson process fit on observed data
(opt_PP_LL_icsmax <- optim(Poisson_process_LL_icsmax, par = c(0.1, 0), data = gron_eq_cat, covariates = covariates,
                           threshold_obs = gron_eq_cat$threshold_A2, covariates_threshold = covariates$best_threshold ,
                           thresh_fit = thresh_fit_A2, control=list(fnscale=-1), hessian=T))


# future_covariates$intensity_above_0 <- resulting_intensity_icsmax(opt_PP_LL_icsmax, future_covariates, 0, thresh_fit_A2)
# future_covariates$sigma_0 <- thresh_fit_A2$par[1] + thresh_fit_A2$par[2] * future_covariates$ICS_max

# write.csv(future_covariates, "Data/covariates/future_covariates_2024-2055.csv", row.names = FALSE)

# These values will mostly be zero or close to zero due to the stresses not increasing

grid_box_E <- (max(unique(future_covariates$Easting))-min(unique(future_covariates$Easting)))/length(unique(future_covariates$Easting))
grid_box_N <- (max(unique(future_covariates$Northing))-min(unique(future_covariates$Northing)))/length(unique(future_covariates$Northing))
grid_box_area <- grid_box_E/1000* grid_box_N/1000

agg_intensity <- sum(future_covariates$intensity_above_0, na.rm = TRUE)

future_covariates <- future_covariates %>% group_by(Year) %>%
  mutate(agg_intensity_per_year = sum(intensity_above_0, na.rm = TRUE))%>%
  ungroup()

# future_covariates$normalised_intensity <- future_covariates$intensity_above_0 / future_covariates$agg_intensity_per_year    
# 
# future_covariates$endpoint <- -(thresh_fit_A2$par[1] + thresh_fit_A2$par[2] * future_covariates$ICS_max)/thresh_fit_A2$par[3]

endpoint_by_year <- future_covariates %>% group_by(Year) %>%
  summarise(endpoint_wm = sum(endpoint*normalised_intensity, na.rm = TRUE), agg_intensity_per_year = sum(intensity_above_0, na.rm = TRUE)) 

dev.new(height=5, width=10, noRStudioGD = TRUE)
par(mfrow=c(1,2), bg='transparent')
plot(endpoint_by_year$Year, endpoint_by_year$agg_intensity_per_year/agg_intensity, type = "l", 
     xlab = "Year", ylab = "Ratio of integrated intensities", lwd=2, ylim=c(0,0.08))
plot(endpoint_by_year$Year, endpoint_by_year$endpoint_wm, type = "l", 
     xlab = "Year", ylab = "Weighted mean endpoint", lwd=2)

  # Functions to evaluate endpoints
evaluate_endpoints <- function(future_covariates, intensity_par, gpd_par, grid_box_area = 0.2445286) {
  
  endpoints <- -(gpd_par[1] + gpd_par[2] * future_covariates$ICS_max)/gpd_par[3]
  
  endpoint_max <- max(endpoints, na.rm = TRUE)
  
  # weighted mean
  intensity_above0 <- future_covariates$dsmaxdt * exp(intensity_par[1] + intensity_par[2] * future_covariates$ICS_max) 
  
  agg_intensity <- sum(intensity_above0, na.rm = TRUE) 
  
  # Calculate the weighted mean endpoint
  endpoint_wm <- sum(endpoints * intensity_above0/agg_intensity , na.rm = TRUE)
  
  return(c(endpoint_max, endpoint_wm))
}


# Functions for evaluated v level  ----------------------------------------

# v-level extreme event estimation
v_level_extreme_event <- function(v, prob, future_covariates, intensity_par, gpd_par, grid_box_area = 0.2445286) {
  
  intensity_above0 <- future_covariates$dsmaxdt * exp(intensity_par[1] + intensity_par[2] * future_covariates$ICS_max)
  sigma0 <- gpd_par[1] + gpd_par[2] * future_covariates$ICS_max
  intensity_above_v <-  intensity_above0 * (1 + gpd_par[3] * (v) /sigma0 )^(-1 / gpd_par[3])
  
  intensity_above_v[1 + gpd_par[3] * (v) / sigma0 < 0] <- 0
  intensity_above_v[sigma0 <= 0] <- 0
  # Calculate aggregated intensity
  agg_intensity_v <- sum(intensity_above_v, na.rm = TRUE) * grid_box_area
  
  function_to_solve <- (exp(-agg_intensity_v) - prob)^2
  
  return(function_to_solve)
}

# v_level solver function
v_level_solver <- function(prob, future_covariates,  intensity_par, gpd_par, grid_box_area = 0.2445286) {
  
  if(gpd_par[3] >= 0) {
    upper_limit <- 13
  }
  else{
    # Calculate the endpoints
    endpoints <- -(gpd_par[1] + gpd_par[2] * future_covariates$ICS_max)/gpd_par[3]
    
    endpoint_max <- max(endpoints, na.rm = TRUE)
    if(endpoint_max > 13){
      upper_limit <- 13
    }
    else{
      upper_limit <- endpoint_max
    }
  }
  v_solution <- optimize(v_level_extreme_event, interval = c(0, upper_limit), prob = 0.9387,
                         future_covariates = future_covariates, intensity_par = intensity_par, gpd_par=gpd_par)$minimum
  return(v_solution)
}

v_level_solver(prob = 0.9387, future_covariates = future_covariates, 
               intensity_par=opt_PP_LL_icsmax$par, 
               gpd_par=thresh_fit_A2$par)
  


# Future inference with uncertainty ---------------------------------------


# Alg 1

bootstrap_estimates_Alg1 <- readRDS("uncertainty/bootstrap_fits_Alg1.rds")

v_levels_Alg1 <- unlist(lapply(bootstrap_estimates_Alg1, function(x) {v_level_solver(prob = 0.9387, 
                                                                  future_covariates = future_covariates, 
                                                                  intensity_par = x$intensity_fit, 
                                                                  gpd_par = x$fit)$minimum}))
total_sum <- 0

for(jj in 1:200){
  sum_current <- sum(bootstrap_estimates_Alg1[[jj]]$fit[3] >  0)
  if(bootstrap_estimates_Alg1[[jj]]$fit[3] >  0){
    print(bootstrap_estimates_Alg1[[jj]])
    stop()
  }
total_sum  <- total_sum + sum_current
}

quantile(sapply(bootstrap_estimates_Alg1, function(x){x$fit[3]}), c(0.025, 0.975))

sd(sapply(bootstrap_estimates_Alg1, function(x){x$fit[3]}))
(CI_v_levels <- quantile(v_levels_Alg1, c(0.025, 0.975)))


endpoints_Alg1 <- matrix(unlist(lapply(bootstrap_estimates_Alg1, function(x) {
  evaluate_endpoints(future_covariates, x$intensity_fit, x$fit)})), ncol = 2, byrow = TRUE)

(CI_endpoints_Alg1 <- apply(endpoints_Alg1, 2, quantile, c(0.025, 0.975)))

# Alg 2 and 3 done in parallel on storm

future_inferences_Alg2 <- readRDS("uncertainty/future_inferences_Alg2.rds")
future_inferences_Alg3 <- readRDS("uncertainty/future_inferences_Alg3.rds")

future_inferences_Alg2[is.na(future_inferences_Alg2[,2]),2] <- Inf
future_inferences_Alg3[is.na(future_inferences_Alg3[,2]),2] <- Inf

future_inferences_Alg2[is.na(future_inferences_Alg2[,3]),3] <- Inf
future_inferences_Alg3[is.na(future_inferences_Alg3[,3]),3] <- Inf

apply(future_inferences_Alg2, 2, quantile, c(0.025, 0.975)) 
apply(future_inferences_Alg3, 2, quantile, c(0.025, 0.975))

sum(is.na(future_inferences_Alg2[,2]))/40000

shape_vec <- numeric(40000)
for (ii in 1:200) {
  boot_ests <- bootstrap_fits_Alg1and3[[ii]]
  # Extract the shape parameter from the fit
  # Assuming the shape parameter is the third element in the fit vector
  shape_vec[(ii - 1) * 200 + 1:(200*ii)] <- do.call(rbind, lapply(boot_ests, function(x) x$fit[3]))
}


quantile(shape_vec, c(0.025, 0.975), na.rm = TRUE)
hist(shape_vec)
num_excess <- numeric(200)
for (ii in 1:200) {
  boot_ests <- bootstrap_model_selection_results_Alg3[[ii]]
  # Extract the shape parameter from the fit
  # Assuming the shape parameter is the third element in the fit vector
  num_excess[ii] <- boot_ests$model_results$num_excess
}

min(num_excess)
max(num_excess)
mean(num_excess)
sd(num_excess)

threshold_values_uncertainty_results_Alg2[[1]]
quantile(endpoint_sorted, c(0.25, 0.75))

hist(future_inferences_Alg2[,3])

endpoint_sorted <- sort(future_inferences_Alg2[,3], decreasing = TRUE)

hist(endpoint_sorted[(0.025*40000):40000])

range(shape_vec[which(future_inferences_Alg2[,2] > 10)])

shape_vec <- numeric(200)
for (ii in 1:200) {
  boot_ests <- bootstrap_fits_Alg1[[ii]]
  # Extract the shape parameter from the fit
  # Assuming the shape parameter is the third element in the fit vector
  shape_vec[ii] <- boot_ests$fit[3]
}
sd(shape_vec)
hist(shape_vec, breaks = 50, main = "Shape parameter distribution", xlab = "Shape parameter")
sum(shape_vec > 0)/200

quantile(future_inferences_Alg3[,2], c(0.25, 0.75), na.rm = TRUE)

# Proporition of model forms chosen
chosen_models <- do.call(rbind, lapply(bootstrap_model_selection_results_Alg3, function(x) x$chosen_form))

sum(chosen_models %in% c(1,2,3,4))/ length(chosen_models)
sum(chosen_models %in% c(5,6,7,8))/ length(chosen_models)
sum(chosen_models %in% c(9,10,11,12))/ length(chosen_models)

sum(chosen_models %in% c(1,5,9))/ length(chosen_models)
sum(chosen_models %in% c(2,6,10))/ length(chosen_models)
sum(chosen_models %in% c(3,7,11))/ length(chosen_models)
sum(chosen_models %in% c(4,8,12))/ length(chosen_models)

sum(chosen_models == 2)/ length(chosen_models)
sum(chosen_models == 5)/ length(chosen_models)
sum(chosen_models == 10)/ length(chosen_models)

# Conservative future inference
future_inferences_Alg1_conservative <- readRDS("uncertainty/future_inferences_Alg1_conservative.rds")
bootstrap_estimates_Alg1_conservative <- readRDS("uncertainty/bootstrap_fits_Alg1_conservative.rds")

future_inferences_Alg1_conservative[is.na(future_inferences_Alg1_conservative[,2]),2] <- Inf
future_inferences_Alg1_conservative[is.na(future_inferences_Alg1_conservative[,3]),3] <- Inf

apply(future_inferences_Alg1_conservative, 2, quantile, c(0.025, 0.975))

gpd_par_mat <- matrix(unlist(lapply(bootstrap_estimates_Alg1_conservative, function(x) x$fit)), ncol = 3, byrow = TRUE)
intensity_par_mat <- matrix(unlist(lapply(bootstrap_estimates_Alg1_conservative, function(x) x$intensity_fit)), ncol = 2, byrow = TRUE)

weird_indexes <- which(future_inferences_Alg1_conservative[,2] > 20)
future_inferences_Alg1_conservative[which(future_inferences_Alg1_conservative[,2] > 20),2 ]

hist(future_inferences_Alg1_conservative[which(future_inferences_Alg1_conservative[,1] < 20),1])

weird_indexes
gpd_par_mat[weird_indexes[8], ]
intensity_par_mat[weird_indexes[8], ]
future_inferences_Alg1_conservative[weird_indexes[8], ]

x <- seq(0, 20, length.out = 50)
test_v_levels <- numeric(length(x))
for (ii in 1:length(x)) {
  test_v_levels[ii] <- v_level_extreme_event(x[ii], prob = 0.9387, 
                                              future_covariates = future_covariates, 
                                              intensity_par = intensity_par_mat[weird_indexes[8], ], 
                                              gpd_par = gpd_par_mat[weird_indexes[8], ])
}

plot(x, test_v_levels, type = "l", xlab = "v-level", ylab = "Function value")


