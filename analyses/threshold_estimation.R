#Loading data and relevant files

library(dplyr)
source("src/eqd.R")
source("src/eqd_geo.R")
source("src/distance_to_nearest_geo.R")

gron_eq_cat <- read.csv("Data/Events/unrounded_after_1995_in_polygon_with_covariates.csv", header=T)
covariates <- read.csv("Data/covariates/covariates_1995-2024.csv", header=T)
geophones_deepest <- read.csv("Data/Geophones/Geophones_processed_03-07-2024_deepest_only.csv")

mags <- gron_eq_cat$Magnitude
num_boot <- 200


# Geophone-based thresholds
intercepts <- seq(0, 1.5, by=0.02)
slopes <- seq(0,0.8, by=0.02)
threshold_matrix <- as.matrix(expand.grid(intercepts, slopes))

nearest_dist_matrix <- matrix(c(gron_eq_cat$V_1, gron_eq_cat$V_2, gron_eq_cat$V_3, gron_eq_cat$V_4), ncol=4, byrow=F)

for (ii in c(1:4)) {
  nearest_dist <- nearest_dist_matrix[,ii]
  log_nearest_dist <- log(nearest_dist)
  sqrt_nearest_dist <- sqrt(nearest_dist)
  
  # V
  cat(paste0("Starting selection for V_", ii, " threshold"))
  set.seed(11111)
  geo_thresh_fit <- eqd_geo_ics(data=mags, thresh = threshold_matrix, distance_to_geo = nearest_dist, k=num_boot, ics = gron_eq_cat$ICS_max)
  saveRDS(geo_thresh_fit, paste0("threshold_results/geo_thresh_fit_V",ii,".rds"))
  
  # log(V)
  cat(paste0("Starting selection for log(V_", ii, ") threshold"))
  set.seed(11111)
  log_geo_thresh_fit <- eqd_geo_ics(data=mags, thresh = threshold_matrix, distance_to_geo =  log_nearest_dist, k=num_boot, ics = gron_eq_cat$ICS_max)
  saveRDS(log_geo_thresh_fit, paste0("threshold_results/geo_thresh_fit_logV",ii,".rds"))
  
  # sqrt(V)
  cat(paste0("Starting selection for sqrt(V_", ii, ") threshold"))
  set.seed(11111)
  sqrt_geo_thresh_fit <- eqd_geo_ics(data=mags, thresh = threshold_matrix, distance_to_geo = sqrt_nearest_dist, k=num_boot, ics = gron_eq_cat$ICS_max)
  saveRDS(sqrt_geo_thresh_fit, paste0("threshold_results/geo_thresh_fit_sqrtV",ii,".rds"))
}
  
