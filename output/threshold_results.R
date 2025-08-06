library(ggplot2)
library(dplyr)
source("src/helper_functions.R")


gron_eq_cat <- read.csv("Data/Events/unrounded_after_1995_in_polygon_with_covariates.csv", header=T)
mags <- gron_eq_cat$Magnitude
nearest_dist_matrix <- matrix(c(gron_eq_cat$V_1, gron_eq_cat$V_2, gron_eq_cat$V_3, gron_eq_cat$V_4), byrow=F, ncol=4)
covariates <- read.csv("Data/covariates/covariates_1995-2024.csv", header=T)

# Threshold results

#ICS max (B=200) - Results shown in paper
thresh_fit_V1_ics <- readRDS("threshold_results/geo_thresh_fit_V1.rds")
thresh_fit_V2_ics <- readRDS("threshold_results/geo_thresh_fit_V2.rds")
thresh_fit_V3_ics <- readRDS("threshold_results/geo_thresh_fit_V3.rds")
thresh_fit_V4_ics <- readRDS("threshold_results/geo_thresh_fit_V4.rds")

thresh_fit_logV1_ics <- readRDS("threshold_results/geo_thresh_fit_logV1.rds")
thresh_fit_logV2_ics <- readRDS("threshold_results/geo_thresh_fit_logV2.rds")
thresh_fit_logV3_ics <- readRDS("threshold_results/geo_thresh_fit_logV3.rds")
thresh_fit_logV4_ics <- readRDS("threshold_results/geo_thresh_fit_logV4.rds")

thresh_fit_sqrtV1_ics <- readRDS("threshold_results/geo_thresh_fit_sqrtV1.rds")
thresh_fit_sqrtV2_ics <- readRDS("threshold_results/geo_thresh_fit_sqrtV2.rds")
thresh_fit_sqrtV3_ics <- readRDS("threshold_results/geo_thresh_fit_sqrtV3.rds")
thresh_fit_sqrtV4_ics <- readRDS("threshold_results/geo_thresh_fit_sqrtV4.rds")

# B=1000
# thresh_fit_V1_ics <- readRDS("threshold_results/boot1000/geo_thresh_fit_V1.rds")
# thresh_fit_V2_ics <- readRDS("threshold_results/boot1000/geo_thresh_fit_V2.rds")
# thresh_fit_V3_ics <- readRDS("threshold_results/boot1000/geo_thresh_fit_V3.rds")
# thresh_fit_V4_ics <- readRDS("threshold_results/boot1000/geo_thresh_fit_V4.rds")
# 
# thresh_fit_logV1_ics <- readRDS("threshold_results/boot1000/geo_thresh_fit_logV1.rds")
# thresh_fit_logV2_ics <- readRDS("threshold_results/boot1000/geo_thresh_fit_logV2.rds")
# thresh_fit_logV3_ics <- readRDS("threshold_results/boot1000/geo_thresh_fit_logV3.rds")
# thresh_fit_logV4_ics <- readRDS("threshold_results/boot1000/geo_thresh_fit_logV4.rds")
# 
# thresh_fit_sqrtV1_ics <- readRDS("threshold_results/boot1000/geo_thresh_fit_sqrtV1.rds")
# thresh_fit_sqrtV2_ics <- readRDS("threshold_results/boot1000/geo_thresh_fit_sqrtV2.rds")
# thresh_fit_sqrtV3_ics <- readRDS("threshold_results/boot1000/geo_thresh_fit_sqrtV3.rds")
# thresh_fit_sqrtV4_ics <- readRDS("threshold_results/boot1000/geo_thresh_fit_sqrtV4.rds")



# Table 1 -----------------------------------------------------------------


geo_thresh_fit <- vector("list", 4)
log_geo_thresh_fit <- vector("list", 4)
sqrt_geo_thresh_fit <- vector("list", 4)
for(i in 1:4){
  geo_thresh_fit[[i]] <- readRDS(paste0("threshold_results/geo_thresh_fit_V", i, ".rds"))
  log_geo_thresh_fit[[i]] <- readRDS(paste0("threshold_results/geo_thresh_fit_logV", i, ".rds"))
  sqrt_geo_thresh_fit[[i]] <- readRDS(paste0("threshold_results/geo_thresh_fit_sqrtV", i, ".rds"))
}

# EQD results for the geophone-based thresholds

# EQD values (min distances) from each threshold fit
eqd_vals <- numeric(12)
for(i in 1:4){
  eqd_vals[i] <- min(geo_thresh_fit[[i]]$dists, na.rm = TRUE)
  eqd_vals[i+4] <- min(log_geo_thresh_fit[[i]]$dists, na.rm = TRUE)
  eqd_vals[i+8] <- min(sqrt_geo_thresh_fit[[i]]$dists, na.rm = TRUE)
}
print(eqd_vals)


                               
