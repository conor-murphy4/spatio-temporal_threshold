#Loading data and relevant files
#library(quantreg)
library(dplyr)
source("src/eqd.R")
source("src/eqd_geo.R")
source("src/distance_to_nearest_geo.R")

gron_eq_cat <- read.csv("Data/Events/unrounded_after_1995_in_polygon_with_covariates.csv", header=T)
covariates <- read.csv("Data/covariates/covariates_1995-2024.csv", header=T)
geophones_deepest <- read.csv("Data/Geophones/Geophones_processed_03-07-2024_deepest_only.csv")

mags <- gron_eq_cat$Magnitude
num_boot <- 200

################################################################################
#covariates$V2 <- distance_to_nearest(covariates,geophones_deepest, index=2)
#covariates$V3 <- distance_to_nearest(covariates,geophones_deepest, index=3)
#covariates$V4 <- distance_to_nearest(covariates,geophones_deepest, index=4)

#covariates$V2
#min(covariates$Date)
#sum(is.na(covariates$V4))
################################################################################
################################################################################

# EQD fit
thresholds_eqd <- quantile(mags, seq(0,0.95, by=0.01))
cat("Starting selection for EQD threshold")
set.seed(11111)
eqd_thresh_fit <- eqd_exp(mags, thresh = thresholds_eqd, k = 200)
saveRDS(eqd_thresh_fit, "threshold_results/eqd_thresh_fit.rds")

# Stepped threshold fit
thresh_before_change <- seq(0.5, 2, by=0.02)
thresh_after_change <- seq(0.4, 1.5, by=0.02)
thresholds_pc <- as.matrix(expand.grid(thresh_before_change, thresh_after_change))
cat("Starting selection for stepped threshold")
set.seed(11111)
pc_thresh_fit <- eqd_stepped(mags, thresholds = thresholds_pc, change_index = change_index, k = 200)
saveRDS(pc_thresh_fit, "threshold_results/pc_thresh_fit.rds")

# V-based thresholds
intercepts <- seq(0, 1.5, by=0.02)
slopes <- seq(0,0.8, by=0.02)
threshold_matrix <- as.matrix(expand.grid(intercepts, slopes))

nearest_dist_matrix <- matrix(c(gron_eq_cat$V_1, gron_eq_cat$V_2, gron_eq_cat$V_3, gron_eq_cat$V_4), ncol=4, byrow=F)
current_covariates <- covariates %>% filter(Year == max(covariates$Year))
max_dists <- c(max(current_covariates$V1), max(current_covariates$V2), 
               max(current_covariates$V3), max(current_covariates$V4, na.rm = T))
min_ics <- min(covariates$ICS_max)
for (ii in c(1:4)) {
  nearest_dist <- nearest_dist_matrix[,ii]
  log_nearest_dist <- log(nearest_dist)
  sqrt_nearest_dist <- sqrt(nearest_dist)
  
  # V
  cat(paste0("Starting selection for V_", ii, " threshold"))
  set.seed(11111)
  geo_thresh_fit <- eqd_geo_ics(data=mags, thresh = threshold_matrix, distance_to_geo = nearest_dist, k=num_boot, ics = gron_eq_cat$ICS_max, max_dist = max_dists[ii], min_ics = min_ics)
  saveRDS(geo_thresh_fit, paste0("threshold_results/icsmax/geo_thresh_fit_V",ii,".rds"))
  
  # log(V)
  cat(paste0("Starting selection for log(V_", ii, ") threshold"))
  set.seed(11111)
  log_geo_thresh_fit <- eqd_geo_ics(data=mags, thresh = threshold_matrix, distance_to_geo =  log_nearest_dist, k=num_boot, ics = gron_eq_cat$ICS_max, max_dist = log(max_dists[ii]), min_ics = min_ics)
  saveRDS(log_geo_thresh_fit, paste0("threshold_results/icsmax/geo_thresh_fit_logV",ii,".rds"))
  
  # sqrt(V)
  cat(paste0("Starting selection for sqrt(V_", ii, ") threshold"))
  set.seed(11111)
  sqrt_geo_thresh_fit <- eqd_geo_ics(data=mags, thresh = threshold_matrix, distance_to_geo = sqrt_nearest_dist, k=num_boot, ics = gron_eq_cat$ICS_max, max_dist = sqrt(max_dists[ii]), min_ics = min_ics)
  saveRDS(sqrt_geo_thresh_fit, paste0("threshold_results/icsmax/geo_thresh_fit_sqrtV",ii,".rds"))
}
  
# Constructing splits in data ---------------------------------------------

# Split before and after first event in 2015
change_index <- which(gron_eq_cat$Date == "2015-01-06")[1]
gron_eq_cat_before <- gron_eq_cat[1:change_index,]
gron_eq_cat_after <- gron_eq_cat[(change_index+1):nrow(gron_eq_cat),]
mags_before <- gron_eq_cat_before$Magnitude
mags_after <- gron_eq_cat_after$Magnitude
third_nearest_dist_2d_before <- gron_eq_cat_before$V
third_nearest_dist_2d_after <- gron_eq_cat_after$V

third_nearest_dist_before <- gron_eq_cat_before$V
log_third_nearest_dist_before <- log(third_nearest_dist_before)
sqrt_third_nearest_dist_before <- sqrt(third_nearest_dist_before)

third_nearest_dist_after <- gron_eq_cat_after$V
log_third_nearest_dist_after <- log(third_nearest_dist_after)
sqrt_third_nearest_dist_after <- sqrt(third_nearest_dist_after)

# Split below and above V=5
V_cutoff <- 5
gron_eq_cat_belowV_2d <- gron_eq_cat[gron_eq_cat$V <= V_cutoff,]
gron_eq_cat_aboveV_2d <- gron_eq_cat[gron_eq_cat$V > V_cutoff,]
mags_belowV_2d <- gron_eq_cat_belowV_2d$Magnitude
mags_aboveV_2d <- gron_eq_cat_aboveV_2d$Magnitude
third_nearest_dist_belowV_2d <- gron_eq_cat_belowV_2d$V
third_nearest_dist_aboveV_2d <- gron_eq_cat_aboveV_2d$V

gron_eq_cat_belowV <- gron_eq_cat[gron_eq_cat$V <= V_cutoff,]
gron_eq_cat_aboveV <- gron_eq_cat[gron_eq_cat$V > V_cutoff,]
mags_belowV <- gron_eq_cat_belowV$Magnitude
mags_aboveV <- gron_eq_cat_aboveV$Magnitude

third_nearest_dist_belowV <- gron_eq_cat_belowV$V
log_third_nearest_dist_belowV <- log(third_nearest_dist_belowV)
sqrt_third_nearest_dist_belowV <- sqrt(third_nearest_dist_belowV)

third_nearest_dist_aboveV <- gron_eq_cat_aboveV$V
log_third_nearest_dist_aboveV <- log(third_nearest_dist_aboveV)
sqrt_third_nearest_dist_aboveV <- sqrt(third_nearest_dist_aboveV)


# Estimating thresholds for whole catalogue -------------------------------




# V_2D
cat("Starting selection for V_2D threshold")
set.seed(11111)
geo_thresh_fit_2d <- eqd_geo(data = mags, thresh = threshold_matrix, third_nearest_distance = third_nearest_dist_2d, k=500, min_dist = min(third_nearest_dist_2d), max_dist = max(third_nearest_dist_2d))
saveRDS(geo_thresh_fit_2d, "threshold_results/geo_thresh_fit_2d_k500.rds")


# Estimating thresholds before and after change in 2015 -------------------

# Before
cat("Starting selection for V_2D threshold before")
set.seed(11111)
geo_thresh_fit_2d_before <- eqd_geo(data=mags_before, thresh = threshold_matrix, third_nearest_distance = third_nearest_dist_2d_before, k=200, min_dist = min(third_nearest_dist_2d_before), max_dist = max(third_nearest_dist_2d_before))
saveRDS(geo_thresh_fit_2d_before, "threshold_results/geo_thresh_fit_2d_before.rds")

cat("Starting selection for V threshold before")
set.seed(11111)
geo_thresh_fit_before <- eqd_geo(data=mags_before, thresh = threshold_matrix, third_nearest_distance = third_nearest_dist_before, k=200, min_dist = min(third_nearest_dist_before), max_dist = max(third_nearest_dist_before))
saveRDS(geo_thresh_fit_before, "threshold_results/geo_thresh_fit_before.rds")

cat("Starting selection for log(V) threshold before")
set.seed(11111)
log_geo_thresh_fit_before <- eqd_geo(data=mags_before, thresh = threshold_matrix, third_nearest_distance = log_third_nearest_dist_before, k=200, min_dist = min(log_third_nearest_dist_before), max_dist = max(log_third_nearest_dist_before))
saveRDS(log_geo_thresh_fit_before, "threshold_results/log_geo_thresh_fit_before.rds")

cat("Starting selection for sqrt(V) threshold before")
set.seed(11111)
sqrt_geo_thresh_fit_before <- eqd_geo(data=mags_before, thresh = threshold_matrix, third_nearest_distance = sqrt_third_nearest_dist_before, k=200, min_dist = min(sqrt_third_nearest_dist_before), max_dist = max(sqrt_third_nearest_dist_before))
saveRDS(sqrt_geo_thresh_fit_before, "threshold_results/sqrt_geo_thresh_fit_before.rds")

# After
cat("Starting selection for V_2D threshold after")
set.seed(11111)
geo_thresh_fit_2d_after <- eqd_geo(data=mags_after, thresh = threshold_matrix, third_nearest_distance = third_nearest_dist_2d_after, k=200, min_dist = min(third_nearest_dist_2d_after), max_dist = max(third_nearest_dist_2d_after))
saveRDS(geo_thresh_fit_2d_after, "threshold_results/geo_thresh_fit_2d_after.rds")

cat("Starting selection for V threshold after")
set.seed(11111)
geo_thresh_fit_after <- eqd_geo(data=mags_after, thresh = threshold_matrix, third_nearest_distance = third_nearest_dist_after, k=200, min_dist = min(third_nearest_dist_after), max_dist = max(third_nearest_dist_after))
saveRDS(geo_thresh_fit_after, "threshold_results/geo_thresh_fit_after.rds")

cat("Starting selection for log(V) threshold after")
set.seed(11111)
log_geo_thresh_fit_after <- eqd_geo(data=mags_after, thresh = threshold_matrix, third_nearest_distance = log_third_nearest_dist_after, k=200, min_dist = min(log_third_nearest_dist_after), max_dist = max(log_third_nearest_dist_after))
saveRDS(log_geo_thresh_fit_after, "threshold_results/log_geo_thresh_fit_after.rds")

cat("Starting selection for log(V) threshold after")
set.seed(11111)
sqrt_geo_thresh_fit_after <- eqd_geo(data=mags_after, thresh = threshold_matrix, third_nearest_distance = sqrt_third_nearest_dist_after, k=200, min_dist = min(sqrt_third_nearest_dist_after), max_dist = max(sqrt_third_nearest_dist_after))
saveRDS(sqrt_geo_thresh_fit_after, "threshold_results/sqrt_geo_thresh_fit_after.rds")

# Estimating thresholds below and above V=5 -------------------------------

# Below
cat("Starting selection for V_2D threshold below")
set.seed(11111)
geo_thresh_fit_belowV_2d <- eqd_geo(data=mags_belowV_2d, thresh = threshold_matrix, third_nearest_distance = third_nearest_dist_belowV_2d, k=200, min_dist = min(third_nearest_dist_belowV_2d), max_dist = max(third_nearest_dist_belowV_2d))
saveRDS(geo_thresh_fit_belowV_2d, "threshold_results/geo_thresh_fit_belowV_2d.rds")

cat("Starting selection for V threshold below")
set.seed(11111)
geo_thresh_fit_belowV <- eqd_geo(data=mags_belowV, thresh = threshold_matrix, third_nearest_distance = third_nearest_dist_belowV, k=200, min_dist = min(third_nearest_dist_belowV), max_dist = max(third_nearest_dist_belowV))
saveRDS(geo_thresh_fit_belowV, "threshold_results/geo_thresh_fit_belowV.rds")

cat("Starting selection for log(V) threshold below")
set.seed(11111)
log_geo_thresh_fit_belowV <- eqd_geo(data=mags_belowV, thresh = threshold_matrix, third_nearest_distance = log_third_nearest_dist_belowV, k=200, min_dist = min(log_third_nearest_dist_belowV), max_dist = max(log_third_nearest_dist_belowV))
saveRDS(log_geo_thresh_fit_belowV, "threshold_results/log_geo_thresh_fit_belowV.rds")

cat("Starting selection for sqrt(V) threshold below")
set.seed(11111)
sqrt_geo_thresh_fit_belowV <- eqd_geo(data=mags_belowV, thresh = threshold_matrix, third_nearest_distance = sqrt_third_nearest_dist_belowV, k=200, min_dist = min(sqrt_third_nearest_dist_belowV), max_dist = max(sqrt_third_nearest_dist_belowV))
saveRDS(sqrt_geo_thresh_fit_belowV, "threshold_results/sqrt_geo_thresh_fit_belowV.rds")

# Above
cat("Starting selection for V_2D threshold above")
set.seed(11111)
geo_thresh_fit_aboveV_2d <- eqd_geo(data=mags_aboveV_2d, thresh = threshold_matrix, third_nearest_distance = third_nearest_dist_aboveV_2d, k=200, min_dist = min(third_nearest_dist_aboveV_2d), max_dist = max(third_nearest_dist_aboveV_2d))
saveRDS(geo_thresh_fit_aboveV_2d, "threshold_results/geo_thresh_fit_aboveV_2d.rds")

cat("Starting selection for V threshold above")
set.seed(11111)
geo_thresh_fit_aboveV <- eqd_geo(data=mags_aboveV, thresh = threshold_matrix, third_nearest_distance = third_nearest_dist_aboveV, k=200, min_dist = min(third_nearest_dist_aboveV), max_dist = max(third_nearest_dist_aboveV))
saveRDS(geo_thresh_fit_aboveV, "threshold_results/geo_thresh_fit_aboveV.rds")

cat("Starting selection for log(V) threshold above")
set.seed(11111)
log_geo_thresh_fit_aboveV <- eqd_geo(data=mags_aboveV, thresh = threshold_matrix, third_nearest_distance = log_third_nearest_dist_aboveV, k=200, min_dist = min(log_third_nearest_dist_aboveV), max_dist = max(log_third_nearest_dist_aboveV))
saveRDS(log_geo_thresh_fit_aboveV, "threshold_results/log_geo_thresh_fit_aboveV.rds")

cat("Starting selection for sqrt(V) threshold above")
set.seed(11111)
sqrt_geo_thresh_fit_aboveV <- eqd_geo(data=mags_aboveV, thresh = threshold_matrix, third_nearest_distance = sqrt_third_nearest_dist_aboveV, k=200, min_dist = min(sqrt_third_nearest_dist_aboveV), max_dist = max(sqrt_third_nearest_dist_aboveV))
saveRDS(sqrt_geo_thresh_fit_aboveV, "threshold_results/sqrt_geo_thresh_fit_aboveV.rds")