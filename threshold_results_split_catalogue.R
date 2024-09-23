source("src/helper_functions.R")

gron_eq_cat <- read.csv("Data/Events/unrounded_after_1995_in_polygon.csv", header=T)
mags <- gron_eq_cat$Magnitude
third_nearest_dist_3d <- gron_eq_cat$V_3d
log_third_nearest_dist_3d <- log(third_nearest_dist_3d)
sqrt_third_nearest_dist_3d <- sqrt(third_nearest_dist_3d)

#Before and after changepoint
u_h_length <- which(gron_eq_cat$Date == as.Date("2015-01-06"))[1]
gron_eq_cat_before <- gron_eq_cat[1:u_h_length,]
mags_before <- gron_eq_cat_before$Magnitude
third_nearest_dist_3d_before <- gron_eq_cat_before$V_3d
log_third_nearest_dist_3d_before <- log(third_nearest_dist_3d_before)
sqrt_third_nearest_dist_3d_before <- sqrt(third_nearest_dist_3d_before)

gron_eq_cat_after <- gron_eq_cat[(u_h_length+1):nrow(gron_eq_cat),]
mags_after <- gron_eq_cat_after$Magnitude
third_nearest_dist_3d_after <- gron_eq_cat_after$V_3d
log_third_nearest_dist_3d_after <- log(third_nearest_dist_3d_after)
sqrt_third_nearest_dist_3d_after <- sqrt(third_nearest_dist_3d_after)

#Below and above V=5
gron_eq_cat_belowV <- gron_eq_cat[gron_eq_cat$V_3d <= 5,]
mags_belowV <- gron_eq_cat_belowV$Magnitude
third_nearest_dist_3d_belowV <- gron_eq_cat_belowV$V_3d
log_third_nearest_dist_3d_belowV <- log(third_nearest_dist_3d_belowV)
sqrt_third_nearest_dist_3d_belowV <- sqrt(third_nearest_dist_3d_belowV)

gron_eq_cat_aboveV <- gron_eq_cat[gron_eq_cat$V_3d > 5,]
mags_aboveV <- gron_eq_cat_aboveV$Magnitude
third_nearest_dist_3d_aboveV <- gron_eq_cat_aboveV$V_3d
log_third_nearest_dist_3d_aboveV <- log(third_nearest_dist_3d_aboveV)
sqrt_third_nearest_dist_3d_aboveV <- sqrt(third_nearest_dist_3d_aboveV)

# Fitted thresholds -------------------------------------------------------

# EQD threshold
eqd_thresh_fit <- readRDS("threshold_results/const_thresh_fit.rds")
eqd_threshold <- rep(eqd_thresh_fit$thresh, length(mags))

eqd_thresh_before <- eqd_threshold[1:u_h_length]
eqd_thresh_after <- eqd_threshold[(u_h_length+1):length(mags)]

eqd_belowV <- eqd_threshold[gron_eq_cat$V_3d <= 5]
eqd_aboveV <- eqd_threshold[gron_eq_cat$V_3d > 5]

# Conservative constant threshold from Zak's work -------------------------
conservative_threshold <- rep(1.45, length(mags))

conserv_thresh_before <- conservative_threshold[1:u_h_length]
conserv_thresh_after <- conservative_threshold[(u_h_length+1):length(mags)]

conserv_belowV <- conservative_threshold[gron_eq_cat$V_3d <= 5]
conserv_aboveV <- conservative_threshold[gron_eq_cat$V_3d > 5]

# Piece-wise constant threshold according to index ------------------------

u_h_length <- which(gron_eq_cat$Date == as.Date("2015-01-06"))[1]
piecewise_const_thresh <- c(rep(1.15, u_h_length), rep(0.76, length(mags) - u_h_length))

pc_thresh_before <- piecewise_const_thresh[1:u_h_length]
pc_thresh_after <- piecewise_const_thresh[(u_h_length+1):length(mags)]

pc_belowV <- piecewise_const_thresh[gron_eq_cat$V_3d <= 5]
pc_aboveV <- piecewise_const_thresh[gron_eq_cat$V_3d > 5]

# V-based thresholds
geo_thresh_fit_3d_before <- readRDS("threshold_results/geo_thresh_fit_3d_before.rds")
log_geo_thresh_fit_3d_before <- readRDS("threshold_results/log_geo_thresh_fit_3d_before.rds")
sqrt_geo_thresh_fit_3d_before <- readRDS("threshold_results/sqrt_geo_thresh_fit_3d_before.rds")

geo_thresh_fit_3d_after <- readRDS("threshold_results/geo_thresh_fit_3d_after.rds")
log_geo_thresh_fit_3d_after <- readRDS("threshold_results/log_geo_thresh_fit_3d_after.rds")
sqrt_geo_thresh_fit_3d_after <- readRDS("threshold_results/sqrt_geo_thresh_fit_3d_after.rds")

geo_thresh_fit_3d_belowV <- readRDS("threshold_results/geo_thresh_fit_3d_belowV.rds")
log_geo_thresh_fit_3d_belowV <- readRDS("threshold_results/log_geo_thresh_fit_3d_belowV.rds")
sqrt_geo_thresh_fit_3d_belowV <- readRDS("threshold_results/sqrt_geo_thresh_fit_3d_belowV.rds")

geo_thresh_fit_3d_aboveV <- readRDS("threshold_results/geo_thresh_fit_3d_aboveV.rds")
log_geo_thresh_fit_3d_aboveV <- readRDS("threshold_results/log_geo_thresh_fit_3d_aboveV.rds")
sqrt_geo_thresh_fit_3d_aboveV <- readRDS("threshold_results/sqrt_geo_thresh_fit_3d_aboveV.rds")

#Confidence intervals for parameters

get_par_ests_step(mags_before, conserv_thresh_before)
get_par_ests_step(mags_before, eqd_thresh_before)
get_par_ests_step(mags_before, pc_thresh_before)

get_par_ests_geo(mags_before, geo_thresh_fit_3d_before, third_nearest_dist_3d_before)
get_par_ests_geo(mags_before, log_geo_thresh_fit_3d_before, log_third_nearest_dist_3d_before)
get_par_ests_geo(mags_before, sqrt_geo_thresh_fit_3d_before, sqrt_third_nearest_dist_3d_before)

get_par_ests_step(mags_after, conserv_thresh_after)
get_par_ests_step(mags_after, eqd_thresh_after)
get_par_ests_step(mags_after, pc_thresh_after)

get_par_ests_geo(mags_after, geo_thresh_fit_3d_after, third_nearest_dist_3d_after)
get_par_ests_geo(mags_after, log_geo_thresh_fit_3d_after, log_third_nearest_dist_3d_after)
get_par_ests_geo(mags_after, sqrt_geo_thresh_fit_3d_after, sqrt_third_nearest_dist_3d_after)

# QQ plots

get_qq_plot(mags_before, geo_thresh_fit_3d_before, third_nearest_dist_3d_before)

plot_threshold <- function(gron_eq_cat, thresh_fit, third_nearest_dist){
  threshold <- thresh_fit$thresh[1] + thresh_fit$thresh[2]*third_nearest_dist
  plot(third_nearest_dist, gron_eq_cat$Magnitude, xlab = "Third nearest distance", ylab = "Magnitude")
  points(third_nearest_dist, threshold, col="red", pch=19)

}

plot_threshold(gron_eq_cat_before, geo_thresh_fit_3d_before, third_nearest_dist_3d_before)



# After changepoint

get_par_ests(mags_after, geo_thresh_fit_3d_after, third_nearest_dist_3d_after)
get_par_ests(mags_after, log_geo_thresh_fit_3d_after, log_third_nearest_dist_3d_after)
get_par_ests(mags_after, sqrt_geo_thresh_fit_3d_after, sqrt_third_nearest_dist_3d_after)


get_par_ests_geo(mags_belowV, geo_thresh_fit_3d_belowV, third_nearest_dist_3d_belowV)
get_par_ests_geo(mags_belowV, log_geo_thresh_fit_3d_belowV, log_third_nearest_dist_3d_belowV)
get_par_ests_geo(mags_belowV, sqrt_geo_thresh_fit_3d_belowV, sqrt_third_nearest_dist_3d_belowV)

get_par_ests_geo(mags_aboveV, geo_thresh_fit_3d_aboveV, third_nearest_dist_3d_aboveV)
get_par_ests_geo(mags_aboveV, log_geo_thresh_fit_3d_aboveV, log_third_nearest_dist_3d_aboveV)
get_par_ests_geo(mags_aboveV, sqrt_geo_thresh_fit_3d_aboveV, sqrt_third_nearest_dist_3d_aboveV)

