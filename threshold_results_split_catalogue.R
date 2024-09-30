source("src/helper_functions.R")

gron_eq_cat <- read.csv("Data/Events/unrounded_after_1995_in_polygon.csv", header=T)
mags <- gron_eq_cat$Magnitude
third_nearest_dist_2d <- gron_eq_cat$V
third_nearest_dist_3d <- gron_eq_cat$V_3d
log_third_nearest_dist_3d <- log(third_nearest_dist_3d)
sqrt_third_nearest_dist_3d <- sqrt(third_nearest_dist_3d)

#Before and after changepoint
u_h_length <- which(gron_eq_cat$Date == as.Date("2015-01-06"))[1]
gron_eq_cat_before <- gron_eq_cat[1:u_h_length,]
mags_before <- gron_eq_cat_before$Magnitude
third_nearest_dist_2d_before <- gron_eq_cat_before$V
third_nearest_dist_3d_before <- gron_eq_cat_before$V_3d
log_third_nearest_dist_3d_before <- log(third_nearest_dist_3d_before)
sqrt_third_nearest_dist_3d_before <- sqrt(third_nearest_dist_3d_before)

gron_eq_cat_after <- gron_eq_cat[(u_h_length+1):nrow(gron_eq_cat),]
mags_after <- gron_eq_cat_after$Magnitude
third_nearest_dist_2d_after <- gron_eq_cat_after$V
third_nearest_dist_3d_after <- gron_eq_cat_after$V_3d
log_third_nearest_dist_3d_after <- log(third_nearest_dist_3d_after)
sqrt_third_nearest_dist_3d_after <- sqrt(third_nearest_dist_3d_after)

#Below and above V=5
V_cutoff <- 5
gron_eq_cat_belowV <- gron_eq_cat[gron_eq_cat$V_3d <= V_cutoff,]
mags_belowV <- gron_eq_cat_belowV$Magnitude
third_nearest_dist_3d_belowV <- gron_eq_cat_belowV$V_3d
log_third_nearest_dist_3d_belowV <- log(third_nearest_dist_3d_belowV)
sqrt_third_nearest_dist_3d_belowV <- sqrt(third_nearest_dist_3d_belowV)

gron_eq_cat_belowV_2d <- gron_eq_cat[gron_eq_cat$V <= V_cutoff,]
mags_belowV_2d <- gron_eq_cat_belowV_2d$Magnitude
third_nearest_dist_2d_belowV <- gron_eq_cat_belowV_2d$V

gron_eq_cat_aboveV <- gron_eq_cat[gron_eq_cat$V_3d > V_cutoff,]
mags_aboveV <- gron_eq_cat_aboveV$Magnitude
third_nearest_dist_2d_aboveV <- gron_eq_cat_aboveV$V
third_nearest_dist_3d_aboveV <- gron_eq_cat_aboveV$V_3d
log_third_nearest_dist_3d_aboveV <- log(third_nearest_dist_3d_aboveV)
sqrt_third_nearest_dist_3d_aboveV <- sqrt(third_nearest_dist_3d_aboveV)

gron_eq_cat_aboveV_2d <- gron_eq_cat[gron_eq_cat$V > V_cutoff,]
mags_aboveV_2d <- gron_eq_cat_aboveV_2d$Magnitude
third_nearest_dist_2d_aboveV <- gron_eq_cat_aboveV_2d$V
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

# V-based thresholds ------------------------------------------------------
geo_thresh_fit_2d_before <- readRDS("threshold_results/geo_thresh_fit_2d_before.rds")
geo_thresh_fit_3d_before <- readRDS("threshold_results/geo_thresh_fit_3d_before.rds")
log_geo_thresh_fit_3d_before <- readRDS("threshold_results/log_geo_thresh_fit_3d_before.rds")
sqrt_geo_thresh_fit_3d_before <- readRDS("threshold_results/sqrt_geo_thresh_fit_3d_before.rds")

geo_thresh_fit_2d_after <- readRDS("threshold_results/geo_thresh_fit_2d_after.rds")
geo_thresh_fit_3d_after <- readRDS("threshold_results/geo_thresh_fit_3d_after.rds")
log_geo_thresh_fit_3d_after <- readRDS("threshold_results/log_geo_thresh_fit_3d_after.rds")
sqrt_geo_thresh_fit_3d_after <- readRDS("threshold_results/sqrt_geo_thresh_fit_3d_after.rds")

geo_thresh_fit_2d_belowV <- readRDS("threshold_results/geo_thresh_fit_2d_belowV.rds")
geo_thresh_fit_3d_belowV <- readRDS("threshold_results/geo_thresh_fit_3d_belowV.rds")
log_geo_thresh_fit_3d_belowV <- readRDS("threshold_results/log_geo_thresh_fit_3d_belowV.rds")
sqrt_geo_thresh_fit_3d_belowV <- readRDS("threshold_results/sqrt_geo_thresh_fit_3d_belowV.rds")

geo_thresh_fit_2d_aboveV <- readRDS("threshold_results/geo_thresh_fit_2d_aboveV.rds")
geo_thresh_fit_3d_aboveV <- readRDS("threshold_results/geo_thresh_fit_3d_aboveV.rds")
log_geo_thresh_fit_3d_aboveV <- readRDS("threshold_results/log_geo_thresh_fit_3d_aboveV.rds")
sqrt_geo_thresh_fit_3d_aboveV <- readRDS("threshold_results/sqrt_geo_thresh_fit_3d_aboveV.rds")

#Confidence intervals for parameters

#Before and after 06-01-2015
get_par_ests_step(mags_before, conserv_thresh_before)
get_par_ests_step(mags_before, eqd_thresh_before)
get_par_ests_step(mags_before, pc_thresh_before)

get_par_ests_geo(mags_before, geo_thresh_fit_2d_before, third_nearest_dist_2d_before)
get_par_ests_geo(mags_before, geo_thresh_fit_3d_before, third_nearest_dist_3d_before)
get_par_ests_geo(mags_before, log_geo_thresh_fit_3d_before, log_third_nearest_dist_3d_before)
get_par_ests_geo(mags_before, sqrt_geo_thresh_fit_3d_before, sqrt_third_nearest_dist_3d_before)

get_par_ests_step(mags_after, conserv_thresh_after)
get_par_ests_step(mags_after, eqd_thresh_after)
get_par_ests_step(mags_after, pc_thresh_after)

get_par_ests_geo(mags_after, geo_thresh_fit_2d_after, third_nearest_dist_2d_after)
get_par_ests_geo(mags_after, geo_thresh_fit_3d_after, third_nearest_dist_3d_after)
get_par_ests_geo(mags_after, log_geo_thresh_fit_3d_after, log_third_nearest_dist_3d_after)
get_par_ests_geo(mags_after, sqrt_geo_thresh_fit_3d_after, sqrt_third_nearest_dist_3d_after)

#Below and above V=5
get_par_ests_step(mags_belowV, conserv_belowV)
get_par_ests_step(mags_belowV, eqd_belowV)
get_par_ests_step(mags_belowV, pc_belowV)

get_par_ests_geo(mags_belowV_2d, geo_thresh_fit_2d_belowV, third_nearest_dist_2d_belowV)

get_par_ests_geo(mags_belowV, geo_thresh_fit_3d_belowV, third_nearest_dist_3d_belowV)
get_par_ests_geo(mags_belowV, log_geo_thresh_fit_3d_belowV, log_third_nearest_dist_3d_belowV)
get_par_ests_geo(mags_belowV, sqrt_geo_thresh_fit_3d_belowV, sqrt_third_nearest_dist_3d_belowV)

get_par_ests_step(mags_aboveV, conserv_aboveV)
get_par_ests_step(mags_aboveV, eqd_aboveV)
get_par_ests_step(mags_aboveV, pc_aboveV)

get_par_ests_geo(mags_aboveV_2d, geo_thresh_fit_2d_aboveV, third_nearest_dist_2d_aboveV)

get_par_ests_geo(mags_aboveV, geo_thresh_fit_3d_aboveV, third_nearest_dist_3d_aboveV)
get_par_ests_geo(mags_aboveV, log_geo_thresh_fit_3d_aboveV, log_third_nearest_dist_3d_aboveV)
get_par_ests_geo(mags_aboveV, sqrt_geo_thresh_fit_3d_aboveV, sqrt_third_nearest_dist_3d_aboveV)


# QQ plots

get_qq_plot(mags_before, geo_thresh_fit_3d_before, third_nearest_dist_3d_before)

# Visualising thresholds
plot_threshold <- function(gron_eq_cat, thresh_fit, third_nearest_dist, xlab=V){
  threshold <- thresh_fit$thresh[1] + thresh_fit$thresh[2]*third_nearest_dist
  plot(third_nearest_dist, gron_eq_cat$Magnitude, ylab = "Magnitude", xlab=xlab)
  points(third_nearest_dist, threshold, col="red", pch=19)

}

# Against V, logV and sqrt(V)
dev.new(height=12, width=30, noRStudioGD = TRUE)
par(mfrow=c(1,3), bg="transparent")
plot_threshold(gron_eq_cat, geo_thresh_fit_3d, third_nearest_dist_3d, xlab = "V_3d")
points(third_nearest_dist_3d_before, geo_thresh_fit_3d_before$thresh[1] + geo_thresh_fit_3d_before$thresh[2]*third_nearest_dist_3d_before, col="blue", pch=19)
points(third_nearest_dist_3d_after, geo_thresh_fit_3d_after$thresh[1] + geo_thresh_fit_3d_after$thresh[2]*third_nearest_dist_3d_after, col="green", pch=19)
points(third_nearest_dist_3d_belowV, geo_thresh_fit_3d_belowV$thresh[1] + geo_thresh_fit_3d_belowV$thresh[2]*third_nearest_dist_3d_belowV, col="purple", pch=19)
points(third_nearest_dist_3d_aboveV, geo_thresh_fit_3d_aboveV$thresh[1] + geo_thresh_fit_3d_aboveV$thresh[2]*third_nearest_dist_3d_aboveV, col="orange", pch=19)
legend("bottomright", legend=c("All", "Before", "After", "Below V=5", "Above V=5"), col=c("red", "blue", "green", "purple", "orange"), pch=19)

plot_threshold(gron_eq_cat, log_geo_thresh_fit_3d, log_third_nearest_dist_3d, xlab = "log(V_3d)")
points(log_third_nearest_dist_3d_before, log_geo_thresh_fit_3d_before$thresh[1] + log_geo_thresh_fit_3d_before$thresh[2]*log_third_nearest_dist_3d_before, col="blue", pch=19)
points(log_third_nearest_dist_3d_after, log_geo_thresh_fit_3d_after$thresh[1] + log_geo_thresh_fit_3d_after$thresh[2]*log_third_nearest_dist_3d_after, col="green", pch=19)
points(log_third_nearest_dist_3d_belowV, log_geo_thresh_fit_3d_belowV$thresh[1] + log_geo_thresh_fit_3d_belowV$thresh[2]*log_third_nearest_dist_3d_belowV, col="purple", pch=19)
points(log_third_nearest_dist_3d_aboveV, log_geo_thresh_fit_3d_aboveV$thresh[1] + log_geo_thresh_fit_3d_aboveV$thresh[2]*log_third_nearest_dist_3d_aboveV, col="orange", pch=19)
legend("bottomright", legend=c("All", "Before", "After", "Below V=5", "Above V=5"), col=c("red", "blue", "green", "purple", "orange"), pch=19)

plot_threshold(gron_eq_cat, sqrt_geo_thresh_fit_3d, sqrt_third_nearest_dist_3d, xlab = "sqrt(V_3d)")
points(sqrt_third_nearest_dist_3d_before, sqrt_geo_thresh_fit_3d_before$thresh[1] + sqrt_geo_thresh_fit_3d_before$thresh[2]*sqrt_third_nearest_dist_3d_before, col="blue", pch=19)
points(sqrt_third_nearest_dist_3d_after, sqrt_geo_thresh_fit_3d_after$thresh[1] + sqrt_geo_thresh_fit_3d_after$thresh[2]*sqrt_third_nearest_dist_3d_after, col="green", pch=19)
points(sqrt_third_nearest_dist_3d_belowV, sqrt_geo_thresh_fit_3d_belowV$thresh[1] + sqrt_geo_thresh_fit_3d_belowV$thresh[2]*sqrt_third_nearest_dist_3d_belowV, col="purple", pch=19)
points(sqrt_third_nearest_dist_3d_aboveV, sqrt_geo_thresh_fit_3d_aboveV$thresh[1] + sqrt_geo_thresh_fit_3d_aboveV$thresh[2]*sqrt_third_nearest_dist_3d_aboveV, col="orange", pch=19)
legend("bottomright", legend=c("All", "Before", "After", "Below V=5", "Above V=5"), col=c("red", "blue", "green", "purple", "orange"), pch=19)

#Against index
dev.new(height=12, width=30, noRStudioGD = TRUE)
par(mfrow=c(1,3), bg="transparent")

plot(gron_eq_cat$Magnitude, ylab = "Magnitude", xlab="Index", main="V_3d")
points(geo_thresh_fit_3d$thresh[1] + geo_thresh_fit_3d$thresh[2]*third_nearest_dist_3d, col="red", pch=19)
points(c(1:u_h_length),geo_thresh_fit_3d_before$thresh[1] + geo_thresh_fit_3d_before$thresh[2]*third_nearest_dist_3d_before, col="blue", pch=19)
points(c((u_h_length+1):nrow(gron_eq_cat)),geo_thresh_fit_3d_after$thresh[1] + geo_thresh_fit_3d_after$thresh[2]*third_nearest_dist_3d_after, col="green", pch=19)
points(which(gron_eq_cat$V_3d <= V_cutoff), geo_thresh_fit_3d_belowV$thresh[1] + geo_thresh_fit_3d_belowV$thresh[2]*third_nearest_dist_3d_belowV, col="purple", pch=19)
points(which(gron_eq_cat$V_3d > V_cutoff), geo_thresh_fit_3d_aboveV$thresh[1] + geo_thresh_fit_3d_aboveV$thresh[2]*third_nearest_dist_3d_aboveV, col="orange", pch=19)
legend("bottomright", legend=c("All", "Before", "After", "Below V=5", "Above V=5"), col=c("red", "blue", "green", "purple", "orange"), pch=19)

plot(gron_eq_cat$Magnitude, ylab = "Magnitude", xlab="Index", main="log(V_3d)")
points(log_geo_thresh_fit_3d$thresh[1] + log_geo_thresh_fit_3d$thresh[2]*log_third_nearest_dist_3d, col="red", pch=19)
points(c(1:u_h_length),log_geo_thresh_fit_3d_before$thresh[1] + log_geo_thresh_fit_3d_before$thresh[2]*log_third_nearest_dist_3d_before, col="blue", pch=19)
points(c((u_h_length+1):nrow(gron_eq_cat)),log_geo_thresh_fit_3d_after$thresh[1] + log_geo_thresh_fit_3d_after$thresh[2]*log_third_nearest_dist_3d_after, col="green", pch=19)
points(which(gron_eq_cat$V_3d <= V_cutoff), log_geo_thresh_fit_3d_belowV$thresh[1] + log_geo_thresh_fit_3d_belowV$thresh[2]*log_third_nearest_dist_3d_belowV, col="purple", pch=19)
points(which(gron_eq_cat$V_3d > V_cutoff), log_geo_thresh_fit_3d_aboveV$thresh[1] + log_geo_thresh_fit_3d_aboveV$thresh[2]*log_third_nearest_dist_3d_aboveV, col="orange", pch=19)
legend("bottomright", legend=c("All", "Before", "After", "Below V=5", "Above V=5"), col=c("red", "blue", "green", "purple", "orange"), pch=19)

plot(gron_eq_cat$Magnitude, ylab = "Magnitude", xlab="Index", main="sqrt(V_3d)")
points(sqrt_geo_thresh_fit_3d$thresh[1] + sqrt_geo_thresh_fit_3d$thresh[2]*sqrt_third_nearest_dist_3d, col="red", pch=19)
points(c(1:u_h_length),sqrt_geo_thresh_fit_3d_before$thresh[1] + sqrt_geo_thresh_fit_3d_before$thresh[2]*sqrt_third_nearest_dist_3d_before, col="blue", pch=19)
points(c((u_h_length+1):nrow(gron_eq_cat)), sqrt_geo_thresh_fit_3d_after$thresh[1] + sqrt_geo_thresh_fit_3d_after$thresh[2]*sqrt_third_nearest_dist_3d_after, col="green", pch=19)
points(which(gron_eq_cat$V_3d <= V_cutoff), sqrt_geo_thresh_fit_3d_belowV$thresh[1] + sqrt_geo_thresh_fit_3d_belowV$thresh[2]*sqrt_third_nearest_dist_3d_belowV, col="purple", pch=19)
points(which(gron_eq_cat$V_3d > V_cutoff), sqrt_geo_thresh_fit_3d_aboveV$thresh[1] + sqrt_geo_thresh_fit_3d_aboveV$thresh[2]*sqrt_third_nearest_dist_3d_aboveV, col="orange", pch=19)
legend("bottomright", legend=c("All", "Before", "After", "Below V=5", "Above V=5"), col=c("red", "blue", "green", "purple", "orange"), pch=19)

#V_2d
dev.new(height=12, width=10, noRStudioGD = TRUE)
par(mfrow=c(1,1), bg="transparent")
plot_threshold(gron_eq_cat, geo_thresh_fit_2d, third_nearest_dist_2d, xlab = "V_2d")
points(third_nearest_dist_2d_before, geo_thresh_fit_2d_before$thresh[1] + geo_thresh_fit_2d_before$thresh[2]*third_nearest_dist_2d_before, col="blue", pch=19)
points(third_nearest_dist_2d_after, geo_thresh_fit_2d_after$thresh[1] + geo_thresh_fit_2d_after$thresh[2]*third_nearest_dist_2d_after, col="green", pch=19)
points(third_nearest_dist_2d_belowV, geo_thresh_fit_2d_belowV$thresh[1] + geo_thresh_fit_2d_belowV$thresh[2]*third_nearest_dist_2d_belowV, col="purple", pch=19)
points(third_nearest_dist_2d_aboveV, geo_thresh_fit_2d_aboveV$thresh[1] + geo_thresh_fit_2d_aboveV$thresh[2]*third_nearest_dist_2d_aboveV, col="orange", pch=19)
legend("bottomright", legend=c("All", "Before", "After", "Below V=5", "Above V=5"), col=c("red", "blue", "green", "purple", "orange"), pch=19)


#QQplots 
dev.new(height=20, width=30, noRStudioGD = TRUE)
par(mfrow=c(2,3), bg="transparent")

get_qq_plot(mags_before, geo_thresh_fit_3d_before, third_nearest_dist_3d_before, main="V_3d (before)")
get_qq_plot(mags_before, log_geo_thresh_fit_3d_before, log_third_nearest_dist_3d_before, main="log(V_3d) (before)")
get_qq_plot(mags_before, sqrt_geo_thresh_fit_3d_before, sqrt_third_nearest_dist_3d_before, main="sqrt(V_3d) (before)")

get_qq_plot(mags_after, geo_thresh_fit_3d_after, third_nearest_dist_3d_after, main="V_3d (after)")
get_qq_plot(mags_after, log_geo_thresh_fit_3d_after, log_third_nearest_dist_3d_after, main="log(V_3d) (after)")
get_qq_plot(mags_after, sqrt_geo_thresh_fit_3d_after, sqrt_third_nearest_dist_3d_after, main="sqrt(V_3d) (after)")

dev.new(height=20, width=30, noRStudioGD = TRUE)
par(mfrow=c(2,3), bg="transparent")

get_qq_plot(mags_belowV, geo_thresh_fit_3d_belowV, third_nearest_dist_3d_belowV, main="V_3d (below V=5)")
get_qq_plot(mags_belowV, log_geo_thresh_fit_3d_belowV, log_third_nearest_dist_3d_belowV, main="log(V_3d) (below V=5)")
get_qq_plot(mags_belowV, sqrt_geo_thresh_fit_3d_belowV, sqrt_third_nearest_dist_3d_belowV, main="sqrt(V_3d) (below V=5)")

get_qq_plot(mags_aboveV, geo_thresh_fit_3d_aboveV, third_nearest_dist_3d_aboveV, main="V_3d (above V=5)")
get_qq_plot(mags_aboveV, log_geo_thresh_fit_3d_aboveV, log_third_nearest_dist_3d_aboveV, main="log(V_3d) (above V=5)")
get_qq_plot(mags_aboveV, sqrt_geo_thresh_fit_3d_aboveV, sqrt_third_nearest_dist_3d_aboveV, main="sqrt(V_3d) (above V=5)")

