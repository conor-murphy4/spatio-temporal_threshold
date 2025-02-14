source("src/helper_functions.R")

library(ggplot2)

gron_eq_cat <- read.csv("Data/Events/unrounded_after_1995_in_polygon.csv", header=T)

#Excluding largest V
gron_eq_cat <- gron_eq_cat[gron_eq_cat$V_3d <= 250,]

write.csv(gron_eq_cat, "Data/Events/unrounded_after_1995_in_polygon.csv", row.names=FALSE)

mags <- gron_eq_cat$Magnitude
third_nearest_dist_2d <- gron_eq_cat$V
third_nearest_dist_3d <- gron_eq_cat$V_3d
log_third_nearest_dist_3d <- log(third_nearest_dist_3d)
sqrt_third_nearest_dist_3d <- sqrt(third_nearest_dist_3d)


# Before and after changepoint --------------------------------------------

change_index <- which(gron_eq_cat$Date == "2015-01-06")[1]
gron_eq_cat_before <- gron_eq_cat[1:change_index,]
mags_before <- gron_eq_cat_before$Magnitude
third_nearest_dist_2d_before <- gron_eq_cat_before$V
third_nearest_dist_3d_before <- gron_eq_cat_before$V_3d
log_third_nearest_dist_3d_before <- log(third_nearest_dist_3d_before)
sqrt_third_nearest_dist_3d_before <- sqrt(third_nearest_dist_3d_before)

gron_eq_cat_after <- gron_eq_cat[(change_index+1):nrow(gron_eq_cat),]
mags_after <- gron_eq_cat_after$Magnitude
third_nearest_dist_2d_after <- gron_eq_cat_after$V
third_nearest_dist_3d_after <- gron_eq_cat_after$V_3d
log_third_nearest_dist_3d_after <- log(third_nearest_dist_3d_after)
sqrt_third_nearest_dist_3d_after <- sqrt(third_nearest_dist_3d_after)

# Fitted thresholds -------------------------------------------------------

# EQD threshold
eqd_thresh_fit <- readRDS("threshold_results/eqd_thresh_fit_without_largeV.rds")
eqd_threshold <- rep(eqd_thresh_fit$thresh, length(mags))

eqd_thresh_before <- eqd_threshold[1:change_index]
eqd_thresh_after <- eqd_threshold[(change_index+1):length(mags)]

# Conservative constant threshold from Zak's work -------------------------
conservative_threshold <- rep(1.45, length(mags))

conserv_thresh_before <- conservative_threshold[1:change_index]
conserv_thresh_after <- conservative_threshold[(change_index+1):length(mags)]

# Piece-wise constant threshold according to index ------------------------

piecewise_const_thresh <- c(rep(1.15, change_index), rep(0.76, length(mags) - change_index))

pc_thresh_before <- piecewise_const_thresh[1:change_index]
pc_thresh_after <- piecewise_const_thresh[(change_index+1):length(mags)]


# Fitted stepped threshold ------------------------------------------------
pc_thresh_fit <- readRDS("threshold_results/pc_thresh_fit_without_largeV.rds")
pc_threshold <- c(rep(pc_thresh_fit$thresh_par[1], change_index), rep(pc_thresh_fit$thresh_par[2], length(mags) - change_index))

pc_thresh_fitted_before <- pc_threshold[1:change_index]
pc_thresh_fitted_after <- pc_threshold[(change_index+1):length(mags)]


# V-based thresholds ------------------------------------------------------
geo_thresh_fit_2d_before <- readRDS("threshold_results/geo_thresh_fit_2d_before_without_largeV.rds")
geo_thresh_fit_3d_before <- readRDS("threshold_results/geo_thresh_fit_3d_before_without_largeV.rds")
log_geo_thresh_fit_3d_before <- readRDS("threshold_results/log_geo_thresh_fit_3d_before_without_largeV.rds")
sqrt_geo_thresh_fit_3d_before <- readRDS("threshold_results/sqrt_geo_thresh_fit_3d_before_without_largeV.rds")

geo_thresh_fit_2d_after <- readRDS("threshold_results/geo_thresh_fit_2d_after_without_largeV.rds")
geo_thresh_fit_3d_after <- readRDS("threshold_results/geo_thresh_fit_3d_after_without_largeV.rds")
log_geo_thresh_fit_3d_after <- readRDS("threshold_results/log_geo_thresh_fit_3d_after_without_largeV.rds")
sqrt_geo_thresh_fit_3d_after <- readRDS("threshold_results/sqrt_geo_thresh_fit_3d_after_without_largeV.rds")


#Confidence intervals for parameters

get_par_ests_step(mags_before, conserv_thresh_before)
get_par_ests_step(mags_before, eqd_thresh_before)
get_par_ests_step(mags_before, pc_thresh_before)
get_par_ests_step(mags_before, pc_thresh_fitted_before)

get_par_ests_geo(mags_before, geo_thresh_fit_2d_before, third_nearest_dist_2d_before)
get_par_ests_geo(mags_before, geo_thresh_fit_3d_before, third_nearest_dist_3d_before)
get_par_ests_geo(mags_before, log_geo_thresh_fit_3d_before, log_third_nearest_dist_3d_before)
get_par_ests_geo(mags_before, sqrt_geo_thresh_fit_3d_before, sqrt_third_nearest_dist_3d_before)

get_par_ests_step(mags_after, conserv_thresh_after)
get_par_ests_step(mags_after, eqd_thresh_after)
get_par_ests_step(mags_after, pc_thresh_after)
get_par_ests_step(mags_after, pc_thresh_fitted_after)

get_par_ests_geo(mags_after, geo_thresh_fit_2d_after, third_nearest_dist_2d_after)
get_par_ests_geo(mags_after, geo_thresh_fit_3d_after, third_nearest_dist_3d_after)
get_par_ests_geo(mags_after, log_geo_thresh_fit_3d_after, log_third_nearest_dist_3d_after)
get_par_ests_geo(mags_after, sqrt_geo_thresh_fit_3d_after, sqrt_third_nearest_dist_3d_after)


# EQD values --------------------------------------------------------------

#Before
set.seed(11111)
get_eqd_value(mags_before, conserv_thresh_before, par = get_par_ests_step(mags_before, conserv_thresh_before)$par, step = TRUE)

set.seed(11111)
get_eqd_value(mags_before, eqd_thresh_before, par = get_par_ests_step(mags_before, eqd_thresh_before)$par, step = TRUE)

set.seed(11111)
get_eqd_value(mags_before, pc_thresh_before, par = get_par_ests_step(mags_before, pc_thresh_before)$par, step = TRUE)

set.seed(11111)
get_eqd_value(mags_before, pc_thresh_fitted_before, par = get_par_ests_step(mags_before, pc_thresh_fitted_before)$par, step = TRUE)

min(geo_thresh_fit_2d_before$dists, na.rm=TRUE)
min(geo_thresh_fit_3d_before$dists, na.rm=TRUE)
min(log_geo_thresh_fit_3d_before$dists, na.rm=TRUE)
min(sqrt_geo_thresh_fit_3d_before$dists, na.rm=TRUE)

#After
set.seed(11111)
get_eqd_value(mags_after, conserv_thresh_after, par = get_par_ests_step(mags_after, conserv_thresh_after)$par, step = TRUE)

set.seed(11111)
get_eqd_value(mags_after, eqd_thresh_after, par = get_par_ests_step(mags_after, eqd_thresh_after)$par, step = TRUE)

set.seed(11111)
get_eqd_value(mags_after, pc_thresh_after, par = get_par_ests_step(mags_after, pc_thresh_after)$par, step = TRUE)

set.seed(11111)
get_eqd_value(mags_after, pc_thresh_fitted_after, par = get_par_ests_step(mags_after, pc_thresh_fitted_after)$par, step = TRUE)

min(geo_thresh_fit_2d_after$dists, na.rm=TRUE)
min(geo_thresh_fit_3d_after$dists, na.rm=TRUE)
min(log_geo_thresh_fit_3d_after$dists, na.rm=TRUE)
min(sqrt_geo_thresh_fit_3d_after$dists, na.rm=TRUE)

# QQ plots

#QQplots
dev.new(height=20, width=40, noRStudioGD = TRUE)
par(mfrow=c(2,4), bg='transparent')
get_qq_plot_const(mags_before, eqd_thresh_before, main="EQD thresh (before)")
get_qq_plot_const(mags_before, conserv_thresh_before, main="Conservative (before)")
get_qq_plot_const(mags_before, pc_thresh_before, main="Zak's stepped (before)")
get_qq_plot_const(mags_before, pc_thresh_fitted_before, main="Fitted stepped (before)")

get_qq_plot_geo(mags_before, geo_thresh_fit_2d_before, third_nearest_dist_2d_before, main="V_2d (before)")
get_qq_plot_geo(mags_before, geo_thresh_fit_3d_before, third_nearest_dist_3d_before, main="V_3d (before)" )
get_qq_plot_geo(mags_before, log_geo_thresh_fit_3d_before, log_third_nearest_dist_3d_before, main="log(V_3d) (before)")
get_qq_plot_geo(mags_before, sqrt_geo_thresh_fit_3d_before, sqrt_third_nearest_dist_3d_before, main="sqrt(V_3d) (before)")

dev.new(height=20, width=40, noRStudioGD = TRUE)
par(mfrow=c(2,4), bg='transparent')
get_qq_plot_const(mags_after, eqd_thresh_after, main="EQD thresh (after)")
get_qq_plot_const(mags_after, conserv_thresh_after, main="Conservative (after)")
get_qq_plot_const(mags_after, pc_thresh_after, main="Zak's stepped (after)")
get_qq_plot_const(mags_after, pc_thresh_fitted_after, main="Fitted stepped (after)")

get_qq_plot_geo(mags_after, geo_thresh_fit_2d_after, third_nearest_dist_2d_after, main="V_2d (after)")
get_qq_plot_geo(mags_after, geo_thresh_fit_3d_after, third_nearest_dist_3d_after, main="V_3d (after)" )
get_qq_plot_geo(mags_after, log_geo_thresh_fit_3d_after, log_third_nearest_dist_3d_after, main="log(V_3d) (after)")
get_qq_plot_geo(mags_after, sqrt_geo_thresh_fit_3d_after, sqrt_third_nearest_dist_3d_after, main="sqrt(V_3d) (after)")


#Image plots of dists

dev.new(height=20, width=20, noRStudioGD = TRUE)
par(mfrow=c(1,1), bg='transparent')

intercepts <- seq(0, 1.5, by=0.02)
slopes <- seq(0,2.5, by=0.02)
threshold_matrix <- as.matrix(expand.grid(intercepts, slopes))


image_plot_df <- data.frame(a=threshold_matrix[,1], b=threshold_matrix[,2], v_dists=geo_thresh_fit_3d$dists, logv_dists=log_geo_thresh_fit_3d$dists, sqrtv_dists=sqrt_geo_thresh_fit_3d$dists)

ggplot(image_plot_df, aes(x=a, y=b, fill=v_dists)) + geom_tile() + scale_fill_viridis_c() + ggtitle("EQD V") + theme_minimal()
ggplot(image_plot_df, aes(x=a, y=b, fill=logv_dists)) + geom_tile() + scale_fill_viridis_c() + ggtitle("EQD logV") + theme_minimal()
ggplot(image_plot_df, aes(x=a, y=b, fill=sqrtv_dists)) + geom_tile() + scale_fill_viridis_c() + ggtitle("EQD sqrtV") + theme_minimal()


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
points(third_nearest_dist_belowV_3d, geo_thresh_fit_belowV_3d$thresh[1] + geo_thresh_fit_belowV_3d$thresh[2]*third_nearest_dist_belowV_3d, col="purple", pch=19)
points(third_nearest_dist_aboveV_3d, geo_thresh_fit_aboveV_3d$thresh[1] + geo_thresh_fit_aboveV_3d$thresh[2]*third_nearest_dist_aboveV_3d, col="orange", pch=19)
legend("bottomright", legend=c("All", "Before", "After", "Below V=V_cutoff", "Above V=V_cutoff"), col=c("red", "blue", "green", "purple", "orange"), pch=19)

plot_threshold(gron_eq_cat, log_geo_thresh_fit_3d, log_third_nearest_dist_3d, xlab = "log(V_3d)")
points(log_third_nearest_dist_3d_before, log_geo_thresh_fit_3d_before$thresh[1] + log_geo_thresh_fit_3d_before$thresh[2]*log_third_nearest_dist_3d_before, col="blue", pch=19)
points(log_third_nearest_dist_3d_after, log_geo_thresh_fit_3d_after$thresh[1] + log_geo_thresh_fit_3d_after$thresh[2]*log_third_nearest_dist_3d_after, col="green", pch=19)
points(log_third_nearest_dist_belowV_3d, log_geo_thresh_fit_belowV_3d$thresh[1] + log_geo_thresh_fit_belowV_3d$thresh[2]*log_third_nearest_dist_belowV_3d, col="purple", pch=19)
points(log_third_nearest_dist_aboveV_3d, log_geo_thresh_fit_aboveV_3d$thresh[1] + log_geo_thresh_fit_aboveV_3d$thresh[2]*log_third_nearest_dist_aboveV_3d, col="orange", pch=19)
legend("bottomright", legend=c("All", "Before", "After", "Below V=V_cutoff", "Above V=V_cutoff"), col=c("red", "blue", "green", "purple", "orange"), pch=19)

plot_threshold(gron_eq_cat, sqrt_geo_thresh_fit_3d, sqrt_third_nearest_dist_3d, xlab = "sqrt(V_3d)")
points(sqrt_third_nearest_dist_3d_before, sqrt_geo_thresh_fit_3d_before$thresh[1] + sqrt_geo_thresh_fit_3d_before$thresh[2]*sqrt_third_nearest_dist_3d_before, col="blue", pch=19)
points(sqrt_third_nearest_dist_3d_after, sqrt_geo_thresh_fit_3d_after$thresh[1] + sqrt_geo_thresh_fit_3d_after$thresh[2]*sqrt_third_nearest_dist_3d_after, col="green", pch=19)
points(sqrt_third_nearest_dist_belowV_3d, sqrt_geo_thresh_fit_belowV_3d$thresh[1] + sqrt_geo_thresh_fit_belowV_3d$thresh[2]*sqrt_third_nearest_dist_belowV_3d, col="purple", pch=19)
points(sqrt_third_nearest_dist_aboveV_3d, sqrt_geo_thresh_fit_aboveV_3d$thresh[1] + sqrt_geo_thresh_fit_aboveV_3d$thresh[2]*sqrt_third_nearest_dist_aboveV_3d, col="orange", pch=19)
legend("bottomright", legend=c("All", "Before", "After", "Below V=V_cutoff", "Above V=V_cutoff"), col=c("red", "blue", "green", "purple", "orange"), pch=19)

#Against index
dev.new(height=12, width=30, noRStudioGD = TRUE)
par(mfrow=c(1,3), bg="transparent")

plot(gron_eq_cat$Magnitude, ylab = "Magnitude", xlab="Index", main="V_3d")
points(geo_thresh_fit_3d$thresh[1] + geo_thresh_fit_3d$thresh[2]*third_nearest_dist_3d, col="red", pch=19)
points(c(1:change_index),geo_thresh_fit_3d_before$thresh[1] + geo_thresh_fit_3d_before$thresh[2]*third_nearest_dist_3d_before, col="blue", pch=19)
points(c((change_index+1):nrow(gron_eq_cat)),geo_thresh_fit_3d_after$thresh[1] + geo_thresh_fit_3d_after$thresh[2]*third_nearest_dist_3d_after, col="green", pch=19)
points(which(gron_eq_cat$V_3d <= V_cutoff), geo_thresh_fit_belowV_3d$thresh[1] + geo_thresh_fit_belowV_3d$thresh[2]*third_nearest_dist_belowV_3d, col="purple", pch=19)
points(which(gron_eq_cat$V_3d > V_cutoff), geo_thresh_fit_aboveV_3d$thresh[1] + geo_thresh_fit_aboveV_3d$thresh[2]*third_nearest_dist_aboveV_3d, col="orange", pch=19)
legend("bottomright", legend=c("All", "Before", "After", "Below V=V_cutoff", "Above V=V_cutoff"), col=c("red", "blue", "green", "purple", "orange"), pch=19)

plot(gron_eq_cat$Magnitude, ylab = "Magnitude", xlab="Index", main="log(V_3d)")
points(log_geo_thresh_fit_3d$thresh[1] + log_geo_thresh_fit_3d$thresh[2]*log_third_nearest_dist_3d, col="red", pch=19)
points(c(1:change_index),log_geo_thresh_fit_3d_before$thresh[1] + log_geo_thresh_fit_3d_before$thresh[2]*log_third_nearest_dist_3d_before, col="blue", pch=19)
points(c((change_index+1):nrow(gron_eq_cat)),log_geo_thresh_fit_3d_after$thresh[1] + log_geo_thresh_fit_3d_after$thresh[2]*log_third_nearest_dist_3d_after, col="green", pch=19)
points(which(gron_eq_cat$V_3d <= V_cutoff), log_geo_thresh_fit_belowV_3d$thresh[1] + log_geo_thresh_fit_belowV_3d$thresh[2]*log_third_nearest_dist_belowV_3d, col="purple", pch=19)
points(which(gron_eq_cat$V_3d > V_cutoff), log_geo_thresh_fit_aboveV_3d$thresh[1] + log_geo_thresh_fit_aboveV_3d$thresh[2]*log_third_nearest_dist_aboveV_3d, col="orange", pch=19)
legend("bottomright", legend=c("All", "Before", "After", "Below V=V_cutoff", "Above V=V_cutoff"), col=c("red", "blue", "green", "purple", "orange"), pch=19)

plot(gron_eq_cat$Magnitude, ylab = "Magnitude", xlab="Index", main="sqrt(V_3d)")
points(sqrt_geo_thresh_fit_3d$thresh[1] + sqrt_geo_thresh_fit_3d$thresh[2]*sqrt_third_nearest_dist_3d, col="red", pch=19)
points(c(1:change_index),sqrt_geo_thresh_fit_3d_before$thresh[1] + sqrt_geo_thresh_fit_3d_before$thresh[2]*sqrt_third_nearest_dist_3d_before, col="blue", pch=19)
points(c((change_index+1):nrow(gron_eq_cat)), sqrt_geo_thresh_fit_3d_after$thresh[1] + sqrt_geo_thresh_fit_3d_after$thresh[2]*sqrt_third_nearest_dist_3d_after, col="green", pch=19)
points(which(gron_eq_cat$V_3d <= V_cutoff), sqrt_geo_thresh_fit_belowV_3d$thresh[1] + sqrt_geo_thresh_fit_belowV_3d$thresh[2]*sqrt_third_nearest_dist_belowV_3d, col="purple", pch=19)
points(which(gron_eq_cat$V_3d > V_cutoff), sqrt_geo_thresh_fit_aboveV_3d$thresh[1] + sqrt_geo_thresh_fit_aboveV_3d$thresh[2]*sqrt_third_nearest_dist_aboveV_3d, col="orange", pch=19)
legend("bottomright", legend=c("All", "Before", "After", "Below V=V_cutoff", "Above V=V_cutoff"), col=c("red", "blue", "green", "purple", "orange"), pch=19)

#V_2d
dev.new(height=12, width=10, noRStudioGD = TRUE)
par(mfrow=c(1,1), bg="transparent")
plot_threshold(gron_eq_cat, geo_thresh_fit_2d, third_nearest_dist_2d, xlab = "V_2d")
points(third_nearest_dist_2d_before, geo_thresh_fit_2d_before$thresh[1] + geo_thresh_fit_2d_before$thresh[2]*third_nearest_dist_2d_before, col="blue", pch=19)
points(third_nearest_dist_2d_after, geo_thresh_fit_2d_after$thresh[1] + geo_thresh_fit_2d_after$thresh[2]*third_nearest_dist_2d_after, col="green", pch=19)
points(third_nearest_dist_belowV_2d, geo_thresh_fit_belowV_2d$thresh[1] + geo_thresh_fit_belowV_2d$thresh[2]*third_nearest_dist_belowV_2d, col="purple", pch=19)
points(third_nearest_dist_aboveV_2d, geo_thresh_fit_aboveV_2d$thresh[1] + geo_thresh_fit_aboveV_2d$thresh[2]*third_nearest_dist_aboveV_2d, col="orange", pch=19)
legend("bottomright", legend=c("All", "Before", "After", "Below V=V_cutoff", "Above V=V_cutoff"), col=c("red", "blue", "green", "purple", "orange"), pch=19)




