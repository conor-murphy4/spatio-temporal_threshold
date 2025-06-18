library(ggplot2)
library(dplyr)
source("src/helper_functions.R")


gron_eq_cat <- read.csv("Data/Events/unrounded_after_1995_in_polygon_with_covariates.csv", header=T)
mags <- gron_eq_cat$Magnitude
nearest_dist_matrix <- matrix(c(gron_eq_cat$V_1, gron_eq_cat$V_2, gron_eq_cat$V_3, gron_eq_cat$V_4), byrow=F, ncol=4)
covariates <- read.csv("Data/covariates/covariates_1995-2024.csv", header=T)

# Most recent threshold results

#ICS max (B=200) - Results shown in paper
thresh_fit_V1_ics <- readRDS("threshold_results/icsmax/geo_thresh_fit_V1.rds")
thresh_fit_V2_ics <- readRDS("threshold_results/icsmax/geo_thresh_fit_V2.rds")
thresh_fit_V3_ics <- readRDS("threshold_results/icsmax/geo_thresh_fit_V3.rds")
thresh_fit_V4_ics <- readRDS("threshold_results/icsmax/geo_thresh_fit_V4.rds")

thresh_fit_logV1_ics <- readRDS("threshold_results/icsmax/geo_thresh_fit_logV1.rds")
thresh_fit_logV2_ics <- readRDS("threshold_results/icsmax/geo_thresh_fit_logV2.rds")
thresh_fit_logV3_ics <- readRDS("threshold_results/icsmax/geo_thresh_fit_logV3.rds")
thresh_fit_logV4_ics <- readRDS("threshold_results/icsmax/geo_thresh_fit_logV4.rds")

thresh_fit_sqrtV1_ics <- readRDS("threshold_results/icsmax/geo_thresh_fit_sqrtV1.rds")
thresh_fit_sqrtV2_ics <- readRDS("threshold_results/icsmax/geo_thresh_fit_sqrtV2.rds")
thresh_fit_sqrtV3_ics <- readRDS("threshold_results/icsmax/geo_thresh_fit_sqrtV3.rds")
thresh_fit_sqrtV4_ics <- readRDS("threshold_results/icsmax/geo_thresh_fit_sqrtV4.rds")

# B=1000
# thresh_fit_V1_ics <- readRDS("threshold_results/icsmax/boot1000/geo_thresh_fit_V1.rds")
# thresh_fit_V2_ics <- readRDS("threshold_results/icsmax/boot1000/geo_thresh_fit_V2.rds")
# thresh_fit_V3_ics <- readRDS("threshold_results/icsmax/boot1000/geo_thresh_fit_V3.rds")
# thresh_fit_V4_ics <- readRDS("threshold_results/icsmax/boot1000/geo_thresh_fit_V4.rds")
# 
# thresh_fit_logV1_ics <- readRDS("threshold_results/icsmax/boot1000/geo_thresh_fit_logV1.rds")
# thresh_fit_logV2_ics <- readRDS("threshold_results/icsmax/boot1000/geo_thresh_fit_logV2.rds")
# thresh_fit_logV3_ics <- readRDS("threshold_results/icsmax/boot1000/geo_thresh_fit_logV3.rds")
# thresh_fit_logV4_ics <- readRDS("threshold_results/icsmax/boot1000/geo_thresh_fit_logV4.rds")
# 
# thresh_fit_sqrtV1_ics <- readRDS("threshold_results/icsmax/boot1000/geo_thresh_fit_sqrtV1.rds")
# thresh_fit_sqrtV2_ics <- readRDS("threshold_results/icsmax/boot1000/geo_thresh_fit_sqrtV2.rds")
# thresh_fit_sqrtV3_ics <- readRDS("threshold_results/icsmax/boot1000/geo_thresh_fit_sqrtV3.rds")
# thresh_fit_sqrtV4_ics <- readRDS("threshold_results/icsmax/boot1000/geo_thresh_fit_sqrtV4.rds")

# Minimum dists values
min(thresh_fit_V1_ics$dists, na.rm = TRUE)
min(thresh_fit_V2_ics$dists, na.rm = TRUE)
min(thresh_fit_V3_ics$dists, na.rm = TRUE)
min(thresh_fit_V4_ics$dists, na.rm = TRUE)

min(thresh_fit_logV1_ics$dists, na.rm = TRUE)
min(thresh_fit_logV2_ics$dists, na.rm = TRUE)
min(thresh_fit_logV3_ics$dists, na.rm = TRUE)
min(thresh_fit_logV4_ics$dists, na.rm = TRUE)

min(thresh_fit_sqrtV1_ics$dists, na.rm = TRUE)
min(thresh_fit_sqrtV2_ics$dists, na.rm = TRUE)
min(thresh_fit_sqrtV3_ics$dists, na.rm = TRUE)
min(thresh_fit_sqrtV4_ics$dists, na.rm = TRUE)

#Fitted thresholds
eqd_thresh_fit <- readRDS("threshold_results/eqd_thresh_fit.rds")

eqd_threshold <- rep(eqd_thresh_fit$thresh, length(mags))

conservative_threshold <- rep(1.45, length(mags))

u_h_length <- which(gron_eq_cat$Date == as.Date("2015-01-06"))[1]
piecewise_const_thresh <- c(rep(1.15, u_h_length), rep(0.76, length(mags) - u_h_length))

pc_thresh_fit <- readRDS("threshold_results/pc_thresh_fit.rds")
pc_fitted_threshold <- c(rep(pc_thresh_fit$thresh[1], u_h_length), rep(pc_thresh_fit$thresh[2], length(mags) - u_h_length))

geo_thresh_fit <- vector("list", 4)
log_geo_thresh_fit <- vector("list", 4)
sqrt_geo_thresh_fit <- vector("list", 4)
for(i in 1:4){
  geo_thresh_fit[[i]] <- readRDS(paste0("threshold_results/icsmax/geo_thresh_fit_V", i, ".rds"))
  log_geo_thresh_fit[[i]] <- readRDS(paste0("threshold_results/icsmax/geo_thresh_fit_logV", i, ".rds"))
  sqrt_geo_thresh_fit[[i]] <- readRDS(paste0("threshold_results/icsmax/geo_thresh_fit_sqrtV", i, ".rds"))
}

# Results for the constant and stepped thresholds

get_table_of_results(mags, conservative_threshold)
get_table_of_results(mags, eqd_threshold)
get_table_of_results(mags, piecewise_const_thresh)
get_table_of_results(mags, pc_fitted_threshold)

#Results for the geophone-based thresholds

min_ics <- min(covariates$ICS_max, na.rm = TRUE)
current_covariates <- covariates %>% filter(Year == max(covariates$Year))
max_dists <- c(max(current_covariates$V1), max(current_covariates$V2), 
               max(current_covariates$V3), max(current_covariates$V4, na.rm=T))

get_table_of_results_geo_ics(mags, geo_thresh_fit, nearest_dist_matrix, gron_eq_cat$ICS_max, 
                            min_ics, max_dists)
get_table_of_results_geo_ics(mags, log_geo_thresh_fit, log(nearest_dist_matrix), gron_eq_cat$ICS_max, 
                            min_ics, log(max_dists))
get_table_of_results_geo_ics(mags, sqrt_geo_thresh_fit, sqrt(nearest_dist_matrix), gron_eq_cat$ICS_max,
                            min_ics, sqrt(max_dists))



                               
#QQplots
dev.new(height=40, width=40, noRStudioGD = TRUE)
par(mfrow=c(4,4), bg='transparent')
get_qq_plot_const(mags, eqd_threshold, main="EQD threshold")
get_qq_plot_const(mags, conservative_threshold, main="Conservative threshold")
get_qq_plot_const(mags, piecewise_const_thresh, main="Zak's stepped threshold")
get_qq_plot_const(mags, pc_fitted_threshold, main="Fitted stepped threshold")

#Change the below for each variation of V-based threshold, i.e. V, logV, sqrtV
list_of_fits <- sqrt_geo_thresh_fit
dist_matrix <- sqrt(third_nearest_dist_matrix)
for(i in 1:4){
  get_qq_plot_geo(mags, list_of_fits[[i]], dist_matrix[,i], main=paste0("Threshold on sqrtV", i))
}

# Visualising thresholds

dev.new(width=20, height=20,noRStudioGD = TRUE)
par(mfrow=c(2,2),bg='transparent')

plot(gron_eq_cat$Magnitude, xlab = "Index", ylab = "Magnitude", main = "Constant & Stepped", col="grey", pch=19)
lines(piecewise_const_thresh, col="red", lw=2)
lines(pc_fitted_threshold, col="orange", lw=2)
lines(conservative_threshold, col="blue", lw=2)
lines(eqd_threshold, col="green", lw=2)
legend("topright", legend=c("Zak's stepped", "Fitted stepped", "Conservative", "EQD constant"), col=c("red","orange", "blue", "green"), lty=1, cex=0.8)

# Visualising the geophone-based thresholds #Try to automate below code so no manual changes are needed...
list_of_fits <- sqrt_geo_thresh_fit
dist_matrix <- sqrt(third_nearest_dist_matrix)
plot(gron_eq_cat$Magnitude, xlab = "Index", ylab = "Magnitude", main = "sqrtV-based", col="grey", pch=19)
titles <- vector("character", 4)
for(i in 1:4){
  lines(list_of_fits[[i]]$thresh_par[1] + list_of_fits[[i]]$thresh_par[2]*dist_matrix[,i], col=i, lw=2)
  titles[i] <- paste0("Threshold on sqrtV", i)
}
legend("topright", legend=titles, col=c(1,2,3,4), lty=1, cex=0.8)

#Image plot for EQD values across thresh_matrix

dev.new(height=20, width=20, noRStudioGD = TRUE)
par(mfrow=c(1,1), bg='transparent')

get_image_plots <- function(threshold_matrix, list_of_fits){
  for(i in 1:4){
    dev.new(height=20, width=20, noRStudioGD = TRUE)
    par(mfrow=c(1,1), bg='transparent')
    dists <- list_of_fits[[i]]$dists
    image_plot_df <- data.frame(a=threshold_matrix[,1], b=threshold_matrix[,2], dists=dists)
    print(ggplot(image_plot_df, aes(x=a, y=b, fill=dists)) + geom_tile() + scale_fill_viridis_c() + ggtitle(paste("EQD V", i)) + theme_minimal())
  }
}

# NOTE: Adjust the below threshold matrix as currently, the slopes go too high
intercepts <- seq(0, 1.5, by=0.02)
slopes <- seq(0,1, by=0.02)
threshold_matrix <- as.matrix(expand.grid(intercepts, slopes))

get_image_plots(threshold_matrix, geo_thresh_fit)

image_plot_df <- data.frame(a=threshold_matrix[,1], b=threshold_matrix[,2], dists=thresh_fit_geo_ics$dists)
ggplot(image_plot_df, aes(x=a, y=b, fill=dists)) + geom_tile() + scale_fill_viridis_c() + ggtitle("EQD V1 with ICS") + theme_minimal()

# Plot of profile EQD plots foe each intercept

get_profile_eqd <- function(threshold_matrix, list_of_fits, margin = 1){
  prof_dists <- matrix(NA, nrow=length(unique(threshold_matrix[,margin])), ncol=4)
  par_vals <- matrix(NA, nrow=length(unique(threshold_matrix[,margin])), ncol=4)
  for(i in 1:4){
    dists <- list_of_fits[[i]]$dists
    for(j in 1:length(unique(threshold_matrix[,margin]))){
      prof_axis_val <- unique(threshold_matrix[,margin])[j]
      corresponding_dists <- dists[threshold_matrix[,margin] == prof_axis_val]
      if(all(is.na(corresponding_dists))){
        prof_dists[j, i] <- NA
        par_vals[j, i] <- NA
      }
      else{
        corresponding_other_par <- threshold_matrix[threshold_matrix[,margin] == prof_axis_val, -margin]
        min_index <- which.min(corresponding_dists)
        prof_dists[j, i] <- corresponding_dists[min_index]
        par_vals[j, i] <- corresponding_other_par[min_index]
      }
    }
  }
  return(list(profile_axis = unique(threshold_matrix[,margin]), eqd_dists = prof_dists, corresponding_pars = par_vals))
}

profile_over_a_geo <- get_profile_eqd(threshold_matrix, geo_thresh_fit)
profile_over_b_geo <- get_profile_eqd(threshold_matrix, geo_thresh_fit, margin=2)

profile_over_a_log_geo <- get_profile_eqd(threshold_matrix, log_geo_thresh_fit)
profile_over_b_log_geo <- get_profile_eqd(threshold_matrix, log_geo_thresh_fit, margin=2)

profile_over_a_sqrt_geo <- get_profile_eqd(threshold_matrix, sqrt_geo_thresh_fit)
profile_over_b_sqrt_geo <- get_profile_eqd(threshold_matrix, sqrt_geo_thresh_fit, margin=2)

# Plot the results

# Profile over intercept
dev.new(height=30, width=20, noRStudioGD = TRUE)
par(mfrow=c(3,2), bg='transparent')

lo1 <- loess(profile_over_a_geo$eqd_dists[,1] ~ profile_over_a_geo$profile_axis)
lo2 <- loess(profile_over_a_geo$eqd_dists[,2] ~ profile_over_a_geo$profile_axis)
lo3 <- loess(profile_over_a_geo$eqd_dists[,3] ~ profile_over_a_geo$profile_axis)
lo4 <- loess(profile_over_a_geo$eqd_dists[,4] ~ profile_over_a_geo$profile_axis)
plot(profile_over_a_geo$profile_axis, predict(lo1), type="l", col="red", xlab="Intercept", ylab="EQD", main = "V")
lines(profile_over_a_geo$profile_axis, predict(lo2), col="blue")
lines(profile_over_a_geo$profile_axis, predict(lo3), col="green")
lines(profile_over_a_geo$profile_axis, predict(lo4), col="black")
legend("topright", legend=c("V1", "V2", "V3", "V4"), col=c("red", "blue", "green", "black"), lty=1, cex=0.8)
abline(h=0.0576, col = "grey", lty=2, lwd=2)
abline(h=0.0365, col = "grey", lty=2, lwd=2)
abline(h=0.0322, col = "grey", lty=2, lwd=2)
abline(h=0.0304, col = "grey", lty=2, lwd=2)

plot(profile_over_a_geo$profile_axis, profile_over_a_geo$corresponding_pars[,1], type="l", col="red", xlab="Intercept", ylab="Slope", main = "V")
lines(profile_over_a_geo$profile_axis, profile_over_a_geo$corresponding_pars[,2], col="blue")
lines(profile_over_a_geo$profile_axis, profile_over_a_geo$corresponding_pars[,3], col="green")
lines(profile_over_a_geo$profile_axis, profile_over_a_geo$corresponding_pars[,4], col="black")
legend("topright", legend=c("V1", "V2", "V3", "V4"), col=c("red", "blue", "green", "black"), lty=1, cex=0.8)

lo1 <- loess(profile_over_a_log_geo$eqd_dists[,1] ~ profile_over_a_log_geo$profile_axis)
lo2 <- loess(profile_over_a_log_geo$eqd_dists[,2] ~ profile_over_a_log_geo$profile_axis)
lo3 <- loess(profile_over_a_log_geo$eqd_dists[,3] ~ profile_over_a_log_geo$profile_axis)
lo4 <- loess(profile_over_a_log_geo$eqd_dists[,4] ~ profile_over_a_log_geo$profile_axis)
plot(profile_over_a_log_geo$profile_axis, predict(lo1), type="l", col="red", xlab="Intercept", ylab="EQD", main = "logV")
lines(profile_over_a_log_geo$profile_axis, predict(lo2), col="blue")
lines(profile_over_a_log_geo$profile_axis, predict(lo3), col="green")
lines(profile_over_a_log_geo$profile_axis, predict(lo4), col="black")
legend("topright", legend=c("V1", "V2", "V3", "V4"), col=c("red", "blue", "green", "black"), lty=1, cex=0.8)
abline(h=0.0576, col = "grey", lty=2, lwd=2)
abline(h=0.0365, col = "grey", lty=2, lwd=2)
abline(h=0.0322, col = "grey", lty=2, lwd=2)
abline(h=0.0304, col = "grey", lty=2, lwd=2)

plot(profile_over_a_log_geo$profile_axis, profile_over_a_log_geo$corresponding_pars[,1], type="l", col="red", xlab="Intercept", ylab="Slope", main = "V")
lines(profile_over_a_log_geo$profile_axis, profile_over_a_log_geo$corresponding_pars[,2], col="blue")
lines(profile_over_a_log_geo$profile_axis, profile_over_a_log_geo$corresponding_pars[,3], col="green")
lines(profile_over_a_log_geo$profile_axis, profile_over_a_log_geo$corresponding_pars[,4], col="black")
legend("topright", legend=c("V1", "V2", "V3", "V4"), col=c("red", "blue", "green", "black"), lty=1, cex=0.8)

lo1 <- loess(profile_over_a_sqrt_geo$eqd_dists[,1] ~ profile_over_a_sqrt_geo$profile_axis)
lo2 <- loess(profile_over_a_sqrt_geo$eqd_dists[,2] ~ profile_over_a_sqrt_geo$profile_axis)
lo3 <- loess(profile_over_a_sqrt_geo$eqd_dists[,3] ~ profile_over_a_sqrt_geo$profile_axis)
lo4 <- loess(profile_over_a_sqrt_geo$eqd_dists[,4] ~ profile_over_a_sqrt_geo$profile_axis)
plot(profile_over_a_sqrt_geo$profile_axis, predict(lo1), type="l", col="red", xlab="Intercept", ylab="EQD", main = "sqrtV")
lines(profile_over_a_sqrt_geo$profile_axis, predict(lo2), col="blue")
lines(profile_over_a_sqrt_geo$profile_axis, predict(lo3), col="green")
lines(profile_over_a_sqrt_geo$profile_axis, predict(lo4), col="black")
legend("topright", legend=c("V1", "V2", "V3", "V4"), col=c("red", "blue", "green", "black"), lty=1, cex=0.8)
abline(h=0.0576, col = "grey", lty=2, lwd=2)
abline(h=0.0365, col = "grey", lty=2, lwd=2)
abline(h=0.0322, col = "grey", lty=2, lwd=2)
abline(h=0.0304, col = "grey", lty=2, lwd=2)

plot(profile_over_a_sqrt_geo$profile_axis, profile_over_a_sqrt_geo$corresponding_pars[,1], type="l", col="red", xlab="Intercept", ylab="Slope", main = "sqrtV")
lines(profile_over_a_sqrt_geo$profile_axis, profile_over_a_sqrt_geo$corresponding_pars[,2], col="blue")
lines(profile_over_a_sqrt_geo$profile_axis, profile_over_a_sqrt_geo$corresponding_pars[,3], col="green")
lines(profile_over_a_sqrt_geo$profile_axis, profile_over_a_sqrt_geo$corresponding_pars[,4], col="black")
legend("topright", legend=c("V1", "V2", "V3", "V4"), col=c("red", "blue", "green", "black"), lty=1, cex=0.8)

# Profile over slope
lo1 <- loess(profile_over_b_geo$eqd_dists[,1] ~ profile_over_b_geo$profile_axis)
lo2 <- loess(profile_over_b_geo$eqd_dists[,2] ~ profile_over_b_geo$profile_axis)
lo3 <- loess(profile_over_b_geo$eqd_dists[,3] ~ profile_over_b_geo$profile_axis)
lo4 <- loess(profile_over_b_geo$eqd_dists[,4] ~ profile_over_b_geo$profile_axis)
plot(profile_over_b_geo$profile_axis, predict(lo1,profile_over_b_geo$profile_axis ), type="l", col="red", xlab="Slope", ylab="EQD", main = "V")
lines(profile_over_b_geo$profile_axis, predict(lo2, profile_over_b_geo$profile_axis), col="blue")
lines(profile_over_b_geo$profile_axis, predict(lo3, profile_over_b_geo$profile_axis), col="green")
lines(profile_over_b_geo$profile_axis, predict(lo4, profile_over_b_geo$profile_axis), col="black")
legend("topright", legend=c("V1", "V2", "V3", "V4"), col=c("red", "blue", "green", "black"), lty=1, cex=0.8)
abline(h=0.0576, col = "grey", lty=2, lwd=2)
abline(h=0.0365, col = "grey", lty=2, lwd=2)
abline(h=0.0322, col = "grey", lty=2, lwd=2)
abline(h=0.0304, col = "grey", lty=2, lwd=2)

plot(profile_over_b_geo$profile_axis, profile_over_b_geo$corresponding_pars[,1], type="l", col="red", xlab="Slope", ylab="Intercept", main = "V")
lines(profile_over_b_geo$profile_axis, profile_over_b_geo$corresponding_pars[,2], col="blue")
lines(profile_over_b_geo$profile_axis, profile_over_b_geo$corresponding_pars[,3], col="green")
lines(profile_over_b_geo$profile_axis, profile_over_b_geo$corresponding_pars[,4], col="black")
legend("topright", legend=c("V1", "V2", "V3", "V4"), col=c("red", "blue", "green", "black"), lty=1, cex=0.8)

lo1 <- loess(profile_over_b_log_geo$eqd_dists[,1] ~ profile_over_b_log_geo$profile_axis)
lo2 <- loess(profile_over_b_log_geo$eqd_dists[,2] ~ profile_over_b_log_geo$profile_axis)
lo3 <- loess(profile_over_b_log_geo$eqd_dists[,3] ~ profile_over_b_log_geo$profile_axis)
lo4 <- loess(profile_over_b_log_geo$eqd_dists[,4] ~ profile_over_b_log_geo$profile_axis)
plot(profile_over_b_log_geo$profile_axis, predict(lo1, profile_over_b_log_geo$profile_axis), type="l", col="red", xlab="Slope", ylab="EQD", main = "logV")
lines(profile_over_b_log_geo$profile_axis, predict(lo2, profile_over_b_log_geo$profile_axis), col="blue")
lines(profile_over_b_log_geo$profile_axis, predict(lo3, profile_over_b_log_geo$profile_axis), col="green")
lines(profile_over_b_log_geo$profile_axis, predict(lo4, profile_over_b_log_geo$profile_axis), col="black")
legend("topright", legend=c("V1", "V2", "V3", "V4"), col=c("red", "blue", "green", "black"), lty=1, cex=0.8)
abline(h=0.0576, col = "grey", lty=2, lwd=2)
abline(h=0.0365, col = "grey", lty=2, lwd=2)
abline(h=0.0322, col = "grey", lty=2, lwd=2)
abline(h=0.0304, col = "grey", lty=2, lwd=2)

plot(profile_over_b_log_geo$profile_axis, profile_over_b_log_geo$corresponding_pars[,1], type="l", col="red", xlab="Slope", ylab="Intercept", main = "V")
lines(profile_over_b_log_geo$profile_axis, profile_over_b_log_geo$corresponding_pars[,2], col="blue")
lines(profile_over_b_log_geo$profile_axis, profile_over_b_log_geo$corresponding_pars[,3], col="green")
lines(profile_over_b_log_geo$profile_axis, profile_over_b_log_geo$corresponding_pars[,4], col="black")
legend("topright", legend=c("V1", "V2", "V3", "V4"), col=c("red", "blue", "green", "black"), lty=1, cex=0.8)

lo1 <- loess(profile_over_b_sqrt_geo$eqd_dists[,1] ~ profile_over_b_sqrt_geo$profile_axis)
lo2 <- loess(profile_over_b_sqrt_geo$eqd_dists[,2] ~ profile_over_b_sqrt_geo$profile_axis)
lo3 <- loess(profile_over_b_sqrt_geo$eqd_dists[,3] ~ profile_over_b_sqrt_geo$profile_axis)
lo4 <- loess(profile_over_b_sqrt_geo$eqd_dists[,4] ~ profile_over_b_sqrt_geo$profile_axis)
plot(profile_over_b_sqrt_geo$profile_axis, predict(lo1, profile_over_b_sqrt_geo$profile_axis), type="l", col="red", xlab="Slope", ylab="EQD", main = "sqrtV")
lines(profile_over_b_sqrt_geo$profile_axis, predict(lo2, profile_over_b_sqrt_geo$profile_axis), col="blue")
lines(profile_over_b_sqrt_geo$profile_axis, predict(lo3, profile_over_b_sqrt_geo$profile_axis), col="green")
lines(profile_over_b_sqrt_geo$profile_axis, predict(lo4, profile_over_b_sqrt_geo$profile_axis), col="black")
legend("topright", legend=c("V1", "V2", "V3", "V4"), col=c("red", "blue", "green", "black"), lty=1, cex=0.8)
abline(h=0.0576, col = "grey", lty=2, lwd=2)
abline(h=0.0365, col = "grey", lty=2, lwd=2)
abline(h=0.0322, col = "grey", lty=2, lwd=2)
abline(h=0.0304, col = "grey", lty=2, lwd=2)

plot(profile_over_b_sqrt_geo$profile_axis, profile_over_b_sqrt_geo$corresponding_pars[,1], type="l", col="red", xlab="Slope", ylab="Intercept", main = "sqrtV")
lines(profile_over_b_sqrt_geo$profile_axis, profile_over_b_sqrt_geo$corresponding_pars[,2], col="blue")
lines(profile_over_b_sqrt_geo$profile_axis, profile_over_b_sqrt_geo$corresponding_pars[,3], col="green")
lines(profile_over_b_sqrt_geo$profile_axis, profile_over_b_sqrt_geo$corresponding_pars[,4], col="black")
legend("topright", legend=c("V1", "V2", "V3", "V4"), col=c("red", "blue", "green", "black"), lty=1, cex=0.8)


lo1 <- loess(log_geo_thresh_fit[[1]]$dists ~ threshold_matrix)

ggplot(data.frame(a=threshold_matrix[,1], b=threshold_matrix[,2], dists=predict(lo1, threshold_matrix)), aes(x=a, y=b, fill=dists)) + geom_tile() + scale_fill_viridis_c() + ggtitle("EQD V1") + theme_minimal()

threshold_matrix[which.min(predict(lo1, threshold_matrix)),]

# Plot Magnitudes against V1-V4
dev.new(height=20, width=20, noRStudioGD = TRUE)
par(mfrow=c(2,2), bg='transparent')
xlimits <- c(min(gron_eq_cat$V_1, gron_eq_cat$V_2, gron_eq_cat$V_3, gron_eq_cat$V_4), max(gron_eq_cat$V_1, gron_eq_cat$V_2, gron_eq_cat$V_3, gron_eq_cat$V_4))
plot(gron_eq_cat$V_1, gron_eq_cat$Magnitude, xlab = "V1", ylab = "Magnitude", main = "V1", col="grey", pch=19, xlim=xlimits)
plot(gron_eq_cat$V_2, gron_eq_cat$Magnitude, xlab = "V2", ylab = "Magnitude", main = "V2", col="grey", pch=19, xlim=xlimits)
plot(gron_eq_cat$V_3, gron_eq_cat$Magnitude, xlab = "V3", ylab = "Magnitude", main = "V3", col="grey", pch=19, xlim=xlimits)
plot(gron_eq_cat$V_4, gron_eq_cat$Magnitude, xlab = "V4", ylab = "Magnitude", main = "V4", col="grey", pch=19, xlim=xlimits)


# Fitting Inf thresolds
threshold_choice <- threshold_matrix[which(log_geo_thresh_fit[[1]]$dists == Inf),]

fitted_threshold <- threshold_choice[1] + threshold_choice[2]*log(third_nearest_dist_matrix[,1])

excesses <- mags[mags > fitted_threshold]
log_third_ex <- log(third_nearest_dist_matrix[,1])[mags > fitted_threshold]

sampled_indices <- sample(1:length(excesses), length(excesses), replace = TRUE)
excesses <- excesses[sampled_indices]
log_third_ex <- log_third_ex[sampled_indices]

gpd_fit <- optim(GPD_LL_given_third_nearest, par = c(mean(excesses), 0.1), 
                 excess = excesses, thresh_par = threshold_choice, third_nearest = log_third_ex, control=list(fnscale = -1), hessian = TRUE)

sigma_tilde <- gpd_fit$par[2] - gpd_fit$par[1]*(threshold_choice[1] + threshold_choice[2]*log_third_ex)
transformed_excesses <- transform_to_exp(excesses, sigma_tilde, gpd_fit$par[1])

probs <- 1:length(transformed_excesses)/(length(transformed_excesses)+1)
m2 <- mean(abs(quantile(transformed_excesses,probs)- qexp(probs)))


















