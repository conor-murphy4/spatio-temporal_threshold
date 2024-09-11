gron_eq_cat <- read.csv("Data/Events/unrounded_after_geophone_start_in_polygon_with_V_3d.csv", header=T)

#All commented code could possibly be remioved.. This is code related to 2D version of V


# Constant GPD fit with constant threshold estimated using EQD --------------------------------

source("src/eqd.R")
mags <- gron_eq_cat$Magnitude

thresholds <- quantile(mags, seq(0,0.95, by=0.01))
const_thresh_fit <- eqd_exp(mags, thresh = thresholds, k=200)

excesses_const <- mags[mags > const_thresh_fit$thresh] - const_thresh_fit$thresh

# Conservative constant threshold from Zak's work -------------------------
conservative_threshold <- 1.45

excesses_conservative <- mags[mags > conservative_threshold] - conservative_threshold

conservative_thresh_fit <- optim(GPD_LL, par=c(mean(excesses_conservative), 0.1), z = excesses_conservative, control = list(fnscale=-1))

# Piece-wise constant threshold according to index ------------------------

u_h_length <- 1089
piecewise_const_thresh <- c(rep(1.15, u_h_length), rep(0.76, length(mags) - u_h_length))

excesses_piecewise_const <- mags[mags > piecewise_const_thresh] - piecewise_const_thresh[mags > piecewise_const_thresh]

step_thresh_exc <- piecewise_const_thresh[mags > piecewise_const_thresh]

piecewise_const_thresh_fit <- optim(GPD_LL_step, par=c(mean(excesses_piecewise_const), 0.1), excess = excesses_piecewise_const, thresh = step_thresh_exc, control = list(fnscale=-1))

# Non-stationary GPD fit with threshold based on EQD given third_nearest_dist --------

source("src/eqd_geo.R")
# third_nearest_dist <- gron_eq_cat$V
third_nearest_dist_3d <- gron_eq_cat$V_3d

# intercepts <- seq(0, 1.5, by=0.01)
# slopes <- seq(0,2.5, by=0.01)
# threshold_matrix <- as.matrix(expand.grid(intercepts, slopes))
# 
# geo_thresh_fit <- eqd_geo(data=mags, thresh = threshold_matrix, third_nearest_distance = third_nearest_dist, k=200, min_dist = min(third_nearest_dist), max_dist = max(third_nearest_dist))
# 
# geo_chosen_threshold <- geo_thresh_fit$thresh[1] + geo_thresh_fit$thresh[2]*third_nearest_dist
# excesses_geo <- mags[mags > geo_chosen_threshold] - geo_chosen_threshold[mags > geo_chosen_threshold]
# geo_sigma <- geo_thresh_fit$par[1] + geo_thresh_fit$par[2]*geo_chosen_threshold[mags > geo_chosen_threshold]
# 
# # Non-stationary GPD fit with threshold based on EQD given log(third_nearest_dist)--------
# intercepts <- 1.4 #seq(0.3, 2, by=0.1)
# slopes <- seq(0,0.5, by=0.01)
# threshold_matrix <- as.matrix(expand.grid(intercepts, slopes))
# 
# log_third_nearest_dist <- log(third_nearest_dist)
# 
# log_geo_thresh_fit <- eqd_geo(data=mags, thresh = threshold_matrix, third_nearest_distance = log_third_nearest_dist, k=200, min_dist = min(log_third_nearest_dist), max_dist = max(log_third_nearest_dist))
# 
# log_geo_chosen_threshold <- log_geo_thresh_fit$thresh[1] + log_geo_thresh_fit$thresh[2]*log_third_nearest_dist
# excesses_log_geo <- mags[mags > log_geo_chosen_threshold] - log_geo_chosen_threshold[mags > log_geo_chosen_threshold]
# log_geo_sigma <- log_geo_thresh_fit$par[1] + log_geo_thresh_fit$par[2]*log_geo_chosen_threshold[mags > log_geo_chosen_threshold]
# 
# # Non-stationary GPD fit with threshold based on EQD given sqrt(third_nearest_dist)--------
# intercepts <- seq(0, 1, by=0.1)
# slopes <- seq(0.1,2, by=0.1)
# threshold_matrix <- as.matrix(expand.grid(intercepts, slopes))
# 
# sqrt_third_nearest_dist <- sqrt(third_nearest_dist)
# 
# sqrt_geo_thresh_fit <- eqd_geo(data=mags, thresh = threshold_matrix, third_nearest_distance = sqrt_third_nearest_dist, k=200, min_dist = min(sqrt_third_nearest_dist), max_dist = max(sqrt_third_nearest_dist))
# 
# sqrt_geo_chosen_threshold <- sqrt_geo_thresh_fit$thresh[1] + sqrt_geo_thresh_fit$thresh[2]*sqrt_third_nearest_dist
# excesses_sqrt_geo <- mags[mags > sqrt_geo_chosen_threshold] - sqrt_geo_chosen_threshold[mags > sqrt_geo_chosen_threshold]
# sqrt_geo_sigma <- sqrt_geo_thresh_fit$par[1] + sqrt_geo_thresh_fit$par[2]*sqrt_geo_chosen_threshold[mags > sqrt_geo_chosen_threshold]

#3d versions including depth in V function
intercepts <- seq(0, 1.5, by=0.01)
slopes <- seq(0,2.5, by=0.01)
threshold_matrix <- as.matrix(expand.grid(intercepts, slopes))
geo_thresh_fit_3d <- eqd_geo(data=mags, thresh = threshold_matrix, third_nearest_distance = third_nearest_dist_3d, k=200, min_dist = min(third_nearest_dist_3d), max_dist = max(third_nearest_dist_3d))

geo_chosen_threshold_3d <- geo_thresh_fit_3d$thresh[1] + geo_thresh_fit_3d$thresh[2]*third_nearest_dist_3d
excesses_geo_3d <- mags[mags > geo_chosen_threshold_3d] - geo_chosen_threshold_3d[mags > geo_chosen_threshold_3d]
geo_sigma_3d <- geo_thresh_fit_3d$par[1] + geo_thresh_fit_3d$par[2]*geo_chosen_threshold_3d[mags > geo_chosen_threshold_3d]

log_third_nearest_dist_3d <- log(third_nearest_dist_3d)

log_geo_thresh_fit_3d <- eqd_geo(data=mags, thresh = threshold_matrix, third_nearest_distance = log_third_nearest_dist_3d, k=200, min_dist = min(log_third_nearest_dist_3d), max_dist = max(log_third_nearest_dist_3d))

log_geo_chosen_threshold_3d <- log_geo_thresh_fit_3d$thresh[1] + log_geo_thresh_fit_3d$thresh[2]*log_third_nearest_dist_3d
excesses_log_geo_3d <- mags[mags > log_geo_chosen_threshold_3d] - log_geo_chosen_threshold_3d[mags > log_geo_chosen_threshold_3d]
log_geo_sigma_3d <- log_geo_thresh_fit_3d$par[1] + log_geo_thresh_fit_3d$par[2]*log_geo_chosen_threshold_3d[mags > log_geo_chosen_threshold_3d]

sqrt_third_nearest_dist_3d <- sqrt(third_nearest_dist_3d)

sqrt_geo_thresh_fit_3d <- eqd_geo(data=mags, thresh = threshold_matrix, third_nearest_distance = sqrt_third_nearest_dist_3d, k=200, min_dist = min(sqrt_third_nearest_dist_3d), max_dist = max(sqrt_third_nearest_dist_3d))

sqrt_geo_chosen_threshold_3d <- sqrt_geo_thresh_fit_3d$thresh[1] + sqrt_geo_thresh_fit_3d$thresh[2]*sqrt_third_nearest_dist_3d
excesses_sqrt_geo_3d <- mags[mags > sqrt_geo_chosen_threshold_3d] - sqrt_geo_chosen_threshold_3d[mags > sqrt_geo_chosen_threshold_3d]
sqrt_geo_sigma_3d <- sqrt_geo_thresh_fit_3d$par[1] + sqrt_geo_thresh_fit_3d$par[2]*sqrt_geo_chosen_threshold_3d[mags > sqrt_geo_chosen_threshold_3d]

# Visualising thresholds
dev.new(width=30, height=10,noRStudioGD = TRUE)
par(mfrow=c(1,3),bg='transparent')
plot(gron_eq_cat$Magnitude, main = "Magnitudes over index", xlab = "Index", ylab = "Magnitude", col="grey", pch=19)
lines(piecewise_const_thresh, col="red", lw=2)
lines(rep(conservative_threshold, length(mags)), col="blue", lw=2)
lines(rep(const_thresh_fit$thresh, length(mags)), col="green", lw=2)
legend("topright", legend=c("Piecewise constant", "Conservative", "EQD constant"), col=c("red", "blue", "green"), lty=1, cex=0.8)

plot(gron_eq_cat$Magnitude, main = "Magnitudes over index", xlab = "Index", ylab = "Magnitude", col="grey", pch=19)
lines(geo_thresh_fit$thresh[[1]] + geo_thresh_fit$thresh[[2]]*third_nearest_dist, col="purple")
lines(log_geo_thresh_fit$thresh[[1]] + log_geo_thresh_fit$thresh[[2]]*log_third_nearest_dist, col="orange")
lines(sqrt_geo_thresh_fit$thresh[[1]] + sqrt_geo_thresh_fit$thresh[[2]]*sqrt_third_nearest_dist, col="brown")
legend("topright", legend=c("Threshold on V", "Threshold on log(V)", "Threshold on sqrt(V)"), col=c("purple", "orange", "brown"), lty=1, cex=0.8)

plot(gron_eq_cat$Magnitude, main = "Magnitudes over index", xlab = "Index", ylab = "Magnitude", col="grey", pch=19)
lines(geo_thresh_fit_3d$thresh[[1]] + geo_thresh_fit_3d$thresh[[2]]*third_nearest_dist_3d, col="purple")
lines(log_geo_thresh_fit_3d$thresh[[1]] + log_geo_thresh_fit_3d$thresh[[2]]*log_third_nearest_dist_3d, col="orange")
lines(sqrt_geo_thresh_fit_3d$thresh[[1]] + sqrt_geo_thresh_fit_3d$thresh[[2]]*sqrt_third_nearest_dist_3d, col="brown")
legend("topright", legend=c("Threshold on V_3d", "Threshold on log(V_3d)", "Threshold on sqrt(V_3d)"), col=c("purple", "orange", "brown"), lty=1, cex=0.8)

# Comparing QQplots on standard Exponential margins -----------------------

#QQplot on Exp scale for easy comparison to more complex models

exp_transformed_excesses_const <- transform_to_exp(excesses_const, sig = const_thresh_fit$par[1], xi = const_thresh_fit$par[2])
exp_transformed_excesses_conservative <- transform_to_exp(excesses_conservative, sig = conservative_thresh_fit$par[1], xi = conservative_thresh_fit$par[2])
exp_transformed_excesses_piecewise_const <- transform_to_exp(excesses_piecewise_const, sig = piecewise_const_thresh_fit$par[1], xi = piecewise_const_thresh_fit$par[2])
# exp_transformed_excesses_geo <- transform_to_exp(excesses_geo, sig = geo_sigma, xi = geo_thresh_fit$par[2])
# exp_transformed_excesses_log_geo <- transform_to_exp(excesses_log_geo, sig = log_geo_sigma, xi = log_geo_thresh_fit$par[2])
# exp_transformed_excesses_sqrt_geo <- transform_to_exp(excesses_sqrt_geo, sig = sqrt_geo_sigma, xi = sqrt_geo_thresh_fit$par[2])

#3d excesses
exp_transformed_excesses_geo_3d <- transform_to_exp(excesses_geo_3d, sig = geo_sigma_3d, xi = geo_thresh_fit_3d$par[2])
exp_transformed_excesses_log_geo_3d <- transform_to_exp(excesses_log_geo_3d, sig = log_geo_sigma_3d, xi = log_geo_thresh_fit_3d$par[2])
exp_transformed_excesses_sqrt_geo_3d <- transform_to_exp(excesses_sqrt_geo_3d, sig = sqrt_geo_sigma_3d, xi = sqrt_geo_thresh_fit_3d$par[2])

sample_quantiles_exp_const <- quantile(exp_transformed_excesses_const, c(1:length(excesses_const))/(length(excesses_const)+1))
sample_quantiles_exp_conservative <- quantile(exp_transformed_excesses_conservative, c(1:length(excesses_conservative))/(length(excesses_conservative)+1))
sample_quantiles_exp_piecewise_const <- quantile(exp_transformed_excesses_piecewise_const, c(1:length(excesses_piecewise_const))/(length(excesses_piecewise_const)+1))
# sample_quantiles_exp_geo <- quantile(exp_transformed_excesses_geo, probs)
# sample_quantiles_exp_log_geo <- quantile(exp_transformed_excesses_log_geo, probs)
# sample_quantiles_exp_sqrt_geo <- quantile(exp_transformed_excesses_sqrt_geo, probs)

#3d sample quantiles
sample_quantiles_exp_geo_3d <- quantile(exp_transformed_excesses_geo_3d, c(1:length(excesses_geo_3d))/(length(excesses_geo_3d)+1))
sample_quantiles_exp_log_geo_3d <- quantile(exp_transformed_excesses_log_geo_3d, c(1:length(excesses_log_geo_3d))/(length(excesses_log_geo_3d)+1))
sample_quantiles_exp_sqrt_geo_3d <- quantile(exp_transformed_excesses_sqrt_geo_3d, probs)

model_quantiles_exp <- qexp(probs, rate = 1)

n_boot <- 200
bootstrapped_model_quants_const <- matrix(NA, nrow = n_boot, ncol= length(excesses_const))
bootstrapped_model_quants_conservative <- matrix(NA, nrow = n_boot, ncol= length(excesses_conservative))
bootstrapped_model_quants_piecewise_const <- matrix(NA, nrow = n_boot, ncol= length(excesses_piecewise_const))
# bootstrapped_model_quants_geo <- matrix(NA, nrow = n_boot, ncol= length(probs))
# bootstrapped_model_quants_log_geo <- matrix(NA, nrow = n_boot, ncol= length(probs))
# bootstrapped_model_quants_sqrt_geo <- matrix(NA, nrow = n_boot, ncol= length(probs))
#3d bootstrapped model quants
bootstrapped_model_quants_geo_3d <- matrix(NA, nrow = n_boot, ncol= length(excesses_geo_3d))
bootstrapped_model_quants_log_geo_3d <- matrix(NA, nrow = n_boot, ncol= length(excesses_log_geo_3d))
bootstrapped_model_quants_sqrt_geo_3d <- matrix(NA, nrow = n_boot, ncol= length(excesses_sqrt_geo_3d))

for(i in 1:n_boot){
  boot_excess_exp_const <- rexp(length(excesses_const), rate = 1)
  boot_excess_exp_conservative <- rexp(length(excesses_conservative), rate = 1)
  boot_excess_exp_piecewise_const <- rexp(length(excesses_piecewise_const), rate = 1)
  # boot_excess_exp_geo <- rexp(length(excesses_geo), rate = 1)
  # boot_excess_exp_log_geo <- rexp(length(excesses_log_geo), rate = 1)
  # boot_excess_exp_sqrt_geo <- rexp(length(excesses_sqrt_geo), rate = 1)
  
  #3d bootstrapped model quants
  boot_excess_exp_geo_3d <- rexp(length(excesses_geo_3d), rate = 1)
  boot_excess_exp_log_geo_3d <- rexp(length(excesses_log_geo_3d), rate = 1)
  boot_excess_exp_sqrt_geo_3d <- rexp(length(excesses_sqrt_geo_3d), rate = 1)
  
  #Do we not need to refit to incorporate the uncertainty in the model fit into the tolerance bounds?
  
  bootstrapped_model_quants_const[i,] <- quantile(boot_excess_exp_const, c(1:length(excesses_const))/(length(excesses_const)+1))
  bootstrapped_model_quants_conservative[i,] <- quantile(boot_excess_exp_conservative, c(1:length(excesses_conservative))/(length(excesses_conservative)+1))
  bootstrapped_model_quants_piecewise_const[i,] <- quantile(boot_excess_exp_piecewise_const, c(1:length(excesses_piecewise_const))/(length(excesses_piecewise_const)+1))
  # bootstrapped_model_quants_geo[i,] <- quantile(boot_excess_exp_geo, probs)
  # bootstrapped_model_quants_log_geo[i,] <- quantile(boot_excess_exp_log_geo, probs)
  # bootstrapped_model_quants_sqrt_geo[i,] <- quantile(boot_excess_exp_sqrt_geo, probs)
  
  #3d bootstrapped model quants
  bootstrapped_model_quants_geo_3d[i,] <- quantile(boot_excess_exp_geo_3d, c(1:length(excesses_geo_3d))/(length(excesses_geo_3d)+1))
  bootstrapped_model_quants_log_geo_3d[i,] <- quantile(boot_excess_exp_log_geo_3d, c(1:length(excesses_log_geo_3d))/(length(excesses_log_geo_3d)+1))
  bootstrapped_model_quants_sqrt_geo_3d[i,] <- quantile(boot_excess_exp_sqrt_geo_3d, c(1:length(excesses_sqrt_geo_3d))/(length(excesses_sqrt_geo_3d)+1))
}

upper_tolerance_int_const <- apply(bootstrapped_model_quants_const, 2, quantile, prob = 0.975)
upper_tolerance_int_conservative <- apply(bootstrapped_model_quants_conservative, 2, quantile, prob = 0.975)
upper_tolerance_int_piecewise_const <- apply(bootstrapped_model_quants_piecewise_const, 2, quantile, prob = 0.975)
# upper_tolerance_int_geo <- apply(bootstrapped_model_quants_geo, 2, quantile, prob = 0.975)
# upper_tolerance_int_log_geo <- apply(bootstrapped_model_quants_log_geo, 2, quantile, prob = 0.975)
# upper_tolerance_int_sqrt_geo <- apply(bootstrapped_model_quants_sqrt_geo, 2, quantile, prob = 0.975)

#3d upper tolerance intervals
upper_tolerance_int_geo_3d <- apply(bootstrapped_model_quants_geo_3d, 2, quantile, prob = 0.975)
upper_tolerance_int_log_geo_3d <- apply(bootstrapped_model_quants_log_geo_3d, 2, quantile, prob = 0.975)
upper_tolerance_int_sqrt_geo_3d <- apply(bootstrapped_model_quants_sqrt_geo_3d, 2, quantile, prob = 0.975)

lower_tolerance_int_const <- apply(bootstrapped_model_quants_const, 2, quantile, prob = 0.025)
lower_tolerance_int_conservative <- apply(bootstrapped_model_quants_conservative, 2, quantile, prob = 0.025)
lower_tolerance_int_piecewise_const <- apply(bootstrapped_model_quants_piecewise_const, 2, quantile, prob = 0.025)

# lower_tolerance_int_geo <- apply(bootstrapped_model_quants_geo, 2, quantile, prob = 0.025)
# lower_tolerance_int_log_geo <- apply(bootstrapped_model_quants_log_geo, 2, quantile, prob = 0.025)
# lower_tolerance_int_sqrt_geo <- apply(bootstrapped_model_quants_sqrt_geo, 2, quantile, prob = 0.025)

#3d lower tolerance intervals
lower_tolerance_int_geo_3d <- apply(bootstrapped_model_quants_geo_3d, 2, quantile, prob = 0.025)
lower_tolerance_int_log_geo_3d <- apply(bootstrapped_model_quants_log_geo_3d, 2, quantile, prob = 0.025)
lower_tolerance_int_sqrt_geo_3d <- apply(bootstrapped_model_quants_sqrt_geo_3d, 2, quantile, prob = 0.025)

dev.new(height=20, width=30, noRStudioGD = TRUE)
par(mfrow=c(3,2), bg='transparent')

max_plot <- max(c(upper_tolerance_int_const, upper_tolerance_int_geo, upper_tolerance_int_log_geo, upper_tolerance_int_sqrt_geo, upper_tolerance_int_geo_3d, upper_tolerance_int_log_geo_3d, upper_tolerance_int_sqrt_geo_3d))

model_quantiles_exp <- qexp(c(1:length(excesses_const))/(length(excesses_const)+1), rate = 1)
plot(model_quantiles_exp, sample_quantiles_exp_const, main = "EQD threshold", ylab ="Sample quantiles", xlab = "Model quantiles", pch=19, asp=1, ylim = c(0, max_plot))
lines(model_quantiles_exp, upper_tolerance_int_const, lty="dashed", lwd=2, col="red")
lines(model_quantiles_exp, lower_tolerance_int_const, lty="dashed", lwd=2, col="red")
abline(a=0, b=1, lwd=2, col="grey")

model_quantiles_exp <- qexp(c(1:length(excesses_conservative))/(length(excesses_conservative)+1), rate = 1)
plot(model_quantiles_exp, sample_quantiles_exp_conservative, main = "Conservative threshold", ylab ="Sample quantiles", xlab = "Model quantiles", pch=19, asp=1, ylim = c(0, max_plot))
lines(model_quantiles_exp, upper_tolerance_int_conservative, lty="dashed", lwd=2, col="red")
lines(model_quantiles_exp, lower_tolerance_int_conservative, lty="dashed", lwd=2, col="red")
abline(a=0, b=1, lwd=2, col="grey")

model_quantiles_exp <- qexp(c(1:length(excesses_piecewise_const))/(length(excesses_piecewise_const)+1), rate = 1)
plot(model_quantiles_exp, sample_quantiles_exp_piecewise_const, main = "Piecewise constant threshold", ylab ="Sample quantiles", xlab = "Model quantiles", pch=19, asp=1, ylim = c(0, max_plot))
lines(model_quantiles_exp, upper_tolerance_int_piecewise_const, lty="dashed", lwd=2, col="red")
lines(model_quantiles_exp, lower_tolerance_int_piecewise_const, lty="dashed", lwd=2, col="red")
abline(a=0, b=1, lwd=2, col="grey")

# plot(model_quantiles_exp, sample_quantiles_exp_geo, main = "EQD threshold with V", ylab ="Sample quantiles", xlab = "Model quantiles", pch=19, asp=1, ylim = c(0, max_plot))
# lines(model_quantiles_exp, upper_tolerance_int_geo, lty="dashed", lwd=2, col="red")
# lines(model_quantiles_exp, lower_tolerance_int_geo, lty="dashed", lwd=2, col="red")
# abline(a=0, b=1, lwd=2, col="grey")
# 
# plot(model_quantiles_exp, sample_quantiles_exp_log_geo, main = "EQD threshold with logV", ylab ="Sample quantiles", xlab = "Model quantiles", pch=19, asp=1, ylim = c(0, max_plot))
# lines(model_quantiles_exp, upper_tolerance_int_log_geo, lty="dashed", lwd=2, col="red")
# lines(model_quantiles_exp, lower_tolerance_int_log_geo, lty="dashed", lwd=2, col="red")
# abline(a=0, b=1, lwd=2, col="grey")
# 
# plot(model_quantiles_exp, sample_quantiles_exp_sqrt_geo, main = "EQD threshold with sqrtV", ylab ="Sample quantiles", xlab = "Model quantiles", pch=19, asp=1, ylim = c(0, max_plot))
# lines(model_quantiles_exp, upper_tolerance_int_sqrt_geo, lty="dashed", lwd=2, col="red")
# lines(model_quantiles_exp, lower_tolerance_int_sqrt_geo, lty="dashed", lwd=2, col="red")
# abline(a=0, b=1, lwd=2, col="grey")

model_quantiles_exp <- qexp(c(1:length(excesses_geo_3d))/(length(excesses_geo_3d)+1), rate = 1)
plot(model_quantiles_exp, sample_quantiles_exp_geo_3d, main = "EQD threshold with V 3D", ylab ="Sample quantiles", xlab = "Model quantiles", pch=19, asp=1, ylim = c(0, max_plot))
lines(model_quantiles_exp, upper_tolerance_int_geo_3d, lty="dashed", lwd=2, col="red")
lines(model_quantiles_exp, lower_tolerance_int_geo_3d, lty="dashed", lwd=2, col="red")
abline(a=0, b=1, lwd=2, col="grey")

model_quantiles_exp <- qexp(c(1:length(excesses_log_geo_3d))/(length(excesses_log_geo_3d)+1), rate = 1)
plot(model_quantiles_exp, sample_quantiles_exp_log_geo_3d, main = "EQD threshold with logV 3D", ylab ="Sample quantiles", xlab = "Model quantiles", pch=19, asp=1, ylim = c(0, max_plot))
lines(model_quantiles_exp, upper_tolerance_int_log_geo_3d, lty="dashed", lwd=2, col="red")
lines(model_quantiles_exp, lower_tolerance_int_log_geo_3d, lty="dashed", lwd=2, col="red")
abline(a=0, b=1, lwd=2, col="grey")

model_quantiles_exp <- qexp(c(1:length(excesses_sqrt_geo_3d))/(length(excesses_sqrt_geo_3d)+1), rate = 1)
plot(model_quantiles_exp, sample_quantiles_exp_sqrt_geo_3d, main = "EQD threshold with sqrtV 3D", ylab ="Sample quantiles", xlab = "Model quantiles", pch=19, asp=1, ylim = c(0, max_plot))
lines(model_quantiles_exp, upper_tolerance_int_sqrt_geo_3d, lty="dashed", lwd=2, col="red")
lines(model_quantiles_exp, lower_tolerance_int_sqrt_geo_3d, lty="dashed", lwd=2, col="red")
abline(a=0, b=1, lwd=2, col="grey")

#Plot of dists from each model

dev.new(height=20, width=20, noRStudioGD = TRUE)
par(mfrow=c(1,1), bg='transparent')

plot(geo_thresh_fit$dists, main = "EQD dists", ylab = "EQD value", xlab="Index", ylim=c(0, 0.7))
points(log_geo_thresh_fit$dists, col="orange")
points(sqrt_geo_thresh_fit$dists, col="red")
points(const_thresh_fit$dists, col="blue")
#3d dists
points(geo_thresh_fit_3d$dists, col="pink")
points(log_geo_thresh_fit_3d$dists, col="green")
points(sqrt_geo_thresh_fit_3d$dists, col="yellow")


legend("bottomright", legend=c("V", "logV", "sqrtV", "EQD"), col=c("black", "orange", "red", "blue"), pch=19, cex=1.5)


#Fits from each model

geo_thresh_fit_3d$thresh_par
geo_thresh_fit_3d$par
geo_thresh_fit_3d$num_excess

log_geo_thresh_fit_3d$thresh_par
log_geo_thresh_fit_3d$par
log_geo_thresh_fit_3d$num_excess

sqrt_geo_thresh_fit_3d$thresh_par
sqrt_geo_thresh_fit_3d$par
sqrt_geo_thresh_fit_3d$num_excess

const_thresh_fit$thresh
const_thresh_fit$par
const_thresh_fit$num_excess

conservative_threshold
conservative_thresh_fit$par
length(excesses_conservative)

c(1.15, 0.76)
piecewise_const_thresh_fit$par
length(excesses_piecewise_const)

# Visualising GPD fits over time

dates <- gron_eq_cat$Date[mags > geo_chosen_threshold_3d]
dates_log <- gron_eq_cat$Date[mags > log_geo_chosen_threshold_3d]
dates_sqrt <- gron_eq_cat$Date[mags > sqrt_geo_chosen_threshold_3d]

dev.new(height=10, width=20, noRStudioGD = TRUE)
par(mfrow=c(1,2), bg='transparent')

min_plot <- min(c(geo_sigma_3d, log_geo_sigma_3d, sqrt_geo_sigma_3d))
max_plot <- max(c(geo_sigma_3d, log_geo_sigma_3d, sqrt_geo_sigma_3d))
plot(as.Date(dates), geo_sigma_3d, main = "Sigma over time", xlab = "Date", ylab = "Sigma", col="purple", type='l', ylim=c(min_plot, max_plot))
lines(as.Date(dates_log), log_geo_sigma_3d, col="orange")
lines(as.Date(dates_sqrt), sqrt_geo_sigma_3d, col="red")
legend("bottomright", legend=c("V (xi=-0.243)", "logV (xi=-0.317)", "sqrtV (xi=-0.319)"), col=c("purple", "orange", "red"), lty=1, cex=1.2)

plot(which(mags > geo_chosen_threshold_3d), geo_sigma_3d, main = "Sigma over index", xlab = "Index", ylab = "Sigma", col="purple", type='l', ylim=c(min_plot, max_plot), xlim = c(0, length(mags)))
lines(which(mags > log_geo_chosen_threshold_3d), log_geo_sigma_3d, col="orange")
lines(which(mags > sqrt_geo_chosen_threshold_3d), sqrt_geo_sigma_3d, col="red")
legend("bottomright", legend=c("V (xi=-0.243)", "logV (xi=-0.317)", "sqrtV (xi=-0.319)"), col=c("purple", "orange", "red"), lty=1, cex=1.2)

