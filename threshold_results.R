source("src/helper_functions.R")

gron_eq_cat <- read.csv("Data/Events/unrounded_after_geophone_start_in_polygon_with_V_3d.csv", header=T)
mags <- gron_eq_cat$Magnitude
third_nearest_dist_3d <- gron_eq_cat$V_3d
log_third_nearest_dist_3d <- log(third_nearest_dist_3d)
sqrt_third_nearest_dist_3d <- sqrt(third_nearest_dist_3d)

#Fitted thresholds
const_thresh_fit <- readRDS("threshold_results/const_thresh_fit.rds")
geo_thresh_fit_3d <- readRDS("threshold_results/geo_thresh_fit_3d.rds")
log_geo_thresh_fit_3d <- readRDS("threshold_results/log_geo_thresh_fit_3d.rds")
sqrt_geo_thresh_fit_3d <- readRDS("threshold_results/sqrt_geo_thresh_fit_3d.rds")

# Conservative constant threshold from Zak's work -------------------------
conservative_threshold <- 1.45

excesses_conservative <- mags[mags > conservative_threshold] - conservative_threshold

conservative_thresh_fit <- optim(GPD_LL, par=c(mean(excesses_conservative), 0.1), z = excesses_conservative, control = list(fnscale=-1), hessian = TRUE)

(scale_ci_conserv <- conservative_thresh_fit$par[1] + c(-1.96,1.96)*sqrt(diag(solve(-conservative_thresh_fit$hessian)))[1])
(shape_ci_conserv <- conservative_thresh_fit$par[2] + c(-1.96,1.96)*sqrt(diag(solve(-conservative_thresh_fit$hessian)))[2])

# Piece-wise constant threshold according to index ------------------------

u_h_length <- 1089
piecewise_const_thresh <- c(rep(1.15, u_h_length), rep(0.76, length(mags) - u_h_length))

excesses_piecewise_const <- mags[mags > piecewise_const_thresh] - piecewise_const_thresh[mags > piecewise_const_thresh]

step_thresh_exc <- piecewise_const_thresh[mags > piecewise_const_thresh]

piecewise_const_thresh_fit <- optim(GPD_LL_step, par=c(mean(excesses_piecewise_const), 0.1), excess = excesses_piecewise_const, thresh = step_thresh_exc, control = list(fnscale=-1), hessian = T)

(scale_ci_pc <- piecewise_const_thresh_fit$par[1] + c(-1.96,1.96)*sqrt(diag(solve(-piecewise_const_thresh_fit$hessian)))[1])
(shape_ci_pc <- piecewise_const_thresh_fit$par[2] + c(-1.96,1.96)*sqrt(diag(solve(-piecewise_const_thresh_fit$hessian)))[2])

#Resulting thresholds, excesses and sigma(x,t) --------------------------------
excesses_const <- mags[mags > const_thresh_fit$thresh] - const_thresh_fit$thresh

const_opt_fit <- optim(GPD_LL, par=c(mean(excesses_const), 0.1), z = excesses_const, control = list(fnscale=-1), hessian = TRUE)

(scale_ci_const <- const_opt_fit$par[1] + c(-1.96,1.96)*sqrt(diag(solve(-const_opt_fit$hessian)))[1])
(shape_ci_const <- const_opt_fit$par[2] + c(-1.96,1.96)*sqrt(diag(solve(-const_opt_fit$hessian)))[2])

#Threshold based on V
geo_chosen_threshold_3d <- geo_thresh_fit_3d$thresh[1] + geo_thresh_fit_3d$thresh[2]*third_nearest_dist_3d
third_nearest_dist_excess <- third_nearest_dist_3d[mags > geo_chosen_threshold_3d]
excesses_geo_3d <- mags[mags > geo_chosen_threshold_3d] - geo_chosen_threshold_3d[mags > geo_chosen_threshold_3d]
geo_sigma_3d <- geo_thresh_fit_3d$par[1] + geo_thresh_fit_3d$par[2]*geo_chosen_threshold_3d[mags > geo_chosen_threshold_3d]

geo_thresh_opt_fit <- optim(GPD_LL_given_third_nearest, thresh_par=geo_thresh_fit_3d$thresh_par , par=c(mean(excesses_geo_3d), 0.1), excess = excesses_geo_3d, third_nearest = third_nearest_dist_excess, control = list(fnscale=-1), hessian = T, min_dist=min(third_nearest_dist_3d), max_dist=max(third_nearest_dist_3d))

(scale_ci_geo <- geo_thresh_opt_fit$par[1] + c(-1.96,1.96)*sqrt(diag(solve(-geo_thresh_opt_fit$hessian)))[1])
(shape_ci_geo <- geo_thresh_opt_fit$par[2] + c(-1.96,1.96)*sqrt(diag(solve(-geo_thresh_opt_fit$hessian)))[2])

log_geo_chosen_threshold_3d <- log_geo_thresh_fit_3d$thresh_par[1] + log_geo_thresh_fit_3d$thresh_par[2]*log_third_nearest_dist_3d
log_third_nearest_dist_excess <- log_third_nearest_dist_3d[mags > log_geo_chosen_threshold_3d]
excesses_log_geo_3d <- mags[mags > log_geo_chosen_threshold_3d] - log_geo_chosen_threshold_3d[mags > log_geo_chosen_threshold_3d]
log_geo_sigma_3d <- log_geo_thresh_fit_3d$par[1] + log_geo_thresh_fit_3d$par[2]*log_geo_chosen_threshold_3d[mags > log_geo_chosen_threshold_3d]

log_geo_thresh_opt_fit <- optim(GPD_LL_given_third_nearest, thresh_par=log_geo_thresh_fit_3d$thresh_par , par=c(mean(excesses_log_geo_3d), 0.1), excess = excesses_log_geo_3d, third_nearest = log_third_nearest_dist_excess, control = list(fnscale=-1), hessian = T, min_dist=min(log_third_nearest_dist_3d), max_dist=max(log_third_nearest_dist_3d))

(scale_ci_log_geo <- log_geo_thresh_opt_fit$par[1] + c(-1.96,1.96)*sqrt(diag(solve(-log_geo_thresh_opt_fit$hessian)))[1])
(shape_ci_log_geo <- log_geo_thresh_opt_fit$par[2] + c(-1.96,1.96)*sqrt(diag(solve(-log_geo_thresh_opt_fit$hessian)))[2])

sqrt_geo_chosen_threshold_3d <- sqrt_geo_thresh_fit_3d$thresh[1] + sqrt_geo_thresh_fit_3d$thresh[2]*sqrt_third_nearest_dist_3d
sqrt_third_nearest_dist_excess <- sqrt_third_nearest_dist_3d[mags > sqrt_geo_chosen_threshold_3d]
excesses_sqrt_geo_3d <- mags[mags > sqrt_geo_chosen_threshold_3d] - sqrt_geo_chosen_threshold_3d[mags > sqrt_geo_chosen_threshold_3d]
sqrt_geo_sigma_3d <- sqrt_geo_thresh_fit_3d$par[1] + sqrt_geo_thresh_fit_3d$par[2]*sqrt_geo_chosen_threshold_3d[mags > sqrt_geo_chosen_threshold_3d]

sqrt_geo_thresh_opt_fit <- optim(GPD_LL_given_third_nearest, thresh_par=sqrt_geo_thresh_fit_3d$thresh_par , par=c(mean(excesses_sqrt_geo_3d), 0.1), excess = excesses_sqrt_geo_3d, third_nearest = sqrt_third_nearest_dist_excess, control = list(fnscale=-1), hessian = T, min_dist=min(sqrt_third_nearest_dist_3d), max_dist=max(sqrt_third_nearest_dist_3d))

(scale_ci_sqrt_geo <- sqrt_geo_thresh_opt_fit$par[1] + c(-1.96,1.96)*sqrt(diag(solve(-sqrt_geo_thresh_opt_fit$hessian)))[1])
(shape_ci_sqrt_geo <- sqrt_geo_thresh_opt_fit$par[2] + c(-1.96,1.96)*sqrt(diag(solve(-sqrt_geo_thresh_opt_fit$hessian)))[2])


# Visualising thresholds

dev.new(width=20, height=10,noRStudioGD = TRUE)
par(mfrow=c(1,2),bg='transparent')

plot(gron_eq_cat$Magnitude, main = "Magnitudes over index", xlab = "Index", ylab = "Magnitude", col="grey", pch=19)
lines(piecewise_const_thresh, col="red", lw=2)
lines(rep(conservative_threshold, length(mags)), col="blue", lw=2)
lines(rep(const_thresh_fit$thresh, length(mags)), col="green", lw=2)
legend("topright", legend=c("Piecewise constant", "Conservative", "EQD constant"), col=c("red", "blue", "green"), lty=1, cex=0.8)

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
exp_transformed_excesses_geo_3d <- transform_to_exp(excesses_geo_3d, sig = geo_sigma_3d, xi = geo_thresh_fit_3d$par[2])
exp_transformed_excesses_log_geo_3d <- transform_to_exp(excesses_log_geo_3d, sig = log_geo_sigma_3d, xi = log_geo_thresh_fit_3d$par[2])
exp_transformed_excesses_sqrt_geo_3d <- transform_to_exp(excesses_sqrt_geo_3d, sig = sqrt_geo_sigma_3d, xi = sqrt_geo_thresh_fit_3d$par[2])

sample_quantiles_exp_const <- quantile(exp_transformed_excesses_const, c(1:length(excesses_const))/(length(excesses_const)+1))
sample_quantiles_exp_conservative <- quantile(exp_transformed_excesses_conservative, c(1:length(excesses_conservative))/(length(excesses_conservative)+1))
sample_quantiles_exp_piecewise_const <- quantile(exp_transformed_excesses_piecewise_const, c(1:length(excesses_piecewise_const))/(length(excesses_piecewise_const)+1))
sample_quantiles_exp_geo_3d <- quantile(exp_transformed_excesses_geo_3d, c(1:length(excesses_geo_3d))/(length(excesses_geo_3d)+1))
sample_quantiles_exp_log_geo_3d <- quantile(exp_transformed_excesses_log_geo_3d, c(1:length(excesses_log_geo_3d))/(length(excesses_log_geo_3d)+1))
sample_quantiles_exp_sqrt_geo_3d <- quantile(exp_transformed_excesses_sqrt_geo_3d, c(1:length(excesses_sqrt_geo_3d))/(length(excesses_sqrt_geo_3d)+1))

n_boot <- 200
bootstrapped_model_quants_const <- matrix(NA, nrow = n_boot, ncol= length(excesses_const))
bootstrapped_model_quants_conservative <- matrix(NA, nrow = n_boot, ncol= length(excesses_conservative))
bootstrapped_model_quants_piecewise_const <- matrix(NA, nrow = n_boot, ncol= length(excesses_piecewise_const))
bootstrapped_model_quants_geo_3d <- matrix(NA, nrow = n_boot, ncol= length(excesses_geo_3d))
bootstrapped_model_quants_log_geo_3d <- matrix(NA, nrow = n_boot, ncol= length(excesses_log_geo_3d))
bootstrapped_model_quants_sqrt_geo_3d <- matrix(NA, nrow = n_boot, ncol= length(excesses_sqrt_geo_3d))

for(i in 1:n_boot){
  boot_excess_exp_const <- rexp(length(excesses_const), rate = 1)
  boot_excess_exp_conservative <- rexp(length(excesses_conservative), rate = 1)
  boot_excess_exp_piecewise_const <- rexp(length(excesses_piecewise_const), rate = 1)
  boot_excess_exp_geo_3d <- rexp(length(excesses_geo_3d), rate = 1)
  boot_excess_exp_log_geo_3d <- rexp(length(excesses_log_geo_3d), rate = 1)
  boot_excess_exp_sqrt_geo_3d <- rexp(length(excesses_sqrt_geo_3d), rate = 1)
  
  bootstrapped_model_quants_const[i,] <- quantile(boot_excess_exp_const, c(1:length(excesses_const))/(length(excesses_const)+1))
  bootstrapped_model_quants_conservative[i,] <- quantile(boot_excess_exp_conservative, c(1:length(excesses_conservative))/(length(excesses_conservative)+1))
  bootstrapped_model_quants_piecewise_const[i,] <- quantile(boot_excess_exp_piecewise_const, c(1:length(excesses_piecewise_const))/(length(excesses_piecewise_const)+1))
  bootstrapped_model_quants_geo_3d[i,] <- quantile(boot_excess_exp_geo_3d, c(1:length(excesses_geo_3d))/(length(excesses_geo_3d)+1))
  bootstrapped_model_quants_log_geo_3d[i,] <- quantile(boot_excess_exp_log_geo_3d, c(1:length(excesses_log_geo_3d))/(length(excesses_log_geo_3d)+1))
  bootstrapped_model_quants_sqrt_geo_3d[i,] <- quantile(boot_excess_exp_sqrt_geo_3d, c(1:length(excesses_sqrt_geo_3d))/(length(excesses_sqrt_geo_3d)+1))
}

upper_tolerance_int_const <- apply(bootstrapped_model_quants_const, 2, quantile, prob = 0.975)
upper_tolerance_int_conservative <- apply(bootstrapped_model_quants_conservative, 2, quantile, prob = 0.975)
upper_tolerance_int_piecewise_const <- apply(bootstrapped_model_quants_piecewise_const, 2, quantile, prob = 0.975)
upper_tolerance_int_geo_3d <- apply(bootstrapped_model_quants_geo_3d, 2, quantile, prob = 0.975)
upper_tolerance_int_log_geo_3d <- apply(bootstrapped_model_quants_log_geo_3d, 2, quantile, prob = 0.975)
upper_tolerance_int_sqrt_geo_3d <- apply(bootstrapped_model_quants_sqrt_geo_3d, 2, quantile, prob = 0.975)

lower_tolerance_int_const <- apply(bootstrapped_model_quants_const, 2, quantile, prob = 0.025)
lower_tolerance_int_conservative <- apply(bootstrapped_model_quants_conservative, 2, quantile, prob = 0.025)
lower_tolerance_int_piecewise_const <- apply(bootstrapped_model_quants_piecewise_const, 2, quantile, prob = 0.025)
lower_tolerance_int_geo_3d <- apply(bootstrapped_model_quants_geo_3d, 2, quantile, prob = 0.025)
lower_tolerance_int_log_geo_3d <- apply(bootstrapped_model_quants_log_geo_3d, 2, quantile, prob = 0.025)
lower_tolerance_int_sqrt_geo_3d <- apply(bootstrapped_model_quants_sqrt_geo_3d, 2, quantile, prob = 0.025)

dev.new(height=20, width=30, noRStudioGD = TRUE)
par(mfrow=c(2,3), bg='transparent')

max_plot <- max(c(upper_tolerance_int_const, upper_tolerance_int_geo_3d, upper_tolerance_int_log_geo_3d, upper_tolerance_int_sqrt_geo_3d))

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

# Focusing on before/after changepoint and redoing QQplots

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

pc_thresh_before <- piecewise_const_thresh[1:u_h_length]
pc_thresh_after <- piecewise_const_thresh[(u_h_length+1):nrow(gron_eq_cat)]

geo_chosen_threshold_3d_before <- geo_thresh_fit_3d$thresh[1] + geo_thresh_fit_3d$thresh[2]*third_nearest_dist_3d_before
geo_chosen_threshold_3d_after <- geo_thresh_fit_3d$thresh[1] + geo_thresh_fit_3d$thresh[2]*third_nearest_dist_3d_after

log_geo_chosen_threshold_3d_before <- log_geo_thresh_fit_3d$thresh_par[1] + log_geo_thresh_fit_3d$thresh_par[2]*log_third_nearest_dist_3d_before
log_geo_chosen_threshold_3d_after <- log_geo_thresh_fit_3d$thresh_par[1] + log_geo_thresh_fit_3d$thresh_par[2]*log_third_nearest_dist_3d_after

sqrt_geo_chosen_threshold_3d_before <- sqrt_geo_thresh_fit_3d$thresh[1] + sqrt_geo_thresh_fit_3d$thresh[2]*sqrt_third_nearest_dist_3d_before
sqrt_geo_chosen_threshold_3d_after <- sqrt_geo_thresh_fit_3d$thresh[1] + sqrt_geo_thresh_fit_3d$thresh[2]*sqrt_third_nearest_dist_3d_after

geo_sigma_3d_before <- geo_thresh_fit_3d$par[1] + geo_thresh_fit_3d$par[2]*geo_chosen_threshold_3d_before[mags_before > geo_chosen_threshold_3d_before]
geo_sigma_3d_after <- geo_thresh_fit_3d$par[1] + geo_thresh_fit_3d$par[2]*geo_chosen_threshold_3d_after[mags_after > geo_chosen_threshold_3d_after]

log_geo_sigma_3d_before <- log_geo_thresh_fit_3d$par[1] + log_geo_thresh_fit_3d$par[2]*log_geo_chosen_threshold_3d_before[mags_before > log_geo_chosen_threshold_3d_before]
log_geo_sigma_3d_after <- log_geo_thresh_fit_3d$par[1] + log_geo_thresh_fit_3d$par[2]*log_geo_chosen_threshold_3d_after[mags_after > log_geo_chosen_threshold_3d_after]

sqrt_geo_sigma_3d_before <- sqrt_geo_thresh_fit_3d$par[1] + sqrt_geo_thresh_fit_3d$par[2]*sqrt_geo_chosen_threshold_3d_before[mags_before > sqrt_geo_chosen_threshold_3d_before]
sqrt_geo_sigma_3d_after <- sqrt_geo_thresh_fit_3d$par[1] + sqrt_geo_thresh_fit_3d$par[2]*sqrt_geo_chosen_threshold_3d_after[mags_after > sqrt_geo_chosen_threshold_3d_after]

#Excesses
excesses_conservative_before <- mags_before[mags_before > conservative_threshold] - conservative_threshold
excesses_conservative_after <- mags_after[mags_after > conservative_threshold] - conservative_threshold

excesses_pc_before <- mags_before[mags_before > pc_thresh_before] - pc_thresh_before[mags_before > pc_thresh_before]
excesses_pc_after <- mags_after[mags_after > pc_thresh_after] - pc_thresh_after[mags_after > pc_thresh_after]

excesses_const_before <- mags_before[mags_before > const_thresh_fit$thresh] - const_thresh_fit$thresh
excesses_const_after <- mags_after[mags_after > const_thresh_fit$thresh] - const_thresh_fit$thresh

excesses_geo_3d_before <- mags_before[mags_before > geo_chosen_threshold_3d_before] - geo_chosen_threshold_3d_before[mags_before > geo_chosen_threshold_3d_before]
excesses_geo_3d_after <- mags_after[mags_after > geo_chosen_threshold_3d_after] - geo_chosen_threshold_3d_after[mags_after > geo_chosen_threshold_3d_after]

excesses_log_geo_3d_before <- mags_before[mags_before > log_geo_chosen_threshold_3d_before] - log_geo_chosen_threshold_3d_before[mags_before > log_geo_chosen_threshold_3d_before]
excesses_log_geo_3d_after <- mags_after[mags_after > log_geo_chosen_threshold_3d_after] - log_geo_chosen_threshold_3d_after[mags_after > log_geo_chosen_threshold_3d_after]

excesses_sqrt_geo_3d_before <- mags_before[mags_before > sqrt_geo_chosen_threshold_3d_before] - sqrt_geo_chosen_threshold_3d_before[mags_before > sqrt_geo_chosen_threshold_3d_before]
excesses_sqrt_geo_3d_after <- mags_after[mags_after > sqrt_geo_chosen_threshold_3d_after] - sqrt_geo_chosen_threshold_3d_after[mags_after > sqrt_geo_chosen_threshold_3d_after]

#QQplots
exp_transformed_excesses_const_before <- transform_to_exp(excesses_const_before, sig = const_thresh_fit$par[1], xi = const_thresh_fit$par[2])
exp_transformed_excesses_conservative_before <- transform_to_exp(excesses_conservative_before, sig = conservative_thresh_fit$par[1], xi = conservative_thresh_fit$par[2])
exp_transformed_excesses_pc_before <- transform_to_exp(excesses_pc_before, sig = piecewise_const_thresh_fit$par[1], xi = piecewise_const_thresh_fit$par[2])
exp_transformed_excesses_geo_3d_before <- transform_to_exp(excesses_geo_3d_before, sig = geo_sigma_3d_before, xi = geo_thresh_fit_3d$par[2])
exp_transformed_excesses_log_geo_3d_before <- transform_to_exp(excesses_log_geo_3d_before, sig = log_geo_sigma_3d_before, xi = log_geo_thresh_fit_3d$par[2])
exp_transformed_excesses_sqrt_geo_3d_before <- transform_to_exp(excesses_sqrt_geo_3d_before, sig = sqrt_geo_sigma_3d_before, xi = sqrt_geo_thresh_fit_3d$par[2])

exp_transformed_excesses_const_after <- transform_to_exp(excesses_const_after, sig = const_thresh_fit$par[1], xi = const_thresh_fit$par[2])
exp_transformed_excesses_conservative_after <- transform_to_exp(excesses_conservative_after, sig = conservative_thresh_fit$par[1], xi = conservative_thresh_fit$par[2])
exp_transformed_excesses_pc_after <- transform_to_exp(excesses_pc_after, sig = piecewise_const_thresh_fit$par[1], xi = piecewise_const_thresh_fit$par[2])
exp_transformed_excesses_geo_3d_after <- transform_to_exp(excesses_geo_3d_after, sig = geo_sigma_3d_after, xi = geo_thresh_fit_3d$par[2])
exp_transformed_excesses_log_geo_3d_after <- transform_to_exp(excesses_log_geo_3d_after, sig = log_geo_sigma_3d_after, xi = log_geo_thresh_fit_3d$par[2])
exp_transformed_excesses_sqrt_geo_3d_after <- transform_to_exp(excesses_sqrt_geo_3d_after, sig = sqrt_geo_sigma_3d_after, xi = sqrt_geo_thresh_fit_3d$par[2])

sample_quantiles_exp_const_before <- quantile(exp_transformed_excesses_const_before, c(1:length(excesses_const_before))/(length(excesses_const_before)+1))
sample_quantiles_exp_conservative_before <- quantile(exp_transformed_excesses_conservative_before, c(1:length(excesses_conservative_before))/(length(excesses_conservative_before)+1))
sample_quantiles_exp_pc_before <- quantile(exp_transformed_excesses_pc_before, c(1:length(excesses_pc_before))/(length(excesses_pc_before)+1))
sample_quantiles_exp_geo_3d_before <- quantile(exp_transformed_excesses_geo_3d_before, c(1:length(excesses_geo_3d_before))/(length(excesses_geo_3d_before)+1))
sample_quantiles_exp_log_geo_3d_before <- quantile(exp_transformed_excesses_log_geo_3d_before, c(1:length(excesses_log_geo_3d_before))/(length(excesses_log_geo_3d_before)+1))
sample_quantiles_exp_sqrt_geo_3d_before <- quantile(exp_transformed_excesses_sqrt_geo_3d_before, c(1:length(excesses_sqrt_geo_3d_before))/(length(excesses_sqrt_geo_3d_before)+1))

sample_quantiles_exp_const_after <- quantile(exp_transformed_excesses_const_after, c(1:length(excesses_const_after))/(length(excesses_const_after)+1))
sample_quantiles_exp_conservative_after <- quantile(exp_transformed_excesses_conservative_after, c(1:length(excesses_conservative_after))/(length(excesses_conservative_after)+1))
sample_quantiles_exp_pc_after <- quantile(exp_transformed_excesses_pc_after, c(1:length(excesses_pc_after))/(length(excesses_pc_after)+1))
sample_quantiles_exp_geo_3d_after <- quantile(exp_transformed_excesses_geo_3d_after, c(1:length(excesses_geo_3d_after))/(length(excesses_geo_3d_after)+1))
sample_quantiles_exp_log_geo_3d_after <- quantile(exp_transformed_excesses_log_geo_3d_after, c(1:length(excesses_log_geo_3d_after))/(length(excesses_log_geo_3d_after)+1))
sample_quantiles_exp_sqrt_geo_3d_after <- quantile(exp_transformed_excesses_sqrt_geo_3d_after, c(1:length(excesses_sqrt_geo_3d_after))/(length(excesses_sqrt_geo_3d_after)+1))

n_boot <- 200
bootstrapped_model_quants_const_before <- matrix(NA, nrow = n_boot, ncol= length(excesses_const_before))
bootstrapped_model_quants_conservative_before <- matrix(NA, nrow = n_boot, ncol= length(excesses_conservative_before))
bootstrapped_model_quants_pc_before <- matrix(NA, nrow = n_boot, ncol= length(excesses_pc_before))
bootstrapped_model_quants_geo_3d_before <- matrix(NA, nrow = n_boot, ncol= length(excesses_geo_3d_before))
bootstrapped_model_quants_log_geo_3d_before <- matrix(NA, nrow = n_boot, ncol= length(excesses_log_geo_3d_before))
bootstrapped_model_quants_sqrt_geo_3d_before <- matrix(NA, nrow = n_boot, ncol= length(excesses_sqrt_geo_3d_before))

bootstrapped_model_quants_const_after <- matrix(NA, nrow = n_boot, ncol= length(excesses_const_after))
bootstrapped_model_quants_conservative_after <- matrix(NA, nrow = n_boot, ncol= length(excesses_conservative_after))
bootstrapped_model_quants_pc_after <- matrix(NA, nrow = n_boot, ncol= length(excesses_pc_after))
bootstrapped_model_quants_geo_3d_after <- matrix(NA, nrow = n_boot, ncol= length(excesses_geo_3d_after))
bootstrapped_model_quants_log_geo_3d_after <- matrix(NA, nrow = n_boot, ncol= length(excesses_log_geo_3d_after))
bootstrapped_model_quants_sqrt_geo_3d_after <- matrix(NA, nrow = n_boot, ncol= length(excesses_sqrt_geo_3d_after))

for(i in 1:n_boot){
  boot_excess_exp_const_before <- rexp(length(excesses_const_before), rate = 1)
  boot_excess_exp_conservative_before <- rexp(length(excesses_conservative_before), rate = 1)
  boot_excess_exp_pc_before <- rexp(length(excesses_pc_before), rate = 1)
  boot_excess_exp_geo_3d_before <- rexp(length(excesses_geo_3d_before), rate = 1)
  boot_excess_exp_log_geo_3d_before <- rexp(length(excesses_log_geo_3d_before), rate = 1)
  boot_excess_exp_sqrt_geo_3d_before <- rexp(length(excesses_sqrt_geo_3d_before), rate = 1)
  
  boot_excess_exp_const_after <- rexp(length(excesses_const_after), rate = 1)
  boot_excess_exp_conservative_after <- rexp(length(excesses_conservative_after), rate = 1)
  boot_excess_exp_pc_after <- rexp(length(excesses_pc_after), rate = 1)
  boot_excess_exp_geo_3d_after <- rexp(length(excesses_geo_3d_after), rate = 1)
  boot_excess_exp_log_geo_3d_after <- rexp(length(excesses_log_geo_3d_after), rate = 1)
  boot_excess_exp_sqrt_geo_3d_after <- rexp(length(excesses_sqrt_geo_3d_after), rate = 1)
  
  bootstrapped_model_quants_const_before[i,] <- quantile(boot_excess_exp_const_before, c(1:length(excesses_const_before))/(length(excesses_const_before)+1))
  bootstrapped_model_quants_conservative_before[i,] <- quantile(boot_excess_exp_conservative_before, c(1:length(excesses_conservative_before))/(length(excesses_conservative_before)+1))
  bootstrapped_model_quants_pc_before[i,] <- quantile(boot_excess_exp_pc_before, c(1:length(excesses_pc_before))/(length(excesses_pc_before)+1))
  bootstrapped_model_quants_geo_3d_before[i,] <- quantile(boot_excess_exp_geo_3d_before, c(1:length(excesses_geo_3d_before))/(length(excesses_geo_3d_before)+1))
  bootstrapped_model_quants_log_geo_3d_before[i,] <- quantile(boot_excess_exp_log_geo_3d_before, c(1:length(excesses_log_geo_3d_before))/(length(excesses_log_geo_3d_before)+1))
  bootstrapped_model_quants_sqrt_geo_3d_before[i,] <- quantile(boot_excess_exp_sqrt_geo_3d_before, c(1:length(excesses_sqrt_geo_3d_before))/(length(excesses_sqrt_geo_3d_before)+1))
  
  bootstrapped_model_quants_const_after[i,] <- quantile(boot_excess_exp_const_after, c(1:length(excesses_const_after))/(length(excesses_const_after)+1))
  bootstrapped_model_quants_conservative_after[i,] <- quantile(boot_excess_exp_conservative_after, c(1:length(excesses_conservative_after))/(length(excesses_conservative_after)+1))
  bootstrapped_model_quants_pc_after[i,] <- quantile(boot_excess_exp_pc_after, c(1:length(excesses_pc_after))/(length(excesses_pc_after)+1))
  bootstrapped_model_quants_geo_3d_after[i,] <- quantile(boot_excess_exp_geo_3d_after, c(1:length(excesses_geo_3d_after))/(length(excesses_geo_3d_after)+1))
  bootstrapped_model_quants_log_geo_3d_after[i,] <- quantile(boot_excess_exp_log_geo_3d_after, c(1:length(excesses_log_geo_3d_after))/(length(excesses_log_geo_3d_after)+1))
  bootstrapped_model_quants_sqrt_geo_3d_after[i,] <- quantile(boot_excess_exp_sqrt_geo_3d_after, c(1:length(excesses_sqrt_geo_3d_after))/(length(excesses_sqrt_geo_3d_after)+1))
}

upper_tolerance_int_const_before <- apply(bootstrapped_model_quants_const_before, 2, quantile, prob = 0.975)
upper_tolerance_int_conservative_before <- apply(bootstrapped_model_quants_conservative_before, 2, quantile, prob = 0.975)
upper_tolerance_int_pc_before <- apply(bootstrapped_model_quants_pc_before, 2, quantile, prob = 0.975)
upper_tolerance_int_geo_3d_before <- apply(bootstrapped_model_quants_geo_3d_before, 2, quantile, prob = 0.975)
upper_tolerance_int_log_geo_3d_before <- apply(bootstrapped_model_quants_log_geo_3d_before, 2, quantile, prob = 0.975)
upper_tolerance_int_sqrt_geo_3d_before <- apply(bootstrapped_model_quants_sqrt_geo_3d_before, 2, quantile, prob = 0.975)

lower_tolerance_int_const_before <- apply(bootstrapped_model_quants_const_before, 2, quantile, prob = 0.025)
lower_tolerance_int_conservative_before <- apply(bootstrapped_model_quants_conservative_before, 2, quantile, prob = 0.025)
lower_tolerance_int_pc_before <- apply(bootstrapped_model_quants_pc_before, 2, quantile, prob = 0.025)
lower_tolerance_int_geo_3d_before <- apply(bootstrapped_model_quants_geo_3d_before, 2, quantile, prob = 0.025)
lower_tolerance_int_log_geo_3d_before <- apply(bootstrapped_model_quants_log_geo_3d_before, 2, quantile, prob = 0.025)
lower_tolerance_int_sqrt_geo_3d_before <- apply(bootstrapped_model_quants_sqrt_geo_3d_before, 2, quantile, prob = 0.025)

upper_tolerance_int_const_after <- apply(bootstrapped_model_quants_const_after, 2, quantile, prob = 0.975)
upper_tolerance_int_conservative_after <- apply(bootstrapped_model_quants_conservative_after, 2, quantile, prob = 0.975)
upper_tolerance_int_pc_after <- apply(bootstrapped_model_quants_pc_after, 2, quantile, prob = 0.975)
upper_tolerance_int_geo_3d_after <- apply(bootstrapped_model_quants_geo_3d_after, 2, quantile, prob = 0.975)
upper_tolerance_int_log_geo_3d_after <- apply(bootstrapped_model_quants_log_geo_3d_after, 2, quantile, prob = 0.975)
upper_tolerance_int_sqrt_geo_3d_after <- apply(bootstrapped_model_quants_sqrt_geo_3d_after, 2, quantile, prob = 0.975)

lower_tolerance_int_const_after <- apply(bootstrapped_model_quants_const_after, 2, quantile, prob = 0.025)
lower_tolerance_int_conservative_after <- apply(bootstrapped_model_quants_conservative_after, 2, quantile, prob = 0.025)
lower_tolerance_int_pc_after <- apply(bootstrapped_model_quants_pc_after, 2, quantile, prob = 0.025)
lower_tolerance_int_geo_3d_after <- apply(bootstrapped_model_quants_geo_3d_after, 2, quantile, prob = 0.025)
lower_tolerance_int_log_geo_3d_after <- apply(bootstrapped_model_quants_log_geo_3d_after, 2, quantile, prob = 0.025)
lower_tolerance_int_sqrt_geo_3d_after <- apply(bootstrapped_model_quants_sqrt_geo_3d_after, 2, quantile, prob = 0.025)

dev.new(height=20, width=30, noRStudioGD = TRUE)
par(mfrow=c(2,3), bg='transparent')

max_plot <- max(c(upper_tolerance_int_const_before, upper_tolerance_int_geo_3d_before, upper_tolerance_int_log_geo_3d_before, upper_tolerance_int_sqrt_geo_3d_before))

model_quantiles_exp <- qexp(c(1:length(excesses_const_before))/(length(excesses_const_before)+1), rate = 1)
plot(model_quantiles_exp, sample_quantiles_exp_const_before, main = "EQD threshold before changepoint", ylab ="Sample quantiles", xlab = "Model quantiles", pch=19, asp=1, ylim = c(0, max_plot))
lines(model_quantiles_exp, upper_tolerance_int_const_before, lty="dashed", lwd=2, col="red")
lines(model_quantiles_exp, lower_tolerance_int_const_before, lty="dashed", lwd=2, col="red")
abline(a=0, b=1, lwd=2, col="grey")

model_quantiles_exp <- qexp(c(1:length(excesses_conservative_before))/(length(excesses_conservative_before)+1), rate = 1)
plot(model_quantiles_exp, sample_quantiles_exp_conservative_before, main = "Conservative threshold before changepoint", ylab ="Sample quantiles", xlab = "Model quantiles", pch=19, asp=1, ylim = c(0, max_plot))
lines(model_quantiles_exp, upper_tolerance_int_conservative_before, lty="dashed", lwd=2, col="red")
lines(model_quantiles_exp, lower_tolerance_int_conservative_before, lty="dashed", lwd=2, col="red")
abline(a=0, b=1, lwd=2, col="grey")

model_quantiles_exp <- qexp(c(1:length(excesses_pc_before))/(length(excesses_pc_before)+1), rate = 1)
plot(model_quantiles_exp, sample_quantiles_exp_pc_before, main = "PC threshold before changepoint", ylab ="Sample quantiles", xlab = "Model quantiles", pch=19, asp=1, ylim = c(0, max_plot))
lines(model_quantiles_exp, upper_tolerance_int_pc_before, lty="dashed", lwd=2, col="red")
lines(model_quantiles_exp, lower_tolerance_int_pc_before, lty="dashed", lwd=2, col="red")
abline(a=0, b=1, lwd=2, col="grey")

model_quantiles_exp <- qexp(c(1:length(excesses_geo_3d_before))/(length(excesses_geo_3d_before)+1), rate = 1)
plot(model_quantiles_exp, sample_quantiles_exp_geo_3d_before, main = "EQD V before changepoint", ylab ="Sample quantiles", xlab = "Model quantiles", pch=19, asp=1, ylim = c(0, max_plot))
lines(model_quantiles_exp, upper_tolerance_int_geo_3d_before, lty="dashed", lwd=2, col="red")
lines(model_quantiles_exp, lower_tolerance_int_geo_3d_before, lty="dashed", lwd=2, col="red")
abline(a=0, b=1, lwd=2, col="grey")

model_quantiles_exp <- qexp(c(1:length(excesses_log_geo_3d_before))/(length(excesses_log_geo_3d_before)+1), rate = 1)
plot(model_quantiles_exp, sample_quantiles_exp_log_geo_3d_before, main = "EQD logV before changepoint", ylab ="Sample quantiles", xlab = "Model quantiles", pch=19, asp=1, ylim = c(0, max_plot))
lines(model_quantiles_exp, upper_tolerance_int_log_geo_3d_before, lty="dashed", lwd=2, col="red")
lines(model_quantiles_exp, lower_tolerance_int_log_geo_3d_before, lty="dashed", lwd=2, col="red")
abline(a=0, b=1, lwd=2, col="grey")

model_quantiles_exp <- qexp(c(1:length(excesses_sqrt_geo_3d_before))/(length(excesses_sqrt_geo_3d_before)+1), rate = 1)
plot(model_quantiles_exp, sample_quantiles_exp_sqrt_geo_3d_before, main = "EQD sqrtV before changepoint", ylab ="Sample quantiles", xlab = "Model quantiles", pch=19, asp=1, ylim = c(0, max_plot))
lines(model_quantiles_exp, upper_tolerance_int_sqrt_geo_3d_before, lty="dashed", lwd=2, col="red")
lines(model_quantiles_exp, lower_tolerance_int_sqrt_geo_3d_before, lty="dashed", lwd=2, col="red")
abline(a=0, b=1, lwd=2, col="grey")

max_plot <- max(c(upper_tolerance_int_const_after, upper_tolerance_int_geo_3d_after, upper_tolerance_int_log_geo_3d_after, upper_tolerance_int_sqrt_geo_3d_after))

model_quantiles_exp <- qexp(c(1:length(excesses_const_after))/(length(excesses_const_after)+1), rate = 1)
plot(model_quantiles_exp, sample_quantiles_exp_const_after, main = "EQD threshold after changepoint", ylab ="Sample quantiles", xlab = "Model quantiles", pch=19, asp=1, ylim = c(0, max_plot))
lines(model_quantiles_exp, upper_tolerance_int_const_after, lty="dashed", lwd=2, col="red")
lines(model_quantiles_exp, lower_tolerance_int_const_after, lty="dashed", lwd=2, col="red")
abline(a=0, b=1, lwd=2, col="grey")

model_quantiles_exp <- qexp(c(1:length(excesses_conservative_after))/(length(excesses_conservative_after)+1), rate = 1)
plot(model_quantiles_exp, sample_quantiles_exp_conservative_after, main = "Conservative threshold after changepoint", ylab ="Sample quantiles", xlab = "Model quantiles", pch=19, asp=1, ylim = c(0, max_plot))
lines(model_quantiles_exp, upper_tolerance_int_conservative_after, lty="dashed", lwd=2, col="red")
lines(model_quantiles_exp, lower_tolerance_int_conservative_after, lty="dashed", lwd=2, col="red")
abline(a=0, b=1, lwd=2, col="grey")

model_quantiles_exp <- qexp(c(1:length(excesses_pc_after))/(length(excesses_pc_after)+1), rate = 1)
plot(model_quantiles_exp, sample_quantiles_exp_pc_after, main = "PC threshold after changepoint", ylab ="Sample quantiles", xlab = "Model quantiles", pch=19, asp=1, ylim = c(0, max_plot))
lines(model_quantiles_exp, upper_tolerance_int_pc_after, lty="dashed", lwd=2, col="red")
lines(model_quantiles_exp, lower_tolerance_int_pc_after, lty="dashed", lwd=2, col="red")
abline(a=0, b=1, lwd=2, col="grey")

model_quantiles_exp <- qexp(c(1:length(excesses_geo_3d_after))/(length(excesses_geo_3d_after)+1), rate = 1)
plot(model_quantiles_exp, sample_quantiles_exp_geo_3d_after, main = "EQD V after changepoint", ylab ="Sample quantiles", xlab = "Model quantiles", pch=19, asp=1, ylim = c(0, max_plot))
lines(model_quantiles_exp, upper_tolerance_int_geo_3d_after, lty="dashed", lwd=2, col="red")
lines(model_quantiles_exp, lower_tolerance_int_geo_3d_after, lty="dashed", lwd=2, col="red")
abline(a=0, b=1, lwd=2, col="grey")

model_quantiles_exp <- qexp(c(1:length(excesses_log_geo_3d_after))/(length(excesses_log_geo_3d_after)+1), rate = 1)
plot(model_quantiles_exp, sample_quantiles_exp_log_geo_3d_after, main = "EQD logV after changepoint", ylab ="Sample quantiles", xlab = "Model quantiles", pch=19, asp=1, ylim = c(0, max_plot))
lines(model_quantiles_exp, upper_tolerance_int_log_geo_3d_after, lty="dashed", lwd=2, col="red")
lines(model_quantiles_exp, lower_tolerance_int_log_geo_3d_after, lty="dashed", lwd=2, col="red")
abline(a=0, b=1, lwd=2, col="grey")

model_quantiles_exp <- qexp(c(1:length(excesses_sqrt_geo_3d_after))/(length(excesses_sqrt_geo_3d_after)+1), rate = 1)
plot(model_quantiles_exp, sample_quantiles_exp_sqrt_geo_3d_after, main = "EQD sqrtV after changepoint", ylab ="Sample quantiles", xlab = "Model quantiles", pch=19, asp=1, ylim = c(0, max_plot))
lines(model_quantiles_exp, upper_tolerance_int_sqrt_geo_3d_after, lty="dashed", lwd=2, col="red")
lines(model_quantiles_exp, lower_tolerance_int_sqrt_geo_3d_after, lty="dashed", lwd=2, col="red")
abline(a=0, b=1, lwd=2, col="grey")






# Helper function to calculate quantiles and bootstrapped confidence intervals
get_bootstrapped_quants <- function(excesses_before, excesses_after, n_boot) {
  
  # Function to bootstrap quantiles for a given set of excesses
  bootstrap_quants <- function(excesses, n_boot) {
    bootstrapped_quants <- matrix(NA, nrow = n_boot, ncol = length(excesses))
    for(i in 1:n_boot) {
      boot_excess <- rexp(length(excesses), rate = 1)
      bootstrapped_quants[i,] <- quantile(boot_excess, c(1:length(excesses)) / (length(excesses) + 1))
    }
    list(
      upper = apply(bootstrapped_quants, 2, quantile, prob = 0.975),
      lower = apply(bootstrapped_quants, 2, quantile, prob = 0.025)
    )
  }
  
  list(
    before = bootstrap_quants(excesses_before, n_boot),
    after = bootstrap_quants(excesses_after, n_boot)
  )
}

# Consolidated function to get excesses
get_excesses <- function(mags, threshold) {
  excesses <- mags[mags > threshold] - threshold[mags > threshold]
  return(excesses)
}

# Data preparation
prep_data <- function(gron_eq_cat, start, end, thresh_fit, dist_col) {
  data <- gron_eq_cat[start:end,]
  mags <- data$Magnitude
  third_nearest_dist <- data[[dist_col]]
  log_dist <- log(third_nearest_dist)
  sqrt_dist <- sqrt(third_nearest_dist)
  
  chosen_threshold <- thresh_fit$thresh[1] + thresh_fit$thresh[2] * third_nearest_dist
  sigma <- thresh_fit$par[1] + thresh_fit$par[2] * chosen_threshold[mags > chosen_threshold]
  
  excesses <- get_excesses(mags, chosen_threshold)
  
  return(list(
    mags = mags, 
    log_dist = log_dist, 
    sqrt_dist = sqrt_dist, 
    sigma = sigma, 
    excesses = excesses
  ))
}

# Prepare data for before and after changepoint
data_before <- prep_data(gron_eq_cat, 1, u_h_length, geo_thresh_fit_3d, "V_3d")
data_after <- prep_data(gron_eq_cat, u_h_length + 1, nrow(gron_eq_cat), geo_thresh_fit_3d, "V_3d")

# Bootstrapped quantiles
n_boot <- 200
boot_results <- get_bootstrapped_quants(data_before$excesses, data_after$excesses, n_boot)

# Sample quantiles for QQ plots
sample_quantiles_before <- quantile(transform_to_exp(data_before$excesses, sig = data_before$sigma, xi = geo_thresh_fit_3d$par[2]), 
                                    c(1:length(data_before$excesses)) / (length(data_before$excesses) + 1))
sample_quantiles_after <- quantile(transform_to_exp(data_after$excesses, sig = data_after$sigma, xi = geo_thresh_fit_3d$par[2]), 
                                   c(1:length(data_after$excesses)) / (length(data_after$excesses) + 1))

# Plotting function for both before and after
plot_qq <- function(model_quants, sample_quants, upper, lower, title) {
  plot(model_quants, sample_quants, main = title, ylab = "Sample quantiles", xlab = "Model quantiles", pch = 19, asp = 1)
  lines(model_quants, upper, lty = "dashed", lwd = 2, col = "red")
  lines(model_quants, lower, lty = "dashed", lwd = 2, col = "red")
  abline(a = 0, b = 1, lwd = 2, col = "grey")
}

# Model quantiles for QQ plot
model_quantiles_before <- qexp(c(1:length(data_before$excesses)) / (length(data_before$excesses) + 1), rate = 1)
model_quantiles_after <- qexp(c(1:length(data_after$excesses)) / (length(data_after$excesses) + 1), rate = 1)

# Plot for before and after
dev.new(height=20, width=30, noRStudioGD = TRUE)
par(mfrow=c(2,3), bg='transparent')

plot_qq(model_quantiles_before, sample_quantiles_before, boot_results$before$upper, boot_results$before$lower, "Before changepoint")

plot_qq(model_quantiles_after, sample_quantiles_after, boot_results$after$upper, boot_results$after$lower, "After changepoint")

#Image plot for EQD values across thresh_matrix

dev.new(height=20, width=20, noRStudioGD = TRUE)
par(mfrow=c(1,1), bg='transparent')

image_plot_df <- data.frame(a=threshold_matrix[,1], b=threshold_matrix[,2], v_dists=geo_thresh_fit_3d$dists, logv_dists=log_geo_thresh_fit_3d$dists, sqrtv_dists=sqrt_geo_thresh_fit_3d$dists)

ggplot(image_plot_df, aes(x=a, y=b, fill=v_dists)) + geom_tile() + scale_fill_viridis_c() + ggtitle("EQD V") + theme_minimal()
ggplot(image_plot_df, aes(x=a, y=b, fill=logv_dists)) + geom_tile() + scale_fill_viridis_c() + ggtitle("EQD logV") + theme_minimal()
ggplot(image_plot_df, aes(x=a, y=b, fill=sqrtv_dists)) + geom_tile() + scale_fill_viridis_c() + ggtitle("EQD sqrtV") + theme_minimal()

library(ggplot2)












