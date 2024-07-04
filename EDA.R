Geophones <- read.csv("Data/Geophones/Geophones_processed_03-07-2024.csv")
gron_outline <- read.csv('Data/Geophones/Groningen_Field_outline.csv', header=T)

min(Geophones$Start_date)

#Rounded (to 1 d.p.) catalogue as of 12-04-2022
rounded_cat <- read.csv("Data/Events/2022-04-12_15-09-25_cat.csv", header=T)

#New catalogue with rounding to 2 d.p. as of 04-01-2024
gron_eq_cat <- read.csv("Data/Events/gr_earthquakes_after_geophone_start_20240104.csv", header=T, row.names = 1)

# Plotting both catalogues
dev.new(height=5, width=9, noRStudioGD = TRUE)
par(mfrow=c(1,2), bg='transparent')

#Rounded to 2 d.p
plot(gron_eq_cat$Date, gron_eq_cat$Magnitude, main = "Magnitudes over time", xlab="Event time", ylab = "Magnitude")
plot(gron_eq_cat$Magnitude, main = "Magnitudes over index", xlab = "Index", ylab = "Magnitude")

#Rounded to 1 d.p.
plot(as.Date(rounded_cat$date), rounded_cat$mag, main = "Magnitudes over time", xlab="Event time", ylab = "Magnitude")
plot(rev(rounded_cat$mag), main = "Magnitudes over index", xlab = "Index", ylab = "Magnitude")

#ACF and PACF plots checking for serial dependence

acf(gron_eq_cat$Magnitude)
pacf(gron_eq_cat$Magnitude)

#Fitting a GPD with constant threshold

source("src/eqd.R")
mags <- gron_eq_cat$Magnitude

thresholds <- quantile(mags, seq(0,0.95, by=0.01))
const_thresh_fit <- eqd(mags, thresh = thresholds, k=1000)

#QQplot on Exp scale for easy comparison to more complex models

probs <- c(1:length(mags))/(length(mags)+1)
excesses <- mags[mags > const_thresh_fit$thresh] - const_thresh_fit$thresh
exp_transformed_excesses <- transform_to_exp(excesses, sig = const_thresh_fit$par[1], xi = const_thresh_fit$par[2])
sample_quantiles_exp <- quantile(exp_transformed_excesses, probs)
model_quantiles_exp <- qexp(probs, rate = 1)

n_boot <- 200
bootstrapped_model_quants <- matrix(NA, nrow = n_boot, ncol= length(probs))
for(i in 1:n_boot){
  boot_excess <- rgpd(length(excesses),scale = const_thresh_fit$par[1], shape = const_thresh_fit$par[2] )
  boot_gpd_fit <- optim(GPD_LL, par=c(mean(boot_excess), 0.1), z = boot_excess, control = list(fnscale=-1))
  exp_boot_excess <- transform_to_exp(boot_excess, sig = boot_gpd_fit$par[1], xi = boot_gpd_fit$par[2])
  bootstrapped_model_quants[i,] <- quantile(exp_boot_excess, probs)
}

upper_tolerance_int <- apply(bootstrapped_model_quants, 2, quantile, prob = 0.975)
lower_tolerance_int <- apply(bootstrapped_model_quants, 2, quantile, prob = 0.025)
dev.new(height=5, width=5, noRStudioGD = TRUE)
par(mfrow=c(1,1), bg='transparent')

plot(sample_quantiles_exp, model_quantiles_exp, main = "QQ-plot for Std GPD fit on Exp scale", xlab ="Sample quantiles", ylab = "Model quantiles", pch=19)
lines(model_quantiles_exp, upper_tolerance_int, lty="dashed", lwd=2)
lines(model_quantiles_exp, lower_tolerance_int, lty="dashed", lwd=2)
abline(a=0, b=1, lwd=2)
