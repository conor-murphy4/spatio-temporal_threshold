setwd("C:/Users/murphyc4/OneDrive - Lancaster University/STOR-i/PhD/Projects/Spatio-temporal")

library(Metrics)
library(plotrix)
library(devtools)

#Proposed set of theta (0,2), true theta (0,1.9)
Est_scale_gpd <- readRDS("Comparison_GPD_EXP/Est_scale_gpd.rds")
Est_shape_gpd <- readRDS("Comparison_GPD_EXP/Est_shape_gpd.rds")
Est_theta_gpd <- readRDS("Comparison_GPD_EXP/Est_theta_gpd.rds")

Est_scale_exp <- readRDS("Comparison_GPD_EXP/Est_scale_exp.rds")
Est_shape_exp <- readRDS("Comparison_GPD_EXP/Est_shape_exp.rds")
Est_theta_exp <- readRDS("Comparison_GPD_EXP/Est_theta_exp.rds")

#Results with higher range of proposed theta. proposed (0,2.5)

Est_scale_gpd <- readRDS("STORM/Est_scale_gpd.rds")
Est_shape_gpd <- readRDS("STORM/Est_shape_gpd.rds")
Est_theta_gpd <- readRDS("STORM/Est_theta_gpd.rds")

Est_scale_exp <- readRDS("STORM/Est_scale_exp.rds")
Est_shape_exp <- readRDS("STORM/Est_shape_exp.rds")
Est_theta_exp <- readRDS("STORM/Est_theta_exp.rds")


scale <- 20
shape <- 0.1
n_samples <- 200
n_theta <- 20
num_events <- 1960

thetas <- seq(0,1.9, by=0.1)
thetas_rep <- rep(seq(0,1.9, by=0.1), each=n_samples)

#--------------------------Estimated theta against true theta-------------------
dev.new(width=9.17, height=5.53,noRStudioGD = T)
par(mfrow=c(1,1),bg='transparent')
#GPD
sizeplot(thetas_rep, Est_theta_gpd, xlab=expression(theta), ylab=expression(hat(theta)))
abline(a=0,b=1, col="red", lwd=2.5)
plot(thetas_rep, Est_theta_gpd, xlab=expression(theta), ylab=expression(hat(theta)))
abline(a=0,b=1, col="red", lwd=2.5)
#EXP
sizeplot(thetas_rep, Est_theta_exp, xlab=expression(theta), ylab=expression(hat(theta)))
abline(a=0,b=1, col="red", lwd=2.5)
plot(thetas_rep, Est_theta_exp, xlab=expression(theta), ylab=expression(hat(theta)))
abline(a=0,b=1, col="red", lwd=2.5)




#---------------------------Hypthesis Test for parameters-----------------------

#h0: b = 0, h1: b /= 0

##-------GPD-----------
#Scale parameter
b_scale_gpd <- apply(Est_scale_gpd, 2, function(x,z){x-z}, z= scale)
sd_scale_gpd <- apply(b_scale_gpd, 2, sd)
h0_upper <- 1.96*sd_scale_gpd/sqrt(n_samples)
h0_lower <- -1.96*sd_scale_gpd/sqrt(n_samples)

E_bias_scale_gpd <- apply(Est_scale_gpd, 2, bias, actual=scale)

dev.new(width=9.17, height=5.53,noRStudioGD = TRUE)
par(mfrow=c(1,1),bg='transparent')
plot(thetas, E_bias_scale_gpd, ylab = "Bias", xlab = expression(theta), ylim=c(-3, 1.5))
abline(h=0)
lines(thetas, h0_lower, lty="dashed", col="blue")
lines(thetas, h0_upper, lty="dashed", col="blue")

#Shape parameter
b_shape_gpd <- apply(Est_shape_gpd, 2, function(x,z){x-z}, z= shape)
sd_shape_gpd <- apply(b_shape_gpd, 2, sd)
h0_upper <- 1.96*sd_shape_gpd/sqrt(n_samples)
h0_lower <- -1.96*sd_shape_gpd/sqrt(n_samples)

E_bias_shape_gpd <- apply(Est_shape_gpd, 2, bias, actual=shape)

plot(thetas, E_bias_shape_gpd, ylab = "Bias", xlab = expression(theta), ylim=c(-0.1,0.1))
abline(h=0)
lines(thetas, h0_lower, lty="dashed", col="blue")
lines(thetas, h0_upper, lty="dashed", col="blue")

##------EXP------------

#Scale parameter
b_scale_exp <- apply(Est_scale_exp, 2, function(x,z){x-z}, z= scale)
sd_scale_exp <- apply(b_scale_exp, 2, sd)
h0_upper <- 1.96*sd_scale_exp/sqrt(n_samples)
h0_lower <- -1.96*sd_scale_exp/sqrt(n_samples)

E_bias_scale_exp <- apply(Est_scale_exp, 2, bias, actual=scale)

plot(thetas, E_bias_scale_exp, ylab = "Bias", xlab = expression(theta), ylim=c(-3, 1.5))
abline(h=0)
lines(thetas, h0_lower, lty="dashed", col="blue")
lines(thetas, h0_upper, lty="dashed", col="blue")

#Shape parameter
b_shape_exp <- apply(Est_shape_exp, 2, function(x,z){x-z}, z= shape)
sd_shape_exp <- apply(b_shape_exp, 2, sd)
h0_upper <- 1.96*sd_shape_exp/sqrt(n_samples)
h0_lower <- -1.96*sd_shape_exp/sqrt(n_samples)

E_bias_shape_exp <- apply(Est_shape_exp, 2, bias, actual=shape)

plot(thetas, E_bias_shape_exp, ylab = "Bias", xlab = expression(theta), ylim=c(-0.1,0.1))
abline(h=0)
lines(thetas, h0_lower, lty="dashed", col="blue")
lines(thetas, h0_upper, lty="dashed", col="blue")


#-------------95% CIs based on distribution of estimates------------------------
#Do the CIs cover 0?

#---------GPD-----------

#Scale parameter
Est_scale_gpd_sorted <- apply(Est_scale_gpd, 2, sort)
Est_scale_gpd_l <- Est_scale_gpd_sorted[0.025*n_samples,]
Est_scale_gpd_u <- Est_scale_gpd_sorted[0.975*n_samples,]

dev.new(width=9.17, height=5.53,noRStudioGD = TRUE)
par(mfrow=c(1,1),bg='transparent')
plot(thetas, Est_scale_gpd_l-scale, ylim=c(-10, 30),type='l', ylab = expression(hat(sigma)-sigma), xlab=expression(theta))
lines(thetas, Est_scale_gpd_u-scale)
abline(h=0, col="blue")

#Shape parameter
Est_shape_gpd_sorted <- apply(Est_shape_gpd, 2, sort)
Est_shape_gpd_l <- Est_shape_gpd_sorted[0.025*n_samples,]
Est_shape_gpd_u <- Est_shape_gpd_sorted[0.975*n_samples,]


plot(thetas, Est_shape_gpd_l-shape,type='l',ylim=c(-0.5, 0.5), ylab = expression(hat(xi)-xi), xlab=expression(theta))
lines(thetas, Est_shape_gpd_u-shape)
abline(h=0, col="blue")


#----------EXP-----------

#Scale parameter
Est_scale_exp_sorted <- apply(Est_scale_exp, 2, sort)
Est_scale_exp_l <- Est_scale_exp_sorted[0.025*n_samples,]
Est_scale_exp_u <- Est_scale_exp_sorted[0.975*n_samples,]


plot(thetas, Est_scale_exp_l-scale, ylim=c(-10, 30),type='l', ylab = expression(hat(sigma)-sigma), xlab=expression(theta))
lines(thetas, Est_scale_exp_u-scale)
abline(h=0, col="blue")

#Shape parameter
Est_shape_exp_sorted <- apply(Est_shape_exp, 2, sort)
Est_shape_exp_l <- Est_shape_exp_sorted[0.025*n_samples,]
Est_shape_exp_u <- Est_shape_exp_sorted[0.975*n_samples,]

plot(thetas, Est_shape_exp_l-shape,type='l',ylim=c(-0.5, 0.5), ylab = expression(hat(xi)-xi), xlab=expression(theta))
lines(thetas, Est_shape_exp_u-shape)
abline(h=0, col="blue")
