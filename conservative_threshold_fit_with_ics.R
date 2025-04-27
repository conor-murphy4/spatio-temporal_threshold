source("src/helper_functions.R")


gron_eq_cat <- read.csv("Data/Events/unrounded_after_1995_in_polygon_with_covariates.csv", header=T)
mags <- gron_eq_cat$Magnitude

conservative_threshold <- rep(1.45, length(mags))

get_table_of_results <- function(mags, threshold, ics){
  thresh_par <- unique(threshold)
  par_ests <- get_par_ests_const_ics(mags, threshold, ics)
  num_excess <- sum(mags > threshold)
  return(list(thresh_par, par_ests, num_excess))
}

get_table_of_results(mags, conservative_threshold, ics = gron_eq_cat$ICS)

conserv_thresh_fit <- get_par_ests_const_ics(mags, conservative_threshold, ics = gron_eq_cat$ICS)

# Intensity estimation

(opt_PP <- optim(Poisson_process_LL, par = c(0.1, 0), data = gron_eq_cat, covariates = covariates, 
                 threshold_obs = conservative_threshold, thresh_fit = conserv_thresh_fit, control=list(fnscale=-1), hessian = T))

cov_matrix <- solve(-1*opt_PP$hessian)
se <- sqrt(diag(cov_matrix))
z <- qnorm(0.975)
lower <- opt_PP$par - z*se
upper <- opt_PP$par + z*se
(CI_gamma0 <- round(c(lower[1], upper[1]), 3) )
(CI_gamma1 <- round(c(lower[2], upper[2]),3) )


# Results from V1 threshold
get_par_ests_geo_ics(mags, thresh_fit = thresh_fit_V1_ics, V = gron_eq_cat$V_1,
                     ics = gron_eq_cat$ICS, min_ics = min(covariates$ICS), max_dist = max(covariates$V1), 
                     show_fit = F)

(opt_PP <- optim(Poisson_process_LL, par = c(0.1, 0), data = gron_eq_cat, covariates = covariates, 
                 threshold_obs = threshold_obs_V1, thresh_fit = thresh_fit_V1_ics, control=list(fnscale=-1), hessian = T))

cov_matrix <- solve(-1*opt_PP$hessian)
se <- sqrt(diag(cov_matrix))
z <- qnorm(0.975)
lower <- opt_PP$par - z*se
upper <- opt_PP$par + z*se
(CI_gamma0 <- round(c(lower[1], upper[1]), 3) )
(CI_gamma1 <- round(c(lower[2], upper[2]),3) )
