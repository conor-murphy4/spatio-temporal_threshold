## Poisson process moedl using linear ICS and optimal smoothing

print("This is S4_fitting.R.")

##
## 00: Set up envrionment. Load packages, scripts and data. --------------------
##

PLOT_PATH  <- "Point processes/Covariates/Output/plots/"
DATA_PATH  <- "Point processes/Covariates/Output/data/"
MAKE_PLOTS <- TRUE
SAVE_PLOTS <- TRUE
RUN_ALL <- TRUE    # run all code or load results for slow parts?

## Install required packages ------------
library(lubridate)
library(sp)
library(MASS)
library(fields)
library(mgcv)
library(purrr)

## Source required scripts --------------
source("Point processes/Covariates/S4_ics_exponential/S4_functions.R")

## Read in relevant data ----------------

# Groningen field outline
SPgfo <- readRDS("Point processes/Covariates/00_data/derived/field_outline/SPgfo.RDS")
inPolyVec <- readRDS("Point processes/Covariates/00_data/derived/inPoly/inPolyVec.RDS")
inPolyMat <- readRDS("Point processes/Covariates/00_data/derived/inPoly/inPolyMat.RDS")

# incremental Coulomb Stress covariate
ics_mats <- readRDS("Point processes/Covariates/00_data/derived/covariates/incremental_Coulomb_stress_mats.RDS")
ics_mats <- subset_to_years_cov_mats(ics_mats, 1995, 2017)

#ics_rate_mats <- temporal_difference_cov_mats(ics_mats, "ics_rate")
#ics_rate_mats <- replace_negatives_cov_mats(ics_rate_mats,replace = TRUE, replacement_value = 1e-6)
#plot_cov_mats(ics_rate_mats, zlim = c(0,0.4))

# Earthquake catalogue
EqCat <- read.csv('Data/Events/2022-04-12_15-09-25_cat.csv', header = TRUE)
colnames(EqCat)[1:5] <- c("Date", "Time_UT", "Easting", "Northing", "Magnitude" )
Cat <- Catalogue(EC = EqCat, minMag = 1.5, maxDate = '2016-12-31')
rm(EqCat)



##
## 01: Model Fitting -----------------------------------------------------------
##
debugonce(llik_S4)
# params := c(beta_0, beta_1, sigma_cumulative)
llik_S4(params = c(0.1,0.01, 500),
            catalogue = Cat,
            cov_mats = ics_mats,
            inpolymat = inPolyMat,
            negative = TRUE)


# Minimise negative log-likelihood ------------
Opt<- optim(par = c(0.1, 0.01, 500),
             fn = llik_S4,
             #method = 'L-BFGS-B',
             #lower = c(0.3, 1000),
             #upper = c(0.5, 4000),
             catalogue = Cat,
             cov_mats = ics_mats,
             inpolymat = inPolyMat,
             negative = TRUE,
             hessian = TRUE
)

rds_path <- paste0(DATA_PATH, "Opt_S4.RDS")
saveRDS(object = Opt, file = rds_path)

# Point estimates
Opt$par
# (very) approximate standard errors
sqrt(diag(solve(Opt$hessian)))
# Log-likelhood value at the mle
-Opt$value


# Testing profile likelihood
prof_llik_S4(
  params = c(0.1,0.01),
  catalogue = Cat,
  cov_mats = ics_mats,
  inpolymat = inPolyMat,
  negative = FALSE)

# EDITED TO HERE 2021-04-19

# Deviance based CI   -------------------------
alpha = 0.05
critical_value <- 0.5 * qchisq(p = 1 - alpha, df = 2)

b0 <- seq(0.0001,0.1, length.out = 21)
b1 <- seq(-0.085, 1.5, length.out = 23)
pars <- expand.grid(b0, b1)

prof_llikvals <- matrix(NA, nrow = length(b0), ncol = length(b1))

rds_path <- paste0(DATA_PATH, "prof_llikvals_S2.rds")
if (RUN_ALL) {

  for (i in seq_along(b0)) {
    for (j in seq_along(b1)) {
      prof_llikvals[i,j] <- prof_llik_S2(
        params = c(b0[i], b1[j]),
        catalogue = Cat,
        cov_mats = ics_mats,
        inpolymat = inPolyMat,
        negative = FALSE,
        detailed = FALSE
      )
      print(paste("(",i,",",j,")"))
    }
  }

  saveRDS(object = list(beta_0 = b0, beta_1 = b1, prof_llikvals = prof_llikvals),
          file = rds_path)
} else {
  prof_llikvals <- readRDS(rds_path)[[3]]
}
## UP TO HERE!! ---------------------------------------------------------------------
pdf_path = paste0(PLOT_PATH, 'S2_likelihood.pdf')
if( MAKE_PLOTS) {

  if(SAVE_PLOTS){ pdf(pdf_path, width = 7, height = 5) }

  contour(
    x = b0,
    y = b1,
    z = prof_llikvals,
    xlab = "beta_0",
    ylab = "beta_1",
    main = "Poisson Process log-likelihood (S2 ICS linear)",
    levels = - (Opt$value + seq(0, 10, by = 1)),
    xlim = c(0.025,0.055),
    ylim = c(-0.1,0))

  points(x = Opt$par[1], y = Opt$par[2])

  contour(
    x = b0,
    y = b1,
    z = prof_llikvals,
    levels = - (Opt$value + critical_value),
    col = 2,
    add = TRUE)

  abline(h = -1 / max(ics_mats$ics, na.rm = TRUE), col = 2, lty = 2)
  if(SAVE_PLOTS){ dev.off() }
}

#
# 02: Model evaluation ---------------------------------------------------------
#

ics_rate_mats_opt <- smooth_cov_mats(
  cov_mats = ics_rate_mats,
  lambda = Opt$par[2],
  inpolymat = inPolyMat)

debugonce(llik_S1)
llik_S1(params = Opt$par,
        catalogue = Cat,
        cov_rate_mats = ics_rate_mats,
        inpolymat = inPolyMat,
        negative = FALSE)

#array of expected event counts
Ecountmats <- Opt$par[1] * ics_rate_mats_opt[[1]] * 0.5^2
Lambda <- list(Ecountmats,
               ics_rate_mats_opt$xgrid,
               ics_rate_mats_opt$ygrid,
               ics_rate_mats_opt$tgrid)
names(Lambda) <- c('Emats', 'xgrid', 'ygrid', 'tgrid')
rm(Ecountmats)


# Plot annual expected event count per box -------------------

pdf_path <- paste0(PLOT_PATH, "S1_annual_expected_count_maps.pdf")

if(MAKE_PLOTS){

  if(SAVE_PLOTS){
    pdf(file = pdf_path, width = 7 , height = 7)
    par(mfrow = c(1,1), omi = c(0.5,0.5,0.5,0.5))
  }

  for (yr in 1:22){

    YEAR = yr + 1994
    date_condition <- year(Cat$ld) == YEAR
    mag_condition  <- Cat$mag >= 1.5

    events_E <- Cat$x[  date_condition & mag_condition]
    events_N <- Cat$y[  date_condition & mag_condition]
    events_M <- Cat$mag[date_condition & mag_condition]

    expected_count <- sum(Lambda$Emats[,,yr], na.rm = TRUE)
    fitted_CI <- qpois(lambda = expected_count, p =  c(0.025,0.975))
    observed_count <- length(events_E)

    YEAR
    rounded_fitted_count <- round(expected_count, 2)
    rounded_fitted_CI <- round(fitted_CI, 3)

    image.plot(
      x = Lambda$xgrid,
      y = Lambda$ygrid,
      z = Lambda$Emats[,,yr] / (0.5 ^ 2),
      asp= 1,
      xlab = 'Easting',
      ylab = 'Northing',
      #zlim = range(Lambda$Emats[,,yr] / (0.5 ^ 2), na.rm = TRUE), #high contrast
      zlim = range(Lambda$Emats / (0.5 ^ 2), na.rm = TRUE), #comparable years
      main = paste0('Expected EQ count ',
                  YEAR,
                  ". Obs = ",
                  observed_count,
                  ", Exp = ",
                  rounded_fitted_count
                  #"(",
                  #rounded_fitted_CI[1],
                  #",",
                  #rounded_fitted_CI[2],
                  #")."
                  )
    )
    points(x = events_E, y = events_N, cex = sqrt(events_M / pi))
    plot(SPgfo, add= TRUE, asp =1)
  }

  if(SAVE_PLOTS){ dev.off() }
}

# Compare annual predicted and observed event counts ------------
# NB: Now takes into account parameter uncertainty

pred.annual <- apply(X = Lambda$Emats, MARGIN = 3, FUN = sum, na.rm = TRUE)
eventsinyear <- function(year){sum(Cat$Year == year)}
obsv.annual <- sapply(X = 1995:2016, FUN = eventsinyear)


Eyr <- function(param, cov_rate_mats = ics_rate_mats, inpolymat = inPolyMat){
  #for a sampled parameter, create a vector of expected counts
  cov_rate_mats_smoothed <- smooth_cov_mats(
    cov_mats = cov_rate_mats,
    lambda = param[2],
    inpolymat = inPolyMat)

  Emats <- (param[1] * cov_rate_mats_smoothed[[1]] ) * 0.5^2
  apply(X = Emats, MARGIN = 3, FUN = sum, na.rm=TRUE)
}

CI<- function(params_mle, params_hessian, nsim = 1000, cov_rate_mats = ics_rate_mats, inpolymat = inPolyMat){
  #vector of sampled parameter values
  Sigma <- solve(params_hessian)
  param_sim <- mgcv::rmvn(nsim, mu = params_mle, V = Sigma)
  param_sim <- param_sim[param_sim[,2]>0, ]
  print(paste0("Using ",nrow(param_sim), " valid parameter values from ", nsim, " simulated values."))
  # for each sampled parameter, create a vector of expected counts
  Eyear <- apply(X = param_sim, MARGIN = 1, FUN = Eyr, cov_rate_mats = cov_rate_mats, inpolymat = inpolymat)
  #for each vector, simulate counts from poisson
  sim_year<- apply(X = Eyear, MARGIN = 2, FUN = rpois, n=dim(Eyear)[1])
  #for each year calculate quantiles
  CI_year <-apply(X =sim_year, MARGIN = 1,FUN = quantile, probs=c(0.025,0.975))
  return(CI_year)
}

# This takes a little while - look to optimise when have the time.
rds_path <- paste0(DATA_PATH, "pred_ci.rds")
if (RUN_ALL) {
  pred_ci <- CI(
    params_mle = Opt$par,
    params_hessian =  Opt$hessian,
    nsim = 1000,
    cov_rate_mats = ics_rate_mats,
    inpolymat = inPolyMat
  )

  saveRDS(pred_ci, file = rds_path)
} else {
  pred_ci <- readRDS(rds_path)
}


# plot annual predicted and annual observed event counts
pdf_path = paste0(PLOT_PATH,"S1_annual_counts_fitted.pdf")
if (MAKE_PLOTS){
  if(SAVE_PLOTS){
    pdf(file = pdf_path, width = 8, height = 5)
    par(mfrow = c(1,1), omi = c(0,0,0,0), mar = c(4.5,4.5,1,1))
  }
  plot(y = obsv.annual, x = 1995:2016, xlab = 'year', ylab = 'Event Count', ylim = c(0,35), pch = 16)
  lines(y = pred_ci[1,], x = 1995:2017, lty = 2)
  lines(y = pred_ci[2,], x = 1995:2017, lty = 2)
  lines(y = pred.annual, x = 1995:2017)
  if(SAVE_PLOTS){ dev.off()}
}


# Simulate a catalogue from the fitted model ----------------

#expected and simulated total number of events
Eevents <- sum(Lambda$Emats, na.rm = TRUE)
nevents <- rpois(lambda = Eevents, n = 1)
#assign each to box prop to its contribution
box.probs <- as.vector(Lambda$Emats)/Eevents
box.probs[is.na(box.probs)] <- 0
box.no <- length(box.probs)
box.smpl <- sample(x = 1:box.no,
                   size = nevents,
                   prob =  box.probs, replace = TRUE)
#convert ot space-time position of box
sim.cat <- box.no2xyt(box.smpl,minX = min(ics_rate_mats$xgrid),
                      minY = min(ics_rate_mats$ygrid), minYr = 1995)
sim.cat<- as.data.frame(sim.cat)
#loacte each event randomly within assigned box
xJtr  <- (runif(nevents) - 0.5)*500
yJtr  <- (runif(nevents) - 0.5)*500
tJtr  <- runif(nevents)

sim.cat$x <- sim.cat$x + xJtr
sim.cat$y <- sim.cat$y + yJtr
sim.cat$t <- sim.cat$t + tJtr

# Plot simulated catalogue  ---------------------------
pdf_path <- paste0(PLOT_PATH, "S1_simulated_catalogue.pdf")
if (MAKE_PLOTS){
  if(SAVE_PLOTS){
    pdf(file = pdf_path, width = 8, height = 5)
    par(mfrow=c(1,2), mar = c(1,1,4.5,1))
  }
  plot(SPgfo, asp =1, main = 'Simulated')
  points(sim.cat[,1:2])
  plot(SPgfo, asp =1, main = 'Observed')
  points(Cat$x, Cat$y)
  if(SAVE_PLOTS){ dev.off() }
}


## EOF ---------------
