#=====================================================================
# Functions for Generalised Pareto Distribution.
# Added option to use nu parameterisation
# checks that param values are valid
#=====================================================================
# pgpd
# qgpd
# dgpd
# rgpd

# GPD_LL
#=====================================================================
#' Generalised Pareto Distribution
#'
#' Cumulative density function of the GPD specified in terms of (sigma,xi) or (nu,xi).
#' Improvement over evir function as it returns an error if the (implied) shape
#' parameter is non-positive. Also properly handles cases where and xi=0 or p<mu.
#'
#' @author Zak Varty
#'
#' @param q vector of quantiles.
#' @param shape shape parameter (xi)
#' @param scale scale parameter (sigma)
#' @param nu  alternative scale parameter: nu = sigma/(1+xi)
#' @param mu  location parameter
#' @param skip_checks logical. Speed up evaluation by skipping checks on inputs? (Beware!)
#' @return Probability of the GPD X<=q
#' @importFrom stats pexp
#' @examples
#' pgpd(q = c(-1,1.5,3), shape = 1, scale = 1)
#' pgpd(q = 1.5, shape = c(0,-1), scale = c(0.1,1))
#' @export

pgpd <- function(q, shape, scale = NULL, nu = NULL, mu = 0, skip_checks = FALSE){
  
  if (!skip_checks) {
    # one and only one of {nu, scale} may be specified
    if (is.null(scale) & is.null(nu)) {
      stop('Define one of the parameters nu or scale.')
    }
    if (!is.null(scale) & !is.null(nu)) {
      stop('Define only one of the parameters nu and scale.')
    }
    # Calculate scale from nu if required
    if (!is.null(nu) & is.null(scale)) {
      scale <- nu / (1 + shape)
      if (any(scale <= 0)) {
        stop('Implied scale parameter(s) must be positive.')
      }
      
    }
    # Check that scale value(s) are positive
    if (any(scale <= 0)) {
      stop('Scale parameter(s) must be positive.')
    }
    
    # Ensure q, scale, shape and mu are of same length.
    if (length(scale) == 1 & length(q) > 1) {
      scale <- rep(scale, length(q))
    }
    if (length(shape) == 1 & length(q) > 1) {
      shape <- rep(shape, length(q))
    }
    if (length(mu) == 1 & length(q) > 1) {
      mu <- rep(mu, length(q))
    }
  } else {
    if (!is.null(nu) & is.null(scale)) {
      scale <- nu / (1 + shape)
    }
  }
  #calculate probabilities
  p <- (1 - (1 + (shape * (q - mu))/scale)^(-1/shape))
  #correct probabilities below mu or above upper end point
  p[q < mu] <- 0
  p[(shape < 0) & (q >= (mu - scale/shape))] <- 1
  
  #correct probabilities where xi = 0
  if (any(abs(shape) < 1e-10)) {
    #ex <- which(shape ==0)
    ex <- which(abs(shape) < 1e-10)
    p[ex] <- pexp(q = q[ex] - mu[ex], rate = 1 / scale[ex])
  }
  
  return(p)
}

#' Generalised Pareto Distribution
#'
#' Cumulative density function of the GPD specified in terms of (sigma,xi) or (nu,xi).
#' Improvement over evir function as it returns an error if the (implied) shape
#' parameter is non-positive. Also properly handles cases where and xi=0 or p is not a valid
#' probability.
#'
#' @author Zak Varty
#'
#' @param p vector of quantiles.
#' @param shape shape parameter (xi)
#' @param scale scale parameter (sigma)
#' @param nu  alternative scale parameter: nu = sigma/(1+xi)
#' @param mu  location parameter
#' @return Probability of the GPD X<=x
#' @examples
#' qgpd(p = 0.5, shape = 0.5, scale = 0.5)
#' \dontrun{ qgpd(p = -0.1, shape = 0, scale = 1, mu = 0.1) }
#' @export
qgpd <- function(p, shape, scale = NULL, nu = NULL, mu = 0){
  # one and only one of {nu, scale} may be specified
  if (is.null(scale) & is.null(nu)) {
    stop('Define one of the parameters nu or scale.')
  }
  if (!is.null(scale) & !is.null(nu)) {
    stop('Define only one of the parameters nu and scale.')
  }
  
  # Probabilities must all be positive
  if (!all((p >= 0) & (p <= 1))) {
    stop('Probabilities p must be in the range [0,1].')
  }
  
  # Calculate scale from nu if required
  if (!is.null(nu) & is.null(scale)) {
    scale <- nu / (1 + shape)
    if (any(scale <= 0)) {
      stop('Implied scale parameter(s) must be positive.')
    }
    
  }
  
  # Check that scale value(s) are positive
  if (any(scale <= 0)) {
    stop('Scale parameter(s) must be positive.')
  }
  # Ensure p, scale, shape and mu are of same length.
  if (length(scale) == 1 & length(p) > 1) {
    scale <- rep(scale, length(p))
  }
  if (length(shape) == 1 & length(p) > 1) {
    shape <- rep(shape, length(p))
  }
  if (length(mu) == 1 & length(p) > 1) {
    mu <- rep(mu, length(p))
  }
  
  #calculate quantiles
  q <- mu + (scale/shape) * ((1 - p)^(-shape) - 1)
  
  #correct quantiles where xi = 0
  #ex <- which(shape ==0)
  if (any(abs(shape) < 1e-10)) {
    ex <- which(abs(shape) < 1e-10)
    q[ex] <- mu[ex] + stats::qexp(p = p[ex],rate = 1/scale[ex])
  }
  return(q)
}



#' Generalised Pareto Distribution
#'
#' Density function of the GPD specified in terms of (sigma,xi) or (nu,xi).
#' Improvement over evir function as it returns an error if the (implied) shape
#' parameter is non-positive. Also properly handles cases where and xi=0 or x is
#' outside of the domain of the given distribution.
#'
#' @author Zak Varty
#'
#' @param x vector of values as which to evaluate density.
#' @param shape shape parameter (xi)
#' @param scale scale parameter (sigma)
#' @param nu  alternative scale parameter
#' @param mu  location parameter
#' @param log  locical. Return log
#' @return density of the GPD at x
#' @examples
#' dgpd(x = c(-1,0.5,1,1.9,5),shape = -0.5, scale = 1)
#' @export
#'
dgpd <- function(x, shape, scale = NULL, nu = NULL, mu = 0, log = FALSE){
  # one and only one of {nu, scale} may be specified
  if (is.null(scale) & is.null(nu)) {
    stop('Define one of the parameters nu or scale.')
  }
  if (!is.null(scale) & !is.null(nu)) {
    stop('Define only one of the parameters nu and scale.')
  }
  
  # Calculate scale from nu if required
  if (!is.null(nu) & is.null(scale)) {
    scale <- nu / (1 + shape)
    if (any(scale <= 0)) {
      stop('Implied scale parameter(s) must be positive.')
    }
  }
  
  # Check that scale value(s) are positive
  if (any(scale <= 0)) {
    stop('Scale parameter(s) must be positive.')
  }
  # Ensure x, scale, shape and mu are of same length.
  if (length(scale) == 1 & length(x) > 1) {
    scale <- rep(scale, length(x))
  }
  if (length(shape) == 1 & length(x) > 1) {
    shape <- rep(shape, length(x))
  }
  if (length(mu) == 1 & length(x) > 1) {
    mu <- rep(mu, length(x))
  }
  
  if (log == FALSE) {
    out <- (scale^(-1)) * pmax((1 + shape * (x - mu)/scale),0)^((-1/shape) - 1)
    # amend values below threshold
    out[which(x < mu)] <- 0
    # amend values above upper endpoint (if it exists)
    out[which((shape < 0) & (x >= (mu - scale/shape)))] <- 0
    # amend values where xi = 0 (if they exist)
    if (any(abs(shape < 1e-10))) {
      ex <- which(abs(shape) < 1e-10)
      out[ex] <- stats::dexp(x = x[ex] - mu[ex], rate = 1/scale[ex])
    }
  } else {
    out <-  -log(scale) + ((-1/shape) - 1)*log(pmax((1 + shape * (x - mu)/scale),0))
    # amend values below threshold
    out[which(x < mu)] <- -Inf
    # amend values above upper endpoint (if it exists)
    out[which((shape < 0) & (x >= (mu - scale/shape)))] <- -Inf
    # amend values where xi = 0 (if they exist)
    if (any(abs(shape) < 1e-10)) {
      ex <- which(abs(shape) < 1e-10)
      out[ex] <- stats::dexp(x = x[ex] - mu[ex], rate = 1 / scale[ex],log = TRUE)
    }
  }
  return(out)
}

#' Generalised Pareto Distribution
#'
#' Sample the GPD specified in terms of (sigma,xi) or (nu,xi).
#' Improvement over evir function as it returns an error if the (implied) shape
#' parameter is non-positive. Also properly handles cases where and xi=0.
#'
#' @author Zak Varty
#'
#' @param n sample size.
#' @param shape shape parameter (xi).
#' @param scale scale parameter (sigma).
#' @param nu  alternative scale parameter.
#' @param mu  location parameter.
#' @return Random sample from generalised pareto distirbution.
#'
#' @examples
#' rgpd(n = 100, shape = 0, scale = 1:100)
#' @export
rgpd <- function(n, shape, scale = NULL, nu = NULL, mu = 0){
  ## Input checks
  # one and only one of {nu, scale} may be specified
  if (is.null(scale) & is.null(nu)) {
    stop('Define one of the parameters nu or scale.')
  }
  if (!is.null(scale) & !is.null(nu)) {
    stop('Define only one of the parameters nu and scale.')
  }
  # Calculate scale from nu if required
  if (!is.null(nu) & is.null(scale)) {
    scale <- nu / (1 + shape)
    if (any(scale <= 0)) {
      stop('Implied scale parameter(s) must be positive.')
    }
  }
  # Check that scale value(s) are positive
  if (any(scale <= 0)) {
    stop('Scale parameter(s) must be positive.')
  }
  # Ensure q, scale, shape and mu are of same length.
  if ((length(scale) == 1) & (n > 1)) {
    scale <- rep(scale, n)
  }
  if ((length(shape) == 1) & (n > 1)) {
    shape <- rep(shape, n)
  }
  if ((length(mu) == 1) & (n > 1)) {
    mu <- rep(mu, n)
  }
  
  #simulate sample
  sample <- mu + (scale/shape) * ((1 - stats::runif(n))^(-shape) - 1)
  #correct sample values where xi = 0
  #ex <- which(shape ==0)
  if (any(abs(shape) < 1e-10)) {
    ex <- which(abs(shape) < 1e-10)
    sample[ex] <- mu[ex] +
      stats::rexp(n = length(ex),rate = 1/scale[ex])
  }
  return(sample)
}


#' Generalised Pareto log-likelihood
#'
#' @author Conor Murphy
#'
#' @param par A numeric vector of parameter values of length 2.
#' @param z A numeric vector of excesses of some threshold.
#'
#' @returns A numeric value of the log-likeihood.
#'
#' @examples
#' test1 <- rgpd(1000, shape = 0.1, scale=0.5, mu=1)
#' excess <- test1[test1>1.5] - 1.5
#' GPD_LL(par=c(1,0.4), z=excess)
#' 
GPD_LL <- function(par, z){
  sigma <- par[1]
  xi <- par[2]
  if (sigma < 0) {return(-1e7)}
  
  if (abs(xi) < 1e-10) {
    return(-length(z) * log(sigma) - ((1 / sigma) * sum(z)))
  }
  
  if (any(1 + (xi * z) / sigma <= 0)) {return(-1e6)}
  
  LL1 <- -(length(z) * log(sigma))
  LL2 <- -((1 + 1 / xi) * (sum(log(1 + (xi * z) / sigma))))
  return(LL1 + LL2)    
}

#' Transform GPD excesses to standard exponential
#' 
#' @author Conor Murphy
#' 
#' @param y A numeric vector of excesses of some threshold.
#' @param sig A numeric value or vector of values of the scale parameter.
#' @param xi A numeric value or vector of values of the shape parameter.
#' 
#' @returns A numeric vector of transformed excesses.
transform_to_exp <- function(y, sig, xi){
  
  if (length(sig) != length(y) & length(sig) != 1) {stop("sig must be of length 1 or length(y)")}
  if (length(xi) != length(y) & length(xi) != 1) {stop("xi must be of length 1 or length(y)")}
  
  std_exp <- (1 / xi) * log( 1 + xi * (y  /sig))  
  return(std_exp)
}

#' Obtain MLEs, SEs and CIs for standard GPD parameters
#'
#' @author Conor Murphy
#'
#' @param mags A numeric vector of magnitudes.
#' @param threshold A numeric value of the threshold.
#' @param show_fit Logical. If TRUE, the full fit object is returned.
#'
#' @returns A list with MLEs, SEs and CIs for the scale and shape parameters.
get_par_ests <- function(mags, threshold, show_fit = TRUE){
  excesses <- mags[mags > threshold] - threshold
  mle0 <- mean(excesses)
  
  gpd_fit <- optim(GPD_LL,
                   z = excesses,
                   par = c(mle0,0.1),
                   control = list(fnscale = -1),
                   hessian = TRUE)
  
  hessian <- gpd_fit$hessian
  cov_matrix <- solve(-1*hessian)
  se <- sqrt(diag(cov_matrix))
  z <- qnorm(0.975)
  lower <- gpd_fit$par - z * se
  upper <- gpd_fit$par + z * se
  
  output <- list(
    par = gpd_fit$par,
    SE = se,
    CI_scale = round(c(lower[1], upper[1]), 3),
    CI_shape = round(c(lower[2], upper[2]),3))
  
  if (show_fit) {output$fit <- gpd_fit}
  
  return(output)
}

#' GPD log-likelihood with covariates
#' 
#' @author Conor Murphy
#' 
#' @param par A numeric vector of parameter values of length 3.
#' @param excess A numeric vector of excesses of some threshold.
#' @param thresh_par A numeric vector of threshold parameter values of length 2.
#' @param V A numeric vector of distance values corresponding to each excess.
#' @param ics A numeric vector of stress values corresponding to each excess.
#' 
#' @returns A numeric value of the log-likeihood.
GPD_LL_given_V_ICS <- function(par, excess, thresh_par, V, ics){
  
  # input checks
  stopifnot(length(par) == 3)
  stopifnot(length(thresh_par) == 2)
  if (!is.numeric(excess)) stop("excess must be a vector")
  if (!is.numeric(V)) stop("V must be vector")
  if (!is.numeric(ics)) stop("ics must be vector")
  stopifnot(length(excess) == length(V))
  stopifnot(length(excess) == length(ics))
  
  sigma_par <- par[1:2]
  xi <- par[3]
  
  sigma_ics <- sigma_par[1] + sigma_par[2]*ics
  sigma_tilde <- sigma_ics + xi*(thresh_par[[1]] + thresh_par[[2]]*V)
  sigma_check <- c(sigma_ics, sigma_tilde)
  
  if (any(sigma_check < 0) || xi <= -0.9 || sigma_par[2] < 0) {return(-1e7)}
  
  if (abs(xi) < 1e-10) {return(-sum(log(sigma_tilde)) - sum(excess/sigma_tilde))}
  
  if (any(1 + (xi * excess) / sigma_tilde <= 0)) {return(-1e6)}
  
  LL1 <- -sum(log(sigma_tilde))
  LL2 <- -(1 + 1 / xi) * (sum(log(1 + (xi * excess) / sigma_tilde)))
  return(LL1 + LL2)
  
}

#' Obtain MLEs, SEs and CIs for parameters from a GPD with covariates
#' 
#' @author Conor Murphy
#' 
#' @param mags A numeric vector of magnitudes.
#' @param thresh_fit A list containing the threshold parameters.
#' @param V A numeric vector of distance values corresponding to each magnitude.
#' @param ics A numeric vector of stress values corresponding to each magnitude.
#' @param show_fit Logical. If TRUE, the full fit object is returned.
#' 
#' @returns A list with MLEs, SEs and CIs for the parameters of the 
#' covariate dependent GPD model.
get_par_ests_geo_ics <- function(mags, thresh_fit, V, ics, show_fit = TRUE){
  
  threshold <- thresh_fit$thresh_par[1] + thresh_fit$thresh_par[2]*V
  excesses <- mags[mags > threshold] - threshold[mags > threshold]
  V_excess <- V[mags > threshold]
  ics_excess <- ics[mags > threshold]
  
  mle0 <- mean(excesses)
  gpd_fit <- optim(GPD_LL_given_V_ICS,
                   excess = excesses,
                   par = c(mle0, 0, 0.1),
                   control = list(fnscale = -1),
                   thresh_par = thresh_fit$thresh_par, 
                   V = V_excess,
                   ics = ics_excess,
                   hessian = TRUE)
  
  check <- gpd_fit$par[1] == thresh_fit$par[1] &
    gpd_fit$par[2] == thresh_fit$par[2] &
    gpd_fit$par[3] == thresh_fit$par[3]
  if (!check) { print(gpd_fit$par); print("Parameter estimates don't agree") }
  
  hessian <- gpd_fit$hessian
  cov_matrix <- solve(-1*hessian)
  se <- sqrt(diag(cov_matrix))
  z <- qnorm(0.975)
  lower <- gpd_fit$par - z*se
  upper <- gpd_fit$par + z*se
  
  output <- list(
    par = gpd_fit$par,
    SE = se,
    CI_beta0 = round(c(lower[1], upper[1]), 3),
    CI_beta1 = round(c(lower[2], upper[2]), 3), 
    CI_shape = round(c(lower[3], upper[3]), 3))
  
  if (show_fit) {output$fit <- gpd_fit}
  
  return(output)
}

#' Generate Q-Q plot on Exponential margins for GPD with covariates
#' 
#' @author Conor Murphy
#' 
#' @param mags A numeric vector of magnitudes.
#' @param thresh_fit A list containing the threshold parameters.
#' @param V A numeric vector of distance values corresponding to each magnitude.
#' @param ics A numeric vector of stress values corresponding to each magnitude.
#' @param n_boot Number of bootstrap samples to use for tolerance intervals.
#' @param main Title for the plot.
#' @param xlim x-axis limits for the plot.
#' @param ylim y-axis limits for the plot.
#' 
#' @returns A Q-Q plot with 95% tolerance intervals.
get_qq_plot_geo_ics <- function(mags, thresh_fit, V, ics, n_boot = 200, main = "Q-Q plot", xlim=c(0,6), ylim=c(0,6)){
  
  threshold <- thresh_fit$thresh_par[1] + thresh_fit$thresh_par[2] * V
  excesses <- mags[mags > threshold] - threshold[mags > threshold]
  ics_excess <- ics[mags > threshold]
  sigma_tilde <- thresh_fit$par[1] + 
    thresh_fit$par[2] * ics_excess + 
    thresh_fit$par[3] * threshold[mags > threshold]
  
  transformed_excess <- transform_to_exp(y = excesses,
                                         sig = sigma_tilde,
                                         xi = thresh_fit$par[3])
  probs <- seq_along(excesses)/(length(excesses) + 1)
  sample_quantiles <- quantile(transformed_excess, probs = probs)
  model_quantiles <- qexp(probs, rate = 1)
  
  bootstrapped_quantiles <- matrix(NA, nrow = n_boot, ncol = length(probs))
  for (i in 1:n_boot) {
    excess_boot <- rexp(length(excesses), rate = 1)
    bootstrapped_quantiles[i,] <- quantile(excess_boot, probs = probs)
  }
  upper <- apply(bootstrapped_quantiles, 2, quantile, prob = 0.975)
  lower <- apply(bootstrapped_quantiles, 2, quantile, prob = 0.025)
  plot(
    x = model_quantiles,
    y = sample_quantiles,
    xlab = "Theoretical quantiles",
    ylab = "Sample quantiles",
    pch = 19,
    asp = 1,
    main = main,
    xlim = xlim,
    ylim = ylim)
  abline(a = 0, b = 1, col = "grey")
  lines(model_quantiles, upper, col = "red", lwd = 2, lty = "dashed")
  lines(model_quantiles, lower, col = "red", lwd = 2, lty = "dashed")
}


#' Generate Q-Q plot on Exponential margins for standard GPD
#' 
#' @author Conor Murphy
#' 
#' @param mags A numeric vector of magnitudes.
#' @param threshold A numeric value of the threshold.
#' @param n_boot Number of bootstrap samples to use for tolerance intervals.
#' @param main Title for the plot.
#' @param xlim x-axis limits for the plot.
#' @param ylim y-axis limits for the plot.
#' 
#' @returns A Q-Q plot with 95% tolerance intervals.
get_qq_plot_const <- function(mags, threshold, n_boot = 200, main = "Q-Q plot", xlim=c(0,6), ylim=c(0,6)){
  
  excesses <- mags[mags > threshold] - threshold
  par_ests <- get_par_ests(mags = mags, threshold = threshold)$par
  transformed_excess <- transform_to_exp(y = excesses, sig = par_ests[1], xi = par_ests[2])
  
  probs <- seq_along(excesses)/(length(excesses) + 1)
  sample_quantiles <- quantile(transformed_excess, probs = probs)
  model_quantiles <- qexp(probs, rate = 1)
  
  bootstrapped_quantiles <- matrix(NA, nrow = n_boot, ncol = length(probs))
  for (i in 1:n_boot) {
    excess_boot <- rexp(length(excesses), rate = 1)
    bootstrapped_quantiles[i,] <- quantile(excess_boot, probs = probs)
  }
  
  upper <- apply(bootstrapped_quantiles, 2, quantile, prob = 0.975)
  lower <- apply(bootstrapped_quantiles, 2, quantile, prob = 0.025)
  plot(
    x = model_quantiles,
    y = sample_quantiles,
    xlab = "Theoretical quantiles",
    ylab = "Sample quantiles",
    pch = 19,
    asp = 1,
    main = main,
    xlim = xlim,
    ylim = ylim)
  abline(a = 0, b = 1, col = "grey")
  lines(model_quantiles, upper, col = "red", lwd = 2, lty = "dashed")
  lines(model_quantiles, lower, col = "red", lwd = 2, lty = "dashed")
}

#' Generate P-P plot for standard GPD model
#' 
#' @author Conor Murphy
#' 
#' @param mags A numeric vector of magnitudes.
#' @param threshold A numeric value of the threshold.
#' @param n_boot Number of bootstrap samples to use for tolerance intervals.
#' @param main Title for the plot.
#' @param xlim x-axis limits for the plot.
#' @param ylim y-axis limits for the plot.
#' 
#' @returns A P-P plot with 95% tolerance intervals.
get_pp_plot_const <- function(mags, threshold, n_boot = 200, main = "P-P plot", xlim=c(0,1), ylim=c(0,1)){
  
  excesses <- mags[mags > threshold] - threshold
  par_ests <- get_par_ests(mags = mags, threshold = threshold)$par
  transformed_excess <- transform_to_exp(y = excesses, sig = par_ests[1], xi = par_ests[2])
  
  transformed_excess_sorted <- sort(transformed_excess)
  sample_probs <- seq_along(transformed_excess) / (length(transformed_excess) + 1)
  model_probs <- pexp(transformed_excess_sorted, rate = 1)
  
  bootstrapped_probs <- matrix(NA, nrow = n_boot, ncol = length(excesses))
  
  for (i in 1:n_boot) {bootstrapped_probs[i,] <- sort(runif(length(excesses)))}
  
  upper <- apply(bootstrapped_probs, 2, quantile, prob = 0.975)
  lower <- apply(bootstrapped_probs, 2, quantile, prob = 0.025)
  
  plot(
    x = sample_probs, 
    y = model_probs, 
    xlab = "Sample probabilities", 
    ylab = "Model probabilities", 
    pch = 19,
    asp = 1,
    main = main,
    xlim = xlim,
    ylim = ylim)
  abline(a = 0, b = 1, col = "grey")
  lines(sample_probs, upper, col = "red", lwd = 2)
  lines(sample_probs, lower, col = "red", lwd = 2)
}

#' Generate P-P plot for GPD with covariates
#' 
#' @author Conor Murphy
#' 
#' @param mags A numeric vector of magnitudes.
#' @param thresh_fit A list containing the threshold parameters.
#' @param V A numeric vector of distance values corresponding to each magnitude.
#' @param ics A numeric vector of stress values corresponding to each magnitude.
#' @param n_boot Number of bootstrap samples to use for tolerance intervals.
#' @param main Title for the plot.
#' @param xlim x-axis limits for the plot.
#' @param ylim y-axis limits for the plot.
#' 
#' @returns A P-P plot with 95% tolerance intervals.
get_pp_plot_geo_ics <- function(mags, thresh_fit, V, ics, n_boot = 200, main = "P-P plot", xlim=c(0,1), ylim=c(0,1)){

  threshold <- thresh_fit$thresh_par[1] + thresh_fit$thresh_par[2]*V
  excesses <- mags[mags > threshold] - threshold[mags > threshold]
  ics_excess <- ics[mags > threshold]
  sigma_tilde <- thresh_fit$par[1] + 
    thresh_fit$par[2] * ics_excess +
    thresh_fit$par[3] * threshold[mags > threshold]
  
  transformed_excess <- transform_to_exp(y = excesses,
                                         sig = sigma_tilde,
                                         xi = thresh_fit$par[3])
  transformed_excess_sorted <- sort(transformed_excess)
  
  sample_probs <- seq_along(transformed_excess) / (length(transformed_excess) + 1)
  model_probs <- pexp(transformed_excess_sorted, rate = 1)
  
  bootstrapped_probs <- matrix(NA, nrow = n_boot, ncol = length(excesses))
  
  for (i in 1:n_boot) {bootstrapped_probs[i,] <- sort(runif(length(excesses)))}
  
  upper <- apply(bootstrapped_probs, 2, quantile, prob = 0.975)
  lower <- apply(bootstrapped_probs, 2, quantile, prob = 0.025)
  
  plot(
    x = sample_probs, 
    y = model_probs, 
    xlab = "Sample probabilities", 
    ylab = "Model probabilities", 
    pch = 19,
    asp = 1,
    main = main,
    xlim = xlim,
    ylim = ylim)
  abline(a = 0, b = 1, col = "grey")
  lines(sample_probs, upper, col = "red", lwd = 2)
  lines(sample_probs, lower, col = "red", lwd = 2)
}

