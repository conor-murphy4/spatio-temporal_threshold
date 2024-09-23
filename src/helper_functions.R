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

#rgpd_rd <- function(n, sig, xi, mu, to_nearest, mu_latent = NULL){
#  if(is.null(mu_latent)) mu_latent = mu - 0.5 * to_nearest
#  x <- rgpd(n = n, shape = xi, scale = sig, mu = mu_latent)
#  x <- round_to_nearest(x, to_nearest)
#  return(x)
#}

#' evalue probability mass function of rounded generalised Pareto distribution
#'
#' @author Zak Varty
#'
#' @param x Vector values at which to evaluate mass function
#' @param u Vector of latent threshold values
#' @param sig_u Vector of latent scale parameters (for exceedances of u)
#' @param xi Latent shape parameter
#' @param to_nearest Level of rounding
#'
#' @return pmf evaluated at x. NOTE: does not check validity of x values.
#'
#' @examples
#' dgpd_rd(x = seq(0.1, 1.5, by = 0.1), to_nearest = 0.1, u = 0, sig_u = 1, xi = 0)
#' dgpd_rd(x = seq(0.1, 1.5, by = 0.1), to_nearest = 0.1, u = seq(0, 1.4, by = 0.1), sig_u = 1, xi = 0)
#' # CAUTION:
#' gpd_rd(x = 0.15, to_nearest = 0.1,  u = 0, sig_u = 1, xi = 0)
dgpd_rd <- function(x, u, sig_u, xi, to_nearest){

  # If (Y_i - u_i | Y_i > u_i) ~ GPD(sig_i, xi)
  # then Z_i  = [(Y_i - u_i)/sig_i | Y_i - u_i > 0] ~ GPD(1, xi)

  # range of z values that lead to observing x
  x_low <- pmax(x - to_nearest/2, u)
  z_low <- (x_low - u) / sig_u
  z_high <- (x - u + to_nearest/2)/sig_u

  # calculate probability of z in that range
  p_high <- pgpd(q = z_high, scale = 1, shape = xi, mu = 0)
  p_low  <- pgpd(q = z_low,  scale = 1, shape = xi, mu = 0)
  p <- p_high - p_low

  return(p)
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


GPD_LL <- function(par, z){
  sigma <- par[1]
  xi <- par[2]
  if (sigma > 0) {
    if (abs(xi) < 1e-10) {
      return(-length(z) * log(sigma) - ((1 / sigma) * sum(z)))
    }
    else {
      if (all(1 + (xi * z) / sigma > 0)) {
        return(-(length(z) * log(sigma)) - ((1 + 1 / xi)*(sum(log(1 + (xi * z) / sigma)))))
      }
      else{
        return(-1e6)
      }
    }
  }
  else{
    return(-1e7)
  }
}


transform_to_exp <- function (y, sig, xi){
  std_exp <- (1 / xi) * log( 1 + xi * (y/sig))  
  return(std_exp)
}

GPD_LL_given_third_nearest <- function(par, excess, thresh_par, third_nearest, min_dist = 0, max_dist = 33.782){
  #Current min and max distances are conservative and should be updated based on Stephen's guidance
  
  if(length(par)!=2) stop("par must be a vector of length 2")
  if(length(thresh_par)!=2) stop("thresh must be a vector of length 2")
  if (!is.numeric(excess)) stop("excess must be a vector")
  if (!is.numeric(third_nearest)) stop("third_nearest must be vector")
  if(length(excess) != length(third_nearest)) stop("excess and third_nearest must be the same length")
  if(length(min_dist)!=1) stop("min_dist must be a scalar")
  if(length(max_dist)!=1) stop("max_dist must be a scalar")
  
  sigma<-par[1]
  xi<-par[2]
  
  sigma_tilde <- sigma + xi*(thresh_par[[1]] + thresh_par[[2]]*third_nearest)
  
  #Extra conditions to check for spatio-temporal model
  sigma_max <- sigma + xi*(thresh_par[[1]] + thresh_par[[2]]*max_dist)
  sigma_min <- sigma + xi*(thresh_par[[1]] + thresh_par[[2]]*min_dist)
  
  sigma_check <- c(sigma, sigma_min, sigma_max, sigma_tilde)
  
  if(all(sigma_check>0)){
    if(abs(xi)<1e-10){
      return(-sum(log(sigma_tilde))-sum(excess/sigma_tilde))
    }
    else{
      if(all(1+(xi*excess)/sigma_tilde >0)){
        return(-sum(log(sigma_tilde))-(1+1/xi)*(sum(log(1+(xi*excess)/sigma_tilde))))
      }
      else{
        return(-1e6)
      }
    }
  }
  else{
    return(-1e7)
  }
}

GPD_LL_given_third_nearest_unconstrained <- function(par, excess, thresh_par, third_nearest){
  
  if(length(par)!=2) stop("par must be a vector of length 2")
  if(length(thresh_par)!=2) stop("thresh must be a vector of length 2")
  if (!is.numeric(excess)) stop("excess must be a vector")
  if (!is.numeric(third_nearest)) stop("third_nearest must be vector")
  if(length(excess) != length(third_nearest)) stop("excess and third_nearest must be the same length")
  #if(length(min_dist)!=1) stop("min_dist must be a scalar")
  #if(length(max_dist)!=1) stop("max_dist must be a scalar")
  
  sigma<-par[1]
  xi<-par[2]
  
  sigma_tilde <- sigma + xi*(thresh_par[[1]] + thresh_par[[2]]*third_nearest)
  
  sigma_check <- c(sigma, sigma_tilde)
  
  if(all(sigma_check>0)){
    if(abs(xi)<1e-10){
      return(-sum(log(sigma_tilde))-sum(excess/sigma_tilde))
    }
    else{
      if(all(1+(xi*excess)/sigma_tilde >0)){
        return(-sum(log(sigma_tilde))-(1+1/xi)*(sum(log(1+(xi*excess)/sigma_tilde))))
      }
      else{
        return(-1e6)
      }
    }
  }
  else{
    return(-1e7)
  }
}

GPD_LL_step <- function(par, excess, thresh){
  #Current min and max distances are conservative and should be updated based on Stephen's guidance
  
  if(length(par)!=2) stop("par must be a vector of length 2")
  if(length(thresh)!=length(excess)) stop("excess and thresh must be the same length")
  if (!is.numeric(excess)) stop("excess must be a vector")
 
  sigma<-par[1]
  xi<-par[2]
  
  sigma_tilde <- sigma + xi*(thresh)
  
  sigma_check <- c(sigma, sigma_tilde)
  
  if(all(sigma_check>0)){
    if(abs(xi)<1e-10){
      return(-sum(log(sigma_tilde))-sum(excess/sigma_tilde))
    }
    else{
      if(all(1+(xi*excess)/sigma_tilde >0)){
        return(-sum(log(sigma_tilde))-(1+1/xi)*(sum(log(1+(xi*excess)/sigma_tilde))))
      }
      else{
        return(-1e6)
      }
    }
  }
  else{
    return(-1e7)
  }
}

#Function to get par ests and confidence intervals
get_par_ests_geo <- function(mags, thresh_fit, third_nearest_dist){
  threshold <- thresh_fit$thresh[1] + thresh_fit$thresh[2]*third_nearest_dist
  excesses <- mags[mags > threshold] - threshold[mags > threshold]
  third_nearest_dist_excess <- third_nearest_dist[mags > threshold]
  mle0 <- mean(excesses)
  gpd_fit <- optim(GPD_LL_given_third_nearest, excess = excesses, par = c(mle0,0.1), control = list(fnscale = -1), thresh_par=thresh_fit$thresh, 
                   third_nearest = third_nearest_dist_excess, min_dist = min(third_nearest_dist), max_dist = max(third_nearest_dist), hessian = TRUE)
  check <- gpd_fit$par[1] == thresh_fit$par[1] & gpd_fit$par[2] == thresh_fit$par[2]
  if(!check){
    stop("Parameter estimates don't agree")
  }
  hessian <- gpd_fit$hessian
  cov_matrix <- solve(-1*hessian)
  se <- sqrt(diag(cov_matrix))
  z <- qnorm(0.975)
  lower <- gpd_fit$par - z*se
  upper <- gpd_fit$par + z*se
  return(list(par=gpd_fit$par, lower=lower, upper=upper))
}

get_par_ests_step <- function(mags, threshold){
  threshold_excess <- threshold[mags > threshold]
  excesses <- mags[mags > threshold] - threshold_excess
  mle0 <- mean(excesses)
  gpd_fit <- optim(GPD_LL_step, excess = excesses, par = c(mle0,0.1), control = list(fnscale = -1), thresh=threshold_excess, hessian = TRUE)
  hessian <- gpd_fit$hessian
  cov_matrix <- solve(-1*hessian)
  se <- sqrt(diag(cov_matrix))
  z <- qnorm(0.975)
  lower <- gpd_fit$par - z*se
  upper <- gpd_fit$par + z*se
  return(list(par=gpd_fit$par, lower=lower, upper=upper))
}

get_qq_plot <- function(mags, thresh_fit, third_nearest_dist, n_boot = 200){
  threshold <- thresh_fit$thresh[1] + thresh_fit$thresh[2]*third_nearest_dist
  excesses <- mags[mags > threshold] - threshold[mags > threshold]
  third_nearest_dist_excess <- third_nearest_dist[mags > threshold]
  sigma_tilde <- thresh_fit$par[1] + thresh_fit$par[2]*threshold[mags > threshold]
  transformed_excess <- transform_to_exp(y = excesses, sig = sigma_tilde, xi = thresh_fit$par[2])
  probs <- c(1:length(excesses))/(length(excesses)+1)
  sample_quantiles <- quantile(transformed_excess, probs = probs)
  model_quantiles <- qexp(probs, rate = 1)
  
  bootstrapped_quantiles <- matrix(NA, nrow = n_boot, ncol = length(probs))
  for(i in 1:n_boot){
    excess_boot <- rexp(length(excesses), rate = 1)
    bootstrapped_quantiles[i,] <- quantile(excess_boot, probs = probs)
  }
  upper <- apply(bootstrapped_quantiles, 2, quantile, prob = 0.975)
  lower <- apply(bootstrapped_quantiles, 2, quantile, prob = 0.025)
  plot(model_quantiles, sample_quantiles, xlab = "Theoretical quantiles", ylab = "Sample quantiles", pch=19, asp=1)
  abline(a = 0, b = 1, col="grey")
  lines(model_quantiles, upper, col = "red", lwd=2, lty="dashed")
  lines(model_quantiles, lower, col = "red", lwd=2, lty="dashed")
}
