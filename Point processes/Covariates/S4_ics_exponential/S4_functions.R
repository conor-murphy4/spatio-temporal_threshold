## functions for Poisson_Topograd.R

# Covariate <- function(CovDF){
#
#   Grd<- CovDF[,1:2]
#   Cov <- CovDF[,3]
#
#   out <- list(Grd,Cov)
#   names(out) <- c('Grd','Cov')
#   return(out)
# }

Catalogue<- function(EC, minMag = 1.5, minDate = '1995-01-01' , maxDate = '2014-12-31'){
  # Format EC columnns as dates
  EC$Date <- dmy(EC$Date)
  EC$Time_UT <- hms(EC$Time_UT)

  #Trim to range
  EC <- EC[EC$Magnitude > minMag,]
  EC <- EC[EC$Date >= ymd(minDate),]
  EC <- EC[EC$Date <= ymd(maxDate),]

  #List detailing EC properties

  x <- EC$Easting
  y <- EC$Northing
  time <- EC$Time_UT
  date <- EC$Date
  mag  <- EC$Magnitude
  ld   <- date + time
  nEQ  <- nrow(EC)
  Mc   <- minMag
  minDate <- ymd(minDate)
  maxDate <- ymd(maxDate)
  Year    <- year(EC$Date)
  nYr     <- max(Year) - min(Year)
  TG      <- c()

  Cat <- list(x, y, time, date, mag, ld,
              nEQ, Mc, minDate, maxDate, Year, nYr,
              TG)
  names(Cat) <- c('x', 'y', 'time', 'date', 'mag',
                  'ld', 'nEQ', 'Mc', 'minDate','maxDate',
                  'Year', 'nYr', 'TG')
  # output
  return(Cat)
}


cov_lookup <- function(x, y, t, cov_mats){

  tFloor <- floor(t)
  xRound <- round(x/500,digits = 0)*500
  yRound <- round(y/500,digits = 0)*500
  locRound <- cbind(xRound,yRound)

  lookups_early <- cbind(
    match(xRound, cov_mats$xgrid),
    match(yRound, cov_mats$ygrid),
    match(tFloor, cov_mats$tgrid))

  lookups_late <- cbind(
    match(xRound, cov_mats$xgrid),
    match(yRound, cov_mats$ygrid),
    match(tFloor + 1, cov_mats$tgrid))

  vals_early <- cov_mats[[1]][lookups_early]
  vals_late <- cov_mats[[1]][lookups_late]

  vals <- vals_early + (vals_late - vals_early) * ( t %% 1 )
  out <- vals
  return(out)
}

cov_rate_lookup <- function(x, y, t, cov_mats){

  tFloor <- floor(t)
  xRound <- round(x/500,digits = 0)*500
  yRound <- round(y/500,digits = 0)*500
  locRound <- cbind(xRound,yRound)

  lookups_early <- cbind(
    match(xRound, cov_mats$xgrid),
    match(yRound, cov_mats$ygrid),
    match(tFloor, cov_mats$tgrid))

  lookups_late <- cbind(
    match(xRound, cov_mats$xgrid),
    match(yRound, cov_mats$ygrid),
    match(tFloor + 1, cov_mats$tgrid))

  vals_early <- cov_mats[[1]][lookups_early]
  vals_late  <- cov_mats[[1]][lookups_late]
  out <- vals_late - vals_early
  return(out)
}
# example <- smooth_cov_lookup(x = Cat$x, y = Cat$y, smoothed_cov_mats = test)

#
# llik <- function(params, EC){
#   if(params[1]<=0 ){
#     out <- -100000
#   }
#   else{
#     n <- EC$nEQ
#     out <- -params[1] * sum(TG$Cov * inPolyVec) * EC$nYr * 0.5^2
#     out <- out + n*log(params[1])  + sum(log(EC$TG))
#   }
#   return(out)
#
# }
#
#
llik_S4 <- function(params, catalogue, cov_mats, inpolymat = inPolyMat, negative = FALSE){
  #Check parameters are valid
  beta_0 <- params[1]
  beta_1 <- params[2]
  sigma <- params[3]
  beta_1_lower_limit <- 0 # -1 / max(cov_mats[[1]], na.rm = TRUE)

  if(beta_0 <= 0 | beta_1 <= beta_1_lower_limit | sigma <= 0){

    llh <- -1e6

  } else {

    # apply smoothing to covariate rate array - :checkmark:
    smoothed_cov_mats <- smooth_cov_mats(
      cov_mats = cov_mats,
      lambda =  sigma,
      inpolymat =  inpolymat)

    # look up covariate values at each event - :checkmark:
    cov_at_eqs <- cov_lookup(
      x = catalogue$x,
      y = catalogue$y,
      t = decimal_date(catalogue$ld),
      cov_mats = smoothed_cov_mats)

    # look up covariate rate at each event - :checkmark:
    cov_rate_at_eqs <- cov_rate_lookup(
      x = catalogue$x,
      y = catalogue$y,
      t = decimal_date(catalogue$ld),
      cov_mats = smoothed_cov_mats)

    # Calculate log-likelihood
    n_slices <- dim(smoothed_cov_mats[[1]])[3]
    n_eqs <- catalogue$nEQ
    final_slice <- smoothed_cov_mats[[1]][,,n_slices]
    initial_slice <- smoothed_cov_mats[[1]][,,1]

    Lambda_mat <- (beta_1)^(-1) * beta_0 * (0.25) * as.numeric(inpolymat) * (exp(beta_1 * final_slice) - exp(beta_1 * initial_slice))

    llh <- - sum(Lambda_mat, na.rm = TRUE)
    llh <- llh + n_eqs * log(beta_0)  + sum(log(cov_rate_at_eqs), na.rm = TRUE) + sum(beta_1 * cov_at_eqs, na.rm = TRUE)
  }

  if(negative) llh <- llh * (-1)
  return(llh)
}
# example <- llik_smooth(params = c(1,800), catalogue = Cat, cov_mats = tgmats)

#auxihiary likelhiood for use in profile llh
aux_llik_S4 <- function(sigma, params, catalogue, cov_mats, inpolymat = inPolyMat, neg = FALSE){
  #Check parameters are valid
  beta_0 <- params[1]
  beta_1 <- params[2]
  beta_1_lower_limit <- 0 #-1 / max(cov_mats[[1]], na.rm = TRUE)

  if(beta_0 <= 0 | beta_1 <= beta_1_lower_limit | sigma <= 0){

    llh <- -1e6

  } else {

    # apply smoothing to covariate rate array - :checkmark:
    smoothed_cov_mats <- smooth_cov_mats(
      cov_mats = cov_mats,
      lambda =  sigma,
      inpolymat =  inpolymat)

    # look up covariate values at each event - :checkmark:
    cov_at_eqs <- cov_lookup(
      x = catalogue$x,
      y = catalogue$y,
      t = decimal_date(catalogue$ld),
      cov_mats = smoothed_cov_mats)

    # look up covariate rate at each event - :checkmark:
    cov_rate_at_eqs <- cov_rate_lookup(
      x = catalogue$x,
      y = catalogue$y,
      t = decimal_date(catalogue$ld),
      cov_mats = smoothed_cov_mats)

    # Calculate log-likelihood
    n_slices <- dim(smoothed_cov_mats[[1]])[3]
    n_eqs <- catalogue$nEQ
    final_slice <- smoothed_cov_mats[[1]][,,n_slices]
    initial_slice <- smoothed_cov_mats[[1]][,,1]

    Lambda_mat <- (beta_1)^(-1) * beta_0 * (0.25) * as.numeric(inpolymat) * (exp(beta_1 * final_slice) - exp(beta_1 * initial_slice))

    llh <- - sum(Lambda_mat, na.rm = TRUE)
    llh <- llh + n_eqs * log(beta_0)  + sum(log(cov_rate_at_eqs), na.rm = TRUE) + sum(beta_1 * cov_at_eqs, na.rm = TRUE)
  }

  if(neg) llh <- llh * (-1)
  return(llh)

}

prof_llik_S4 <- function(params, catalogue, cov_mats, inpolymat = inPolyMat, negative = TRUE, detailed = FALSE, ...){
  opt <- optim(par = sigma,
               fn = aux_llik_S4,
               method = "Brent",
               lower = 0,
               upper = 40000,
               params = params,
               catalogue = catalogue,
               cov_mats = cov_mats,
               inpolymat = inpolymat,
               neg = TRUE,
               ... )
  if(detailed) return(opt)

  out <- ifelse(negative, opt$value, -opt$value)
  return(out)
}

smooth_cov_mat <- function(image, lambda, inpolymat){
  # image  = matrix of observed values
  # lambda = smoothing length scale (m)
  im <- spatstat.geom::as.im(X = image)
  smoothed_im <- spatstat.core::blur(
    x = im,
    sigma = lambda / 500,
    kernel = "gaussian",
    normalise = TRUE,
    bleed = FALSE
  )
  in_poly_else_NA <- inpolymat
  in_poly_else_NA[in_poly_else_NA == 0] <- NA
  return(smoothed_im$v * in_poly_else_NA)
}

smooth_cov_mats <- function(cov_mats, lambda = 1e-6, inpolymat = inPolyMat){
  smoothed_cov <- plyr::aaply(
    .data = cov_mats[[1]],
    .margins = 3,
    .fun = smooth_cov_mat,
    lambda = lambda,
    inpolymat = inpolymat)
  smoothed_cov <- aperm(unname(smoothed_cov), c(2,3,1))

  out <- list(
    smoothed_cov = smoothed_cov,
    xgrid = cov_mats$xgrid,
    ygrid = cov_mats$ygrid,
    tgrid = cov_mats$tgrid)
  return(out)
}
#test <- smooth_cov_mats(cov_mats = tgmats, lambda = 1, inpolymat = inPolyMat)

# make.tg_mats <- function(Y1_lim = c(1959,1959)){
#
#   nx   <- length(unique(tg$X))
#   minx <- min(tg$X)
#   maxx <- max(tg$X)
#   ny   <- length(unique(tg$Y))
#   miny <- min(tg$Y)
#   maxy <- max(tg$Y)
#
#   xgrid = seq(minx,maxx,length= nx ) # x values of observation points
#   ygrid = seq(miny,maxy,length= ny) # y values of observation points
#
#   year <- 1:64
#   tg_mats <- array(data = NA, dim = c(nx,ny,length(year)))
#
#   for (j in year){
#     tgs <- tg[[(Y1_lim[2] - 1959 + 2 + 1)]]
#
#     tg_mat <- matrix(NA, nrow = nx, ncol = ny)
#     for( i in 1: length(tg$X)){
#       x <- (tg$X[i] - minx)/500 + 1
#       y <- (tg$Y[i] - miny)/500 + 1
#       tg_mat[x,y] <- tgs[i]
#     }
#
#     tg_mats[,,j] <- tg_mat * inPolyMat
#   }
#
#   out <- list(tg_mats,xgrid,ygrid)
#   names(out) <- c('tg_mats','xgrid','ygrid')
#   return(out)
# }
#
box.no2xyt<- function(n,len.x = 84,len.y = 99, minX = NA, minY = NA, minYr =NA){
  #Input: n box number or vector of box numbers
  #       len.x number of pixels in x direction of image
  #       len.y number of pixels in y direction of image
  #Action: returns (x,y,t) centred spatially but at start of year
  #Output: vector of length 3. (x,y,t)
  a <-  (n %% (len.x*len.y))%%len.x
  b <-  ((n %% (len.x*len.y))-a)/len.x +1
  c <-  (n-a-(b-1)*len.x)/(len.x * len.y) +1

  x<- (a-1)*500 + minX
  y<- (b-1)*500 + minY
  t<- (minYr - 1) + c

  if( any(is.na(c(minX, minY, minYr)))){
    print("returning (a,b,c)")
    return(data.frame(a,b,c))
  }else{
    print("returning (x,y,t)")
    return(data.frame(x,y,t))
  }
}

subset_to_years_cov_mats<- function(cov_mats, year_min, year_max){
  slice_min <- which(cov_mats$tgrid == year_min)
  slice_max <- which(cov_mats$tgrid == year_max)

  out <- cov_mats
  out[[1]] <- out[[1]][,,slice_min:slice_max]
  out$tgrid <- out$tgrid[slice_min:slice_max]

  return(out)
}

#' Calculate the change in a covariate array over each epoch
#'
#' @param cov_mats A covariate array object
#' @param new_name (char) Name for array component of new covariate array object
#'
#' @return A covariate array object of the same dimensions as cov_mats. Note the last slice will be filled with NAs
temporal_difference_cov_mats <- function(cov_mats, new_name= NULL){
  if(is.null(new_name)){ new_name = paste0(names(cov_mats)[1],"_rate") }

  out <- cov_mats
  n_slices <- dim(out[[1]])[3]
  out[[1]] <- out[[1]][,,-1] - out[[1]][,,-n_slices]
  out[[1]] <- abind::abind(out[[1]], out[[1]][,,1]*NA, along = 3)
  out
}

replace_negatives_cov_mats <- function(cov_mats, replace,replacement_value){
  n_to_change <- sum(cov_mats[[1]] <= 0, na.rm = TRUE)
  if(replace){
    cov_mats[[1]][cov_mats[[1]] <= 0] <- replacement_value
  }
  print(paste(n_to_change, "values changed."))
  cov_mats
}
