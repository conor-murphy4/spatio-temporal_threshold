#source("_gpd.R")

library(dplyr)

#-------------Spatio-temporal distance functions-------------------------

sq_eu_dist<- function(Eq_loc, GeoCoord){
  dist <- apply(Eq_loc, 1, function(x,z){(x[[1]]/1000-z[[1]]/1000)^2+(x[[2]]/1000-z[[2]]/1000)^2}, z=GeoCoord)
  return(dist)
}

Geo_current <- function(Geo, year){
  ind <- which(Geo$Start_Year <= year & Geo$End_Year >= year)
  Geo_current <- Geo[ind,]
  return(Geo_current)
}

V_2d_time <- function(Eq.df, GeoCoord){
  V <- numeric(length(Eq.df$X))
  for(i in 1:length(Eq.df$X)){
    Eq_current <- Eq.df[i,]
    Current_coords <- Geo_current(GeoCoord,Eq_current$Year)
    Eq_loc_current <- data.frame(X=Eq_current$X, Y=Eq_current$Y)
    dist <- sq_eu_dist(Eq_loc_current, Current_coords)
    V[i] <- sort(dist)[3]
  }
  return(V)
}

#---------Estimating parameters with censored data-----------------------

Spat_temp_GPD_LL<-function(par,z,V){
  sigma<-par[1]
  xi<-par[2]
  theta <- par[3]
  sigma_tilde <- sigma + xi*theta*V
  if(theta <= 0 || sigma <=0){
    return(-1e8)
  }
  else{
    if(all(sigma_tilde>0)){
      if(abs(xi)<1e-10){
        return(-sum(log(sigma_tilde))-sum(z/sigma_tilde))
      }
      else{
        if( all(1+(xi*z)/sigma_tilde >0)){
          return(-sum(log(sigma_tilde))-(1+1/xi)*(sum(log(1+(xi*z)/sigma_tilde))))
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
}

#------------------Threshold Selection------------------------------------

#Exponential transformation function

transform_to_exp <- function (y,sig, xi, theta, V){
  std_exp <- (1 / xi) * log( 1 + xi * (y/(sig+xi*theta*V)))  
  return(std_exp)
}

#Stationary GPD transformation

transform_to_gpd <- function(y, sig, xi, theta, V){
  std_gpd <- y/(sig + theta*xi*V)
  return(std_gpd)
}

#Spatio-temporal threshold selection technique
#With EXP transformation

#df needs to contain EQmags, loc and V sorted by loc
Spat_temp_thr_exp <- function(df,thresholds,k=100,m=500){
  # if(length(dim(data)) != 0){
  #   stop("Data must be a vector")
  # }
  # if(length(dim(thresholds)) != 0){
  #   stop("u to be tested needs to be a vector")
  # }
  # if(k <= 0 | k%%1 != 0){
  #   stop("Number of bootstrapped samples must be a positive integer")
  # }
  # if(m <= 0 | m%%1 != 0){
  #   stop("Number of equally spaced probabilities must be a positive integer")
  # }
  meandistances <- xis <- sigmas <- thetas <- num_excess <- numeric(length(thresholds))
  for(i in 1:length(thresholds)){
    df$thresh <- thresholds[i]*df$V
    subdf <- subset(df, df$M > df$thresh)
    subdf$excess <- subdf$M - subdf$thresh
    num_excess[i] <- nrow(subdf)
    if(num_excess[i]>10){
      mle0 <- mean(subdf$excess)
      init.fit <- optim(Spat_temp_GPD_LL,z=subdf$excess,par=c(mle0,0.1,1),V=subdf$V,control=list(fnscale=-1)) 
      sigmas[i] <- init.fit$par[[1]]
      xis[i] <- init.fit$par[[2]]
      thetas[i] <- init.fit$par[[3]]
      distances <- numeric(k)
      for(j in 1:k){
        sample.df <- sample_n(subdf, num_excess[i], replace=TRUE)
        mle <- mean(sample.df$excess)
        ifelse(xis[i] < 0, pars_init <-  c(mle, 0.1, 1) ,pars_init <- c(sigmas[i], xis[i], thetas[i]) )
        gpd.fit<-optim(Spat_temp_GPD_LL,z=sample.df$excess,V=sample.df$V,par=pars_init,control=list(fnscale=-1))
        std_exp <- transform_to_exp(sample.df$excess, sig= gpd.fit$par[1], xi= gpd.fit$par[2], theta=gpd.fit$par[3], V=sample.df$V)
        quants<-qexp(c(1:m)/(m+1),rate=1)
        distances[j] <- (1/m)*sum(abs(quantile(std_exp,probs=c(1:m)/(m+1))-quants))
      }
      meandistances[i] <- mean(distances)
    }
    else{
      meandistances[i] <- NA
    }
  }
  chosen_ind <- which.min(meandistances)
  choice <- thresholds[chosen_ind]
  xi <- xis[chosen_ind]
  sigma <- sigmas[chosen_ind]
  len <- num_excess[chosen_ind]
  result <- list(thresh=choice, par=c(sigma,xi), num=len, dists=meandistances)
  return(result)
}

#With GPD transformation
Spat_temp_thr_gpd <- function(df,thresholds,k=100,m=500){
  # if(length(dim(data)) != 0){
  #   stop("Data must be a vector")
  # }
  # if(length(dim(thresholds)) != 0){
  #   stop("u to be tested needs to be a vector")
  # }
  # if(k <= 0 | k%%1 != 0){
  #   stop("Number of bootstrapped samples must be a positive integer")
  # }
  # if(m <= 0 | m%%1 != 0){
  #   stop("Number of equally spaced probabilities must be a positive integer")
  # }
  meandistances <- xis <- sigmas <- thetas <- num_excess <- numeric(length(thresholds))
  for(i in 1:length(thresholds)){
    df$thresh <- thresholds[i]*df$V
    subdf <- subset(df, df$M > df$thresh)
    subdf$excess <- subdf$M - subdf$thresh
    num_excess[i] <- nrow(subdf)
    if(num_excess[i]>10){
      mle0 <- mean(subdf$excess)
      init.fit <- optim(Spat_temp_GPD_LL,z=subdf$excess,par=c(mle0,0.1,1),V=subdf$V,control=list(fnscale=-1)) 
      sigmas[i] <- init.fit$par[[1]]
      xis[i] <- init.fit$par[[2]]
      thetas[i] <- init.fit$par[[3]]
      distances <- numeric(k)
      for(j in 1:k){
        sample.df <- sample_n(subdf, num_excess[i], replace=TRUE)
        mle <- mean(sample.df$excess)
        ifelse(xis[i] < 0, pars_init <-  c(mle, 0.1, 1) ,pars_init <- c(sigmas[i], xis[i], thetas[i]) )
        gpd.fit<-optim(Spat_temp_GPD_LL,z=sample.df$excess,V=sample.df$V,par=pars_init,control=list(fnscale=-1))
        std_gpd <- transform_to_gpd(sample.df$excess, sig= gpd.fit$par[1], xi= gpd.fit$par[2],  theta=gpd.fit$par[3], V=sample.df$V)
        quants<-qgpd(c(1:m)/(m+1),scale = 1,shape=gpd.fit$par[[2]])
        distances[j] <- (1/m)*sum(abs(quantile(std_gpd,probs=c(1:m)/(m+1))-quants))
      }
      meandistances[i] <- mean(distances)
    }
    else{
      meandistances[i] <- NA
    }
  }
  chosen_ind <- which.min(meandistances)
  choice <- thresholds[chosen_ind]
  xi <- xis[chosen_ind]
  sigma <- sigmas[chosen_ind]
  len <- num_excess[chosen_ind]
  result <- list(thresh=choice, par=c(sigma,xi), num=len, dists=meandistances)
  return(result)
}