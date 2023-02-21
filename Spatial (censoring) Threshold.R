setwd("C:/Users/murphyc4/OneDrive - Lancaster University") 
source('STOR-i/PhD/R/threshold_paper_code-main/00_src/_gpd.R')


#---------Estimating parameters with censored data-----------------------

GPD_LL_adj<-function(par,z,V){
  xi<-par[1]
  sigma<-par[2]
  theta <- par[3]
  sigma_tilde <- sigma + xi*theta*V
  if(theta <= 0 || sigma <=0){
    return(-1e8)
  }
  if(all(sigma_tilde>0)){
    if(xi==0){
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

#------------------Threshold Selection------------------------------------

#Exponential transformation function
standardise_gpd_ns_sample <- function (y, xi, sig, theta, V){
  std_exp <- c()
  for(i in 1:length(y)){
    std_exp[i] <- (1 / xi) * log( 1 + xi * (y[i]/(sig+xi*theta*V[i])))  
  }
  return(std_exp)
}

#NOTE: Make function take whole data frame so locations are included so when sampling excesses, locations and V at locations are available for tranforming to exponential

#Adjusted Selection technique
#df needs to contain EQmags, loc and V sorted by loc
spat_distmetric <- function(df,k=100,thresholds,mm=500){
  meandistances <- NULL
  xis <- NULL
  sigmas <- NULL
  thetas <- NULL
  lens <- NULL
  meandists <<- NULL
  for(i in 1:length(thresholds)){
    df$thresh <- thresholds[i]*df$V
    m <- nrow(subset(df,df$mags>df$thresh))
    lens <- c(lens,m)
    if(m>10){
      subdf <- subset(df, df$mags>df$thresh)
      subdf$mags <- subdf$mags - subdf$thresh
      mle0 <- mean(subdf$mags)
      init.fit <- optim(GPD_LL_adj,z=subdf$mags,V=subdf$V,par=c(0.1,mle0,thresholds[i]),control=list(fnscale=-1)) 
      xi <- init.fit$par[[1]]
      sigma <- init.fit$par[[2]]
      theta <- init.fit$par[[3]]
      xis <- c(xis,xi)
      sigmas <- c(sigmas,sigma)
      thetas <- c(thetas,theta)
      distances <- NULL
      for(j in 1:k){
        sample.df <- sample_n(subdf, nrow(subdf), replace=TRUE)
        if(xi < 0){
          xi0 <- 0.1
          sigma0 <-  mean(sample.df$mags)
        }
        else{
          xi0 <- xi
          sigma0 <- sigma
        }
        gpd.fit<-optim(GPD_LL_adj,z=sample.df$mags,V=sample.df$V,par=c(xi0,sigma0,theta),control=list(fnscale=-1))
        std_exp <- standardise_gpd_ns_sample(sample.df$mags, xi= gpd.fit$par[1], sig= gpd.fit$par[2], theta=gpd.fit$par[3], V=sample.df$V)
        quants<-qexp(c(1:mm)/(mm+1),rate=1)
        dist <- (1/mm)*sum(abs(quantile(std_exp,probs=c(1:mm)/(mm+1))-quants))
        distances <- c(distances,dist)
      }
      meandist <- mean(distances)
      meandistances <- c(meandistances,meandist)
    }
    else{
      meandistances <- c(meandistances,NA)
    }
  }
  choice <- thresholds[which.min(meandistances)]
  xi <- xis[which.min(meandistances)]
  sigma <- sigmas[which.min(meandistances)]
  len <- lens[which.min(meandistances)]
  meandists <<- meandistances
  result <- c(choice, xi, sigma, len)
  return(result)
}

EQ.df <- subset(EQ_sorted, EQ_sorted$keep==1)
spat_distmetric(df=EQ.df,thresholds = seq(0.5,3.0, by=0.1))



