#----------Generating EQ locations and calculating V function----------------
n <- 1000
eq_loc <- runif(n,0,10)
eq_mags <- rgpd(length(eq_loc), shape=0.1, scale=1)
geo <- runif(20,0,10)

third_dist <- function(x){
  min_ind <- which.min(x)
  x_2 <- x[-min_ind]
  min_ind_2 <- which.min(x_2)
  x_3 <- x_2[-min_ind_2]
  min_ind_3 <- which.min(x_3)
  return(x_3[min_ind_3])
}


V <- function(eq_loc, geo){
  dist <- sapply(eq_loc, function(x){sapply(geo, function(z){(x-z)^2})})
  v <- apply(dist, 2, third_dist) 
  return(v)
}

plot(eq_loc,V(eq_loc, geo))
points(geo,rep(0,length(geo)), col="red", pch=4)

#-----------Generating EQ mags and censoring----------------------------
setwd("C:/Users/murphyc4/OneDrive - Lancaster University") #Windows
source('STOR-i/PhD/R/threshold_paper_code-main/00_src/_gpd.R')
library(evir)

plot(eq_loc,eq_mags)

#Constant threshold
cons <- 1.0
abline(h=cons, col="red")
abline(v=geo, col="blue")

#Probabilistic censoring with constant threshold
cens_thr<-cons*rbeta(length(eq_mags),1,1)
keep <- eq_mags>cens_thr
u<-cons*rbeta(length(eq_mags),1,1)
keep <- eq_mags>u
plot(EQ$loc,EQ$mags)
points(EQ$loc[keep], EQ$mags[keep], pch=20)
abline(h=cons)

EQ <- data.frame(mags=eq_mags,loc=eq_loc,thr=cens_thr,keep=keep)
hist(EQ$mags[keep], breaks=30, freq = FALSE)
lines(density(EQ$mags[keep]))
abline(v=cons, lwd=2, col="red")


#Probabilistic censoring with a varying threshold
dist <- V(eq_loc,geo)
plot(eq_loc,V1)
theta <- 1.0
cens_thr <- theta*dist
u<-cens_thr*rbeta(length(eq_mags),1,1)
keep <- eq_mags>u
EQ <- data.frame(mags=eq_mags,loc=eq_loc,thr=cens_thr,keep=keep, V=dist)
plot(EQ$loc,EQ$mags)
#points(eq_loc,u , col="green")
points(EQ$loc,cens_thr, col="red")
points(EQ$loc[EQ$keep], EQ$mags[EQ$keep], pch=20)
abline(v=geo, col="lightblue")


hist(EQ$mags[EQ$keep], breaks=30, freq = FALSE)
lines(density(EQ$mags[EQ$keep]))
abline(v=theta, lwd=2, col="red")

EQ_sorted <- EQ[order(EQ$loc),]
EQ_kept <- subset(EQ_sorted, EQ_sorted$keep)

#Threshold selection
library(dplyr)
GPD_LL_adj<-function(par,z,V){
  xi<-par[1]
  sigma<-par[2]
  theta <- par[3]
  sigma_tilde <- sigma + xi*theta*V
  if(theta <= 0 || sigma <=0){
    return(-1e6)
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
    return(-1e6)
  }
}
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
      # gpddata <- gpd(data,thresh, method = "pwm")
      # par1 <- gpddata$par.ests[[1]]
      # par2 <- gpddata$par.ests[[2]]
      #print(c("pars=",par1,par2))
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
      #print(c("init=", xi,sigma))
      distances <- NULL
      for(j in 1:k){
        #X <- rgpd(n=m,shape = xi,scale = sigma)
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
        #print(c("fit=",gpd.fit$par[[1]],gpd.fit$par[[2]]))
        # xis <- c(xis, gpd.fit$par[[1]])
        # sigmas <- c(sigmas, gpd.fit$par[[2]])
        std_exp <- standardise_gpd_ns_sample(sample.df$mags, xi= gpd.fit$par[1], sig= gpd.fit$par[2], theta=gpd.fit$par[3], V=sample.df$V)
        quants<-qexp(c(1:mm)/(mm+1),rate=1)
        #X_sorted <- sort(X)
        #plot(quants, X_sorted)
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
  #result <- data.frame(thr=choice,xi=xi,sigma=sigma,n=len)
  result <- c(choice, xi, sigma, len)
  return(result)
}

thresholds <- seq(0.4, 3, by=0.05)
spat_distmetric(df=EQ_kept, thresholds = thresholds)
