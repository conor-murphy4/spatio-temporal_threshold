rm(list=ls())
source("helper_functions.R")
source("thresh_qq_metric.R")
library(dplyr)
#Events
cat <- read.csv("Data/Events/2022-04-12_15-09-25_cat.csv")

#Groningen outline
gron_outline <- read.csv("Data/Geophones/Geophone data/Groningen_Field_outline.csv")

#Plotting
plot(cat$RD_X, cat$RD_Y, asp=1, xlab="Easting (m)", ylab = "Northing (m)")

points(gron_outline$X, gron_outline$Y, pch=19, col="blue")
y_div <- 735000 - 0.6*seq(220000, 280000, by=1)
points(seq(220000, 280000, by=1), y_div)

#Average mags across fields
library(pracma)
mean_grid <- function(start_x, end_x, start_y, end_y, cell_dim, cat){
  length_x <- (end_x-start_x)/cell_dim
  length_y <- (end_y-start_y)/cell_dim
  mean_mat <- matrix(NA, nrow = length_y, ncol = length_x)
  for(i in 1:length_x){
    for(j in 1:length_y){
      current_x <- start_x + (i-1)*cell_dim
      current_y <- start_y + (j-1)*cell_dim
      poly_x <- c(current_x, current_x + cell_dim,  current_x + cell_dim, current_x)
      poly_y <- c(current_y, current_y, current_y + cell_dim, current_y + cell_dim)
      current_loc <- inpolygon(cat$RD_X, cat$RD_Y, poly_x, poly_y)
      mean_mat[j,i] <- mean(cat$excess[current_loc])
    }
  }
  return(mean_mat)
}

cat$thresh <- 1.2
cat$thresh[cat$RD_Y <= (735000 - 0.6*cat$RD_X)] <- 0.876
extreme_whole <- cat %>% filter(cat$mag > cat$thresh)
extreme_whole$excess<- extreme_whole$mag - extreme_whole$thresh


mean_mat <- mean_grid(230000, 270000, 565000, 615000, 1000, extreme_whole)

x <- seq(230000,270000, by=1000)
y <- seq(565000, 615000, by=1000)
fields::image.plot(x,y,t(mean_mat), asp=1)
lines(gron_outline$X, gron_outline$Y, pch=19)
y_div <- 735000 - 0.6*seq(220000, 280000, by=1)
lines(seq(220000, 280000, by=1), y_div, lwd=2)

#Check
poly_x <- c(240000, 245000,  245000, 240000)
poly_y <- c(610000, 610000, 615000, 615000)
current_loc <- inpolygon(cat$RD_X, cat$RD_Y, poly_x, poly_y)
(mean_mat <- mean(cat$mag[current_loc]))


cat_upper <- cat[cat$RD_Y > (735000 - 0.6*cat$RD_X),]
cat_lower <- cat[cat$RD_Y <= (735000 - 0.6*cat$RD_X),]

cat_upper$index <- seq(length(cat_upper$mag), 1, by=-1)
cat_lower$index <- seq(length(cat_lower$mag), 1, by=-1)

par(mfrow=c(1,2))

plot(cat_lower$index,cat_lower$mag, ylab="Magnitude", xlab="Index", main="Lower region", ylim=c(0,3))
abline(h=1.07, col="red", lwd=2)
plot(cat_upper$index,cat_upper$mag, ylab = "Magnitude", xlab="Index", main = "Upper region", ylim=c(0,3))
abline(h=1.07, col="red", lwd=2)

#Threshold selection
thresh_whole <- quantile(cat$mag, seq(0,0.95, by=0.01))
thr_selection_whole <- thresh_qq_metric(cat$mag, thresh_whole, k=200)

thresh_upper <- quantile(cat_upper$mag, seq(0,0.95, by=0.01))
thr_selection_upper <- thresh_qq_metric(cat_upper$mag, thresh_upper, k=200)

thresh_lower <- quantile(cat_lower$mag, seq(0,0.95, by=0.01))
thr_selection_lower <- thresh_qq_metric(cat_lower$mag, thresh_lower, k=200)

#-------------------------------------- LRTs for variation in parameters across two regions----------------------------------------------
#LRT
##H0: Parameters constant across space
#H1: Varying scale with constant shape

#Constant threshold (1.07) with constant GPD across gas-field
extreme_whole <- cat$mag[cat$mag > 1.07] - 1.07
(gpd.fit_whole <- optim(GPD_LL, par=c(mean(extreme_whole), 0.1), z=extreme_whole, control = list(fnscale=-1)))
(loglik_constant <- gpd.fit_whole$value)

#Constant threshold (1.07) with GPD scale varying and constant shape across gas-field

extreme_whole <- cat %>% filter(cat$mag > 1.07)
extreme_whole$mag<- extreme_whole$mag - 1.07
extreme_upper <- cat_upper$mag[cat_upper$mag > 1.07] - 1.07
extreme_lower <- cat_lower$mag[cat_lower$mag > 1.07] - 1.07
gpd_whole_sigR <- optim(GPD_LL_sigR, par=c(mean(extreme_lower), mean(extreme_upper),0.1), cat=extreme_whole, control = list(fnscale=-1))
(loglik_sigR <- gpd_whole_sigR$value)

(test_stat <- -2*(loglik_constant - loglik_sigR))

(p.val <- pchisq(test_stat, df = 1, lower.tail = FALSE))

#Constant threshold (1.318) with constant GPD across gas-field
extreme_whole <- cat$mag[cat$mag > 1.318] - 1.318
(gpd.fit_whole <- optim(GPD_LL, par=c(mean(extreme_whole), 0.1), z=extreme_whole, control = list(fnscale=-1)))
(loglik_constant <- gpd.fit_whole$value)

#Constant threshold (1.318) with GPD scale varying and constant shape across gas-field

extreme_whole <- cat %>% filter(cat$mag > 1.318)
extreme_whole$mag<- extreme_whole$mag - 1.318
extreme_upper <- cat_upper$mag[cat_upper$mag > 1.318] - 1.318
extreme_lower <- cat_lower$mag[cat_lower$mag > 1.318] - 1.318
gpd_whole_sigR <- optim(GPD_LL_sigR, par=c(mean(extreme_lower), mean(extreme_upper),0.1), cat=extreme_whole, control = list(fnscale=-1))
(loglik_sigR <- gpd_whole_sigR$value)

(test_stat <- -2*(loglik_constant - loglik_sigR))

(p.val <- pchisq(test_stat, df = 1, lower.tail = FALSE))

#Stepped threshold (1.2, 0.876) with constant GPD across gas-field
extreme_upper <- cat_upper$mag[cat_upper$mag > 1.2] - 1.2
extreme_lower <- cat_lower$mag[cat_lower$mag > 0.876] - 0.876
extreme_whole <- c(extreme_upper, extreme_lower)
(gpd.fit_whole <- optim(GPD_LL, par=c(mean(extreme_whole), 0.1), z=extreme_whole, control = list(fnscale=-1)))
(loglik_constant <- gpd.fit_whole$value)


#Stepped threshold (1.2, 0.876) with GPD scale varying and constant shape across gas-field
cat$thr <- 1.2
cat$thr[cat$RD_Y <= (735000 - 0.6*cat$RD_X)] <- 0.876
extreme_whole <- cat %>% filter(cat$mag > cat$thr)
extreme_whole$mag<- extreme_whole$mag - extreme_whole$thr
extreme_upper <- cat_upper$mag[cat_upper$mag > 1.2] - 1.2
extreme_lower <- cat_lower$mag[cat_lower$mag > 0.876] - 0.876
gpd_whole_sigR <- optim(GPD_LL_sigR, par=c(mean(extreme_lower), mean(extreme_upper),0.1), cat=extreme_whole, control = list(fnscale=-1))
(loglik_sigR <- gpd_whole_sigR$value)


(test_stat <- -2*(loglik_constant - loglik_sigR))

(p.val <- pchisq(test_stat, df = 1, lower.tail = FALSE))


#LRT 
#H0: Parameters constant across space
#H1: Parameters significantly different in two regions

#Constant threshold (1.07) with constant GPD across gas-field
extreme_whole <- cat$mag[cat$mag > 1.07] - 1.07
(gpd.fit_whole <- optim(GPD_LL, par=c(mean(extreme_whole), 0.1), z=extreme_whole, control = list(fnscale=-1)))
(loglik_constant <- gpd.fit_whole$value)

#Constant threshold (1.07) with GPD parameters verying across regions

extreme_upper <- cat_upper$mag[cat_upper$mag > 1.07] - 1.07
(gpd_fit_upper <- optim(GPD_LL, par=c(mean(extreme_upper), 0.1), z=extreme_upper, control = list(fnscale=-1)))

extreme_lower <- cat_lower$mag[cat_lower$mag > 1.07] - 1.07
(gpd_fit_lower <- optim(GPD_LL, par=c(mean(extreme_lower), 0.1), z=extreme_lower, control = list(fnscale=-1)))

(loglik_spatial <- gpd_fit_lower$value + gpd_fit_upper$value)  

(test_stat <- -2*(loglik_constant - loglik_spatial))

(p.val <- pchisq(test_stat, df = 2, lower.tail = FALSE))

#Constant threshold (1.318) with constant GPD across gas-field
extreme_whole <- cat$mag[cat$mag > 1.318] - 1.318
(gpd.fit_whole <- optim(GPD_LL, par=c(mean(extreme_whole), 0.1), z=extreme_whole, control = list(fnscale=-1)))
(loglik_constant <- gpd.fit_whole$value)

#Constant threshold (1.318) with GPD parameters verying across regions

extreme_upper <- cat_upper$mag[cat_upper$mag > 1.318] - 1.318
(gpd_fit_upper <- optim(GPD_LL, par=c(mean(extreme_upper), 0.1), z=extreme_upper, control = list(fnscale=-1)))

extreme_lower <- cat_lower$mag[cat_lower$mag > 1.318] - 1.318
(gpd_fit_lower <- optim(GPD_LL, par=c(mean(extreme_lower), 0.1), z=extreme_lower, control = list(fnscale=-1)))

(loglik_spatial <- gpd_fit_lower$value + gpd_fit_upper$value)  

(test_stat <- -2*(loglik_constant - loglik_spatial))

(p.val <- pchisq(test_stat, df = 2, lower.tail = FALSE))

#Stepped threshold (1.2, 0.876) with constant GPD across gas-field
extreme_upper <- cat_upper$mag[cat_upper$mag > 1.2] - 1.2
extreme_lower <- cat_lower$mag[cat_lower$mag > 0.876] - 0.876
extreme_whole <- c(extreme_upper, extreme_lower)
(gpd.fit_whole <- optim(GPD_LL, par=c(mean(extreme_whole), 0.1), z=extreme_whole, control = list(fnscale=-1)))
(loglik_constant <- gpd.fit_whole$value)

#Stepped threshold (1.2, 0.876) with  GPD parameters verying across regions

extreme_upper <- cat_upper$mag[cat_upper$mag > 1.2] - 1.2
(gpd_fit_upper <- optim(GPD_LL, par=c(mean(extreme_upper), 0.1), z=extreme_upper, control = list(fnscale=-1)))

extreme_lower <- cat_lower$mag[cat_lower$mag > 0.876] - 0.876
(gpd_fit_lower <- optim(GPD_LL, par=c(mean(extreme_lower), 0.1), z=extreme_lower, control = list(fnscale=-1)))

(loglik_spatial <- gpd_fit_lower$value + gpd_fit_upper$value)  

(test_stat <- -2*(loglik_constant - loglik_spatial))

(p.val <- pchisq(test_stat, df = 2, lower.tail = FALSE))


#LRT
##H0: Varying scale with constant shape
#H1: Varying scale and shape

#Constant threshold (1.07) with GPD scale varying and constant shape across gas-field

extreme_whole <- cat %>% filter(cat$mag > 1.07)
extreme_whole$mag<- extreme_whole$mag - 1.07
extreme_upper <- cat_upper$mag[cat_upper$mag > 1.07] - 1.07
extreme_lower <- cat_lower$mag[cat_lower$mag > 1.07] - 1.07
gpd_whole_sigR <- optim(GPD_LL_sigR, par=c(mean(extreme_lower), mean(extreme_upper),0.1), cat=extreme_whole, control = list(fnscale=-1))
(loglik_sigR <- gpd_whole_sigR$value)

#Constant threshold (1.07) with GPD parameters verying across regions

extreme_upper <- cat_upper$mag[cat_upper$mag > 1.07] - 1.07
(gpd_fit_upper <- optim(GPD_LL, par=c(mean(extreme_upper), 0.1), z=extreme_upper, control = list(fnscale=-1)))

extreme_lower <- cat_lower$mag[cat_lower$mag > 1.07] - 1.07
(gpd_fit_lower <- optim(GPD_LL, par=c(mean(extreme_lower), 0.1), z=extreme_lower, control = list(fnscale=-1)))

(loglik_spatial <- gpd_fit_lower$value + gpd_fit_upper$value)

(test_stat <- -2*(loglik_sigR - loglik_spatial))

(p.val <- pchisq(test_stat, df = 1, lower.tail = FALSE))

#Constant threshold (1.318) with GPD scale varying and constant shape across gas-field

extreme_whole <- cat %>% filter(cat$mag > 1.318)
extreme_whole$mag<- extreme_whole$mag - 1.318
extreme_upper <- cat_upper$mag[cat_upper$mag > 1.318] - 1.318
extreme_lower <- cat_lower$mag[cat_lower$mag > 1.318] - 1.318
gpd_whole_sigR <- optim(GPD_LL_sigR, par=c(mean(extreme_lower), mean(extreme_upper),0.1), cat=extreme_whole, control = list(fnscale=-1))
(loglik_sigR <- gpd_whole_sigR$value)

#Constant threshold (1.318) with GPD parameters verying across regions

extreme_upper <- cat_upper$mag[cat_upper$mag > 1.318] - 1.318
(gpd_fit_upper <- optim(GPD_LL, par=c(mean(extreme_upper), 0.1), z=extreme_upper, control = list(fnscale=-1)))

extreme_lower <- cat_lower$mag[cat_lower$mag > 1.318] - 1.318
(gpd_fit_lower <- optim(GPD_LL, par=c(mean(extreme_lower), 0.1), z=extreme_lower, control = list(fnscale=-1)))

(loglik_spatial <- gpd_fit_lower$value + gpd_fit_upper$value)

(test_stat <- -2*(loglik_sigR - loglik_spatial))

(p.val <- pchisq(test_stat, df = 1, lower.tail = FALSE))

#Stepped threshold (1.2, 0.876) with GPD scale varying and constant shape across gas-field
cat$thr <- 1.2
cat$thr[cat$RD_Y <= (735000 - 0.6*cat$RD_X)] <- 0.876
extreme_whole <- cat %>% filter(cat$mag > cat$thr)
extreme_whole$mag<- extreme_whole$mag - extreme_whole$thr
extreme_upper <- cat_upper$mag[cat_upper$mag > 1.2] - 1.2
extreme_lower <- cat_lower$mag[cat_lower$mag > 0.876] - 0.876
gpd_whole_sigR <- optim(GPD_LL_sigR, par=c(mean(extreme_lower), mean(extreme_upper),0.1), cat=extreme_whole, control = list(fnscale=-1))
(loglik_sigR <- gpd_whole_sigR$value)

#Stepped threshold (1.2, 0.876) with GPD parameters verying across regions

extreme_upper <- cat_upper$mag[cat_upper$mag > 1.2] - 1.2
(gpd_fit_upper <- optim(GPD_LL, par=c(mean(extreme_upper), 0.1), z=extreme_upper, control = list(fnscale=-1)))

extreme_lower <- cat_lower$mag[cat_lower$mag > 0.876] - 0.876
(gpd_fit_lower <- optim(GPD_LL, par=c(mean(extreme_lower), 0.1), z=extreme_lower, control = list(fnscale=-1)))

(loglik_spatial <- gpd_fit_lower$value + gpd_fit_upper$value)

(test_stat <- -2*(loglik_sigR - loglik_spatial))

(p.val <- pchisq(test_stat, df = 1, lower.tail = FALSE))

#--------------------------------LRTs for variation in parameters across two regions allowing distance d between regions------------------------------
d=4000
#Plotting
dev.new(width=7.5, height=6,noRStudioGD = TRUE)
par(mfrow=c(1,1),bg='transparent')
plot(cat$RD_X, cat$RD_Y, asp=1, xlab="Easting (m)", ylab = "Northing (m)")
lines(gron_outline$X, gron_outline$Y, pch=19, cex=0.5, col="red")
y <- 735000 - 0.6*seq(220000, 280000, by=1) 
y_up <- 735000 - 0.6*seq(220000, 280000, by=1) + d
y_down <- 735000 - 0.6*seq(220000, 280000, by=1) - d
lines(seq(220000, 280000, by=1), y, col="blue")
lines(seq(220000, 280000, by=1), y_up, col="blue", lty="dashed")
lines(seq(220000, 280000, by=1), y_down, col="blue", lty="dashed")

cat_upper <- cat[cat$RD_Y > (735000 - 0.6*cat$RD_X) + d,]
cat_lower <- cat[cat$RD_Y <= (735000 - 0.6*cat$RD_X) - d,]

cat_upper$index <- seq(length(cat_upper$mag), 1, by=-1)
cat_lower$index <- seq(length(cat_lower$mag), 1, by=-1)

cat <- rbind(cat_upper, cat_lower)

#LRT
##H0: Parameters constant across space
#H1: Varying scale with constant shape

#Constant threshold (1.07) with constant GPD across gas-field
extreme_whole <- cat$mag[cat$mag > 1.07] - 1.07
(gpd.fit_whole <- optim(GPD_LL, par=c(mean(extreme_whole), 0.1), z=extreme_whole, control = list(fnscale=-1)))
(loglik_constant <- gpd.fit_whole$value)

#Constant threshold (1.07) with GPD scale varying and constant shape across gas-field

extreme_whole <- cat %>% filter(cat$mag > 1.07)
extreme_whole$mag<- extreme_whole$mag - 1.07
extreme_upper <- cat_upper$mag[cat_upper$mag > 1.07] - 1.07
extreme_lower <- cat_lower$mag[cat_lower$mag > 1.07] - 1.07
gpd_whole_sigR <- optim(GPD_LL_sigR, par=c(mean(extreme_lower), mean(extreme_upper),0.1), cat=extreme_whole, control = list(fnscale=-1))
(loglik_sigR <- gpd_whole_sigR$value)

(test_stat <- -2*(loglik_constant - loglik_sigR))

(p.val <- pchisq(test_stat, df = 1, lower.tail = FALSE))

#Constant threshold (1.318) with constant GPD across gas-field
extreme_whole <- cat$mag[cat$mag > 1.318] - 1.318
(gpd.fit_whole <- optim(GPD_LL, par=c(mean(extreme_whole), 0.1), z=extreme_whole, control = list(fnscale=-1)))
(loglik_constant <- gpd.fit_whole$value)

#Constant threshold (1.318) with GPD scale varying and constant shape across gas-field

extreme_whole <- cat %>% filter(cat$mag > 1.318)
extreme_whole$mag<- extreme_whole$mag - 1.318
extreme_upper <- cat_upper$mag[cat_upper$mag > 1.318] - 1.318
extreme_lower <- cat_lower$mag[cat_lower$mag > 1.318] - 1.318
gpd_whole_sigR <- optim(GPD_LL_sigR, par=c(mean(extreme_lower), mean(extreme_upper),0.1), cat=extreme_whole, control = list(fnscale=-1))
(loglik_sigR <- gpd_whole_sigR$value)

(test_stat <- -2*(loglik_constant - loglik_sigR))

(p.val <- pchisq(test_stat, df = 1, lower.tail = FALSE))

#Stepped threshold (1.2, 0.876) with constant GPD across gas-field
extreme_upper <- cat_upper$mag[cat_upper$mag > 1.2] - 1.2
extreme_lower <- cat_lower$mag[cat_lower$mag > 0.876] - 0.876
extreme_whole <- c(extreme_upper, extreme_lower)
(gpd.fit_whole <- optim(GPD_LL, par=c(mean(extreme_whole), 0.1), z=extreme_whole, control = list(fnscale=-1)))
(loglik_constant <- gpd.fit_whole$value)


#Stepped threshold (1.2, 0.876) with GPD scale varying and constant shape across gas-field
cat$thr <- 1.2
cat$thr[cat$RD_Y <= (735000 - 0.6*cat$RD_X)] <- 0.876
extreme_whole <- cat %>% filter(cat$mag > cat$thr)
extreme_whole$mag<- extreme_whole$mag - extreme_whole$thr
extreme_upper <- cat_upper$mag[cat_upper$mag > 1.2] - 1.2
extreme_lower <- cat_lower$mag[cat_lower$mag > 0.876] - 0.876
gpd_whole_sigR <- optim(GPD_LL_sigR, par=c(mean(extreme_lower), mean(extreme_upper),0.1), cat=extreme_whole, control = list(fnscale=-1))
(loglik_sigR <- gpd_whole_sigR$value)


(test_stat <- -2*(loglik_constant - loglik_sigR))

(p.val <- pchisq(test_stat, df = 1, lower.tail = FALSE))


#LRT 
#H0: Parameters constant across space
#H1: Parameters significantly different in two regions

#Constant threshold (1.07) with constant GPD across gas-field
extreme_whole <- cat$mag[cat$mag > 1.07] - 1.07
(gpd.fit_whole <- optim(GPD_LL, par=c(mean(extreme_whole), 0.1), z=extreme_whole, control = list(fnscale=-1)))
(loglik_constant <- gpd.fit_whole$value)

#Constant threshold (1.07) with GPD parameters verying across regions

extreme_upper <- cat_upper$mag[cat_upper$mag > 1.07] - 1.07
(gpd_fit_upper <- optim(GPD_LL, par=c(mean(extreme_upper), 0.1), z=extreme_upper, control = list(fnscale=-1)))

extreme_lower <- cat_lower$mag[cat_lower$mag > 1.07] - 1.07
(gpd_fit_lower <- optim(GPD_LL, par=c(mean(extreme_lower), 0.1), z=extreme_lower, control = list(fnscale=-1)))

(loglik_spatial <- gpd_fit_lower$value + gpd_fit_upper$value)  

(test_stat <- -2*(loglik_constant - loglik_spatial))

(p.val <- pchisq(test_stat, df = 2, lower.tail = FALSE))

#Constant threshold (1.318) with constant GPD across gas-field
extreme_whole <- cat$mag[cat$mag > 1.318] - 1.318
(gpd.fit_whole <- optim(GPD_LL, par=c(mean(extreme_whole), 0.1), z=extreme_whole, control = list(fnscale=-1)))
(loglik_constant <- gpd.fit_whole$value)

#Constant threshold (1.318) with GPD parameters verying across regions

extreme_upper <- cat_upper$mag[cat_upper$mag > 1.318] - 1.318
(gpd_fit_upper <- optim(GPD_LL, par=c(mean(extreme_upper), 0.1), z=extreme_upper, control = list(fnscale=-1)))

extreme_lower <- cat_lower$mag[cat_lower$mag > 1.318] - 1.318
(gpd_fit_lower <- optim(GPD_LL, par=c(mean(extreme_lower), 0.1), z=extreme_lower, control = list(fnscale=-1)))

(loglik_spatial <- gpd_fit_lower$value + gpd_fit_upper$value)  

(test_stat <- -2*(loglik_constant - loglik_spatial))

(p.val <- pchisq(test_stat, df = 2, lower.tail = FALSE))

#Stepped threshold (1.2, 0.876) with constant GPD across gas-field
extreme_upper <- cat_upper$mag[cat_upper$mag > 1.2] - 1.2
extreme_lower <- cat_lower$mag[cat_lower$mag > 0.876] - 0.876
extreme_whole <- c(extreme_upper, extreme_lower)
(gpd.fit_whole <- optim(GPD_LL, par=c(mean(extreme_whole), 0.1), z=extreme_whole, control = list(fnscale=-1)))
(loglik_constant <- gpd.fit_whole$value)

#Stepped threshold (1.2, 0.876) with  GPD parameters verying across regions

extreme_upper <- cat_upper$mag[cat_upper$mag > 1.2] - 1.2
(gpd_fit_upper <- optim(GPD_LL, par=c(mean(extreme_upper), 0.1), z=extreme_upper, control = list(fnscale=-1)))

extreme_lower <- cat_lower$mag[cat_lower$mag > 0.876] - 0.876
(gpd_fit_lower <- optim(GPD_LL, par=c(mean(extreme_lower), 0.1), z=extreme_lower, control = list(fnscale=-1)))

(loglik_spatial <- gpd_fit_lower$value + gpd_fit_upper$value)  

(test_stat <- -2*(loglik_constant - loglik_spatial))

(p.val <- pchisq(test_stat, df = 2, lower.tail = FALSE))


#LRT
##H0: Varying scale with constant shape
#H1: Varying scale and shape

#Constant threshold (1.07) with GPD scale varying and constant shape across gas-field

extreme_whole <- cat %>% filter(cat$mag > 1.07)
extreme_whole$mag<- extreme_whole$mag - 1.07
extreme_upper <- cat_upper$mag[cat_upper$mag > 1.07] - 1.07
extreme_lower <- cat_lower$mag[cat_lower$mag > 1.07] - 1.07
gpd_whole_sigR <- optim(GPD_LL_sigR, par=c(mean(extreme_lower), mean(extreme_upper),0.1), cat=extreme_whole, control = list(fnscale=-1))
(loglik_sigR <- gpd_whole_sigR$value)

#Constant threshold (1.07) with GPD parameters verying across regions

extreme_upper <- cat_upper$mag[cat_upper$mag > 1.07] - 1.07
(gpd_fit_upper <- optim(GPD_LL, par=c(mean(extreme_upper), 0.1), z=extreme_upper, control = list(fnscale=-1)))

extreme_lower <- cat_lower$mag[cat_lower$mag > 1.07] - 1.07
(gpd_fit_lower <- optim(GPD_LL, par=c(mean(extreme_lower), 0.1), z=extreme_lower, control = list(fnscale=-1)))

(loglik_spatial <- gpd_fit_lower$value + gpd_fit_upper$value)

(test_stat <- -2*(loglik_sigR - loglik_spatial))

(p.val <- pchisq(test_stat, df = 1, lower.tail = FALSE))

#Constant threshold (1.318) with GPD scale varying and constant shape across gas-field

extreme_whole <- cat %>% filter(cat$mag > 1.318)
extreme_whole$mag<- extreme_whole$mag - 1.318
extreme_upper <- cat_upper$mag[cat_upper$mag > 1.318] - 1.318
extreme_lower <- cat_lower$mag[cat_lower$mag > 1.318] - 1.318
gpd_whole_sigR <- optim(GPD_LL_sigR, par=c(mean(extreme_lower), mean(extreme_upper),0.1), cat=extreme_whole, control = list(fnscale=-1))
(loglik_sigR <- gpd_whole_sigR$value)

#Constant threshold (1.318) with GPD parameters verying across regions

extreme_upper <- cat_upper$mag[cat_upper$mag > 1.318] - 1.318
(gpd_fit_upper <- optim(GPD_LL, par=c(mean(extreme_upper), 0.1), z=extreme_upper, control = list(fnscale=-1)))

extreme_lower <- cat_lower$mag[cat_lower$mag > 1.318] - 1.318
(gpd_fit_lower <- optim(GPD_LL, par=c(mean(extreme_lower), 0.1), z=extreme_lower, control = list(fnscale=-1)))

(loglik_spatial <- gpd_fit_lower$value + gpd_fit_upper$value)

(test_stat <- -2*(loglik_sigR - loglik_spatial))

(p.val <- pchisq(test_stat, df = 1, lower.tail = FALSE))

#Stepped threshold (1.2, 0.876) with GPD scale varying and constant shape across gas-field
cat$thr <- 1.2
cat$thr[cat$RD_Y <= (735000 - 0.6*cat$RD_X)] <- 0.876
extreme_whole <- cat %>% filter(cat$mag > cat$thr)
extreme_whole$mag<- extreme_whole$mag - extreme_whole$thr
extreme_upper <- cat_upper$mag[cat_upper$mag > 1.2] - 1.2
extreme_lower <- cat_lower$mag[cat_lower$mag > 0.876] - 0.876
gpd_whole_sigR <- optim(GPD_LL_sigR, par=c(mean(extreme_lower), mean(extreme_upper),0.1), cat=extreme_whole, control = list(fnscale=-1))
(loglik_sigR <- gpd_whole_sigR$value)

#Stepped threshold (1.2, 0.876) with GPD parameters verying across regions

extreme_upper <- cat_upper$mag[cat_upper$mag > 1.2] - 1.2
(gpd_fit_upper <- optim(GPD_LL, par=c(mean(extreme_upper), 0.1), z=extreme_upper, control = list(fnscale=-1)))

extreme_lower <- cat_lower$mag[cat_lower$mag > 0.876] - 0.876
(gpd_fit_lower <- optim(GPD_LL, par=c(mean(extreme_lower), 0.1), z=extreme_lower, control = list(fnscale=-1)))

(loglik_spatial <- gpd_fit_lower$value + gpd_fit_upper$value)

(test_stat <- -2*(loglik_sigR - loglik_spatial))

(p.val <- pchisq(test_stat, df = 1, lower.tail = FALSE))





