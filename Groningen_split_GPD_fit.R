rm(list=ls())
source("helper_functions.R")
source("thresh_qq_metric.R")
#Events
cat <- read.csv("Data/Events/2022-04-12_15-09-25_cat.csv")

#Groningen outline
gron_outline <- read.csv("Data/Geophones/Geophone data/Groningen_Field_outline.csv")

plot(cat$RD_X, cat$RD_Y, asp=1)
points(gron_outline$X, gron_outline$Y, pch=19, col="blue")
y <- 735000 - 0.6*seq(220000, 280000, by=1)
points(seq(220000, 280000, by=1), y)


cat_upper <- cat[cat$RD_Y > (735000 - 0.6*cat$RD_X),]
cat_lower <- cat[cat$RD_Y <= (735000 - 0.6*cat$RD_X),]

#Threshold selection
thresh_whole <- quantile(cat$mag, seq(0,0.95, by=0.01))
thr_selection_whole <- thresh_qq_metric(cat$mag, thresh_whole, k=200)

thresh_upper <- quantile(cat_upper$mag, seq(0,0.95, by=0.01))
thr_selection_upper <- thresh_qq_metric(cat_upper$mag, thresh_upper, k=200)

thresh_lower <- quantile(cat_lower$mag, seq(0,0.95, by=0.01))
thr_selection_lower <- thresh_qq_metric(cat_lower$mag, thresh_lower, k=200)

thr_selection_lower$thresh
thr_selection_upper$thresh

#Fitting separate GPDs over separate thresholds
thr_selection_lower$par
thr_selection_upper$par

#Fitting separate GPDs over common thresholds
extreme_upper <- cat_upper$mag[cat_upper$mag > thr_selection_whole$thresh] - thr_selection_whole$thresh
(gpd_fit_upper <- optim(GPD_LL, par=c(mean(extreme_upper), 0.1), z=extreme_upper, control = list(fnscale=-1)))

extreme_lower <- cat_lower$mag[cat_lower$mag > thr_selection_whole$thresh] - thr_selection_whole$thresh
(gpd_fit_lower <- optim(GPD_LL, par=c(mean(extreme_lower), 0.1), z=extreme_lower, control = list(fnscale=-1)))

#Common threshold of 1.07 (Zak's thesis)
extreme_upper <- cat_upper$mag[cat_upper$mag > 1.07] - 1.07
(gpd_fit_upper <- optim(GPD_LL, par=c(mean(extreme_upper), 0.1), z=extreme_upper, control = list(fnscale=-1)))

extreme_lower <- cat_lower$mag[cat_lower$mag > 1.07] - 1.07
(gpd_fit_lower <- optim(GPD_LL, par=c(mean(extreme_lower), 0.1), z=extreme_lower, control = list(fnscale=-1)))


#LRT 
#H0: Parameters constant across space
#H1: Parameters significantly differently in two regions

#Constant threshold selected above
extreme_whole <- cat$mag[cat$mag > thr_selection_whole$thresh] - thr_selection_whole$thresh
(gpd.fit_whole <- optim(GPD_LL, par=c(mean(extreme_whole), 0.1), z=extreme_whole, control = list(fnscale=-1)))
(loglik_constant <- gpd.fit_whole$value)

(loglik_spatial <- gpd_fit_lower$value + gpd_fit_upper$value)  

(test_stat <- -2*(loglik_constant - loglik_spatial))

(p.val <- pchisq(test_stat, df = 2, lower.tail = FALSE))

#Threshold of 1.07 (Zak's paper)  
extreme_whole <- cat$mag[cat$mag > 1.07] - 1.07
(gpd.fit_whole <- optim(GPD_LL, par=c(mean(extreme_whole), 0.1), z=extreme_whole, control = list(fnscale=-1)))
(loglik_constant <- gpd.fit_whole$value)

(loglik_spatial <- gpd_fit_lower$value + gpd_fit_upper$value)  

(test_stat <- -2*(loglik_constant - loglik_spatial))

(p.val <- pchisq(test_stat, df = 2, lower.tail = FALSE))
