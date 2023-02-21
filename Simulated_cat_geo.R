setwd("C:/Users/murphyc4/OneDrive - Lancaster University/STOR-i/PhD/Projects/Spatio-temporal") 
source("Spatio-temporal threshold.R")

#Simple simulation of space-time varying earthquakes/geophones
num_events <- 1960
M <- rgpd(num_events, shape = 0.1, scale=20)
#New limits of rectangle!!!!!!!!!!!!!!!!!!!!!!!!!!!!!---------------------------------------------------------
X <- runif(num_events, 0, 45000)
Y <- runif(num_events, 0, 79000)
GeoX <- runif(60, 0, 45000)
GeoY <- runif(60, 0, 79000)
plot(X,Y)
points(GeoX, GeoY, col="red", pch=19, cex=2)

EQ.df <- data.frame(X=X, Y=Y, M=M)
EQ.df$Year <- rep(seq(1995, 2022, by=1), times=1, each=70)
Geo.df <- data.frame(X=GeoX, Y=GeoY)
Geo.df$Start_Year <- c(rep(1995, 35), rep(1996,4), rep(2000, 5), rep(2005, 2), rep(2015, 14)) 
Geo.df$End_Year <- c(rep(2005, 5 ), rep(2015, 3), rep(2022, 9), rep(2016, 4), rep(2022, 39) )

#Visualising individual years
year<-1995
plot(EQ.df$X[EQ.df$Year ==year],EQ.df$Y[EQ.df$Year ==year])
Current_Geo <- Geo_current(Geo.df, year)
points(Current_Geo$X, Current_Geo$Y, col="red", cex=2, pch=19)
year<-1996
points(EQ.df$X[EQ.df$Year ==year],EQ.df$Y[EQ.df$Year ==year])
Current_Geo <- Geo_current(Geo.df, year)
points(Current_Geo$X, Current_Geo$Y, col="red", cex=2, pch=19)
year<-1997
points(EQ.df$X[EQ.df$Year ==year],EQ.df$Y[EQ.df$Year ==year])
Current_Geo <- Geo_current(Geo.df, year)
points(Current_Geo$X, Current_Geo$Y, col="red", cex=2, pch=19)
year<-1998
points(EQ.df$X[EQ.df$Year ==year],EQ.df$Y[EQ.df$Year ==year])
Current_Geo <- Geo_current(Geo.df, year)
points(Current_Geo$X, Current_Geo$Y, col="red", cex=2, pch=19)
year<-1999
points(EQ.df$X[EQ.df$Year ==year],EQ.df$Y[EQ.df$Year ==year])
Current_Geo <- Geo_current(Geo.df, year)
points(Current_Geo$X, Current_Geo$Y, col="red", cex=2, pch=19)
year<-2000
points(EQ.df$X[EQ.df$Year ==year],EQ.df$Y[EQ.df$Year ==year])
Current_Geo <- Geo_current(Geo.df, year)
points(Current_Geo$X, Current_Geo$Y, col="red", cex=2, pch=19)
year<-2015
points(EQ.df$X[EQ.df$Year <= year],EQ.df$Y[EQ.df$Year <= year])
Current_Geo <- Geo_current(Geo.df, year)
points(Current_Geo$X, Current_Geo$Y, col="red", cex=2, pch=19)
year<-2022
points(EQ.df$X[EQ.df$Year <= year],EQ.df$Y[EQ.df$Year <= year])
Current_Geo <- Geo_current(Geo.df, year)
points(Current_Geo$X, Current_Geo$Y, col="red", cex=2, pch=19)

#Calculating V(x,t)
V_xt <- V_2d_time(EQ.df, Geo.df)

#Censoring
theta <- 0.2 #Chosen scaling factor to estimate
cens_thr <- theta*V_xt #Space-time varying threshold below which observations are censored
keep <- M > cens_thr 
Sim_EQ.df <- data.frame(M=M,X=X, Y=Y, V=V_xt)
Sim_EQ.df <- subset(Sim_EQ.df, keep) #Subsetting dataset according to keep

points(Sim_EQ.df$X, Sim_EQ.df$Y, pch=19, col = "green") #Visualise observed points over original

#Estimating parameters on above simulated dataset
thresholds <- seq(0,1.0, by=0.1)
Spat_temp_thr(Sim_EQ.df, thresholds = thresholds)

#-----------------------------Poisson process simulation---------------------------------

#library(spatstat)


#Region size
x_l <- 227000 
x_u <- 272000
y_l <- 535000
y_u <- 614000
x_delta <- (x_u - x_l)
y_delta <- (y_u - y_l)
area <- x_delta*y_delta/1000000

#Spatially-varying intensity
alpha <-1
beta <-1
scale <- 10
x_init <- (x_u+x_l)/2
y_init <- (y_u+y_l)/2
xy_init <- c(x_init, y_init)
lambda_fun <- function(x,y){
  return(10*exp(-(alpha*((x-x_init)/1000)^2 + beta*((y-y_init)/1000)^2)/(2*scale^2)))
}

#Max value of intensity
neg_lambda <- function(x){
  return(-lambda_fun(x[1],x[2]))
}
opt_max <- optim(xy_init,neg_lambda,method = "L-BFGS-B",lower = c(x_l,y_l), upper = c(x_u, y_u))
max_lambda <- -opt_max$value

#Thinning prob function
thin <- function(x,y){
  return(lambda_fun(x,y)/max_lambda)
}

#Simulate Poisson process
num_points <- rpois(1,max_lambda*area)
xx <- runif(num_points, 0, x_delta) + x_l
yy <- runif(num_points, 0, y_delta) + y_l

#Spatially-dependent thinning probs
p <- thin(xx,yy)

keep <- runif(num_points, 0, 1) < p

xx_kept <- xx[keep]
yy_kept <- yy[keep]

plot(xx_kept, yy_kept, ylim=c(y_l,y_u), xlim=c(x_l,x_u))

