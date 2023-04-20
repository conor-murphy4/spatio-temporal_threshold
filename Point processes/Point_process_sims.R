library(tidyverse)
library(dplyr)

source("Point processes/Covariates/making_covariate_arrays/covariate_array_functions.R")

Opt  <- readRDS("Point processes/Covariates/Output/data/Opt_S4.RDS")
Opt
B0 <- Opt$par[1]
B1 <- Opt$par[2]
sigma <- Opt$par[3]

#Covariates
ics <- readRDS("Point processes/Covariates/00_data/derived/covariates/incremental_Coulomb_stress_mats.RDS")

#Events
cat <- read.csv("Data/Events/2022-04-12_15-09-25_cat.csv")

#Groningen outline
gron_outline <- read.csv("Data/Geophones/Geophone data/Groningen_Field_outline.csv")

#Geophones
geo <- read.csv("Data/Geophones/Geophone data/Geo_orig_post_process.csv")

plot(cat$RD_X, cat$RD_Y)
points(gron_outline$X, gron_outline$Y, pch=19, col="blue")
points(geo$Xcoord, geo$Ycoord, pch=4, col="red", cex=1.5)

#Regular grid
a <- seq(227750, 269750, by=500)
b <- seq(565250, 615750, by=500)
grid <- expand.grid(a,b)
points(grid$Var1, grid$Var2, pch="|")

# icsgrid <- expand.grid(ics$xgrid, ics$ygrid)
# points(icsgrid$Var1, icsgrid$Var2, col="red")

#Accessing elements
head(ics)
dim(ics[[1]][,,45]) #dimensions of matrix of ics across space at point in time
print(ics[[1]][40,45,]) #vector of ics over time at one point in space
length(ics[[1]][40,45,]) #length of vector of ics over time at one point in space

#Temporal derivative
ics_rate <- temporal_difference_cov_mats(ics)
dim(ics_rate[[1]])

#Intensity (w/o temporal derivative)
intensity <- ics
intensity[[1]] <- B0*ics_rate[[1]]*exp(B1*ics[[1]]) #This is our space-time varying intensity

#Now, we need to isolate grid cubes, take intensity value, calculate no. of points 
#within grid box by n=intensity[x,y,t]*volume and locate events uniformly across the cube

