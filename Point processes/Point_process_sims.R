rm(list=ls())
library(tidyverse)
library(dplyr)

source("Point processes/Covariates/making_covariate_arrays/covariate_array_functions.R")
source("Point processes/Covariates/S4_ics_exponential/S4_functions.R")
source("helper_functions.R")

Opt  <- readRDS("Point processes/Covariates/Output/data/Opt_S4.RDS")
Opt
B0 <- Opt$par[1]
B1 <- Opt$par[2]
sigma <- Opt$par[3]

#Covariates
ics <- readRDS("Point processes/Covariates/00_data/derived/covariates/incremental_Coulomb_stress_mats.RDS")
ICS_new <- readRDS("Point processes/Covariates/00_data/derived/covariates/ICS_1995_2025.RDS")
par(mfrow=c(1,1))
diff_ics <- ics$ics[,,45] - ICS_new$ics[,,45]
image.plot(diff_ics)
max(diff_ics, na.rm = TRUE)
ICS_new$tgrid
min(ICS_new$ics, na.rm = TRUE)
max(ICS_new$ics, na.rm=TRUE)
par(mfrow=c(2,2))
t_year <- 2015
image.plot(ics$xgrid, ics$ygrid,0.1*ics$ics[,,ics$tgrid==t_year], ylab="Northing (m)", xlab="Easting (m)", main="2015")
image.plot(ICS_new$xgrid, ICS_new$ygrid, ICS_new$ics[,,ICS_new$tgrid==t_year], ylab="Northing (m)", xlab="Easting (m)", main="2015")

#Events
cat <- read.csv("Data/Events/2022-04-12_15-09-25_cat.csv")

#Groningen outline
gron_outline <- read.csv("Data/Geophones/Geophone data/Groningen_Field_outline.csv")

#Geophones
geo <- read.csv("Data/Geophones/Geophone data/Geo_orig_post_process.csv")

plot(cat$RD_X, cat$RD_Y)
points(gron_outline$X, gron_outline$Y, pch=19, col="blue")
points(geo$Xcoord, geo$Ycoord, pch=4, col="red", cex=1.5)
plot(geo$Xcoord, geo$Ycoord)
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

#Temporal derivative ####NEED TO UPDATE FUNCTION!!
ics_rate <- temporal_difference_cov_mats(ics)
dim(ics_rate[[1]])

#Intensity 
intensity <- ics
intensity[[1]] <- B0*ics_rate[[1]]*exp(B1*ics[[1]]) #This is our space-time varying intensity
names(intensity)[1] <- "intensity"

REPLACE_NEGATIVES <- TRUE
REPLACEMENT_VALUE <- 1e-9

intensity<- replace_negatives_cov_mats(
  cov_mats = intensity,
  replace = REPLACE_NEGATIVES,
  replacement_value = REPLACEMENT_VALUE)

ICS_new <- ICS_1995_2025 
ICS_old <- ics_mats
#Now, we need to isolate grid cubes, take intensity value, calculate no. of points 
#within grid box by n=intensity[x,y,t]*volume and locate events uniformly across the cube
library(fields)
PLOT_PATH <- "Point processes/Simulations/Plots/PP_sim1.pdf"
pdf(PLOT_PATH, width = 7,height = 5)
plot_cov_mats_with_sim(intensity,PP_sim1_list,cov_name = "Intensity and Simulation")
dev.off()

#Simulating PP
sim_PP <- function(intensity, grid_box_area, x_length = 250, y_length = 250, t_length = 0.5 ){
  Int_intensity <- sum(intensity[[1]], na.rm = TRUE)*grid_box_area
  n_total <- rpois(1, Int_intensity)
  # x_loc <- y_loc <- t_loc <- numeric(n_total)
  # tgrid <- intensity$tgrid
  # ygrid <- intensity$ygrid
  # xgrid <- intensity$xgrid
  n_elements <- length(c(intensity[[1]]))
  probs <- c(intensity[[1]])
  probs[is.na(probs)] <- 0
  grid_boxes <- sample( c(1:n_elements), n_total, prob=probs, replace = TRUE)
  points <- box.no2xyt(grid_boxes, minX = min(intensity$xgrid), minY = min(intensity$ygrid), minYr = min(intensity$tgrid))
  x_loc <- points[[1]] + x_length*runif(n_total, -1, 1)
  y_loc <- points[[2]] + y_length*runif(n_total, -1, 1)
  t_loc <- points[[3]] + t_length*runif(n_total, -1, 1)
  # for (t in 1:length(tgrid)){
  #   for(y in 1:length(ygrid)){
  #     for(x in 1:length(xgrid)){
  #       int <- intensity[[1]][x,y,t]
  #       if(!is.na(int)){
  #         Lambda <- grid_box_area*int
  #         N[[1]][x,y,t] <- rpois(1, lambda = Lambda)
  #         x_new <- runif(N[[1]][x,y,t], xgrid[x] - x_length, xgrid[x] + x_length)
  #         y_new <- runif(N[[1]][x,y,t], ygrid[y] - y_length, ygrid[y] + y_length)
  #         t_new <- runif(N[[1]][x,y,t], tgrid[t] - t_length, tgrid[t] + t_length)
  #         x_loc <- c(x_loc, x_new)
  #         y_loc <- c(y_loc, y_new)
  #         t_loc <- c(t_loc, t_new)
  #       }
  #       
  #     }
  #   }
  # }
  mags <- rgpd(length(x_loc), scale=1, shape=0.1)
  sim_data <- data.frame(X=x_loc, Y=y_loc, t=t_loc, M=mags)
  return(sim_data)
}
grid_box_area = 0.5^2
PP_sim1 <- sim_PP(intensity, grid_box_area)

#Taking year out of time
PP_sim1$Year <- stringr::str_sub(PP_sim1$t,1,4)

#Convert uniformly located times to dates
library(lubridate)
date <- ymd(PP_sim1$t)

Int_intensity <- sum(intensity[[1]], na.rm = TRUE)*grid_box_area
n_total <- rpois(1, Int_intensity)
max(intensity[[1]], na.rm = TRUE)

PP_sim1[PP_sim1$Year == 1995,]

PP_matrix <- data.matrix(PP_sim1)

PP_sim1_list <- list()
Years <- sort(unique(PP_sim1$Year))
for(i in 1:length(Years)){
  current_year <- Years[i]
  current_PP <- PP_sim1[PP_sim1$Year==current_year,1:2]
  name <- paste0(current_year)
  PP_sim1_list[[i]] <- as.matrix(current_PP)
}
names(PP_sim1_list) <- Years

pdf_path <- "Point processes/Simulations/Plots/PP_sim1.pdf"
pdf(pdf_path, width = 7,height = 5)
for(i in 1:length(Years)){
  plot_data <- PP_sim1_list[[i]]
  plot(gron_outline$X, gron_outline$Y, asp=1, main=  paste0("Simulated events - ", Years[i]), xlab = "Easting(m)", ylab = "Northing(m)", pch=19, cex=0.5)
  points(plot_data[,1], plot_data[,2], pch=19, col="blue", cex=1.5)
}
dev.off()


