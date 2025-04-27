
gron_eq_cat <- read.csv("Data/Events/unrounded_after_1995_in_polygon_with_covariates.csv", header=T)
covariates <- read.csv("Data/covariates/covariates_1995-2024.csv", header=T)
geophones_deepest <- read.csv("Data/Geophones/Geophones_processed_03-07-2024_deepest_only.csv", header=T, row.names=1)
gron_outline <- read.csv('Data/Geophones/Groningen_Field_outline.csv', header=T)
gron_polygon <- read.table('Data/Geophones/polygon_for_groningen_earthquakes.txt', header=T)
gron_rect <- data.frame(X=c(210000,275000, 275000, 210000, 210000), Y=c(560000, 560000, 625000, 625000, 560000))

library(ggplot2)
library(ggspatial)
library(pracma)
library(cowplot)

# Zak's sigmoid threshold
sigmoid_threshold <- function(x, vl = 1.15, vr = 0.76, mu, zeta){
  return(vr + (vl - vr) * pnorm((mu - x) / zeta))
}

change_index <- which(gron_eq_cat$Date == as.Date("2015-12-25"))
pc_threshold <- c(rep(1.15, change_index), rep(0.76, nrow(gron_eq_cat)-change_index))

dev.new(height=5, width=10, noRStudioGD = TRUE)
par(mfrow=c(1,2), bg='transparent')
plot(as.Date(gron_eq_cat$Date),gron_eq_cat$Magnitude, xlab="Event time", ylab = "Magnitude", pch=19, col="grey", cex=0.7)
lines(as.Date(gron_eq_cat$Date),rep(1.45, nrow(gron_eq_cat)), col="red", lty=2, lwd=2)
#gron_eq_cat$zak_sigmoid <- sigmoid_threshold(x=c(1:nrow(gron_eq_cat)), mu = 746, zeta = 1)
lines(as.Date(gron_eq_cat$Date), pc_threshold, col="blue", lwd=2)

plot(gron_polygon$POINT_X, gron_polygon$POINT_Y, col="green", xlab="Easting (m)", ylab = "Northing (m)", type = 'l', lty=2, asp=1, lwd=2)
points(gron_eq_cat$Easting, gron_eq_cat$Northing, pch=19, col="grey", cex=0.7, asp=1)
lines(gron_outline$X, gron_outline$Y, col="black", asp=1)

#Rates of evens below pc threshold
mags_h <- gron_eq_cat$Magnitude[1:change_index]
p_h <- length(mags_h[mags_h < 0.76]) / length(mags_h)
se_h <- sqrt(p_h * (1 - p_h) / length(mags_h))
mags_l <- gron_eq_cat$Magnitude[(change_index+1):nrow(gron_eq_cat)]
p_l <- length(mags_l[mags_l < 0.76]) / length(mags_l)
se_l <- sqrt(p_l * (1 - p_l) / length(mags_l))


#Mags vs ICS
dev.new(height=5, width=10, noRStudioGD = TRUE)
par(mfrow=c(1,1), bg='transparent')
plot(gron_eq_cat$ICS, gron_eq_cat$Magnitude, xlab="ICS", ylab = "Magnitude", pch=19, col="grey", cex=0.7)


#Average stress plots for different years
dev.new(height=5, width=5, noRStudioGD = TRUE)
par(mfrow=c(1,1), bg='transparent')
chosen_years <- c("2020")
for(year in chosen_years){
  current_covariates <- covariates[covariates$Year == year,]
  #current_exceedances <- gron_eq_cat_exceed_V1[gron_eq_cat_exceed_V1$Year == year,]
  #Make an empty dataframe to store the average threshold for each location
  plot_df <- data.frame(Easting = numeric(0), Northing = numeric(0), avg_stress = numeric(0))
  for(easting in unique(current_covariates$Easting)){
    for(northing in unique(current_covariates$Northing)){
      plot_df <- rbind(plot_df, data.frame(Easting=easting, Northing = northing, Average_ICS = mean(current_covariates$ICS[current_covariates$Easting == easting & current_covariates$Northing == northing])))
    }
  }
  plot_df <- plot_df[complete.cases(plot_df),]
  # Edit below code to remove key from plot
  # And plot the key separately
  ggplot(plot_df, aes(x = Easting, y = Northing, fill = Average_ICS)) + geom_tile() + scale_fill_gradient(low = "red", high = "yellow") + fixed_plot_aspect(ratio = 1) + theme_classic() +
    theme(plot.background = element_blank()) + geom_point(data=gron_outline, aes(x = X, y = Y), size=0.5, shape=1, fill = "black") +
    labs(x = "Easting (m)", y = "Northing (m)", fill = "Average ICS") + coord_fixed()
    #theme(legend.position="none") +  coord_fixed()  
  
  
  # print(ggplot(plot_df, aes(x = Easting, y = Northing, fill = Average_ICS)) + geom_tile() + scale_fill_gradient(low = "red", high = "yellow")) + fixed_plot_aspect(ratio = 1) + theme_classic() +
    # theme(plot.background = element_blank()) + geom_point(data=gron_outline, aes(x = X, y = Y), size=0.5, shape=1, fill = "black") 
          #geom_point(data = current_exceedances, aes(x = Easting, y = Northing), size=2, shape=19, fill = "black"))
}

# Boxplots of magnitudes for high/low stress
mags_high_stress <- gron_eq_cat$Magnitude[gron_eq_cat$ICS > median(gron_eq_cat$ICS)]
mags_low_stress <- gron_eq_cat$Magnitude[gron_eq_cat$ICS <= median(gron_eq_cat$ICS)]
boxplot(mags_low_stress, mags_high_stress, names=c("Low stress", "High stress"), ylab = "Magnitude")

# Boxplots of magnitudes for high/low stress fr exceedances of conservative mc
gron_eq_cat_ex_conserv <- gron_eq_cat[gron_eq_cat$Magnitude > 1.45,]
mags_high_stress_conserv <- gron_eq_cat_ex_conserv$Magnitude[gron_eq_cat_ex_conserv$ICS > median(gron_eq_cat_ex_conserv$ICS)]
mags_low_stress_conserv <- gron_eq_cat_ex_conserv$Magnitude[gron_eq_cat_ex_conserv$ICS <= median(gron_eq_cat_ex_conserv$ICS)]
boxplot(mags_low_stress_conserv, mags_high_stress_conserv, names=c("Low stress", "High stress"), ylab = "Magnitudes")

# Boxplots of magnitudes for high/low stress fr exceedances of estimated threshold
threshold_est <- geo_thresh_fit[[1]]$thresh_par[1] + geo_thresh_fit[[1]]$thresh_par[2] * gron_eq_cat$V_1
gron_eq_cat_ex_est <- gron_eq_cat[gron_eq_cat$Magnitude > threshold_est,]
mags_high_stress_est <- gron_eq_cat_ex_est$Magnitude[gron_eq_cat_ex_est$ICS > median(gron_eq_cat_ex_est$ICS)]
mags_low_stress_est <- gron_eq_cat_ex_est$Magnitude[gron_eq_cat_ex_est$ICS <= median(gron_eq_cat_ex_est$ICS)]
boxplot(mags_low_stress_est, mags_high_stress_est, names=c("Low stress", "High stress"), ylab = "Magnitudes")

# Development of geophones
dates <- seq(min(as.Date(gron_eq_cat$Date)), max(as.Date(gron_eq_cat$Date)), by="day")
num_geophones_in_rect <- numeric(length(dates))
num_geophones_in_polygon <- numeric(length(dates))
num_geophones_in_gasfield <- numeric(length(dates))
geophones_in_rect <- geophones_deepest[inpolygon(geophones_deepest$Xcoord, geophones_deepest$Ycoord, gron_rect$X, gron_rect$Y),]
geophones_in_polygon <- geophones_deepest[inpolygon(geophones_deepest$Xcoord, geophones_deepest$Ycoord, gron_polygon$POINT_X, gron_polygon$POINT_Y),]
geophones_in_gasfield <- geophones_deepest[inpolygon(geophones_deepest$Xcoord, geophones_deepest$Ycoord, gron_outline$X, gron_outline$Y),]
for(i in 1:length(dates)){
  date <- dates[i]
  current_geo_in_rect <- geophones_in_rect[geophones_in_rect$Start_date <= date & geophones_in_rect$End_date >= date,]
  current_geo_in_polygon <- geophones_in_polygon[geophones_in_polygon$Start_date <= date & geophones_in_polygon$End_date >= date,]
  current_geo_in_gasfield <- geophones_in_gasfield[geophones_in_gasfield$Start_date <= date & geophones_in_gasfield$End_date >= date,]
  num_geophones_in_rect[i] <- nrow(current_geo_in_rect)
  num_geophones_in_polygon[i] <- nrow(current_geo_in_polygon)
  num_geophones_in_gasfield[i] <- nrow(current_geo_in_gasfield)
}

dev.new(height=5, width=15, noRStudioGD = TRUE)
par(mfrow=c(1,3), bg='transparent')
plot(dates, num_geophones_in_rect, xlab="Time", ylab = "Number of geophones", type='l', col="blue")
lines(dates, num_geophones_in_polygon, col="green")
lines(dates, num_geophones_in_gasfield, col="black")

chosen_years <- c("2010", "2020")
for(year in chosen_years){
  current_geophones <- geophones_deepest[stringr::str_sub(geophones_deepest$Start_date, 1, 4) <= year & stringr::str_sub(geophones_deepest$End_date, 1, 4) >= year,]
  plot(gron_polygon$POINT_X, gron_polygon$POINT_Y, col="green", xlab="Easting (m)", ylab = "Northing (m)", type = 'l', lty=2, lwd=2, asp=1)
  points(current_geophones$Xcoord, current_geophones$Ycoord, pch=4, col="blue", cex=0.7)
  lines(gron_outline$X, gron_outline$Y, col="black")
}

gron_eq_cat[gron_eq_cat$Northing == max(gron_eq_cat$Northing),]
max(covariates$Northing)


# Development of geophones
dates <- seq(as.Date("1995-01-01"), min(as.Date(gron_eq_cat$Date)), by="day")
num_geophones <- numeric(length(dates))
for(i in 1:length(dates)){
  date <- dates[i]
  current_geo <- geophones_deepest[geophones_deepest$Start_date <= date & geophones_deepest$End_date >= date,]
  num_geophones[i] <- nrow(current_geo)
}

dates[num_geophones < 5]

last_date <- dates[num_geophones < 5][1]
nrow(gron_eq_cat[gron_eq_cat$Date <= last_date,])

gron_eq_cat[gron_eq_cat$Date < as.Date("1995-05-01"),]
date <- as.Date("1995-04-26")
current_geo <- geophones_deepest[geophones_deepest$Start_date <= date & geophones_deepest$End_date >= date,]
plot(current_geo$Xcoord, current_geo$Ycoord, pch=4, col="blue", cex=0.7, xlab="Easting (m)", ylab = "Northing (m)", asp=1)
lines(gron_polygon$POINT_X, gron_polygon$POINT_Y, col="green", lwd=2)
lines(gron_outline$X, gron_outline$Y, col="black")
lines(gron_rect$X, gron_rect$Y, col="red", lwd=2)
nrow(current_geo)
