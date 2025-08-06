source("src/distance_to_nearest_geo.R")
library(pracma)
library(ggplot2)
source("src/eqd_geo.R")

gron_outline <- read.csv('Data/Geophones/Groningen_Field_outline.csv', header=T)
gron_polygon <- read.table('Data/Geophones/polygon_for_groningen_earthquakes.txt', header=T)
geophones <- read.csv("Data/Geophones/Geophones_processed_03-07-2024_without_duplicates.csv", header=T, row.names = 1)
geophones_deepest <- read.csv("Data/Geophones/Geophones_processed_03-07-2024_deepest_only.csv", header=T, row.names = 1)
gron_eq_cat_all_netherlands <- read.csv("Data/Events/unrounded_after_geophone_start.csv", header=T, row.names = 1)
gron_eq_cat <- read.csv("Data/Events/unrounded_after_1995_in_polygon_with_covariates.csv", header=T)
gron_eq_cat_all <- read.csv("Data/Events/unrounded_after_geophone_start_in_polygon_with_V_3d.csv", header=T)
gron_eq_cat_old <- read.csv("Data/Events/2022-04-12_15-09-25_cat.csv", header=T)
covariates <- read.csv("Data/covariates/covariates_1995-2024.csv", header=T)
covariates_1957_2055 <- read.csv("Data/covariates/covariates_1957-2055.csv", header=T)
covariates_2055 <- read.csv("Data/covariates/covariates_1995-2055.csv", header=T)
covariates_1957_2055_in_poly <- read.csv("Data/covariates/covariates_1957-2055_in_poly.csv", header=T)

covariates_1957_2055_in_poly <- covariates_1957_2055[inpolygon(covariates_1957_2055$x, covariates_1957_2055$y, gron_polygon$POINT_X, gron_polygon$POINT_Y),]

covariates <- covariates_2055[covariates_2055$Date <= max(gron_eq_cat$Date) & covariates_2055$Date >= "1995-04-01",]

# Creating ICS_max 
covariates_2055$ICS_max <- rep(NA, nrow(covariates_2055))
start_time <- Sys.time()
for(date in unique(covariates_2055$Date)[2]){
  print(date)
  current_covariates <- covariates_2055[covariates_2055$Date == date,]
  covariates_up_to_date <- covariates_2055[covariates_2055$Date <= date,]
  # For each locatiion, find the maximum ICS value up to that date
  # and assign it to the ICS_max column for that date
  for(i in 1:nrow(current_covariates)){
    print(paste(i, "of", nrow(current_covariates), "for date", date))
    current_location <- current_covariates[i,]
    ics_up_to_date_in_location <- covariates_up_to_date$ICS[covariates_up_to_date$Easting == current_location$Easting & covariates_up_to_date$Northing == current_location$Northing]
    max_ics <- max(ics_up_to_date_in_location, na.rm = TRUE)
    covariates_2055$ICS_max[covariates_2055$Date == date & covariates_2055$Easting == current_location$Easting & covariates_2055$Northing == current_location$Northing] <- max_ics
  }
}

#How could I change the above loop to make it run faster?
library(dplyr)


covariates_2055 <- covariates_2055 %>%
  arrange(Easting, Northing, Date) %>%
  group_by(Easting, Northing) %>%
  mutate(ICS_max = cummax(ICS)) %>%
  ungroup()


# Check if ICS_max is correctly calculated

plot(as.Date(covariates_2055$Date), covariates_2055$ICS_max, xlab = "Date", ylab = "ICS max", col="grey", type='l')

write.csv(covariates_2055, "Data/covariates/covariates_1995-2055.csv", row.names = F)

length(covariates_2055$ICS_max[is.na(covariates_2055$ICS_max)])

#Temporal derivative of ICS max
# ICS max values - previous ICS max values

#Calculating max values for first date in observation window
# Filter to relevant time period first
filtered_stress <- covariates_1957_2055_in_poly %>%
  filter(time <= as.Date("1995-03-01"))

# Compute max CoulombStress per location
max_stress_by_location <- filtered_stress %>%
  group_by(x, y) %>%
  summarise(ICS_max = max(CoulombStress, na.rm = TRUE), .groups = "drop") %>%
  rename(Easting = x, Northing = y)

# Get covariates for the earliest date
current_covariates <- covariates_2055 %>%
  filter(Date == min(Date))

# Create previous_covariates by joining with the precomputed max values
previous_covariates <- current_covariates %>%
  select(Easting, Northing) %>%
  distinct() %>%
  left_join(max_stress_by_location, by = c("Easting", "Northing")) %>%
  mutate(Date = as.Date("1995-03-01"))


# Calculating all values for dsmaxdt
covariates_2055 <- covariates_2055 %>%
  arrange(Easting, Northing, Date) %>%
  group_by(Easting, Northing) %>%
  mutate(dsmaxdt = ICS_max - lag(ICS_max)) %>%
  ungroup()

#Adding in values for first date for dsmaxdt
covariates_2055$dsmaxdt[covariates_2055$Date == min(covariates_2055$Date)] <- current_covariates$ICS_max - previous_covariates$ICS_max


#order covariates_2055 by date
covariates_2055 <- covariates_2055[order(covariates_2055$Date),]
#Remove columns from covariates_2055
covariates_2055 <- covariates_2055[,c("Date", "MonthYear","Easting", "Northing", "Depth", "V1", "ICS", "dsdt", "ICS_max", "dsmaxdt")]

# Taking out values within observation window
covariates <- covariates_2055[covariates_2055$Date <= max(covariates$Date) & covariates_2055$Date >= min(covariates$Date),] 

#Checking time development of ICS max and dsmaxdt for a location
dev.new(height=10, width=10, noRStudioGD = TRUE)
par(mfrow=c(1,1), bg = "transparent")
loc <- c(255000, 570000)
covariates_loc <- covariates[covariates$Easting == loc[1] & covariates$Northing == loc[2],]
plot(as.Date(covariates_loc$Date), covariates_loc$dsmaxdt, xlab = "Date", ylab = "ICS max temporal derivative", col="grey", type='l')
plot(as.Date(covariates_loc$Date), covariates_loc$ICS_max, xlab = "Date", ylab = "ICS max", col="grey", type='l')



#Checking dates where less than 4 geophones are present
for(date in unique(covariates$Date)){
  geos <- geo_current(geophones_deepest, date)
  if(nrow(geos) < 4){
    print(date)
  }
}

#Extracting ICS and dsdt for each event
gron_eq_cat$ICS_max <- rep(NA, nrow(gron_eq_cat))
gron_eq_cat$dsmaxdt <- rep(NA, nrow(gron_eq_cat))

for(monthyear in unique(gron_eq_cat$MonthYear)){
  print(monthyear)
  eq_cat_ind <- which(gron_eq_cat$MonthYear == monthyear)
  current_eq_cat <- gron_eq_cat[eq_cat_ind,]
  current_covariates <- covariates[covariates$MonthYear == monthyear,]
  # Find nearest location in ICS data for each earthquake
  extracted_ics_max <- rep(NA, nrow(current_eq_cat))
  extracted_dsmaxdt <- rep(NA, nrow(current_eq_cat))
  for(i in 1:nrow(current_eq_cat)){
    extracted_ics_max[i] <- current_covariates$ICS_max[which.min(sqrt((current_eq_cat$Easting[i] - current_covariates$Easting)^2 + (current_eq_cat$Northing[i] - current_covariates$Northing)^2))]
    extracted_dsmaxdt[i] <- current_covariates$dsmaxdt[which.min(sqrt((current_eq_cat$Easting[i] - current_covariates$Easting)^2 + (current_eq_cat$Northing[i] - current_covariates$Northing)^2))]
  }
  gron_eq_cat$ICS_max[eq_cat_ind] <- extracted_ics_max
  gron_eq_cat$dsmaxdt[eq_cat_ind] <- extracted_dsmaxdt
}

plot(as.Date(gron_eq_cat$Date), gron_eq_cat$dsmaxdt, xlab = "Date", ylab = "ICS temporal derivative", col="grey", type='l')
abline(h=0, col="red", lty=2)
names(gron_eq_cat)

gron_eq_cat <- gron_eq_cat[,c("Datetime", "Date", "MonthYear", "Year", "Easting", "Northing", "Depth", "V_1", "V_2", "V_3", "V_4", "Magnitude", "ICS", "dsdt", "ICS_max", "dsmaxdt")]
write.csv(gron_eq_cat, "Data/Events/unrounded_after_1995_in_polygon_with_covariates.csv", row.names = F)

# Covariate dataset -------------------------------------------------------

# covariates$V1 <- distance_to_nearest(covariates, geophones_deepest, 1)
# covariates$V2 <- distance_to_nearest(covariates, geophones_deepest, 2)
# covariates$V3 <- distance_to_nearest(covariates, geophones_deepest, 3)
# covariates$V4 <- distance_to_nearest(covariates, geophones_deepest, 4)
# 
# 
# covariates <- covariates[covariates$Date >= "1995-04-06",]
# 
# write.csv(covariates, "Data/covariates/covariates_1995-2024.csv", row.names = F)

# ------ Space-time fields of threshold, GPD parameters, summaries
# Using V1 threshold fit as example, will need to be repeated for best choice
thresh_fit_V2_ics <- readRDS("threshold_results/geo_thresh_fit_V2with_ics.rds")
covariates$threshold <- thresh$thresh_par[1] + thresh_fit_geo_ics$thresh_par[2] * covariates$V1
covariates$sigma <- thresh_fit_geo_ics$par[1] + thresh_fit_geo_ics$par[2] * covariates$ICS + thresh_fit_geo_ics$par[3] * covariates$threshold 
covariates$xi <- rep(thresh_fit_geo_ics$par[3], nrow(covariates))

# Theoretical mean excess
covariates$mean_excess <- covariates$threshold + covariates$sigma / (1 - covariates$xi)

# Endpoint
covariates$endpoint <- -(thresh_fit_geo_ics$par[1] + thresh_fit_geo_ics$par[2] * covariates$ICS) / covariates$xi

#write.csv(covariates, "Data/covariates/covariates_1995-2055.csv", row.names = F)

# covariates_obs <- covariates[covariates$Date < "2025-01-01",]
#write.csv(covariates_obs, "Data/covariates/covariates_1995-2024.csv", row.names = F)

# Average threshold across space for each date
avg_threshold_over_space <- numeric(length(unique(covariates$Date)))
for(i in 1:length(unique(covariates$Date))){
  avg_threshold_over_space[i] <- mean(covariates$threshold[covariates$Date == unique(covariates$Date)[i]])
}

# Max endpoint across space for each date
max_endpoint_over_space <- numeric(length(unique(covariates$Date)))
avg_endpoint_over_space <- numeric(length(unique(covariates$Date)))
for(i in 1:length(unique(covariates$Date))){
  max_endpoint_over_space[i] <- max(covariates$endpoint[covariates$Date == unique(covariates$Date)[i]])
  avg_endpoint_over_space[i] <- mean(covariates$endpoint[covariates$Date == unique(covariates$Date)[i]])
}

# Zak's sigmoid threshold
sigmoid_threshold <- function(x, vr = 1.15, vl = 0.76, mu, zeta){
  return(a + (b - a) * pnorm((mu - x) / zeta))
}

# Threshold on observed events
thresh_fit_V1_ics <- readRDS("threshold_results/geo_thresh_fit_V1with_ics.rds")

thresh_fit_V2_ics <- readRDS("threshold_results/geo_thresh_fit_V2with_ics.rds")
thresh_fit_V3_ics <- readRDS("threshold_results/geo_thresh_fit_V3with_ics.rds")
thresh_fit_V4_ics <- readRDS("threshold_results/geo_thresh_fit_V4with_ics.rds")

thresh_fit_logV1_ics <- readRDS("threshold_results/log_geo_thresh_fit_V1with_ics.rds")
thresh_fit_logV2_ics <- readRDS("threshold_results/log_geo_thresh_fit_V2with_ics.rds")
thresh_fit_logV3_ics <- readRDS("threshold_results/log_geo_thresh_fit_V3with_ics.rds")
thresh_fit_logV4_ics <- readRDS("threshold_results/log_geo_thresh_fit_V4with_ics.rds")

thresh_fit_sqrtV1_ics <- readRDS("threshold_results/sqrt_geo_thresh_fit_V1with_ics.rds")
thresh_fit_sqrtV2_ics <- readRDS("threshold_results/sqrt_geo_thresh_fit_V2with_ics.rds")
thresh_fit_sqrtV3_ics <- readRDS("threshold_results/sqrt_geo_thresh_fit_V3with_ics.rds")
thresh_fit_sqrtV4_ics <- readRDS("threshold_results/sqrt_geo_thresh_fit_V4with_ics.rds")

##################################-------------------------------------------------------------------
# Threshold on observed events with num_boot=1000
thresh_fit_V1_ics <- readRDS("threshold_results/k1000/geo_thresh_fit_V1with_ics.rds")
thresh_fit_V2_ics <- readRDS("threshold_results/k1000/geo_thresh_fit_V2with_ics.rds")
thresh_fit_V3_ics <- readRDS("threshold_results/k1000/geo_thresh_fit_V3with_ics.rds")
thresh_fit_V4_ics <- readRDS("threshold_results/k1000/geo_thresh_fit_V4with_ics.rds")

min(thresh_fit_V1_ics$dists, na.rm = T)
min(thresh_fit_V2_ics$dists, na.rm = T)
min(thresh_fit_V3_ics$dists, na.rm = T)
min(thresh_fit_V4_ics$dists, na.rm = T)

thresh_fit_logV1_ics <- readRDS("threshold_results/k1000/log_geo_thresh_fit_V1with_ics.rds")
thresh_fit_logV2_ics <- readRDS("threshold_results/k1000/log_geo_thresh_fit_V2with_ics.rds")
thresh_fit_logV3_ics <- readRDS("threshold_results/k1000/log_geo_thresh_fit_V3with_ics.rds")
thresh_fit_logV4_ics <- readRDS("threshold_results/k1000/log_geo_thresh_fit_V4with_ics.rds")

min(thresh_fit_logV1_ics$dists, na.rm = T)
min(thresh_fit_logV2_ics$dists, na.rm = T)
min(thresh_fit_logV3_ics$dists, na.rm = T)
min(thresh_fit_logV4_ics$dists, na.rm = T)

thresh_fit_sqrtV1_ics <- readRDS("threshold_results/k1000/sqrt_geo_thresh_fit_V1with_ics.rds")
thresh_fit_sqrtV2_ics <- readRDS("threshold_results/k1000/sqrt_geo_thresh_fit_V2with_ics.rds")
thresh_fit_sqrtV3_ics <- readRDS("threshold_results/k1000/sqrt_geo_thresh_fit_V3with_ics.rds")
thresh_fit_sqrtV4_ics <- readRDS("threshold_results/k1000/sqrt_geo_thresh_fit_V4with_ics.rds")

min(thresh_fit_sqrtV1_ics$dists, na.rm = T)
min(thresh_fit_sqrtV2_ics$dists, na.rm = T)
min(thresh_fit_sqrtV3_ics$dists, na.rm = T)
min(thresh_fit_sqrtV4_ics$dists, na.rm = T)

##################################-------------------------------------------------------------------

threshold_obs_V1 <- thresh_fit_V1_ics$thresh_par[1] + thresh_fit_V1_ics$thresh_par[2] * gron_eq_cat$V_1

threshold_obs_V2 <- thresh_fit_V2_ics$thresh_par[1] + thresh_fit_V2_ics$thresh_par[2] * gron_eq_cat$V_2
threshold_obs_V3 <- thresh_fit_V3_ics$thresh_par[1] + thresh_fit_V3_ics$thresh_par[2] * gron_eq_cat$V_3
threshold_obs_V4 <- thresh_fit_V4_ics$thresh_par[1] + thresh_fit_V4_ics$thresh_par[2] * gron_eq_cat$V_4

gron_eq_cat_exceed_V1 <- gron_eq_cat[gron_eq_cat$Magnitude > threshold_obs_V1,] 

gron_eq_cat_exceed_V2 <- gron_eq_cat[gron_eq_cat$Magnitude > threshold_obs_V2,]
gron_eq_cat_exceed_V3 <- gron_eq_cat[gron_eq_cat$Magnitude > threshold_obs_V3,]
gron_eq_cat_exceed_V4 <- gron_eq_cat[gron_eq_cat$Magnitude > threshold_obs_V4,]

# Plotting ---------------------------------------------------------------
ylimits <- c(min(gron_eq_cat$Magnitude), max(max_endpoint_over_space))
plot(as.Date(gron_eq_cat$Date), gron_eq_cat$Magnitude, xlab = "Date", ylab = "Magnitude", col="grey", ylim = ylimits)
lines(as.Date(gron_eq_cat$Date), threshold_obs_V1, col = "red", type='l', lwd=2)
lines(as.Date(unique(covariates$Date)), avg_threshold_over_space, col = "blue", type='l', lwd=2)
lines(as.Date(unique(covariates$Date)), max_endpoint_over_space, col = "black", type='l', lwd=2)
lines(as.Date(unique(covariates$Date)), avg_endpoint_over_space, col = "green", type='l', lwd=2)
lines(as.Date(gron_eq_cat$Date[1:change_index]), rep(1.15, change_index), col = "yellow", type='l', lwd=2)
lines(as.Date(gron_eq_cat$Date[(change_index+60):length(gron_eq_cat$Date)]), rep(0.76, length(gron_eq_cat$Date) - (change_index+59)), col = "yellow", type='l', lwd=2)

plot(as.Date(gron_eq_cat$Date), gron_eq_cat$Magnitude, xlab = "Date", ylab = "Magnitude", col="grey")
lines(as.Date(gron_eq_cat$Date), threshold_obs_V1, col = rgb(1, 0, 0), type='l', lwd=2)
lines(as.Date(gron_eq_cat$Date), threshold_obs_V2, col = rgb(0, 0, 1, alpha = 0.5), type='l', lwd=2)
lines(as.Date(gron_eq_cat$Date), threshold_obs_V3, col = rgb(0, 1, 0, alpha = 0.3), type='l', lwd=2)
lines(as.Date(gron_eq_cat$Date), threshold_obs_V4, col = rgb(0.5, 0.5, 0.5, alpha = 0.5), type='l', lwd=2)
lines(as.Date(gron_eq_cat$Date[1:change_index]), rep(1.15, change_index), col = "yellow", type='l', lwd=2)
lines(as.Date(gron_eq_cat$Date[(change_index+60):length(gron_eq_cat$Date)]), rep(0.76, length(gron_eq_cat$Date) - (change_index+59)), col = "yellow", type='l', lwd=2)


# Spatial plots of mean excess and endpoint for chosen timepoints
chosen_dates <- c("2000-01-01", "2005-01-01", "2010-01-01", "2015-01-01", "2020-01-01")
for(date in chosen_dates){
  current_covariates <- covariates[covariates$Date == date,]
  plot_df <- data.frame(Easting = current_covariates$Easting, Northing = current_covariates$Northing, mean_excess = current_covariates$mean_excess, endpoint = current_covariates$endpoint, threshold = current_covariates$threshold)
  print(ggplot(plot_df, aes(x = Easting, y = Northing, fill = mean_excess)) + geom_tile() + scale_fill_gradient(low = "black", high = "red") + ggtitle(paste("Mean excess on", date)))
  print(ggplot(plot_df, aes(x = Easting, y = Northing, fill = endpoint)) + geom_tile() + scale_fill_gradient(low = "blue", high = "yellow") + ggtitle(paste("Endpoint on", date)))
  print(ggplot(plot_df, aes(x = Easting, y = Northing, fill = threshold)) + geom_tile() + scale_fill_gradient(low = "blue", high = "red") + ggtitle(paste("Threshold on", date)))
}

#Average threshold plots for different years
chosen_years <- c("2000", "2005", "2010", "2015", "2020")
for(year in chosen_years){
  current_covariates <- covariates[covariates$Year == year,]
  current_exceedances <- gron_eq_cat_exceed_V1[gron_eq_cat_exceed_V1$Year == year,]
  #Make an empty dataframe to store the average threshold for each location
  plot_df <- data.frame(Easting = numeric(0), Northing = numeric(0), avg_threshold = numeric(0))
  for(easting in unique(current_covariates$Easting)){
    for(northing in unique(current_covariates$Northing)){
      plot_df <- rbind(plot_df, data.frame(Easting=easting, Northing = northing, avg_threshold = mean(current_covariates$threshold[current_covariates$Easting == easting & current_covariates$Northing == northing])))
    }
  }
  plot_df <- plot_df[complete.cases(plot_df),]
  ggplot(plot_df, aes(x = Easting, y = Northing, fill = avg_threshold)) + geom_tile() + scale_fill_gradient(low = "blue", high = "red") + + fixed_plot_aspect(ratio = 1) + theme_classic() +
          geom_point(data = current_exceedances, aes(x = Easting, y = Northing), size=2, shape=19, fill = "black")
}

#Average stress plots for different years
chosen_years <- c("2000", "2005", "2010", "2015", "2020")
for(year in chosen_years){
  current_covariates <- covariates[covariates$Year == year,]
  current_exceedances <- gron_eq_cat_exceed_V1[gron_eq_cat_exceed_V1$Year == year,]
  #Make an empty dataframe to store the average threshold for each location
  plot_df <- data.frame(Easting = numeric(0), Northing = numeric(0), avg_stress = numeric(0))
  for(easting in unique(current_covariates$Easting)){
    for(northing in unique(current_covariates$Northing)){
      plot_df <- rbind(plot_df, data.frame(Easting=easting, Northing = northing, avg_stress = mean(current_covariates$ICS[current_covariates$Easting == easting & current_covariates$Northing == northing])))
    }
  }
  plot_df <- plot_df[complete.cases(plot_df),]
  print(ggplot(plot_df, aes(x = Easting, y = Northing, fill = avg_stress)) + geom_tile() + scale_fill_gradient(low = "black", high = "red") + ggtitle(paste("Average stress and exceedance locations in", year)) +
          geom_point(data = current_exceedances, aes(x = Easting, y = Northing), size=2, shape=19, fill = "black"))
}


#Average intensity plots for different years
chosen_years <- c("2000", "2005", "2010", "2015", "2020")
for(year in chosen_years){
  current_covariates <- covariates[covariates$Year == year,]
  current_exceedances <- gron_eq_cat_exceed_V1[gron_eq_cat_exceed_V1$Year == year,]
  #Make an empty dataframe to store the average threshold for each location
  plot_df <- data.frame(Easting = numeric(0), Northing = numeric(0), avg_intensity= numeric(0))
  for(easting in unique(current_covariates$Easting)){
    for(northing in unique(current_covariates$Northing)){
      plot_df <- rbind(plot_df, data.frame(Easting=easting, Northing = northing, avg_intensity = mean(current_covariates$intensity_above_threshold_V1_per_km2[current_covariates$Easting == easting & current_covariates$Northing == northing])))
    }
  }
  plot_df <- plot_df[complete.cases(plot_df),]
  print(ggplot(plot_df, aes(x = Easting, y = Northing, fill = avg_intensity)) + geom_tile() + scale_fill_gradient(low = "blue", high = "red") + ggtitle(paste("Average intensity above V1 threshold per km2 and exceedance locations in", year)) +
          geom_point(data = current_exceedances, aes(x = Easting, y = Northing), size=2, shape=19, fill = "black"))
}

#Average prob of exceeding plots for different years
chosen_years <- c("2000", "2005", "2010", "2015", "2020")
for(year in chosen_years){
  current_covariates <- covariates[covariates$Year == year,]
  current_exceedances <- gron_eq_cat_exceed_V1[gron_eq_cat_exceed_V1$Year == year,]
  #Make an empty dataframe to store the average threshold for each location
  plot_df <- data.frame(Easting = numeric(0), Northing = numeric(0), prob_ex= numeric(0))
  for(easting in unique(current_covariates$Easting)){
    for(northing in unique(current_covariates$Northing)){
      plot_df <- rbind(plot_df, data.frame(Easting=easting, Northing = northing, prob_ex = mean(current_covariates$prob_exceedance_V1[current_covariates$Easting == easting & current_covariates$Northing == northing])))
    }
  }
  plot_df <- plot_df[complete.cases(plot_df),]
  print(ggplot(plot_df, aes(x = Easting, y = Northing, fill = prob_ex)) + geom_tile() + scale_fill_gradient(low = "blue", high = "green") + ggtitle(paste("Average probability of exceedance of V1 threshold and exceedance locations in", year)) +
          geom_point(data = current_exceedances, aes(x = Easting, y = Northing), size=2, shape=19, fill = "black"))
}

# Split according to two modes of EQ activity -----------------------------
max_dist <- max(covariates$V1)
min_ics <- min(covariates$ICS)

plot(gron_eq_cat$Easting, gron_eq_cat$Northing, asp=1, xlab="Easting (m)", ylab = "Northing (m)", col="grey")
points(gron_outline$X, gron_outline$Y, pch=19, col="blue", cex=0.3)
y_div <- 735000 - 0.6*seq(220000, 280000, by=1)
points(seq(220000, 280000, by=1), y_div, col = "red", cex=0.5)

# Dummy variable to split catalogue (1= Upper, 0 = Lower)
gron_eq_cat$region <- rep(1, nrow(gron_eq_cat))
gron_eq_cat$region[gron_eq_cat$Northing < 735000 - 0.6*gron_eq_cat$Easting] <- 0

# Threshold
thresh_fit_geo_ics <- readRDS("threshold_results/thresh_fit_V1_ics.rds")
gron_eq_cat$threshold <- thresh_fit_geo_ics$thresh_par[1] + thresh_fit_geo_ics$thresh_par[2] * gron_eq_cat$V_1

# Upper region
gron_eq_cat_upper <- gron_eq_cat[gron_eq_cat$region == 1,]
excess_upper <- gron_eq_cat_upper$Magnitude[gron_eq_cat_upper$Magnitude > gron_eq_cat_upper$threshold] - gron_eq_cat_upper$threshold[gron_eq_cat_upper$Magnitude > gron_eq_cat_upper$threshold]
pars_init <- c(mean(excess_upper), 0, 0.1)
gpd_ics_V1_fit_upper <- optim(GPD_LL_given_V_ICS, excess = excess_upper, par = pars_init, control = list(fnscale = -1), thresh_par=thresh_fit_geo_ics$thresh_par, V = gron_eq_cat_upper$V_1[gron_eq_cat_upper$Magnitude > gron_eq_cat_upper$threshold], 
                              ics = gron_eq_cat_upper$ICS[gron_eq_cat_upper$Magnitude > gron_eq_cat_upper$threshold], max_dist = max_dist, min_ics = min_ics)
gpd_ics_V1_fit_upper

# Lower region
gron_eq_cat_lower <- gron_eq_cat[gron_eq_cat$region == 0,]
excess_lower <- gron_eq_cat_lower$Magnitude[gron_eq_cat_lower$Magnitude > gron_eq_cat_lower$threshold] - gron_eq_cat_lower$threshold[gron_eq_cat_lower$Magnitude > gron_eq_cat_lower$threshold]
pars_init <- c(mean(excess_lower), 0, 0.1)
gpd_ics_V1_fit_lower <- optim(GPD_LL_given_V_ICS, excess = excess_lower, par = pars_init, control = list(fnscale = -1), thresh_par=thresh_fit_geo_ics$thresh_par, V = gron_eq_cat_lower$V_1[gron_eq_cat_lower$Magnitude > gron_eq_cat_lower$threshold], 
                              ics = gron_eq_cat_lower$ICS[gron_eq_cat_lower$Magnitude > gron_eq_cat_lower$threshold], max_dist = max_dist, min_ics = min_ics)
gpd_ics_V1_fit_lower

# Whole field
excess <- gron_eq_cat$Magnitude[gron_eq_cat$Magnitude > gron_eq_cat$threshold] - gron_eq_cat$threshold[gron_eq_cat$Magnitude > gron_eq_cat$threshold]
pars_init <- c(mean(excess), 0, 0.1)
gpd_ics_V1_fit <- optim(GPD_LL_given_V_ICS, excess = excess, par = pars_init, control = list(fnscale = -1), thresh_par=thresh_fit_geo_ics$thresh_par, V = gron_eq_cat$V_1[gron_eq_cat$Magnitude > gron_eq_cat$threshold], 
                        ics = gron_eq_cat$ICS[gron_eq_cat$Magnitude > gron_eq_cat$threshold], max_dist = max_dist, min_ics = min_ics)
gpd_ics_V1_fit

# Likelihood ratio test
LR_test <- 2 * (gpd_ics_V1_fit_upper$value + gpd_ics_V1_fit_lower$value - gpd_ics_V1_fit$value)
p_value <- 1 - pchisq(LR_test, df = 3, lower.tail = F)
p_value




# Testing for inclusion of ICS --------------------------------------------
# Fitting with ICS above conservative threshold
mags <- gron_eq_cat$Magnitude
conservative_threshold <- rep(1.45, length(mags))
excess_threshold <- conservative_threshold[mags > conservative_threshold]
excesses_conserv <- mags[mags > conservative_threshold] - excess_threshold
ics_conserv <- gron_eq_cat$ICS[mags > conservative_threshold]

(fit_conserv <- optim(GPD_LL_ICS_constant_thresh, par=c(mean(excesses_conserv), 0, 0.1),
                      excess=excesses_conserv, ics=ics_conserv, thresh = excess_threshold, control=list(fnscale=-1)))
# (fit_conserv <- optim(GPD_LL_ICS_constant_thresh_2, par=c(mean(excesses_conserv), 0, 0.1), 
#                      excess=excesses_conserv, ics=ics_conserv, control=list(fnscale=-1)))
sigma_var <- fit_conserv$par[1] + fit_conserv$par[2]*ics_conserv + fit_conserv$par[3]*excess_threshold
transformed_excesses <- transform_to_exp(excesses_conserv, sigma_var, fit_conserv$par[3])
# Quantile regression
# Exp(1) excesses against ics
taus <- seq(0.1, 0.9, by=0.1)
cols <- rainbow(length(taus))
qr_fits <- lapply(taus, function(tau) rq(transformed_excesses ~ ics_conserv, tau=tau))
plot(ics_conserv, transformed_excesses, main = "Conservative excesses on Exp(1) against ICS", xlab="ICS", ylab="Excesses")
for(i in 1:length(qr_fits)) {
  qr_fit <- qr_fits[[i]]
  abline(qr_fit, col=cols[i], lwd=2)
}
legend("topleft", legend=taus, col=cols, lty=1)

# Exp(1) excesses against time
qr_fits <- lapply(taus, function(tau) rq(transformed_excesses ~ as.Date(gron_eq_cat$Date[mags > conservative_threshold]), tau=tau))
plot(as.Date(gron_eq_cat$Date[mags > conservative_threshold]), transformed_excesses, main = "Conservative excesses on Exp(1) against time", xlab="Time", ylab="Excesses")
for(i in 1:length(qr_fits)) {
  qr_fit <- qr_fits[[i]]
  abline(qr_fit, col=cols[i], lwd=2)
}
legend("topleft", legend=taus, col=cols, lty=1)

# QQ-plot
probs <- (1:length(transformed_excesses))/(length(transformed_excesses)+1)
emp_quants <- quantile(transformed_excesses, probs)
theoretical_quants <- qexp(probs, rate=1)

num_boot <- 200
boot_theoretical_quants <- matrix(NA, nrow=num_boot, ncol=length(probs))
for(ii in 1:num_boot) {
  boot_excesses <- rexp(length(excesses_conserv), rate=1)
  boot_theoretical_quants[ii,] <- quantile(boot_excesses, probs)
}
upper_tol <- apply(boot_theoretical_quants, 2, quantile, probs=0.975)
lower_tol <- apply(boot_theoretical_quants, 2, quantile, probs=0.025)

dev.new(height=10, width=10, noRStudioGD = TRUE)
par(mfrow=c(1,1), bg = "transparent")
limits <- range(c(emp_quants, theoretical_quants, upper_tol, lower_tol))
plot(theoretical_quants, emp_quants, xlab="Empirical quantiles", 
     ylab="Theoretical quantiles", xlim=limits, ylim=limits, main = "Conservative threshold")
abline(0,1)
lines(theoretical_quants, upper_tol, lty=2)
lines(theoretical_quants, lower_tol, lty=2)


# Fitting with ICS above stepped threshold
u_h_length <- which(gron_eq_cat$Date == as.Date("2015-01-06"))[1]
piecewise_const_thresh <- c(rep(1.15, u_h_length), rep(0.76, length(mags) - u_h_length))
excess_threshold <- piecewise_const_thresh[mags > piecewise_const_thresh]
excesses_pc <- mags[mags > piecewise_const_thresh] - excess_threshold
ics_pc <- gron_eq_cat$ICS[mags > piecewise_const_thresh]

(fit_pc <- optim(GPD_LL_ICS_constant_thresh, par=c(mean(excesses_pc), 0, 0.1),
                 excess=excesses_pc, ics=ics_pc, thresh = excess_threshold, control=list(fnscale=-1), hessian = T))

sqrt(diag(solve(-1*fit_pc$hessian)))

# (fit_pc <- optim(GPD_LL_ICS_constant_thresh_2, par=c(mean(excesses_pc), 0, 0.1), 
#                       excess=excesses_pc, ics=ics_pc, control=list(fnscale=-1)))

ns_sigma <- fit_pc$par[1] + fit_pc$par[2]*ics_pc + fit_pc$par[3]*excess_threshold
transformed_excesses <- transform_to_exp(excesses_pc, ns_sigma, fit_pc$par[3])

# Quantile regression
# Exp(1) excesses against ics
taus <- seq(0.1, 0.9, by=0.1)
cols <- rainbow(length(taus))
qr_fits <- lapply(taus, function(tau) rq(transformed_excesses ~ ics_pc, tau=tau))
plot(ics_pc, transformed_excesses, main = "Stepped threshold excesses on Exp(1) against ICS", xlab="ICS", ylab="Excesses")
for(i in 1:length(qr_fits)) {
  qr_fit <- qr_fits[[i]]
  abline(qr_fit, col=cols[i], lwd=2)
}
legend("topleft", legend=taus, col=cols, lty=1)

# Exp(1) excesses against time
qr_fits <- lapply(taus, function(tau) rq(transformed_excesses ~ as.Date(gron_eq_cat$Date[mags > piecewise_const_thresh]), tau=tau))
plot(as.Date(gron_eq_cat$Date[mags > piecewise_const_thresh]), transformed_excesses, main = "Stepped threshold excesses on Exp(1) against time", xlab="Time", ylab="Excesses")
for(i in 1:length(qr_fits)) {
  qr_fit <- qr_fits[[i]]
  abline(qr_fit, col=cols[i], lwd=2)
}
legend("topleft", legend=taus, col=cols, lty=1)

# QQ-plot
probs <- (1:length(transformed_excesses))/(length(transformed_excesses)+1)
emp_quants <- quantile(transformed_excesses, probs)
theoretical_quants <- qexp(probs, rate=1)

num_boot <- 200
boot_theoretical_quants <- matrix(NA, nrow=num_boot, ncol=length(probs))
for(ii in 1:num_boot) {
  boot_excesses <- rexp(length(excesses_pc), rate=1)
  boot_theoretical_quants[ii,] <- quantile(boot_excesses, probs)
}
upper_tol <- apply(boot_theoretical_quants, 2, quantile, probs=0.975)
lower_tol <- apply(boot_theoretical_quants, 2, quantile, probs=0.025)

dev.new(height=10, width=30, noRStudioGD = TRUE)
par(mfrow=c(1,3), bg = "transparent")
limits <- range(c(emp_quants, theoretical_quants, upper_tol, lower_tol))
plot(theoretical_quants, emp_quants, xlab="Empirical quantiles", 
     ylab="Theoretical quantiles", xlim=limits, ylim=limits, main = "Stepped threshold (Whole)")
abline(0,1)
lines(theoretical_quants, upper_tol, lty=2)
lines(theoretical_quants, lower_tol, lty=2)

# Fitting separately in each time region using stepped threshold
mags_before <- gron_eq_cat$Magnitude[1:u_h_length]
mags_after <- gron_eq_cat$Magnitude[(u_h_length+1):length(mags)]
ics_before <- gron_eq_cat$ICS[1:u_h_length]
ics_after <- gron_eq_cat$ICS[(u_h_length+1):length(mags)]

threshold_before <- piecewise_const_thresh[1]
threshold_after <- piecewise_const_thresh[u_h_length+1]

excesses_before <- mags_before[mags_before > threshold_before] - threshold_before
excesses_after <- mags_after[mags_after > threshold_after] - threshold_after

ics_before_excess <- ics_before[mags_before > threshold_before]
ics_after_excess <- ics_after[mags_after > threshold_after]

(fit_before <- optim(GPD_LL_ICS_constant_thresh, par=c(mean(excesses_before), 0, 0.1),
                     excess=excesses_before, ics=ics_before_excess, thresh = threshold_before, control=list(fnscale=-1), hessian = T))

sqrt(diag(solve(-1*fit_before$hessian)))

(fit_after <- optim(GPD_LL_ICS_constant_thresh, par=c(mean(excesses_after), 0, 0.1),
                    excess=excesses_after, ics=ics_after_excess, thresh = threshold_after, control=list(fnscale=-1), hessian = T))

sqrt(diag(solve(-1*fit_after$hessian)))

# Log-likelihood ratio test
LL_before_after <- -1*(fit_before$value + fit_after$value)
LL_all <- -1*(fit_pc$value)
test_stat <- 2*(LL_all - LL_before_after)
(p_val <- 1-pchisq(test_stat, 3))

# QQ-plots
#Before
ns_sigma_before <- fit_before$par[1] + fit_before$par[2]*ics_before_excess + fit_before$par[3]*threshold_before
transformed_excesses_before <- transform_to_exp(excesses_before, ns_sigma_before, fit_before$par[3])
probs <- (1:length(transformed_excesses_before))/(length(transformed_excesses_before)+1)
emp_quants_before <- quantile(transformed_excesses_before, probs)
theoretical_quants_before <- qexp(probs, rate=1)

num_boot <- 200
boot_theoretical_quants_before <- matrix(NA, nrow=num_boot, ncol=length(probs))
for(ii in 1:num_boot) {
  boot_excesses <- rexp(length(excesses_before), rate=1)
  boot_theoretical_quants_before[ii,] <- quantile(boot_excesses, probs)
}
upper_tol_before <- apply(boot_theoretical_quants_before, 2, quantile, probs=0.975)
lower_tol_before <- apply(boot_theoretical_quants_before, 2, quantile, probs=0.025)

limits <- range(c(emp_quants_before, theoretical_quants_before, upper_tol_before, lower_tol_before))
plot(theoretical_quants_before, emp_quants_before, xlab="Empirical quantiles", 
     ylab="Theoretical quantiles", xlim=limits, ylim=limits, main = "Stepped Threshold (Before)")
abline(0,1)
lines(theoretical_quants_before, upper_tol_before, lty=2)
lines(theoretical_quants_before, lower_tol_before, lty=2)

#After
ns_sigma_after <- fit_after$par[1] + fit_after$par[2]*ics_after_excess + fit_after$par[3]*threshold_after
transformed_excesses_after <- transform_to_exp(excesses_after, ns_sigma_after, fit_after$par[3])
probs <- (1:length(transformed_excesses_after))/(length(transformed_excesses_after)+1)
emp_quants_after <- quantile(transformed_excesses_after, probs)
theoretical_quants_after <- qexp(probs, rate=1)

num_boot <- 200
boot_theoretical_quants_after <- matrix(NA, nrow=num_boot, ncol=length(probs))
for(ii in 1:num_boot) {
  boot_excesses <- rexp(length(excesses_after), rate=1)
  boot_theoretical_quants_after[ii,] <- quantile(boot_excesses, probs)
}
upper_tol_after <- apply(boot_theoretical_quants_after, 2, quantile, probs=0.975)
lower_tol_after <- apply(boot_theoretical_quants_after, 2, quantile, probs=0.025)

limits <- range(c(emp_quants_after, theoretical_quants_after, upper_tol_after, lower_tol_after))
plot(theoretical_quants_after, emp_quants_after, xlab="Empirical quantiles", 
     ylab="Theoretical quantiles", xlim=limits, ylim=limits, main = "Stepped Threshold (After)")
abline(0,1)
lines(theoretical_quants_after, upper_tol_after, lty=2)
lines(theoretical_quants_after, lower_tol_after, lty=2)



# Repeat above with shape constrained to be constant
(fit_before <- optim(GPD_LL_ICS_constant_thresh_const_shape, par=c(mean(excesses_before), 0),
                     excess=excesses_before, ics=ics_before_excess, thresh = threshold_before, shape=-0.07, control=list(fnscale=-1), hessian = T))
sqrt(diag(solve(-1*fit_before$hessian)))

(fit_after <- optim(GPD_LL_ICS_constant_thresh_const_shape, par=c(mean(excesses_after), 0),
                    excess=excesses_after, ics=ics_after_excess, thresh = threshold_after, shape=-0.07, control=list(fnscale=-1)))



# Investigating Infs
lens <- numeric(1108)
shapes <- numeric(1108)
for(i in 1:1108){
  lens[i] <- length(thresh_fit_geo_ics$inf_list_mean[[i]]$excess)
  shapes[i] <- thresh_fit_geo_ics$inf_list_mean[[i]]$par_ests[3]
}
hist(lens)
hist(shapes)

transform_to_exp(y = thresh_fit_geo_ics$inf_list_boot[[1]]$excess_boot, sig = thresh_fit_geo_ics$inf_list_boot[[1]]$sigma_var, xi = thresh_fit_geo_ics$inf_list_boot[[1]]$est[3])
(1/thresh_fit_geo_ics$inf_list_boot[[1]]$est[3])*log(1 + thresh_fit_geo_ics$inf_list_boot[[1]]$est[3] * thresh_fit_geo_ics$inf_list_boot[[1]]$excess_boot / thresh_fit_geo_ics$inf_list_boot[[1]]$sigma_var)


# Threshold selection with ICS and V1
# Choosing max and min values for distance and ics
covariates <- read.csv("Data/covariates/covariates_1995-2024.csv", header=T)
max_dist <- max(covariates$V1)
min_ics <- min(covariates$ICS)

set.seed(11111)
thresh_fit_geo_ics <- eqd_geo_ics(data=mags, thresh = threshold_matrix, distance_to_geo = gron_eq_cat$V_1, k=num_boot, ics = gron_eq_cat$ICS, max_dist = max_dist, min_ics = min_ics)
saveRDS(thresh_fit_geo_ics, paste0("threshold_results/thresh_fit_V1_ics.rds"))

thresh_fit_geo_ics <- readRDS("threshold_results/thresh_fit_V1_ics.rds")


# Resulting fit

chosen_threshold <- thresh_fit_geo_ics$thresh[1] + thresh_fit_geo_ics$thresh[2]*gron_eq_cat$V_1
excesses <- mags[mags > chosen_threshold] - chosen_threshold[mags > chosen_threshold]
ics_excess <- gron_eq_cat$ICS[mags > chosen_threshold]
ns_sigma <- thresh_fit_geo_ics$par[1] + thresh_fit_geo_ics$par[2]*ics_excess + thresh_fit_geo_ics$par[3]*chosen_threshold[mags > chosen_threshold]
transformed_excesses <- transform_to_exp(excesses, ns_sigma, thresh_fit_geo_ics$par[3])

# Quantile regression
# Exp(1) excesses against ics
taus <- seq(0.1, 0.9, by=0.1)
cols <- rainbow(length(taus))
qr_fits <- lapply(taus, function(tau) rq(transformed_excesses ~ ics_excess, tau=tau))
plot(ics_excess, transformed_excesses, main = "V1 excesses on Exp(1) against ICS", xlab="ICS", ylab="Excesses")
for(i in 1:length(qr_fits)) {
  qr_fit <- qr_fits[[i]]
  abline(qr_fit, col=cols[i], lwd=2)
}
legend("topleft", legend=taus, col=cols, lty=1)

# Exp(1) excesses against time
qr_fits <- lapply(taus, function(tau) rq(transformed_excesses ~ as.Date(gron_eq_cat$Date[mags > chosen_threshold]), tau=tau))
plot(as.Date(gron_eq_cat$Date[mags > chosen_threshold]), transformed_excesses, main = "V1 excesses on Exp(1) against time", xlab="Time", ylab="Excesses")
for(i in 1:length(qr_fits)) {
  qr_fit <- qr_fits[[i]]
  abline(qr_fit, col=cols[i], lwd=2)
}
legend("topleft", legend=taus, col=cols, lty=1)

# QQ-plot
probs <- (1:length(transformed_excesses))/(length(transformed_excesses)+1)
emp_quants <- quantile(transformed_excesses, probs)
theoretical_quants <- qexp(probs, rate=1)

num_boot <- 200
boot_theoretical_quants <- matrix(NA, nrow=num_boot, ncol=length(probs))
for(ii in 1:num_boot) {
  boot_excesses <- rexp(length(excesses), rate=1)
  boot_theoretical_quants[ii,] <- quantile(boot_excesses, probs)
}
upper <- apply(boot_theoretical_quants, 2, quantile, probs=0.975)
lower <- apply(boot_theoretical_quants, 2, quantile, probs=0.025)

dev.new(height=10, width=10, noRStudioGD = TRUE)
par(mfrow=c(1,1), bg = "transparent")
limits <- range(c(emp_quants, theoretical_quants, upper, lower))
plot(theoretical_quants, emp_quants, xlab="Empirical quantiles", 
     ylab="Theoretical quantiles", xlim=limits, ylim=limits, main = "V1 and ICS")
abline(0,1)
lines(theoretical_quants, upper, lty=2)
lines(theoretical_quants, lower, lty=2)

