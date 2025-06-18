
gron_eq_cat <- read.csv("Data/Events/unrounded_after_1995_in_polygon_with_covariates.csv", header=T)
covariates <- read.csv("Data/covariates/covariates_1995-2024.csv", header=T)
covariates_2055 <- read.csv("Data/covariates/covariates_1995-2055.csv", header=T)
geophones_deepest <- read.csv("Data/Geophones/Geophones_processed_03-07-2024_deepest_only.csv", header=T, row.names=1)
gron_outline <- read.csv('Data/Geophones/Groningen_Field_outline.csv', header=T)
gron_polygon <- read.table('Data/Geophones/polygon_for_groningen_earthquakes.txt', header=T)
gron_rect <- data.frame(X=c(210000,275000, 275000, 210000, 210000), Y=c(560000, 560000, 625000, 625000, 560000))

library(ggplot2)
library(ggspatial)
library(pracma)
library(cowplot)
library(dplyr)
library(purrr)
library(patchwork)

source("src/helper_functions.R")
source("intensity_estimation.R")

# TODO: Need to add resulting thresholds to covariates for ICS max

# Max dists to use in any fits
max_dists <- c(max(covariates$V1), max(covariates$V2), max(covariates$V3), max(covariates$V4, na.rm = TRUE))

# Zak's sigmoid threshold
sigmoid_threshold <- function(x, vl = 1.15, vr = 0.76, mu, zeta){
  return(vr + (vl - vr) * pnorm((mu - x) / zeta))
}

change_index <- which(gron_eq_cat$Date == as.Date("2015-12-25"))
pc_threshold <- c(rep(1.15, change_index), rep(0.76, nrow(gron_eq_cat)-change_index))

dev.new(height=5, width=10, noRStudioGD = TRUE)
par(mfrow=c(1,2), bg='transparent')
plot(as.Date(gron_eq_cat$Date),gron_eq_cat$Magnitude, xlab="Time", ylab = expression(Magnitude~(M[L])), pch=19, col="grey", cex=0.7)
lines(as.Date(gron_eq_cat$Date),rep(1.45, nrow(gron_eq_cat)), col="red", lty=2, lwd=2)
#gron_eq_cat$zak_sigmoid <- sigmoid_threshold(x=c(1:nrow(gron_eq_cat)), mu = 746, zeta = 1)
lines(as.Date(gron_eq_cat$Date), pc_threshold, col="blue", lwd=2)

plot(gron_polygon$POINT_X, gron_polygon$POINT_Y, col="green", xlab="Easting (m)", ylab = "Northing (m)", type = 'l', lty=2, asp=1, lwd=2)
points(gron_eq_cat$Easting, gron_eq_cat$Northing, pch=19, col="grey", cex=0.7, asp=1)
lines(gron_outline$X, gron_outline$Y, col="black", asp=1)

#Rates of evets below conservtiave
mags_h <- gron_eq_cat$Magnitude[1:change_index]
(p_h <- length(mags_h[mags_h <= 1.45]) / length(mags_h))
(se_h <- sqrt(p_h * (1 - p_h) / length(mags_h)))
mags_l <- gron_eq_cat$Magnitude[(change_index+1):nrow(gron_eq_cat)]
(p_l <- length(mags_l[mags_l <= 1.45]) / length(mags_l))
(se_l <- sqrt(p_l * (1 - p_l) / length(mags_l)))

(rate_per_year_h <- length(mags_h)/19.75)
(rate_per_year_l <- length(mags_l)/8)
(num_above_h <- length(mags_h[mags_h > 1.45]))
(num_above_l <- length(mags_l[mags_l > 1.45]))


# Average max stress over time (including future)
average_ICS_df <- covariates_2055 %>%
  group_by(Date) %>%
  summarise(average_ICS = mean(ICS_max, na.rm = TRUE), .groups = "drop")

saveRDS(average_ICS_df, file = "Data/covariates/average_ICS_max_1995-2055.rds")

#Adding three locations
locations <- data.frame(Easting = c(250000, 250000, 250000), Northing = c(575000, 590000, 605000))
ics_at_location_1 <- covariates_2055$ICS_max[covariates_2055$Easting == locations$Easting[1] & covariates_2055$Northing == locations$Northing[1]]
ics_at_location_2 <- covariates_2055$ICS_max[covariates_2055$Easting == locations$Easting[2] & covariates_2055$Northing == locations$Northing[2]]
ics_at_location_3 <- covariates_2055$ICS_max[covariates_2055$Easting == locations$Easting[3] & covariates_2055$Northing == locations$Northing[3]]

dev.new(height=5, width=5, noRStudioGD = TRUE)
par(mfrow=c(1,1), bg='transparent')
plot(dates, average_ICS_df$average_ICS, xlab = "Time", ylab = "Average KS (MPa)", type = 'l', lwd=2, ylim = c(0, 0.4), col="black")
lines(dates, ics_at_location_1, col="blue", lwd=2)
lines(dates, ics_at_location_2, col="green", lwd=2)
lines(dates, ics_at_location_3, col="red", lwd=2)
abline(v=as.Date("2015-12-25"), col="grey", lty=2)
abline(v=as.Date(max(gron_eq_cat$Date)), col="grey", lty=2)
# Average Kaiser stress plots for different years -------------------------

dev.new(height=5, width=5, noRStudioGD = TRUE)
par(mfrow=c(1,1), bg='transparent')
year <- c("2020")
current_covariates <- covariates[covariates$Year == year,]
plot_df <- data.frame(Easting = numeric(0), Northing = numeric(0), avg_stress = numeric(0))
  for(easting in unique(current_covariates$Easting)){
    for(northing in unique(current_covariates$Northing)){
      plot_df <- rbind(plot_df, data.frame(Easting=easting, Northing = northing, Average_ICS = mean(current_covariates$ICS_max[current_covariates$Easting == easting & current_covariates$Northing == northing])))
    }
  }
plot_df <- plot_df[complete.cases(plot_df),]

ggplot(plot_df, aes(x = Easting, y = Northing, fill = Average_ICS)) + geom_tile() + scale_fill_gradient(low = "red", high = "yellow") + fixed_plot_aspect(ratio = 1) + theme_classic() +
    theme(plot.background = element_blank()) + geom_point(data=gron_outline, aes(x = X, y = Y), size=0.5, shape=1, fill = "black")  +
    geom_point(data=locations, aes(x = Easting, y = Northing), size=2, shape=19, fill = "black") +
    labs(x = "Easting (m)", y = "Northing (m)", fill = "Average KS") + coord_fixed() 


# Filter data for the selected year
current_covariates <- covariates %>%
  filter(Year == "2020")  # Assuming Year is character; if it's numeric, use 2020

# Compute average ICS_max per grid cell
plot_df <- current_covariates %>%
  group_by(Easting, Northing) %>%
  summarise(Average_ICS = mean(ICS_max, na.rm = TRUE), .groups = "drop") %>%
  filter(complete.cases(.))

# Plot
ggplot(plot_df, aes(x = Easting, y = Northing, fill = Average_ICS)) +
  geom_tile() +
  scale_fill_gradient(low = "red", high = "yellow") +
  coord_fixed() +
  theme_classic() +
  theme(plot.background = element_blank()) +
  geom_point(data = gron_outline, aes(x = X, y = Y), size = 0.5, shape = 1, fill = "black") +
  geom_point(data = locations, aes(x = Easting, y = Northing), size = 2, shape = 19, fill = "black") +
  labs(x = "Easting (m)", y = "Northing (m)", fill = "Average KS")

# Boxplots of magnitudes for high/low Kaiser stress for exceedances of conservative mc
gron_eq_cat_ex_conserv <- gron_eq_cat[gron_eq_cat$Magnitude > 1.45,]
mags_high_stress_conserv <- gron_eq_cat_ex_conserv$Magnitude[gron_eq_cat_ex_conserv$ICS_max > median(gron_eq_cat_ex_conserv$ICS_max)]
mags_low_stress_conserv <- gron_eq_cat_ex_conserv$Magnitude[gron_eq_cat_ex_conserv$ICS_max <= median(gron_eq_cat_ex_conserv$ICS_max)]
boxplot(mags_low_stress_conserv, mags_high_stress_conserv, names=c("Low KS", "High KS"), ylab = expression(Magnitude~(M[L])))


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



# Plots of best thresholds ------------------------------------------------

# ICS max
thresh_fit_A2 <- readRDS("threshold_results/icsmax/geo_thresh_fit_V2.rds")

thresh_fit_B1 <- readRDS("threshold_results/icsmax/geo_thresh_fit_logV1.rds")

thresh_fit_C2 <- readRDS("threshold_results/icsmax/geo_thresh_fit_sqrtV2.rds")


# Observed thresholds
gron_eq_cat$threshold_A2 <- thresh_fit_A2$thresh_par[1] + thresh_fit_A2$thresh_par[2] * gron_eq_cat$V_2
gron_eq_cat$threshold_B1 <- thresh_fit_B1$thresh_par[1] + thresh_fit_B1$thresh_par[2] * log(gron_eq_cat$V_1)
gron_eq_cat$threshold_C2 <- thresh_fit_C2$thresh_par[1] + thresh_fit_C2$thresh_par[2] * sqrt(gron_eq_cat$V_2)

gron_eq_cat$endpoint_A2 <- -(thresh_fit_A2$par[1] + thresh_fit_A2$par[2] * gron_eq_cat$ICS_max) / thresh_fit_A2$par[3]
gron_eq_cat$endpoint_B1 <- -(thresh_fit_B1$par[1] + thresh_fit_B1$par[2] * gron_eq_cat$ICS_max) / thresh_fit_B1$par[3]
gron_eq_cat$endpoint_C2 <- -(thresh_fit_C2$par[1] + thresh_fit_C2$par[2] * gron_eq_cat$ICS_max) / thresh_fit_C2$par[3]

gron_eq_cat[gron_eq_cat$Magnitude > 3.5,]

gron_eq_cat_exceed_A2 <- gron_eq_cat[gron_eq_cat$Magnitude > gron_eq_cat$threshold_A2,] 

sum(gron_eq_cat$Magnitude > gron_eq_cat$threshold_A2)
sum(gron_eq_cat$Magnitude > gron_eq_cat$threshold_B1)
sum(gron_eq_cat$Magnitude > gron_eq_cat$threshold_C2)

dev.new(height=5, width=5, noRStudioGD = TRUE)
par(mfrow=c(1,1), bg='transparent')
plot(as.Date(gron_eq_cat$Date), gron_eq_cat$Magnitude, xlab="Event time", ylab = "Magnitude", pch=19, col="grey", cex=0.7)
lines(as.Date(gron_eq_cat$Date), gron_eq_cat$threshold_A2, lwd=2)
lines(as.Date(gron_eq_cat$Date),rep(1.45, nrow(gron_eq_cat)), col="red", lty=2, lwd=2)
#gron_eq_cat$zak_sigmoid <- sigmoid_threshold(x=c(1:nrow(gron_eq_cat)), mu = 746, zeta = 1)
lines(as.Date(gron_eq_cat$Date), pc_threshold, col="blue", lwd=2)

diff_thr_B1 <- gron_eq_cat$threshold_A2 - gron_eq_cat$threshold_B1
diff_thr_C2 <- gron_eq_cat$threshold_A2 - gron_eq_cat$threshold_C2

plot(as.Date(gron_eq_cat$Date), diff_thr_B1, xlab="Event time", ylab = "Difference in threshold", col="green", type='l', 
     ylim=c(min(diff_thr_B1, diff_thr_C2), max(diff_thr_B1, diff_thr_C2)))
lines(as.Date(gron_eq_cat$Date), diff_thr_C2, col="red")
abline(h=0, col="black", lty=2)


# Goodness of fit ---------------------------------------------------------

#QQplot
dev.new(height=5, width=5, noRStudioGD = TRUE)
par(mfrow=c(1,1), bg='transparent')
get_qq_plot_const(gron_eq_cat$Magnitude, 1.45, main="" )
get_qq_plot_geo_ics(gron_eq_cat$Magnitude, thresh_fit_A2, gron_eq_cat$V_2, gron_eq_cat$ICS_max, main="" )

#PPplots
dev.new(height=5, width=5, noRStudioGD = TRUE)
par(mfrow=c(1,1), bg='transparent')

get_pp_plot_const(gron_eq_cat$Magnitude, 1.45, main="", n_boot=1000 )
get_pp_plot_geo_ics(gron_eq_cat$Magnitude, thresh_fit_A2, gron_eq_cat$V_2, gron_eq_cat$ICS_max, main="", n_boot=1000 )

# Spatial threshold plots --------------------------------------------------

dev.new(height=5, width=10, noRStudioGD = TRUE)
par(mfrow=c(1,1), bg='transparent')

# Ensure proper date conversion
covariates$Date <- as.Date(covariates$Date)
geophones_deepest$Start_date <- as.Date(geophones_deepest$Start_date)
geophones_deepest$End_date <- as.Date(geophones_deepest$End_date)
chosen_dates <- as.Date(c("2010-01-01", "2020-01-01"))  

covariates_for_dates <- covariates %>% filter(Date %in% chosen_dates)
# Ensure fill scale range is consistent
fill_limits <- range(covariates_for_dates$best_threshold, na.rm = TRUE)

# Function to create plot for a given date
plot_threshold_for_date <- function(date) {
  current_covariates <- covariates %>% filter(Date == date)
  current_geophones <- geophones_deepest %>%
    filter(Start_date <= date, End_date >= date)
  current_geo_in_polygon <- current_geophones %>%
    filter(inpolygon(Xcoord, Ycoord, gron_polygon$POINT_X, gron_polygon$POINT_Y))
  
  ggplot(current_covariates, aes(x = Easting, y = Northing, fill = best_threshold)) +
    geom_tile() + fixed_plot_aspect(ratio = 1) +
    scale_fill_gradient(low = "blue", high = "red", limits = fill_limits) +
    coord_fixed() +
    theme_classic() +
    geom_point(data = current_geo_in_polygon, aes(x = Xcoord, y = Ycoord), size = 1, shape = 19, fill = "black") +
    labs(fill = "Threshold")
}

# Generate plots
plots <- lapply(chosen_dates, plot_threshold_for_date)

# Confirm that both are valid ggplot objects
if (all(sapply(plots, inherits, "ggplot"))) {
  combined_plot <- wrap_plots(plots, guides = "collect") & theme(legend.position = "right")
  print(combined_plot)
} else {
  stop("One or more plots are not ggplot objects.")
}


# Poisson process intensity fit -------------------------------------------

gron_eq_cat_exceed_A2 <- gron_eq_cat[gron_eq_cat$Magnitude > gron_eq_cat$threshold_A2,] 

sum(gron_eq_cat_exceed_A2$dsmaxdt == 0)
# sum(gron_eq_cat_exceed_A2$dsmaxdt <= 1e-5)
# sum(gron_eq_cat_exceed_A2$dsmaxdt <= 1e-4)
# sum(gron_eq_cat_exceed_A2$dsmaxdt <= 5.5e-5)

# Likelihood was adjusted each time for above limits and rerun below to obtain 
# corresponding par ests

(opt_PP_LL_icsmax <- optim(Poisson_process_LL_icsmax, par = c(0.1, 0), data = gron_eq_cat, covariates = covariates, 
                           threshold_obs = gron_eq_cat$threshold_A2, covariates_threshold = covariates$best_threshold , 
                           thresh_fit = thresh_fit_A2, control=list(fnscale=-1), hessian=T))

(se <- sqrt(diag(solve(-opt_PP_LL_icsmax$hessian))))
(CIs <- cbind(opt_PP_LL_icsmax$par - qnorm(0.975) * se, opt_PP_LL_icsmax$par + qnorm(0.975) * se))

# Corresponding par ests if using the intensity above conservative threshold
# as in Bourne 2018
#Checking how many events are removed by dsmaxdt condition
gron_eq_cat_exceed_1.45 <- gron_eq_cat[gron_eq_cat$Magnitude > 1.45,]
sum(gron_eq_cat_exceed_1.45$dsmaxdt == 0)

(opt_PP_LL_icsmax_conserv <- optim(Poisson_process_LL_const_thresh, par = c(0.1, 0), data = gron_eq_cat, covariates = covariates, 
                                     threshold = 1.45, control=list(fnscale=-1), hessian=T))

(se <- sqrt(diag(solve(-opt_PP_LL_icsmax_conserv$hessian))))
(CIs <- cbind(opt_PP_LL_icsmax_conserv$par - qnorm(0.975) * se, opt_PP_LL_icsmax_conserv$par + qnorm(0.975) * se))


# Resulting intensities
covariates$intensity_above_threshold_A2 <- resulting_intensity_icsmax(opt_PP_LL_icsmax, covariates, covariates$best_threshold, thresh_fit_A2)
covariates$intensity_above_0 <- resulting_intensity_icsmax(opt_PP_LL_icsmax, covariates, 0, thresh_fit_A2)

grid_box_E <- (max(unique(covariates$Easting))-min(unique(covariates$Easting)))/length(unique(covariates$Easting))
grid_box_N <- (max(unique(covariates$Northing))-min(unique(covariates$Northing)))/length(unique(covariates$Northing))
grid_box_area <- grid_box_E/1000* grid_box_N/1000

covariates$lambda_A2_per_km2 <- covariates$intensity_above_threshold_A2/grid_box_area


# Aggregated intensity plots for different years -----------------------------

# NOTE: Below filters out exceedances corresponding to dsmaxdt = 0 to make sure plots
# are comparing observed values which correspond to how intensity was estimated

dev.new(height=5, width=10, noRStudioGD = TRUE)
par(mfrow=c(1,1), bg='transparent')
chosen_years <- c("2010", "2020")

covariates_for_2010 <- covariates %>% filter(Year == "2010") %>% group_by(Easting, Northing) %>%
  summarise(agg_intensity = sum(lambda_A2_per_km2, na.rm = TRUE), .groups = "drop") %>%
  filter(!is.na(agg_intensity))

covariates_for_2020 <- covariates %>% filter(Year == "2020") %>% group_by(Easting, Northing) %>%
  summarise(agg_intensity = sum(lambda_A2_per_km2, na.rm = TRUE), .groups = "drop") %>%
  filter(!is.na(agg_intensity))

fill_limits <- range(c(covariates_for_2010$agg_intensity, covariates_for_2020$agg_intensity), na.rm = TRUE)

plot_intensity_for_year <- function(year) {
  current_covariates <- covariates %>% filter(Year == year)
  current_exceedances <- gron_eq_cat_exceed_A2 %>% filter(Year == year) %>% filter(dsmaxdt > 0)
  
  plot_df <- current_covariates %>%
    group_by(Easting, Northing) %>%
    summarise(agg_intensity = sum(lambda_A2_per_km2, na.rm = TRUE), .groups = "drop") %>%
    filter(!is.na(agg_intensity))
  
  ggplot(plot_df, aes(x = Easting, y = Northing, fill = agg_intensity)) +
    geom_tile() + fixed_plot_aspect(ratio = 1) + theme_classic() +
    theme(plot.background = element_blank()) + scale_fill_gradient(low = "blue", high = "red", limits = fill_limits) +
    geom_point(data = current_exceedances, aes(x = Easting, y = Northing), size=1, shape=19, fill = "black")  + 
    labs(x = "Easting (m)", y = "Northing (m)", fill = expression(Lambda[u])) + coord_fixed() 
}

# Generate plots
plots <- lapply(chosen_years, plot_intensity_for_year)

# Confirm that both are valid ggplot objects
if (all(sapply(plots, inherits, "ggplot"))) {
  combined_plot <- wrap_plots(plots, guides = "collect") & theme(legend.position = "right")
  print(combined_plot)
} else {
  stop("One or more plots are not ggplot objects.")
}
# Aggregated intensity over space for all years ----------------------------------------

dev.new(height=5, width=5, noRStudioGD = TRUE)
par(mfrow=c(1,1), bg='transparent')

#Aggregated intensity by year
agg_intensity_df <- covariates %>%
  group_by(Year) %>%
  summarise(agg_intensity = sum(intensity_above_threshold_A2, na.rm = TRUE), .groups = "drop")

# Remove last year
agg_intensity_df <- agg_intensity_df[agg_intensity_df$Year != max(agg_intensity_df$Year),]

# Compute number of exceedances for each year
# NOTE: Below filters out exceedances corresponding to dsmaxdt =0 to make sure plots
# are comparing observed values which correspond to how intensity was estimated

num_exceedances_per_year <- gron_eq_cat_exceed_A2 %>% filter(dsmaxdt > 0) %>%
  group_by(Year) %>%  summarise(num_exceedances = n(), .groups = "drop")

# Convert Year to numeric for merging
agg_intensity_df$Year <- as.numeric(as.character(agg_intensity_df$Year))

# Merge the two data frames
agg_intensity_df <- agg_intensity_df %>%
  left_join(num_exceedances_per_year, by = "Year")
#Adding zeros for number of exceedances in first two years
agg_intensity_df$num_exceedances[is.na(agg_intensity_df$num_exceedances)] <- 0

# Plot agg intensity and number of exceedances against year
ggplot(agg_intensity_df, aes(x = Year)) +
  geom_point(aes(y = agg_intensity*grid_box_area, color = "Aggregated Intensity"), size = 2) +
  geom_line(aes(y = agg_intensity*grid_box_area, color = "Aggregated Intensity"), size = 1) +
  geom_point(aes(y = num_exceedances , color = "Number of Exceedances"), size = 2) +
  geom_line(aes(y = num_exceedances, color = "Number of Exceedances"), size = 1) +
  labs(x = "Year", y = "Number per year") +
  theme_classic() +
  theme(plot.background = element_blank()) +
  scale_color_manual(values = c("blue", "red"), guide="none")


# Aggregated intensity above 0 for all years ----------------------------------------
dev.new(height=5, width=5, noRStudioGD = TRUE)
par(mfrow=c(1,1), bg='transparent')

#Aggregated intensity by year
agg_intensity_df <- covariates %>%
  group_by(Year) %>%
  summarise(agg_intensity = sum(intensity_above_0, na.rm = TRUE), .groups = "drop")

# Remove last year
agg_intensity_df <- agg_intensity_df[agg_intensity_df$Year != max(agg_intensity_df$Year),]

# Compute number of exceedances for each year
# NOTE: Below filters out exceedances corresponding to dsmaxdt =0 to make sure plots
# are comparing observed values which correspond to how intensity was estimated

gron_eq_cat_exceed_0 <- gron_eq_cat[gron_eq_cat$Magnitude > 0,]
num_exceedances_per_year <- gron_eq_cat_exceed_0 %>% filter(dsmaxdt > 0) %>%
  group_by(Year) %>%
  summarise(num_exceedances = n(), .groups = "drop")

# Convert Year to numeric for merging
agg_intensity_df$Year <- as.numeric(as.character(agg_intensity_df$Year))

# Merge the two data frames
agg_intensity_df <- agg_intensity_df %>%
  left_join(num_exceedances_per_year, by = "Year")

# Plot agg intensity and number of exceedances against year
ggplot(agg_intensity_df, aes(x = Year)) +
  geom_point(aes(y = agg_intensity*grid_box_area, color = "Aggregated Intensity"), size = 2) +
  geom_line(aes(y = agg_intensity*grid_box_area, color = "Aggregated Intensity"), size = 1) +
  geom_point(aes(y = num_exceedances , color = "Number of Exceedances"), size = 2) +
  geom_line(aes(y = num_exceedances, color = "Number of Exceedances"), size = 1) +
  labs(x = "Year", y = "Number per year") +
  theme_classic() +
  theme(plot.background = element_blank()) +
  scale_color_manual(values = c("blue", "red"), guide="none")

# Proportion observed relative to expected (above 0)
ggplot(agg_intensity_df, aes(x = Year)) +
  geom_point(aes(y = num_exceedances/(agg_intensity*grid_box_area), color = "Proportion observed"), size = 2) +
  geom_line(aes(y = num_exceedances/(agg_intensity*grid_box_area), color = "Proportion observed"), size = 1) +
  labs(x = "Year", y = "Proportion observed") +
  theme_classic() +
  theme(plot.background = element_blank()) +
  scale_color_discrete(guide="none")


