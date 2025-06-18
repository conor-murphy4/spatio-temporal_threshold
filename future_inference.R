
source("src/distance_to_nearest_geo.R")
source("intensity_estimation.R")

library(dplyr)
library(pracma)
library(ggplot2)
library(ggspatial)
library(pracma)
library(cowplot)
library(dplyr)
library(purrr)
library(patchwork)

gron_eq_cat <- read.csv("Data/Events/unrounded_after_1995_in_polygon_with_covariates.csv", header=T)
geophones_deepest <- read.csv("Data/Geophones/Geophones_processed_03-07-2024_deepest_only.csv", header=T, row.names = 1)
covariates_2055 <- read.csv("Data/covariates/covariates_1995-2055.csv", header=T)
covariates <- read.csv("Data/covariates/covariates_1995-2024.csv", header=T)
future_covariates <- read.csv("Data/covariates/future_covariates_2024-2055.csv", header=T)
gron_outline <- read.csv('Data/Geophones/Groningen_Field_outline.csv', header=T)

# Future period of 30 years from 01-01-2025 to 01-01-2055
# future_covariates <- future_covariates %>%
#  filter(Date >= "2025-01-01" ) 
# 
# future_covariates$Year <- stringr::str_sub(future_covariates$Date, 1, 4) %>% as.numeric()

# future_covariates$V2 <- distance_to_nearest(future_covariates, geophones_deepest,2)
# future_covariates$V3 <- distance_to_nearest(future_covariates, geophones_deepest,3)
# Best threshold fit
thresh_fit_A2 <- readRDS("threshold_results/icsmax/geo_thresh_fit_V2.rds")
gron_eq_cat$threshold_A2 <- thresh_fit_A2$thresh_par[1] + thresh_fit_A2$thresh_par[2] * gron_eq_cat$V_2


# Poisson process fit on observed data
(opt_PP_LL_icsmax <- optim(Poisson_process_LL_icsmax, par = c(0.1, 0), data = gron_eq_cat, covariates = covariates, 
                           threshold_obs = gron_eq_cat$threshold_A2, covariates_threshold = covariates$best_threshold , 
                           thresh_fit = thresh_fit_A2, control=list(fnscale=-1), hessian=T))


# future_covariates$intensity_above_0 <- resulting_intensity_icsmax(opt_PP_LL_icsmax, future_covariates, 0, thresh_fit_A2)
# future_covariates$sigma_0 <- thresh_fit_A2$par[1] + thresh_fit_A2$par[2] * future_covariates$ICS_max

# write.csv(future_covariates, "Data/covariates/future_covariates_2024-2055.csv", row.names = FALSE)

# These values will mostly be zero or close to zero due to the stresses not increasing

grid_box_E <- (max(unique(covariates$Easting))-min(unique(covariates$Easting)))/length(unique(covariates$Easting))
grid_box_N <- (max(unique(covariates$Northing))-min(unique(covariates$Northing)))/length(unique(covariates$Northing))
grid_box_area <- grid_box_E/1000* grid_box_N/1000

agg_intensity <- sum(future_covariates$intensity_above_0, na.rm = TRUE)

future_covariates <- future_covariates %>% group_by(Year) %>%
  mutate(agg_intensity_per_year = sum(intensity_above_0, na.rm = TRUE))%>%
  ungroup()

future_covariates$normalised_intensity <- future_covariates$intensity_above_0 / future_covariates$agg_intensity_per_year    
future_covariates$endpoint <- 1
future_covariates$endpoint <- -(thresh_fit_A2$par[1] + thresh_fit_A2$par[2] * future_covariates$ICS_max)/thresh_fit_A2$par[3]

(endpoint_max <- max(future_covariates$endpoint, na.rm = TRUE))

endpoint_by_year <- future_covariates %>% group_by(Year) %>%
  summarise(endpoint_wm = sum(endpoint*normalised_intensity, na.rm = TRUE), agg_intensity_per_year = sum(intensity_above_0, na.rm = TRUE)) 

endpoint_by_year <- endpoint_by_year[-nrow(endpoint_by_year),]

sum(endpoint_by_year$endpoint_wm * endpoint_by_year$agg_intensity_per_year/agg_intensity)

dev.new(height=5, width=10, noRStudioGD = TRUE)
par(mfrow=c(1,2), bg='transparent')
plot(endpoint_by_year$Year, endpoint_by_year$agg_intensity_per_year/agg_intensity, type = "l", 
     xlab = "Year", ylab = "Ratio of integrated intensities", lwd=2, ylim=c(0,0.08))
plot(endpoint_by_year$Year, endpoint_by_year$endpoint_wm, type = "l", 
     xlab = "Year", ylab = "Weighted mean endpoint", lwd=2)


chosen_dates <- c("2025-01-01", "2040-01-01", "2055-01-01")
future_covariates_selected <- future_covariates %>% filter(Date %in% chosen_dates)
# Plot of normalised intensity for 2025, 2040 and 2055


endpoint_wm_by_year <- future_covariates %>% group_by(Year) %>%
  summarise(endpoint_wm_year = mean(endpoint_wm, na.rm = TRUE)) 

intensity_by_year <- future_covariates %>% group_by(Year) %>%
  summarise(agg_intensity_year = mean(intensity_above_0, na.rm = TRUE))
plot(endpoint_wm_by_year$Year,endpoint_wm_by_year$endpoint_wm_year , type = "l", 
     xlab = "Date", ylab = "Weighted mean endpoint", lwd=2)
plot(stress_by_year$Year, stress_by_year$stress_year, type = "l", 
     xlab = "Date", ylab = "Mean stress (ICS_max)", lwd=2)
plot(agg_intensity_by_year$Year, agg_intensity_by_year$agg_intensity_year, type = "l",
     xlab = "Date", ylab = "Mean aggregated intensity", lwd=2)

# Plot of normalised intensity for 2025, 2040 and 2055
library(ggplot2)
library(dplyr)
future_covariates_2025 <- future_covariates %>% filter(Date == "2025-01-01") 
future_covariates_2040 <- future_covariates %>% filter(Date == "2040-01-01")
future_covariates_2055 <- future_covariates %>% filter(Date == "2055-01-01")

change_2025_2040 <- future_covariates_2040 %>%
    left_join(future_covariates_2025, by = c("Easting", "Northing"), suffix = c("_2040", "_2025")) %>%
    mutate(change_normalised_intensity = normalised_intensity_2040 - normalised_intensity_2025,
           change_endpoint = endpoint_2040 - endpoint_2025)
change_2040_2055 <- future_covariates_2055 %>%
    left_join(future_covariates_2040, by = c("Easting", "Northing"), suffix = c("_2055", "_2040")) %>%
    mutate(change_normalised_intensity = normalised_intensity_2055 - normalised_intensity_2040,
           change_endpoint = endpoint_2055 - endpoint_2040)

fill_limits <- range(c(future_covariates_2025$normalised_intensity,
                       future_covariates_2040$normalised_intensity,
                       future_covariates_2055$normalised_intensity), na.rm = TRUE)

# plot1 <- ggplot(future_covariates_2025, aes(x = Easting, y = Northing, fill = normalised_intensity)) +
#     geom_tile() + fixed_plot_aspect(ratio = 1) + theme_classic() +
#     theme(plot.background = element_blank()) + scale_fill_gradient(low = "blue", high = "red", limits = fill_limits) +
#     labs(x = "Easting (m)", y = "Northing (m)", fill = expression(lambda[0]/Lambda[0])) + coord_fixed()
# plot2 <- ggplot(future_covariates_2040, aes(x = Easting, y = Northing, fill = normalised_intensity)) +
#     geom_tile() + fixed_plot_aspect(ratio = 1) + theme_classic() +
#     theme(plot.background = element_blank()) + scale_fill_gradient(low = "blue", high = "red", limits = fill_limits) +
#     labs(x = "Easting (m)", y = "Northing (m)", fill = expression(lambda[0]/Lambda[0])) + coord_fixed()
# plot3 <- ggplot(future_covariates_2055, aes(x = Easting, y = Northing, fill = normalised_intensity)) +
#     geom_tile() + fixed_plot_aspect(ratio = 1) + theme_classic() +
#     theme(plot.background = element_blank()) + scale_fill_gradient(low = "blue", high = "red", limits = fill_limits) +
#     labs(x = "Easting (m)", y = "Northing (m)", fill = expression(lambda[0]/Lambda[0])) + coord_fixed()

# fill_limits <- range(c(future_covariates_2025$endpoint,
#                        future_covariates_2040$endpoint,
#                        future_covariates_2055$endpoint), na.rm = TRUE)
# plot1 <- ggplot(future_covariates_2025, aes(x = Easting, y = Northing, fill = endpoint)) +
#     geom_tile() + fixed_plot_aspect(ratio = 1) + theme_classic() +
#     theme(plot.background = element_blank()) + scale_fill_gradient(low = "blue", high = "red", limits = fill_limits) +
#     labs(x = "Easting (m)", y = "Northing (m)", fill = "endpoint") + coord_fixed()
# plot2 <- ggplot(future_covariates_2040, aes(x = Easting, y = Northing, fill = endpoint)) +
#     geom_tile() + fixed_plot_aspect(ratio = 1) + theme_classic() +
#     theme(plot.background = element_blank()) + scale_fill_gradient(low = "blue", high = "red", limits = fill_limits) +
#     labs(x = "Easting (m)", y = "Northing (m)", fill = "endpoint") + coord_fixed()
# plot3 <- ggplot(future_covariates_2055, aes(x = Easting, y = Northing, fill = endpoint)) +
#     geom_tile() + fixed_plot_aspect(ratio = 1) + theme_classic() +
#     theme(plot.background = element_blank()) + scale_fill_gradient(low = "blue", high = "red", limits = fill_limits) +
#     labs(x = "Easting (m)", y = "Northing (m)", fill = "endpoint") + coord_fixed()

fill_limits <- range(c(change_2025_2040$change_endpoint,
                       change_2040_2055$change_endpoint), na.rm = TRUE)
plot1 <- ggplot(future_covariates_2025, aes(x = Easting, y = Northing, fill = endpoint)) +
  geom_tile() + fixed_plot_aspect(ratio = 1) + theme_classic() +
  theme(plot.background = element_blank()) + scale_fill_gradient(low = "blue", high = "red") +
  labs(x = "Easting (m)", y = "Northing (m)", fill = "Endpoint") + coord_fixed()
plot2 <- ggplot(change_2025_2040, aes(x = Easting, y = Northing, fill = change_endpoint)) +
  geom_tile() + fixed_plot_aspect(ratio = 1) + theme_classic() +
  theme(plot.background = element_blank()) + scale_fill_gradient(low = "blue", high = "red", limits = fill_limits) +
  labs(x = "Easting (m)", y = "Northing (m)", fill = "Change in Endpoint") + coord_fixed()
plot3 <- ggplot(change_2040_2055, aes(x = Easting, y = Northing, fill = change_endpoint)) +
  geom_tile() + fixed_plot_aspect(ratio = 1) + theme_classic() +
  theme(plot.background = element_blank()) + scale_fill_gradient(low = "blue", high = "red", limits = fill_limits) +
  labs(x = "Easting (m)", y = "Northing (m)", fill = "Change in Endpoint") + coord_fixed()

# Combine the plots using patchwork
plots <- list(plot2, plot3)
# Confirm that both are valid ggplot objects
dev.new(height=5, width=10, noRStudioGD = TRUE)
par(mfrow=c(1,1), bg='transparent')
if (all(sapply(plots, inherits, "ggplot"))) {
  combined_plot <- wrap_plots(plots, guides = "collect") & theme(legend.position = "right")
  print(combined_plot)
} else {
  stop("One or more plots are not ggplot objects.")
}

unique(future_covariates_2025$endpoint_wm2)
unique(future_covariates_2040$endpoint_wm2)
unique(future_covariates_2055$endpoint_wm2)
# aggregated intensity for region
aggregated_intensity_v <- function(v, future_covariates, thresh_fit, region = NULL, grid_box_area = 0.2445286) {
  intensity_above_v <-  future_covariates$intensity_above_0 * (1 + thresh_fit$par[3] * (v) / future_covariates$sigma_0)^(-1 / thresh_fit$par[3])
  # Do we need below?
  intensity_above_v[1 + thresh_fit$par[3] * (v) / future_covariates$sigma_0 < 0] <- 0
  # Filter by region
  if (!is.null(region)) {
    intensity_above_v <- intensity_above_v[inpolygon(future_covariates$Easting, future_covariates$Northing, 
                                                      region$Easting, region$Northing),]
  }
  # Calculate aggregated intensity
  agg_intensity_v <- sum(intensity_above_v, na.rm = TRUE) * grid_box_area
  return(agg_intensity_v)
}

#aggregated_intensity_v(0, future_covariates, thresh_fit_A2)

# v-level extreme event estimation
v_level_extreme_event <- function(v, prob, future_covariates, thresh_fit, grid_box_area = 0.2445286) {
  # Calculate aggregated intensity
  intensity_above_v <-  future_covariates$intensity_above_0 * (1 + thresh_fit$par[3] * (v) / future_covariates$sigma_0)^(-1 / thresh_fit$par[3])
  # Do we need below?
  intensity_above_v[1 + thresh_fit$par[3] * (v) / future_covariates$sigma_0 < 0] <- 0
  #Calculate aggregated intensity
  agg_intensity_v <- sum(intensity_above_v, na.rm = TRUE) * grid_box_area
  
  function_to_solve <- (agg_intensity_v + log(prob))^2
  
  return(function_to_solve)
}

# v_level solver function
v_level_solver <- function(prob, future_covariates, thresh_fit, region = NULL, grid_box_area = 0.2445286) {
  # Define the function to solver
  if(!is.null(region)) {
    future_covariates_in_region <- future_covariates[inpolygon(future_covariates$Easting, future_covariates$Northing, 
                                                      region$Easting, region$Northing),] 
  } else {
    future_covariates_in_region <- future_covariates
  }
  # Use uniroot to find the root
  v_solution <- optimize(v_level_extreme_event, interval = c(0, 10), prob = 0.9387, 
                         future_covariates = future_covariates_in_region, thresh_fit = thresh_fit_A2)
  return(v_solution)
}

v_level_solver(prob = 0.9387, future_covariates = future_covariates, 
                thresh_fit = thresh_fit_A2, region = gron_outline)
  
