# Code to reproduce figures in the main text 

# Load required datasets -------------------------------------------------------
# TODO: Should we also add file paths for threshold results here?
file_paths <- list(
  gron_eq_cat = "data/events/unrounded_after_1995_in_polygon_with_covariates.csv",
  covariates = "C:/Users/murphyc4/OneDrive/OneDrive - Lancaster/STOR-i/PhD/Projects/Induced-seismicity/Other code files/Messy versions/Data/covariates/covariates_1995-2024.csv",
  covariates_full_period = "data/covariates/covariates_1995-2055.csv",
  future_covariates = "C:/Users/murphyc4/OneDrive/OneDrive - Lancaster/STOR-i/PhD/Projects/Induced-seismicity/Other code files/Messy versions/Data/covariates/covariates_2024-2055.csv",
  geophones_deepest = "data/geophones/Geophones_processed_03-07-2024_deepest_only.csv",
  gron_outline = "data/geophones/Groningen_Field_outline.csv",
  gron_polygon = "data/geophones/polygon_for_groningen_earthquakes.txt", 
  covariates_in_G = "data/covariates/covariates_in_gasfield_1995-2024.csv",
  alg3_results = "in_development/uncertainty/bootstrap_model_selection_results_Alg3.rds"
)

# TODO: Remove unnecessary file paths?
# TODO: Change file paths when moved out of in_development
# TODO: Set up file same as main_figures.R and include relevant output paths
# TODO: Edit file path for covariates once set up properly in repo

output_paths <- list(
  fig_1 = "outputs/figures/main/fig_1_eq_cat.pdf",
  fig_2 = "outputs/figures/main/fig_2_geophone_network.pdf", 
  fig_3a = "outputs/figures/main/fig_3a_average_kaiser_stress_2020.pdf",
  fig_3b = "outputs/figures/main/fig_3b_temporal_kaiser_stress.pdf", 
  data_3 = "Data/covariates/average_ICS_max_1995-2055.rds",
  fig_S2 = "outputs/figures/supp/fig_S2_spatial_bootstrap_SE.pdf",
  fig_S3 = "outputs/figures/supp/fig_S3_sigma_0_variation.pdf",
  fig_S4 = "outputs/figures/supp/fig_S4_ppplots.pdf",
  fig_S5 = "outputs/figures/supp/fig_S6_future_eq_properties.pdf",
  fig_S6a = "outputs/figures/supp/fig_S6a_future_spatial_endpoint_2025.pdf",
  fig_S6b = "outputs/figures/supp/fig_S6b_future_spatial_changes_endpoint.pdf"
)

chosen_times <- c("1995-04-01", "2005-01-01", "2015-01-01", "2024-01-01")

alg3_results <- readRDS(file_paths$alg3_results)


gron_eq_cat <- read.csv(file_paths$gron_eq_cat, header = TRUE)
covariates <- read.csv(file_paths$covariates, header = TRUE)
covariates_2055 <- read.csv(file_paths$covariates_2025, header = TRUE)
future_covariates <- read.csv(file_paths$future_covariates, header = TRUE)
geophones_deepest <- read.csv(file_paths$geophones_deepest, header = TRUE, row.names = 1)
gron_outline <- read.csv(file_paths$gron_outline, header = TRUE)
gron_polygon <- read.table(file_paths$gron_polygon, header = TRUE)

covariates_in_G <- read.csv(file_paths$covariates_in_G, header = TRUE)

gron_rect <- data.frame(X = c(210000, 275000, 275000, 210000, 210000),
                        Y = c(560000, 560000, 625000, 625000, 560000))
locations <- data.frame(Easting  = c(250000, 250000, 250000),
                        Northing = c(575000, 590000, 605000))


# load required libraries and functions ----------------------------------------
library(ggplot2)
library(ggspatial)
library(pracma)
library(cowplot)
library(dplyr)
library(purrr)
library(patchwork)

source("src/helper_functions.R")
source("src/intensity_estimation.R")

# Table S1 -----------------------------------------------------------------


# Proportion of models chosen
chosen_models <- numeric(200)
for(i in 1:200){
  chosen_models[i] <- alg3_results[[i]]$chosen_form
}
table(chosen_models)
prop.table(table(chosen_models))

# Figure S.2 --------------------------------------------------------------

# Spatial bootstrap SE plots from Alg 2 and 3

threshold_values_uncertainty_results_Alg2 <- readRDS("in_development/uncertainty/threshold_values_uncertainty_results_Alg2.rds")
bootstrap_model_selection_results_Alg3 <- readRDS("in_development/uncertainty/bootstrap_model_selection_results_Alg3.rds")
thresh_fit_A2 <- readRDS("outputs/threshold_results/geo_thresh_fit_V2.rds")

# Ensure proper date conversion
covariates$Date <- as.Date(covariates$Date)
geophones_deepest$Start_date <- as.Date(geophones_deepest$Start_date)
geophones_deepest$End_date <- as.Date(geophones_deepest$End_date)
chosen_dates <- as.Date(c("2010-01-01", "2020-01-01"))  

covariates_for_dates <- covariates %>% filter(Date %in% chosen_dates)

# Algorithm 2 -----------------------
thresh_par_Alg2 <- do.call(rbind,lapply(threshold_values_uncertainty_results_Alg2, function(x){
  x$thresh_par
}))

# Calculate boostrapped SE for each spatial location and date in covariates_for_dates

covariates_for_dates <- covariates_for_dates %>%
  rowwise() %>%
  mutate(boot_SE_Alg2 = sd(thresh_par_Alg2[, 1] + thresh_par_Alg2[, 2] * V2))

# Algorithm 3 -----------------------
thresh_par_Alg3 <- do.call(rbind,lapply(bootstrap_model_selection_results_Alg3, function(x){
  x$model_results$thresh_par
}))

chosen_models <- do.call(rbind, lapply(bootstrap_model_selection_results_Alg3, function(x){
  x$chosen_form
}))

covariates_distances <- cbind(covariates_for_dates$V1, covariates_for_dates$V2, covariates_for_dates$V3, covariates_for_dates$V4,
                              log(covariates_for_dates$V1), log(covariates_for_dates$V2), log(covariates_for_dates$V3), log(covariates_for_dates$V4),
                              sqrt(covariates_for_dates$V1), sqrt(covariates_for_dates$V2), sqrt(covariates_for_dates$V3), sqrt(covariates_for_dates$V4))


# Calculate boostrapped SE for each spatial location and date in covariates_for_dates
dist_boot <- t(
  covariates_distances[, chosen_models]
)

threshold_boot <- thresh_par_Alg3[, 1] +
  thresh_par_Alg3[, 2] * dist_boot

covariates_for_dates$boot_SE_Alg3 <- apply(threshold_boot, 2, sd)

# Ensure fill scale range is consistent
fill_limits <- range(covariates_for_dates$boot_SE_Alg2, 
                     covariates_for_dates$boot_SE_Alg3, na.rm = TRUE)

# Function to create plot for a given date
plot_SE_for_date_Alg2 <- function(date) {
  current_covariates <- covariates_for_dates %>% filter(Date == date)
  current_geophones <- geophones_deepest %>%
    filter(Start_date <= date, End_date >= date)
  current_geo_in_polygon <- current_geophones %>%
    filter(inpolygon(Xcoord, Ycoord, gron_polygon$POINT_X, gron_polygon$POINT_Y))
  
  ggplot(current_covariates, aes(x = Easting, y = Northing, fill = boot_SE_Alg2)) +
    geom_tile() + 
    fixed_plot_aspect(ratio = 1) +
    scale_fill_gradient(low = "blue", high = "red", limits = fill_limits) +
    coord_fixed() +
    geom_point(data = current_geo_in_polygon, 
               aes(x = Xcoord, y = Ycoord),
               size = 1,
               shape = 19,
               fill = "black") +
    labs(fill = "Std Error", x = "Easting (km)", y = "Northing (km)") +
    scale_x_continuous(labels = function(x) x / 1000) +
    scale_y_continuous(labels = function(x) x / 1000) +
    theme_classic()
}

# Function to create plot for a given date
plot_SE_for_date_Alg3 <- function(date) {
  current_covariates <- covariates_for_dates %>% filter(Date == date)
  current_geophones <- geophones_deepest %>%
    filter(Start_date <= date, End_date >= date)
  current_geo_in_polygon <- current_geophones %>%
    filter(inpolygon(Xcoord, Ycoord, gron_polygon$POINT_X, gron_polygon$POINT_Y))
  
  ggplot(current_covariates, aes(x = Easting, y = Northing, fill = boot_SE_Alg3)) +
    geom_tile() + 
    fixed_plot_aspect(ratio = 1) +
    scale_fill_gradient(low = "blue", high = "red", limits = fill_limits) +
    coord_fixed() +
    geom_point(data = current_geo_in_polygon, 
               aes(x = Xcoord, y = Ycoord),
               size = 1,
               shape = 19,
               fill = "black") +
    labs(fill = "Std Error", x = "Easting (km)", y = "Northing (km)") +
    scale_x_continuous(labels = function(x) x / 1000) +
    scale_y_continuous(labels = function(x) x / 1000) +
    theme_classic()
}


plots_Alg2 <- lapply(chosen_dates, plot_SE_for_date_Alg2)

plots_Alg3 <- lapply(chosen_dates, plot_SE_for_date_Alg3)

# Combine both lists into one
plots <- c(plots_Alg2, plots_Alg3)


path = output_paths$fig_S2
pdf(file = path, height = 10, width = 10)
par(mfrow = c(1,1), bg = 'transparent')

# Confirm that both are valid ggplot objects
if (all(sapply(plots, inherits, "ggplot"))) {
  combined_plot <- wrap_plots(plots, guides = "collect") & theme(legend.position = "right")
  print(combined_plot)
} else {
  stop("One or more plots are not ggplot objects.")
}

dev.off()


# Figure S.3 --------------------------------------------------------------

# Plots of sigma_0 for different time points ------------------------------

# Results for best-performing threshold
thresh_fit_A2 <- readRDS("outputs/threshold_results/geo_thresh_fit_V2.rds")

sig_pars <- thresh_fit_A2$par[1:2]

covariates_for_chosen_times <- filter(covariates, Date %in% chosen_times)

sig_0_values <- covariates_for_chosen_times %>%
  mutate(sig_0 = sig_pars[1] + sig_pars[2] * ICS_max) %>%
  pull(sig_0)

fill_limits <- range(sig_0_values, na.rm = TRUE)

# Function to create plot for chosen date
plot_sigma0_for_date <- function(date) {
  current_covariates <- filter(covariates, Date == date)
  
  plot_df <- current_covariates %>%
    mutate(sig_0 = sig_pars[1] + sig_pars[2] * ICS_max)
  
  ggplot(
    data = plot_df,
    aes(x = Easting, y = Northing, fill = sig_0)) +
    geom_tile() +
    scale_fill_gradient(low = "red", high = "yellow", limits=fill_limits) +
    geom_point(data = gron_outline,
               aes(x = Easting, y = Northing),
               size = 0.5,
               shape = 1,
               fill = "black") +
    coord_fixed() +
    theme_classic() +
    theme(plot.background = element_blank()) +
    scale_x_continuous(
      labels = function(x) x / 1000
    ) +
    scale_y_continuous(
      labels = function(y) y / 1000
    ) +
    labs(x = "Easting (km)", y = "Northing (km)", fill = expression(sigma[0]))
  
}

plots <- lapply(chosen_times, plot_sigma0_for_date)

path <- output_paths$fig_S3 
pdf(file = path, height = 10, width = 10)
par(mfrow = c(2,2), bg = 'transparent')

# Confirm that both are valid ggplot objects
if (all(sapply(plots, inherits, "ggplot"))) {
  combined_plot <- wrap_plots(plots, guides = "collect") & theme(legend.position = "right")
  print(combined_plot)
} else {
  stop("One or more plots are not ggplot objects.")
}

dev.off()

# Figure S.4 --------------------------------------------------------------

# PPplots 

threshold <- 1.45
excess_data <- filter(gron_eq_cat, Magnitude > threshold)
excesses <- excess_data$Magnitude - threshold

fit_obs_without_KS <- optim(
  fn = GPD_LL, 
  par = c(mean(excesses), 0.1), 
  z = excesses, 
  control = list(fnscale = -1))

fit_obs_with_KS <- optim(
  fn = GPD_LL_given_V_ICS,
  par = c(0.1, 0, 0.1),
  excess = excesses,
  thresh_par = c(1.45, 0),
  V = excess_data$V_1,
  ics = excess_data$ICS_max,
  control = list(fnscale = -1))

thresh_fit_with_KS <- list(thresh_par = c(1.45, 0), par = fit_obs_with_KS$par)

path <- output_paths$fig_S4
pdf(file = path, height=5, width=15)  
par(mfrow=c(1,3), bg='transparent')

get_pp_plot_const(gron_eq_cat$Magnitude, threshold, main="" )
get_pp_plot_geo_ics(gron_eq_cat$Magnitude, thresh_fit_with_KS, gron_eq_cat$V_1, gron_eq_cat$ICS_max, main="" )
get_pp_plot_geo_ics(gron_eq_cat$Magnitude, thresh_fit_A2, gron_eq_cat$V_2, gron_eq_cat$ICS_max, main="" )

dev.off()


# Figure S.5 --------------------------------------------------------

# Future endpoint plot over time
(endpoint_max <- max(future_covariates$endpoint, na.rm = TRUE))

endpoint_by_year <- future_covariates %>% group_by(Year) %>%
  summarise(endpoint_wm = sum(endpoint*normalised_intensity, na.rm = TRUE), agg_intensity_per_year = sum(intensity_above_0, na.rm = TRUE)) 

endpoint_by_year <- endpoint_by_year[-nrow(endpoint_by_year),]

sum(endpoint_by_year$endpoint_wm * endpoint_by_year$agg_intensity_per_year/agg_intensity)

path <- output_paths$fig_S5
pdf(file = path, height=5, width=10)  
par(mfrow=c(1,2), bg='transparent')

plot(endpoint_by_year$Year, endpoint_by_year$agg_intensity_per_year/agg_intensity, type = "l", 
     xlab = "Year", ylab = "Ratio of integrated intensities", lwd=2, ylim=c(0,0.08))
plot(endpoint_by_year$Year, endpoint_by_year$endpoint_wm, type = "l", 
     xlab = "Year", ylab = "Weighted mean endpoint", lwd=2)


# Spatial Endpoint plots

future_covariates_2025 <- future_covariates %>% filter(Date == "2025-01-01") 
future_covariates_2040 <- future_covariates %>% filter(Date == "2040-01-01")
future_covariates_2055 <- future_covariates %>% filter(Date == "2055-01-01")

# Endpoint plot for 2025
path <- output_paths$fig_S6a
pdf(file = path, height=5, width=5)  
par(mfrow=c(1,1), bg='transparent')

ggplot(future_covariates_2025, aes(x = Easting, y = Northing, fill = endpoint)) +
  geom_tile() + fixed_plot_aspect(ratio = 1) + theme_classic() +
  theme(plot.background = element_blank()) + scale_fill_gradient(low = "blue", high = "red") +
  labs(x = "Easting (km)", y = "Northing (km)", fill = "Endpoint") + coord_fixed() +
  scale_x_continuous(labels = function(x) x / 1000) +
  scale_y_continuous(labels = function(y) y / 1000)


dev.off()

# Changes in endpoint between 2025-2040 and 2040-2055

change_2025_2040 <- future_covariates_2040 %>%
  left_join(future_covariates_2025, by = c("Easting", "Northing"), suffix = c("_2040", "_2025")) %>%
  mutate(change_normalised_intensity = normalised_intensity_2040 - normalised_intensity_2025,
         change_endpoint = endpoint_2040 - endpoint_2025)
change_2040_2055 <- future_covariates_2055 %>%
  left_join(future_covariates_2040, by = c("Easting", "Northing"), suffix = c("_2055", "_2040")) %>%
  mutate(change_normalised_intensity = normalised_intensity_2055 - normalised_intensity_2040,
         change_endpoint = endpoint_2055 - endpoint_2040)

fill_limits <- range(c(change_2025_2040$change_endpoint,
                       change_2040_2055$change_endpoint), na.rm = TRUE)

plot2 <- ggplot(change_2025_2040, aes(x = Easting, y = Northing, fill = change_endpoint)) +
  geom_tile() + fixed_plot_aspect(ratio = 1) + theme_classic() +
  theme(plot.background = element_blank()) + scale_fill_gradient(low = "blue", high = "red", limits = fill_limits) +
  labs(x = "Easting (km)", y = "Northing (km)", fill = "Change in Endpoint") + coord_fixed() +
  scale_x_continuous(labels = function(x) x / 1000) +
  scale_y_continuous(labels = function(y) y / 1000)

plot3 <- ggplot(change_2040_2055, aes(x = Easting, y = Northing, fill = change_endpoint)) +
  geom_tile() + fixed_plot_aspect(ratio = 1) + theme_classic() +
  theme(plot.background = element_blank()) + scale_fill_gradient(low = "blue", high = "red", limits = fill_limits) +
  labs(x = "Easting (km)", y = "Northing (km)", fill = "Change in Endpoint") + coord_fixed() +
  scale_x_continuous(labels = function(x) x / 1000) +
  scale_y_continuous(labels = function(y) y / 1000)

# Combine the plots using patchwork
plots <- list(plot2, plot3)
# Confirm that both are valid ggplot objects

path <- output_paths$fig_S6b
pdf(file = path, height=5, width=10)  
par(mfrow=c(1,2), bg='transparent')


if (all(sapply(plots, inherits, "ggplot"))) {
  combined_plot <- wrap_plots(plots, guides = "collect") & theme(legend.position = "right")
  print(combined_plot)
} else {
  stop("One or more plots are not ggplot objects.")
}
dev.off()
