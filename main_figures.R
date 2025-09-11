# Code to reproduce figures in the main text 

# Load required datasets -------------------------------------------------------
file_paths <- list(
  gron_eq_cat = "data/events/unrounded_after_1995_in_polygon_with_covariates.csv",
  covariates = "data/covariates/covariates_1995-2024.csv",
  covariates_2025 = "data/covariates/covariates_1995-2055.csv",
  geophones_deepest = "data/geophones/Geophones_processed_03-07-2024_deepest_only.csv",
  gron_outline = "data/geophones/Groningen_Field_outline.csv",
  gron_polygon = "data/geophones/polygon_for_groningen_earthquakes.txt", 
  covariates_in_G = "data/covariates/covariates_in_gasfield_1995-2024.csv"
)

# TODO: CM to replace with informative names that match LaTeX source
output_paths <- list(
  fig_1 = "outputs/figures/fig_1.pdf",
  fig_2 = "outputs/figures/fig_2.pdf", 
  fig_3a = "outputs/figures/fig_3a.pdf",
  fig_3b = "outputs/figures/fig_3b.pdf", 
  fig_4 = "outputs/figures/fig_4.pdf",
  data_3 = "Data/covariates/average_ICS_max_1995-2055.rds"
)


chosen_years <- c("2010", "2020") # in Figure 2, plot geophone networks in these years


gron_eq_cat <- read.csv(file_paths$gron_eq_cat, header = TRUE)
covariates <- read.csv(file_paths$covariates, header = TRUE)
covariates_2055 <- read.csv(file_paths$covariates_2025, header = TRUE)
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


# Figure 1 ----------------------------------------------------------------

# changepoint threshold and earthquake locations
change_index <- which(gron_eq_cat$Date == as.Date("2015-12-25"))
pc_threshold <- c(rep(1.15, change_index), rep(0.76, nrow(gron_eq_cat) - change_index))

path <- output_paths$fig_1
pdf(file = path, height = 5, width = 10)
par(mfrow = c(1, 2), bg = 'transparent')
plot(
  x = as.Date(gron_eq_cat$Date),
  y = gron_eq_cat$Magnitude,
  xlab = "Time",
  ylab = expression(Magnitude~(M[L])),
  pch = 19,
  col = "grey",
  cex = 0.7)
lines(
  x = as.Date(gron_eq_cat$Date),
  y = rep(1.45, nrow(gron_eq_cat)),
  col = "red",
  lty = 2,
  lwd = 2)
lines(x = as.Date(gron_eq_cat$Date), y = pc_threshold, col = "blue", lwd = 2)

plot( 
  x = gron_polygon$POINT_X,
  y = gron_polygon$POINT_Y,
  col = "green",
  xlab = "Easting (m)",
  ylab = "Northing (m)",
  type = 'l',
  lty = 2,
  asp = 1,
  lwd = 2)
points(x = gron_eq_cat$Easting, y = gron_eq_cat$Northing, pch = 19, col = "grey", cex = 0.7, asp = 1)
lines(x = gron_outline$Easting, y = gron_outline$Northing, col = "black", asp = 1)
dev.off()

# Figure 2 ----------------------------------------------------------------

# Number of geophones over time
dates <- seq(from = min(as.Date(gron_eq_cat$Date)),
             to = max(as.Date(gron_eq_cat$Date)),
             by = "day")
num_geophones_in_rect <- numeric(length(dates))
num_geophones_in_polygon <- numeric(length(dates))
num_geophones_in_gasfield <- numeric(length(dates))

# subset of geophones belonging to each spatial region

geophones_deepest$in_rect  <- inpolygon(geophones_deepest$Xcoord,
                                        geophones_deepest$Ycoord,
                                        gron_rect$X,
                                        gron_rect$Y)
geophones_deepest$in_poly  <- inpolygon(geophones_deepest$Xcoord,
                                        geophones_deepest$Ycoord,
                                        gron_polygon$POINT_X,
                                        gron_polygon$POINT_Y)
geophones_deepest$in_field <- inpolygon(geophones_deepest$Xcoord,
                                        geophones_deepest$Ycoord,
                                        gron_outline$Easting,
                                        gron_outline$Northing)

geophones_in_rect <- filter(geophones_deepest, in_rect)
geophones_in_polygon <- filter(geophones_deepest, in_poly)
geophones_in_gasfield <- filter(geophones_deepest, in_field)

filter_to_date <- function(df, date){df[df$Start_date <= date & df$End_date >= date,]}

for (i in seq_along(dates)) {
  date <- dates[i]
  num_geophones_in_rect[i] <- filter_to_date(geophones_in_rect, date) |> nrow()
  num_geophones_in_polygon[i] <- filter_to_date(geophones_in_polygon, date) |> nrow()
  num_geophones_in_gasfield[i] <- filter_to_date(geophones_in_gasfield, date) |> nrow() 
  if (i %% 250 == 0) cat(i, " of ", length(dates), "\n")
}

# construct plots
path <- output_paths$fig_2
pdf(file = path, height = 5, width = 15)
par(mfrow = c(1, 3), bg = 'transparent')

plot(
  x = dates,
  y = num_geophones_in_rect,
  xlab = "Time",
  ylab = "Number of geophones",
  type = 'l',
  col = "blue")
lines(x = dates, y = num_geophones_in_polygon, col = "green")
lines(x = dates, y = num_geophones_in_gasfield, col = "black")

# Spatial plots of geophones active in particular years
for (year in chosen_years) {
  is_in_year <- stringr::str_sub(geophones_deepest$Start_date, 1, 4) <= year &
                stringr::str_sub(geophones_deepest$End_date, 1, 4) >= year
  current_geophones <- geophones_deepest[is_in_year,]
  plot(
    x = gron_polygon$POINT_X, 
    y = gron_polygon$POINT_Y,
    col = "green",
    xlab = "Easting (m)",
    ylab = "Northing (m)",
    type = 'l',
    lty = 2,
    lwd = 2,
    asp = 1)
  points(x = current_geophones$Xcoord,
         y = current_geophones$Ycoord,
         pch = 4,
         col = "blue",
         cex = 0.7)
  lines(x = gron_outline$Easting, y = gron_outline$Northing, col = "black")
}
dev.off()

# Figure 3 ----------------------------------------------------------------

# Average Kaiser stress plots for different years -------------------------

# Filter data for the selected year
current_covariates <- filter(covariates, Year == "2020")

# Compute average ICS_max per grid cell
plot_df <- current_covariates %>%
  group_by(Easting, Northing) %>%
  summarise(Average_ICS = mean(ICS_max, na.rm = TRUE), .groups = "drop") %>%
  filter(complete.cases(.))

# Plot

#dev.new(height = 5, width = 5, noRStudioGD = TRUE)
#par(mfrow = c(1,1), bg = 'transparent')
path <- output_paths$fig_3a 
pdf(file = path, height = 5, width = 5)
ggplot(
  data = plot_df,
  aes(x = Easting, y = Northing, fill = Average_ICS)) +
  geom_tile() +
  scale_fill_gradient(low = "red", high = "yellow") +
  geom_point(data = gron_outline,
             aes(x = Easting, y = Northing),
             size = 0.5,
             shape = 1,
             fill = "black") +
  geom_point(data = locations,
             aes(x = Easting, y = Northing),
             size = 2,
             shape = 19,
             fill = "black") +
  coord_fixed() +
  theme_classic() +
  theme(plot.background = element_blank()) +
  labs(x = "Easting (m)", y = "Northing (m)", fill = "Average KS")
dev.off()

# Average max stress over time (including future)
average_ICS_df <- covariates_2055 %>%
  group_by(Date) %>%
  summarise(average_ICS = mean(ICS_max, na.rm = TRUE), .groups = "drop") %>% 
  mutate(Date = as.Date(Date))

if (!file.exists(output_paths$data_3)) {
  saveRDS(average_ICS_df, file = output_paths$data_3)
}

ics_at_locations <- vector(mode = "list", length = nrow(locations))
for (i in 1:nrow(locations)) {
  ics_at_locations[[i]] <- covariates_2055 %>%  
    filter(Easting == locations$Easting[i] & Northing == locations$Northing[i]) %>% 
    pull(ICS_max)
}

path <- output_paths$fig_3b
pdf(file = path, height = 5, width = 5)
par(mfrow = c(1,1), bg = 'transparent')
plot(x = average_ICS_df$Date,
     y = average_ICS_df$average_ICS,
     xlab = "Time",
     ylab = "Average KS (MPa)",
     type = 'l',
     lwd = 2,
     ylim = c(0, 0.4),
     col = "black")
lines(average_ICS_df$Date, ics_at_locations[[1]], lwd = 2, col = "blue")
lines(average_ICS_df$Date, ics_at_locations[[2]], lwd = 2, col = "green",)
lines(average_ICS_df$Date, ics_at_locations[[3]], lwd = 2, col = "red")
dev.off()


# Figure 4 ----------------------------------------------------------------

# Loasbest performing thresholds (based on EQD)
thresh_fit_A2 <- readRDS("outputs/threshold_results/geo_thresh_fit_V2.rds")
thresh_fit_B1 <- readRDS("outputs/threshold_results/geo_thresh_fit_logV1.rds")
thresh_fit_C2 <- readRDS("outputs/threshold_results/geo_thresh_fit_sqrtV2.rds")

# add threshold values (and differences) at each earthquake location
gron_eq_cat <- gron_eq_cat %>% mutate(
  Date = as.Date(Date),
  threshold_A2 = thresh_fit_A2$thresh_par[1] + thresh_fit_A2$thresh_par[2] * V_2,
  threshold_B1 = thresh_fit_B1$thresh_par[1] + thresh_fit_B1$thresh_par[2] * log(V_2),
  threshold_C2 = thresh_fit_C2$thresh_par[1] + thresh_fit_C2$thresh_par[2] * sqrt(V_2),
  threshold_PC = if_else(as.Date(Date) > "2015-12-25", true = 0.76, false = 1.15),
  threshold_con = 1.45, 
  diff_thr_B1 = threshold_A2 - threshold_B1, 
  diff_thr_C2 = threshold_A2 - threshold_C2
)

# Evaluate threshold A2 threshold for entire observation window
covariates$threshold_A2 <- thresh_fit_A2$thresh_par[1] + thresh_fit_A2$thresh_par[2] * covariates$V2

# Add indicator of which spatio-temporal grid cells are within the field outline
# (for efficiency, do first time point and repeat)
cov_slice <- covariates %>% filter(MonthYear == "1995-04")
in_field <- inpolygon(cov_slice$Easting,
                      cov_slice$Northing,
                      gron_outline$Easting,
                      gron_outline$Northing)
n_MonthYear <- length(unique(covariates$MonthYear))
covariates$in_field <- rep(x = in_field, times = n_MonthYear)
rm(cov_slice, n_MonthYear, in_field)


average_A2_threshold_over_field <- covariates %>% 
  filter(in_field) %>%
  group_by(Date) %>%
  summarise(threshold = mean(threshold_A2, na.rm = TRUE), .groups = "drop")

thresh_A2_at_locations <- vector("list", length = nrow(locations))
for (i in 1:3) {
  ei <- locations$Easting[i]
  ni <- locations$Northing[i]
  
  thresh_A2_at_locations[[i]] <- covariates %>% 
    filter(Easting == ei & Northing == ni) %>% 
    pull(threshold_A2)
}
rm(ei, ni)

# Create figures ---
path = output_paths$fig_4
pdf(file = path, height = 5, width = 5)
par(mfrow = c(1,1), bg = 'transparent')

plot(
  x = gron_eq_cat$Date,
  y = gron_eq_cat$Magnitude,
  xlab = "Event time",
  ylab = expression(Magnitude~(M[L])),
  pch = 19,
  col = "grey",
  cex = 0.7)
with(gron_eq_cat, lines(x = Date, y = threshold_A2, lwd = 2))
with(gron_eq_cat, lines(x = Date, y = threshold_con, col = "red", lty = 2, lwd = 2))
with(gron_eq_cat, lines(x = Date, y = pc_threshold, col = "blue", lwd = 2))
with(average_A2_threshold_over_field[-1,], lines(x = as.Date(Date), y = threshold, lwd = 2, col = "green"))

thresh_df <- data.frame(
  date = as.Date(unique(covariates$Date)[-1]),
  loc_1 = thresh_A2_at_locations[[1]][-1],
  loc_2 = thresh_A2_at_locations[[2]][-1],
  loc_3 = thresh_A2_at_locations[[3]][-1])
plot(
  x = gron_eq_cat$Date,
  y = gron_eq_cat$Magnitude,
  xlab = "Event time",
  ylab = expression(Magnitude~(M[L])),
  pch = 19,
  col = "grey",
  cex = 0.7)
lines(x = thresh_df$date, y = thresh_df$loc_1, lwd = 2, col = "purple")
lines(x = thresh_df$date, y = thresh_df$loc_2, lwd = 2, col = "darkgreen")
lines(x = thresh_df$date, y = thresh_df$loc_3, lwd = 2, col = "orange")
lines(x = gron_eq_cat$Date, y = gron_eq_cat$threshold_con, lty = 2, lwd = 2, col = "red")
lines(x = gron_eq_cat$Date, y = gron_eq_cat$threshold_PC, lwd = 2, col = "blue")
rm(thresh_df)

y_min <- min(gron_eq_cat$diff_thr_B1, gron_eq_cat$diff_thr_C2)
y_max <- max(gron_eq_cat$diff_thr_B1, gron_eq_cat$diff_thr_C2)
plot(
  x = gron_eq_cat$Date,
  y = gron_eq_cat$diff_thr_B1,
  xlab = "Event time",
  ylab = "Difference in threshold",
  col = "yellow",
  type = 'l', 
  ylim = c(y_min, y_max))
lines(x = gron_eq_cat$Date, y = gron_eq_cat$diff_thr_C2, col = "red")
abline(h = 0, col = "black", lty = 2)
rm(y_min, y_max)
dev.off()

# Figure 5 ----------------------------------------------------------------

# Spatial threshold plots
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
    labs(fill = "Threshold", x = "Easting (m)", y = "Northing (m)") 
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



# Figure 6 ----------------------------------------------------------------

#QQplots
dev.new(height=5, width=5, noRStudioGD = TRUE)
par(mfrow=c(1,1), bg='transparent')

threshold <- 1.45

excess_data <- gron_eq_cat[gron_eq_cat$Magnitude > threshold,]

fit_obs <- optim(GPD_LL_given_V_ICS, par = c(0.1, 0, 0.1), excess = excess_data$Magnitude - threshold,
                 thresh_par = c(1.45, 0), V = excess_data$V_1, ics = excess_data$ICS_max,
                 control = list(fnscale = -1))


thresh_fit <- list(thresh_par = c(1.45, 0), par = fit_obs$par)
get_qq_plot_geo_ics(gron_eq_cat$Magnitude, thresh_fit, gron_eq_cat$V_1, gron_eq_cat$ICS_max, main="" )
get_qq_plot_geo_ics(gron_eq_cat$Magnitude, thresh_fit_A2, gron_eq_cat$V_2, gron_eq_cat$ICS_max, main="" )

#PPplots (in supp)
dev.new(height=5, width=5, noRStudioGD = TRUE)
par(mfrow=c(1,1), bg='transparent')

get_pp_plot_geo_ics(gron_eq_cat$Magnitude, thresh_fit, gron_eq_cat$V_1, gron_eq_cat$ICS_max, main="", n_boot=1000 )
get_pp_plot_geo_ics(gron_eq_cat$Magnitude, thresh_fit_A2, gron_eq_cat$V_2, gron_eq_cat$ICS_max, main="", n_boot=1000 )



# Poisson process intensity fit -------------------------------------------

gron_eq_cat_exceed_A2 <- gron_eq_cat[gron_eq_cat$Magnitude > gron_eq_cat$threshold_A2,] 
#Checking how many events are removed by dsmaxdt condition
sum(gron_eq_cat_exceed_A2$dsmaxdt == 0)


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

# Figure 7 ----------------------------------------------------------------

# Aggregated intensity over space for all years ------------

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


# Aggregated intensity above 0 for all years ----------
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

# Proportion observed relative to expected (above 0)---------
ggplot(agg_intensity_df, aes(x = Year)) +
  geom_point(aes(y = num_exceedances/(agg_intensity*grid_box_area), color = "Proportion observed"), size = 2) +
  geom_line(aes(y = num_exceedances/(agg_intensity*grid_box_area), color = "Proportion observed"), size = 1) +
  labs(x = "Year", y = "Proportion observed") +
  theme_classic() +
  theme(plot.background = element_blank()) +
  scale_color_discrete(guide="none")




# Figure 8 ----------------------------------------------------------------

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


# Figure 9 ----------------------------------------------------------------

threshold_values_uncertainty_results_Alg2 <- readRDS("uncertainty/threshold_values_uncertainty_results_Alg2.rds")
bootstrap_model_selection_results_Alg3 <- readRDS("uncertainty/bootstrap_model_selection_results_Alg3.rds")
thresh_fit_A2 <- readRDS("threshold_results/geo_thresh_fit_V2.rds")

# Algorithm 2 -----------------------
thresh_par_Alg2 <- do.call(rbind,lapply(threshold_values_uncertainty_results_Alg2, function(x){
  x$thresh_par
}))

threshold_A2_over_G <- covariates_in_G %>% group_by(Date) %>%
  summarise(mean_thresh = mean(thresh_fit_A2$thresh_par[1] + thresh_fit_A2$thresh_par[2] * V2, na.rm = TRUE), .groups = "drop")

dev.new(height=5, width=5, noRStudioGD = TRUE)
par(mfrow=c(1,1), bg='transparent')
plot(as.Date(gron_eq_cat$Date), gron_eq_cat$Magnitude, xlab = "Time", ylab = "Average threshold across gasfield", 
     ylim=c(min(gron_eq_cat$Magnitude),4.5),pch=19, col="grey", cex=0.7)

mean_mat <- matrix(0, nrow = nrow(thresh_par_Alg2), ncol = length(unique(covariates_in_G$Date)))
for(i in 1:nrow(thresh_par_Alg2)) {
  average_boot_thresh <- covariates_in_G %>%
    group_by(Date) %>%
    summarise(mean_thresh = mean(thresh_par_Alg2[i, 1] + thresh_par_Alg2[i, 2] * V2, na.rm = TRUE), .groups = "drop")
  mean_mat[i, ] <- average_boot_thresh$mean_thresh
  lines(as.Date(average_boot_thresh$Date[-1]), average_boot_thresh$mean_thresh[-1], col=rgb(0,0,0.8,0.1))
}
lines(as.Date(threshold_A2_over_G$Date[-1]), threshold_A2_over_G$mean_thresh[-1], col="green", lwd=3)
upper_bound <- apply(mean_mat, 2, quantile, probs = 0.975)
lower_bound <- apply(mean_mat, 2, quantile, probs = 0.025)
lines(as.Date(average_boot_thresh$Date[-1]), upper_bound[-1], col="darkorange", lwd=2)
lines(as.Date(average_boot_thresh$Date[-1]), lower_bound[-1], col="darkorange", lwd=2)

# Algorithm 3 -----------------------
thresh_par_Alg3 <- do.call(rbind,lapply(bootstrap_model_selection_results_Alg3, function(x){
  x$model_results$thresh_par
}))

chosen_models <- do.call(rbind, lapply(bootstrap_model_selection_results_Alg3, function(x){
  x$chosen_form
}))

threshold_A2_over_G <- covariates_in_G %>% group_by(Date) %>%
  summarise(mean_thresh = mean(thresh_fit_A2$thresh_par[1] + thresh_fit_A2$thresh_par[2] * V2, na.rm = TRUE), .groups = "drop")

covariates_distances_in_G <- cbind(covariates_in_G$V1, covariates_in_G$V2, covariates_in_G$V3, covariates_in_G$V4,
                                   log(covariates_in_G$V1), log(covariates_in_G$V2), log(covariates_in_G$V3), log(covariates_in_G$V4),
                                   sqrt(covariates_in_G$V1), sqrt(covariates_in_G$V2), sqrt(covariates_in_G$V3), sqrt(covariates_in_G$V4))

dev.new(height=5, width=5, noRStudioGD = TRUE)
par(mfrow=c(1,1), bg='transparent')
plot(as.Date(gron_eq_cat$Date), gron_eq_cat$Magnitude, xlab = "Time", ylab = "Average threshold across gasfield", 
     ylim=c(min(gron_eq_cat$Magnitude),4.5),pch=19, col="grey", cex=0.7)

mean_mat <- matrix(0, nrow = nrow(thresh_par_Alg3), ncol = length(unique(covariates_in_G$Date)))
for(i in 1:nrow(thresh_par_Alg3)) {
  covariates_in_G$current_dist <- covariates_distances_in_G[,chosen_models[i]]
  average_boot_thresh <- covariates_in_G %>%
    group_by(Date) %>%
    summarise(mean_thresh = mean(thresh_par_Alg3[i, 1] + thresh_par_Alg3[i, 2] * current_dist, na.rm = TRUE), .groups = "drop")
  mean_mat[i, ] <- average_boot_thresh$mean_thresh
  lines(as.Date(average_boot_thresh$Date[-1]), average_boot_thresh$mean_thresh[-1], col=rgb(0,0,0.8,0.1))
}
lines(as.Date(threshold_A2_over_G$Date[-1]), threshold_A2_over_G$mean_thresh[-1], col="green", lwd=3)
upper_bound <- apply(mean_mat, 2, quantile, probs = 0.975)
lower_bound <- apply(mean_mat, 2, quantile, probs = 0.025)
lines(as.Date(average_boot_thresh$Date[-1]), upper_bound[-1], col="darkorange", lwd=2)
lines(as.Date(average_boot_thresh$Date[-1]), lower_bound[-1], col="darkorange", lwd=2)




