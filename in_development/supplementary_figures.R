# Code to reproduce figures in the main text 

# Load required datasets -------------------------------------------------------
file_paths <- list(
  gron_eq_cat = "data/events/unrounded_after_1995_in_polygon_with_covariates.csv",
  covariates = "data/covariates/covariates_1995-2024.csv",
  covariates_2025 = "data/covariates/covariates_1995-2055.csv",
  geophones_deepest = "data/geophones/Geophones_processed_03-07-2024_deepest_only.csv",
  gron_outline = "data/geophones/Groningen_Field_outline.csv",
  gron_polygon = "data/geophones/polygon_for_groningen_earthquakes.txt", 
  covariates_in_G = "data/covariates/covariates_in_gasfield_1995-2024.csv",
  alg3_results = "in_development/uncertainty/bootstrap_model_selection_results_Alg3.rds"
)

# TODO: Remove unnecessary file paths and change file paths when moved out of in_development
# TODO: Set up file same as main_figures.R and include relevant output paths
output_paths <- list(
  fig_1 = "outputs/figures/fig_1_eq_cat.pdf",
  fig_2 = "outputs/figures/fig_2_geophone_network.pdf", 
  fig_3a = "outputs/figures/fig_3a_average_kaiser_stress_2020.pdf",
  fig_3b = "outputs/figures/fig_3b_temporal_kaiser_stress.pdf", 
  data_3 = "Data/covariates/average_ICS_max_1995-2055.rds"
)

alg3_results <- readRDS(file_paths$alg3_results)

# Table S1 -----------------------------------------------------------------


# Proportion of models chosen
chosen_models <- numeric(200)
for(i in 1:200){
  chosen_models[i] <- alg3_results[[i]]$chosen_form
}
table(chosen_models)
prop.table(table(chosen_models))

# Future inference (supp) --------------------------------------------------------
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
