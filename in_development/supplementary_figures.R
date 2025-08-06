

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
