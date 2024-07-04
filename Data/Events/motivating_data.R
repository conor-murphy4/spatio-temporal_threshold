# Groningnen catalogue plots:
# 0) Filter full catalogue to date range [1995-01-01 to 2020-01-01)
# 1) Plot on natural time scale with loess mean
# 2) Plot on index time scale with loess mean

## SOURCE ----------------------------------------------------------------------
library(dplyr)
source('helper_functions.R')

## LOAD DATA -------------------------------------------------------------------
cat_full <- readr::read_csv(file = "Data/Events/2022-04-12_15-09-25_cat.csv")
cat_full <- cat_full %>%
  mutate(index = order(julian2)) %>%
  filter(date >= as.Date("1995-01-01")) %>%
  filter(date < as.Date("2022-01-01")) %>%
  select(mag,index, date)

## PLOT CATALOGUES -------------------------------------------------------------

pdf(file = 'C:/Users/murphyc4/OneDrive - Lancaster University/STOR-i/PhD/Projects/spatio-temporal_threshold/Data/Events/groningen_catalgoue_thresholds.pdf', width = 8, height = 5)
## PLOT 1 -- Natural timescale
dev.new(width=9.17, height=5.53,noRStudioGD = TRUE)
par(mfrow=c(1,1),bg='transparent')
ggplot2::ggplot(cat_full, ggplot2::aes(x = date, y = mag)) +
  ggplot2:: geom_point(col = 'grey60', size = 2.5) +
  ggplot2::geom_hline(yintercept=0.5, color='red', size=1.2) +
  ggplot2::geom_hline(yintercept=1.5, color='blue', size=1.2) +
  ggplot2::geom_smooth(method = 'loess', col = 'black',) +
  ggplot2::ylab('magnitude') +
  ggplot2::xlab('event time') +
  ggplot2::theme_minimal(base_size = 22)

## PLOT 2 -- Index timescale
ggplot2::ggplot(cat_full, ggplot2::aes(x = index, y = mag)) +
  ggplot2:: geom_point(col = 'grey60', size = 2.5) +
  ggplot2::geom_smooth(method = 'loess', col = 'black',) +
  ggplot2::ylab('magnitude') +
  ggplot2::xlab('event index') +
  ggplot2::theme_minimal(base_size = 22)
dev.off()

exceed <- cat_full[cat_full$mag>1.4,]
exceed$gpd_density <- gpd_density
gpd_fit <- optim(GPD_LL, par=c(mean(excess),0.1), z=excess, control=list(fnscale=-1))
gpd_density <- dgpd(excess, shape = gpd_fit$par[2], scale = gpd_fit$par[1])

ex_plot <- data.frame(mag=cat_full$mag, y=c(rep(0, length(cat_full$mag[cat_full$mag<=1.4])),gpd_density))
ggplot2::ggplot(cat_full, ggplot2::aes(x=mag)) + 
  ggplot2::geom_histogram() +
  ggplot2::xlab('magnitude') +
  ggplot2::geom_vline(xintercept = 1.4, color="red", lwd=1) +
  ggplot2::geom_point(data = ex_plot, ggplot2::aes(x=mag, y=y))
  

ggplot2::ggplot(ggplot2::aes(x=cat_full$mag[cat_full$mag>1.4], y=gpd_density))
## SAME PLOTS BUT PNG ----------------------------------------------------------

hist(cat_full$mag, breaks=30, main='', xlab="magnitude")

## PLOT 1 -- Natural timescale
ggplot2::ggplot(cat_full, ggplot2::aes(x = date, y = mag)) +
  ggplot2:: geom_point(col = 'grey60', size = 2.5) +
  ggplot2::geom_smooth(method = 'loess', col = 'black',) +
  ggplot2::ylab('magnitude') +
  ggplot2::xlab('event time') +
  ggplot2::theme_minimal(base_size = 22)
ggplot2::ggsave("./output/plots/groningen_catalogue_natural.jpg", width = 8, height = 5)

## PLOT 2 -- Index timescale
ggplot2::ggplot(cat_full, ggplot2::aes(x = index, y = mag)) +
  ggplot2:: geom_point(col = 'grey60', size = 2.5) +
  ggplot2::geom_smooth(method = 'loess', col = 'black',) +
  ggplot2::ylab('magnitude') +
  ggplot2::xlab('event index') +
  ggplot2::theme_minimal(base_size = 22)
ggplot2::ggsave("./output/plots/groningen_catalogue_index.jpg", width = 8, height = 5)

