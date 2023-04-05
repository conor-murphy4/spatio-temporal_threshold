# Groningnen catalogue plots:
# 0) Filter full catalogue to date range [1995-01-01 to 2020-01-01)
# 1) Plot on natural time scale with loess mean
# 2) Plot on index time scale with loess mean

## SOURCE ----------------------------------------------------------------------
library(dplyr)

## LOAD DATA -------------------------------------------------------------------
cat_full <- readr::read_csv(file = "C:/Users/murphyc4/OneDrive - Lancaster University/STOR-i/PhD/Projects/Spatio-temporal/Data/Events/2022-04-12_15-09-25_cat.csv")
cat_full <- cat_full %>%
  mutate(index = order(julian2)) %>%
  filter(date >= as.Date("1995-01-01")) %>%
  filter(date < as.Date("2022-01-01")) %>%
  select(mag,index, date)

## PLOT CATALOGUES -------------------------------------------------------------

pdf(file = 'C:/Users/murphyc4/OneDrive - Lancaster University/STOR-i/PhD/Projects/Spatio-temporal/Data/Events/groningen_catalgoue.pdf', width = 8, height = 5)
## PLOT 1 -- Natural timescale
dev.new(width=9.17, height=5.53,noRStudioGD = TRUE)
par(mfrow=c(1,1),bg='transparent')
ggplot2::ggplot(cat_full, ggplot2::aes(x = date, y = mag)) +
  ggplot2:: geom_point(col = 'grey60', size = 2.5) +
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

## SAME PLOTS BUT PNG ----------------------------------------------------------


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

