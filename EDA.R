library(pracma)

source("src/distance_to_third_nearest_geo.R")

gron_outline <- read.csv('Data/Geophones/Groningen_Field_outline.csv', header=T)
gron_polygon <- read.table('Data/Geophones/polygon_for_groningen_earthquakes.txt', header=T)
geophones <- read.csv("Data/Geophones/Geophones_processed_03-07-2024_without_duplicates.csv", header=T, row.names = 1)
geophones_deepest <- read.csv("Data/Geophones/Geophones_processed_03-07-2024_deepest_only.csv", header=T, row.names = 1)
gron_eq_cat_all_netherlands <- read.csv("Data/Events/unrounded_after_geophone_start.csv", header=T, row.names = 1)
gron_eq_cat <- read.csv("Data/Events/unrounded_after_1995_in_polygon.csv", header=T)
gron_eq_cat_all <- read.csv("Data/Events/unrounded_after_geophone_start_in_polygon_with_V_3d.csv", header=T)
gron_eq_cat_old <- read.csv("Data/Events/2022-04-12_15-09-25_cat.csv", header=T)

# included_geos <- what_geos(gron_eq_cat, geophones)

#Restricting EQs to within the groningen polygon
# gron_eq_cat_in_polygon <- gron_eq_cat[inpolygon(gron_eq_cat$Easting, gron_eq_cat$Northing, gron_polygon$POINT_X, gron_polygon$POINT_Y),]
# write.csv(gron_eq_cat_in_polygon, "Data/Events/unrounded_after_geophone_start_in_polygon.csv")

#Restricting EQs past 1995
# gron_eq_cat_past_1995 <- gron_eq_cat[as.Date(gron_eq_cat$Date) > as.Date("1995-01-01"),]
# write.csv(gron_eq_cat_past_1995, "Data/Events/unrounded_after_1995_in_polygon.csv")

#Comparing old and new datasets

# Plotting both catalogues
dev.new(height=5, width=9, noRStudioGD = TRUE)
par(mfrow=c(1,2), bg='transparent')

#Rounded to 2 d.p
plot(as.Date(gron_eq_cat$Date), gron_eq_cat$Magnitude, main = "Magnitudes over time", xlab="Event time", ylab = "Magnitude")
plot(gron_eq_cat$Magnitude, main = "Magnitudes over index", xlab = "Index", ylab = "Magnitude")

#Rounded to 1 d.p.
plot(as.Date(gron_eq_cat_old$date),gron_eq_cat_old$mag, main = "Magnitudes over time", xlab="Event time", ylab = "Magnitude")
plot(rev(gron_eq_cat_old$mag), main = "Magnitudes over index", xlab = "Index", ylab = "Magnitude")

#Plotting event locations of both datasets
dev.new(width=15, height=10,noRStudioGD = TRUE)
par(mfrow=c(1,1), bg='transparent')
plot(gron_eq_cat$Easting, gron_eq_cat$Northing, pch=19, col="blue",asp=1, xlab="Easting (m)", ylab = "Northing (m)", ylim=c(560000, 630000), xlim=c(200000, 300000) )
points(gron_outline$X, gron_outline$Y, pch=19, cex=0.2, col="red")
points(geophones$Xcoord, geophones$Ycoord, pch=4, col="green")
lines(gron_polygon$POINT_X, gron_polygon$POINT_Y, col="grey", lty="dashed", lw=2)
legend("topright", legend=c("Event locations" , "Geophones", "Field outline", "Polygon"), col=c("blue", "green" , "red", "grey"), pch=c(19,4,NA,NA), lty = c(NA,NA,1,2) )

#Superimposing magnitudes on time for both datasets
plot(as.Date(gron_eq_cat$Date),gron_eq_cat$Magnitude, main = "Magnitudes over time", xlab="Event time", ylab = "Magnitude", pch=19, col="blue")
points(as.Date(gron_eq_cat_old$date), gron_eq_cat_old$mag, pch=19, col="red")

#find largest length of gron_outline
max(gron_outline$X) - min(gron_outline$X)

#Finding distance to third nearest geophone for observed events in gron_eq_cat

third_nearest_dist <- distance_to_third_nearest(gron_eq_cat, geophones_deepest)
gron_eq_cat$V <- third_nearest_dist

third_nearest_dist_3d <- distance_to_third_nearest_3d(gron_eq_cat, geophones_deepest)
third_nearest_dist_3d_all <- distance_to_third_nearest_3d(gron_eq_cat, geophones) 
gron_eq_cat$V_3d <- third_nearest_dist_3d
write.csv(gron_eq_cat, "Data/Events/unrounded_after_1995_in_polygon.csv",row.names = FALSE)

#For each year in gron_eq_cat, find the minimum V
gron_eq_cat$Year <- as.numeric(format(as.Date(gron_eq_cat$Date), "%Y"))

min_V_per_year <- aggregate(gron_eq_cat$V, by=list(gron_eq_cat$Year), FUN=min)

#Calculate minimum V per year using all events which occurred in and before that year
min_V_per_year$x <- sapply(min_V_per_year$Group.1, function(i) min(gron_eq_cat$V[gron_eq_cat$Year <= i]))

plot(min_V_per_year$Group.1, min_V_per_year$x, xlab="Year", ylab="Minimum V", main="Minimum V per year", pch=19)
points(min_V_per_year$Group.1, min_V_per_year$x, col="red",pch=19)

#Exploratory plots of relationship between magnitudes and V
dev.new(width=30, height=10,noRStudioGD = TRUE)
par(mfrow=c(1,3),bg='transparent')
plot(gron_eq_cat$V, gron_eq_cat$Magnitude, xlab = "V_2d", ylab = "Magnitude")
plot(log(gron_eq_cat$V), gron_eq_cat$Magnitude, xlab = "log(V_2d)", ylab = "Magnitude")
plot(sqrt(gron_eq_cat$V), gron_eq_cat$Magnitude, xlab = "sqrt(V_2d)", ylab = "Magnitude")

#Relationship between V_3d and Magnitude
plot(gron_eq_cat$V_3d, gron_eq_cat$Magnitude, xlab = "V_3d", ylab = "Magnitude")
points(gron_eq_cat$V_3d, geo_chosen_threshold_3d, col="red", pch=19)
plot(log(gron_eq_cat$V_3d), gron_eq_cat$Magnitude, xlab = "log(V_3d)", ylab = "Magnitude")
points(log(gron_eq_cat$V_3d), log_geo_chosen_threshold_3d, col="red", pch=19)
plot(sqrt(gron_eq_cat$V_3d), gron_eq_cat$Magnitude, xlab = "sqrt(V_3d)", ylab = "Magnitude")
points(sqrt(gron_eq_cat$V_3d), sqrt_geo_chosen_threshold_3d, col="red", pch=19)

#Splitting gas field in two regions
#Plotting
# dev.new(width=10, height=10,noRStudioGD = TRUE)
# par(mfrow=c(1,1), bg='transparent')
# plot(gron_outline$X, gron_outline$Y, pch=19, col="blue",asp=1, xlab="Easting (m)", ylab = "Northing (m)")
# points(gron_eq_cat$Easting, gron_eq_cat$Northing )
# y_div <- 735000 - 0.6*seq(220000, 280000, by=1)
# points(seq(220000, 280000, by=1), y_div)
# 
# # Relationship with V and Magnitude in the two regions
# gron_eq_cat_upper <- gron_eq_cat[gron_eq_cat$Northing > (735000 - 0.6*gron_eq_cat$Easting),]
# dev.new(width=30, height=10,noRStudioGD = TRUE)
# par(mfrow=c(1,3),bg='transparent')
# plot(gron_eq_cat_upper$V, gron_eq_cat_upper$Magnitude, xlab = "V", ylab = "Magnitude", main = "Upper region")
# plot(log(gron_eq_cat_upper$V), gron_eq_cat_upper$Magnitude, xlab = "log(V)", ylab = "Magnitude", main = "Upper region")
# plot(sqrt(gron_eq_cat_upper$V), gron_eq_cat_upper$Magnitude, xlab = "sqrt(V)", ylab = "Magnitude", main = "Upper region")
# 
# gron_eq_cat_lower <- gron_eq_cat[gron_eq_cat$Northing <= (735000 - 0.6*gron_eq_cat$Easting),]
# dev.new(width=30, height=10,noRStudioGD = TRUE)
# par(mfrow=c(1,3),bg='transparent')
# plot(gron_eq_cat_lower$V, gron_eq_cat_lower$Magnitude, xlab = "V", ylab = "Magnitude", main = "Lower region")
# plot(log(gron_eq_cat_lower$V), gron_eq_cat_lower$Magnitude, xlab = "log(V)", ylab = "Magnitude", main = "Lower region")
# plot(sqrt(gron_eq_cat_lower$V), gron_eq_cat_lower$Magnitude, xlab = "sqrt(V)", ylab = "Magnitude", main = "Lower region")

#ACF and PACF plots checking for serial dependence

acf(gron_eq_cat$Magnitude)
pacf(gron_eq_cat$Magnitude)


# Exploratory analysis on third_nearest_dist ------------------------------

# Construct a grid within the groningen polygon

num_grid_points <- 100
grid <- expand.grid(X=seq(min(gron_polygon$POINT_X), max(gron_polygon$POINT_X), length.out=num_grid_points), Y=seq(min(gron_polygon$POINT_Y), max(gron_polygon$POINT_Y), length.out=num_grid_points))
grid_3d <- expand.grid(X=seq(min(gron_polygon$POINT_X), max(gron_polygon$POINT_X), length.out=num_grid_points), Y=seq(min(gron_polygon$POINT_Y), max(gron_polygon$POINT_Y), length.out=num_grid_points), Z=3)

#Check which grid points are in the polygon
grid_in_polygon <- grid[inpolygon(grid$X, grid$Y, gron_polygon$POINT_X, gron_polygon$POINT_Y),]
grid_in_outline <- grid[inpolygon(grid$X, grid$Y, gron_outline$X, gron_outline$Y),]

grid_in_polygon_3d <- grid_3d[inpolygon(grid_3d$X, grid_3d$Y, gron_polygon$POINT_X, gron_polygon$POINT_Y),]
grid_in_outline_3d <- grid_3d[inpolygon(grid_3d$X, grid_3d$Y, gron_outline$X, gron_outline$Y),]

# Find the distance to the third nearest geophone for each grid point and day
date_sequence_old <- sort(unique(c(geophones$Start_date, geophones$End_date)))
date_sequence <- sort(unique(geophones_deepest$Start_date, geophones_deepest$End_date))

min(as.Date(date_sequence_old))
max(as.Date(date_sequence_old))

mean_third_nearest_dist_polygon <- mean_third_distance_grid(grid_in_polygon, geophones, date_sequence)
mean_third_nearest_dist_outline <- mean_third_distance_grid(grid_in_outline, geophones, date_sequence)

loc1 <- data.frame(X=245000, Y=600000)
loc2 <- data.frame(X=255000, Y=595000)
loc3 <- data.frame(X=260000, Y=580000)
loc4 <- data.frame(X=245000, Y=580000)
loc5 <- data.frame(X=235000, Y=615000)

dist1 <- loc_third_distance(loc1, geophones, date_sequence)
dist2 <- loc_third_distance(loc2, geophones, date_sequence)
dist3 <- loc_third_distance(loc3, geophones, date_sequence)
dist4 <- loc_third_distance(loc4, geophones, date_sequence)
dist5 <- loc_third_distance(loc5, geophones, date_sequence)

third_nearest_dist_summaries_polygon <- mean_min_max_third_distance_grid_3d(grid_in_polygon_3d, geophones, date_sequence)
third_nearest_dist_summaries_outline <- mean_min_max_third_distance_grid_3d(grid_in_outline_3d, geophones, date_sequence)

loc1 <- data.frame(X=245000, Y=600000, Z=3)
loc2 <- data.frame(X=255000, Y=595000, Z=3)
loc3 <- data.frame(X=260000, Y=580000, Z=3)
loc4 <- data.frame(X=245000, Y=580000, Z=3)
loc5 <- data.frame(X=235000, Y=615000, Z=3)

dist1_3d <- loc_third_distance_3d(loc1, geophones, date_sequence)
dist2_3d <- loc_third_distance_3d(loc2, geophones, date_sequence)
dist3_3d <- loc_third_distance_3d(loc3, geophones, date_sequence)
dist4_3d <- loc_third_distance_3d(loc4, geophones, date_sequence)
dist5_3d <- loc_third_distance_3d(loc5, geophones, date_sequence)

dev.new(width=20, height=10,noRStudioGD = TRUE)
par(mfrow=c(1,2),bg='transparent')

plot(gron_polygon$POINT_X, gron_polygon$POINT_Y, type='l',col="grey",asp=1, xlab="Easting (m)", ylab = "Northing (m)")
lines(gron_outline$X, gron_outline$Y, col="blue")
text(loc1$X, loc1$Y, "1", cex=2, col="red")
text(loc2$X, loc2$Y, "2", cex=2, col="green")
text(loc3$X, loc3$Y, "3", cex=2, col="purple")
text(loc4$X, loc4$Y, "4", cex=2, col="orange")
text(loc5$X, loc5$Y, "5", cex=2, col="yellow")
legend("topright", legend=c("Groningen polygon", "Groningen outline"), col=c("grey", "blue"), lty=1)

# min_plot <- min(c(mean_third_nearest_dist_polygon, mean_third_nearest_dist_outline, dist1, dist2, dist3, dist4, dist5))
# max_plot <- max(c(mean_third_nearest_dist_polygon, mean_third_nearest_dist_outline, dist1, dist2, dist3, dist4, dist5))
# plot(as.Date(date_sequence), mean_third_nearest_dist_polygon, type='l', lwd=2, xlab="Date", ylab="Distance to third nearest geophone", ylim=c(min_plot, max_plot))
# lines(as.Date(date_sequence), mean_third_nearest_dist_outline, col="blue", lwd=2)
# lines(as.Date(date_sequence), dist1, col="red")
# lines(as.Date(date_sequence), dist2, col="green")
# lines(as.Date(date_sequence), dist3, col="purple")
# lines(as.Date(date_sequence), dist4, col="orange")
# lines(as.Date(date_sequence), dist5, col="yellow")
# legend("topright", legend=c("Mean distance (100x100 grid) in polygon", "Mean distance (100x100 grid) in outline"), col=c("black", "blue"), lty=1, lwd=2)

min_plot <- min(c(third_nearest_dist_summaries_polygon$min, third_nearest_dist_summaries_outline$min, dist1_3d, dist2_3d, dist3_3d, dist4_3d, dist5_3d))
max_plot <- max(c(third_nearest_dist_summaries_polygon$max, third_nearest_dist_summaries_outline$max, dist1_3d, dist2_3d, dist3_3d, dist4_3d, dist5_3d))
plot(as.Date(date_sequence), third_nearest_dist_summaries_polygon$mean, type='l', lwd=2, xlab="Date", ylab="Distance to third nearest geophone", ylim=c(min_plot, max_plot))
# lines(as.Date(date_sequence), third_nearest_dist_summaries_polygon$min, lwd=2, lty="dashed")
# lines(as.Date(date_sequence), third_nearest_dist_summaries_polygon$max, lwd=2, lty="dashed")
lines(as.Date(date_sequence), third_nearest_dist_summaries_outline$mean, col="blue", lwd=2)
# lines(as.Date(date_sequence), third_nearest_dist_summaries_outline$min, col="blue", lwd=2, lty="dashed")
# lines(as.Date(date_sequence), third_nearest_dist_summaries_outline$max, col="blue", lwd=2, lty="dashed")
lines(as.Date(date_sequence), dist1_3d, col="red")
lines(as.Date(date_sequence), dist2_3d, col="green")
lines(as.Date(date_sequence), dist3_3d, col="purple")
lines(as.Date(date_sequence), dist4_3d, col="orange")
lines(as.Date(date_sequence), dist5_3d, col="yellow")
legend("topright", legend=c("Mean distance in polygon", "Mean distance in outline"), col=c("black", "blue"), lty=1, lwd=2)
