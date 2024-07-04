geophones <- read.csv("Data/Geophones/Geophones_processed_03-07-2024.csv", header=T, row.names = 1)
gron_eq_cat <- read.csv("Data/Events/gr_earthquakes_after_geophone_start_20240104.csv", header=T, row.names = 1)

#Functions for third distance calculations

euclid_dist<- function(eq_locations, geo_coords){
  dist <- apply(eq_locations, 1, function(x,z){sqrt(((x[[1]]-z[[1]])/1000)^2+((x[[2]]-z[[2]])/1000)^2)}, z=geo_coords)
  return(dist)
}

geo_current <- function(geophones, date){
  indices_current <- which(geophones$Start_date <= date & geophones$End_date >= date)
  current_geophones <- geophones[indices_current,]
  current_geo_coords <- data.frame(X=current_geophones$Xcoord, Y=current_geophones$Ycoord)
  return(current_geo_coords)
}

distance_to_third_nearest <- function(eq_cat, geophones){
  third_distances <- numeric(nrow(eq_cat))
  for(i in 1:nrow(eq_cat)){
    eq_cat_current <- eq_cat[i,]
    current_geo_coords <- geo_current(geophones,eq_cat_current$Date)
    eq_loc_current <- data.frame(X=eq_cat_current$Easting, Y=eq_cat_current$Northing)
    distances <- euclid_dist(eq_loc_current, current_geo_coords)
    third_distances[i] <- sort(distances)[3]
  }
  return(third_distances)
}

#####---- Move below to EDA??? ---#######

gron_eq_cat$V <- distance_to_third_nearest(gron_eq_cat, geophones)

dev.new(width=30, height=10,noRStudioGD = TRUE)
par(mfrow=c(1,3),bg='transparent')
plot(gron_eq_cat$V, gron_eq_cat$Magnitude, xlab = "V", ylab = "Magnitude", main = "Whole catalogue")
plot(log(gron_eq_cat$V), gron_eq_cat$Magnitude, xlab = "log(V)", ylab = "Magnitude", main = "Whole catalogue")
plot(sqrt(gron_eq_cat$V), gron_eq_cat$Magnitude, xlab = "sqrt(V)", ylab = "Magnitude", main = "Whole catalogue")


#Splitting gas field in two
#Plotting
dev.new(width=10, height=10,noRStudioGD = TRUE)
par(mfrow=c(1,1), bg='transparent')
plot(gron_outline$X, gron_outline$Y, pch=19, col="blue",asp=1, xlab="Easting (m)", ylab = "Northing (m)")
points(gron_eq_cat$Easting, gron_eq_cat$Northing, )
y_div <- 735000 - 0.6*seq(220000, 280000, by=1)
points(seq(220000, 280000, by=1), y_div)

gron_eq_cat_upper <- gron_eq_cat[gron_eq_cat$Northing > (735000 - 0.6*gron_eq_cat$Easting),]
dev.new(width=30, height=10,noRStudioGD = TRUE)
par(mfrow=c(1,3),bg='transparent')
plot(gron_eq_cat_upper$V, gron_eq_cat_upper$Magnitude, xlab = "V", ylab = "Magnitude", main = "Upper region")
plot(log(gron_eq_cat_upper$V), gron_eq_cat_upper$Magnitude, xlab = "log(V)", ylab = "Magnitude", main = "Upper region")
plot(sqrt(gron_eq_cat_upper$V), gron_eq_cat_upper$Magnitude, xlab = "sqrt(V)", ylab = "Magnitude", main = "Upper region")

gron_eq_cat_lower <- gron_eq_cat[gron_eq_cat$Northing <= (735000 - 0.6*gron_eq_cat$Easting),]
dev.new(width=30, height=10,noRStudioGD = TRUE)
par(mfrow=c(1,3),bg='transparent')
plot(gron_eq_cat_lower$V, gron_eq_cat_lower$Magnitude, xlab = "V", ylab = "Magnitude", main = "Lower region")
plot(log(gron_eq_cat_lower$V), gron_eq_cat_lower$Magnitude, xlab = "log(V)", ylab = "Magnitude", main = "Lower region")
plot(sqrt(gron_eq_cat_lower$V), gron_eq_cat_lower$Magnitude, xlab = "sqrt(V)", ylab = "Magnitude", main = "Lower region")
