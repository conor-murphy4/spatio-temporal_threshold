
#Functions for third distance calculations

euclid_dist<- function(eq_locations, geo_coords){
  dist <- apply(eq_locations, 1, function(x,z){sqrt(((x[[1]]-z[[1]])/1000)^2+((x[[2]]-z[[2]])/1000)^2)}, z=geo_coords)
  return(dist)
}

euclid_dist_3d <- function(eq_locations, geo_coords){
  dist <- apply(eq_locations, 1, function(x,z){sqrt(((x[[1]]-z[[1]])/1000)^2+((x[[2]]-z[[2]])/1000)^2+(x[[3]]-z[[3]]/1000)^2)}, z=geo_coords)
  return(dist)
}

locs <- data.frame(X=runif(10,250000,270000), Y=runif(10,550000,610000), Z=runif(10, 0,3))
coords <- data.frame(X=runif(5,230000,275000), Y=runif(5,520000,625000), Z=runif(5, 0,3))
euclid_dist_3d(locs, coords)

geo_current <- function(geophones, date){
  indices_current <- which(geophones$Start_date <= date & geophones$End_date >= date)
  current_geophones <- geophones[indices_current,]
  current_geo_coords <- data.frame(X=current_geophones$Xcoord, Y=current_geophones$Ycoord, Z=current_geophones$Depth, Start_date=current_geophones$Start_date, End_date=current_geophones$End_date)
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

distance_to_third_nearest_3d <- function(eq_cat, geophones){
  third_distances <- numeric(nrow(eq_cat))
  for(i in 1:nrow(eq_cat)){
    eq_cat_current <- eq_cat[i,]
    current_geo_coords <- geo_current(geophones,eq_cat_current$Date)
    eq_loc_current <- data.frame(X=eq_cat_current$Easting, Y=eq_cat_current$Northing, Z=eq_cat_current$Depth)
    distances <- euclid_dist_3d(eq_loc_current, current_geo_coords)
    third_distances[i] <- sort(distances)[3]
  }
  return(third_distances)
}

mean_third_distance_grid <- function(grid, geophones, dates){
  mean_third_distances <- numeric(length(dates))
  for(i in 1:length(dates)){
    current_date <- dates[i]
    current_geo_coords <- geo_current(geophones,current_date)
    distances <- euclid_dist(grid, current_geo_coords)
    third_distances <- apply(distances, 2, function(x){sort(x)[3]})
    mean_third_distances[i] <- mean(third_distances)
  }
  return(mean_third_distances)
}

mean_min_max_third_distance_grid_3d <- function(grid, geophones, dates){
  #Check grid has 3 coordinates
  if(ncol(grid)!=3) stop("Grid must have 3 coordinates")
  
  mean_third_distances <- numeric(length(dates))
  min_third_distances <- numeric(length(dates))
  max_third_distances <- numeric(length(dates))
  for(i in 1:length(dates)){
    current_date <- dates[i]
    current_geo_coords <- geo_current(geophones,current_date)
    distances <- euclid_dist_3d(grid, current_geo_coords)
    third_distances <- apply(distances, 2, function(x){sort(x)[3]})
    mean_third_distances[i] <- mean(third_distances)
    min_third_distances[i] <- min(third_distances)
    max_third_distances[i] <- max(third_distances)
  }
  return(list(mean=mean_third_distances, min=min_third_distances, max=max_third_distances))
}

loc_third_distance <- function(loc, geophones, dates){
  third_distances <- numeric(length(dates))
  for(i in 1:length(dates)){
    current_date <- dates[i]
    current_geo_coords <- geo_current(geophones,current_date)
    distances <- euclid_dist(loc, current_geo_coords)
    third_distances[i] <- apply(distances, 2, function(x){sort(x)[3]})
  }
  return(third_distances)
}

loc_third_distance_3d <- function(loc, geophones, dates){
  #Check locations has 3 coordinates
  if(ncol(loc)!=3) stop("Location must have 3 coordinates")
  
  third_distances <- numeric(length(dates))
  for(i in 1:length(dates)){
    current_date <- dates[i]
    current_geo_coords <- geo_current(geophones,current_date)
    distances <- euclid_dist_3d(loc, current_geo_coords)
    third_distances[i] <- apply(distances, 2, function(x){sort(x)[3]})
  }
  return(third_distances)
}

What_geos <- function(Eq.df, GeoCoord){
  Geo_final <- data.frame(X   = numeric(length(Eq.df$X)), Y = numeric(length(Eq.df$X)), Start_date= character(length(Eq.df$X)), End_date= character(length(Eq.df$X)))
  for(i in 1:length(Eq.df$X)){
    Eq_current <- Eq.df[i,]
    Current_coords <- Geo_current(GeoCoord,Eq_current$Date)
    Eq_loc_current <- data.frame(X=Eq_current$X, Y=Eq_current$Y)
    Geo_loc_current <- data.frame(X=Current_coords$Xcoord, Y=Current_coords$Ycoord)
    dist <- sq_eu_dist(Eq_loc_current,Geo_loc_current)
    df <- data.frame(D=dist,CX=Current_coords$Xcoord, CY=Current_coords$Ycoord, Start_date = Current_coords$Start_date, End_date = Current_coords$End_date)
    df <- df[order(df$D),]
    Geo_final[i,1] <- df$CX[3]
    Geo_final[i,2] <- df$CY[3]
    Geo_final[i,3] <- df$Start_date[3]
    Geo_final[i,4] <- df$End_date[3]
  }
  return(Geo_final)
}

what_geos <- function(eq_cat, geophones){
  included_geos <- data.frame(X=numeric(nrow(eq_cat)), Y=numeric(nrow(eq_cat)), Z=numeric(nrow(eq_cat)), Start_date=character(nrow(eq_cat)), End_date=character(nrow(eq_cat)))
  print(dim(included_geos))
  for(i in 1:nrow(eq_cat)){
    print(i)
    eq_cat_current <- eq_cat[i,]
    current_geo_coords <- geo_current(geophones,eq_cat_current$Date)
    eq_loc_current <- data.frame(X=eq_cat_current$Easting, Y=eq_cat_current$Northing, Z=eq_cat_current$Depth)
    distances <- euclid_dist_3d(eq_loc_current, current_geo_coords)
    current_geo_ordered <- current_geo_coords[order(distances),]
    print(dim(current_geo_ordered[3,]))
    included_geos[i,] <- current_geo_ordered[3,]
  }
  return(included_geos)
}

