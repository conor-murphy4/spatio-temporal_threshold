
# Calculate 3-dimensional Euclidean distance between earthquakes and geophones
euclid_dist <- function(eq_locations, geo_coords){
  
  euclid_km_2d <- function(x,z){
    sqrt( ((x[[1]] - z[[1]])/1000)^2 + ((x[[2]] - z[[2]])/1000)^2 )
    }
  
  dist <- apply(eq_locations, MARGIN =  1, FUN = euclid_km_2d , z = geo_coords)
  return(dist)
}

# Calculate 3-dimensional Euclidean distance between earthquakes and geophones
euclid_dist_3d <- function(eq_locations, geo_coords){
  
  euclid_km_3d <- function(x,z){
    sqrt(((x[[1]]-z[[1]])/1000)^2 + ((x[[2]]-z[[2]])/1000)^2 + (x[[3]]-z[[3]]/1000)^2)
    }
  
  dist <- apply(eq_locations, MARGIN = 1, FUN = euclid_km_3d, z=geo_coords)
  return(dist)
}


# Helper function to provide full network of geophones at a given date
geo_current_full_output <- function(geophones, date){
  indices_current <- which(geophones$Start_date <= date & geophones$End_date >= date)
  current_geophones <- geophones[indices_current,]
  return(current_geophones)
}

# Look up location, depth, start- and end-dates for geophones active on given date
geo_current <- function(geophones, date){
  
  current_geophones <- geo_current_full_output(geophones, date)
  
  current_geo_coords <- data.frame(
    X = current_geophones$Xcoord,
    Y = current_geophones$Ycoord,
    Z = current_geophones$Depth,
    Start_date = current_geophones$Start_date,
    End_date = current_geophones$End_date)
  
  return(current_geo_coords)
}


# Calculate distance to ith nearest geophone for earthquake locations in eq_cat 
# (index defines i)
distance_to_nearest <- function(eq_cat, geophones, index = 3){
  
  index_distances <- numeric(nrow(eq_cat))
  
  for (i in 1:nrow(eq_cat)) {
    eq_cat_current <- eq_cat[i,]
    current_geo_coords <- geo_current(geophones, eq_cat_current$Date)
    eq_loc_current <- data.frame(X = eq_cat_current$Easting,
                                 Y = eq_cat_current$Northing,
                                 Z = eq_cat_current$Depth)
    distances <- euclid_dist_3d(eq_loc_current, current_geo_coords)
    index_distances[i] <- sort(distances)[index]
  }
  
  return(index_distances)
}

