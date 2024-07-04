Geo <- read.csv("Data/Geophones/Geophone data/Geophone_rect.csv", header=T)

library(pracma)

#Simple Polygon outline
Gron_polygon <- read.csv('Data/Geophones/Geophone data/groningen_field_simple_polygon.csv', header=T)

#Gron outline
Gron_outline <- read.csv('Data/Geophones/Geophone data/Groningen_Field_outline.csv', header=T)

Geo_rect <- data.frame(X=Geo$Xcoord, Y=Geo$Ycoord)
Geo_polygon<-inpolygon(Geo$Xcoord, Geo$Ycoord, Gron_outline$X, Gron_outline$Y)
Geo_out <- Geo[Geo_polygon,]
Geo_poly <- data.frame(X=Geo_out$Xcoord, Y=Geo_out$Ycoord)


#Distance functions
sq_eu_dist<- function(Eq_loc, GeoCoord){
  dist <- apply(Eq_loc, 1, function(x,z){(x[[1]]-z[[1]])^2+(x[[2]]-z[[2]])^2}, z=GeoCoord)
  return(dist)
}

V_2d <- function(Eq_loc, GeoCoord){
  dist <- sq_eu_dist(Eq_loc, GeoCoord)
  V <- apply(dist, 2, function(x){ sort(x)[3]}) 
  return(V)
}

#Approximating integrating over whole rectangle
l<-100
x <- seq(232000, 270000, length=380)
y <- seq(570000, 615000, length=450)
grid <- expand.grid(X=x, Y=y)
length(unique(grid$X))

grid_poly <- inpolygon(grid$X,grid$Y, Gron_outline$X, Gron_outline$Y)
grid_polygon <- grid[grid_poly,]

V_rect <- V_2d(grid, Geo_rect)
V_df <- data.frame(V=V_rect, X=grid$X, Y=grid$Y)

#grid size
grid_size <- (615000-570000)*(270000-232000)/10000
app_int_rect <- grid_size*sum(V_df$V) 

#Approximating over whole gas field
V_whole <- V_2d(grid_polygon, Geo_poly)

grid_size <- (max(Gron_outline$Y)-min(Gron_outline$Y))*(max(Gron_outline$X)-min(Gron_outline$X))/(length(unique(grid_polygon$X))*length(unique(grid_polygon$Y))) 
app_int_whole <- grid_size*sum(V_whole)

#Approximating integral of V(x) for each year

#Working out grid size


grid_size <- 100 # (max(Gron_outline$Y)-min(Gron_outline$Y))*(max(Gron_outline$X)-min(Gron_outline$X))/(length(unique(grid_polygon$X))*length(unique(grid_polygon$Y))) 

approx_int <- function(year, grid_size){
  ind <- which(Geo_out$Start_Year <= year & Geo_out$End_Year >= year)
  Geo_current <- Geo_out[ind,]
  Current_coords <- data.frame(X=Geo_current$Xcoord, Y=Geo_current$Ycoord)
  V_current <- V_2d(grid_polygon, Current_coords)
  int <- grid_size*sum(V_current)
  return(int)
}

year <- seq(1995, 2022, by=1)
V_int <- numeric(length(year))
for(i in 1:length(year)){
  current_year <- year[i]
  V_int[i] <- approx_int(current_year, grid_size)
}

plot(year, log(V_int), type='l', main="Integrated V", lwd=2)

#Fixing points in space and assessing changes over time

#First loc: (260000,580000)
loc1 <- data.frame(X=245000, Y=600000)
loc2 <- data.frame(X=255000, Y=595000)
loc3 <- data.frame(X=260000, Y=580000)
loc4 <- data.frame(X=245000, Y=580000)
V_time <- function(year, loc){
  ind <- which(Geo_out$Start_Year <= year & Geo_out$End_Year >= year)
  Geo_current <- Geo_out[ind,]
  Current_coords <- data.frame(X=Geo_current$Xcoord, Y=Geo_current$Ycoord)
  V_current <- V_2d(loc, Current_coords)
  return(V_current)
}

year <- seq(1995, 2022, by=1)
V_loc1 <- numeric(length(year))
V_loc2 <- numeric(length(year))
V_loc3 <- numeric(length(year))
V_loc4 <- numeric(length(year))
for(i in 1:length(year)){
  current_year <- year[i]
  V_loc1[i] <- V_time(current_year, loc1)
  V_loc2[i] <- V_time(current_year, loc2)
  V_loc3[i] <- V_time(current_year, loc3)
  V_loc4[i] <- V_time(current_year, loc4)
}

dev.new(width=9.17, height=5.53,noRStudioGD = TRUE)
par(mfrow=c(1,2),bg='transparent')

plot(Gron_outline$X, Gron_outline$Y, pch=19, cex=0.3, col="darkgrey",xlab="Easting (m)", ylab="Northing (m)")
points(Geo_poly$X, Geo_poly$Y, pch=19)
text(loc1$X, loc1$Y, "1", cex=2)
text(loc2$X, loc2$Y, "2", cex=2, col=2)
text(loc3$X, loc3$Y, "3", cex=2, col=3)
text(loc4$X, loc4$Y, "4", cex=2, col=4)

plot(year, log(V_int), type='l', main="Integrated V", lwd=2)
par(mfrow=c(1,1))
plot(year, log(V_loc1),  lwd=2, type='l', ylim = c(12, 20), xlab="Year", ylab="log(V(x,t))")
lines(year, log(V_loc2),  lwd=2, col=2)
lines(year, log(V_loc3),  lwd=2, col=3)
lines(year, log(V_loc4), lwd=2, col=4)
