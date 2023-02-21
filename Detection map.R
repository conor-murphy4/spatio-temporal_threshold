sq_eu_dist<- function(Eq_loc, GeoCoord){
  dist <- apply(Eq_loc, 1, function(x,z){(x[[1]]-z[[1]])^2+(x[[2]]-z[[2]])^2}, z=GeoCoord)
  #Dist <- matrix(dist,byrow=F, length(Eq_loc[[1]]))
  return(dist)
}
third_dist <- function(x){
  min_ind <- which.min(x)
  x_2 <- x[-min_ind]
  min_ind_2 <- which.min(x_2)
  x_3 <- x_2[-min_ind_2]
  min_ind_3 <- which.min(x_3)
  return(x_3[min_ind_3])
}
V_2d <- function(Eq_loc, GeoCoord){
  dist <- sq_eu_dist(Eq_loc, GeoCoord)
  V <- apply(dist, 2, third_dist) 
  return(V)
}

#Reading in Geophones within rectangle data
setwd("C:/Users/murphyc4/OneDrive - Lancaster University") 
Geo <- read.csv("STOR-i/PhD/Data/Geophones/Geophone data/Geophone_rect.csv", header=T)
plot(Geo$Xcoord, Geo$Ycoord)

#Simple Polygon outline
Gron_polygon <- read.csv('STOR-i/PhD/Data/Geophones/Geophone data/groningen_field_simple_polygon.csv', header=T)

#Gron outline
Gron_outline <- read.csv('STOR-i/PhD/Data/Geophones/Geophone data/Groningen_Field_outline.csv', header=T)

#Logical check for polygon
library(pracma)
Geo_polygon<-inpolygon(Geo$Xcoord, Geo$Ycoord, Gron_outline$X, Gron_outline$Y)
Geo_out <- Geo[Geo_polygon,]
Geo_poly <- data.frame(X=Geo_out$Xcoord, Y=Geo_out$Ycoord)

#Grid across rectangle
l<-50
x <- seq(232000, 270000, length=l)
y <- seq(570000, 615000, length=l)
grid <- expand.grid(X=x, Y=y)
#head(sq_eu_dist(grid, Geo_poly))

#Filled contour 
V_rect <- V_2d(grid, Geo_poly)
V_df <- data.frame(V=V_rect, X=grid$X, Y=grid$Y)
V_df$poly<-inpolygon(V_df$X,V_df$Y, Gron_outline$X, Gron_outline$Y)
V_df$V[!V_df$poly] <- NA
V_mat <- matrix(V_df$V, byrow=T, nrow=l)
filled.contour(sort(unique(V_df$X)), sort(unique(V_df$Y)), t(V_mat),xlim=c(230000, 270000), ylim=c(565000, 615000), plot.axes={points(Geo_poly$Xcoord, Geo_poly$Ycoord)})
#lines(Gron_outline$X-3500, Gron_outline$Y)
#points(Geo_poly$Xcoord, Geo_poly$Ycoord, pch=19, col="black")

#Snapshots in time
#2013
ind_2013 <- which(Geo_out$Start_Year <= 2013 & Geo_out$End_Year >= 2013)
Geo_2013 <- Geo_out[ind_2013,]
Coords_2013 <- data.frame(X=Geo_2013$Xcoord, Y=Geo_2013$Ycoord)
V_2013 <- V_2d(grid, Coords_2013)
V_df_2013 <- data.frame(V=V_2013, X=grid$X, Y=grid$Y)
V_df_2013$outline<-inpolygon(V_df_2013$X,V_df_2013$Y, Gron_outline$X, Gron_outline$Y)
V_df_2013$V[!V_df_2013$outline] <- NA
V_mat_2013 <- matrix(V_df_2013$V, byrow=T, nrow=l)

dev.new(width=9.17, height=5.53,noRStudioGD = TRUE)
par(mfrow=c(1,1),bg='transparent')
filled.contour(sort(unique(V_df_2013$X)), sort(unique(V_df_2013$Y)), t(V_mat_2013)/1000,xlim=c(230000, 270000), ylim=c(565000, 615000), plot.axes={points(Geo_2013$Xcoord, Geo_2013$Ycoord)})
# lines(Gron_outline$X-3500, Gron_outline$Y)

#2014
ind_2014 <- which(Geo_out$Start_Year == 2014 & Geo_out$End_Year >= 2014)
Geo_14 <- Geo_out[ind_2014,]
Geo_2014 <- rbind(Geo_2013, Geo_14)
Coords_2014 <- data.frame(X=Geo_2014$Xcoord, Y=Geo_2014$Ycoord)
V_2014 <- V_2d(grid, Coords_2014)
V_df_2014 <- data.frame(V=V_2014, X=grid$X, Y=grid$Y)
V_df_2014$outline<-inpolygon(V_df_2014$X,V_df_2014$Y, Gron_outline$X, Gron_outline$Y)
V_df_2014$V[!V_df_2014$outline] <- NA
V_mat_2014 <- matrix(V_df_2014$V, byrow=T, nrow=l)
filled.contour(sort(unique(V_df_2014$X)), sort(unique(V_df_2014$Y)), t(V_mat_2014),xlim=c(230000, 270000), ylim=c(565000, 615000), plot.axes={points(Geo_2014$Xcoord, Geo_2014$Ycoord)})

#2015
ind_2015 <- which(Geo_out$Start_Year == 2015 & Geo_out$End_Year >= 2015)
Geo_15 <- Geo_out[ind_2015,]
Geo_2015 <- rbind(Geo_2013,Geo_14, Geo_15)
Coords_2015 <- data.frame(X=Geo_2015$Xcoord, Y=Geo_2015$Ycoord)
V_2015 <- V_2d(grid, Coords_2015)
V_df_2015 <- data.frame(V=V_2015, X=grid$X, Y=grid$Y)
V_df_2015$outline<-inpolygon(V_df_2015$X,V_df_2015$Y, Gron_outline$X, Gron_outline$Y)
V_df_2015$V[!V_df_2015$outline] <- NA
V_mat_2015 <- matrix(V_df_2015$V, byrow=T, nrow=l)
filled.contour(sort(unique(V_df_2015$X)), sort(unique(V_df_2015$Y)), t(V_mat_2015)/1000,xlim=c(230000, 270000), ylim=c(565000, 615000), plot.axes={points(Geo_2015$Xcoord, Geo_2015$Ycoord)})

