setwd("C:/Users/murphyc4/OneDrive - Lancaster University/STOR-i/PhD/Projects/Spatio-temporal") 

library(pracma) # Use inpolygon function if needed
library(lubridate)
library(stringr)

#Geophones
Geophones <- read.csv("Data/Geophones/Geophone data/Geo_orig_post_process.csv")

#Adding dates to Geophone csv

Geophones$Start_date <- str_sub(Geophones$StartTime, 1,10)
Geophones$End_date <- str_sub(Geophones$EndTime,1,10)

#Gron outline
Gron_outline <- read.csv('Data/Geophones/Geophone data/Groningen_Field_outline.csv', header=T)
plot(Gron_outline$X, Gron_outline$Y, xlim=c(230000, 275000), ylim = c(560000, 620000))

#squared distance in km
sq_eu_dist<- function(Eq_loc, GeoCoord){
  dist <- apply(Eq_loc, 1, function(x,z){(x[[1]]/1000-z[[1]]/1000)^2+(x[[2]]/1000-z[[2]]/1000)^2}, z=GeoCoord)
  return(dist)
}

#Current set of geophones according to date
Geo_current <- function(Geo, date){
  ind <- which(Geo$Start_date <= date & Geo$End_date >= date)
  Geo_current <- Geo[ind,]
  #Current_coords <- data.frame(X=Geo_current$Xcoord, Y=Geo_current$Ycoord)
  return(Geo_current)#Current_coords)
}

#Spatio-temporal V function
V_2d_time <- function(Eq.df, GeoCoord){
  V <- numeric(length(Eq.df$X))
  for(i in 1:length(Eq.df$X)){
    Eq_current <- Eq.df[i,]
    Current_coords <- Geo_current(GeoCoord,Eq_current$Date)
    Eq_loc_current <- data.frame(X=Eq_current$X, Y=Eq_current$Y)
    dist <- sq_eu_dist(Eq_loc_current, Current_coords)
    V[i] <- sort(dist)[3]
  }
  return(V)
}

#Earthquake catalogue dataframe (now has V in km^2 already included)
load("Data/Catalogue.Rda")

#Calculating V for Groningen EQs
Cat.df$V <- V_2d_time(Cat.df, Geophones)
save(Cat.df, file="Data/Catalogue.Rda")

#Groningen dataframe
load("Data/Catalogue.Rda")

###------------How big is the region including all geophones used in V???---------###

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

Geo_final <- What_geos(Cat.df, Geophones)
write.csv(Geo_final, "Data/Geophones/Geophone data/Geo_final.csv")
plot(Geo_final[,1], Geo_final[,2], col="blue", pch=19, ylim = c(535000, 620000))
points(Gron_outline$X, Gron_outline$Y, col="red", pch=19)
points(Cat.df$X, Cat.df$Y)

#Can do some checks to make sure this is correct? Maybe do a quick simulated version??
#Need to estimate size of rectangle now...

rect_X <- c(227000,227000,272000,272000,227000)
rect_Y <- c(535000,614000,614000,535000,535000)
lines(rect_X, rect_Y,lwd=2)

#Rectangle should therefore be 45km x 79km

small_rect_X <- c(255000,255000,260000,260000,255000)
small_rect_Y <- c(590000,595000,595000,590000,590000)
lines(small_rect_X, small_rect_Y,lwd=2, col="green")


ind <- inpolygon(Cat.df$X, Cat.df$Y, small_rect_X, small_rect_Y)
Cat.df[ind,]
Geos<-What_geos(Cat.df[ind,], Geophones)
Geos[which(Geos$Y < 545000),]

sq_eu_dist(Cat.df[1419,3:4], Geos[which(Geos$Y < 545000),])

points(Cat.df$X[1419], Cat.df$Y[1419], col="orange", cex=3, pch=19)
points(Geos$X[which(Geos$Y < 545000)], Geos$Y[which(Geos$Y < 545000)], col="orange", cex=4, lwd=2)

plot(Geo_final$X[which(Geo_final$Start_date <= "1995-01-24" & Geo_final$End_date < "1995-01-24")],Geo_final$Y[which(Geo_final$Start_date <= "1995-01-24" & Geo_final$End_date < "1995-01-24")], ylab="Y", xlab="X")


