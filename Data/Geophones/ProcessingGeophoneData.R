library("RColorBrewer")
library("sp")
require("jsonlite")
library("lubridate")
library("matrixStats")
library("mgcv")
library("KScorrect")
library("stringr")

#specifying coordinate systems
latlong = "+init=epsg:4326"
google = "+init=epsg:3857"
RD = "+init=epsg:28992"

GeophoneData<-read.csv("Data/Geophones/Geophones_unprocessed_03-07-2024.csv")
idObs=!is.nan(GeophoneData$Longitude)
tpoints <- cbind(as.numeric(GeophoneData$Longitude[idObs]), as.numeric(GeophoneData$Latitude[idObs]))
tpoints2 <- SpatialPointsDataFrame(tpoints, proj4string = CRS(latlong), data=data.frame(ID=1:nrow(tpoints)))
tpoints3 <- spTransform(tpoints2, CRS(RD))

# creating coordinate system
RD_coord=as.data.frame(matrix(NA,nrow(GeophoneData),ncol=2))
RD_coord[idObs,1]=tpoints3@coords[,1]
RD_coord[idObs,2]=tpoints3@coords[,2]
names(RD_coord)=c("Xcoord","Ycoord")

#Adding RD coords
Geo<-as.data.frame(cbind(GeophoneData,RD_coord))

write.csv(Geo,"Data/Geophones/Geophones_with_RDcoord_03-07-2024.csv")

#Cleaning geophone data
Geo <- Geo[,-c(3, seq(9,15,by=1))]

Geo$Start_date <- dmy(str_sub(Geo$StartTime, 1,10))
Geo$End_date <- dmy(str_sub(Geo$EndTime,1,10))

Geo$Channel3 <- str_sub(Geo$Channel,3,3)
E<- Geo[Geo$Channel3=="E",]
N<- Geo[Geo$Channel3=="N",]
Z<- Geo[Geo$Channel3=="Z",]
H1 <- Geo[Geo$Channel=="H1",]
H2 <- Geo[Geo$Channel=="H2",]
HG1 <- Geo[Geo$Channel=="HG1",]
HG2 <- Geo[Geo$Channel=="HG2",]

Geo_fin <- rbind(E,N,Z,H1,H2,HG1,HG2)

Geo_fin <- Geo_fin[,-ncol(Geo_fin)]

write.csv(Geo_fin, "Data/Geophones/Geophones_processed_03-07-2024.csv")

#Removing duplicate geophones
indices_to_remove <- c()
for(i in 1:(nrow(Geo_fin)-1)){
  for(j in (i+1):nrow(Geo_fin)){
    if(i != j){
      check_coords <- c(Geo_fin$Xcoord[i] == Geo_fin$Xcoord[j], Geo_fin$Ycoord[i] == Geo_fin$Ycoord[j])
      check_channel <- Geo_fin$Channel[i] == Geo_fin$Channel[j]
      check_station <- Geo_fin$Station[i] == Geo_fin$Station[j]
      check_depth <- Geo_fin$Depth[i] == Geo_fin$Depth[j]
      check_elevation <- Geo_fin$Elevation[i] == Geo_fin$Elevation[j]
      check_dates <- Geo_fin$Start_date[i] == Geo_fin$Start_date[j] & Geo_fin$End_date[i] == Geo_fin$End_date[j]
      if(all(check_coords) & check_channel & check_station & check_depth & check_elevation & check_dates){
        indices_to_remove <- c(indices_to_remove, j)
      }
    }
  }
}

indices_to_remove <- unique(indices_to_remove)

Geo_fin_without_duplicates <- Geo_fin[-indices_to_remove,]

write.csv(Geo_fin_without_duplicates, "Data/Geophones/Geophones_processed_03-07-2024_without_duplicates.csv")

Geo_fin_without_duplicates <- read.csv("Data/Geophones/Geophones_processed_03-07-2024_without_duplicates.csv")
#Taking only deepest of each station
Geo_fin_deepest <- data.frame(matrix(ncol=ncol(Geo_fin_without_duplicates),nrow=0))
colnames(Geo_fin_deepest) <- colnames(Geo_fin_without_duplicates)

for(x in unique(Geo_fin_without_duplicates$Xcoord)){
  indices_x <- which(Geo_fin_without_duplicates$Xcoord == x)
  Geo_fin_x <- Geo_fin_without_duplicates[indices_x,]
  print(paste("Number of geophones at x=",x," is ",nrow(Geo_fin_x)))
  for(y in unique(Geo_fin_x$Ycoord)){
    indices_y <- which(Geo_fin_x$Ycoord == y)
    Geo_fin_y <- Geo_fin_x[indices_y,]
    print(paste("Number of geophones at x= ", x, "and y=",y," is ",nrow(Geo_fin_y)))
    for(date in unique(Geo_fin_y$Start_date)){
      indices_date <- which(Geo_fin_y$Start_date == date)
      Geo_fin_date <- Geo_fin_y[indices_date,]
      print(paste("Number of geophones at x= ", x, "and y=",y," and date=",date," is ",nrow(Geo_fin_date)))
      depths <- Geo_fin_date$Depth
      max_depth <- max(depths)
      ind_deepest <- which(depths == max_depth)
      if(length(ind_deepest) > 1){
        ind_deepest <- ind_deepest[1]
      }
      Geo_fin_deepest <- rbind(Geo_fin_deepest,Geo_fin_date[ind_deepest,])
    }
  }  
}

write.csv(Geo_fin_deepest, "Data/Geophones/Geophones_processed_03-07-2024_deepest_only.csv")

Gron_outline <- read.csv('Data/Geophones/Groningen_Field_outline.csv', header=T)

plot(Geo$Xcoord,Geo$Ycoord)
points(Gron_outline$X, Gron_outline$Y, col="red", cex=0.4)
points(Geo_fin_deepest$Xcoord,Geo_fin_deepest$Ycoord, pch=19, col="blue")


