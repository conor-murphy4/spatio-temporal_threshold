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

Gron_outline <- read.csv('Data/Geophones/Groningen_Field_outline.csv', header=T)

plot(Geo$Xcoord,Geo$Ycoord)
points(Gron_outline$X, Gron_outline$Y, col="red", cex=0.4)
points(Geo_fin$Xcoord,Geo_fin$Ycoord, pch=19, col="blue")
write.csv(Geo_fin, "Data/Geophones/Geophones_processed_03-07-2024.csv")

