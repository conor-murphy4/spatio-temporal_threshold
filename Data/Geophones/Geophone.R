
#------------------Reading in Geophone data-----------------------------
setwd("C:/Users/murphyc4/OneDrive - Lancaster University") #Windows

Geo <- read.csv('STOR-i/PhD/Projects/Spatio-temporal/Data/Geophones/Geophone data/Geophones_orig.csv', header=T)

Geo1 <- Geo[,-c(4,seq(8,16,by=1))]
Geo2 <- Geo1[-which(Geo1$Station=="<undefined>"),]

library(plotrix)
plot(Geo2$Xcoord,Geo2$Ycoord )#,ylim=c(550000,610000), xlim=c(230000,260000))
points(Gron_outline$X, Gron_outline$Y, col="red", cex=0.4)

Gron_outline <- read.csv('STOR-i/PhD/Projects/Spatio-temporal/Data/Geophones/Geophone data/Groningen_Field_outline.csv', header=T)
Gron_polygon <- read.csv('STOR-i/PhD/Projects/Spatio-temporal/Data/Geophones/Geophone data/groningen_field_simple_polygon.csv', header=T)


plot(Gron_outline$X,Gron_outline$Y, col="red")#, xlim=c(230000,270000), ylim=c(560000,620000))
points(Geo2$Xcoord,Geo2$Ycoord)
lines(Gron_polygon$X, Gron_polygon$Y, col="blue")


X <- c(232000,232000,268000,268000,232000)
Y <- c(568000,614000,614000,568000,568000)

Gron_rect <- data.frame(X=X, Y=Y)
lines(Gron_rect$X, Gron_rect$Y, col="black", lwd=2)

# Geo_gron_rect <- Geo2[which(Geo2$Xcoord<max(Gron_rect$X) & Geo2$Xcoord > min(Gron_rect$X) & Geo2$Ycoord < max(Gron_rect$Y) & Geo2$Ycoord>min(Gron_rect$Y)),]
# points(Geo_gron_rect$Xcoord,Geo_gron_rect$Ycoord, pch=19)



#which(Geo2$Xcoord<max(Gron_rect$X) & Geo2$Xcoord > min(Gron_rect$X) & Geo2$Ycoord < max(Gron_rect$Y) & Geo2$Ycoord>min(Gron_rect$Y))


unique(HG1$Channel)
library(stringr)
Geo2$Channel3 <- str_sub(Geo2$Channel,3,3)
Geo2$Start_Year <- str_sub(Geo2$StartTime, 1,4)
Geo2$End_Year <- str_sub(Geo2$EndTime,1,4)

library(dplyr)
E<-filter(Geo2, grepl("E", Channel3))
N<-filter(Geo2, grepl("N", Channel3))
Z<-filter(Geo2, grepl("Z", Channel3))
H1 <- filter(Geo2, grepl("H1",Channel))
H2 <- filter(Geo2, grepl("H2",Channel))
HG1 <-  filter(Geo2, grepl("HG1",Channel))
HG2 <- filter(Geo2, grepl("HG2",Channel))

Geo_fin <- rbind(E,N,Z,H1,H2,HG1,HG2)

points(Geo_fin$Xcoord,Geo_fin$Ycoord, pch=19, col="blue")

GeoCoord <- data.frame(X=Geo_fin$Xcoord,Y=Geo_fin$Ycoord)
#write.csv(Geo_fin, "Geophone_rect.csv")
write.csv(Geo_fin, "STOR-i/PhD/Projects/Spatio-temporal/Data/Geophones/Geophone data/Geo_orig_post_process.csv")
#Snapshots in time
head(Geo_fin)
unique(Geo_fin$Start_Year)

dev.new(width=9.17, height=5.53,noRStudioGD = TRUE)
par(mfrow=c(1,1),bg='transparent')

plot(Gron_outline$X,Gron_outline$Y, col="red", xlab='', ylab='', pch=19, cex=0.5, cex.axis=1.2)#, xlim=c(230000,270000), ylim=c(560000,620000))
points(Geo_fin$Xcoord[Geo_fin$Start_Year <= 2013 & Geo_fin$End_Year >= 2013], Geo_fin$Ycoord[Geo_fin$Start_Year <= 2013 & Geo_fin$End_Year>=2013], col="black", pch=4)
points(Geo_fin$Xcoord[Geo_fin$Start_Year == 2014 & Geo_fin$End_Year >= 2014], Geo_fin$Ycoord[Geo_fin$Start_Year == 2014 & Geo_fin$End_Year>=2014], col="blue", pch=4)
points(Geo_fin$Xcoord[Geo_fin$Start_Year == 2015 & Geo_fin$End_Year >= 2015], Geo_fin$Ycoord[Geo_fin$Start_Year == 2015 & Geo_fin$End_Year>=2015], col="green", pch=4)
points(Geo_fin$Xcoord[Geo_fin$End_Year<=2015], Geo_fin$Ycoord[Geo_fin$End_Year<=2015], col="red", cex=3)
lines(Gron_rect$X, Gron_rect$Y, col="black", lwd=2)
mtext('X-coordinate',side=1, cex=1.5,line=2.5)
mtext("Y-coordinate", side=2, cex=1.5, line=2.5)
legend("topright", c("<2013", "2014", "2015", "Removed"), col=c("black","blue","green","red"), pch=c(4,4,4,1), cex=1.3 )



####--------- Cross-referencing geophones and current list -----------------------------
setwd("C:/Users/murphyc4/OneDrive - Lancaster University/STOR-i/PhD/Projects/Spatio-temporal/Data/Geophones/Geophone data/")

Station_list <- read.csv("Station list.csv")
KNMI_list <- read.csv("KNMI Station list.csv", header=F)
#Geo_iris <- read.csv("Geophone_iris_KNMI.csv")
iris_all <- read.csv("query.csv")
iris_all$End_date <- stringr::str_sub(iris_all$EndTime,1,10)

#Geo_iris$End_date <- str_sub(Geo_iris$EndTime,1,10)
#Geo_ref <- Geo_iris[Geo_iris$End_date > "2022-12-13",]

iris_all_current <- iris_all[iris_all$End_date > "2022-12-13",]

# library(dplyr)
# Geo_ref$Channel3 <- str_sub(Geo_ref$Channel,3,3)
# E<-filter(Geo_ref, grepl("E", Channel3))
# N<-filter(Geo_ref, grepl("N", Channel3))
# Z<-filter(Geo_ref, grepl("Z", Channel3))

# Geo_ref <- rbind(E,N,Z)
# write.csv(Geo_ref, "Geo_ref.csv")


ind <- KNMI_list[,2] %in% unique(iris_all_current$Station)
(KNMI_missed <- KNMI_list[!ind,])
KNMI_missed[,2]

ind_2 <- Station_list$STATION %in% unique(iris_all_current$Station)
(Station_list_missed <- Station_list[!ind_2,])
Station_list_missed[,7]

ind_3 <- KNMI_missed[,2] %in% Station_list_missed[,7]
(KNMI_missed_from_Station <- KNMI_missed[!ind_3,])
KNMI_missed_from_Station[,2]


#Plotting station list stations which are not included
Station_list$RD_X <- as.numeric(gsub(",", "", Station_list$RD_X))
Station_list$RD_Y <- as.numeric(gsub(",", "", Station_list$RD_Y))

Gron_outline <- read.csv('Groningen_Field_outline.csv', header=T)

plot(Station_list[!ind_2,4], Station_list[!ind_2,5], ylab="Y", xlab="X")
points(Gron_outline$X, Gron_outline$Y, col="red", cex=0.5)
points(Station_list_missed$RD_X[19], Station_list_missed$RD_Y[19], cex=2.5, lwd=1.5, col="blue")
points(Station_list_missed$RD_X[20], Station_list_missed$RD_Y[20], cex=2.5, lwd=1.5, col="blue")
points(Station_list_missed$RD_X[64], Station_list_missed$RD_Y[64], cex=2.5, lwd=1.5, col="blue")

write.csv(KNMI_missed, "KNMI_missed.csv")
