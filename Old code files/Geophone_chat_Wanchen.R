#Full geophone set after some processing (removing irrelevant geophones, inlcuding ones with E,N,Z, etc.)
Geo <- read.csv('Data/Geophones/Geophone data/Geo_orig_post_process.csv', header=T)

#Outline of gasfield
Gron_outline <- read.csv('STOR-i/PhD/Projects/Spatio-temporal/Data/Geophones/Geophone data/Groningen_Field_outline.csv', header=T)

#Rectangle which includes all relevant geophones
X <- c(227000,227000,272000,272000,227000)
Y <- c(535000,614000,614000,535000,535000)
Gron_rect <- data.frame(X=X, Y=Y)

#Geophones actually used in our analysis based on V(x,t)
Geo_used <- read.csv('Data/Geophones/Geophone data/Geo_final.csv', header=T)

#Across Netherlands
plot(Geo$Xcoord, Geo$Ycoord, asp=1, xlab = "Easting(m)", ylab = "Northing(m)")
points(Gron_outline$X, Gron_outline$Y, col="blue", cex=0.5)
points(Geo_used$X, Geo_used$Y, pch=19)
lines(Gron_rect$X, Gron_rect$Y, col="red")

#Zoomed into rectangle
plot(Gron_rect$X, Gron_rect$Y, col="red", asp=1, type='l', xlab = "Easting(m)", ylab = "Northing(m)")
points(Gron_outline$X, Gron_outline$Y, col="blue", cex=0.5)
points(Geo$Xcoord, Geo$Ycoord)
points(Geo_used$X, Geo_used$Y, pch=19)

#Zoomed with just relevant geophones included
plot(Gron_rect$X, Gron_rect$Y, col="red", asp=1, type='l', xlab = "Easting(m)", ylab = "Northing(m)")
points(Gron_outline$X, Gron_outline$Y, col="blue", cex=0.5)
points(Geo_used$X, Geo_used$Y, pch=19)

#Time-development of relevant geophone network 

#Geophones on 01-01-2013
ind_2013_01_01 <- which(Geo_used$Start_date <= "2013-01-01" & Geo_used$End_date > "2013-01-01")
Geo_2013_01_01 <- Geo_used[ind_2013_01_01,]
#Added Geophones in 2013
ind_add_2013 <- which(Geo_used$Start_date > "2013-01-01" & Geo_used$Start_date <= "2014-01-01")
Geo_add_2013 <- Geo_used[ind_add_2013,]

Geo_current <- rbind(Geo_2013_01_01, Geo_add_2013)

#Removed Geophones in 2013
ind_rem_2013 <- which(Geo_current$End_date >= "2013-01-01" & Geo_current$End_date < "2014-01-01")
Geo_rem_2013 <- Geo_2013_01_01[ind_rem_2013,]
#####No geophones removed in 2013

#Added Geophones in 2014
ind_add_2014 <- which(Geo_used$Start_date > "2014-01-01" & Geo_used$Start_date <= "2015-01-01" & Geo_used$End_date > "2015-01-01")
Geo_add_2014 <- Geo_used[ind_add_2014,]

Geo_current <- rbind(Geo_current, Geo_add_2014)

#Removed Geophones in 2014
ind_rem_2014 <- which(Geo_current$End_date >= "2014-01-01" & Geo_current$End_date < "2015-01-01")
Geo_rem_2014 <- Geo_current[ind_rem_2013,]
#####No geophones removed in 2014

#Geophones on 01-01-2015
ind_2015_01_01 <- which(Geo_used$Start_date <= "2015-01-01" & Geo_used$End_date > "2015-01-01")
Geo_2015_01_01 <- Geo_used[ind_2015_01_01,]

#Plotting of time-development of geophones from 2013-2015
plot(Gron_rect$X, Gron_rect$Y, col="red", asp=1, type='l', xlab = "Easting(m)", ylab = "Northing(m)")
points(Gron_outline$X, Gron_outline$Y, col="blue", cex=0.5)
points(Geo_2013_01_01$X, Geo_2013_01_01$Y, pch=19)
points(Geo_add_2013$X, Geo_add_2013$Y, pch=19, col="red")
points(Geo_add_2014$X, Geo_add_2014$Y, pch=19, col="green")

#Time-development of full geophone set

#Geophones on 01-01-2013
ind_2013_01_01 <- which(Geo$Start_date <= "2013-01-01" & Geo$End_date > "2013-01-01")
Geo_2013_01_01 <- Geo[ind_2013_01_01,]
#Added Geophones in 2013
ind_add_2013 <- which(Geo$Start_date > "2013-01-01" & Geo$Start_date <= "2014-01-01")
Geo_add_2013 <- Geo[ind_add_2013,]

Geo_current <- rbind(Geo_2013_01_01, Geo_add_2013)

#Removed Geophones in 2013
ind_rem_2013 <- which(Geo_current$End_date >= "2013-01-01" & Geo_current$End_date < "2014-01-01")
Geo_rem_2013 <- Geo_2013_01_01[ind_rem_2013,]
#####No geophones removed in 2013

#Added Geophones in 2014
ind_add_2014 <- which(Geo$Start_date > "2014-01-01" & Geo$Start_date <= "2015-01-01" & Geo$End_date > "2015-01-01")
Geo_add_2014 <- Geo[ind_add_2014,]

Geo_current <- rbind(Geo_current, Geo_add_2014)

#Removed Geophones in 2014
ind_rem_2014 <- which(Geo_current$End_date >= "2014-01-01" & Geo_current$End_date < "2015-01-01")
Geo_rem_2014 <- Geo_current[ind_rem_2014,]
#####No geophones removed in 2014

#Geophones on 01-01-2015
ind_2015_01_01 <- which(Geo$Start_date <= "2015-01-01" & Geo$End_date > "2015-01-01")
Geo_2015_01_01 <- Geo[ind_2015_01_01,]

#Plotting of time-development of geophones from 2013-2015
plot(Gron_rect$X, Gron_rect$Y, col="red", asp=1, type='l', xlab = "Easting(m)", ylab = "Northing(m)")
points(Gron_outline$X, Gron_outline$Y, col="blue", cex=0.5)
points(Geo_2013_01_01$Xcoord, Geo_2013_01_01$Ycoord, pch=19)
points(Geo_add_2013$Xcoord, Geo_add_2013$Ycoord, pch=19, col="red")
points(Geo_add_2014$Xcoord, Geo_add_2014$Ycoord, pch=19, col="green")
points(Geo_rem_2013$Xcoord, Geo_rem_2013$Ycoord, pch=4, cex=1.5, col="red")
points(Geo_rem_2014$Xcoord, Geo_rem_2014$Ycoord, pch=4, cex=1.5, col="red")

