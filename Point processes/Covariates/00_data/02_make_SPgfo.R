## This file will make a spatial polygons data frame for the Groniningen field.

print("This is 02_make_SPgfo.R")

Outline <- read.csv('../data/Field_outline.csv', header = TRUE)
pol<-Polygon(coords = Outline[,1:2],hole = FALSE)
pols<- Polygons(list(pol),ID="1")
SpPols<- SpatialPolygons(list(pols), proj4string = CRS(as.character(NA)))
SPgfo<- SpatialPolygonsDataFrame(SpPols, data = data.frame(name = "GFO"),match.ID = TRUE)
rm(pol,pols,SpPols, Outline)

saveRDS(SPgfo, file = "./derived/field_outline/SPgfo.RDS")
rm(SPgfo)

## EOF
