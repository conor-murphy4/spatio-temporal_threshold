### Constructs vector and matrix indicating gridcells falling within polygon
print("This is 01_make_inPoly.R")
#---------------------------------------------------------
#install required packages
library(sp)

#---------------------------------------------------------
# Read Field Outline
Outline <- read.csv('Point processes/Covariates/00_data/raw/Field_outline.csv', header = TRUE)
Cov <- read.csv('Point processes/Covariates/00_data//raw/ReservoirModel_sm01_topograds.csv')

pol<-Polygon(coords = Outline[,1:2],hole = FALSE)
pols<- Polygons(list(pol),ID="1")
SpPols<- SpatialPolygons(list(pols), proj4string = CRS(as.character(NA)))
SPgfo<- SpatialPolygonsDataFrame(SpPols, data = data.frame(name = "GFO"),match.ID = TRUE)


polyx<- coordinates(pol)[,1]
polyy<- coordinates(pol)[,2]
boxx <- Cov$X
boxy <- Cov$Y

#---------------------------------------------------------
# Vector and matrix stating if covariate location is in gas field

inPolyVec <- point.in.polygon(point.x = boxx,point.y = boxy,
                              pol.x = polyx, pol.y = polyy)

inPolyMat <- matrix(inPolyVec,nrow = 84,ncol = 99,byrow = TRUE)

rm(Outline,pol,pols,SpPols,SPgfo,polyx, polyy, boxx,boxy, Cov)

#---------------------------------------------------------
# Save logical vector and matrix for use in analyses
saveRDS(object = inPolyVec,file =  "Point processes/Covariates/00_data/derived/inPoly/inPolyVec.RDS")
saveRDS(object = inPolyMat,file =  "Point processes/Covariates/00_data/derived/inPoly/inPolyMat.RDS")

## EOF ----
