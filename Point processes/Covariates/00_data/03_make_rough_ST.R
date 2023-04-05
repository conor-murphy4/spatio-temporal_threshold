# This script will make an unsomoothed version of the strain thickness covariate

## 0.1: Set script arguments and load packages ---------------------------------

# Script arguments
MAKE_PLOTS = FALSE
SAVE_PLOTS = FALSE
PLOT_PATH  = "99_output/99_covariate_plots/ST_annual_rough_and_smooth.pdf"

CSV_PATH = "00_data/derived/covariates/strainThicknessRough_1.csv"

# Load packages
library(fields)
library(MASS)
# library(lubridate)    # soft requirement (plots only, lubridate::year() )

## 0.2: Read data --------------------------------------------------------------

## Read strain thickness data
csv_path <- "00_data/raw/ReservoirModel_rm03_s01_24bcm_strainthickness_linear.csv"
st <- read.csv(csv_path,header = TRUE)

## Read compaction data
csv_path <- "00_data/raw/ReservoirModel_rm02_s01_24bcm_compaction_linear.csv"
compaction <- read.csv(csv_path,header = TRUE)

## Read topographic gradient data
csv_path <- "00_data/raw/ReservoirModel_sm01_topograds.csv"
topo_grads <- read.csv(csv_path,header = TRUE)

if (MAKE_PLOTS){
## Read field outline
rds_path <- "00_data/derived/field_outline/SPgfo.RDS"
SPgfo <- readRDS(rds_path)

## Read event catalogue
csv_path <- "00_data/raw/EventCatalogue_20170428.csv"
EventCatalogue <- read.csv(csv_path)
}

## 1.1: Functions to make covariate arrays   -----------------------------------


# Function returning list of:
# array of gridded (topological graident * compaction) each year
# x grid values
# y grid values

make.st_mats <- function(Y1_lim = c(1959,1959)){

  nx   <- length(unique(st$X))
  minx <- min(st$X)
  maxx <- max(st$X)
  ny   <- length(unique(st$Y))
  miny <- min(st$Y)
  maxy <- max(st$Y)

  xgrid = seq(minx,maxx,length= nx ) # x values of observation points
  ygrid = seq(miny,maxy,length= ny) # y values of observation points

  year <- 1:64
  st_mats <- array(data = NA, dim = c(nx,ny,length(year)))

  for (j in year){
    sts <- st[[(Y1_lim[2] - 1959 + 1 + j + 2)]] - st[[(Y1_lim[1] - 1959 + j + 2)]]

    st_mat <- matrix(NA, nrow = nx, ncol = ny)
    for( i in 1: length(st$X)){
      x <- (st$X[i] - minx)/500 + 1
      y <- (st$Y[i] - miny)/500 + 1
      st_mat[x,y] <- sts[i]
    }

    st_mats[,,j] <- st_mat
  }

  out <- list(st_mats,xgrid,ygrid)
  names(out) <- c('st_mats','xgrid','ygrid')
  return(out)
}


# Function returning list of:
# array of gridded (topological graident * compaction) each year
# x grid values
# y grid values

make.gradcomp_mats <- function(Y1Comp_lim = c(1959,1959)){

  nx   <- length(unique(compaction$X))
  minx <- min(compaction$X)
  maxx <- max(compaction$X)
  ny   <- length(unique(compaction$Y))
  miny <- min(compaction$Y)
  maxy <- max(compaction$Y)

  xgrid = seq(minx,maxx,length= nx ) # x values of observation points
  ygrid = seq(miny,maxy,length= ny) # y values of observation points

  topo_mat <- matrix(NA, nrow = nx, ncol = ny)
  for( i in 1: length(topo_grads$X)){
    x <- (topo_grads$X[i] - minx)/500 + 1
    y <- (topo_grads$Y[i] - miny)/500 + 1
    topo_mat[x,y] <- topo_grads$topo_gradients_smoothed[i]
  }

  year <- 1:63
  gradcomp_mats <- array(data = NA, dim = c(nx,ny,length(year)))

  for (j in year){
    compactions <- compaction[[(Y1Comp_lim[2] - 1959 + 1 + j +2 -1)]]# - Compaction[[(Y1Comp_lim[1] - 1959 + j +2)]]

    comp_mat <- matrix(NA, nrow = nx, ncol = ny)
    for( i in 1: length(compaction$X)){
      x <- (compaction$X[i] - minx)/500 + 1
      y <- (compaction$Y[i] - miny)/500 + 1
      comp_mat[x,y] <- compactions[i]
    }

    gradcomp_mats[,,j] <- comp_mat * topo_mat * 100
  }

  out <- list(gradcomp_mats, xgrid, ygrid)
  names(out) <- c('gradcomp_mats','xgrid','ygrid')
  return(out)
}

## 2.0: Make covariate objects -------------------------------------------------
X <- make.st_mats(Y1_lim = c(1959,1959))
Y <- make.gradcomp_mats(Y1Comp_lim = c(1959,1959))

## 3.0: (optional) plot annual covariates, overlay events and outline ----------

if(MAKE_PLOTS){

  if(SAVE_PLOTS){
    pdf_path <- PLOT_PATH
    pdf(file = pdf_path, width = 12, height = 5)
  }

  par(mfrow = c(1,2),omi = c(0.5,0.5,0.5,0.5))

  # Plot compaction*topograd and ST maps for each year.
  for (yr in 1:63){
    image.plot(x = X$xgrid,
               y = X$ygrid,
               z = X$st_mats[,,yr],
               asp=1,
               xlab = 'Easting',
               ylab = 'Northing',
               #zlim = c(0,0.013),
               main= paste('Strain Thickness',
                           1959+yr,
                           '-',
                           1959+yr,
                           '   ',
                           'Events in',
                           1959+yr)
               )

    date_condition <- lubridate::year(EventCatalogue$Date) == yr +1959
    mag_condition  <- EventCatalogue$Magnitude >= 1.5
    events_E <- EventCatalogue$Easting[  date_condition & mag_condition]
    events_N <- EventCatalogue$Northing[ date_condition & mag_condition]
    events_M <- EventCatalogue$Magnitude[date_condition & mag_condition]

    points(x = events_E,y = events_N, cex = sqrt(events_M/pi))
    plot(SPgfo,add= TRUE)

    image.plot(x = Y$xgrid,
               y = Y$ygrid,
               z = Y$gradcomp_mats[,,yr],
               asp=1,
               xlab = 'Easting',
               ylab = 'Northing',
               #zlim = c(0,0.013),
               main= paste('Comp * Topo',
                           1959 + yr,
                           '-',
                           1959 + yr,
                           '   ',
                           'Events in',
                           1959 + yr)
               )

    date_condition <- lubridate::year(EventCatalogue$Date) == yr + 1959
    mag_condition  <- EventCatalogue$Magnitude >= 1.5
    events_E <- EventCatalogue$Easting[  date_condition & mag_condition]
    events_N <- EventCatalogue$Northing[ date_condition & mag_condition]
    events_M <- EventCatalogue$Magnitude[date_condition & mag_condition]

    points(x = events_E,y = events_N, cex = sqrt(events_M/pi))
    plot(SPgfo,add= TRUE)
  }
}

if(SAVE_PLOTS){dev.off()}

## 4.0: Function to convert covariate array to savable data frame --------------

arrayToTable <- function(A, YrInit){
  #function to convert a covariate array to
  #a table of the form covarites are supplied in.

  # A = covariate array
  # Xgrid = x locations at which recorded
  # Ygrid = y locations at which recorded
  # YrRange = c(minYear, maxYear) for covariate

  nx <-  dim(A[[1]])[1]
  ny <-  dim(A[[1]])[2]
  nYr <- dim(A[[1]])[3]

  out<- cbind( rep(A[[2]],each = ny) , rep(A[[3]], nx))

  for(k in 1: nYr){
    values <- c()
    for (i in 1: nx){
      values <-  c(values,A[[1]][i,,k])
    }

    out <- cbind(out,values)
  }

  out <- data.frame(out)
  names(out)[1:2] <- c('X','Y')
  for(i in 3:(nYr+2)){
    names(out)[i] <- paste('dV', YrInit + i -3,'_01_01', sep = '')
  }

  return(out)
}

## 4.0: Convert and save rough strain thickness --------------------------------
STrough<- arrayToTable(Y,YrInit = 1959)
write.csv(x = STrough,file = CSV_PATH,row.names = FALSE)

## EOF
