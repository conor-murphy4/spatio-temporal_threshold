## This script makes unsmoothed covariate arrays for  earthquake location models

## Script parameters -----------------
SAVE_COVS <- TRUE
DATA_PATH <- "00_data/derived/covariates/"

MAKE_PLOTS <- TRUE
SAVE_PLOTS <- TRUE
PLOT_PATH <- "99_output/00_data/covariates/"

REPLACE_NEGATIVES <- TRUE
REPLACEMENT_VALUE <- 1e-9

## Load required packages ------------
#library(lubridate)
library(sp)     # for field outline
#library(MASS)
library(fields) # Only required for plotting
library(abind)  # for temporal differencing of covariates
#library(mgcv)
#library(purrr)


## Source required scripts --------------
#source("../31_topographic_gradient/topograd_opt_smooth_functions.R")
source("making_covarite_arrays/covariate_array_functions.R")
## Read in relevant data ----------------

# Groningen field outline
SPgfo<- readRDS("00_data/derived/field_outline/SPgfo.RDS")
inPolyVec <- readRDS("00_data/derived/inPoly/inPolyVec.RDS")
inPolyMat <- readRDS("00_data/derived/inPoly/inPolyMat.RDS")

# topographic gradient ---------------------------------------------------------
tg <- read.csv('00_data/raw/ReservoirModel_sm01_topograds.csv',header = TRUE)
tg_mats <- make_cov_mats(
  cov_df = tg,
  df_t0 = 1957,
  array_tlim = c(1957, 2020),
  inpolymat = inPolyMat,
  cov_name = 'tg',
  is_constant = TRUE)

tg_mats<- replace_negatives_cov_mats(
  cov_mats = tg_mats,
  replace = REPLACE_NEGATIVES,
  replacement_value = REPLACEMENT_VALUE)

if(SAVE_COVS){
  rds_path <- paste0(DATA_PATH, "tg_mats.RDS")
  saveRDS(object = tg_mats, file = rds_path)
}

if(MAKE_PLOTS){
  if(SAVE_PLOTS){
    pdf_path <- paste0(PLOT_PATH, "topograd_annual.pdf")
    pdf(pdf_path, width = 7,height = 5)
  }
  plot_cov_mats(tg_mats, cov_name = "togographic gradient")
  if(SAVE_PLOTS){ dev.off()}

}

# Compressibility --------------------------------------------------------------
csv_path <- "C:/Users/murphyc4/OneDrive - Lancaster University/STOR-i/PhD/Projects/zv-cm-code/point-process-simulation/seis-mod-dev-v6/seis-mod-dev-v6/inputs/reservoirs/ReservoirModel_compressibility_20180122.csv"
##----------STOPPED------------------
compress_df <- read.csv(csv_path)

compressibility_mats <- make_cov_mats(
  cov_df = compress_df,
  df_t0 = 1957,
  array_tlim = c(1957, 2020),
  inpolymat = inPolyMat,
  cov_name = 'compressibility',
  is_constant = TRUE)

compressibility_mats <- replace_negatives_cov_mats(
  cov_mats = compressibility_mats,
  replace = REPLACE_NEGATIVES,
  replacement_value = REPLACEMENT_VALUE)

quantile(compressibility_mats[[1]], probs = (0:10)/10, na.rm = TRUE)

if(SAVE_COVS){
  rds_path <- paste0(DATA_PATH, "compressibility_mats.RDS")
  saveRDS(object = tg_mats, file = rds_path)
}

if(MAKE_PLOTS){
  if(SAVE_PLOTS){
    pdf_path <- paste0(PLOT_PATH, "compressibility_annual.pdf")
    pdf(pdf_path, width = 7,height = 5)
  }
  plot_cov_mats(compressibility_mats, cov_name = "Compressibility")
  if(SAVE_PLOTS){ dev.off() }
}

# Poroelastic constant ---------------------------------------------------------

H_s <- 1e10
H <- 1 / (compressibility_mats[[1]] + H_s^(-1))

H_mats <- compressibility_mats
H_mats[[1]] <- H
names(H_mats)[1] <- "H_mats"

# Check range of values
quantile(H_mats[[1]], probs = (0:10)/10, na.rm = TRUE)

H_mats <- replace_negatives_cov_mats(
  cov_mats = H_mats,
  replace = REPLACE_NEGATIVES,
  replacement_value = REPLACEMENT_VALUE)

if(SAVE_COVS){
  rds_path <- paste0(DATA_PATH, "poroelestic_modulus_mats.RDS")
  saveRDS(object = tg_mats, file = rds_path)
}

if(MAKE_PLOTS){
  if(SAVE_PLOTS){
    pdf_path <- paste0(PLOT_PATH, "poroelastic_modulus_annual.pdf")
    pdf(pdf_path, width = 7,height = 5)
  }
  plot_cov_mats(H_mats, cov_name = "Poroelastic modulus", zlim = c(0, 3.5e4))
  if(SAVE_PLOTS){ dev.off() }
}


## Compaction (Shell) ---------------------------------------------------------
csv_path <- "00_data/raw/ReservoirModel_rm02_s01_24bcm_compaction_linear.csv"
compaction_df <- read.csv(csv_path)

compaction_mats <- make_cov_mats(
  cov_df = compaction_df,
  df_t0 = 1959,
  array_tlim = c(1959,2021),
  inpolymat = inPolyMat,
  cov_name = "compaction",
  is_constant = FALSE)

quantile(compaction_mats[[1]], probs = (0:10)/10, na.rm = TRUE)

compaction_mats <- replace_negatives_cov_mats(
  cov_mats = compaction_mats,
  replace = REPLACE_NEGATIVES,
  replacement_value = REPLACEMENT_VALUE)

if(SAVE_COVS){
  rds_path <- paste0(DATA_PATH, "compaction_shell_mats.RDS")
  saveRDS(object = compaction_mats, file = rds_path)
}

if(MAKE_PLOTS){
  if(SAVE_PLOTS){
    pdf_path <- paste0(PLOT_PATH, "compaction_shell_annual.pdf")
    pdf(pdf_path, width = 7,height = 5)
  }
  plot_cov_mats(compaction_mats, cov_name = "Compaction (m)", zlim = c(0, 0.36))
  if(SAVE_PLOTS){ dev.off() }
}

# Pressure depletion -----------------------------------------------------------

pressure <- read.csv("C:/Users/murphyc4/OneDrive - Lancaster University/STOR-i/PhD/Projects/zv-cm-code/point-process-simulation/seis-mod-dev-v6/seis-mod-dev-v6/inputs/reservoirs/reservoir_gf_2021_pressure.csv")

# note: pressure is measured at days rather than over years so there is an extra
# column as compared to compaction. To adjust for this calculate pressure
# depletion over each year and remove first covariate column (3).

depl <- pressure
for (col in NCOL(depl):3){
  depl[,col] <- depl[,3] - depl[,col]
}
names(depl) <- c("X","Y", paste0("dP",as.character(1957:2020),"_01_01"))

# Note 2: Pressure (and therefore depletion) have an extra X and Y value as
# compared to topograd and compressibility. These are the largest X and Y values
# respectively. Remove observations where X = 270,000 or Y = 566,500.
depl_subsetted <- depl[!((depl$X == 270000) | (depl$Y == 566000)), ]


depl_mats <- make_cov_mats(
  cov_df = depl_subsetted,
  df_t0 = 1957,
  array_tlim = c(1957, 2020),
  inpolymat = inPolyMat,
  cov_name = "depletion",
  is_constant = FALSE)

# Check range of depletion values
quantile(depl_mats[[1]], probs = (0:10)/10, na.rm = TRUE)

depl_mats <- replace_negatives_cov_mats(
  cov_mats = depl_mats,
  replace = REPLACE_NEGATIVES,
  replacement_value = REPLACEMENT_VALUE)

if(SAVE_COVS){
  rds_path <- paste0(DATA_PATH, "pressure_depletion_mats.RDS")
  saveRDS(object = depl_mats, file = rds_path)
}

if(MAKE_PLOTS){
  if(SAVE_PLOTS){
    pdf_path <- paste0(PLOT_PATH, "pressure_depletion_annual.pdf")
    pdf(pdf_path, width = 7,height = 5)
  }
  plot_cov_mats(depl_mats, zlim = c(-1, 310), cov_name = "depletion (MPa)")
  if(SAVE_PLOTS){ dev.off() }
}

# compaction (unsmoothed) ------------------------------------------------------

compaction <- multiply_cov_mats(
  cov_mats_1 = depl_mats,
  cov_mats_2 = compressibility_mats,
  new_name = "compaction")

compaction <- replace_negatives_cov_mats(
  cov_mats = compaction,
  replace = REPLACE_NEGATIVES,
  replacement_value = REPLACEMENT_VALUE)


if(SAVE_COVS){
  rds_path <- paste0(DATA_PATH, "compaction_mats.RDS")
  saveRDS(object = compaction, file = rds_path)
}

if(MAKE_PLOTS){
  if(SAVE_PLOTS){
    pdf_path <- paste0(PLOT_PATH, "compaction_annual.pdf")
    pdf(pdf_path, width = 7,height = 5)
  }
  plot_cov_mats(compaction, zlim = c(0, 0.03), cov_name = "Compaction (m)")
  if(SAVE_PLOTS){ dev.off() }
}

# strain thickness -------------------------------------------------------------
strain_thickness <- multiply_cov_mats(
  cov_mats_1 = compaction,
  cov_mats_2 =  tg_mats,
  new_name = "strain_thickness")

quantile(strain_thickness[[1]], probs = (0:10)/10, na.rm = TRUE)

strain_thickness <- replace_negatives_cov_mats(
  cov_mats = strain_thickness,
  replace = REPLACE_NEGATIVES,
  replacement_value = REPLACEMENT_VALUE)

if(SAVE_COVS){
  rds_path <- paste0(DATA_PATH, "strain_thickness_mats.RDS")
  saveRDS(object = strain_thickness, file = rds_path)
}

if(MAKE_PLOTS){
  if(SAVE_PLOTS){
    pdf_path <- paste0(PLOT_PATH, "strain_thickness_annual.pdf")
    pdf(pdf_path, width = 7,height = 5)
  }
  plot_cov_mats(strain_thickness, zlim = c(0, 7.2e-4), cov_name = "strain thickness")
  if(SAVE_PLOTS){ dev.off() }
}

# incremental Coulomb stress ---------------------------------------------------
ics <- multiply_cov_mats(strain_thickness, H_mats, new_name = "ics")

quantile(ics[[1]], probs = (0:10)/10, na.rm = TRUE)

if(SAVE_COVS){
  rds_path <- paste0(DATA_PATH, "incremental_Coulomb_stress_mats.RDS")
  saveRDS(object = ics, file = rds_path)
}

if(MAKE_PLOTS){
  if(SAVE_PLOTS){
    pdf_path <- paste0(PLOT_PATH, "incremental_Coulomb_stress_annual.pdf")
    pdf(pdf_path, width = 7,height = 5)
  }
  plot_cov_mats(ics, cov_name = "ICS")
  if(SAVE_PLOTS){ dev.off() }
}


## CHECKS ON NEGATIVE VALUES - DOES NOT HAVE TO BE RUN ----------------------------------------------

# # Plot indicator of non-postive covariate values
# is_negative_ics <- ics
# is_negative_ics[[1]] <- is_negative_ics[[1]] <= 0
# plot_cov_mats(is_negative_ics)
#       # all before start of modelling period apart from a few tricky pixels
#
# is_negative_st <- st_zv
# is_negative_st[[1]] <- is_negative_st[[1]] <= 0
# plot_cov_mats(is_negative_st)
#
#
# # Investigate non-positive covariate values in each covariate over time
# non_positive_count <- function(x){sum(x <= 0, na.rm = TRUE)}
#
# npc <- data.frame(year = compaction_zv$tgrid)
# npc$depletion <- apply(depl_mats[[1]], MARGIN = 3, FUN = non_positive_count)
# npc$compaction <- apply(compaction_zv[[1]], MARGIN = 3, FUN = non_positive_count)
# npc$compressibility <- apply(compressibility_mats[[1]], MARGIN = 3, FUN = non_positive_count)
# npc$topograd <- apply(tg_mats[[1]], MARGIN = 3, FUN = non_positive_count)
# npc$H <- apply(H_mats[[1]], MARGIN = 3, FUN = non_positive_count)
# npc$st <- apply(st_zv[[1]], MARGIN = 3, FUN = non_positive_count)
# npc$ics <- apply(ics_zv[[1]], MARGIN = 3, FUN = non_positive_count)
#
# # Plot number of non-positive covariate values in each time slice
# plot(x = npc$year, y = npc$depletion, data = npc, type = 'l', ylim = c(0,20)) # none throughout study period, more before
# lines(x = npc$year, y = npc$compaction,  col = 2) # 2 throughout study period, higer before
# lines(x = npc$year, y = npc$compressibility,  col = 3) # 4 constant
# lines(x = npc$year, y = npc$topograd,  col = 4) # none
# lines(x = npc$year, y = npc$H,  col = 5) # none
# lines(x = npc$year,  y = npc$st, col = "orange") # as compaction: 2 in study period
# lines(x = npc$year, y = npc$ics, col = 'hotpink') # as compaciton

## Conclusion: setting negative values to small positive ones will be okay.

