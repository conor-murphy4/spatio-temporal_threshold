setwd("/home/murphyc4/Thr_Selection") #on STORM only (need normal working directory if running file in R)
source("STORM_input/_gpd.R")
source("STORM_input/Spatio-temporal threshold.R")

#----------Sim study--------------------------------------------

n_samples <- 200
n_theta <- 20
shape <- 0.1
scale <- 20
num_events <- 1960
# Est_scale_gpd <- Est_shape_gpd <- Est_theta_gpd <- matrix(NA, nrow = n_samples, ncol = n_theta)
Est_scale_exp <- Est_shape_exp <- Est_theta_exp <-  matrix(NA, nrow = n_samples, ncol = n_theta)
thresholds <- seq(0, 2.5, by=0.1)
set.seed(12345)
for(j in 1:n_samples){
  saveRDS(j, "STORM_output/iter_j.rds")
  #data generation
  num_events <- 1960
  M <- rgpd(num_events, shape = shape, scale=scale)
  X <- runif(num_events, 0, 40000)
  Y <- runif(num_events, 0, 60000)
  GeoX <- runif(60, 0, 40000)
  GeoY <- runif(60, 0, 60000)
  EQ.df <- data.frame(X=X, Y=Y, M=M)
  EQ.df$Year <- rep(seq(1995, 2022, by=1), times=1, each=70)
  Geo.df <- data.frame(X=GeoX, Y=GeoY)
  Geo.df$Start_Year <- c(rep(1995, 35), rep(1996,4), rep(2000, 5), rep(2005, 2), rep(2015, 14)) 
  Geo.df$End_Year <- c(rep(2005, 5 ), rep(2015, 3), rep(2022, 9), rep(2016, 4), rep(2022, 39) )
  V_sim <- V_2d_time(EQ.df, Geo.df)
  Sim_EQ.df <- data.frame(M=M,X=X, Y=Y, V=V_sim)
  theta <- 0
  for(i in 1:n_theta){ 
    saveRDS(i, "STORM_output/iter_i.rds")
    #Censoring
    cens_thr <- theta*V_sim
    keep <- M > cens_thr
    Sim_EQ_sub <- subset(Sim_EQ.df, keep)
    est_par_exp <- Spat_temp_thr_exp(Sim_EQ_sub, thresholds)
    # est_par_gpd <- Spat_temp_thr_gpd(Sim_EQ_sub, thresholds)
    Est_scale_exp[j,i] <- est_par_exp$par[1]
    Est_shape_exp[j,i] <- est_par_exp$par[2]
    Est_theta_exp[j,i] <- est_par_exp$thresh
    # Est_scale_gpd[j,i] <- est_par_gpd$par[1]
    # Est_shape_gpd[j,i] <- est_par_gpd$par[2]
    # Est_theta_gpd[j,i] <- est_par_gpd$thresh
    theta <- theta + 0.1
  }
}


saveRDS(Est_scale_exp, "STORM_output/Est_scale_exp.rds")
saveRDS(Est_shape_exp, "STORM_output/Est_shape_exp.rds")
saveRDS(Est_theta_exp, "STORM_output/Est_theta_exp.rds")

# saveRDS(Est_scale_gpd, "STORM_output/Est_scale_gpd.rds")
# saveRDS(Est_shape_gpd, "STORM_output/Est_shape_gpd.rds")
# saveRDS(Est_theta_gpd, "STORM_output/Est_theta_gpd.rds")

