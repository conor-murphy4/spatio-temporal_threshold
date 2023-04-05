Opt  <- readRDS("Point processes/Covariates/Output/data/Opt_S4.RDS")

B0 <- Opt$par[1]
B1 <- Opt$par[2]
sigma <- Opt$par[3]


intensity <- B0*
