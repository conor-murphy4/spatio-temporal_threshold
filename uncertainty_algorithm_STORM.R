
# Functions

transform_to_exp <- function (y, sig, xi){
  std_exp <- (1 / xi) * log( 1 + xi * (y/sig))  
  return(std_exp)
}

GPD_LL_given_V_ICS <- function(par, excess, thresh_par, V, ics){
  
  if(length(par)!=3) stop("par must be a vector of length 3")
  if(length(thresh_par)!=2) stop("thresh must be a vector of length 2")
  if (!is.numeric(excess)) stop("excess must be a vector")
  if (!is.numeric(V)) stop("V must be vector")
  if (!is.numeric(ics)) stop("ics must be vector")
  if(length(excess) != length(V)) stop("excess and V must be the same length")
  if(length(excess) != length(ics)) stop("excess and ics must be the same length")
  
  sigma_par <- par[1:2]
  xi <- par[3]
  
  sigma_ics <- sigma_par[1] + sigma_par[2]*ics
  
  sigma_tilde <- sigma_ics + xi*(thresh_par[[1]] + thresh_par[[2]]*V)
  
  sigma_check <- c(sigma_ics, sigma_tilde)
  
  if(all(sigma_check > 0) & xi > -0.9){
    if(abs(xi) < 1e-10){
      return(-sum(log(sigma_tilde)) - sum(excess/sigma_tilde))
    }
    else{
      if(all(1+(xi*excess)/sigma_tilde > 0)){
        return(-sum(log(sigma_tilde)) - (1+1/xi)*(sum(log(1+(xi*excess)/sigma_tilde))))
      }
      else{
        return(-1e6)
      }
    }
  }
  else{
    return(-1e7)
  }
}


eqd_geo_ics <- function(data, thresh, distance_to_geo, ics, k = 100, m = 500, underlying_thresh = 0){
  
  # Check inputs are valid
  if (!is.numeric(data)) stop("Data must be a vector")
  if (!is.matrix(thresh)) stop("Thresholds must be a matrix")
  if (!is.numeric(distance_to_geo)) stop("Third nearest distances must be a vector")
  if (length(distance_to_geo) != length(data)) stop("Length of third nearest distances must match length of data")
  if (k <= 0 | k %% 1 != 0) stop("Number of bootstrapped samples must be a positive integer")
  if (m <= 0 | m %% 1 != 0) stop("Number of equally spaced probabilities must be a positive integer")
  
  meandistances <- xis <- num_excess <- numeric(nrow(thresh))
  sigma_par <- matrix(NA, nrow = nrow(thresh), ncol = 2)
  ########################## extra code for Inf investigation ##################
  # inf_list_mean <- list()
  # inf_list_boot <- list()
  ########################## extra code for Inf investigation ##################
  for (i in 1:nrow(thresh)) {
    #print(i)
    u <- thresh[i,1] + thresh[i,2] * distance_to_geo
    #if(any(u < underlying_thresh)) stop(paste("Candidate thresholds must be above the underlying threshold of ", underlying_thresh))
    if(any(u < underlying_thresh)) {
      meandistances[i] <- NA
      next
      }
    excess <- data[data > u] - u[data > u]
    V_excess <- distance_to_geo[data > u]
    ICS_excess <- ics[data > u]
    num_excess[i] <- length(excess)
    if (num_excess[i] > 20) {
      mle0 <- mean(excess)
      init.fit <- optim(GPD_LL_given_V_ICS, excess = excess, par = c(mle0, 0, 0.1), control = list(fnscale = -1), 
                        thresh_par=thresh[i,], V = V_excess, ics=ICS_excess)
      xis[i] <- init.fit$par[3]
      sigma_par[i,] <- init.fit$par[1:2]
      distances <- numeric(k)
      j <- 1
      while(j  <= k) {
        #possible way to bootstrap
        sampled_indices <- sample(1:num_excess[i], num_excess[i], replace = TRUE)
        excess_boot <- excess[sampled_indices]
        V_excess_boot <- V_excess[sampled_indices]
        ICS_excess_boot <- ICS_excess[sampled_indices]
        init_mle <- mean(excess_boot)
        ifelse(xis[i] < 0, pars_init <-  c(init_mle, 0, 0.1) ,pars_init <- c(sigma_par[i,], xis[i]) )
        gpd.fit <- optim(GPD_LL_given_V_ICS, excess = excess_boot, par = pars_init, control = list(fnscale = -1), 
                         thresh_par=thresh[i,], V = V_excess_boot, ics = ICS_excess_boot)
        sigma_given_V <- gpd.fit$par[1] + gpd.fit$par[2]*ICS_excess_boot + gpd.fit$par[3] * (thresh[i,1] + thresh[i,2]*V_excess_boot)
        if(any((1 + gpd.fit$par[3] * excess_boot/sigma_given_V) < 1e-323)) next
        else{
          transformed_excess_boot <- transform_to_exp(y = excess_boot, sig = sigma_given_V, xi = gpd.fit$par[3])
          sample_quants <- quantile(transformed_excess_boot, probs = (1:m) / (m+1))
          model_quants <- qexp((1:m) / (m+1), rate = 1)
          distances[j] <- mean(abs(sample_quants - model_quants))
          ########################## extra code for Inf investigation ##########
          # if(distances[j] == Inf){
          #   inf_list_current <- list(i = i, j = j, excess_boot = excess_boot, transformed_excess = transformed_excess_boot, est = gpd.fit$par, sigma_var = sigma_given_V)
          #   inf_list_boot <- c(inf_list_boot, list(inf_list_current))
          # }
          ########################## extra code for Inf investigation ##########
          j <- j + 1
        }
      }
      meandistances[i] <- mean(distances)
      ########################## extra code for Inf investigation ##############
      # if(meandistances[i] == Inf){
      #   inf_list_mean <- c(inf_list_mean, list(list(i = i, excess = excess, par_ests = init.fit$par, dists = distances)))
      # }
      ########################## extra code for Inf investigation ##############
    }
    else{
      meandistances[i] <- NA
    }
  }
  chosen_index <- which.min(meandistances)
  chosen_threshold_par <- thresh[chosen_index,]
  xi <- xis[chosen_index]
  sigmas <- sigma_par[chosen_index,]
  len <- num_excess[chosen_index]
  result <- list(thresh_par = chosen_threshold_par, par = c(sigmas,xi), num_excess = len, dists = meandistances)
  ########################## extra code for Inf investigation ##################
  # result <- list(thresh_par = chosen_threshold_par, par = c(sigmas,xi), num_excess = len, dists = meandistances, inf_list_mean = inf_list_mean, inf_list_boot = inf_list_boot)
  ########################## extra code for Inf investigation ##################
  return(result)
}

threshold_selection_varying_formulation <- function(data, threshold_matrix, num_boot=200, SEED=11111, output_all = FALSE) {
  message(Sys.time(),"Starting selection")
  nearest_dist_matrix <- matrix(c(data$V_1, data$V_2, data$V_3, data$V_4), ncol=4, byrow=F)
  results <- vector("list", 12)
  for (ii in c(1:4)) {
    nearest_dist <- nearest_dist_matrix[,ii]
    log_nearest_dist <- log(nearest_dist)
    sqrt_nearest_dist <- sqrt(nearest_dist)
    
    # V
    set.seed(SEED)
    results[[ii]] <- eqd_geo_ics(data=data$Magnitude, thresh = threshold_matrix, distance_to_geo = nearest_dist, k=num_boot, ics = data$ICS_max)
    
    # log(V)
    set.seed(SEED)
    results[[ii+4]] <- eqd_geo_ics(data=data$Magnitude, thresh = threshold_matrix, distance_to_geo =  log_nearest_dist, k=num_boot, ics = data$ICS_max)
    
    # sqrt(V)
    set.seed(SEED)
    results[[ii+8]] <- eqd_geo_ics(data=data$Magnitude, thresh = threshold_matrix, distance_to_geo = sqrt_nearest_dist, k=num_boot, ics = data$ICS_max)
  }
  chosen_model_index <- which.min(lapply(results, function(x){
    if (is.null(x) || is.null(x$dists) || all(is.na(x$dists))) return(Inf)
    min(x$dists, na.rm = T)
    } ))
  chosen_model <- results[[chosen_model_index]]
  if (output_all) {
    return(list(chosen_form = chosen_model_index, model_results=results))
  }
  else {
    return(list(chosen_form = chosen_model_index, model_results=chosen_model))
  }
}

non_parametric_samples <- readRDS("STORM_input/non_parametric_samples.rds")


intercepts <- seq(0, 1.5, by=0.02)
slopes <- seq(0,0.8, by=0.02)
threshold_matrix <- as.matrix(expand.grid(intercepts, slopes))

library(parallel)

cl <- makeCluster(100)
clusterExport(cl, varlist = c("threshold_selection_varying_formulation", "threshold_matrix", "eqd_geo_ics", "GPD_LL_given_V_ICS", "transform_to_exp"))
clusterSetRNGStream(cl, 11111)

bootstrap_results <- parLapply(cl, non_parametric_samples, function(sample) {
  threshold_selection_varying_formulation(
    data = sample,
    threshold_matrix = threshold_matrix,
    num_boot = 100,
    SEED = 11111,
    output_all = FALSE
  )
})
stopCluster(cl)

saveRDS(bootstrap_results, file="STORM_output/model_uncertainty_results.rds")
