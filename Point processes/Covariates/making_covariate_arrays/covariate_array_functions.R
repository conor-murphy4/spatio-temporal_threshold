#' Make covariate array from covariate data frame
#'
#' @param cov_df data frame of covariate values. Columns 1 and 2 give X and Y coordinates, subsequent columns give covarite value each year. For static covariates, this has only three columns.
#' @param df_t0 Repeated observations assumed to be on Jan 1st each year. This gives the year of first observation.
#' @param array_tlim integer vecot of length 2. Gives the year of the first and last observation to be included in the covariate array.
#' @param inpolymat Logical matrix indicating whether each pixel on the observation grid is within the gas field.
#' @param cov_name Character string. Name of list element containing covariate array. Default value of "cov".
#' @param is_constant Logical. Is this covariate constant over time? (If so only first three colums used)
#'
#' @return A listobject containing a covatiate array, and the grid values for x, y and t for each observation.
#'
make_cov_mats <- function(cov_df, df_t0, array_tlim = NULL, inpolymat = inPolyMat, cov_name = "cov", is_constant){

  if(is.null(array_tlim)){
    n_yrs_df <- ncol(cov_df) - 2
    array_ylim <- df_y0 + c(0, n_yrs_df)
  }

  nx   <- length(unique(cov_df$X))
  minx <- min(cov_df$X)
  maxx <- max(cov_df$X)
  ny   <- length(unique(cov_df$Y))
  miny <- min(cov_df$Y)
  maxy <- max(cov_df$Y)

  xgrid = seq(minx,maxx,length= nx) # x values of observation points
  ygrid = seq(miny,maxy,length= ny) # y values of observation points
  tgrid = seq(array_tlim[1], array_tlim[2]) # Measured on Jan 1st this year

  year <- seq_along(tgrid)
  cov_mats <- array(data = NA, dim = c(nx,ny,length(year)))

  for (j in year){
    if (is_constant){
    covs <- cov_df[[3]]
    } else {
    covs <- cov_df[[2 + j]]
    }

    cov_mat <- matrix(NA, nrow = nx, ncol = ny)
    for( i in seq_along(cov_df$X)){
      x <- (cov_df$X[i] - minx) / 500 + 1
      y <- (cov_df$Y[i] - miny) / 500 + 1
      cov_mat[x,y] <- covs[i]
    }

    notinpoly_NA <- inpolymat
    notinpoly_NA[notinpoly_NA == 0] <- NA

    cov_mats[,,j] <- cov_mat * notinpoly_NA
  }

  out <- list(cov_mats, xgrid, ygrid, tgrid)
  names(out) <- c(paste0(cov_name,"_mats"),'xgrid','ygrid','tgrid')
  return(out)
}

#' PLot one or more temporal slices of a covariate array
#'
#' @param cov_mats covariate array object
#' @param xlab (char) x axis label
#' @param ylab (char)y axis label
#' @param cov_name (char) Optional name of covariate in plot titles
#' @param main (char) Optional plot title
#' @param asp (dbl) aspect ratio of fields plot
#' @param plot_slices (int) optional vector of tine "slices" to plot
#' @param zlim (dbl, length 2) optional vector of (shared) limits on colour scale
#' @param ... other parameters to pass to fields::image.plot()
#'
#' @return invisibly returns NULL
plot_cov_mats <- function(cov_mats, xlab = "Easting (m)" ,ylab = "Northing (m)", cov_name = "covariate values", main = NULL, asp = 1, plot_slices = seq_along(cov_mats$tgrid), zlim = NULL, ...){

  if(is.null(zlim)){zlim = range(cov_mats[[1]][,,plot_slices], na.rm = TRUE)}

  for (slice in plot_slices ){
    if(is.null(main)){
      title_text <- paste0(cov_name, " 01-01-", cov_mats$tgrid[slice])
    }

    fields::image.plot(
        x = sort(cov_mats$xgrid),
        y = sort(cov_mats$ygrid),
        z = cov_mats[[1]][,,slice],
        main = title_text,
        xlab = xlab,
        ylab = ylab,
        asp = asp,
        zlim = zlim,
        ...
    )
  }
  invisible(NULL)
}


#' Function to multiply two covariate arrays. Note that coordinates are inherited from the first argument.
#'
#' @param cov_mats_1 A covariate array object
#' @param cov_mats_2 A covariate array object
#' @param new_name (char) name for array element of new covariate array object
#'
#' @return A covariate array object, the product of two similar objects
multiply_cov_mats <- function(cov_mats_1, cov_mats_2, new_name){
  print("Remember to check conforming grids.")

  product <- cov_mats_1
  product[[1]] <- cov_mats_1[[1]] * cov_mats_2[[1]]
  names(product)[1] <- new_name

  product
}

# divide_cov_mats <- function(cov_mats_1, cov_mats_2, new_name){
#   print("Remember to check conforming grids.")
#
#   product <- cov_mats_1
#   product[[1]] <- cov_mats_1[[1]] / cov_mats_2[[1]]
#   names(product)[1] <- new_name
#
#   product
# }

#' Calculate the change in a covariate array over each epoch
#'
#' @param cov_mats A covariate array object
#' @param new_name (char) Name for array component of new covariate array object
#'
#' @return A covariate array object of the same dimensions as cov_mats. Note the last slice will be filled with NAs
temporal_difference_cov_mats <- function(cov_mats, new_name= NULL){
  if(is.null(new_name)){ new_name = paste0(names(cov_mats)[1],"_rate") }

  out <- cov_mats
  n_slices <- dim(out[[1]])[3]
  out[[1]] <- out[[1]][,,-1] - out[[1]][,,-n_slices]
  abind::abind(out[[1]], out[[1]][,,1]*NA, along = 3)
  out
}

replace_negatives_cov_mats <- function(cov_mats, replace,replacement_value){
  n_to_change <- sum(cov_mats[[1]] <= 0, na.rm = TRUE)
  if(replace){
    cov_mats[[1]][cov_mats[[1]] <= 0] <- replacement_value
  }
  print(paste(n_to_change, "values changed."))
  cov_mats
}

subset_to_years_cov_mats<- function(cov_mats, year_min, year_max){
  slice_min <- which(cov_mats$tgrid == year_min)
  slice_max <- which(cov_mats$tgrid == year_max)

  out <- cov_mats
  out[[1]] <- out[[1]][,,slice_min:slice_max]
  out$tgrid <- out$tgrid[slice_min:slice_max]

  return(out)
}
