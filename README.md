
<!-- README.md is generated from README.Rmd. Please edit that file -->

# spatio-temporal_threshold_selection

<!-- badges: start -->
<!-- badges: end -->

R code used to output figures and tables in the paper “Spatio-temporal
modelling of extreme induced seismicity in the presence of an evolving
measurement network”.

## Dependencies

To run this code, several R packages are required which may be installed
using the following code:

``` r
    required_pkgs <- c(
          "ggplot2",
          "ggspatial",
          "pracma",
          "cowplot",
          "dplyr",
          "purrr",
          "patchwork",
          "stringr"
        )
        
    install.packages(required_pkgs)
```

## Repository Overview

The structure of the repository is outlined below.

### `/src`

`distance_to_nearest_geo.R` contains the function to calculate the
distance to the ith closest geophone for a set of supplied earthquake
locations. The resulting distance function is used as the covariate in
the estimated spatio-temporal threshold.

`eqd_geo.R` contains R code to estimate a spatio-temporal threshold
function by choosing the most appropriate combination of parameters from
the supplied candidates such that the excesses of the resulting
threshold function can be closely modelled by a covariate-dependent
Generalised Pareto distribution (GPD). The selection is done by
utilising the expected quantile discrepancy (EQD) method.

`helper_functions.R` contains functions for the GPD which feed into
`eqd_geo.R`.

`intensity_estimation.R` contains functions to fit a stress-dependent
intensity function (based on a Poisson process likelihood) to the
earthquake locations.

### `/data`

This section contains the datasets used in our analysis, namely the
Groningen earthquake catalogue, the geophone dataset and the stress
covariate data. The covariate date is split into yearly files to allow
upload to this repository.

### `/analyses`

`main_figures.R` contains code to reproduce Figures 1-9 of the main
text.

`supplementary_figures` contains code to generate the tables and figures
contained in the supplementary material.

`threshold_estimation.R` contains the code to utilise `eqd_geo.R` to
select the appropriate threshold parameter values for each of the 12
threshold function formulations discussed in Section 4.1 of the main
text.

`threshold_results.R` provides code to extract the EQD values from each
threshold fit. These values appear in Table 1 in Section 6.1 of the main
text.

### `/outputs`

`/figures` contains the generated figures from the main text and
supplementary material.

`/threshold_results` contains the chosen thresholds, estimated GPD
parameter values, number of excesses and EQD distances for each of the
12 threshold function formulations. It also contains the same run with
B=1000 bootstraps in the `/boot1000` folder.

### `/in_development`

These code files were used to generate output for the paper, however the
code still needs cleaning and may be difficult to understand currently.

`future_inference.R` provides functions used in estimating future
endpoints and design levels.

`gron_eq_cat_covariate_inclusion.R` contains code to process the stress
data and extract stress values for each observed earthquake location.

`uncertainty_algorithms.R` contains code relevant for utilising the
uncertainty algorithms developed in Section 5 of the paper.

`uncertainty_algorithms_STORM.R` contains code to run the uncertainty
algorithms on a computing cluster.

`/uncertainty` contains all the results from running the different
uncertainty algorithms.

### Contact

If you have questions, please contact <murphyconor148@gmail.com> .
Please include “Induced seismicity code” in the subject of the email.
