#' @param generator 
#' Function which generates data from the underlying base distribution. It is
#' assumend it takes the number of simulated observations `n_obs` as first
#'  argument, as all random generation functions in the \pkg{stats} and 
#' \pkg{extraDistr} do. See details.
#' @param transform_initial
#' Function which specifies the transformation of the underlying
#' dataset `Z` to final dataset `X`. See details.
#' @param names_final
#' NULL or character vector with variable names for final dataset `X`. 
#' Length needs to equal the number of columns of `X`. Overrides other naming options. 
#' @param prefix_final
#' NULL or prefix attached to variables in final dataset `X`. Overriden
#' by `names_final` argument. Set to NULL if no prefixes should 
#' be added. 
#' @param process_final
#' List of lists specifying post-processing functions applied to final
#' datamatrix `X` before returning it. See `\link{process_data}`.
