#' @param generator 
#' Function which generates data from the underlying base distribution. It is
#' assumed it takes the number of simulated observations `n_obs` as first
#'  argument, as all random generation functions in the \pkg{stats} and 
#' \pkg{extraDistr} do. Furthermore, it is expected to return a two-dimensional
#' array as output (matrix or data.frame). See details.
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
#' datamatrix `X` before returning it. See \code{\link{do_processing}}.
