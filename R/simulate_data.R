# Data Simulation ###################################################
#' @title Simulate design matrix
#'
#' @description
#' Generate simulated dataset based on transformation of
#' an underlying base distribution.
#' 
#' @template simulate_data_template
#' @param n_obs
#' Number of simulated observations.
#' @param seed
#' Set random seed to ensure reproducibility of results.
#' @param ...
#' Further arguments passed to `generator` function.
#'
#' @details
#' Data is generated using the following procedure:
#' \enumerate{
#' \item An underlying dataset `Z` is sampled from some distribution. This is 
#' done by a call to the `generator` function. 
#' \item `Z` is then transformed into the final dataset `X` by applying the
#' `transform` function to `Z`.
#' \item `X` is post-processed if specified (e.g. truncation to avoid
#' outliers).
#' }
#' 
#' @section Generators:
#' The `generator` function is assumend to provide the same interface
#' as the random generation functions in the R \pkg{stats} and \pkg{extraDistr}
#' packages. Specifically, that means it takes the number of observations as
#' first argument. All further arguments can be set via passing them as
#' named argument to this function.
#'
#' @section Transformations:
#' Transformations should be applicable to the output of the `generator`
#' function (i.e. take a data.frame or matrix as input) and output another
#' data.frame or matrix. A convenience function `\link{function_list}` is 
#' provided by this package to specify transformations as a list of functions,
#'  which take the whole datamatrix `Z` as single argument and can be used to
#'  apply specific transformations to the columns of that matrix. See the 
#'  documentation for `\link{function_list}` for details.
#'  
#' @section Post-processing:
#' Post-processing the datamatrix is based on `\link{process_data}`.
#'
#' @note
#' This function is best used in conjunction with the `\link{simdesign}`
#' S3 class or any template based upon it, which facilitates further data 
#' visualization and conveniently stores information as a template for simulation tasks.
#' 
#' @return
#' Data.frame or matrix with `n_obs` rows for simulated dataset `X`.
#'
#' @seealso
#' `\link{simdesign}`, 
#' `\link{mvtnorm_simdesign}`, 
#' `\link{simulate_data_conditional}`,
#' `\link{process_data}`
#'
#' @export
simulate_data <- function(generator, ...) {
    UseMethod("simulate_data", generator)
}

#' @describeIn simulate_data Function to be used if no `\link{simdesign}`
#' S3 class is used.
#'
#' @export
#' @method simulate_data default
simulate_data.default <- function(generator, 
                                  n_obs = 1, 
                                  transform_initial = base::identity,
                                  names_final = NULL,
                                  prefix_final = NULL,
                                  process_final = list(),
                                  seed = NULL, 
                                  ...) {

    if (!is.null(seed))
        set.seed(seed)

    # generate initial dataset
    x = generator(n_obs, ...)
    
    # transform to final dataset
    x = transform_initial(x)

    # apply post-processing functions
    x = process_data(x, process_final)

    # rename columns
    if (!is.null(names_final)) {
        colnames(x) = names_final
    } else if (!is.null(prefix_final)) {
        colnames(x) = paste0(prefix_final, 1:ncol(x))
    }

    x
}

#' @describeIn simulate_data Function to be used with `\link{simdesign}` S3 class.
#' 
#' @param apply_transformation
#' If `generator` is a `simdesign` object, then this argument can be set to 
#' FALSE to override the stored information and not transform and process data. 
#' Thus, the raw data from the design generator is returned. This can be useful
#' for debugging purposes, but is usually unnecessary to be changed.
#' @param apply_processing
#' If `generator` is a `simdesign` object, then this argument can be set to 
#' FALSE to override the stored information and not process the data after
#' the initial data is transformed. This can be useful for debugging purposes, 
#' but is usually unnecessary to be changed.
#'
#' @export
#' @method simulate_data simdesign
simulate_data.simdesign <- function(generator,
                                    n_obs,
                                    seed = NULL,
                                    apply_transformation = TRUE, 
                                    apply_processing = TRUE) {
    
    transform_initial = design$transform_initial
    names_final = design$names_final
    process_final = design$process_final
    
    if (!apply_transformation) {
        transform_initial = base::identity
    }
    
    if (!all(apply_processing, apply_transformation)) {
        names_final = NULL
        process_final = list()
    }
    
    simulate_data(generator = design$generator,
                  n_obs = n_obs,
                  transform_initial = transform_initial,
                  names_final = names_final,
                  process_final = process_final,
                  seed = seed)
}

# TODO test where design is used to transform data to contain constant column
# or collinear!!
# Conditional Data Simulation #######################################
#' @title Simulate data which satisfies certain conditions
#'
#' @description
#' Generate simulated dataset based on transformation of
#' an underlying base distribution while checking that certain
#' conditions are met.
#'
#' @template simulate_data_template
#' @param n_obs
#' Number of simulated observations.
#' @param reject
#' Function which takes a matrix or data.frame `X` as single input and outputs 
#' TRUE or FALSE. Specifies when a simulated final datamatrix `X` should
#' be rejected. Functions must output TRUE if condition IS NOT met / FALSE if 
#' condition IS met and matrix can be accepted. Intended to be used with 
#' `\link{function_list}`. See details.
#' @param reject_max_iter
#' Intger > 0. In case of rejection, how many times should a new datamatrix be
#' simulated until the conditions in `reject` are met?
#' @param on_reject
#' If "stop", an error is returned if after `reject_max_iter` times no
#' suitable datamatrix X could be found. If "current", the current datamatrix
#' is returned, regardless of the conditions in `reject`.
#' Otherwise, NULL is returned. In each case a warning is reported.
#' @param ...
#' All further parameters are passed to `\link{simulate_data}`.
#'
#' @details
#' For details on generating, transforming and post-processing datasets, see
#' `\link{simulate_data}`. This function simulates data conditional
#' on certain requirements that must be met by the final datamatrix `X`.
#' This checking is conducted on the output of `simulate_data` (i.e.
#' also includes possible post-processing steps).
#'
#' @section Rejecting Datasets:
#' Examples for restrictions include
#' variance restrictions (e.g. no constant columns which could happen due
#' to extreme transformations of the initial gaussian distribution `Z`), ensuring
#' a sufficient number of observations in a given class (e.g. certain
#' binary variables should have at least x\% events) or preventing
#' multicollinearity (e.g. `X` must have full column rank). If `reject` 
#' evaluates to FALSE, the current datamatrix `X` is rejected.
#' In case of rejection, new datasets can be simulated until the conditions
#' are met or a given maximum iteration limit is hit (`reject_max_iter`),
#' after which the latest datamatrix is returned or an error is reported.
#'
#' @section Rejection Function:
#' The `reject` function should take a single input (a data.frame or matrix) and
#' output TRUE if the dataset is to be rejected or FALSE if it is to be accepted.
#' This package provides the `\link{function_list}` convenience function which
#' allows to easily create a rejection function which assesses several 
#' conditions on the input dataset by simply passing individual test functions
#' to `function_list`. Such test function templates are found in 
#' `\link{is_collinear}` and `\link{contains_constant}`. See the example
#' below.
#' 
#' @return
#' Data.frame or matrix with `n_obs` rows for simulated dataset `X` if all
#' conditions are met within the iteration limit. Otherwise NULL.
#'
#' @seealso
#' `\link{simdesign}`, 
#' `\link{simulate_data}`,
#' `\link{function_list}`,
#' `\link{is_collinear}`, 
#' `\link{contains_constant}`
#' 
#' @examples
#' \dontrun{
#' dsgn = mvtnorm_simdesign(diag(5))
#' simulate_data_conditional(dsgn, 100, 
#'    reject = function_list(is_collinear, contains_constant))
#' }
#'
#' @export
simulate_data_conditional <- function(generator, 
                                      n_obs = 1,
                                      reject = function(x) TRUE,
                                      reject_max_iter = 10,
                                      on_reject = "ignore", 
                                      ...) {
    simulation_success = FALSE

    while (!simulation_success) {
        reject_max_iter = reject_max_iter - 1
        
        x = simulate_data(generator, n_obs, ...)
        
        # check rejections
        simulation_success = !any(unlist(reject(x)))

        if (reject_max_iter == 0) {
            if (!simulation_success) {
                if (on_reject == "stop") {
                    stop("No suitable datamatrix found within iteration limit.\n")
                } else if (on_reject == "current") {
                    warning("No suitable datamatrix found within iteration limit. Returning current candidate matrix.\n")
                    simulation_success = TRUE
                } else {
                    warning("No suitable datamatrix found within iteration limit. Returning NULL.\n")
                    simulation_success = TRUE
                    x = NULL
                }
            }
        }
    }

    x
}
