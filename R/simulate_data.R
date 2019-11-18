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
#' Further arguments passed to `generator`` function.
#' 
#' @return
#' Data.frame or matrix with `n_obs` rows for simulated dataset `X`.
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
#' @seealso
#' `\link{simdesign}`, 
#' `\link{mvtnorm_simdesign}`, 
#' `\link{conditional_simulate_data}`,
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

    # generate initial data matrix
    x = generator(n_obs, ...)
    
    # transform to final data matrix
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

#' @describeIn `simulate_data` Function to be used with `\link{simdesign}` S3 class.
#'
#' @export
simulate_data.simdesign <- function(design,
                                    n_obs,
                                    seed = NULL) {
    simulate_data(generator = design$generator,
                  n_obs = n_obs,
                  transform = design$transform,
                  names_final = design$names_final,
                  process_final = design$process_final,
                  seed = seed)
}


# TODO
# Conditional Data Simulation #######################################
#' @title Simulate data which satisfies certain conditions
#'
#' @description
#' Generate simulated dataset based on transformation of
#' multivariate gaussian distribution while checking certain
#' conditions are met.
#'
#' @param reject
#' List of functions. Specifies when a simulated final datamatrix X should
#' be rejected. Functions must take X as single input and output TRUE if
#' condition IS NOT met / FALSE if condition IS met and matrix can be
#' accepted. See details.
#' @param reject_max_iter
#' In case of rejection, how many times should a new datamatrix be simulated
#' until the conditions in `reject` are met?
#' @param on_reject
#' If "stop", an error is returned if after `reject_max_iter` times no
#' suitable datamatrix X could be found. If "current", the current datamatrix
#' is returned, regardless of the conditions in `reject`.
#' Otherwise, NULL is returned. In each case a warning is reported.
#' @param ...
#' All further parameters are passed to `\link{simulate_data}`.
#'
#' @return
#' Data.frame with simulated datamatrix X if a suitable dataset can be found
#' or the iteration limit is hit.
#'
#' @details
#' For details on generating and post-processing datasets, see
#' `\link{simulate_data}`. This function simulates data conditional
#' on certain requirements that must be met by the final datamatrix X.
#' This checking is conducted on the output of `simulate_data` (i.e.
#' also includes possible post-processing steps).
#'
#' @section Rejecting Datasets:
#' Examples for restrictions include
#' variance restrictions (e.g. no constant columns which could happen due
#' to extreme transformations of the initial gaussian distribution Z), ensuring
#' a sufficient number of observations in a given class (e.g. certain
#' binary variables should have at least x\% events) or preventing
#' multicollinearity (e.g. X must have full column rank). If one of the
#' functions in `reject` evaluates to FALSE, the current datamatrix
#' X is rejected.
#' In case of rejection, new datasets can be simulated until the conditions
#' are met or a given maximum iteration limit is hit (`reject_max_iter`),
#' after which the last datamatrix is returned or an error is reported.
#'
#' @section Rejection Functions:
#' Rejection function templates are found in `\link{is_collinear}` and
#' `\link{contains_constant}`.
#'
#' @note
#' Note that `relations` specifies the correlation / covariance
#' of the underlying gaussian data Z and thus does not directly translate into
#' correlations between the variables of the final datamatrix X.
#'
#' This function is best used in conjunction with the `design` S3 class,
#' which facilitates further data visualization and conveniently stores
#' information as a template for simulation tasks.
#'
#' @seealso
#' `\link{design}`, 
#' `\link{simulate_data}`,
#' `\link{is_collinear}`, 
#' `\link{contains_constant}`
#'
#' @export
conditional_simulate_data <- function(relations, n_obs,
                                      reject = NULL,
                                      reject_max_iter = 10,
                                      on_reject = "ignore", ...) {
    simulation_success = FALSE

    while (!simulation_success) {
        x = simulate_data(relations, n_obs, ...)

        # check rejections
        if (is.null(reject)) {
            simulation_success = TRUE
        } else {
            simulation_success = all(!sapply(reject, function(f, x) f(x), x))
        }

        reject_max_iter = reject_max_iter - 1

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

    return(x)
}
