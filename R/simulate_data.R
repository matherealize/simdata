# TODO: Update doc

# Data Simulation ###################################################
#' @title Simulate design matrix
#'
#' @description
#' Generate simulated dataset based on transformation of
#' multivariate gaussian distribution.
#'
#' @param relations
#' Correlation / Covariance matrix of the initial multivariate
#' gaussian distribution Z or object of class \code{\link{design}}.
#' @param n_obs
#' Number of simulated observations.
#' @param transform
#' List of functions. Specifies transformation of underlying
#' datamatrix Z to final datamatrix X. If not specified, initial
#' data Z is returned. If named list, \code{names(transform)} will be used
#' for columnnames of X. See details for more information.
#' @param mean_initial
#' Vector of mean values of the initial multivariate gaussian
#' distribution Z. Dimension needs to correspond to dimension
#' of \code{relations}.
#' @param sd_initial
#' Vector of standard deviations of the initial multivariate
#' gaussian distribution Z. Dimension needs to correspond to dimension
#' of \code{relations}.
#' @param is_correlation
#' If TRUE, then \code{relations} specifies a correlation matrix (default,
#' this type of specification is usually more natural than specifying
#' a covariance matrix). Otherwise, \code{relations} specifies a
#' covariance matrix.
#' @param names_final
#' Variable names for final datamatrix X. Needs to have same
#' dimension as \code{transform}. Overrides other naming options.
#' @param prefix_final
#' Prefix attached to variables in final datamatrix X. Overriden
#' by \code{names} argument or using a named list for \code{transform}.
#' @param process_final
#' List of lists specifying post-processing functions applied to final
#' datamatrix X before returning it. See \code{\link{process_data}}.
#' @param seed
#' Set random seed to ensure reproducibility of results.
#'
#' @return
#' Data.frame with simulated datamatrix X.
#'
#' @details
#' Data is generated using the following procedure:
#' \enumerate{
#' \item The underlying data matrix Z is sampled from a
#' multivariate gaussian distribution (number of dimensions specified by
#' dimensions of \code{relations}).
#' \item Z is then transformed into the final data matrix X by applying
#' functions from \code{transform} to the columns of Z (final number of
#' variables given by length of \code{transform}).
#' \item X is post-processed if specified (truncation to avoid
#' outliers).
#' }
#'
#' Transformations are specified as a list of functions, which take
#' the whole datamatrix Z as single argument. For example, to multiply
#' column 2 of Z by 2, use function(Z) Z[,2] * 2.
#'
#' Post-processing the datamatrix is based on \code{\link{process_data}}.
#'
#' @note
#' Note that \code{relations} specifies the correlation / covariance
#' of the underlying gaussian data Z and thus does not directly translate into
#' correlations between the variables of the final datamatrix X.
#'
#' This function is best used in conjunction with the \code{design} S3 class,
#' which facilitates further data visualization and conveniently stores
#' information as a template for simulation tasks.
#'
#' @seealso
#' \code{\link{design}}, \code{\link{conditional_simulate_data}},
#' \code{\link{process_data}}
#'
#' @export
simulate_data <- function(generator, ...) {
    UseMethod("simulate_data", generator)
}

#' @describeIn simulate_data Workhorse function to be used if no \code{design}
#' S3 class is used.
#'
#' @export
simulate_data.default <- function(generator, 
                                  n_obs = 1, 
                                  transform = base::identity,
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
    x = transform(x)

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

#' @describeIn simulate_data Function to be used with \code{design} S3 class.
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
#' until the conditions in \code{reject} are met?
#' @param on_reject
#' If "stop", an error is returned if after \code{reject_max_iter} times no
#' suitable datamatrix X could be found. If "current", the current datamatrix
#' is returned, regardless of the conditions in \code{reject}.
#' Otherwise, \code{NULL} is returned. In each case a warning is reported.
#' @param ...
#' All further parameters are passed to \code{\link{simulate_data}}.
#'
#' @return
#' Data.frame with simulated datamatrix X if a suitable dataset can be found
#' or the iteration limit is hit.
#'
#' @details
#' For details on generating and post-processing datasets, see
#' \code{\link{simulate_data}}. This function simulates data conditional
#' on certain requirements that must be met by the final datamatrix X.
#' This checking is conducted on the output of \code{simulate_data} (i.e.
#' also includes possible post-processing steps).
#'
#' @section Rejecting Datasets:
#' Examples for restrictions include
#' variance restrictions (e.g. no constant columns which could happen due
#' to extreme transformations of the initial gaussian distribution Z), ensuring
#' a sufficient number of observations in a given class (e.g. certain
#' binary variables should have at least x\% events) or preventing
#' multicollinearity (e.g. X must have full column rank). If one of the
#' functions in \code{reject} evaluates to \code{FALSE}, the current datamatrix
#' X is rejected.
#' In case of rejection, new datasets can be simulated until the conditions
#' are met or a given maximum iteration limit is hit (\code{reject_max_iter}),
#' after which the last datamatrix is returned or an error is reported.
#'
#' @section Rejection Functions:
#' Rejection function templates are found in \code{\link{is_collinear}} and
#' \code{\link{contains_constant}}.
#'
#' @note
#' Note that \code{relations} specifies the correlation / covariance
#' of the underlying gaussian data Z and thus does not directly translate into
#' correlations between the variables of the final datamatrix X.
#'
#' This function is best used in conjunction with the \code{design} S3 class,
#' which facilitates further data visualization and conveniently stores
#' information as a template for simulation tasks.
#'
#' @seealso
#' \code{\link{design}}, \code{\link{simulate_data}},
#' \code{\link{is_collinear}}, \code{\link{contains_constant}}
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
