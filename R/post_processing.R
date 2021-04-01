# Post-processing Automation ########################################
#' @title Post-processing of datamatrix
#'
#' @description
#' Applies functions to a matrix or data.frame.
#'
#' @param X
#' Matrix or Data.frame.
#' @param functions
#' List of lists, specifying functions to be applied as well as their
#' arguments. See details.
#'
#' @details
#' Functions are passed into the post-processor as a named list. The name
#' `f` of the list entry is the function to be applied via 
#' \code{\link[base:do.call]{base::do.call}}.
#' The list entry itself is another named list, specifying the arguments
#' to the function `f` as named arguments.
#'
#' The functions must take a matrix or data.frame as first argument and
#' return another matrix or data.frame of the same dimensions as
#' single output.
#'
#' Examples of post-processing steps are truncation
#' (\code{\link{process_truncate}}) or centering / standardizing data 
#' (via \code{\link{scale}}, see example section below).
#'
#' Can be useful to apply on simulated datasets, even outside of the
#' simulation function (e.g. when standardization is only required at the
#' modeling step).
#'
#' @note
#' Use with caution - no error checking is done for now so the user has
#' to take care of everything themselves! Furthermore, output of the
#' functions is not checked either.
#' 
#' @return 
#' Matrix or data.frame with post-processing applied.
#'
#' @examples
#' \dontrun{
#' do_processing(diag(5), 
#'   functions = list(scale = list(center = TRUE, scale = FALSE)))
#' }
#' 
#' @seealso 
#' \code{\link{process_truncate}}
#'
#' @export
do_processing <- function(X, functions = list()) {

    # apply post-processing functions
    for (f in names(functions)) {
        X = do.call(f, utils::modifyList(list(X), functions[[f]]))
    }

    X
}

# Post-processing Functions #########################################
#' @title Truncate columns of datamatrix at datamatrix specific thresholds
#' 
#' @description 
#' Truncation based on the interquartile range to be applied to a dataset.
#'
#' @param X
#' Matrix or Data.frame.
#' @param truncate_multipliers
#' Vector of truncation parameters. Either a single value which is
#' replicated as necessary or of same dimension as `ncol(X)`.
#' If any vector entry is NA, the corresponding column will not be
#' truncated. If named, then the names must correspond to columnnames in `X`, 
#' and only specified columns will be processed. See details.
#' @param only_numeric
#' If TRUE and if `X` is a data.frame, then only columns of type `numeric` will
#' be processed. Otherwise all columns will be processed (e.g. also in the 
#' case that `X` is a matrix).
#'
#' @details
#' Truncation is processed as follows:
#' \enumerate{
#' \item Compute the 1st and 3rd quartile q1 / q3 of variables in `X`.
#' \item Multiply these quantities by values in `truncate_multipliers` to obtain
#' _L_ and _U_. If a value is NA, the corresponding variable is not truncated.
#' \item Set any value smaller / larger than _L_ / _U_ to _L_ / _U_.
#' }
#' 
#' Truncation multipliers can be specified in three ways (note that whenever
#' `only_numeric` is set to TRUE, then only numeric columns are affected): 
#' 
#' \itemize{
#' \item A single numeric - then all columns will be processed in the same way
#' \item A numeric vector without names - it is assumed that the length can be
#' replicated to the number of columns in `X`, each column is processed by the 
#' corresponding value in the vector
#' \item A numeric vector with names - length can differ from the columns in 
#' `X` and only the columns for which the names occur in the vector are 
#' processed
#' }
#'
#' @return
#' Matrix or data.frame of same dimensions as input.
#'
#' @export
process_truncate_by_iqr <- function(X, 
                                    truncate_multipliers = NA,
                                    only_numeric = TRUE) {
    
    # if input vector has names -> remember them
    truncate_names = NULL
    if (!is.null(names(truncate_multipliers)) & !is.null(names(X)))
        truncate_names = names(truncate_multipliers)

    # truncation vector
    truncate_vector = rep_len(NA, ncol(X))
    names(truncate_vector) = names(X)
    
    # set values in vector if available
    if (!is.null(truncate_names)) {
        truncate_vector[truncate_names] = truncate_multipliers
    } else truncate_vector = rep_len(truncate_multipliers, ncol(X))
    
    if (is.data.frame(X) & only_numeric) {
        ind_numeric = sapply(X, class) == "numeric"
        truncate_vector[!ind_numeric] = NA
    }

    # which columns should be truncated
    do_trunc = which(!sapply(truncate_vector, is.na))

    if (length(do_trunc) > 0) {
        # determine thresholds to truncate at
        quantiles = apply(X[, do_trunc, drop = FALSE], 2,
                          stats::quantile, c(0.25, 0.75))
        iqr = apply(quantiles, 2, diff)
        truncate_upper = quantiles[2, ] +
            truncate_vector[do_trunc] * iqr
        truncate_lower = quantiles[1, ] -
            truncate_vector[do_trunc] * iqr

        # do truncation
        for (i in seq_along(do_trunc)) {
            X[X[, do_trunc[i]] > truncate_upper[i], do_trunc[i]] =
                truncate_upper[i]
            X[X[, do_trunc[i]] < truncate_lower[i], do_trunc[i]] =
                truncate_lower[i]
        }
    }

    X
}

#' @title Truncate columns of datamatrix at specified thresholds
#' 
#' @description 
#' Truncation based on fixed thresholds to be applied to a dataset. Allows
#' to implement truncation by measures derived from the overall data generating
#' mechanism.
#'
#' @param X
#' Matrix or Data.frame.
#' @param truncate_lower,truncate_upper
#' Vectors of truncation parameters, i.e. lower and upper tresholds for 
#' truncation. 
#' Either a single value which is replicated as necessary or of same dimension
#' as `ncol(X)`. If any vector entry is NA, the corresponding column will not be
#' truncated. Truncation at lower and upper thresholds is treated independently.
#' If named, then the names must correspond to columnnames in `X`, 
#' and only specified columns will be processed. See details.
#' @param only_numeric
#' If TRUE and if `X` is a data.frame, then only columns of type `numeric` will
#' be processed. Otherwise all columns will be processed (e.g. also in the 
#' case that `X` is a matrix).
#'
#' @details
#' Truncation is defined by setting all values below or above the truncation
#' threshold to the truncation threshold.
#' 
#' Truncation parameters can be specified in three ways (note that whenever
#' `only_numeric` is set to TRUE, then only numeric columns are affected): 
#' 
#' \itemize{
#' \item A single numeric - then all columns will be processed in the same way
#' \item A numeric vector without names - it is assumed that the length can be
#' replicated to the number of columns in `X`, each column is processed by the 
#' corresponding value in the vector
#' \item A numeric vector with names - length can differ from the columns in 
#' `X` and only the columns for which the names occur in the vector are 
#' processed
#' }
#'
#' @return
#' Matrix or data.frame of same dimensions as input.
#'
#' @export
process_truncate_by_threshold <- function(X, 
                                          truncate_lower = NA,
                                          truncate_upper = NA,
                                          only_numeric = TRUE) {
    
    # if input vectors have names -> remember them
    truncate_lower_names = NULL
    if (!is.null(names(truncate_lower)) & !is.null(names(X)))
        truncate_lower_names = names(truncate_lower)
    truncate_upper_names = NULL
    if (!is.null(names(truncate_upper)) & !is.null(names(X)))
        truncate_upper_names = names(truncate_upper)
    
    # truncation vectors
    truncate_lower_vector = rep_len(NA, ncol(X))
    names(truncate_lower_vector) = names(X)
    truncate_upper_vector = rep_len(NA, ncol(X))
    names(truncate_upper_vector) = names(X)
    
    # set values in vector if available
    if (!is.null(truncate_lower_names)) {
        truncate_lower_vector[truncate_lower_names] = truncate_lower
    } else truncate_lower_vector = rep_len(truncate_lower, ncol(X))
    if (!is.null(truncate_upper_names)) {
        truncate_upper_vector[truncate_upper_names] = truncate_upper
    } else truncate_upper_vector = rep_len(truncate_upper, ncol(X))
    
    if (is.data.frame(X) & only_numeric) {
        ind_numeric = sapply(X, class) == "numeric"
        truncate_lower_vector[!ind_numeric] = NA
        truncate_upper_vector[!ind_numeric] = NA
    }
    
    # which columns should be truncated
    do_trunc = which(!sapply(truncate_lower_vector, is.na))
    if (length(do_trunc) > 0) {
        # do truncation
        for (i in do_trunc) {
            X[X[, i] < truncate_lower_vector[i], i] =
                truncate_lower_vector[i]
        }
    }
    
    do_trunc = which(!sapply(truncate_upper_vector, is.na))
    if (length(do_trunc) > 0) {
        # do truncation
        for (i in do_trunc) {
            X[X[, i] > truncate_upper_vector[i], i] =
                truncate_upper_vector[i]
        }
    }
    
    X
}