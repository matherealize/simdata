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
#' \code{f} of the list entry is the function to be applied via \code{\link{do.call}}.
#' The list entry itself is another named list, specifying the arguments
#' to the function \code{f} as named arguments.
#'
#' The functions must take a matrix or data.frame as first argument and
#' return another matrix or data.frame of the same dimensions as
#' single output.
#'
#' Examples of post-processing steps are truncation (\code{\link{process_truncate}})
#' or centering / standardizing data (via \code{\link{scale}}, see example section
#' below).
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
#' @examples
#' \dontrun{
#' process_data(diag(5), 
#'   functions = list(scale = list(center = TRUE, scale = FALSE)))
#' }
#' 
#' @seealso 
#' \code{\link{process_truncate}}
#'
#' @export
process_data <- function(X, functions = list()) {

    # apply post-processing functions
    for (f in names(functions)) {
        X = do.call(f, modifyList(list(X), functions[[f]]))
    }

    X
}

# Post-processing Functions #########################################
#' @title Truncate columns of datamatrix
#'
#' @param X
#' Matrix or Data.frame.
#' @param truncate_multipliers
#' Vector of truncation parameters. Either a single value which is
#' replicated as necessaary or of same dimension as \code{ncol(X)}.
#' If any vector entry is NA, the corresponding column will not be
#' truncated. See details.
#'
#' @details
#' Truncation is processed as follows:
#' \enumerate{
#' \item Compute the 1st and 3rd quartile q1 / q3 of variables in _X_.
#' \item Multiply these quantities by values in \code{truncate_final} to obtain
#' L and U. If a value is NA, the corresponding variable is not truncated.
#' \item Set any value smaller / larger than L / U to L / U.
#' }
#'
#' @return
#' Matrix or Data.frame of same dimensions as input.
#'
#' @export
process_truncate <- function(X, truncate_multipliers = NA) {

    # if single input -> replicate to fit X
    truncate_multipliers = rep_len(truncate_multipliers, ncol(X))

    # which columns should be truncated
    do_trunc = which(!sapply(truncate_multipliers, is.na))

    if (length(do_trunc) > 0) {
        # determine thresholds to truncate at
        quantiles = apply(X[, do_trunc, drop = FALSE], 2,
                          quantile, c(0.25, 0.75))
        iqr = apply(quantiles, 2, diff)
        truncate_upper = quantiles[2, ] +
            truncate_multipliers[do_trunc] * iqr
        truncate_lower = quantiles[1, ] -
            truncate_multipliers[do_trunc] * iqr

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

