# Templates for datamatrix rejection functions ######################
#' @title Check if matrix contains constant column(s)
#'
#' @param x
#' Matrix or Data.frame.
#' @param eps
#' Threshold for standard deviation below which a column is considered to be constant. 
#'
#' @return
#' TRUE if one of the columns has standard deviation of below `eps``, else FALSE.
#'
#' @note
#' Prints a warning if constant is found.
#'
#' @export
contains_constant <- function(x, eps = .Machine$double.eps) {
    x = as.matrix(x)
    val = any(apply(x, 2, sd) < eps)
    
    if (val)
        warning("contains_constant: Matrix contains constant column.\n")
    
    val
}

#' @title Check if matrix is collinear
#'
#' @param x
#' Matrix or Data.frame.
#'
#' @return
#' TRUE if matrix is collinear, else FALSE.
#'
#' @note
#' Prints a warning if collinear.
#'
#' @export
is_collinear <- function(x) {
    val = qr(x)$rank < ncol(x)
    
    if (val)
        warning("is_collinear: Matrix is not full rank.\n")
    
    val
}

# Correlation matrix utilities ######################################
#' @title Build correlation matrix
#'
#' @description Use to specify correlation matrix in convenient way
#' by giving entries of the upper triangular part.
#'
#' @param n_var
#' Integer, number of variables (= rows = columns of matrix).
#' @param entries
#' Matrix of correlation entries. Consists of 3 columns
#' (variable_1, variable_2, correlation) that specify both variables
#' and corresponding correlation in the upper triangular part of
#' the matrix (i.e. variable_1 < variable_2) .
#'
#' @return
#' Matrix with user supplied entries.
#'
#' @seealso
#' `\link{cor_to_upper}`
#'
#' @examples
#' \dontrun{
#' cor_from_upper(2, rbind(c(1,2,0.8)))
#' }
#'
#' @export
cor_from_upper <- function(n_var, entries = NULL) {
    # ones in the diagonal
    cor_mat = diag(n_var)

    if (!is.null(entries)) {
        # set correlation entries
        entries = matrix(entries, ncol = 3)
        entries[, c(1,2)] = t(apply(matrix(entries[, c(1,2), drop = FALSE], ncol = 2),
                                   1, sort))
        cor_mat[matrix(entries[, c(1,2)], ncol = 2)] = entries[, 3]
    }

    # symmetrize
    cor_mat[lower.tri(cor_mat)] = t(cor_mat)[lower.tri(cor_mat)]

    cor_mat
}

#' @title Convert correlation matrix to specification used by
#' `cor_from_upper`
#'
#' @param m
#' Symmetric correlation matrix.
#'
#' @return
#' Matrix with 3 columns (variable_1, variable_2, correlation), where
#' correlation gives the entry at position (variable_1, variable_2) of the
#' input correlation matrix.
#'
#' @seealso
#' `\link{cor_from_upper}`
#'
#' @export
cor_to_upper <- function(m) {
    res = matrix(
        data = c(
            unlist(lapply(1:(nrow(m) - 1), function(x) rep(x, ncol(m) - x))),
            unlist(lapply(2:ncol(m), function(x) x:ncol(m)))
        ),
        nrow = nrow(m) * (nrow(m) - 1) / 2, ncol = 2
    )
    res = cbind(res, m[res[, c(1,2), drop = FALSE]]) 
    res = res[abs(res[,3]) > 0, , drop = FALSE]
    
    res
}

#' @title Convert covariance matrix to correlation matrix
#' 
#' @description 
#' Wrapper for `\link[stats:cov2cor]{stats::cov2cor}`.
#' 
#' @param m
#' Symmetric covariance matrix.
#' 
#' @return 
#' Symmetric correlation matrix.
#' 
#' @seealso 
#' `\link{cor_to_cov}`
#' 
#' @export
cov_to_cor <- function(m) {
    cov2cor(m)
}

#' @title Convert correlation matrix to covariance matrix
#' 
#' @description 
#' Rescale correlation matrix by variable standard deviations to yield a
#' covariance matrix.
#' 
#' @param m
#' Symmetric correlation matrix.
#' 
#' @return 
#' Symmetric covariance matrix.
#' 
#' @seealso 
#' `\link{cov_to_cor}`
#' 
#' @export
cor_to_cov <- function(m, sds = NULL) {
    if (is.null(sds))
        sds = rep(1, nrow(m))
    
    diag(sds) %*% m %*% diag(sds)
}

# Function helper ###################################################
#' @title Apply list of functions to input
#' 
#' @param ...
#' Named or unnamed arguments, each of which is a function taking exactly 
#' one input. See details.
#' @param stringsAsFactors,check.names
#' Arguments of `\link[base:data.frame]{data.frame}`.
#' 
#' @details 
#' This is a convenience function which takes a number of functions and returns
#' another function which applies all of the user specified functions to a new 
#' input, and collects the results as list or data.frame.
#' This is useful to e.g. transform columns of a data.frame or check 
#' the validity of a matrix during simulations. See the example here and 
#' in `\link{simulate_data_conditional}`.
#' 
#' The assumptions for the individual functions are: 
#' 
#' \itemize{
#' \item Each function is expected to take a single input. 
#' \item Each function is expected to output a result consistent with the 
#' other functions (i.e. same output length) to ensure that the results can be
#' summarized as a data.frame.
#' } 
#' 
#' @note 
#' This function works fine without naming the input arguments, but the 
#' resulting data.frames have empty column names if that is the case. Thus, 
#' it is recommended to only pass named function arguments. 
#' 
#' @examples 
#' \dontrun{
#' f = function_list(
#'    v1 = function(x) x[,1] * 2, 
#'    v2 = function(x) x[,2] + 10)
#'    
#' f(diag(2))
#' } 
#' 
#' @return 
#' Function with a single input which outputs a data.frame.
#' 
#' @seealso 
#' `\link[base:data.frame]{data.frame}`
#' 
#' @export
function_list <- function(..., 
                          stringsAsFactors = default.stringsAsFactors(), 
                          check.names = TRUE) {
    fct_list = list(...)
    
    function(m) {
        do.call(data.frame, 
                append(
                    lapply(fct_list, function(f, x) f(x), m), 
                    list(
                        stringsAsFactors = stringsAsFactors, 
                        check.names = check.names
                        )
                    )
        )
    }
}

#' @title Extract individual functions from `function_list`
#' 
#' @description 
#' Extract individual function objects from environment of a `function_list` 
#' object.
#' 
#' @param flist
#' `function_list` object.
#' 
#' @return 
#' List with named or unnamed entries corresponding to individual function
#' objects that were passed to the `function_list` object.
get_from_function_list <- function(flist) {
    get("fct_list", envir = environment(flist))
}

#' @title Create `function_list` object from list of functions
#' 
#' @description 
#' Create a `function_list` object from a list of functions. This is useful
#' if such a list is created programmatically.
#' 
#' @param flist
#' List in which each entry is a function object. Can be named or unnamed.
#' @param ...
#' Passed to `\link{function_list}`.
#' 
#' @inherit function_list return
as_function_list <- function(flist, ...) {
    do.call(function_list, append(flist, list(...)))
}

#' @title Helper to apply functions
#' 
#' @description 
#' Used to make use of apply-like operations, regardless of wether the input 
#' is a matrix or a data.frame
apply_array <- function(obj, dim, fun) {
    if (is.matrix(obj)) {
        return(apply(obj, dim, fun))
    }
    if (is.data.frame(obj)) {
        if (dim == 2) {
            return(sapply(obj, fun))            
        } else {
            return(apply(obj, dim, fun))
        }
    }
    
    NULL
}