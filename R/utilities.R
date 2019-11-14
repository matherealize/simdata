# Templates for datamatrix rejection functions ######################
#' @title Check if matrix contains constant column(s)
#'
#' @param x
#' Matrix or Data.frame.
#'
#' @return
#' TRUE if one of the columns has standard deviation of 0, else FALSE.
#'
#' @note
#' Prints a warning if constant is found.
#'
#' @export
contains_constant <- function(x) {
    x = as.matrix(x)
    val = any(apply(x, 2, sd) == 0)
    
    if (val)
        warning("Matrix contains constant column.\n")
    
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
        warning("Matrix is not full rank.\n")
    
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
#' \code{\link{cor_to_upper}}
#'
#' @examples
#' cor_from_upper(2, rbind(c(1,2,0.8)))
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
#' \code{cor_from_upper}
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
#' \code{\link{cor_from_upper}}
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

# TODO: doc both
cov_to_cor <- function(m) {
    cov2cor(m)
}

cor_to_cov <- function(m, sds = NULL) {
    if (is.null(sds))
        sds = rep(1, nrow(m))
    
    diag(sds) %*% m %*% diag(sds)
}

# Functions to work with distributions ##############################
# TODO: document assumption to look like r function for rng
draw_from_distribution <- function(dist_fct, n, ...) {
    do.call(dist_fct,
            modifyList(list(n), list(...))
    )
}

# Transformation functions ##########################################
# TODO DOC 
# make aware that its best to name the arguments
# more general than just transforming column by column...
function_list <- function(...) {
    fct_list = list(...)
    
    function(m) {
        do.call(data.frame, lapply(fct_list, function(f, x) f(x), m))
    }
}