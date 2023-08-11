# Templates for datamatrix rejection functions ######################
#' @title Check if matrix contains constant column(s)
#'
#' @param x
#' Matrix or Data.frame.
#' @param eps
#' Threshold for standard deviation below which a column is considered to be
#' constant.
#'
#' @return
#' TRUE if one of the columns has standard deviation of below `eps``, else
#' FALSE.
#'
#' @note
#' Prints a warning if constant is found.
#'
#' @export
contains_constant <- function(x, eps = .Machine$double.eps) {
    x <- as.matrix(x)
    val <- any(apply(x, 2, stats::sd) < eps)

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
    x <- as.matrix(x)
    val <- qr(x)$rank < ncol(x)

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
#' \code{\link{cor_to_upper}}
#'
#' @examples
#' \dontrun{
#' cor_from_upper(2, rbind(c(1, 2, 0.8)))
#' }
#'
#' @export
cor_from_upper <- function(n_var, entries = NULL) {
    if (any(entries[, 1:2] < 0))
        stop("Indices in 'entries[, 1:2]' must be positive.")
    if (any(entries[, 3] < -1) | any(entries[, 3] > 1))
        stop("All correlations in 'entries[, 3]' must be in [-1,1].")

    # ones in the diagonal
    cor_mat <- diag(n_var)

    if (!is.null(entries)) {
        # set correlation entries
        # ensure the conversion to a matrix works for vectors and data.frames
        entries <- matrix(as.matrix(entries), ncol = 3)

        # ensure that entries are in upper triangular part
        entries[, 1:2] <- cbind(
            pmin.int(entries[, 1], entries[, 2]),
            pmax.int(entries[, 1], entries[, 2])
        )
        cor_mat[entries[, 1:2, drop = FALSE]] <- entries[, 3]

        # symmetrize
        cor_mat[lower.tri(cor_mat)] <- t(cor_mat)[lower.tri(cor_mat)]
    }

    # make sure diagonal elements are 1
    diag(cor_mat) <- 1

    cor_mat
}

#' @title Convert correlation matrix to specification used by
#' `cor_from_upper`
#'
#' @param m
#' Symmetric correlation matrix.
#' @param remove_below
#' Threshold for absolute correlation values below which they are removed from
#' the returned matrix. If NULL then no filtering is applied.
#'
#' @return
#' Matrix with 3 columns (variable_1, variable_2, correlation), where
#' correlation gives the entry at position (variable_1, variable_2) of the
#' input correlation matrix. Note that variable_1 < variable_2 holds for all
#' entries.
#'
#' @seealso
#' \code{\link{cor_from_upper}}
#'
#' @export
cor_to_upper <- function(m, remove_below = .Machine$double.eps) {
    m <- as.matrix(m)
    if (!is_cor_matrix(m))
        stop("'m' is not a proper correlation matrix.")

    ind <- which(upper.tri(m), arr.ind = TRUE)
    res <- cbind(ind, m[ind])

    if (!is.null(remove_below))
        res <- res[abs(res[, 3]) >= remove_below, , drop = FALSE]

    res
}

#' @title Convert correlation matrix to covariance matrix
#'
#' @description
#' Rescale correlation matrix by variable standard deviations to yield a
#' covariance matrix.
#'
#' @param m
#' Symmetric correlation matrix.
#' @param sds
#' Standard deviations of the variables. Set to 1 for all varirables by default.
#'
#' @return
#' Symmetric covariance matrix.
#'
#' @export
cor_to_cov <- function(m, sds = NULL) {
    if (is.null(sds))
        sds <- rep(1, nrow(m))

    if (any(sds < 0))
        warning("All entries in 'sds' must be positive.")

    m <- as.matrix(m)
    if (!is_cor_matrix(m))
        warning("'m' is not a proper correlation matrix.")

    diag(sds) %*% m %*% diag(sds)
}

#' @title Check if matrix is a correlation matrix
#'
#' @description
#' Checks if matrix is numeric, symmetric, has diagonal elements of one,
#' has only entries in `[-1, 1]`, and is positive definite. Prints a warning
#' if a problem was found.
#'
#' @param m
#' Matrix.
#'
#' @return
#' TRUE if matrix is a correlation matrix, else FALSE.
#'
#' @export
is_cor_matrix <- function(m, tol=10e-10) {
    ok <- TRUE
    if (!is.numeric(m)) {
        warning("'m' must be numeric.")
        ok <- FALSE
    }
    if (!isSymmetric(m)) {
        warning("'m' must be a symmetric matrix.")
        ok <- FALSE
    }
    if (any(abs(diag(m) - 1) > tol)) {
        warning("The diagonal elements of 'm' must be equal to one.")
        ok <- FALSE
    }
    if (any(m > 1) | any(m < -1)) {
        warning("All the elements of 'm' must be in [-1,1].")
        ok <- FALSE
    }
    m_ev <- eigen(m, symmetric = TRUE, only.values = TRUE)
    if (any(m_ev$values <= 0)) {
        warning("'m' must be a positive definite matrix.")
        ok <- FALSE
    }

    ok
}

#' @title Find pairwise initial correlation for NORTA from target correlation
#'
#' @description
#' This function can be used to find a suitable initial correlation for use
#' in the NORTA procedure for a pair of variables with given marginal
#' distributions and target correlation.
#'
#' @param cor_target
#' Target correlation of variable pair.
#' @param dist1,dist2
#' Marginal distributions of variable pair, given as univariable quantile
#' functions.
#' @param n_obs
#' Number of observations to be used in the numerical optimization procedure.
#' @param seed
#' Seed for generating standard normal random variables in the numerical
#' optimization procedure.
#' @param tol,...
#' Further parameters passed to \code{\link[stats:uniroot]{stats::uniroot}}.
#'
#' @details
#' Uses \code{\link[stats:uniroot]{stats::uniroot}} for actual optimization.
#'
#' @return
#' Output of \code{\link[stats:uniroot]{stats::uniroot}} for the univariable
#' optimization for find the initial correlation.
#'
#' @importFrom stats cor pnorm uniroot
optimize_cor_for_pair <- function(cor_target, dist1, dist2,
                                  n_obs = 100000, seed = NULL,
                                  tol = 0.01, ...) {

    if (cor_target > 1 | cor_target < -1)
        stop("'cor_target' must be a correlation in [-1,1].")

    objective <- function(cor_current, cor_target, n_obs,
                          dist1, dist2) {
        cor_mat <- diag(2)
        cor_mat[1, 2] <- cor_mat[2, 1] <- cor_current
        z <- pnorm(mvtnorm::rmvnorm(n_obs, mean = c(0, 0),
            sigma = cor_mat,
            method = "svd"))
        x <- cbind(dist1(z[, 1]), dist2(z[, 2]))
        cor(x)[1, 2] - cor_target
    }

    if (!is.null(seed))
        set.seed(seed)

    uniroot(objective, interval = c(-1, 1), tol = tol,
        cor_target = cor_target, n_obs = n_obs,
        dist1 = dist1, dist2 = dist2, ...)
}

#' @title Find initial correlation matrix for NORTA from target correlation
#'
#' @description
#' This function can be used to find a suitable correlation matrix to be used
#' for simulating initial multivariate normal data in a NORTA based simulation
#' design (see \code{\link{simdesign_norta}}).
#'
#' @param cor_target
#' Target correlation matrix.
#' @param dist
#' List of functions of marginal distributions for simulated variables.
#' Must have the same length as the specified correlation matrix
#' (`cor_target`), and the order of the entries must correspond to the
#' variables in the correlation matrix. See \code{\link{simdesign_norta}} for
#' details of the specification of the marginal distributions.
#' @param ensure_cor_mat
#' if TRUE, this function ensures that the optimized matrix is a proper
#' correlation matrix by ensuring positive definitiness. If FALSE, the
#' optimized matrix is returned as is.
#' @param conv_norm_type
#' Metric to be used to find closest positive definite matrix to optimal matrix,
#' used if `ensure_cor_mat` is TRUE.
#' Passed to \code{\link[Matrix:nearPD]{Matrix::nearPD}}.
#' @param return_diagnostics
#' TRUE to return additional diagnostics of the optimization procedure, see
#' below.
#' @param ...
#' Additional parameters passed to \code{\link{optimize_cor_for_pair}}.
#'
#' @details
#' This function first finds a suitable correlation matrix for the underlying
#' multivariate normal data used in the NORTA procedure. It does so by
#' solving k*(k-1) univariable optimisation problems (where k is the number
#' of variables). In case the result is not a positive-definite matrix, the
#' nearest positive-definite matrix is found according to the user specified
#' metric using \code{\link[Matrix:nearPD]{Matrix::nearPD}}.
#' See e.g. Ghosh and Henderson (2003) for an overview of the procedure.
#'
#' @return
#' If `return_diagnostics` is FALSE, a correlation matrix to be used in the
#' definition of a \code{\link{simdesign_norta}} object. If TRUE, then a list
#' with two entries: `cor_mat` containing the correlation matrix, and
#' `convergence` containing a list of objects returned by the individual
#' optimisation problems from \code{\link[stats:uniroot]{stats::uniroot}}.
#'
#' @references Ghosh, S. and Henderson, S. G. (2003) \emph{Behavior of the
#' NORTA method for correlated random vector generation as the dimension
#' increases}. ACM Transactions on Modeling and Computer Simulation.
#'
#' @seealso
#' \code{\link{simdesign_norta}}
#'
#' @export
optimize_cor_mat <- function(cor_target, dist,
                             ensure_cor_mat = TRUE,
                             conv_norm_type = "O",
                             return_diagnostics = FALSE, ...) {

    if (!is_cor_matrix(cor_target)) {
        stop("'cor_target' must be a proper correlation matrix.")
    }
    if (length(dist) != ncol(cor_target)) {
        stop("Number of marginal distributions must be equal to size of ",
            "correlation matrix.")
    }

    cor_mat <- diag(ncol(cor_target))
    conv_res <- list()

    # TODO
    # if cor_target is 0, can the underlying cor be anything different than 0?
    # to speed up adding uncorrelated vars
    # maybe specific tol below which it is considered 0

    for (row in seq(1, nrow(cor_target) - 1)) {
        conv_res[[row]] <- list()
        for (col in seq(row + 1, ncol(cor_target))) {
            conv_res[[row]][[col]] <- optimize_cor_for_pair(
                cor_target = cor_target[row, col],
                dist1 = dist[[row]], dist2 = dist[[col]], ...
            )
            cor_mat[row, col] <- conv_res[[row]][[col]]$root
        }
    }

    cor_mat[lower.tri(cor_mat)] <- t(cor_mat)[lower.tri(cor_mat)]

    if (ensure_cor_mat) {
        # ensure that matrix found is positive definite correlation matrix
        # by finding the closes pd matrix
        respd <- Matrix::nearPD(cor_mat, corr = TRUE, keepDiag = TRUE,
            conv.norm.type = conv_norm_type, trace = FALSE)
        cor_mat <- as.matrix(respd$mat)
    }

    if (return_diagnostics) {
        return(list(cor_mat = cor_mat, convergence = conv_res))
    } else {
        return(cor_mat)
    }
}

# Function helper ###################################################
#' @title Apply list of functions to input
#'
#' @param ...
#' Named or unnamed arguments, each of which is a function taking exactly
#' one input. See details.
#' @param stringsAsFactors,check.names
#' Arguments of \code{\link[base:data.frame]{data.frame}}.
#'
#' @details
#' This is a convenience function which takes a number of functions and returns
#' another function which applies all of the user specified functions to a new
#' input, and collects the results as list or data.frame.
#' This is useful to e.g. transform columns of a data.frame or check
#' the validity of a matrix during simulations. See the example here and
#' in \code{\link{simulate_data_conditional}}.
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
#' f <- function_list(
#'     v1 = function(x) x[, 1] * 2,
#'     v2 = function(x) x[, 2] + 10)
#'
#' f(diag(2))
#'
#' # function_list can be used to add new columns
#' # naming of columns should be handled separately in such cases
#'
#' f <- function_list(
#'     function(x) x, # return x as it is
#'     X1_X2 = function(x) x[, 2] + 10) # add new column
#'
#' f(diag(2))
#' }
#'
#' @return
#' Function with a single input which outputs a data.frame. Has special
#' 'flist' entry in its environment which stores individual functions as list.
#'
#' @seealso
#' \code{\link[base:data.frame]{data.frame}},
#' \code{\link{get_from_function_list}},
#' \code{\link{get_names_from_function_list}}
#'
#' @export
function_list <- function(...,
                          stringsAsFactors = FALSE,
                          check.names = TRUE) {
    fct_list <- list(...)

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
#' `function_list` or function object.
#'
#' @return
#' List with named or unnamed entries corresponding to individual function
#' objects that were passed to the `function_list` object. If `flist` is a
#' simple function, returns NULL.
#'
#' @seealso
#' \code{\link{function_list}}
#'
#' @export
get_from_function_list <- function(flist) {
    get0("fct_list", envir = environment(flist))
}

#' @title Extract names of individual functions from `function_list`
#'
#' @description
#' Extract names of  individual function objects from environment of a
#' `function_list` object.
#'
#' @param flist
#' `function_list` or function object.
#'
#' @return
#' Names of list corresponding to individual function objects that were passed
#' to the `function_list` object. If `flist` is a simple function, returns NULL.
#'
#' @seealso
#' \code{\link{function_list}}
#'
#' @export
get_names_from_function_list <- function(flist) {
    names(get_from_function_list(flist))
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
#' Passed to \code{\link{function_list}}.
#'
#' @inherit function_list return
#'
#' @seealso
#' \code{\link{function_list}}
#'
#' @export
as_function_list <- function(flist, ...) {
    do.call(function_list, append(flist, list(...)))
}

#' @title Helper to apply functions
#'
#' @description
#' Used to make use of apply-like operations, regardless of wether the input
#' is a matrix or a data.frame
#'
#' @param obj
#' Matrix or data.frame.
#' @param dim
#' Dimension to apply function to.
#' @param fun
#' Function object to apply.
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

#' @title Apply list of functions to column of object
#'
#' @description
#' Helper function to simplify workflow with lists of functions.
#'
#' @param obj
#' 2-dimensional array (matrix or data.frame).
#' @param flist
#' List of functions of length equal to the number of columns of `obj`.
#' Each entry must be a function applicable to a single column of `obj`.
#' The i-th entry of `flist` is applied to the i-th column of `obj`.
#'
#' @return
#' Matrix or data.frame (same type as `obj`) with names taken from `obj`.
colapply_functions <- function(obj, flist) {
    if (length(dim(obj)) != 2)
        stop("'obj' must be a 2-dimensional array.")

    if (length(flist) != dim(obj)[2])
        stop("Number of columns of 'obj' must be equal to number of functions.")

    res <- lapply(1:length(flist), function(col) flist[[col]](obj[, col]))
    names(res) <- colnames(obj)
    res <- do.call(cbind, res)
    if (is.data.frame(obj))
        res <- data.frame(res)

    res
}

#' @title Define partial function
#' 
#' @description 
#' Partial functions are useful to define marginal distributions based on 
#' additional parameters. 
#' 
#' @param f 
#' Function in two or more parameters.
#' @param ... 
#' Parameters to be held fixed for function `f`.
#' 
#' @details 
#' This helper function stores passed arguments in a list, and stores this 
#' list in the environment of the returned function. Thus, it remembers the
#' arguments that should be held fixed, such that the returned partial function
#' now is a function with fewer arguments.
#' 
#' @return 
#' Function object.
#' 
#' @examples
#' marginal <- partial(function(x, meanx) qnorm(x, meanx), meanx = 2)
#' marginal(0.5)
#' 
#' @export
partial <- function(f, ...) {
    f_args <- list(...)
    
    function(...) {
        do.call(f, c(f_args, list(...)))
    }
}