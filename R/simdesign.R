#' @title Design specification for simulating datasets
#'
#' @description
#' Stores information necessary to simulate and visualize datasets based
#' on underlying distribution `Z`.
#'
#' @template simulate_data_template
#' @template simdesign_template
#' @param ...
#' Further arguments are directly stored in the list object to be passed to
#' \code{\link{simulate_data}}.
#'
#' @details
#' The `simdesign` class should be used in the following workflow:
#'
#' \enumerate{
#' \item Specify a design template which will be used in subsequent data
#' generating / visualization steps.
#' \item Sample / visualize datamatrix following template (possibly
#'  multiple times) using \code{\link{simulate_data}}.
#' \item Use sampled datamatrix for simulation study.
#' }
#'
#' For more details on generators and transformations, please see the
#' documentation of \code{\link{simulate_data}}.
#'
#' For details on post-processing, please see the documentation of
#' \code{\link{do_processing}}.
#'
#' @section Naming of variables:
#' If `check_and_infer` is set to TRUE, the following procedure determines
#' the names of the variables:
#'
#' \enumerate{
#' \item use `names_final` if specified and of correct length
#' \item otherwise, use the names of `transform_initial` if present and of
#' correct length
#' \item otherwise, use `prefix_final` to prefix the variable number if
#' not NULL
#' \item otherwise, use names from dataset as generated by the `generator`
#' function
#' }
#'
#' @section Simulation Templates:
#' This class is intended to be used as a template for simulation designs
#' which are based on specific underlying distributions. All such a template
#' needs to define is the `generator` function and its construction and
#' pass it to this function along with the other arguments. See
#' \code{\link{simdesign_mvtnorm}} for an example.
#'
#' @return
#' List object with class attribute "simdesign" (S3 class) containing
#' the following entries (if no further information given, entries are
#' directly saved from user input):
#'
#' \describe{
#' \item{`generator`}{}
#' \item{`name`}{}
#' \item{`transform_initial`}{}
#' \item{`n_var_final`}{}
#' \item{`types_final`}{}
#' \item{`names_final`}{}
#' \item{`process_final`}{}
#' \item{`entries for further information as passed by the user`}{}
#' }
#'
#' @examples
#' generator <- function(n) mvtnorm::rmvnorm(n, mean = 0)
#' sim_design <- simdesign(generator)
#' simulate_data(sim_design, 10, seed = 19)
#'
#' @seealso
#' \code{\link{simdesign_mvtnorm}},
#' \code{\link{simulate_data}},
#' \code{\link{simulate_data_conditional}}
#'
#' @importFrom utils modifyList
#'
#' @export
simdesign <- function(generator,
                      transform_initial = base::identity,
                      n_var_final = -1,
                      types_final = NULL,
                      names_final = NULL,
                      prefix_final = "v",
                      process_final = list(),
                      name = "Simulation design",
                      check_and_infer = TRUE,
                      ...) {
    design <- c(
        list(
            generator = generator,
            name = name,
            transform_initial = transform_initial,
            n_var_final = n_var_final,
            types_final = types_final,
            names_final = names_final,
            process_final = process_final
        ),
        list(...)
    )
    class(design) <- "simdesign"

    if (!check_and_infer) {
        return(design)
    }

    # check if simulation design works by simulating 5 samples
    # do not use names for that to catch errors later on
    res <- tryCatch(
        simulate_data(modifyList(design, list(names_final = NULL)), 5),
        error = function(err) err
    )

    if (inherits(res, "condition")) {
        warning(
            "Unable to simulate from design, returning error: ",
            sprintf("%s", res),
            "Please double check arguments. ",
            "Returning the potentially faulty design object."
        )
        return(design)
    }

    if (is.null(res)) {
        warning(
            "Simulation from design returned NULL. ",
            "Please double check arguments. ",
            "Returning the potentially faulty design object."
        )
        return(design)
    }

    if (length(dim(res)) != 2) {
        warning(
            "Simulation from design did not return a 2-dimensional array ",
            "(matrix or data.frame).",
            "Please double check arguments. ",
            "Returning the potentially faulty design object."
        )
        return(design)
    }

    if (dim(res)[1] != 5) {
        warning(
            "Simulation from design returned unexpected number ",
            sprintf("of observations (%d instead of %d). ", length(res), 5),
            "Please double check arguments. ",
            "Returning the potentially faulty design object."
        )
        return(design)
    }

    n_var_res <- ncol(res)
    names_res <- names(res)

    if (is.null(design$n_var_final) | n_var_final != n_var_res) {
        design$n_var_final <- n_var_res
    }

    # infer data types if necessary
    if (is.null(design$types_final) | length(design$types_final) != n_var_res) {
        design$types_final <- apply_array(res, 2, class)
    }

    # infer names if necessary and take in the following order
    # 1) from names final
    # 2) from transform_initial (could be partial!)
    # 3) from prefix_final (to fill in missing names)

    if (!is.null(names_final) &
        n_var_res != length(design$names_final)) {
        warning(
            "Number of simulated variables differs from length of ",
            "'names_final'. Resetting names in output."
        )
    }

    if (length(names_final) == n_var_res) {
        design$names_final <- names_final
    } else if (length(
        get_names_from_function_list(transform_initial)
    ) == n_var_res) {
        design$names_final <- get_names_from_function_list(transform_initial)
    } else if (!is.null(prefix_final)) {
        design$names_final <- paste0(prefix_final, 1:n_var_res)
    } else {
        # note that the colnames could be empty, therefore
        # we need to modifyList
        design <- utils::modifyList(design,
            list(names_final = colnames(res)),
            keep.null = TRUE
        )
    }

    design
}

# Design templates ############################################################

#' @title NORTA-based design specification
#'
#' @description
#' Stores information necessary to simulate datasets based on the NORTA
#' procedure (Cario and Nelson 1997).
#'
#'
#' @param cor_target_final
#' Target correlation matrix for simulated datasets. At least one of
#' `cor_target_final` or `cor_initial` must be specified.
#' @param cor_initial
#' Correlation matrix for underlying multivariate standard normal distribution
#' on which the final data is based on.  At least one of `cor_target_final` or
#' `cor_initial` must be specified. If NULL, then `cor_initial` will be
#' numerically optimized by simulation for the NORTA procedure using
#' `cor_target_final`.
#' @param dist
#' List of functions of marginal distributions for simulated variables.
#' Must have the same length as the specified correlation matrix
#' (`cor_target_final` and / or `cor_inital`), and the order of the entries
#' must correspond to the variables in the correlation matrix. See details for
#' the specification of the marginal distributions.
#' @param tol_initial
#' If `cor_initial` is numerically optimized, specifies the tolerance for the
#' difference to the target correlation `cor_target_final`. Parameter passed to
#' \code{\link{optimize_cor_for_pair}}.
#' @param n_obs_initial
#' If `cor_initial` is numerically optimized, specifies the number of draws in
#' simulation during optimization used to estimate correlations.
#' Parameter passed to \code{\link{optimize_cor_for_pair}}.
#' @param seed_initial
#' Seed used for draws of the initial distribution used during optimization
#' to estimate correlations.
#' @param conv_norm_type
#' If `cor_initial` is numerically optimized and found not to be a proper
#' correlation matrix (i.e. not positive-definite), specifies the metric used to
#' find the nearest positive-definite correlation matrix.
#' Parameter passed to \code{\link[Matrix:nearPD]{Matrix::nearPD}}
#' (conv.norm.type), see there for details.
#' @param method
#' `method` argument of \code{\link[mvtnorm:Mvnorm]{mvtnorm::rmvnorm}}.
#' @param name
#' Character, optional name of the simulation design.
#' @param ...
#' Further arguments are passed to the \code{\link{simdesign}} constructor.
#'
#' @details
#' This S3 class implements a simulation design based on the
#' NORmal-To-Anything (NORTA) procedure by Cario and Nelson (1997). See the
#' corresponding NORTA vignette for usage examples how to approximate real
#' datasets.
#'
#' @section Data Generation:
#' Data will be generated using the following procedure:
#' \enumerate{
#' \item An underlying data matrix `Z` is sampled from a
#' multivariate standard Normal distribution with correlation structure given by
#' `cor_initial`.
#' \item `Z` is then transformed into a dataset `X` by applying
#' the functions given in `dist` to the columns of `Z`. The resulting dataset
#' `X` will then have the desired marginal distributions, and approximate the
#' target correlation `cor_target_final`, if specified.
#' \item `X` is further transformed by the transformation `transform_initial`
#' (note that this may affect the correlation of the final dataset and is not
#' respected by the optimization procedure), and post-processed if specified.
#' }
#'
#' @section Marginal distributions:
#' A list of functions `dist` is used to define the marginal distributions of
#' the variables. Each entry must be a quantile function, i.e. a function
#' that maps `[0, 1]` to the domain of a probability distribution. Each entry
#' must take a single input vector, and return a single numeric vector.
#' Examples for acceptable entries include all standard quantile functions
#' implemented in R (e.g. `qnorm`, `qbinom`, ...), user defined functions
#' wrapping these (e.g. `function(x) = qnorm(x, mean = 10, sd = 4)`), or
#' empirical quantile functions. The helper function 
#' \code{}\link{quantile_functions_from_data} can be used to automatically 
#' estimate empirical quantile functions from a given data to reproduce it using
#' the NORTA approach.See the example in the NORTA vignette of this package for
#' workflow details.
#'
#' @section Target correlations:
#' Not every valid correlation matrix (i.e. symmetric, positive-definite matrix
#' with elements in `[-1, 1]` and unity diagonal) for a number of variables
#' is feasible for given desired marginal distributions (see e.g.
#' Ghosh and Henderson 2003). Therefore, if `cor_target_final` is specified
#' as target correlation, this class optimises `cor_initial` in such a
#' way, that the final simulated dataset has a correlation which approximates
#' `cor_target_final`. However, the actual correlation in the end may differ
#' if `cor_target_final` is infeasible for the given specification, or the
#' NORTA procedure cannot exactly reproduce the target correlation. In general,
#' however, approximations should be acceptable if target correlations and
#' marginal structures are derived from real datasets.
#' See e.g. Ghosh and Henderson 2003 for the motivation why this works.
#'
#' @return
#' List object with class attribute "simdesign_norta" (S3 class), inheriting
#' from "simdesign". It contains the same entries as a \code{\link{simdesign}}
#' object but in addition the following entries:
#'
#' \describe{
#' \item{`cor_target_final`}{}
#' \item{`cor_initial`}{Initial correlation matrix of multivariate normal
#' distribution}
#' \item{`dist`}{}
#' \item{`tol_initial`}{}
#' \item{`n_obs_initial`}{}
#' \item{`conv_norm_type`}{}
#' \item{`method`}{}
#' }
#'
#' @seealso
#' \code{\link{simdesign}},
#' \code{\link{simulate_data}},
#' \code{\link{simulate_data_conditional}}, 
#' \code{\link{quantile_functions_from_data}}
#'
#' @references Cario, M. C. and Nelson, B. L. (1997) \emph{Modeling and
#' generating random vectors with arbitrary marginal distributions and
#' correlation matrix}. Technical Report, Department of Industrial Engineering
#' and Management Sciences, Northwestern University, Evanston, Illinois.
#'
#' Ghosh, S. and Henderson, S. G. (2003) \emph{Behavior of the NORTA method
#' for correlated random vector generation as the dimension increases}. ACM
#' Transactions on Modeling and Computer Simulation.
#'
#' @importFrom stats pnorm
#'
#' @export
simdesign_norta <- function(cor_target_final = NULL,
                            cor_initial = NULL,
                            dist = list(),
                            tol_initial = 0.001,
                            n_obs_initial = 10000,
                            seed_initial = 1,
                            conv_norm_type = "O",
                            method = "svd",
                            name = "NORTA based simulation design",
                            ...) {
    if (is.null(cor_target_final) & is.null(cor_initial)) {
        stop("One of 'cor_target_final' or 'cor_initial' must be specified.")
    }

    if (is.null(cor_initial) & !is.null(cor_target_final)) {
        if (!is_cor_matrix(cor_target_final)) {
            stop("'cor_target_final' must be a proper correlation matrix.")
        }

        # optimize initial correlation structure
        cor_initial <- optimize_cor_mat(cor_target_final, dist,
            ensure_cor_mat = TRUE,
            return_diagnostics = FALSE,
            conv_norm_type = conv_norm_type,
            tol = tol_initial,
            n_obs = n_obs_initial,
            seed = seed_initial
        )
    }

    if (!is_cor_matrix(cor_initial)) {
        stop("'cor_initial' must be a proper correlation matrix.")
    }

    # generator
    generator <- function(n) {
        W <- mvtnorm::rmvnorm(n,
            mean = rep(0, ncol(cor_initial)),
            sigma = cor_initial,
            method = method
        )
        colapply_functions(pnorm(W), dist)
    }

    # setup simulation design
    dsgn <- simdesign(
        generator = generator,
        name = name,
        cor_target_final = cor_target_final,
        cor_initial = cor_initial,
        tol_initial = tol_initial,
        n_obs_initial = n_obs_initial,
        conv_norm_type = conv_norm_type,
        method = method,
        ...
    )

    class(dsgn) <- c("simdesign_norta", class(dsgn))

    dsgn
}

#' @title Multivariate normal design specification
#'
#' @description
#' Stores information necessary to simulate and visualize datasets based
#' on underlying distribution multivariate normal distribution `Z`.
#'
#' @param relations_initial
#' Correlation / Covariance matrix of the initial multivariate
#' Normal distribution `Z`.
#' @param mean_initial
#' Vector of mean values of the initial multivariate Normal
#' distribution `Z`. Dimension needs to correspond to dimension
#' of `relations`.
#' @param sd_initial
#' Vector of standard deviations of the initial multivariate
#' Normal distribution Z. Dimension needs to correspond to dimension
#' of `relations`. Overriden by suqare root of diagonal elements of
#' `relations` if `is_correlation` is FALSE.
#' @param is_correlation
#' If TRUE, then `relations` specifies a correlation matrix (default,
#' this type of specification is usually more natural than specifying
#' a covariance matrix). Otherwise, `relations` specifies a
#' covariance matrix whose square root diagonal elements override
#' `sd_initial`.
#' @param method
#' `method` argument of \code{\link[mvtnorm:Mvnorm]{mvtnorm::rmvnorm}}.
#' @param name
#' Character, optional name of the simulation design.
#' @param ...
#' Further arguments are passed to the \code{\link{simdesign}} constructor.
#'
#' @details
#' This S3 class implements a simulation design based on an underlying
#' multivariate normal distribution by creating a `generator` function
#' based on \code{\link[mvtnorm:Mvnorm]{mvtnorm::rmvnorm}}.
#'
#' @section Data Generation:
#' Data will be generated by \code{\link{simulate_data}} using the
#' following procedure:
#' \enumerate{
#' \item The underlying data matrix `Z` is sampled from a
#' multivariate Normal distribution (number of dimensions specified by
#' dimensions of `relations`).
#' \item `Z` is then transformed into the final dataset `X` by applying
#' the `transform_initial` function to `Z`.
#' \item `X` is post-processed if specified.
#' }
#'
#' @note
#' Note that `relations` specifies the correlation / covariance
#' of the underlying Normal data `Z` and thus does not directly translate into
#' correlations between the variables of the final datamatrix `X`.
#'
#' @return
#' List object with class attribute "simdesign_mvtnorm" (S3 class), inheriting
#' from "simdesign". It contains the same entries as a \code{\link{simdesign}}
#' object but in addition the following entries:
#'
#' \describe{
#' \item{`mean_initial`}{}
#' \item{`sd_initial`}{}
#' \item{`cor_initial`}{Initial correlation matrix of multivariate normal
#' distribution}
#' }
#'
#' @seealso
#' \code{\link{simdesign}},
#' \code{\link{simulate_data}},
#' \code{\link{simulate_data_conditional}},
#' \code{\link{plot_cor_network.simdesign_mvtnorm}}
#'
#' @export
simdesign_mvtnorm <- function(relations_initial,
                              mean_initial = 0,
                              sd_initial = 1,
                              is_correlation = TRUE,
                              method = "svd",
                              name = "Multivariate-normal based simulation design",
                              ...) {
    # prepare means and standard deviations
    mean_initial <- rep_len(mean_initial, nrow(relations_initial))
    sd_initial <- rep_len(sd_initial, nrow(relations_initial))

    # convert correlation to covariance
    cor_initial <- relations_initial
    if (is_correlation) {
        relations_initial <- cor_to_cov(relations_initial, sd_initial)
    } else {
        # is covariance matrix
        sd_initial <- sqrt(diag(relations_initial))
        cor_initial <- stats::cov2cor(cor_initial)
    }

    # define generator
    generator <- function(n) {
        mvtnorm::rmvnorm(n,
            mean = mean_initial,
            sigma = relations_initial,
            method = method
        )
    }

    # setup simulation design
    dsgn <- simdesign(
        generator = generator,
        mean_initial = mean_initial,
        sd_initial = sd_initial,
        cor_initial = cor_initial,
        name = name,
        ...
    )

    class(dsgn) <- c("simdesign_mvtnorm", class(dsgn))

    dsgn
}

#' @title Uniform disc sampling design specification
#'
#' @description
#' Provides 2-dimensional points, spread uniformly over disc, or partial
#' disc segment (i.e. a circle, or ring, or ring segment). Useful for e.g.
#' building up clustering exercises.
#'
#' @param r_min
#' Minimum radius of points.
#' @param r_max
#' Maximum radius of points.
#' @param angle_min
#' Minimum angle of points (between 0 and 2pi).
#' @param angle_max
#' Maximum angle of points (between 0 and 2pi).
#' @param name
#' Character, optional name of the simulation design.
#' @param ...
#' Further arguments are passed to the \code{\link{simdesign}} constructor.
#'
#' @details
#' The distribution of points on a disk depends on the radius - the farther out,
#' the more area the points need to cover. Thus, simply sampling two uniform
#' values for radius and angle will not work. See references.
#' 
#' @return
#' List object with class attribute "simdesign_discunif" (S3 class), inheriting
#' from "simdesign". It contains the same entries as a \code{\link{simdesign}}
#' object but in addition the following entries:
#'
#' \describe{
#' \item{`r_min`}{}
#' \item{`r_max`}{}
#' \item{`angle_min`}{}
#' \item{`angle_max`}{}
#' }
#'
#' @references
#' \url{https://mathworld.wolfram.com/DiskPointPicking.html}
#'
#' @examples
#' disc_sampler <- simdesign_discunif()
#' plot(simulate_data(disc_sampler, 1000, seed = 19))
#'
#' ring_segment_sampler <- simdesign_discunif(r_min = 0.5, angle_min = 0.5 * pi)
#' plot(simulate_data(ring_segment_sampler, 1000, seed = 19))
#'
#' circle_sampler <- simdesign_discunif(r_min = 1)
#' plot(simulate_data(circle_sampler, 1000, seed = 19))
#'
#' @importFrom stats runif
#'
#' @export
simdesign_discunif <- function(r_min = 0, r_max = 1,
                               angle_min = 0, angle_max = 2 * pi,
                               name = "Uniform circle simulation design",
                               ...) {
    if (r_min < 0 | r_max < 0) {
        stop("r_min and r_max must be positive.")
    }
    if (r_max < r_min) {
        stop("r_max must be larger than r_min.")
    }
    if (angle_max < angle_min) {
        stop("angle_max must be larger than angle_min.")
    }

    # define generator
    # use sqrt of radius instead of radius itself
    generator <- function(n) {
        r <- sqrt(runif(n, min = r_min, max = r_max))
        phi <- runif(n, min = angle_min, max = angle_max)
        matrix(c(r * cos(phi), r * sin(phi)), ncol = 2)
    }

    # setup simulation design
    dsgn <- simdesign(
        generator = generator,
        r_min = r_min,
        r_max = r_max,
        angle_min = angle_min,
        angle_max = angle_max,
        name = name,
        ...
    )

    class(dsgn) <- c("simdesign_discunif", class(dsgn))

    dsgn
}
