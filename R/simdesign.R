# TODO: Tests

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
#' \code{\link{process_data}}.
#' 
#' @section Simulation Templates:
#' This class is intended to be used as a template for simulation designs 
#' which are based on specific underlying distributions. All such a template
#' needs to define is the `generator` function and its construction and 
#' pass it to this function along with the other arguments. See 
#' \code{\link{mvtnorm_simdesign}} for an example.
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
#' @seealso
#' \code{\link{mvtnorm_simdesign}},
#' \code{\link{simulate_data}},
#' \code{\link{simulate_data_conditional}}
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
    design = c(
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
    class(design) = "simdesign"
    
    if (n_var_final > 0 & !is.null(prefix_final))
        design$names_final = paste0(prefix_final, 1:n_var_final)
    
    if (!check_and_infer)
        return(design)
    
    # check if simulation design works by simulating 5 samples
    res = tryCatch(
        simulate_data(design, 5), 
        error = function(err) err
    )
    
    if (inherits(res, "condition")) {
        print(res)
        warning(
            "Unable to simulate from design. ",
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
    
    n_var_res = ncol(res)
    
    if (is.null(design$n_var_final) | n_var_final != n_var_res)
        design$n_var_final = n_var_res
    
    # infer data types if necessary
    if (is.null(design$types_final) | length(design$types_final) != n_var_res) {
        design$types_final = apply_array(res, 2, class)
    }
    
    # infer names if necessary
    if (length(design$names_final) != n_var_res) {
        if (!is.null(prefix_final)) {
            design$names_final = paste0(prefix_final, 1:n_var_res)
        } else {
            # note that the colnames could be empty, therefore
            # we need to modifyList
            design = utils::modifyList(design, 
                                       list(names_final = colnames(res)), 
                                       keep.null = TRUE)
        }
    }
    
    design
}

# Design templates ############################################################

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
#' the `transform` function to `Z`.
#' \item `X` is post-processed if specified.
#' }
#' 
#' @note
#' Note that `relations` specifies the correlation / covariance
#' of the underlying Normal data `Z` and thus does not directly translate into
#' correlations between the variables of the final datamatrix `X`.
#'
#' @return
#' List object with class attribute "mvtnorm_simdesign" (S3 class), inheriting
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
#' \code{\link{plot_cor_network.mvtnorm_simdesign}}
#'
#' @export
mvtnorm_simdesign <- function(relations_initial,
                              mean_initial = 0,
                              sd_initial = 1,
                              is_correlation = TRUE, 
                              method = "svd",
                              name = "Multivariate-normal based simulation design",
                              ...) {
    
    # TODO: assertions, i.e. square matrix, pos def matrix
    
    # prepare means and standard deviations
    mean_initial = rep_len(mean_initial, nrow(relations_initial))
    sd_initial = rep_len(sd_initial, nrow(relations_initial))
    
    # convert correlation to covariance
    cor_initial = relations_initial
    if (is_correlation) {
        relations_initial = cor_to_cov(relations_initial, sd_initial)
    } else {
        # is covariance matrix
        sd_initial = sqrt(diag(relations_initial))
        cor_initial = stats::cov2cor(cor_initial)
    }
    
    # define generator
    generator = function(n) mvtnorm::rmvnorm(n, 
                                             mean = mean_initial, 
                                             sigma = relations_initial, 
                                             method = method)
    
    # setup simulation design
    dsgn = simdesign(
        generator = generator, 
        mean_initial = mean_initial, 
        sd_initial = sd_initial, 
        cor_initial = cor_initial,
        name = name,
        ...
    )
    
    class(dsgn) = c("mvtnorm_simdesign", class(dsgn))
    
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
#' 
#' @details 
#' The distribution of points on a disk depends on the radius - the farther out,
#' the more area the points need to cover. Thus, simply sampling two uniform 
#' values for radius and angle will not work. See references.
#' 
#' @references 
#' \url{https://mathworld.wolfram.com/DiskPointPicking.html}
#' 
#' @examples 
#' \dontrun{
#' disc_sampler = discunif_simdesign()
#' plot(simulate_data(disc_sampler, 1000))
#' 
#' ring_segment_sampler = discunif_simdesign(r_min = 0.5, angle_min = 0.5*pi)
#' plot(simulate_data(ring_segment_sampler, 1000))
#' 
#' circle_sampler = discunif_simdesign(r_min = 1)
#' plot(simulate_data(circle_sampler, 1000))
#' } 
#' 
#' @export
discunif_simdesign <- function(r_min = 0, r_max = 1, 
                               angle_min = 0, angle_max = 2*pi, 
                               name = "Uniform circle simulation design",                              
                               ...) {
    
    # define generator
    # use sqrt of radius instead of radius itself
    generator = function(n) {
        r = sqrt(runif(n, min = r_min, max = r_max) )
        phi = runif(n, min = angle_min, max = angle_max)
        matrix(c(r * cos(phi), r * sin(phi)), ncol = 2)
    }
    
    # setup simulation design
    dsgn = simdesign(
        generator = generator, 
        r_min = r_min,
        r_max = r_max, 
        angle_min = angle_min, 
        angle_max = angle_max,
        name = name,
        ...
    )
    
    class(dsgn) = c("discunif_simdesign", class(dsgn))
    
    dsgn
}