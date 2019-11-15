# TODO: Test if names work, ie finalnames == korrekt!

#' @title Design specification for simulating datasets
#'
#' @description
#' Stores information necessary to simulate and visualize datasets based
#' on underlying distribution *Z*.
#' 
#' @template simulate_data_template
#' @template simdesign_template
#' @param ...
#' Further arguments are directly stored in the list object to be passed to 
#' \code{\link{simulate_data}}.
#'
#' @details
#' The \code{simdesign} class should be used in the following workflow:
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
#' needs to define is the \code{generator} function and its construction and 
#' pass it to this function along with the other arguments. See 
#' \code{\link{mvtnorm_simdesign}} for an example.
#'
#' @return
#' List object with class attribute "simdesign" (S3 class) containing
#' the following entries (if no further information given, entries are
#' directly saved from user input):
#'
#' \describe{
#' \item{\code{generator}}{}
#' \item{\code{name}}{}
#' \item{\code{transform_initial}}{}
#' \item{\code{n_var_final}}{}
#' \item{\code{types_final}}{}
#' \item{\code{names_final}}{}
#' \item{\code{process_final}}{}
#' \item{\code{...}}{Further information as passed by the user.}
#' }
#'
#' @seealso
#' \code{\link{mvtnorm_simdesign}},
#' \code{\link{simulate_data}},
#' \code{\link{conditional_simulate_data}}
#'
#' @export
simdesign <- function(generator,
                      transform_initial = NULL,
                      n_var_final = -1,
                      types_final = NULL,
                      names_final = NULL,
                      prefix_final = "v",
                      process_final = list(), 
                      name = "Simulation design",
                      check_and_infer = TRUE, # TODO document
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
    
    if (is.null(transform_initial)) {
        # default is no transformation at all
        design$transform_initial = base::identity
    }
    
    if (!check_and_infer)
        return(design)
    
    # check if simulation design works by simulating 5 samples
    res = tryCatch(
        simulate_data(design, 5), 
        error = function(err) {
            warning("Unable to simulate from design, please double check arguments. Returning NULL.")
            NULL
        }
    )
    
    if (is.null(res))
        return(NULL)
    
    n_var_res = ncol(res)
    
    if (is.null(design$n_var_final) | n_var_final != n_var_res)
        design$n_var_final = n_var_res
    
    # infer data types if necessary
    if (is.null(design$types_final) | length(design$types_final) != n_var_res) {
        design$types_final = apply(res, 2, class)
    }
    
    # infer names if necessary
    if (length(design$names_final) != n_var_res) {
        if (!is.null(prefix_final)) {
            design$names_final = paste0(prefix_final, 1:n_var_res)
        } else {
            # note that the colnames could be empty, therefore
            # we need to modifyList
            design = modifyList(design, 
                                list(names_final = colnames(res)), 
                                keep.null = TRUE)
        }
    }
    
    design
}

#' @title Design specification for simulating datasets
#'
#' @description
#' Stores information necessary to simulate and visualize datasets based
#' on underlying distribution multivariate normal distribution *Z*.
#' 
#' @param relations_initial
#' Correlation / Covariance matrix of the initial multivariate
#' gaussian distribution *Z*.
#' @param mean_initial
#' Vector of mean values of the initial multivariate gaussian
#' distribution *Z*. Dimension needs to correspond to dimension
#' of \code{relations}.
#' @param sd_initial
#' Vector of standard deviations of the initial multivariate
#' gaussian distribution Z. Dimension needs to correspond to dimension
#' of \code{relations}. Overriden by suqare root of diagonal elements of
#' \code{relations} if \code{is_correlation} is FALSE.
#' @param is_correlation
#' If TRUE, then \code{relations} specifies a correlation matrix (default,
#' this type of specification is usually more natural than specifying
#' a covariance matrix). Otherwise, \code{relations} specifies a
#' covariance matrix whose square root diagonal elements override
#' \code{sd_initial}.
#' @param method
#' \code{method} argument of \code{\link{mvtnorm::rmvnorm}}.
#' @template simulate_data_template
#' @template simdesign_template
#' @param ...
#' Further arguments are directly stored in the list object to be passed to 
#' \code{\link{simulate_data}}.
#'
#' @details
#' This S3 class implements a simulation design based on an underlying
#' multivariate normal distribution by creating a \code{generator} function 
#' based on \code{\link{mvtnorm::rmvnorm}}.
#' 
#' @section Data Generation:
#' Data will be generated by \code{\link{simulate_data}} using the
#' following procedure:
#' \enumerate{
#' \item The underlying data matrix *Z* is sampled from a
#' multivariate gaussian distribution (number of dimensions specified by
#' dimensions of \code{relations}).
#' \item *Z* is then transformed into the final dataset *X* by applying
#' the \code{transform} function to *Z*.
#' \item X is post-processed if specified.
#' }
#' 
#' @note
#' Note that \code{relations} specifies the correlation / covariance
#' of the underlying gaussian data *Z* and thus does not directly translate into
#' correlations between the variables of the final datamatrix *X*.
#'
#' @inherit simdesign return
#'
#' @seealso
#' \code{\link{simdesign}},
#' \code{\link{simulate_data}},
#' \code{\link{conditional_simulate_data}}
#'
#' @export
mvtnorm_simdesign <- function(relations_initial,
                              mean_initial = 0,
                              sd_initial = 1,
                              is_correlation = TRUE, 
                              method = "svd",
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
        cor_initial = cov2cor(cor_initial)
    }
    
    # define generator
    generator = function(n) mvtnorm::rmvnorm(n, 
                                             mean = mean_initial, 
                                             sigma = relations_initial, 
                                             method = method)
    
    # setup simulation design
    simdesign(
        generator = generator, 
        mean_initial = mean_initial, 
        sd_initial = sd_initial, 
        cor_initial = cor_initial,
        ...
    )
}
