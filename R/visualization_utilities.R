# TODO

#' @title Visualize fixed correlation structure as a network
#'
#' @description
#' Useful to visualize e.g. the associations of the initial multivariate
#' gaussian distribution used by `\link{mvtnorm_simdesign}`.
#'
#' @param obj
#' Correlation matrix or S3 class object which has a class method available (see below).
#' @param categorical_indices
#' Vector of indices of variables which should be drawn as rectangles
#' (i.e. represent categorical data).
#' @param decimals
#' Number of decimals.
#' @param cor_cutoff
#' Threshold of absolute correlation below which nodes are not considered
#' as connected. Useful to control complexity of drawn network.
#' Set to NULL to disable.
#' @param vertex_labels
#' Character vector of length `nrow(obj)` of labels for vertices. If not NULL, 
#' overrides the `vertex_label_prefix` argument.
#' @param vertex_label_prefix
#' String which is added as prefix to node labels.
#' @param edge_width_function
#' Function which takes one vector input (absolute correlation values) and
#' outputs transformation of this vector (must be >= 0). Defines edge widths.
#' @param use_edge_weights
#' Logical, if TRUE then the layout will be influenced by the absolute 
#' correlations (i.e. edge weights) such that highly correlated variables will
#' be put closer together.
#' If FALSE, then the layout is independent of the correlation structure.
#' @param edge_weight_function
#' Function which takes one vector input (absolute correlation values) and
#' outputs transformation of this vector (must be >= 0). Defines edge weights.
#' Only relevant if `use_edge_weights` is TRUE.
#' @param seed
#' Set random seed to ensure reproducibility of results. Can be fixed to
#' obtain same layout but vary edge widths, correlation functions etc. Can
#' also be used to obtain nicer looking graph layouts.
#' @param return_network
#' If TRUE, the `igraph` network object is returned and can be plotted by 
#' the user using e.g. the interactive `\link[igraph:tkplot]{igraph::tkplot}` 
#' function.
#' @param mar
#' `mar` argument to the `\link[graphics:par]{par}` function to set 
#' margins of the plot (often required when the axes should be drawn). 
#' A numerical vector of the form c(bottom, left, top, right) which gives the 
#' number of lines of margin to be specified on the four sides of the plot. 
#' The default is c(5, 4, 4, 2) + 0.1. Note that this is not the same argument
#' as the `margin` argument for the `igraph::plot.igraph` function.
#' @param ...
#' Passed to `\link[igraph:plot.igraph]{igraph::plot}`.
#' 
#' @details 
#' For an explanation of all parameters not listed here, please refer to
#' `\link[igraph:plot.igraph]{igraph::plot}`.
#'
#' @seealso
#' `\link{plot_cor_network.mvtnorm_simdesign}`,
#' `\link{plot_estimated_cor_network}`
#'
#' @export
plot_cor_network <- function(obj, ...) {
    UseMethod("plot_cor_network", obj)
}

#' @describeIn plot_cor_network Function to be used for correlation matrix.
#'
#' @export
#' @method plot_cor_network default
plot_cor_network.default <- function(obj, categorical_indices = NULL,
                                     decimals = 2,
                                     cor_cutoff = 0.1,
                                     vertex_labels = NULL,
                                     vertex_label_prefix = "z",
                                     edge_width_function = function(x) x * 10,
                                     use_edge_weights = FALSE,
                                     edge_weight_function = base::identity,
                                     seed = NULL,
                                     return_network = FALSE, 
                                     mar = c(0, 0, 0, 0),
                                     vertex.size = 12, margin = 0, asp = 0,
                                     vertex.color = "#fff7bc",
                                     vertex.frame.color = "#d95f0e",
                                     vertex.label.color = "black",
                                     edge.color = "#fec44f",
                                     edge.label.color = "black",
                                     ...) {
    nodes = data.frame(id = 1:ncol(obj))
    associations = cor_to_upper(obj)

    if (!is.null(cor_cutoff))
        associations = associations[abs(associations[, 3]) > cor_cutoff, ,
                                    drop = FALSE]

    edges = data.frame(from = associations[, 1], to = associations[, 2])

    net = igraph::graph_from_data_frame(d = edges, 
                                        vertices = nodes,
                                        directed = FALSE)

    if (is.null(vertex_labels)) {
        igraph::V(net)$label = paste0(vertex_label_prefix, nodes$id)    
    } else igraph::V(net)$label = vertex_labels
    
    igraph::V(net)$shape = "circle"
    if (!is.null(categorical_indices))
        igraph::V(net)$shape[categorical_indices] = "rectangle"

    connection_indices = cbind(edges$from, edges$to)
    igraph::E(net)$width = edge_width_function(abs(obj[connection_indices]))
    igraph::E(net)$label = round(obj[connection_indices], decimals)
    
    weights = NULL
    if (use_edge_weights) {
        weights = edge_weight_function(abs(obj[connection_indices]))
        igraph::E(net)$weight = weights   
    }
    
    if (return_network)
        return(net)

    if (!is.null(seed))
        set.seed(seed)

    layout = igraph::layout_with_fr(net, start.temp = igraph::vcount(net))

    op = par(mar = mar)
    invisible(on.exit(par(op)))
    igraph::plot.igraph(
        net,
        vertex.size = vertex.size,
        margin = margin,
        layout = layout, rescale = TRUE,
        asp = asp,
        vertex.color = vertex.color,
        vertex.frame.color = vertex.frame.color,
        vertex.label.color = vertex.label.color,
        edge.color = edge.color,
        edge.label.color = edge.label.color, ...)
}

#' @describeIn plot_cor_network Function to be used with `\link{mvtnorm_simdesign}` 
#' S3 class object to visualize initial correlation network of the underlying
#' multivariate normal distribution.
#'
#' @export
#' @method plot_cor_network mvtnorm_simdesign
plot_cor_network.mvtnorm_simdesign <- function(obj, ...) {
    plot_cor_network(obj$cor_initial, ...)
}

#' @title Visualize estimated correlation matrix as a network
#'
#' @description
#' Based on approximation via simulation specified by given simulation design.
#' Convenience wrapper for combining `\link{estimate_final_correlation}` and
#' `\link{plot_cor_network}`.
#'
#' @inheritParams estimate_final_correlation
#' @param show_categorical
#' If TRUE, marks categorical variables differently from numeric ones.
#' Determined by the `types_final` slot of the `obj` argument.
#' @param ...
#' Passed to `\link{plot_cor_network}`.
#'
#' @details 
#' This function is useful to estimate the correlation network of a simulation
#' setup after the initial underlying distribution `Z` has been transformed to 
#' the final dataset `X`.
#'
#' @seealso
#' `\link{plot_cor_network}`,
#' `\link{estimate_final_correlation}`
#'
#' @export
plot_estimated_cor_network <- function(obj, n_obs = 100000,
                                       cor_type = "pearson",
                                       seed = NULL,
                                       show_categorical = TRUE, 
                                       ...) {
    
    if (show_categorical) {
        categorical_indices = which(obj$types_final == "logical" |
                                        obj$types_final == "factor")
    } else categorical_indices = NULL

    relations = estimate_final_correlation(obj, 
                                           n_obs = n_obs, 
                                           cor_type = cor_type, 
                                           seed = seed)
    
    plot_cor_network(relations, 
                     seed = seed,
                     categorical_indices = categorical_indices, 
                     vertex_labels = obj$names_final, ...)
}
