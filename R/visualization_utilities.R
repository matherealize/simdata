#' @title Visualize fixed correlation structure as a network
#'
#' @description
#' Useful to visualize e.g. the associations of the initial multivariate
#' gaussian distribution used by \code{\link{simdesign_mvtnorm}}.
#'
#' @param obj
#' Correlation matrix or S3 class object which has a class method available (see below).
#' @param categorical_indices
#' Vector of indices of variables which should be drawn as rectangles
#' (i.e. represent categorical data).
#' @param decimals
#' Number of decimals, used for default labeling of the network edges.
#' @param cor_cutoff
#' Threshold of absolute correlation below which nodes are not considered
#' as connected. Useful to control complexity of drawn network.
#' Set to NULL to disable.
#' @param vertex.color
#' Argument passed to \code{\link[igraph:plot.igraph]{igraph::plot}}. Usually a
#' character vector with a hex color specification for vertex color. 
#' Alternatively a function that takes as input a data.frame with a column "id" 
#' that gives the column number of the simulated data, and outputs a valid color
#' specification for the corresponding vertices (i.e. a single character hex
#' color or a vector of such hex colors of appropriate length).
#' @param vertex_labels
#' Character vector of length `nrow(obj)` of labels for vertices. If not NULL,
#' overrides the `vertex_label_prefix` argument. If set to NA omits all
#' or some vertex labels.
#' @param vertex_label_prefix
#' String which is added as prefix to node labels.
#' @param edge.color
#' Argument passed to \code{\link[igraph:plot.igraph]{igraph::plot}}. This 
#' package implements some special functionality: if `edge.color = "ramp"` then
#' a colorramp from red (-1) via white (0) to blue (1) is mapped to the
#' correlations and the edges colored accordingly.
#' If `edge.color = "clipped-ramp"` then the ramp is restricted to the 
#' correlation values observed, which may be useful if they are low to increase
#' visibility. 
#' If `edge.color = "red-blue"` then all edges with positive correlation values
#' are colored uniformly red, and all edges with negative correlations are 
#' colored uniformly blue. Alternatively, may be a function that takes as 
#' input the edge correlation values and outputs valid color specifications 
#' (i.e. a single hex color or a vector of hex colors of appropriate length).
#' @param edge_width_function
#' Function which takes one vector input (absolute correlation values) and
#' outputs transformation of this vector (must be >= 0). Defines edge widths.
#' @param edge_label_function
#' Function which takes on vector input (absolute correlation values) and
#' outputs labels for these values as character vector. Defines edges labels.
#' If set to NULL, then no edge labels will be displayed.
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
#' the user using e.g. the interactive \code{\link[igraph:tkplot]{igraph::tkplot}}
#' function.
#' @param mar
#' `mar` argument to the \code{\link[graphics:par]{par}} function to set
#' margins of the plot (often required when the axes should be drawn).
#' A numerical vector of the form c(bottom, left, top, right) which gives the
#' number of lines of margin to be specified on the four sides of the plot.
#' The default is c(5, 4, 4, 2) + 0.1. Note that this is not the same argument
#' as the `margin` argument for the `igraph::plot.igraph` function.
#' @param vertex.size,margin,asp,vertex.frame.color,vertex.label.color,edge.label.color,edge.label.cex
#' Arguments to \code{\link[igraph:plot.igraph]{igraph::plot}}, with sensible
#' defaults for this package's usage.
#' @param ...
#' Passed to \code{\link[igraph:plot.igraph]{igraph::plot}}, with a complete list
#' of arguments and details given in \code{\link[igraph:plot.common]{igraph.plotting}}.
#'
#' @details
#' For an explanation of all parameters not listed here, please refer to
#' \code{\link[igraph:plot.igraph]{igraph::plot}}.
#'
#' @seealso
#' \code{\link{plot_cor_network.simdesign_mvtnorm}},
#' \code{\link{plot_estimated_cor_network}}
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
                                     edge_label_function = function(x) round(x, decimals),
                                     use_edge_weights = FALSE,
                                     edge_weight_function = base::identity,
                                     seed = NULL,
                                     return_network = FALSE,
                                     mar = c(0, 0, 0, 0),
                                     vertex.size = 12,
                                     margin = 0, 
                                     asp = 0,
                                     vertex.color = "#ececec",
                                     vertex.frame.color = "#979797",
                                     vertex.label.color = "black",
                                     edge.color = "ramp",
                                     edge.label.color = "black",
                                     edge.label.cex = 0.8,
                                     ...) {
    
    nodes <- data.frame(id = 1:ncol(obj))
    associations <- cor_to_upper(obj)
    
    if (!is.null(cor_cutoff))
        associations <- associations[abs(associations[, 3]) > cor_cutoff, ,
                                     drop = FALSE]
    
    edges <- data.frame(from = associations[, 1], to = associations[, 2])
    connection_indices <- cbind(edges$from, edges$to)
    edge_correlations <- obj[connection_indices]
    
    net <- igraph::graph_from_data_frame(d = edges,
                                         vertices = nodes,
                                         directed = FALSE)
    
    if (is.null(vertex_labels)) {
        igraph::V(net)$label <- paste0(vertex_label_prefix, nodes$id)
    } else igraph::V(net)$label <- vertex_labels
    
    igraph::V(net)$shape <- "circle"
    if (!is.null(categorical_indices))
        igraph::V(net)$shape[categorical_indices] <- "rectangle"
    
    igraph::E(net)$width <- edge_width_function(abs(edge_correlations))
    if (!is.null(edge_label_function)) {
        igraph::E(net)$label <- edge_label_function(edge_correlations)
    }
    
    weights <- NULL
    if (use_edge_weights) {
        weights <- edge_weight_function(abs(edge_correlations))
        igraph::E(net)$weight <- weights
    }
    
    if (is.function(edge.color)) {
        edge.color <- edge.color(edge_correlations)
    } else if (tolower(edge.color) == "ramp") {
        # palette with 21 segments, white in the middle
        segments <- seq(-1, 1, 2 / 21)
        pal <- grDevices::colorRampPalette(colors = c("#003ed6", "white", "#d60000"))
        pal_ind <- cut(edge_correlations, breaks = segments, 
                       include.lowest = TRUE, 
                       ordered_result = TRUE)
        edge.color <- pal(length(segments) - 1)[as.integer(pal_ind)]
    } else if (tolower(edge.color) == "clipped-ramp") {
        # clip separately for positive and negative correlations
        # else problems when data only has one of those
        lower <- min(min(edge_correlations), -0.001)
        upper <- max(max(edge_correlations), 0.001)
        # scale segments
        segments <- seq(-1, 1, 2 / 21)
        segments[1:11] <- segments[1:11] * abs(lower)
        segments[12:22] <- segments[12:22] * upper
        pal <- grDevices::colorRampPalette(colors = c("#003ed6", "white", "#d60000"))
        pal_ind <- cut(edge_correlations, breaks = segments, 
                       include.lowest = TRUE, 
                       ordered_result = TRUE)
        edge.color <- pal(length(segments) - 1)[as.integer(pal_ind)]
    } else if (tolower(edge.color) == "red-blue") {
        edge.color <- ifelse(edge_correlations > 0, "#d60000", "#003ed6")
    }
    
    if (is.function(vertex.color)) {
        vertex.color <- vertex.color(nodes)
    } 
    
    if (return_network)
        return(net)
    
    if (!is.null(seed))
        set.seed(seed)
    
    layout <- igraph::layout_with_fr(net, start.temp = igraph::vcount(net))
    
    op <- graphics::par(mar = mar)
    invisible(on.exit(graphics::par(op)))
    igraph::plot.igraph(
        net,
        vertex.size = vertex.size,
        margin = margin,
        layout = layout, 
        rescale = TRUE,
        asp = asp,
        vertex.color = vertex.color,
        vertex.frame.color = vertex.frame.color,
        vertex.label.color = vertex.label.color,
        edge.color = edge.color,
        edge.label.color = edge.label.color, 
        edge.label.cex = edge.label.cex,
        ...
    )
}

#' @describeIn plot_cor_network Function to be used with \code{\link{simdesign_mvtnorm}}
#' S3 class object to visualize initial correlation network of the underlying
#' multivariate normal distribution.
#'
#' @export
#' @method plot_cor_network simdesign_mvtnorm
plot_cor_network.simdesign_mvtnorm <- function(obj, ...) {
    plot_cor_network(obj$cor_initial, ...)
}

#' @title Visualize estimated correlation matrix as a network
#'
#' @description
#' Based on approximation via simulation specified by given simulation design.
#' Convenience wrapper for combining \code{\link{estimate_final_correlation}} and
#' \code{\link{plot_cor_network}}.
#'
#' @inheritParams estimate_final_correlation
#' @param show_categorical
#' If TRUE, marks categorical variables differently from numeric ones.
#' Determined by the `types_final` slot of the `obj` argument.
#' @param ...
#' Passed to \code{\link{plot_cor_network}}.
#'
#' @details
#' This function is useful to estimate the correlation network of a simulation
#' setup after the initial underlying distribution `Z` has been transformed to
#' the final dataset `X`.
#'
#' @seealso
#' \code{\link{plot_cor_network}},
#' \code{\link{estimate_final_correlation}}
#'
#' @export
plot_estimated_cor_network <- function(obj, n_obs = 100000,
                                       cor_type = "pearson",
                                       seed = NULL,
                                       show_categorical = TRUE,
                                       ...) {

    if (show_categorical) {
        categorical_indices <- which(obj$types_final == "logical" |
            obj$types_final == "factor")
    } else categorical_indices <- NULL

    relations <- estimate_final_correlation(obj,
        n_obs = n_obs,
        cor_type = cor_type,
        seed = seed)

    plot_cor_network(relations,
        seed = seed,
        categorical_indices = categorical_indices,
        vertex_labels = obj$names_final, ...)
}
