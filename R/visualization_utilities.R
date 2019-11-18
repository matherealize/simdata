# TODO

#' @title Visualize fixed correlation structure as a network
#'
#' @description
#' Useful to visualize e.g. the associations of the initial multivariate
#' gaussian distribution used by \code{\link{design}} and
#' \code{\link{simulate_data}}.
#'
#' @param obj
#' Correlation matrix.
#' @param categorical_indices
#' Vector of indices of variables which should be drawn as rectangles
#' (i.e. represent categorical data).
#' @param decimals
#' Number of decimals.
#' @param cor_cutoff
#' Threshold of absolute correlation below which nodes are not considered
#' as connected. Useful to control complexity of drawn network.
#' Set to \code{NULL} to disable.
#' @param vertex_label_prefix
#' String which is added as prefix to node labels.
#' @param edge_width_function
#' Function which takes one vector input (absolute correlation values) and
#' outputs transformation of this vector (must be >= 0). Defines edge widths.
#' @param seed
#' Set random seed to ensure reproducibility of results. Can be fixed to
#' obtain same layout but vary edge widths, correlation functions etc. Can
#' also be used to obtain nicer looking graph layouts.
#'
#' @details
#' For plotting arguments see \code{igraph::plot}. Arguments via ... are
#' passed directly to \code{igraph::plot}.
#'
#' @seealso
#' \code{\link{plot_initial_cor_network}},
#' \code{\link{plot_transformed_cor_network}}
#'
#' @export
plot_cor_network <- function(obj, categorical_indices = NULL,
                             decimals = 2,
                             cor_cutoff = 0.1,
                             vertex_label_prefix = "z",
                             edge_width_function = function(x) x * 10,
                             vertex.size = 12, margin = 0, asp = 0,
                             vertex.color = "#fff7bc",
                             vertex.frame.color = "#d95f0e",
                             vertex.label.color = "black",
                             edge.color = "#fec44f",
                             edge.label.color = "black",
                             seed = NULL, ...) {
    nodes = data.frame(id = 1:nrow(obj))
    associations = cor_to_upper(obj)

    if (!is.null(cor_cutoff))
        associations = associations[abs(associations[, 3]) > cor_cutoff, ,
                                    drop = FALSE]

    edges = data.frame(from = associations[, 1], to = associations[, 2])

    net = igraph::graph_from_data_frame(d = edges, vertices = nodes,
                                        directed = FALSE)

    igraph::V(net)$label = paste0(vertex_label_prefix, nodes$id)
    igraph::V(net)$shape = "circle"
    if (!is.null(categorical_indices))
        igraph::V(net)$shape[categorical_indices] = "rectangle"

    igraph::E(net)$width = edge_width_function(abs(obj[cbind(edges$from, edges$to)]))
    igraph::E(net)$label = round(obj[cbind(edges$from, edges$to)], decimals)

    if (!is.null(seed))
        set.seed(seed)

    layout = igraph::layout_with_fr(net, start.temp = igraph::vcount(net))

    op = par(mar = c(0, 0, 0, 0))
    plot(net,
         vertex.size = vertex.size,
         margin = margin,
         layout = layout, rescale = TRUE,
         asp = asp,
         vertex.color = vertex.color,
         vertex.frame.color = vertex.frame.color,
         vertex.label.color = vertex.label.color,
         edge.color = edge.color,
         edge.label.color = edge.label.color, ...)
    par(op)
}

#' @title Visualize correlation of initial
#' multivariate gaussian distribution of simulation design
#'
#' @param obj
#' \code{design} S3 class object.
#'
#' @details
#' For explanation of further parameters, see \code{\link{plot_cor_network}}.
#'
#' @seealso
#' \code{\link{plot_cor_network}}
#'
#' @export
plot_initial_cor_network <- function(obj, ...) {
    plot_cor_network(obj$cor_initial, ...)
}

#' @title Visualize estimated correlation of final simulated datamatrix
#'
#' @description
#' Based on approximation via simulation of the final datamatrix specified by
#' given simulation design.
#'
#' @param obj
#' \code{design} S3 class object.
#' @param n_obs
#' Number of simulated observations.
#' @param cor_type
#' Either a character ("p" for Pearson correlation, "s" for Spearman
#' correlation) or a user specified function, taking one input (a data.frame)
#' and returning a matrix.
#' @param show_categorical
#' If TRUE, plots categorical variables specially. Determined by
#' \code{transform_type} slot of \code{design} object.
#'
#' @details
#' For explanation of further parameters, see \code{\link{plot_cor_network}}.
#'
#' @seealso
#' \code{\link{plot_cor_network}}
#'
#' @export
plot_final_cor_network <- function(obj, n_obs = 100000,
                                   cor_type = "p",
                                   seed = NULL,
                                   show_categorical = TRUE, ...) {
    if (!is.null(seed))
        set.seed(seed)

    sim_data = simulate_data(obj, n_obs = n_obs)

    if (show_categorical) {
        categorical_indices = which(obj$transform_type == "logical" |
                                        obj$transform_type == "factor")
    } else categorical_indices = NULL

    f_cor = function(x) cor(x, method = "p")
    if (is.character(cor_type) & cor_type == "s") {
        f_cor = function(x) cor(x, method = "s")
    }
    if (class(cor_type) == "function")
        f_cor = cor_type

    relations = f_cor(sim_data)
    plot_cor_network(relations, seed = seed,
                     vertex_label_prefix = "x",
                     categorical_indices = categorical_indices, ...)
}
