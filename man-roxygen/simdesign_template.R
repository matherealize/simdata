#' @param n_var_final
#' Integer, number of columns in final datamatrix *X*. Can be inferred when 
#' \code{check_and_infer} is TRUE.
#' @param types_final
#' Optional vector of length equal to \code{n_var_final} (set by the user or 
#' inferred) and hence number of columns of final dataset *X*. 
#' Allowed entries are "logical", "factor" and "numeric". 
#' Stores the type of the columns of *X*. 
#' If not specified by, inferred if \code{check_and_infer} is set to TRUE.
#' @param name
#' Character, optional name of the simulation design.
#' @param check_and_infer
#' If TRUE, then the simulation design is tested by simulating 5 observations
#' using \code{\link{simulate_data}}. If everything works without error, 
#' the variables \code{n_var_final} and \code{types_final} will be inferred
#' from the results if not already set correctly by the user.
