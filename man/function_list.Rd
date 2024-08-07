% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{function_list}
\alias{function_list}
\title{Apply list of functions to input}
\usage{
function_list(..., stringsAsFactors = FALSE, check.names = TRUE)
}
\arguments{
\item{...}{Named or unnamed arguments, each of which is a function taking exactly
one input. See details.}

\item{stringsAsFactors, check.names}{Arguments of \code{\link[base:data.frame]{data.frame}}.}
}
\value{
Function with a single input which outputs a data.frame. Has special
'flist' entry in its environment which stores individual functions as list.
}
\description{
Apply list of functions to input
}
\details{
This is a convenience function which takes a number of functions and returns
another function which applies all of the user specified functions to a new
input, and collects the results as list or data.frame.
This is useful to e.g. transform columns of a data.frame or check
the validity of a matrix during simulations. See the example here and
in \code{\link{simulate_data_conditional}}.

The assumptions for the individual functions are:

\itemize{
\item Each function is expected to take a single input.
\item Each function is expected to output a result consistent with the
other functions (i.e. same output length) to ensure that the results can be
summarized as a data.frame.
}
}
\note{
This function works fine without naming the input arguments, but the
resulting data.frames have empty column names if that is the case. Thus,
it is recommended to only pass named function arguments.
}
\examples{
f <- function_list(
    v1 = function(x) x[, 1] * 2,
    v2 = function(x) x[, 2] + 10)

f(diag(2))

# function_list can be used to add new columns
# naming of columns should be handled separately in such cases

f <- function_list(
    function(x) x, # return x as it is
    X1_X2 = function(x) x[, 2] + 10) # add new column

f(diag(2))

}
\seealso{
\code{\link[base:data.frame]{data.frame}},
\code{\link{get_from_function_list}},
\code{\link{get_names_from_function_list}}
}
