% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simdesign.R
\name{simdesign}
\alias{simdesign}
\title{Design specification for simulating datasets}
\usage{
simdesign(
  generator,
  transform_initial = base::identity,
  n_var_final = -1,
  types_final = NULL,
  names_final = NULL,
  prefix_final = "v",
  process_final = list(),
  name = "Simulation design",
  check_and_infer = TRUE,
  ...
)
}
\arguments{
\item{generator}{Function which generates data from the underlying base distribution. It is
assumed it takes the number of simulated observations \code{n_obs} as first
argument, as all random generation functions in the \pkg{stats} and
\pkg{extraDistr} do. Furthermore, it is expected to return a two-dimensional
array as output (matrix or data.frame). Alternatively an R object derived
from the \code{simdata::simdesign} class. See details.}

\item{transform_initial}{Function which specifies the transformation of the underlying
dataset \code{Z} to final dataset \code{X}. See details.}

\item{n_var_final}{Integer, number of columns in final datamatrix \code{X}. Can be inferred when
\code{check_and_infer} is TRUE.}

\item{types_final}{Optional vector of length equal to \code{n_var_final} (set by the user or
inferred) and hence number of columns of final dataset \code{X}.
Allowed entries are "logical", "factor" and "numeric".
Stores the type of the columns of \code{X}.
If not specified by, inferred if \code{check_and_infer} is set to TRUE.}

\item{names_final}{NULL or character vector with variable names for final dataset \code{X}.
Length needs to equal the number of columns of \code{X}.
Overrides other naming options. See details.}

\item{prefix_final}{NULL or prefix attached to variables in final dataset \code{X}. Overriden
by \code{names_final} argument. Set to NULL if no prefixes should
be added. See details.}

\item{process_final}{List of lists specifying post-processing functions applied to final
datamatrix \code{X} before returning it. See \code{\link{do_processing}}.}

\item{name}{Character, optional name of the simulation design.}

\item{check_and_infer}{If TRUE, then the simulation design is tested by simulating 5 observations
using \code{\link{simulate_data}}. If everything works without error,
the variables \code{n_var_final} and \code{types_final} will be inferred
from the results if not already set correctly by the user.}

\item{...}{Further arguments are directly stored in the list object to be passed to
\code{\link{simulate_data}}.}
}
\value{
List object with class attribute "simdesign" (S3 class) containing
the following entries (if no further information given, entries are
directly saved from user input):

\describe{
\item{\code{generator}}{}
\item{\code{name}}{}
\item{\code{transform_initial}}{}
\item{\code{n_var_final}}{}
\item{\code{types_final}}{}
\item{\code{names_final}}{}
\item{\code{process_final}}{}
\item{\verb{entries for further information as passed by the user}}{}
}
}
\description{
Stores information necessary to simulate and visualize datasets based
on underlying distribution \code{Z}.
}
\details{
The \code{simdesign} class should be used in the following workflow:

\enumerate{
\item Specify a design template which will be used in subsequent data
generating / visualization steps.
\item Sample / visualize datamatrix following template (possibly
multiple times) using \code{\link{simulate_data}}.
\item Use sampled datamatrix for simulation study.
}

For more details on generators and transformations, please see the
documentation of \code{\link{simulate_data}}.

For details on post-processing, please see the documentation of
\code{\link{do_processing}}.
}
\section{Naming of variables}{

If \code{check_and_infer} is set to TRUE, the following procedure determines
the names of the variables:

\enumerate{
\item use \code{names_final} if specified and of correct length
\item otherwise, use the names of \code{transform_initial} if present and of
correct length
\item otherwise, use \code{prefix_final} to prefix the variable number if
not NULL
\item otherwise, use names from dataset as generated by the \code{generator}
function
}
}

\section{Simulation Templates}{

This class is intended to be used as a template for simulation designs
which are based on specific underlying distributions. All such a template
needs to define is the \code{generator} function and its construction and
pass it to this function along with the other arguments. See
\code{\link{simdesign_mvtnorm}} for an example.
}

\examples{
generator <- function(n) mvtnorm::rmvnorm(n, mean = 0)
sim_design <- simdesign(generator)
simulate_data(sim_design, 10, seed = 19)

}
\seealso{
\code{\link{simdesign_mvtnorm}},
\code{\link{simulate_data}},
\code{\link{simulate_data_conditional}}
}
