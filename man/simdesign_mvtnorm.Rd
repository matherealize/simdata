% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simdesign.R
\name{simdesign_mvtnorm}
\alias{simdesign_mvtnorm}
\title{Multivariate normal design specification}
\usage{
simdesign_mvtnorm(
  relations_initial,
  mean_initial = 0,
  sd_initial = 1,
  is_correlation = TRUE,
  method = "svd",
  name = "Multivariate-normal based simulation design",
  ...
)
}
\arguments{
\item{relations_initial}{Correlation / Covariance matrix of the initial multivariate
Normal distribution \code{Z}.}

\item{mean_initial}{Vector of mean values of the initial multivariate Normal
distribution \code{Z}. Dimension needs to correspond to dimension
of \code{relations}.}

\item{sd_initial}{Vector of standard deviations of the initial multivariate
Normal distribution Z. Dimension needs to correspond to dimension
of \code{relations}. Overriden by suqare root of diagonal elements of
\code{relations} if \code{is_correlation} is FALSE.}

\item{is_correlation}{If TRUE, then \code{relations} specifies a correlation matrix (default,
this type of specification is usually more natural than specifying
a covariance matrix). Otherwise, \code{relations} specifies a
covariance matrix whose square root diagonal elements override
\code{sd_initial}.}

\item{method}{\code{method} argument of \code{\link[mvtnorm:Mvnorm]{mvtnorm::rmvnorm}}.}

\item{name}{Character, optional name of the simulation design.}

\item{...}{Further arguments are passed to the \code{\link{simdesign}} constructor.}
}
\value{
List object with class attribute "simdesign_mvtnorm" (S3 class), inheriting
from "simdesign". It contains the same entries as a \code{\link{simdesign}}
object but in addition the following entries:

\describe{
\item{\code{mean_initial}}{}
\item{\code{sd_initial}}{}
\item{\code{cor_initial}}{Initial correlation matrix of multivariate normal
distribution}
}
}
\description{
Stores information necessary to simulate and visualize datasets based
on underlying distribution multivariate normal distribution \code{Z}.
}
\details{
This S3 class implements a simulation design based on an underlying
multivariate normal distribution by creating a \code{generator} function
based on \code{\link[mvtnorm:Mvnorm]{mvtnorm::rmvnorm}}.
}
\note{
Note that \code{relations} specifies the correlation / covariance
of the underlying Normal data \code{Z} and thus does not directly translate into
correlations between the variables of the final datamatrix \code{X}.
}
\section{Data Generation}{

Data will be generated by \code{\link{simulate_data}} using the
following procedure:
\enumerate{
\item The underlying data matrix \code{Z} is sampled from a
multivariate Normal distribution (number of dimensions specified by
dimensions of \code{relations}).
\item \code{Z} is then transformed into the final dataset \code{X} by applying
the \code{transform_initial} function to \code{Z}.
\item \code{X} is post-processed if specified.
}
}

\seealso{
\code{\link{simdesign}},
\code{\link{simulate_data}},
\code{\link{simulate_data_conditional}},
\code{\link{plot_cor_network.simdesign_mvtnorm}}
}
