% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{get_from_function_list}
\alias{get_from_function_list}
\title{Extract individual functions from \code{function_list}}
\usage{
get_from_function_list(flist)
}
\arguments{
\item{flist}{\code{function_list} or function object.}
}
\value{
List with named or unnamed entries corresponding to individual function
objects that were passed to the \code{function_list} object. If \code{flist} is a
simple function, returns NULL.
}
\description{
Extract individual function objects from environment of a \code{function_list}
object.
}
\seealso{
\code{\link{function_list}}
}
