% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simdesign.R
\name{discunif_simdesign}
\alias{discunif_simdesign}
\title{Uniform disc sampling design specification}
\usage{
discunif_simdesign(
  r_min = 0,
  r_max = 1,
  angle_min = 0,
  angle_max = 2 * pi,
  name = "Uniform circle simulation design",
  ...
)
}
\arguments{
\item{r_min}{Minimum radius of points.}

\item{r_max}{Maximum radius of points.}

\item{angle_min}{Minimum angle of points (between 0 and 2pi).}

\item{angle_max}{Maximum angle of points (between 0 and 2pi).}
}
\description{
Provides 2-dimensional points, spread uniformly over disc, or partial
disc segment (i.e. a circle, or ring, or ring segment). Useful for e.g.
building up clustering exercises.
}
\details{
The distribution of points on a disk depends on the radius - the farther out,
the more area the points need to cover. Thus, simply sampling two uniform
values for radius and angle will not work. See references.
}
\examples{
\dontrun{
disc_sampler = discunif_simdesign()
plot(simulate_data(disc_sampler, 1000))

ring_segment_sampler = discunif_simdesign(r_min = 0.5, angle_min = 0.5*pi)
plot(simulate_data(ring_segment_sampler, 1000))

circle_sampler = discunif_simdesign(r_min = 1)
plot(simulate_data(circle_sampler, 1000))
} 

}
\references{
\url{https://mathworld.wolfram.com/DiskPointPicking.html}
}
