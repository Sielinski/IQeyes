% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/contours.R
\name{get_contour_levels}
\alias{get_contour_levels}
\title{Get contour levels}
\usage{
get_contour_levels(exam_curvature)
}
\arguments{
\item{exam_curvature}{A data frame of numeric values containing \code{x} and \code{y}, the
Cartesian coordinates of a point, and either \code{z} (the anterior radius of
curvature) or \code{power} (the dioptric power) at each (\code{x}, \code{y}).}
}
\value{
A table. The table's labels identify the contour levels, and the values in
the table identify the number of segments that comprise the contour.
}
\description{
Returns the full range of contours and the number of segments needed to plot
\code{exam_curvature} on a curvature map.
}
\details{
This function returns valid contours for an exam based on the breaks defined
by the color scale of \link{plot_scale}.

Because \link{plot_scale} is designed for the dioptric power of anterior
surfaces, this function works for anterior surfaces only.
}
\examples{
sample_curvature |>
  interpolate_measurements() |>
  get_contour_levels()

}
\seealso{
Other contours: 
\code{\link{get_contour}()},
\code{\link{scale_rotate}()}
}
\concept{contours}
