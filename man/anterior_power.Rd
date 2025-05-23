% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/curvature.R
\name{anterior_power}
\alias{anterior_power}
\title{Anterior power}
\usage{
anterior_power(R)
}
\arguments{
\item{R}{Front refractive index.}
}
\value{
Returns the optical power of the anterior surface of a cornea.
}
\description{
Calculates the dioptric power of the anterior (front) surface of a cornea
based on the radius of curvature.
}
\examples{
anterior_power(7)

}
\references{
\emph{Pentacam Interpretation Guide}, pp 6-7.
}
\seealso{
Other curvature: 
\code{\link{kappa_from_radius}()},
\code{\link{posterior_power}()},
\code{\link{true_net_power}()}
}
\concept{curvature}
