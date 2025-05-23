% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/curvature.R
\name{true_net_power}
\alias{true_net_power}
\title{True net power}
\usage{
true_net_power(R_front, R_back)
}
\arguments{
\item{R_front}{Front refractive index.}

\item{R_rear}{Front refractive index.}
}
\value{
Returns the optical power of the cornea.
}
\description{
Calculates the true net power of a cornea based on the front (anterior)
and rear (posterior) refractive indices.
}
\examples{
true_net_power(7, 5)

}
\references{
\emph{Pentacam Interpretation Guide}, p 9.
}
\seealso{
Other curvature: 
\code{\link{anterior_power}()},
\code{\link{kappa_from_radius}()},
\code{\link{posterior_power}()}
}
\concept{curvature}
