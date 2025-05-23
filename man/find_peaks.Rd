% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/k_next.R
\name{find_peaks}
\alias{find_peaks}
\title{Find the peak points of a numeric matrix.}
\usage{
find_peaks(input_matrix)
}
\arguments{
\item{input_matrix}{A numeric matrix.}
}
\value{
A numeric matrix with the same dimensions as \code{input_matrix}. The values
in the matrix identify the peaks and the points that descend from them. The
highest peak is identified by the point(s) equal to 1, and the points that
descend from that peak are identified by -1. If there is a next highest peak,
it is identified by the point(s) equal to 2, and the points that descend from
it are identified by -2. And so on until all of the points in
\code{input_matrix} have been assigned to a peak.

Any points that was \code{NA} in the \code{input_matrix} will be assigned 0.
}
\description{
Finds the peak points (i.e., max values) and the points that descend from
them in a numeric matrix.
}
\examples{
test_m <- matrix(c(12, 12, 11, 12, 12,
                   14, 14, 12, 13, 13,
                   14, 15, 13, 14, 13,
                   13, 14, 13, 14, 12,
                   12, 12, 12, 12, 12),
                 nrow = 5, byrow = T)
find_peaks(test_m)

}
\seealso{
Other K-next: 
\code{\link{antipodal_power}()},
\code{\link{k_next}()}
}
\concept{K-next}
