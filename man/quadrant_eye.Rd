% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coordinate_systems.R
\name{quadrant_eye}
\alias{quadrant_eye}
\title{Quadrant of an eye}
\usage{
quadrant_eye(x, y, exam_eye)
}
\arguments{
\item{x}{The x-axis coordinate.}

\item{y}{The y-axis coordinate.}

\item{exam_eye}{A character vector containing either "left" or "right".}
}
\value{
A two-character string. See details.
}
\description{
Returns the quadrant of a point based on its \emph{x}- and \emph{y}-axis
position and the \code{exam_eye}.
}
\details{
The quadrant system for an eye is split on the \emph{y}-axis by the superior
("S") and the inferior ("I"), and on the \emph{x}-axis by the nasal ("N") and
the temporal ("T"). The superior is always \code{y > 0}, and the inferior is
always \code{y ≤ 0}. The nasal and temporal depend on the eye: For the right
eye, the nasal is \code{x > 0}; for the left eye, the nasal is \code{x ≤ 0}.
}
\examples{
quadrant_eye(1, 1, 'left')
quadrant_eye(1, 1, 'right')

}
\seealso{
Other coordinate systems: 
\code{\link{cartesian_degrees}()},
\code{\link{cartesian_to_polar}()},
\code{\link{deg_to_rad}()},
\code{\link{euclidean_dist}()},
\code{\link{hausdorff_distance}()},
\code{\link{interpolate_measurements}()},
\code{\link{polar_to_cartesian}()},
\code{\link{rad_to_deg}()},
\code{\link{rotate_shape}()},
\code{\link{signed_smallest_angle}()},
\code{\link{smallest_angle}()},
\code{\link{within_2pi}()},
\code{\link{within_360}()}
}
\concept{coordinate systems}
