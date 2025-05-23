% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coordinate_systems.R
\name{rotate_shape}
\alias{rotate_shape}
\title{Rotate a shape}
\usage{
rotate_shape(points, theta_deg)
}
\arguments{
\item{points}{An \emph{n} x 2 matrix of Cartesian coordinates: \emph{n} rows, one for each
point, and two columns, the \emph{x} and \emph{y} coordinates of a point,
respectively.}

\item{theta_deg}{Angle (in degrees) to rotate the shape.}
}
\value{
An \emph{n}x2 matrix containing the points of original matrix rotated by
\emph{theta} degrees.
}
\description{
Rotate a shape (i.e., a matrix of points) by a specified angle (in degrees).
}
\examples{
matrix(rep(seq(1:5), 2), ncol = 2) |> rotate_shape(180)

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
\code{\link{quadrant_eye}()},
\code{\link{rad_to_deg}()},
\code{\link{signed_smallest_angle}()},
\code{\link{smallest_angle}()},
\code{\link{within_2pi}()},
\code{\link{within_360}()}
}
\concept{coordinate systems}
