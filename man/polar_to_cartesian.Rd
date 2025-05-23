% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coordinate_systems.R
\name{polar_to_cartesian}
\alias{polar_to_cartesian}
\title{Convert the polar coordinates of a point to Cartesian}
\usage{
polar_to_cartesian(radius, angle_deg)
}
\arguments{
\item{radius}{The radius of the polar coordinates.}

\item{angle_deg}{The angle (in degrees) of the polar coordinates.}
}
\value{
A data frame with the Cartesian coordinates, \code{x} and \code{y}.
}
\description{
Converts a point from polar (\emph{radius}, \emph{theta}) to
Cartesian (\emph{x}, \emph{y}) coordinates.
}
\examples{
polar_to_cartesian(1, 45)

}
\seealso{
Other coordinate systems: 
\code{\link{cartesian_degrees}()},
\code{\link{cartesian_to_polar}()},
\code{\link{deg_to_rad}()},
\code{\link{euclidean_dist}()},
\code{\link{hausdorff_distance}()},
\code{\link{interpolate_measurements}()},
\code{\link{quadrant_eye}()},
\code{\link{rad_to_deg}()},
\code{\link{rotate_shape}()},
\code{\link{signed_smallest_angle}()},
\code{\link{smallest_angle}()},
\code{\link{within_2pi}()},
\code{\link{within_360}()}
}
\concept{coordinate systems}
