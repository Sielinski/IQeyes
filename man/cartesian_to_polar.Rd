% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coordinate_systems.R
\name{cartesian_to_polar}
\alias{cartesian_to_polar}
\title{Convert the Cartesian coordinates of a point to polar coordinates}
\usage{
cartesian_to_polar(x, y)
}
\arguments{
\item{x}{The \emph{x}-axis component.}

\item{y}{The \emph{y}-axis component.}
}
\value{
A data frame with the polar coordinates, the radius (\code{r}) and angle
(\code{theta}) in radians.
}
\description{
Converts a point from Cartesian (\emph{x}, \emph{y}) to polar (\emph{r},
\emph{theta}) coordinates.
}
\examples{
cartesian_to_polar(-1, 0)

}
\seealso{
Other coordinate systems: 
\code{\link{cartesian_degrees}()},
\code{\link{deg_to_rad}()},
\code{\link{euclidean_dist}()},
\code{\link{hausdorff_distance}()},
\code{\link{interpolate_measurements}()},
\code{\link{polar_to_cartesian}()},
\code{\link{quadrant_eye}()},
\code{\link{rad_to_deg}()},
\code{\link{rotate_shape}()},
\code{\link{signed_smallest_angle}()},
\code{\link{smallest_angle}()},
\code{\link{within_2pi}()},
\code{\link{within_360}()}
}
\concept{coordinate systems}
