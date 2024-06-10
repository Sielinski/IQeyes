########################
## Coordinate Systems
########################

###########################
## smallest_angle
###########################

#' The smallest angle (in degrees) between two radial lines
#' @description
#' Calculates the smallest angle (in degrees) between two radial lines.
#' @param angle_1
#' The angle (in degrees) of the first radial line.
#' @param angle_2
#' The angle (in degrees) of the second radial line.
#' @return
#' The angle between with the two lines.
#' @details
#' The result is essentially the absolute difference between two angles.
#' @seealso
#' [IQeyes::signed_smallest_angle()]
#' #' @examples
#' smallest_angle(45, 90)
#' smallest_angle(90, 45)
#'
smallest_angle <- function(angle_1, angle_2) {
  if (is.na(angle_1) || is.na(angle_2)) return(NA)

  diff <- abs(angle_1 - angle_2)
  diff <- diff %% 360
  if (diff > 180) diff <- 360 - diff
  return(diff)
}

smallest_angle <- Vectorize(smallest_angle)


###########################
## signed_smallest_angle
###########################

#' The smallest signed angle (in degrees) between two radial lines
#' @description
#' Calculates the smallest angle (in degrees) between two radial lines. If
#' \code{angle_1} > \code{angle_2}, the result is negative. If
#' \code{angle_1} < \code{angle_2}, the result is positive.
#' @param angle_1
#' The angle (in degrees) of the first radial line.
#' @param angle_2
#' The angle (in degrees) of the second radial line.
#' @return
#' The angle (in degrees) between with the two lines.
#' @details
#' The result is essentially the signed difference between two angles.
#' @seealso
#' [IQeyes::smallest_angle()]
#' @examples
#' signed_smallest_angle(45, 90)
#'
signed_smallest_angle <- function(angle_1, angle_2) {
  if (is.na(angle_1) || is.na(angle_2)) return(NA)

  raw_diff <- angle_2 - angle_1
  diff <- abs(raw_diff) %% 360
  if (diff > 180) {
    diff <- 360 - diff
    if (raw_diff > 0) diff <- -diff
  } else {
    if (raw_diff < 0) diff <- -diff
  }

  return(diff)
}

signed_smallest_angle <- Vectorize(signed_smallest_angle)


########################
## polar_to_cartesian
########################

#' Convert the polar coordinates of a point to Cartesian
#' @description
#' Converts a point from polar (\emph{radius}, \emph{theta}) to
#' Cartesian (\emph{x}, \emph{y}) coordinates.
#' @param radius
#' The radius of the polar coordinates.
#' @param angle_deg
#' The angle (in degrees) of the polar coordinates.
#' @return
#' A data frame with the Cartesian coordinates, \code{x} and \code{y}.
#' @examples
#' polar_to_cartesian(1, 45)
#'
polar_to_cartesian <- function(radius, angle_deg) {
  # Convert the angle from degrees to radians
  angle_rad <- angle_deg * (pi / 180)

  # Compute the x and y coordinates
  x <- radius * cos(angle_rad)
  y <- radius * sin(angle_rad)

  return(data.frame(x = x, y = y))
}


########################
## cartesian_to_polar
########################

#' Convert the Cartesian coordinates of a point to polar coordinates
#' @description
#' Converts a point from Cartesian (\emph{x}, \emph{y}) to polar (\emph{r},
#' \emph{theta}) coordinates.
#' @param x
#' The \emph{x}-axis component.
#' @param y
#' The \emph{y}-axis component.
#' @return
#' A data frame with the polar coordinates, the radius (\code{r}) and angle
#' (\code{theta}) in radians.
#' @seealso
#' [IQeyes::cartesian_degrees()]
#' @examples
#' cartesian_to_polar(-1, 0)
#'
cartesian_to_polar <- function(x, y) {
  r <- sqrt(x ^ 2 + y ^ 2) # Calculate radial distance
  theta <- atan2(y, x) # Calculate angle in radians

  return(data.frame(r = r, theta = theta))
}


#######################
## cartesian_degrees
#######################

#' Convert the Cartesian coordinates of a point to an angle (in degrees)
#' @description
#' Calculates the angle (in degrees) formed by the \emph{x}- and
#' \emph{y}-coordinates of a point on a Cartesian plane.
#' @param x
#' The \emph{x}-axis component.
#' @param y
#' The \emph{y}-axis component.
#' @return
#' the angle (\emph{theta}) in degrees.
#' @seealso
#' [IQeyes::cartesian_to_polar()]
#' @examples
#' cartesian_degrees(2, 2)
cartesian_degrees <- function(x, y) {
  theta <- atan2(y, x) # Calculate angle in radians

  # Convert angle from radians to degrees
  theta_deg <- rad_to_deg(theta) |>
    within_360()

  return(theta_deg)
}


################
## deg_to_rad
################

#' Convert an angle from degrees to radians
#' @description
#' Convert an angle from degrees to radians.
#' @param angle_deg
#' The angle in degrees.
#' @return
#' The angle in radians.
#' @seealso
#' [IQeyes::rad_to_deg()]
#' @examples
#' deg_to_rad(180)
#'
deg_to_rad <- function(angle_deg) angle_deg * pi / 180


################
## rad_to_deg
################

#' Convert an angle from radians to degrees
#' @description
#' Convert an angle from radians to degrees.
#' @param angle_rad
#' The angle in radians.
#' @return
#' The angle in degrees.
#' @seealso
#' [IQeyes::deg_to_rad()]
#' @examples
#' rad_to_deg(pi)
#'
rad_to_deg <- function(angle_rad) angle_rad * 180 / pi


####################
## euclidean_dist
####################

#' Euclidean distance
#' @description
#' Calculates the Euclidean distance of a point (\emph{x}, \emph{y}) from the
#' origin (0, 0).
#' @param x
#' The \emph{x}-axis component.
#' @param y
#' The \emph{y}-axis component.
#' @return
#' The distance of the point from the origin.
#' @examples
#' euclidean_dist(2, 2)
#'
euclidean_dist <- function(x, y) {
  sqrt(x ^ 2 +y ^ 2)
}


################
## within_360
################

#' Keep an angle within the 0--360° range
#' @description
#' Converts an angle from a -360--720° range to a 0--360° range.
#' @param x
#' An angle within the -360--720° range.
#' @return
#' The same angle within the 0--360° range.
#' @details
#' Adding or subtracting two angles can sometimes result in a new angle that is
#' outside the 0--360° range. This function converts the result to the same
#' angle within the 0--360° range.
#' @examples
#' within_360(551)
#' within_360(42)
#' within_360(-22)
within_360 <- function(x) {
  if (x >= 360) {
    x - 360
  } else if(x < 0) {
    x + 360
  } else {
    x
  }
}

within_360 <- Vectorize(within_360)


##################
## rotate_shape
##################

#' Rotate a shape
#' @description
#' Rotate a shape (i.e., a matrix of points) by a specified angle (in degrees).
#' @param points
#' An \emph{n} x 2 matrix of Cartesian coordinates: \emph{n} rows, one for each
#' point, and two columns, the \emph{x} and \emph{y} coordinates of a point,
#' respectively.
#' @param theta_deg
#' Angle (in degrees) to rotate the shape.
#' @return
#' An \emph{n}x2 matrix containing the points of orininal matrix rotated by
#' \empha{theta} degrees.
#' @examples
#' matrix(rep(seq(1:5), 2), ncol = 2) |> rotate_shape(180)
#'
rotate_shape <- function(points, theta_deg) {
  theta_rad <- deg_to_rad(theta_deg)

  # Create the rotation matrix
  R <- matrix(c(cos(theta_rad),
                -sin(theta_rad),
                sin(theta_rad),
                cos(theta_rad)
  ), ncol = 2)


  return(points %*% t(R))
}


#############################
## interpolate_measurements
#############################

#' Interpolate measurements over an expanded grid
#' @description
#' Interpolates measurements over an expanded
#' grid of Cartesian coordinates (\emph{x} and \emph{y}).
#' @param source_dat
#' A data frame containing one row for each \code{measurement}. Each row must contain
#' at least three columns: a set of Cartesian coordinates (\code{x} and
#' \code{y}) and a \code{measurement}.
#' @param ...
#' Additional parameters that get passed to \code{akima::interp()}
#' @return
#' A new data frame with the interpolated measurements.
#' @details
#' The \emph{x} and \emph{y} axes of the expanded grid will contain the same
#' number of unique values as original axes, but they will be equally spaced
#' along those axes.
#'
#' The interpolated measurements will \emph{not} necessarily include the
#' original measurements.
#' @examples
#' interpolate_measurements(sample_curvature[1:10, c('x', 'y', 'measurement')])
#'
interpolate_measurements <- function(source_dat, ...) {

  #source_dat <- sample_curvature

  # Create a grid that spans the extents of the measured x and y axes
  x_range <- with(source_dat, seq(min(x), max(x), length.out = length(unique(x))))
  y_range <- with(source_dat, seq(min(y), max(y), length.out = length(unique(y))))
  grid <- expand.grid(x = x_range, y = y_range)

  # Interpolate z values on the grid
  interpolated_dat <-
    with(source_dat,
         akima::interp(
           x,
           y,
           z = measurement,
           xo = x_range,
           yo = y_range,
           duplicate = 'strip'
         ))
  grid$z <- as.vector(interpolated_dat$z)

  grid <- grid |>
    dplyr::filter(!is.na(z))

  return(grid)
}
