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
#' @examples
#' smallest_angle(45, 90)
#' smallest_angle(90, 45)
#'
#' @family Coordinate Systems
#'
#' @export
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
#' @examples
#' signed_smallest_angle(45, 90)
#'
#' @family Coordinate Systems
#'
#' @export
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
#' @family Coordinate Systems
#'
#' @export
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
#' @examples
#' cartesian_to_polar(-1, 0)
#'
#' @family Coordinate Systems
#'
#' @export
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
#' The angle (\emph{theta}) in degrees.
#' @examples
#' cartesian_degrees(2, 2)
#'
#' @family Coordinate Systems
#'
#' @export
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
#' @examples
#' deg_to_rad(180)
#'
#' @family Coordinate Systems
#'
#' @export
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
#' @examples
#' rad_to_deg(pi)
#'
#' @family Coordinate Systems
#'
#' @export
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
#' @family Coordinate Systems
#'
#' @export
euclidean_dist <- function(x, y) {
  sqrt(x ^ 2 +y ^ 2)
}


################
## within_360
################

#' Keep an angle within the 0--360° range
#' @description
#' Converts any angle (in degrees) to the equivalent angle in the 0--360° range.
#' @param angle
#' An angle in degrees.
#' @return
#' The equivalent angle within the 0--360° range.
#' @details
#' Adding or subtracting two angles can sometimes result in a new angle that is
#' outside the 0--360° range. This function converts the result to the equivalent
#' angle within the 0--360° range.
#' @examples
#' within_360(551)
#' within_360(42)
#' within_360(-22)
#'
#' @family Coordinate Systems
#'
#' @export
within_360 <- function(angle) {

  # Normalize the angle to be within 0 - 360 degrees
  normalized_angle <- angle %% 360
  # Handle the case where the angle is negative
  if (normalized_angle < 0) {
    normalized_angle <- normalized_angle + 360
  }
  return(normalized_angle)
}

within_360 <- Vectorize(within_360)


################
## within_2pi
################

#' Keep an angle within the 0--2π range
#' @description
#' Converts any angle (in radians) to the equivalent angle in the 0--2π range.
#' @param angle
#' An angle in radians.
#' @return
#' The equivalent angle within the 0--2π range.
#' @details
#' Adding or subtracting two angles can sometimes result in a new angle that is
#' outside the 0--2π range. This function converts the result to the equivalent
#' angle within the 0--2π range.
#' @examples
#' within_2pi(2 * pi)
#' within_2pi(pi)
#' within_2pi(-pi)
#'
#' @family Coordinate Systems
#'
#' @export
within_2pi <- function(angle) {

  # Normalize the angle to be within 0 - 360 degrees
  normalized_angle <- angle %% (2 * pi)
  # Handle the case where the angle is negative
  if (normalized_angle < 0) {
    normalized_angle <- normalized_angle + (2 * pi)
  }
  return(normalized_angle)
}

within_2pi <- Vectorize(within_2pi)


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
#' An \emph{n}x2 matrix containing the points of original matrix rotated by
#' \emph{theta} degrees.
#' @examples
#' matrix(rep(seq(1:5), 2), ncol = 2) |> rotate_shape(180)
#'
#' @family Coordinate Systems
#'
#' @export
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
#'
#' @description
#' Interpolates measurements over an expanded
#' grid of Cartesian coordinates (\emph{x} and \emph{y}).
#'
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
#' interpolate_measurements(sample_curvature) |>
#'    head()
#'
#' @family Coordinate Systems
#'
#' @importFrom dplyr filter
#' @importFrom akima interp
#'
#' @export
interpolate_measurements <- function(source_dat, ...) {

  #source_dat <- sample_curvature

  source_dat <- source_dat |>
    dplyr::filter(!is.na(measurement))

  # Create a grid that spans the extents of the measured x and y axes
  #x_range <- with(source_dat, seq(min(x), max(x), length.out = length(unique(x))))
  #y_range <- with(source_dat, seq(min(y), max(y), length.out = length(unique(y))))
  y_range <- x_range <- seq(-4.5, 4.5, 0.1)
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


##################
## quadrant_eye
##################

#' Quadrant of an eye
#'
#' @description
#' Returns the quadrant of a point based on its \emph{x}- and \emph{y}-axis
#' position and the \code{exam_eye}.

#' @param x
#' The x-axis coordinate.
#' @param y
#' The y-axis coordinate.
#' @param exam_eye
#' A character vector containing either "left" or "right".
#' @return
#' A two-character string. See details.

#' @details
#' The quadrant system for an eye is split on the \emph{y}-axis by the superior
#' ("S") and the inferior ("I"), and on the \emph{x}-axis by the nasal ("N") and
#' the temporal ("T"). The superior is always \code{y > 0}, and the inferior is
#' always \code{y ≤ 0}. The nasal and temporal depend on the eye: For the right
#' eye, the nasal is \code{x > 0}; for the left eye, the nasal is \code{x ≤ 0}.

#' @examples
#' quadrant_eye(1, 1, 'left')
#' quadrant_eye(1, 1, 'right')
#'
#' @family Coordinate Systems
#'
#' @importFrom dplyr filter
#'
#' @export
quadrant_eye <- function(x, y, exam_eye) {

  if (y > 0) {
    # superior
    if (exam_eye == 'right') {
      quadrant <- dplyr::if_else(x > 0, 'SN', 'ST')
    } else {
      # left
      quadrant <- dplyr::if_else(x >= 0, 'ST', 'SN')
    }
  } else {
    # inferior
    if (exam_eye == 'right') {
      quadrant <- dplyr::if_else(x > 0, 'IN', 'IT')
    } else {
      # left
      quadrant <- dplyr::if_else(x >= 0, 'IT', 'IN')
    }
  }

  return(quadrant)
}

quadrant_eye <- Vectorize(quadrant_eye)


########################
## hausdorff_distance
########################

#' Hausdorff distance
#' @description
#' Calculates the Hausdorff distance between two contours.
#' @param contour_A
#' A data frame containing contour A.
#' @param contour_B
#' A data frame containing contour B.
#' @return
#' The Hausdorff distance between the two contours.
#' @examples
#' hausdorff_distance(
#'   contour_A = get_contour(sample_curvature, contour_power = 44),
#'   contour_B = get_contour(sample_curvature, contour_power = 44.5)
#' )
#'
#' @family Coordinate Systems
#'
#' @importFrom pracma hausdorff_dist
#' @importFrom dplyr select
#'
#' @export
hausdorff_distance <- function(contour_A, contour_B) {

  matrix_A <- contour_A |>
    dplyr::select(x, y) |>
    unique() |>
    as.matrix()

  matrix_B <- contour_B |>
    dplyr::select(x, y) |>
    unique() |>
    as.matrix()

  # The Hausdorff distance is the maximum of these two maximum distances
  pracma::hausdorff_dist(matrix_A, matrix_B)
}
