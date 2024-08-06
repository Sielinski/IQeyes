# generate radial arms
generate_radial_arms <- function(radius, grain, angles) {
  all_points <- data.frame()

  for (theta in angles) {
    theta_rad <- theta * (pi / 180)  # Convert angle to radians

    # Generate points along the radial arm at intervals of 0.1
    r_values <- seq(0, radius, by = grain)
    x_values <- (r_values * cos(theta_rad)) |> round(1)
    y_values <- (r_values * sin(theta_rad)) |> round(1)

    # Create a data frame for the current radial arm
    radial_arm <- data.frame(x = x_values,
                             y = y_values,
                             r = r_values,
                             theta = rep(theta, length(r_values)))

    # Combine with all points
    all_points <- rbind(all_points, radial_arm)
  }

  return(all_points)
}

# define the extent and grain of the radial arm(s)
radius <- 7
grain <- 0.1
# no need to go below 1°: that resolution is already more granular than the
# the matrix of measurements can support
angles <- seq(0, 359, by = 1)

# Generate the data frame
radial_arms <- generate_radial_arms(radius, grain, angles)


#########################
## curvature_meridians
#########################

#' Find an exam's semi-meridians
#' @description
#' Find the semi-meridians for an exam at three diameters (3, 5, and 7 mm).
#'
#' @param exam_curvature
#' A data frame with the same structure as [IQeyes::sample_curvature],
#' containing one row for each curvature (i.e., radius of curvature)
#' \code{measurement}.
#' @param exam_astig
#' A data frame having the same columns as [IQeyes::sample_astig] and containing
#' one row for each of the three ring diameters.
#' @param cornea_surface
#' A string, either 'FRONT' or 'BACK', indicating the anterior or posterior
#' surface of the cornea, respectively.
#' @param from_vertex
#' A Boolean. \code{TRUE} to use the vertex as the starting point for the
#' semi-meridians. \code{FALSE} (the default), to use the end point of the
#' prior semi-meridian. (The 3 mm semi-meridian always uses the vertex as its
#' starting point.)
#' @param interp
#' A Boolean. \code{TRUE} to interpolate measurements.
#'
#' @return
#' A data frame containing one row for each semi-meridian.
#'
#' \describe{
#'  \item{ring}{Ring diameter of the semi-meridian.}
#'  \item{radius}{Ending radius of the semi-meridian segment.}
#'  \item{starting_radius}{Starting radius of the semi-meridian segment.}
#'  \item{axis}{Either the \code{flat} or \code{steep} axis of astigmatism,
#'  identifying the semi-meridian's reference axis.}
#'  \item{between_min_max}{A Boolean identifying whether the semi-meridian was
#'  found within a 180° window "above" or "below" the opposite angle of
#'  astigmatism. See details.}
#'  \item{theta}{Offset angle (in degrees) from the 3 o'clock position.}
#'  \item{power}{Dioptric power of of the semi-meridian. When actual
#'  measurement data are used, \code{power} will be the power at the end radius
#'  of the semi-meridian segment. When interpreted data are used, \code{power}
#'  will be the average power of the semi-meridian segment. }
#' }
#'
#' @details
#' The Pentacam provides the flat axis of astigmatism for the 3, 5, and 7 mm
#' rings. This function looks "above" and "below" each of these axes for the
#' steepest semi-meridians (i.e., the radial arms with the highest average
#' power within these 180° windows.)
#'
#' Presuming that the steep axis of astigmatism is perpendicular to the flat
#' axis, this function likewise looks "above" and "below" the perpendicular axes
#' for the flattest semi-meridians (i.e., the radial arms with the lowest
#' average power).
#'
#' Note that multiple semi-meridians can have the same power. When that happens,
#' this function will return the radial arm at the mid-point of the largest
#' cluster of radial arms with the same power.
#'
#' @examples
#' curvature_meridians(sample_curvature, sample_astig) |>
#'   dplyr::select(-tidyselect::all_of(join_fields))
#'
#' @family meridians
#'
#' @importFrom dplyr mutate
#' @importFrom dplyr inner_join
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom dplyr rename
#' @importFrom dplyr if_else
#' @importFrom dplyr cross_join
#' @importFrom dplyr summarize
#' @importFrom dplyr bind_rows
#' @importFrom dplyr group_by
#' @importFrom dplyr group_split
#' @importFrom dplyr ungroup
#' @importFrom tidyr pivot_longer
#' @importFrom tidyselect all_of
#'
#' @importFrom dplyr slice_min
#' @importFrom dplyr slice_max
#'
#' @export
curvature_meridians <- function(exam_curvature,
                                exam_astig,
                                cornea_surface = 'FRONT',
                                from_vertex = F,
                                interp = T) {

  #exam_curvature <- sample_curvature
  #exam_astig <- sample_astig

  exam_record <- exam_curvature[1, ] |>
    dplyr::select(tidyselect::all_of(join_fields))

  if (interp) {
    # interpolate measurement
    interp_z <- exam_curvature |>
      interpolate_measurements() |>
      dplyr::mutate(surface = cornea_surface)

    exam_curvature <- exam_record  |>
      dplyr::cross_join(interp_z) |>
      dplyr::rename(measurement = z)
  }

  # calculate the power of the points that comprise each radial arm
  radial_power <- radial_arms |>
    dplyr::inner_join(exam_curvature, by = c('x', 'y')) |>
    dplyr::mutate(power = anterior_power(measurement))

  # confirm that we have good astig data
  if (nrow(exam_astig) != 3 || anyNA(exam_astig$axs)) return(invisible(exam_record))

  # create the astig windows: goal is to look above and below both the flat and
  # steep axes (as determined by the Pentacam at the 3, 5, and 7 mm diameters)
  astig_windows <- exam_astig |>
    dplyr::select(-all_of(join_fields)) |>
    # by definition, axs is the flat axis
    dplyr::rename(flat = axs) |>
    # get the radii
    dplyr::mutate(radius = ring_diam / 2,
                  # calculate the starting radii for these ranges (basically, 1 mm closer
                  # except for the 3 mm ring, which starts at 0)
                  starting_radius = radius - dplyr::if_else(radius == 1.5, radius, 1),
                  # the steep axis is perpendicular
                  steep = within_360(flat - 90)) |>
    # if from_vertex is set to T, then use the vertex (0) as the starting point
    # of the radial arm
    dplyr::mutate(starting_radius = (if (from_vertex) 0 else starting_radius)) |>
    tidyr::pivot_longer(cols = c(flat, steep), names_to = 'axis', values_to = 'angle') |>
    # calculate the angle opposite to the axis
    dplyr::mutate(opposite = within_360(angle + 180),
           crosses_zero = angle > opposite)

  # get a set of radial arms for each of the three diameters
  radial_windows <- radial_power |>
    dplyr::cross_join(astig_windows) |>
    # Use between_min_max to define the hemispheres: One hemisphere will contain
    # arms between the min and max, and the other hemisphere will contain arms
    # that are *not* between min and max. However, when the range of angles
    # crosses 0°, we need to adjust the logic, using a combination of
    # crosses_zero and !between_max_min to create the same effect
    dplyr::mutate(between_min_max = theta >= angle & theta < opposite,
           between_max_min = theta >= opposite & theta < angle) |>
    dplyr::mutate(between_min_max = between_min_max | (crosses_zero & !between_max_min)) |>
    # only need segments of arms between the starting and stopping radii for
    # each ring
    dplyr::filter(r >= starting_radius & r <= radius)

  # calculate the mean power of each segment for each range (axis, flip_flop)
  arm_power <- radial_windows |>
    dplyr::summarize(mean_power = mean(power), .by = c(theta, axis, between_min_max, radius))

  # find the min and max arm power (just the power) for each range
  min_max_power <- arm_power |>
    dplyr::summarize(min_power = min(mean_power),
              max_power = max(mean_power),
              .by = c(axis, between_min_max, radius))

  # find the radial arms that match these power values
  candidates_flat <- arm_power |>
    dplyr::inner_join(min_max_power, by = c('axis', 'between_min_max', 'radius', 'mean_power' = 'min_power')) |>
    # the flat semi-meridians will be on the either side of the steep axis
    dplyr::filter(axis == 'steep') |>
    dplyr::mutate(axis = 'flat') |>
    dplyr::select(radius, axis, between_min_max, theta, mean_power)

  candidates_steep <- arm_power |>
    dplyr::inner_join(min_max_power, by = c('axis', 'between_min_max', 'radius', 'mean_power' = 'max_power')) |>
    # the steep semi-meridians will be on the either side of the flat axis
    dplyr::filter(axis == 'flat') |>
    dplyr::mutate(axis = 'steep') |>
    dplyr::select(radius, axis, between_min_max, theta, mean_power)

  candidate_angles <- dplyr::bind_rows(candidates_steep, candidates_flat)

  # cluster candidate angles
  cluster_data <- function(data, group, threshold = 1) {

    if (length(data) == 0) {
      return(data.frame(group = group,
                        cluster = NA))
    }

    # Step 1: Sort the data
    data <- sort(data)
    group <- unique(group)

    # Step 2: Initialize clusters
    clusters <- list()
    current_cluster <- data[1]

    if (length(data) == 1) {
      return(data.frame(group = group,
                        cluster = current_cluster))

    } else {
      # Step 3: Iterate through the data
      for (i in 2:length(data)) {
        if (data[i] - data[i - 1] <= threshold) {
          current_cluster <- c(current_cluster, data[i])
        } else {
          clusters <- c(clusters, list(current_cluster))
          current_cluster <- data[i]
        }
      }
      clusters <- c(clusters, list(current_cluster)) # Add the last cluster

      # Step 4: Prioritize clusters by proximity to other numbers
      min_distance <- function(cluster, all_data) {
        min(sapply(all_data, function(x) min(abs(x - cluster))))
      }

      clusters <- clusters[order(sapply(clusters, min_distance, all_data = data))]

      cluster_mean <- mean(clusters[[1]]) |> round(0)

      return(data.frame(group = group,
                        cluster = cluster_mean))
    }
  }

  # use simple clustering to identify the preferred angles
  preferred_angles <- candidate_angles |>
    dplyr::group_by(radius, axis, between_min_max) |>
    dplyr::group_split() |>
    lapply(function(x) cluster_data(x$theta,
                                    group = paste0(x$radius, '_', x$axis, '_', x$between_min_max))) |>
    dplyr::bind_rows() |>
    dplyr::ungroup()

  if (interp) {
    # if data are interpolated, there's a good chance that we won't have the
    # power at the ends of the arm segments, so use the arms' averages instead
    radial_power <- arm_power |>
      dplyr::select(-axis, -between_min_max) |>
      unique() |>
      dplyr::rename(r = radius,
                    power = mean_power)
  }

  # define the semi-meridians
  meridian_angles <- candidate_angles |>
    dplyr::mutate(group = paste0(radius, '_', axis, '_', between_min_max)) |>
    dplyr::inner_join(preferred_angles, by = c('group', 'theta' = 'cluster')) |>
    # get the power at the ends of the arm segments at the preferred angles
    dplyr::inner_join(radial_power, by = c('theta', 'radius' = 'r')) |>
    # get the starting_radius
    dplyr::inner_join(astig_windows, by = c('radius', 'axis')) |>
    dplyr::mutate(power = round(power, 1)) |>
    dplyr::mutate(ring = radius * 2) |>
    dplyr::select(ring, radius, starting_radius, axis, between_min_max, theta, power) |>
    dplyr::arrange(ring, axis)

  exam_record |>
    dplyr::cross_join(meridian_angles)

}


##########
## srax
##########

#' Calculate the skewed radial axis (SRAX) for an exam
#'
#' @description
#' The SRAX is 180 minus the smallest angle between the steepest radial axes
#' "above" and "below" the flat axis of astigmatism.
#'
#' @param exam_meridians
#' A data frame with the same structure as the output of
#' [IQeyes::curvature_meridians].
#' @param radius
#' A number identifying the radius of the axes to use for the skew calculation.
#'
#' @return
#' A data frame containing the \code{join_fields} and \code{srax}, the skew of
#' radial axes.
#'
#' @details
#' The formal definition of SRAX uses the horizontal meridian to split the
#' cornea into hemispheres. Presuming that [IQeyes::curvature_meridians] is used
#' to identify the steepest radial axes, this function uses the flat axis of
#' astigmatism identified by the Pentacam at the 3, 5, or 7 mm ring.
#'
#' @references
#' Rabinowitz, Yaron S. "Videokeratographic Indices to Aid in Screening for
#' Keratoconus." \emph{Journal of Refractive Surgery}, 1995 Sep-Oct;11(5):371-9.
#' \url{https://doi.org/10.3928/1081-597X-19950901-14}.
#'
#' @examples
#' curvature_meridians(sample_curvature, sample_astig) |>
#'   srax()
#'
#' @family meridians
#'
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom dplyr bind_rows
#' @importFrom tidyselect all_of
#'
#' @export
srax <- function(exam_meridians, radius = 1.5) {

  exam_record <- exam_meridians[1, ] |>
    dplyr::select(tidyselect::all_of(join_fields))

  target_meridians <- exam_meridians |>
    dplyr::filter(radius == radius, axis == 'steep')

  srax <- 180 - smallest_angle(target_meridians$theta[1], target_meridians$theta[2])

  return(dplyr::bind_cols(exam_record, srax = srax))
}
