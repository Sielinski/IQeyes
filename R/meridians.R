## ring_diameters ##

# To find the steepest and flattest meridians on the 3, 5, and 7 mm rings around
# the apex of cornea, create a data frame (ring_diameters) with radial positions
# for every integer angle from 0 - 359° on each of the target rings
ring_diameters <- data.frame(ring_diam = 3, angle = seq(0, 359, 1)) |>
  dplyr::bind_rows(data.frame(ring_diam = 5, angle = seq(0, 359, 1))) |>
  dplyr::bind_rows(data.frame(ring_diam = 7, angle = seq(0, 359, 1)))

# Add Cartesian coordinates to ring_diameters
ring_diameters <- ring_diameters |>
  dplyr::bind_cols(polar_to_cartesian(ring_diameters$ring_diam / 2, ring_diameters$angle)) |>
  dplyr::mutate(ring_diam = round(ring_diam, 5), angle = round(angle, 5))


# Interpolate an exam's measurements data along 3, 5, and 7 mm diameter rings.
# The data frame ring_diameters needs to be created before calling this function
interpolate_rings <- function(source_dat) {

  #source_dat <- sample_curvature

  source_dat <- source_dat |>
    dplyr::filter(!is.na(measurement))

  rings_3d <- ring_diameters

  # target the x and y coordinates of the three ring perimeters
  xo <- sort(unique(rings_3d$x))
  yo <- sort(unique(rings_3d$y))

  # interpolate
  rings_interp <-
    with(source_dat, akima::interp(
      x = x,
      y = y,
      z = measurement,
      linear = F,
      xo = xo,
      yo = yo,
      extrap = T,
      duplicate = 'strip'
    ))

  # extract z-axis values, which is a 2d matrix: we have to find the
  # intersection of each point of interest
  z_dat <- unlist(rings_interp[["z"]])

  # add an empty z-axis
  rings_3d$z <- NA

  # might be able to do this without the for-loop
  for (i in 1:nrow(rings_3d)) {
    #i <- 1
    rings_3d$z[i] <- z_dat[which(xo == rings_3d$x[i]), which(yo == rings_3d$y[i])]
  }

  # remove empty values (although we have extrapolated, so should not be any)
  rings_3d <- rings_3d |>
    dplyr::filter(!is.na(z))

}


#########################
## curvature_meridians
#########################

#' Find an exam's hemi-meridians
#' @description
#' Find the hemi-meridians for an exam at three diameters (3, 5, and 7 mm) and
#' four offset angles (0°, 90°, 180°, and 270°).
#' @param exam_curvature
#' A data frame containing one row for each curvature \code{measurement} and the
#' same columns as \code{sample_curvature}.
#' @param exam_astig
#' A data frame containing one row for each astigmatism \code{measurement} and
#' the same columns as \code{sample_astig}.
#' @param cornea_surface
#' A string, either 'FRONT' or 'BACK', indicating the anterior or posterior
#' surface of the cornea, respectively.
#' @param interp
#' A Boolean. \code{TRUE} to interpolate measurements along the rings.
#' \code{FALSE} is \emph{not} currently implemented.
#' @return
#' A data frame containing one row for each hemi-meridian.
#' @examples
#' curvature_meridians(sample_curvature, sample_astig)
#'
#' @family Meridians
#'
#' @importFrom dplyr cross_join
#' @importFrom dplyr mutate
#' @importFrom dplyr inner_join
#' @importFrom dplyr filter
#' @importFrom dplyr slice_min
#' @importFrom dplyr select
#' @importFrom dplyr slice_max
#' @importFrom dplyr bind_rows
#' @importFrom tidyselect all_of
#'
#' @export
curvature_meridians <- function(exam_curvature, exam_astig, cornea_surface = 'FRONT', interp = T) {

  #exam_curvature <- sample_curvature
  #exam_astig <- sample_astig

  # offset angles are the equivalent of the cardinal directions (N, S, E, W)
  offset_angles <- c(0, 90, 180, 270)

  # window size is the ± region (in degrees) around each offset angle
  window_size <- 45

  # create a data frame that defines a 90° window around each of the offset
  # angles, relative to the axis of astigmatism, for each diameter
  astig_axes <- exam_astig |>
    dplyr::cross_join(data.frame(offset = offset_angles)) |>
    dplyr::mutate(
      meridian = ifelse(offset %in% c(0, 180), 'flat', 'steep'),
      mid = within_360(axs + offset),
      min = mid - window_size,
      max = mid + window_size
    )

  # ensure min is >= 0 and min <= max
  astig_axes <- astig_axes |>
    dplyr::mutate(min = ifelse(min < 0, min + 360, min),
           max = ifelse(min > max, max + 360, max))

  if (interp) {
    # interpolate measurements along the 3, 5, and 7 mm rings
    # the third dimension is the measurement
    rings_3d <- interpolate_rings(exam_curvature)
  } else {
    warning('interp = F is not currently implemented.')
  }

  # for each meridian, ring_diam, and axes, return the range of angles bounded
  # by min and max
  curvature_dat <- astig_axes |>
    dplyr::inner_join(rings_3d, by = 'ring_diam', relationship = "many-to-many") |>
    # ensure that the span of values in the angle column increases from min to max
    dplyr::mutate(angle = ifelse(max > 360, angle + 360, angle)) |>
    dplyr::filter(angle >= min & angle <= max, .by = c(meridian, ring_diam, offset)) |>
    # angle_delta allows slice_min to find the angle closest to the primary axes
    # when two or more measurements have the same value
    dplyr::mutate(angle_delta = abs(angle - mid))

  # find the min radius (i.e., max power) within the axes windows of the two steep
  # meridians for each diameter ring
  power_dat <- curvature_dat |>
    dplyr::filter(meridian == 'steep') |>
    dplyr::slice_min(order_by = data.frame(z, angle_delta), n = 1, by = c(meridian, ring_diam, offset)) |>
    dplyr::select(tidyselect::all_of(join_fields), meridian, offset, ring_diam, angle, mid, x, y, z)

  # find the max radius (i.e., min power) within the axes windows of the two flat
  # meridians for each diameter ring
  power_dat <- curvature_dat |>
    dplyr::filter(meridian == 'flat') |>
    dplyr::slice_max(order_by = data.frame(z, angle_delta), n = 1, by = c(meridian, ring_diam, offset)) |>
    dplyr::select(tidyselect::all_of(join_fields), meridian, offset, ring_diam, angle, mid, x, y, z) |>
    # combine min with max
    dplyr::bind_rows(power_dat)

  # convert curvature radius to power and ensure angle (and axis) is within 360°
  power_dat <- power_dat |>
    dplyr::mutate(measurement = anterior_power(z),
           angle = within_360(angle)) |>
    dplyr::select(tidyselect::all_of(join_fields), meridian, offset, ring_diam, angle, x, y, z, measurement)

  return(power_dat)
}
