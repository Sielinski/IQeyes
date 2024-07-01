##################
## scale_rotate
##################

#' Scale and rotate
#' @description
#' Scales and rotates a shape defined by a set of \emph{x} and \emph{y}
#' coordinates. Can be either a curvature (e.g., [IQeyes::sample_curvature]) or
#' contour data frame (e.g., [IQeyes::sample_contour]).
#' @param shape
#' A data frame containing one row for each point that defines a shape (e.g.,
#' a silhouette or contour). Currently, \code{scale_rotate()} expects the
#' data frame to contain the join fields.
#' @param axs
#' A number. Oftentimes, this should be the angle (in degrees) of the primary
#' axis of astigmatism for an exam.
#' @param r_target
#' A number defining the target radial length of the scaled shape.
#' @return
#' A data frame containing the same columns as \code{shape} with transformed
#' \emph{x} and \emph{y} coordinates. No other values are changed.
#' @details
#' The angle of rotation is determined by the primary axis of astigmatism for
#' the exam. If \code{axs} > 90, its value is shifted by -180°, ensuring that
#' the contour isn't flipped along the horizontal axis.
#'
#' The scaling factor is determined by the max radial length of any point on the
#' shape, not the end-to-end length of the shape, which is especially
#' important for asymmetric shapes.
#' @examples
#' scale_rotate(sample_contour, axs = 33.6) |>
#'   head()
#' scale_rotate(sample_curvature, axs = 33.6) |>
#'   head()
#'
#' @family Contours
#'
#' @importFrom dplyr select
#' @importFrom dplyr inner_join
#' @importFrom dplyr mutate
#' @importFrom dplyr rename
#' @importFrom dplyr group_by
#' @importFrom dplyr ungroup
#' @importFrom dplyr across
#' @importFrom dplyr summarize
#' @importFrom tidyselect all_of
#'
#' @export
scale_rotate <- function(shape, axs = 0, r_target = 4) {

  axs <- ifelse(axs > 90, axs - 180, axs)

  shape_rotated <- shape |>
    dplyr::mutate(angle = cartesian_degrees(x, y), r = euclidean_dist(x, y)) |>
    dplyr::mutate(angle_new = within_360(angle - axs),
           x_new = polar_to_cartesian(r, angle_new)$x,
           y_new = polar_to_cartesian(r, angle_new)$y,
    ) |>
    dplyr::rename(x_orig = x,
           y_orig = y,
           angle_orig = angle,
           x = x_new,
           y = y_new,
           angle = angle_new)

  # calculate the scaling factor
  # can't use the full extent of the contour because, for example, asymmetric
  # shapes might cover the full extent of just one half of the cornea, so look
  # at the max extent of the cornea in the top and bottom halves of the map
  scaling_factors <- shape_rotated |>
    # include join fields to preserve them in the output
    dplyr::group_by(dplyr::across(tidyselect::all_of(join_fields))) |>
    dplyr::summarize(r_max = max(r), .groups = 'keep') |>
    dplyr::mutate(scale_factor = r_target / r_max) |>
    dplyr::ungroup()

  # scale the rotated contours
  shape_transformed <- shape_rotated |>
    dplyr::inner_join(scaling_factors, by = join_fields) |>
    dplyr::mutate(r_new = r * scale_factor) |>
    dplyr::mutate(x_new = polar_to_cartesian(r_new, angle)$x,
           y_new = polar_to_cartesian(r_new, angle)$y,
    ) |>
    dplyr::rename(x_old = x,
           y_old = y,
           r_old = r,
           x = x_new,
           y = y_new,
           r = r_new
    )

  shape_transformed |>
    dplyr::select(tidyselect::all_of(colnames(shape)))

}


##################
## fill_contour
##################

#' Fill a contour with points
#' @description
#' Fills a contour with points, creating a silhouette of its shape.
#' @param exam_curvature
#' A data frame containing one row for each curvature \code{measurement} and the
#' same columns as [IQeyes::sample_curvature].
#' @param contour_power
#' A number defining the dioptric power of the contour to fill.
#' @param interp
#' A Boolean. \code{TRUE} to interpolate points within the contour.
#' @param scale_rotate
#' A Boolean. \code{TRUE} to scale and rotate the contour.
#' @param ...
#' Additional parameters to be passed to \code{scale_rotate()}, specifically
#' \code{axs} and \code{r_target}.
#' @return
#' A data frame with the same columns as [IQeyes::sample_curvature] containing
#' the points of the filled and (if \code{scale_rotate == T}) transformed
#' contour. The data frame also includes \code{pct_coverage}, a column
#' indicating the percentage of the cornea surface contained within the contour.
#' @details
#' Currently only handles the anterior (front) surface of the cornea.
#'
#' Because the function creates data (through interpolation), it cannot return
#' any additional columns that may have been included as part of the incoming
#' \code{exam_curvature}.
#' @examples
#' fill_contour(sample_curvature, contour_power = 45.5, axs = 33.6) |>
#'   dplyr::select(-all_of(join_fields)) |>
#'   head()
#'
#' @family Contours
#'
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom dplyr rename
#' @importFrom dplyr mutate
#' @importFrom dplyr cross_join
#' @importFrom tidyselect all_of
#'
#' @export
fill_contour <- function(exam_curvature,
                         contour_power,
                         interp = T,
                         scale_rotate = T,
                         ...
                         ){

  # cleanse the curvature data (currently limited to the front surface)
  curvature_dat <- exam_curvature |>
    dplyr::filter(surface == 'FRONT') |>
    dplyr::filter(!is.na(measurement))

  if (interp) {
    # Interpolation invalidates column_name and requires recalculation of the
    # power and polar coordinates

    # preserve the exam details (i.e., join fields and surface)
    exam_record <- exam_curvature |>
      dplyr::select(tidyselect::all_of(join_fields), surface) |>
      unique()

    # interpolate points
    interp_dat <- curvature_dat |>
      interpolate_measurements() |>
      dplyr::rename(measurement = z)

    # combine interpolated data with exam_record's join fields (replacing the
    # original curvature_dat)
    curvature_dat <- exam_record |>
      dplyr::cross_join(interp_dat)

    # calculate the the dioptric power of the interpolated measurement and the
    # polar coordinates of the location
    curvature_dat <- curvature_dat |>
      dplyr::mutate(ring_diam = cartesian_to_polar(x, y)$r,
                    angle = cartesian_to_polar(x, y)$theta,
                    column_name = 'interpolated'
                    )
  }

  # calculate dioptric power
  curvature_dat <- curvature_dat |>
    dplyr::mutate(power = anterior_power(measurement))

  n_points_all <- nrow(curvature_dat)

  # any point with a power ≥ the contour will be inside the contour
  curvature_dat <- curvature_dat |>
    dplyr::filter(power >= contour_power)

  pct_coverage <- (nrow(curvature_dat) / n_points_all)

  # scale and rotate (*after* data frame has been reduced to the target contour)
  if (scale_rotate) curvature_dat <- scale_rotate(curvature_dat, ...)

  # put columns back in the same order as a typical curvature data frame
  curvature_dat |>
    dplyr::select(dplyr::all_of(colnames(sample_curvature))) |>
    mutate(pct_coverage = pct_coverage)

}


#####################
## contour_polygon
#####################

#' Convert a contour to a closed polygon
#' @description
#' Converts every segment of a contour to a closed polygon. The difference
#' between a \emph{segment} and a \emph{polygon} is that the latter contains two
#' points at the same location, representing the starting and ending points of
#' the polygon.
#' @param exam_contour
#' A data frame containing one row for each point of the contour and the same
#' columns as [IQeyes::sample_contour]. If the contour needs to be scaled and/or
#' rotated, the necessary transforms must be performed prior to calling
#' \code{contour_polygon}.
#' @return
#' A data frame containing the same columns as \code{exam_contour} with one
#' extra row for each segment of the contour.
#' @details
#' A contour is comprised of one or more segments, uniquely identified by
#' \code{segment_id}. \code{contour_polygon()} takes the first point of each
#' segment (based on its ordinal position in the original data frame) and
#' appends a replica to the end of the the original data frame.
#' @examples
#' contour_polygon(sample_contour) |>
#'    head()
#'
#' @family Contours
#'
#' @importFrom dplyr select
#' @importFrom dplyr group_by
#' @importFrom dplyr across
#' @importFrom dplyr ungroup
#' @importFrom dplyr mutate
#' @importFrom dplyr row_number
#' @importFrom dplyr summarize
#' @importFrom dplyr filter
#' @importFrom dplyr inner_join
#' @importFrom dplyr bind_rows
#' @importFrom dplyr arrange
#' @importFrom tidyselect all_of
#'
#' @export
contour_polygon <- function(exam_contour) {

  # assign a number to every row of every segment
  contour_rn <- exam_contour |>
    dplyr::group_by(dplyr::across(tidyselect::all_of(join_fields)), segment_id) |>
    dplyr::mutate(rn = dplyr::row_number()) |>
    dplyr::select(tidyselect::all_of(colnames(exam_contour)), rn) |>
    dplyr::ungroup()

  # find the first and last row of every segment
  contour_limits <- contour_rn |>
    dplyr::group_by(dplyr::across(tidyselect::all_of(join_fields)), segment_id) |>
    dplyr::summarize(min_rn = min(rn), max_rn = max(rn)) |>
    dplyr::ungroup()

  # connect the ends of every segment by adding a duplicate of the first
  # point after the last point
  contour_polygon <- contour_rn |>
    dplyr::filter(rn == 1) |>
    dplyr::inner_join(contour_limits, by = c(join_fields, 'segment_id')) |>
    dplyr::mutate(rn = max_rn + 1) |>
    dplyr::select(colnames(contour_rn)) |>
    dplyr::bind_rows(contour_rn) |>
    dplyr::arrange(dplyr::across(tidyselect::all_of(join_fields)), segment_id, rn)

  contour_polygon |>
    dplyr::select(tidyselect::all_of(colnames(exam_contour)))

}


########################
## get_contour_levels
########################

#' Get contour levels
#' @description
#' Returns the full range of contours and the number of segments needed to plot
#' \code{exam_curvature} on a curvature map.
#' @param exam_curvature
#' A data frame containing one row for each curvature \code{measurement} and the
#' same columns as [IQeyes::sample_curvature].
#' @return
#' A table. The table's labels identify the contour levels, and the values in
#' the table identify the number of segments that comprise the contour.
#' @details
#' This function doesn't currently work for posterior surfaces because
#' [IQeyes::plot_scale] presumes the color scale is on the "absolute
#' scale", which is designed for the dioptric power of anterior surfaces.
#' @examples
#' get_contour_levels(sample_curvature)
#'
#' @family Contours
#'
#' @importFrom dplyr mutate
#' @importFrom dplyr filter
#' @importFrom akima interp
#'
#' @export
get_contour_levels <- function(exam_curvature) {

  if (nrow(exam_curvature) == 0) warning('No data in exam_curvature.')

  # calculate power
  exam_curvature <- exam_curvature |>
    dplyr::mutate(z = anterior_power(measurement)) |>
    dplyr::filter(!is.na(z))

  # Create a grid that spans the extents of the measured x and y axes
  x_range <- with(exam_curvature, seq(min(x), max(x), length.out = length(unique(x))))
  y_range <- with(exam_curvature, seq(min(y), max(y), length.out = length(unique(y))))

  # Interpolate z values on the grid
  interpolated <-
    with(exam_curvature,
         akima::interp(x, y, z,
                       xo = x_range,
                       yo = y_range,
                       duplicate = 'strip'
         ))

  # identify the corresponding range of breakpoints for a contour plot.
  # note that the interpolated data might not include the highest peak, so
  # determine the break points for the contour lines on the interpolated data
  z_levels <- plot_scale(interpolated$z)

  # otherwise, get the range from the measured data, not the interpolated data
  #z_levels <- plot_scale(plot_dat$z)

  # get contour lines for the identified breakpoints
  # the contour lines will be at the boundaries between bins
  # their names start at the upper bound of the first bin
  contours_lst <- contourLines(interpolated, nlevels = length(z_levels), levels = z_levels)

  # identify the levels for each list item
  contour_levels <- sapply(contours_lst, function(item) item$level)

  # count the number of list items (i.e., the number of contour lines) that comprise each level
  segment_counts <- contour_levels |>
    table()

  return(segment_counts)
}


# This is a helper function for get_contour(): When add_edge = T, get_contour()
# adds the edge of the scanned cornea to the shape. This function returns the
# index of a boundary edge (min or max) of a vector of Boolean values, as
# defined by TRUE values.
find_edge <- function(x, min_max = 'min') {
  idx <- which(x == T)
  if (length(idx) == 0) {
    return(NA)
  } else {
    if (min_max == 'min')
      return(min(idx))
    else
      return(max(idx))
  }
}

# Establish the max ring_diam to use for finding the characteristic contour
contour_max_diameter <- 4.5  # Not currently used

#################
## get_contour
#################

#' Get contour
#' @description
#' Returns the \emph{x} and \emph{y} coordinates of the contour for an exam,
#' \code{exam_curvature}, at the specified power, \code{contour_power}.
#' @param exam_curvature
#' A data frame containing one row for each curvature \code{measurement} and the
#' same columns as [IQeyes::sample_curvature].
#' @param add_edge
#' A Boolean. If the contour extends beyond the edge of the scanned cornea,
#' \code{TRUE} will add the edge of the scanned area, connecting the ends of the
#' contour and fully enclosing it. The \code{segment_id} for the added edge will
#' be \code{0}. If the contour doesn't extend beyond the edge of the cornea, no
#' edge will be added.
#' @param contour_power
#' A number identifying the power of the contour to return. Valid values for an
#' exam can be determined by calling [IQeyes::get_contour_levels].
#' @return
#' A data frame containing one row for each point of the contour and the same
#' columns as [IQeyes::sample_contour].
#' @details
#' A contour is comprised of one or more segments, uniquely identified by
#' \code{segment_id}.
#'
#' This function doesn't currently work for posterior surfaces because
#' [IQeyes::plot_scale] presumes the color scale is on the "absolute
#' scale", which is designed for the dioptric power of anterior surfaces.
#' @examples
#' get_contour(sample_curvature, contour_power = 45.5) |>
#'   head()
#'
#' curvature_plot(sample_curvature, labels = F) +
#'   ggplot2::geom_point(data = get_contour(sample_curvature, contour_power = 45.5))
#'
#' @family Contours
#'
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom dplyr mutate
#' @importFrom dplyr bind_rows
#' @importFrom dplyr bind_cols
#' @importFrom akima interp
#'
#' @export
get_contour <- function(exam_curvature, add_edge = F, contour_power) {

  # preserve the exam details (i.e., join fields and surface)
  exam_record <- exam_curvature |>
    dplyr::select(tidyselect::all_of(join_fields), surface) |>
    unique()

  if (nrow(exam_record) > 1) warning('More than one exam record contained in exam_curvature.')

  # code in the following if-block is the same as get_contour_levels()
  # TODO: refactor to reduce replication
  if (T) {
    if (nrow(exam_curvature) == 0) warning('No data in exam_curvature.')

    # calculate power
    exam_curvature <- exam_curvature |>
      dplyr::mutate(z = anterior_power(measurement)) |>
      dplyr::filter(!is.na(z))

    # Create a grid that spans the extents of the measured x and y axes
    x_range <- with(exam_curvature, seq(min(x), max(x), length.out = length(unique(x))))
    y_range <- with(exam_curvature, seq(min(y), max(y), length.out = length(unique(y))))

    # Interpolate z values on the grid
    interpolated <-
      with(exam_curvature,
           akima::interp(x, y, z,
                         xo = x_range,
                         yo = y_range,
                         duplicate = 'strip'
           ))

    # identify the corresponding range of breakpoints for a contour plot.
    # note that the interpolated data might not include the highest peak, so
    # determine the break points for the contour lines on the interpolated data
    z_levels <- plot_scale(interpolated$z)

    # otherwise, get the range from the measured data, not the interpolated data
    #z_levels <- plot_scale(plot_dat$z)

    # get contour lines for the identified breakpoints
    # the contour lines will be at the boundaries between bins
    # their names start at the upper bound of the first bin
    contours_lst <- contourLines(interpolated, nlevels = length(z_levels), levels = z_levels)

    # identify the levels for each list item
    contour_levels <- sapply(contours_lst, function(item) item$level)

    # count the number of list items (i.e., the number of contour lines) that comprise each level
    segment_counts <- contour_levels |>
      table()
  }

  # if contour_power is specified, use that
  if (is.na(contour_power)) {
    warning('contour_power must be specified.')
  } else {
    if(is.na(segment_counts[as.character(contour_power)])) warning('No contour at specified power.')
    break_name <- as.character(contour_power)
  }

  # get item(s) from the list of contours that match the identified break_name
  target_items <- lapply(contours_lst, function(i) i$level == break_name) |>
    unlist() |>
    which()

  # reduce the list to just the target elements
  contours_lst <- contours_lst[target_items]

  # collate the points that comprise the shape's outline
  outline <-
    lapply(seq_along(contours_lst), function(i) {
      data.frame(x = contours_lst[[i]]$x,
                 y = contours_lst[[i]]$y,
                 # keep track of the individual contours
                 segment_id = i,
                 # remember the value of the characteristic contour
                 contour = contours_lst[[i]]$level)
    }) |>
    dplyr::bind_rows()

  if (add_edge) {
    # TODO: arrange the points in the same way as the other contours
    # get the perimeter of the map

    # TODO: make sure this works for multiple segments

    z_NA <- is.na(interpolated$z)

    edge_bottom <- apply(!z_NA, 1, function(x) find_edge(x, 'min'))
    edge_top <- apply(!z_NA, 1, function(x) find_edge(x, 'max'))
    edge_left <- apply(!z_NA, 2, function(x) find_edge(x, 'min'))
    edge_right <- apply(!z_NA, 2, function(x) find_edge(x, 'max'))

    # x and y are reversed (sort of): the first index of the grid is x, so a
    # "View" of the transpose puts x and y on the Cartesian coordinate system
    #View(t(interpolated$z))

    edge_coords <- data.frame(
      idx_x = c(edge_left, edge_right, seq_along(interpolated$x), seq_along(interpolated$x)),
      idx_y = c(seq_along(interpolated$y), seq_along(interpolated$y), edge_top, edge_bottom)
    ) |>
      dplyr::filter(!is.na(idx_x), !is.na(idx_y))  |>
      dplyr::mutate(x = interpolated$x[idx_x],
             y = interpolated$y[idx_y],
             z = interpolated$z[cbind(idx_x, idx_y)],
             contour = z_levels[findInterval(z, z_levels)])

    outline <- edge_coords |>
      dplyr::filter(contour == as.numeric(break_name)) |>
      dplyr::mutate(segment_id = 0) |>
      dplyr::select(x, y, segment_id, contour) |>
      dplyr::bind_rows(outline)

    #ggplot2::ggplot(outline, ggplot2::aes(x = x, y = y, color = as.factor(segment_id))) + ggplot2::geom_point()
  }

  # return exam_record with outline's x- and y-axis values
  exam_record |>
    dplyr::bind_cols(outline)
}

