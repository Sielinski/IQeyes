##################
## scale_rotate
##################

#' Scale and rotate
#'
#' @description
#' Scales and rotates a shape defined by a set of \emph{x} and \emph{y}
#' coordinates. Can be either a curvature (e.g., [IQeyes::sample_curvature]) or
#' contour data frame (e.g., [IQeyes::sample_contour]).
#'
#' @param shape
#' A data frame containing one row for each point that defines a shape (e.g.,
#' a silhouette or contour). Currently, \code{scale_rotate()} expects the
#' data frame to contain the join fields.
#' @param axs
#' A number. Oftentimes, this should be the angle (in degrees) of the primary
#' axis of astigmatism for an exam.
#' @param r_target
#' A number defining the target radial length of the scaled shape.
#'
#' @return
#' A data frame containing the same columns as \code{shape} with transformed
#' \emph{x} and \emph{y} coordinates. No other values are changed.
#'
#' @details
#' The angle of rotation is determined by the primary axis of astigmatism for
#' the exam. If \code{axs} > 90, its value is shifted by -180°, ensuring that
#' the contour isn't flipped along the horizontal axis.
#'
#' The scaling factor is determined by the max radial length of any point on the
#' shape, not the end-to-end length of the shape, which is especially
#' important for asymmetric shapes.
#'
#' @examples
#' scale_rotate(sample_contour, axs = 33.6) |>
#'   head()
#'
#' scale_rotate(sample_curvature, axs = 33.6) |>
#'   head()
#'
#' @family contours
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


########################
## get_contour_levels
########################

#' Get contour levels
#'
#' @description
#' Returns the full range of contours and the number of segments needed to plot
#' \code{exam_curvature} on a curvature map.
#'
#' @param exam_curvature
#' A data frame of numeric values containing \code{x} and \code{y}, the
#' Cartesian coordinates of a point, and either \code{z} (the anterior radius of
#' curvature) or \code{power} (the dioptric power) at each (\code{x}, \code{y}).
#'
#' @return
#' A table. The table's labels identify the contour levels, and the values in
#' the table identify the number of segments that comprise the contour.
#'
#' @details
#' This function returns valid contours for an exam based on the breaks defined
#' by the color scale of [IQeyes::plot_scale].
#'
#' Because [IQeyes::plot_scale] is designed for the dioptric power of anterior
#' surfaces, this function works for anterior surfaces only.
#'
#' @examples
#' sample_curvature |>
#'   interpolate_measurements() |>
#'   get_contour_levels()
#'
#' @family contours
#'
#' @importFrom dplyr mutate
#' @importFrom dplyr filter
#' @importFrom akima interp
#'
#' @export
get_contour_levels <- function(exam_curvature) {

  #exam_curvature <- sample_curvature |> interpolate_measurements()
  #exam_curvature <- z_dat

  if (nrow(exam_curvature) == 0) warning('No data in exam_curvature.')

  if (!'power' %in% colnames(exam_curvature)) {
    # calculate power
    if ('measurement' %in% colnames(exam_curvature)) {
      exam_curvature <- exam_curvature |>
        dplyr::mutate(power = anterior_power(measurement))
    } else if ('z' %in% colnames(exam_curvature)) {
      exam_curvature <- exam_curvature |>
        dplyr::mutate(power = anterior_power(z))
    } else {
      warning('Power data not found in exam_curvature.')
      return(NA)
    }
  }

  # turn the data frame into a matrix
  curvature_m <- reshape2::acast(data = exam_curvature,
                                 formula = y ~ x,
                                 value.var = 'power')

  # identify the corresponding range of breakpoints for a contour plot.
  # note that the interpolated data might not include the highest peak, so
  # determine the break points for the contour lines on the interpolated data
  contour_levels <- plot_scale(curvature_m)

  # get contour lines for the identified breakpoints
  # the contour lines will be at the boundaries between bins
  # their names start at the upper bound of the first bin
  # for contourLines, x is rows and y is cols: they are *not* Cartesian
  # coordinates, so we have to swap our Cartesian x and y
  contours_lst <- contourLines(x = rownames(curvature_m) |> as.numeric() |> round(1),
                               y = colnames(curvature_m) |> as.numeric() |> round(1),
                               z = curvature_m,
                               levels = contour_levels)

  if (length(contours_lst) == 0) {
    warning('exam_curvature does not have enough points to identify contours.')
    return(NA)
  }

  # identify the levels for each list item
  item_levels <- sapply(contours_lst, function(item) item$level)

  # count the number of list items (i.e., the number of contour lines) that comprise each level
  segment_counts <- item_levels |>
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
#'
#' @description
#' Returns the \emph{x} and \emph{y} coordinates of the contour for an exam,
#' \code{exam_curvature}, at the specified power, \code{contour_power}.
#'
#' @param exam_curvature
#' A data frame with the same structure as [IQeyes::sample_curvature],
#' containing one row for each curvature (i.e., radius of curvature)
#' \code{measurement}.
#' @param interp
#' A Boolean. \code{TRUE} to interpolate measurements. See details.
#' @param add_edge
#' A Boolean. If the contour extends beyond the edge of the scanned cornea,
#' \code{TRUE} will add the edge of the scanned area, connecting the ends of the
#' contour and fully enclosing it. The \code{segment_id} for the added edge will
#' be \code{0}. If the contour doesn't extend beyond the edge of the cornea, no
#' edge will be added.
#' @param contour_power
#' Required. A number identifying the power of the contour to return. Valid
#' values for an exam can be determined by calling [IQeyes::get_contour_levels].
#'
#' @return
#' A data frame containing one row for each point of the contour and the same
#' columns as [IQeyes::sample_contour].
#'
#' @details
#' This function requires a fairly high degree of resolution. If more data
#' are needed, use the parameter \code{interp = T}.
#'
#' The resulting contour will be comprised of one or more segments, uniquely
#' identified by \code{segment_id}.
#'
#' This function doesn't currently work for posterior surfaces because
#' [IQeyes::plot_scale] presumes the color scale is on the "absolute
#' scale", which is designed for the dioptric power of anterior surfaces.
#'
#' @examples
#' get_contour(sample_curvature, interp = T, contour_power = 45.5) |>
#'   head()
#'
#' curvature_plot(sample_curvature, labels = F) +
#' ggplot2::geom_point(data = get_contour(sample_curvature,
#'                                        interp = T,
#'                                        contour_power = 45.5),
#'                     ggplot2::aes(x, y))
#'
#' @family contours
#'
#' @importFrom dplyr rename
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom dplyr mutate
#' @importFrom dplyr bind_rows
#' @importFrom dplyr bind_cols
#'
#' @export
get_contour <- function(exam_curvature, interp = F, add_edge = F, contour_power = NA) {

  # exam_curvature <- sample_curvature

  # preserve the exam details (i.e., join fields and surface)
  exam_record <- exam_curvature[1, ] |>
    dplyr::select(tidyselect::all_of(join_fields), surface)

  if (nrow(exam_record) > 1) {
    warning('More than one exam (or surface) contained in exam_curvature.')
    return(NA)
  }

  if (is.na(contour_power) || length(contour_power) == 0) {
    warning('contour_power must be specified.')
    return(exam_record)
  }

  # This function requires a fairly high degree of resolution. If more data
  # are needed, they can be interpolated
  if (interp) {
    exam_curvature <- exam_curvature |>
      interpolate_measurements() |>
      # the output of interpolate_measurements() is a z column: rename
      dplyr::rename(measurement = z) |>
      #dplyr::mutate(surface = 'FRONT') |>
      dplyr::cross_join(exam_record)
  }

  # the code in this if-block is the same as get_contour_levels
  if (T) {
    if (!'power' %in% colnames(exam_curvature)) {
      if ('measurement' %in% colnames(exam_curvature)) {
        # calculate power
        exam_curvature <- exam_curvature |>
          dplyr::mutate(power = anterior_power(measurement))
      } else if ('z' %in% colnames(exam_curvature)) {
        exam_curvature <- exam_curvature |>
          dplyr::mutate(power = anterior_power(z))
      } else {
        warning('Power data not found in exam_curvature.')
        return(invisible(exam_record))
      }
    }

    # turn the data frame into a matrix
    # formula expects the rows (y) to be on the left-hand side
    # this function inverts the y-axis values
    curvature_m <- reshape2::acast(data = exam_curvature,
                                   formula = y ~ x,
                                   value.var = 'power')

    # identify the corresponding range of breakpoints for a contour plot.
    # note that the interpolated data might not include the highest peak, so
    # determine the break points for the contour lines on the interpolated data
    contour_levels <- plot_scale(curvature_m)

    #if (F) {
    #  ggplot(exam_curvature, aes(x = x, y = y, z = power)) +
    #    geom_contour(breaks = contour_levels)
    #
    #  reshape2::melt(curvature_m, varnames = c('y', 'x'), value.name = "power", na.rm = T) |>
    #    ggplot(aes(x = x, y = y, z = power)) +
    #    geom_contour(breaks = contour_levels)
    #}

    # get contour lines for the identified breakpoints
    # the contour lines will be at the boundaries between bins
    # their names start at the upper bound of the first bin
    # for contourLines, x is rows and y is cols: they are *not* Cartesian
    # coordinates, so we have to swap our Cartesian x and y
    contours_lst <- contourLines(x = rownames(curvature_m) |> as.numeric() |> round(1),
                                 y = colnames(curvature_m) |> as.numeric() |> round(1),
                                 z = curvature_m,
                                 levels = contour_levels)

    if (length(contours_lst) == 0) {
      warning('exam_curvature does not have enough points to identify contours.')
      return(invisible(exam_record))
    }

    # identify the levels for each list item
    item_levels <- sapply(contours_lst, function(item) item$level)

    # count the number of list items (i.e., the number of contour lines) that comprise each level
    segment_counts <- item_levels |>
      table()
  }

  # ensure that specified contour_power is an option
  if(is.na(segment_counts[as.character(contour_power)])) {
    warning('No contour at the specified contour_power.')
    return(invisible(exam_record))
  }
  break_name <- as.character(contour_power)

  # get item(s) from the list of contours that match the identified break_name
  target_items <- lapply(contours_lst, function(i) i$level == break_name) |>
    unlist() |>
    which()

  # reduce the list to just the target elements
  contours_lst <- contours_lst[target_items]

  # collate the points that comprise the shape's outline
  # again, for contourLines, x is rows and y is cols: they are *not* Cartesian
  # coordinates, so we have to swap x and y
  outline <-
    lapply(seq_along(contours_lst), function(i) {
      data.frame(x = contours_lst[[i]]$y,
                 y = contours_lst[[i]]$x,
                 # keep track of the individual segments
                 segment_id = i,
                 # remember the value of the characteristic contour
                 contour = contours_lst[[i]]$level)
    }) |>
    dplyr::bind_rows()

  #ggplot2::ggplot(outline, ggplot2::aes(x = x, y = y, color = as.factor(segment_id))) + ggplot2::geom_point()

  if (add_edge) {
    # TODO: arrange the points in the same way as the other contours
    # get the perimeter of the map

    # TODO: make sure this works for multiple segments

    # establish the boundary of the cornea
    z_NA <- is.na(curvature_m)

    # find the matrix indices of the cornea's boundary from the left and right
    edge_right <- apply(!z_NA, 1, function(x) find_edge(x, 'max'))
    edge_left <- apply(!z_NA, 1, function(x) find_edge(x, 'min'))
    edge_bottom <- apply(!z_NA, 2, function(x) find_edge(x, 'max'))
    edge_top <- apply(!z_NA, 2, function(x) find_edge(x, 'min'))

    x_values <- colnames(curvature_m) |> as.numeric()
    y_values <- rownames(curvature_m) |> as.numeric()

    # turn those coordinates into a data frame
    edge_coords <- data.frame(
      idx_y = c(seq_along(edge_left), seq_along(edge_right), edge_bottom, edge_top),
      idx_x = c(edge_left, edge_right, seq_along(edge_bottom), seq_along(edge_top))
    ) |>
      dplyr::filter(!is.na(idx_x), !is.na(idx_y))  |>
      dplyr::mutate(x = x_values[idx_x],
             y = y_values[idx_y],
             power = curvature_m[cbind(idx_y, idx_x)],
             interval = findInterval(power, contour_levels)) |>
      unique()

    # filter out intervals that are below legitimate levels
    # TODO: Why are there illegitimate levels
    edge_coords <- edge_coords |>
      dplyr::filter(interval > 0) |>
      dplyr::mutate(contour = contour_levels[findInterval(power, contour_levels)])

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
