##############
## contours
##############

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
#' A number defining the target radius of the scaled shape.
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
#' @seealso
#' [IQeyes::fill_contour()], [IQeyes::contour_polygon()]
#' [IQeyes::silhouette_compare()], [IQeyes::silhouette_compare_group()]
#' @examples
#' scale_rotate(sample_contour, axs = 36.2) |> head()
#' scale_rotate(sample_curvature, axs = 36.2) |> head()
#'
scale_rotate <- function(shape, axs = 0, r_target = 4) {

  axs = ifelse(axs > 90, axs - 180, axs)

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
#' Fills a characteristic contour with points, creating a silhouette of its
#' shape.
#' @param exam_curvature
#' A data frame containing one row for each curvature \code{measurement} and the
#' same columns as [IQeyes::sample_curvature].
#' @param exam_contour
#' A data frame containing one row for each point of the contour and the same
#' columns as [IQeyes::sample_contour].
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
#' the points of the filled and transformed contour.
#' @details
#' Currently only handles the anterior (front) surface of the cornea.
#'
#' Because the function creates data (through interpolation), it cannot return
#' any additional columns that may have been included as part of the incoming
#' \code{exam_curvature}.
#' @seealso
#' [IQeyes::scale_rotate()], [IQeyes::contour_polygon()],
#' [IQeyes::silhouette_compare()], [IQeyes::silhouette_compare_group()]
#' @examples
#' fill_contour(sample_curvature, contour_power = 43, axs = 36.2) |> head()
#'
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
      dplyr::mutate(power = anterior_power(measurement),
                    ring_diam = cartesian_to_polar(x, y)$r,
                    angle = cartesian_to_polar(x, y)$theta,
                    column_name = 'interpolated'
                    )
  }

  # any point with a power ≥ the contour will be inside the contour
  curvature_dat <- curvature_dat |>
    dplyr::filter(power >= contour_power)

  # scale and rotate (*after* data frame has been reduced to the target contour)
  if (scale_rotate) curvature_dat <- scale_rotate(curvature_dat, ...)

  # put columns back in the same order as a typical curvature data frame
  curvature_dat |>
    dplyr::select(dplyr::all_of(colnames(sample_curvature)))

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
#' A contour is comprised of one or more segments. \code{contour_polygon}
#' takes the first point of each segment (based on its ordinal
#' position in the original data frame) and appends a replica to the end of the
#' the original data frame.
#' @seealso
#' [IQeyes::fill_contour()], [IQeyes::scale_rotate()]
#' [IQeyes::silhouette_compare()], [IQeyes::silhouette_compare_group()]
#' @examples
#' contour_polygon(sample_contour) |> head()
#'
contour_polygon <- function(exam_contour) {

  # assign a number to every row of every segment
  contour_rn <- exam_contour |>
    dplyr::group_by(dplyr::across(tidyselect::all_of(join_fields)), contour_id) |>
    dplyr::mutate(rn = dplyr::row_number()) |>
    dplyr::select(tidyselect::all_of(colnames(exam_contour)), rn)

  # find the first and last row of every segment
  contour_limits <- contour_rn |>
    dplyr::group_by(dplyr::across(tidyselect::all_of(join_fields)), contour_id) |>
    dplyr::summarize(min_rn = min(rn), max_rn = max(rn)) |>
    dplyr::ungroup()

  # connect the ends of every segment by adding a duplicate of the first
  # point after the last point
  contour_polygon <- contour_rn |>
    dplyr::filter(rn == 1) |>
    dplyr::inner_join(contour_limits, by = c(join_fields, 'contour_id')) |>
    dplyr::mutate(rn = max_rn + 1) |>
    dplyr::select(colnames(contour_rn)) |>
    dplyr::bind_rows(contour_rn) |>
    dplyr::arrange(dplyr::across(tidyselect::all_of(join_fields)), contour_id, rn)

  contour_polygon |>
    dplyr::select(tidyselect::all_of(colnames(exam_contour)))

}


########################
## silhouette_compare
########################

#' Compare the silhouettes of two contours
#' @description
#' Each of silhouettes is defined by two objects, a filled shape and a polygon.
#' These objects must be created before calling \code{silhouette_compare()}. The
#' filled shapes are created by calling \code{fill_contour()}, and polygon
#' objects are created by calling \code{contour_polygon()}.
#' @param shape_A
#' A data frame containing the filled shape for contour A.
#' @param polygon_A
#' A data frame containing the polygon object for contour A.
#' @param shape_B
#' A data frame containing the filled shape for contour B.
#' @param polygon_B
#' A data frame containing the polygon object for contour B.
#' @param show_plots
#' A Boolean. \code{TRUE} to render two plots that illustrate comparisons being
#' made.
#' @return
#' The percentage of overlap between the two silhouettes, averaged over both
#' directions.
#' @seealso
#' [IQeyes::fill_contour()], [IQeyes::contour_polygon()],
#' [IQeyes::scale_rotate()], [IQeyes::silhouette_compare_group()]
#' @examples
#' silhouette_compare(
#'   shape_A = fill_contour(sample_curvature, sample_contour),
#'   polygon_A = contour_polygon(sample_contour),
#'   shape_B = fill_contour(sample_curvature, sample_contour),
#'   polygon_B = contour_polygon(sample_contour)
#' )
#'
silhouette_compare <- function(shape_A, polygon_A,
                               shape_B, polygon_B,
                               show_plots = F) {

  # identify the points of the shape A that lie within the outline of shape B
  A_to_B <- sp::point.in.polygon(
    point.x = shape_A$x,
    point.y = shape_A$y,
    pol.x = polygon_B$x,
    pol.y = polygon_B$y
  )

  # identify the points of the shape B that lie within the outline of shape A
  B_to_A <- sp::point.in.polygon(
    point.x = shape_B$x,
    point.y = shape_B$y,
    pol.x = polygon_A$x,
    pol.y = polygon_A$y
  )

  A_to_B_pct <- sum((as.logical(A_to_B))) / length(A_to_B)
  B_to_A_pct <- sum((as.logical(B_to_A))) / length(B_to_A)

  if (show_plots) {
    # plot A to B
    p <- shape_A |>
      dplyr::bind_cols(inside = as.logical(A_to_B)) |>
      ggplot2::ggplot(ggplot2::aes(x, y)) +
      ggplot2::geom_point(ggplot2::aes(color = inside)) +
      ggplot2::geom_point(data = polygon_B, color = 'grey') +
      ggplot2::labs(title = 'A\'s shape with a shadow cast by B\'s contour')
    print(p)

    # plot B to A
    p <- shape_B |>
      dplyr::bind_cols(inside = as.logical(B_to_A)) |>
      ggplot2::ggplot(ggplot2::aes(x, y)) +
      ggplot2::geom_point(ggplot2::aes(color = inside)) +
      ggplot2::geom_point(data = polygon_A, color = 'grey') +
      ggplot2::labs(title = 'B\'s shape with a shadow cast by A\'s contour')
    print(p)
  }

  (A_to_B_pct + B_to_A_pct) / 2
}


##############################
## silhouette_compare_group
##############################

#' Compare the silhouette of one contour with a group of silhouettes
#' @description
#' Every silhouette is defined by two objects, a filled shape and a polygon.
#' These objects must be created before calling
#' \code{silhouette_compare_group()}. The filled shapes are created by calling
#' \code{fill_contour()}, and the polygon objects are created by calling
#' \code{contour_polygon()}.
#'
#' Each of the comparisons is made in both directions.
#' @param ref_shape
#' A data frame containing the filled shape for the reference contour
#' @param ref_polygon
#' A data frame containing the polygon object for the reference contour
#' @param group_shapes
#' A data frame containing filled shapes for the group of contours.
#' \emph{Requires} a \code{cluster} column that uniquely identifies
#' each member of the group. The cluster IDs for individual contours (i.e.,
#' members of the group) must be the same in both \code{group_shapes} and
#' \code{group_polygons}.
#' @param group_polygons
#' A data frame containing polygon objects for the group of contours.
#' \emph{Requires} a \code{cluster} column that uniquely identifies
#' each member of the group. The cluster IDs for individual contours (i.e.,
#' members of the group) must be the same in both \code{group_shapes} and
#' \code{group_polygons}.
#' @param show_plots
#' A Boolean. \code{TRUE} to render two plots that illustrate comparisons being
#' made.
#' @param return_detail
#' A Boolean. \code{TRUE} to return a list object with two elements:
#' \code{closest_fit}, the index of the group member that has the highest
#' percentage of overlap with the reference silhouette, and
#' \code{average_distance}, a vector containing the average percentage of
#' overlap between the reference silhouette and every member of the group.
#' @details
#' Both \code{group_shapes} and \code{group_polygons} require a \code{cluster}
#' column. Although the cluster column can be passed into and returned by
#' \code{contour_polygon()}, the same cannot be done by \code{fill_contour()},
#' which means the process for creating \code{group_shapes} requires the
#' additional step of appending the \code{cluster} column.
#' @return
#' The index of the group member that has the highest percentage of overlap with
#' the reference silhouette.
#' @seealso
#' [IQeyes::fill_contour()], [IQeyes::contour_polygon()],
#' [IQeyes::scale_rotate()], [IQeyes::silhouette_compare()]
#'
# compare a shape's silhouette to a dataframe of polygons in both directions
silhouette_compare_group <- function(ref_shape, ref_polygon,
                                     group_shapes, group_polygons,
                                     show_plots = F, return_detail = F, ...) {

  # identify the points of the shape that lie within the boundary of every member of the group
  in_out <- group_polygons |>
    dplyr::group_by(cluster) |>
    dplyr::group_split() |>
    lapply(function(x) {
      df <- data.frame(
        inside = sp::point.in.polygon(
          point.x = ref_shape$x,
          point.y = ref_shape$y,
          pol.x = x$x,
          pol.y = x$y
        ))
      df$inside <- as.logical(df$inside)
      df$cluster <- min(x$cluster)
      df
    })

  # identify the points of every member of the group that lie within the boundary of the shape
  out_in <- group_shapes |>
    dplyr::group_by(cluster) |>
    dplyr::group_split() |>
    lapply(function(x) {
      df <- data.frame(
        inside = sp::point.in.polygon(
          point.x = x$x,
          point.y = x$y,
          pol.x = ref_polygon$x,
          pol.y = ref_polygon$y
        ))
      df$inside <- as.logical(df$inside)
      df$cluster <- min(x$cluster)
      df
    })

  in_out_pct <- lapply(in_out, function(x) sum(x$inside) / length(x$inside)) |> unlist()
  out_in_pct <- lapply(out_in, function(x) sum(x$inside) / length(x$inside)) |> unlist()

  if (show_plots) {
    # plot out_in
    out_in_dat <- lapply(out_in, function(x) {
      item_cluster <- min(x$cluster)
      # add a column to each cluster's fill, identifying whether the fill's
      # point is inside the shadow of the reference polygon
      group_shapes |>
        dplyr::filter(cluster == item_cluster) |>
        dplyr::bind_cols(inside = x$inside)
    }) |>
      dplyr::bind_rows()

    p <- ggplot2::ggplot(out_in_dat, ggplot2::aes(x, y)) +
      ggplot2::geom_point(ggplot2::aes(color = inside)) +
      ggplot2::geom_point(data = ref_polygon, color = 'grey') +
      ggplot2::facet_wrap(. ~ cluster) +
      ggplot2::labs(title = 'The group\'s shapes with a shadow cast by the reference contour')

    print(p)

    # plot in_out
    in_out_dat <- lapply(in_out, function(x) {
      # add two columns to the reference fill, identifying the cluster of
      # comparison and whether the fill's point is inside the shadow of that
      # cluster's polygon
      ref_shape |>
        dplyr::bind_cols(inside = x$inside, cluster = x$cluster)
    }) |>
      dplyr::bind_rows()

    p <- ggplot2::ggplot(in_out_dat, ggplot2::aes(x, y)) +
      ggplot2::geom_point(ggplot2::aes(color = inside)) +
      ggplot2::geom_point(data = group_polygons, color = 'grey') +
      ggplot2::facet_wrap(. ~ cluster) +
      ggplot2::labs(title = 'The reference shape with shadows cast by the group\'s contours')

    print(p)

  }

  average_distance <- (in_out_pct + out_in_pct) / 2

  # identify the item with the greatest percentage coverage in both directions
  closest_fit <- which.max(average_distance)

  if (return_detail) {
    return(list(
      closest_fit = closest_fit,
      average_distance = average_distance
    ))

  } else {
    return(closest_fit)
  }
}

