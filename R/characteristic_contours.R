########################
## candidate_contours
########################

#' Identify candidates for the characteristic contour
#' @description
#' Identifies the contours of an exam that are candidates for the
#' characteristic contour.
#' @param exam_curvature
#' A data frame containing one row for each curvature \code{measurement} and the
#' same columns as [IQeyes::sample_curvature].
#' @return
#' An emph{n}-row data frame containing one row for each candidate contour. Each
#' row will contain the [IQeyes::join_fields], the cornea \code{surface}, and
#' the dioptric power of the candidate \code{contour}. In addition, each row
#' will contain the number of \code{segments} that comprise the contour, and the
#' radial distance of K-max(\code{r_max}) and K-next (\code{r_next}) from the
#' apex of the curvature map.
#' @details
#' The underlying algorithm operates at the segment level. If any contour has a
#' segment that contains both K-max and K-next, only contours with segments
#' containing both points are considered candidates. If no contour has a segment
#' that contains both points (or K-next doesn't exist), contours with segments
#' containing a single point are considered candidates. In the unlikely case
#' that no contours have any segments containing any points, all contours are
#' considered candidates.
#' @examples
#' candidate_contours(sample_curvature)
#'
#' @family Characteristic Contours
#'
#' @importFrom dplyr select
#' @importFrom dplyr mutate
#' @importFrom dplyr rename
#' @importFrom dplyr inner_join
#' @importFrom dplyr bind_rows
#' @importFrom dplyr bind_cols
#' @importFrom dplyr filter
#' @importFrom dplyr case_when
#' @importFrom tidyselect all_of
#' @importFrom sp point.in.polygon
#'
#' @export
candidate_contours <- function(exam_curvature) {

  #exam_curvature <- sample_curvature

  # preserve the exam details (i.e., join fields and surface)
  exam_record <- exam_curvature |>
    dplyr::select(tidyselect::all_of(join_fields), surface) |>
    unique()

  # get k_next and its x- and y-coordinates
  point_dat <- k_next(exam_curvature) |>
    dplyr::mutate(x = polar_to_cartesian(ring_diam, angle)$x,
           y = polar_to_cartesian(ring_diam, angle)$y) |>
    dplyr::mutate(point_name = 'k_next') |>
    dplyr::rename(power = k_next,
           radius = ring_diam) |>
    dplyr::select(point_name, x, y, radius, power)

  # add k_max and its x- and y-coordinates
  point_dat <- exam_curvature |>
    dplyr::filter(column_name == 'R_min') |>
    dplyr::mutate(power = anterior_power(measurement),
                  point_name = 'k_max') |>
    dplyr::rename(radius = ring_diam) |>
    dplyr::select(colnames(point_dat)) |>
    dplyr::bind_rows(point_dat)

  # is k_next within the plotted region of the curvature map?
  # is k_max within the plotted region of the curvature map?
  r_max = point_dat$radius[which(point_dat$point_name == 'k_max')]
  r_next = point_dat$radius[which(point_dat$point_name == 'k_next')]

  # remove the row for k_next if value is NA
  point_dat <- point_dat |>
    dplyr::filter(!is.na(power))

  # get complete list of contour levels for this exam
  segments <- get_contour_levels(exam_curvature)

  # create copies of the segments table to track the whether the segments that
  # comprise a contour enclose k_max and k_next
  partial_levels <- contour_levels <- segments

  # Contours are comprised of one or more segments. Evaluate each segment to see
  # if it contains k_max and/or k_next. We're interested in contours whose
  # segments contain both points, even if the points are in separate segments.
  # Evaluating the individual segments is a higher bar than evaluating the
  # segments of the contour all at once, ensuring that an individual segment
  # is complete enough to enclose at least one the key points.

  for (i in seq_along(contour_levels)) {
    #i <- 1
    polygon_dat <- get_contour(exam_curvature, add_edge = F, contour_power = names(contour_levels)[i])

    in_out_combined <- c(0, 0)

    # evaluate each individual segment to see if it contains k_max or k_next
    for (j in unique(polygon_dat$segment_id)) {
      in_out <- with(
        dplyr::filter(polygon_dat, segment_id == j),
        sp::point.in.polygon(point_dat$x, point_dat$y, x, y)
      )
      # OR the result for each segment with a combined result for the contour
      in_out_combined <- in_out_combined | in_out
    }

    # T if *both* k_max and k_next are enclosed
    contour_levels[i] <- in_out_combined[1] & in_out_combined[2]
    # T if *either* k_max or k_next is enclosed
    partial_levels[i] <- in_out_combined[1] | in_out_combined[2]
  }

  # convert vector of candidate level to Boolean, prioritizing AND or OR
  contour_levels <- dplyr::case_when (
    # at least one level contains both k_max and k_next, so only complete levels
    # are candidates
    sum(contour_levels) > 0 ~ as.logical(contour_levels),
    # at least one level contains one, so the partial levels are candidates
    sum(partial_levels) > 0 ~ as.logical(partial_levels),
    # no level contains either point, so every contour is a candidate
    .default = rep(T, length(contour_levels))
  )

  # reduce both segments and contour_levels to just candidates
  segments <- segments[contour_levels == T]
  #contour_levels <- contour_levels[contour_levels == T]

  # create a data frame of the candidates
  candidates <- data.frame(contour = as.numeric(names(segments)),
                           segments = as.numeric(segments),
                           r_max = r_max,
                           r_next = r_next
  )

  return(dplyr::bind_cols(exam_record, candidates))
}


#####################
## contour_context
#####################

#' Identify candidates for characteristic contour and calculate their distance
#' from K-max and K-next
#' @description
#' Identifies the contours of an exam that are candidates for the
#' characteristic contour, and calculates the number of \emph{other}
#' contours that exist between the candidate contour and K-max and K-next.
#' @param exam_curvature
#' A data frame containing one row for each curvature \code{measurement} and the
#' same columns as [IQeyes::sample_curvature].
#' @return
#' A data frame, containing the output of [IQeyes::candidate_contours] and three
#' additional columns: \code{k_max_distance}, the number of other contours that
#' exist between a candidate and K-max, \code{k_next_distance}, the number of
#' contours that exist between a candidate and K-next, and
#' \code{k_max_next_distance}, the difference between \code{k_max_distance} and
#' \code{k_next_distance}.
#' @details
#' When trying to identify an exam's characteristic contour, the actual values
#' of K-max and K-next don't matter, but the relative context that they create
#' \emph{is} important. This function quantifies that context by counting the
#' number of \emph{other} contours that exist between a candidate contour and
#' K-max and K-next.
#' @examples
#' contour_context(sample_curvature)
#'
#' @family Characteristic Contours
#'
#' @importFrom dplyr select
#' @importFrom tidyselect all_of
#'
#' @export
contour_context <- function(exam_curvature) {

  #exam_curvature <- sample_curvature

  # get the values of k_max and k_next
  k_values <- k_next(exam_curvature) |>
    dplyr::select(tidyselect::all_of(join_fields), k_max, k_next)

  # get the candidate contours
  candidates <- candidate_contours(exam_curvature)

  # find contour that contains k_max (based on the absolute scale)
  k_max_interval <- findInterval(k_values$k_max, absolute_scale)
  # determine number of contours between every contour and k_max contour
  k_max_distance <- k_max_interval - seq_along(absolute_scale)
  # extract distances for just the candidates
  k_max_distance <- k_max_distance[which(absolute_scale %in% candidates$contour)]
  # add to the data frame
  candidates$k_max_distance <- k_max_distance

  # repeat for k_next (if it exists)
  if (!is.na(k_values$k_next)) {
    k_next_interval <- findInterval(k_values$k_next, absolute_scale)
    k_next_distance <- k_next_interval - seq_along(absolute_scale)
    k_next_distance <- k_next_distance[which(absolute_scale %in% candidates$contour)]

    candidates$k_next_distance <- k_next_distance

    # number of contours between k_max and k_next
    candidates$k_max_next_distance <- k_max_interval - k_next_interval

  } else {
    candidates$k_next_distance <- NA
    candidates$k_max_next_distance <- NA

  }

  return(candidates)
}


########################
## closest_silhouette
########################

#' Identify the closest canonical silhouette
#' @description
#' Compares the silhouette of an exam's contour to the silhouettes of the
#' canonical shapes and identifies the canonical shape with the highest
#' percentage of overlap.
#'
#' This is essentially a wrapper function for
#' [IQeyes::silhouette_compare_group], so it is no less expensive from a
#' time/compute perspective. See details below.
#'
#' @param exam_curvature
#' A data frame containing one row for each curvature \code{measurement} and the
#' same columns as [IQeyes::sample_curvature].
#' @param exam_power
#' A number identifying the power of the exam's contour. Valid values for an
#' exam can be determined by calling [IQeyes::get_contour_levels].
#' @param exam_astig
#' A number, the angle (in degrees) of the exam's primary axis of astigmatism.
#' @param return_detail
#' A Boolean. \code{FALSE} (the default) to return just the index of the
#' canonical shape that has the highest percentage of overlap with the reference
#' silhouette \emph{and} the percentage of overlap. \code{TRUE} to return the
#' percentage overlap between the reference silhouette and every canonical
#' shape.
#'
#' @return
#' A data frame containing the [IQeyes::join_fields], corneal \code{surface},
#' and three additional columns: \code{contour}, the dioptric power of the
#' exam's contour, \code{closest_fit}, the index of the canonical silhouette
#' that has the highest percentage overlap with the exam's silhouette, and
#' \code{overlap}, the average percentage of overlap between the exam's
#' silhouette and the closest fitting member of the group.
#'
#' If \code{return_detail = T}, the data frame will include an additional
#' set of columns, one for each canonical shape, containing the percentage of
#' overlap between the exam and that shape.
#'
#' @details
#' This function creates all of the necessary objects to perform the comparison.
#' Calling [IQeyes:silhouette_compare] or [IQeyes:silhouette_group_compare]
#' directly requires the creation of two objects for every contour, a silhouette
#' (filled shape) and a polygon. This function creates those objects for the
#' exam, but it requires the power of the contour and the angle of astigmatism
#' to do so. The silhouettes [IQeyes:canonical_silhouettes] and polygons
#' [IQeyes:canonical_polygons] for the canonical shapes are included as objects
#' in this package.
#'
#' @examples
#' closest_silhouette(sample_curvature, exam_power = 45.5, exam_astig = 33.6, return_detail = T) |>
#'   dplyr::select(-all_of(join_fields))
#'
#' @family Characteristic Contours
#'
#' @importFrom dplyr select
#' @importFrom dplyr bind_cols
#' @importFrom tidyselect all_of
#'
#' @export
closest_silhouette <- function(exam_curvature, exam_power, exam_astig, return_detail = F) {

  # extract the exam details (i.e., join fields and surface)
  exam_record <- exam_curvature |>
    dplyr::select(tidyselect::all_of(join_fields), surface) |>
    unique()

  # get the filled shape for the exam at the specified power, scaled and rotated
  # at the specified astigmatism
  ref_shape <- fill_contour(exam_curvature,
                            scale_rotate = T,
                            axs = exam_astig,
                            contour_power = exam_power)

  # get the polygon for the exam at the specified power, scaled and rotated at
  # the specified astigmatism
  ref_polygon <- get_contour(exam_curvature, contour_power = exam_power) |>
    scale_rotate(axs = exam_astig) |>
    contour_polygon()

  # find the closest member of the group
  compare_results <- silhouette_compare_group(
    ref_shape = ref_shape,
    ref_polygon = ref_polygon,
    group_shapes = canonical_silhouettes,
    group_polygons = canonical_polygons,
    show_plots = F,
    return_detail = T
  )

  # create a single-row data frame for the results
  base_result <-
    dplyr::bind_cols(
      exam_record,
      data.frame(
        contour = exam_power,
        closest_fit = compare_results$closest_fit,
        overlap = compare_results$overlap[compare_results$closest_fit]
      )
    )

  if (return_detail) {
    # add the shape distances to each of the canonical contours
    names(compare_results$overlap) <- canonical_shapes$shape

    base_result |>
      dplyr::bind_cols(t(compare_results$overlap))

  } else {
    base_result
  }

}
