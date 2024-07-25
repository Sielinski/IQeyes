########################
## candidate_contours
########################

#' Identify candidates for characteristic contour
#' @description
#' Identifies the contours of an exam that are candidates for the
#' characteristic contour.
#'
#' @param exam_curvature
#' A data frame with the same structure as [IQeyes::sample_curvature],
#' containing one row for each curvature (i.e., radius of curvature)
#' \code{measurement}.
#' @param exam_peaks
#' A data frame with the same structure as the output of [IQeyes::k_next],
#' containing one row for each peak.
#'
#' @return
#' An \emph{n}-row data frame containing one row for each candidate contour. Each
#' row will contain the [IQeyes::join_fields] and the following fields:
#'
#' \describe{
#'   \item{surface}{The corneal surface. Either \code{front} or \code{back}.}
#'   \item{contour}{The dioptric power of the candidate contour.}
#'   \item{segments}{The number of segments (polygons) that comprise the contour.}
#'   \item{r_max}{The radial distance of K-max from the apex.}
#'   \item{r_next}{The radial distance of K-next from the apex.}
#' }
#'
#' @details
#' This function is generally called by [IQeyes::contour_context], which
#' pre-processes the incoming \code{exam_curvature} and calls [IQeyes::k_next]
#' using that pre-processed data. If calling this function directly, ensure that
#' \code{exam_curvature} and \code{exam_peaks} are based on the same data.
#'
#' The underlying algorithm operates at the segment level. If any contour has a
#' segment that contains both K-max and K-next, only contours with segments
#' containing both points are considered candidates. If no contour has a segment
#' that contains both points (or K-next doesn't exist), contours with segments
#' containing a single point are considered candidates. In the unlikely case
#' that no contours have segments containing points, all contours are
#' considered candidates.
#'
#' @examples
#' exam_interp <- interpolate_measurements(sample_curvature)
#'
#' exam_peaks <- exam_interp |>
#'   dplyr::mutate(power = anterior_power(z)) |>
#'   dplyr::select(x, y, power) |>
#'   reshape2::acast(y ~ x, value.var = 'power') |>
#'   k_next()
#'
#' exam_curvature <- sample_curvature[1, ] |>
#'   dplyr::select(tidyselect::all_of(join_fields), 'surface') |>
#'   dplyr::cross_join(exam_interp)
#'
#' candidate_contours(exam_curvature, exam_peaks)
#'
#' @family characteristic contours
#'
#' @importFrom dplyr select
#' @importFrom dplyr mutate
#' @importFrom dplyr rename
#' @importFrom dplyr inner_join
#' @importFrom dplyr bind_rows
#' @importFrom dplyr bind_cols
#' @importFrom dplyr filter
#' @importFrom dplyr if_else
#' @importFrom dplyr case_when
#' @importFrom tidyselect all_of
#' @importFrom sp point.in.polygon
#'
#' @export
candidate_contours <- function(exam_curvature, exam_peaks) {

  #exam_curvature <- sample_curvature

  # preserve the exam details (i.e., join fields and surface)
  exam_record <- exam_curvature[1, ] |>
    dplyr::select(tidyselect::all_of(join_fields), surface)

  # get just k_max and k_next and their coordinates
  point_dat <- exam_peaks |>
    dplyr::filter(peak_id <= 2) |>
    dplyr::mutate(point_name = dplyr::if_else(peak_id == 1, 'k_max', 'k_next')) |>
    dplyr::select(point_name, x, y, r, power) |>
    dplyr::rename(radius = r)

  # get k_next and its x- and y-coordinates
  #point_dat <- k_next_old(exam_curvature) |>
  #  dplyr::mutate(x = polar_to_cartesian(ring_diam, angle)$x,
  #         y = polar_to_cartesian(ring_diam, angle)$y) |>
  #  dplyr::mutate(point_name = 'k_next') |>
  #  dplyr::rename(power = k_next,
  #         radius = ring_diam) |>
  #  dplyr::select(point_name, x, y, radius, power)
  #
  ## add k_max and its x- and y-coordinates
  #point_dat <- exam_curvature |>
  #  dplyr::filter(column_name == 'R_min') |>
  #  dplyr::mutate(power = anterior_power(measurement),
  #                point_name = 'k_max') |>
  #  dplyr::rename(radius = ring_diam) |>
  #  dplyr::select(colnames(point_dat)) |>
  #  dplyr::bind_rows(point_dat)

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
                           r_next = ifelse(length(r_next) == 0, NA, r_next)
  )

  return(dplyr::bind_cols(exam_record, candidates))
}


#####################
## contour_context
#####################

#' Identify candidates for characteristic contour
#'
#' @description
#' Identifies the contours of an exam that are candidates for the
#' characteristic contour, and calculates the number of \emph{other}
#' contours that exist between a candidate contour and both K-max and K-next.
#'
#' @param exam_curvature
#' A data frame with the same structure as [IQeyes::sample_curvature],
#' containing one row for each curvature (i.e., radius of curvature)
#' \code{measurement}.
#'
#' @return
#' A data frame, containing the output of [IQeyes::candidate_contours] and three
#' additional columns:
#'
#' \describe{
#'  \item{k_max_distance}{The number of other contours between a candidate and
#'  K-max.}
#'  \item{k_next_distance}{The number of other contours between a candidate and
#'  K-next.}
#'  \item{k_max_next_distance}{The difference between \code{k_max_distance} and
#' \code{k_next_distance}.}
#' }
#'
#' @details
#' When trying to identify an exam's characteristic contour, the actual values
#' of K-max and K-next don't matter, but the relative context that they create
#' is important. This function quantifies that context by counting the number of
#' \emph{other} contours that exist between a candidate contour and both K-max
#' and K-next.
#'
#' Internally, this function calls [IQeyes::k_next] and
#' [IQeyes::candidate_contours].
#'
#' @examples
#' contour_context(sample_curvature, interp = T) |>
#'   dplyr::select(-tidyselect::all_of(join_fields))
#'
#' @family characteristic contours
#'
#' @importFrom dplyr select
#' @importFrom tidyselect all_of
#'
#' @export
contour_context <- function(exam_curvature, interp = F) {

  #exam_curvature <- sample_curvature

  # preserve the exam details (i.e., join fields and surface)
  exam_record <- exam_curvature[1, ] |>
    dplyr::select(tidyselect::all_of(join_fields), surface)

  if (interp) {
    z_dat <- exam_curvature |>
      interpolate_measurements()
  } else {
    z_dat <- exam_curvature |>
      dplyr::select(x, y, measurement) |>
      dplyr::rename(z = measurement)
  }

  # get exam's curvature peaks
  exam_peaks <- z_dat |>
    dplyr::mutate(power = anterior_power(z)) |>
    reshape2::acast(y ~ x, value.var = 'power') |>
    k_next()

  k_values <- exam_peaks |>
    dplyr::filter(peak_id <= 2) |>
    dplyr::mutate(point_name = dplyr::if_else(peak_id == 1, 'k_max', 'k_next')) |>
    dplyr::select(point_name, power) |>
    tidyr::pivot_wider(names_from = point_name, values_from = power) |>
    dplyr::cross_join(exam_record)

  # get the values of k_max and k_next
  #k_values <- k_next_old(exam_curvature) |>
  #  dplyr::select(tidyselect::all_of(join_fields), k_max, k_next)

  # get the candidate contours
  candidates <- z_dat |>
    dplyr::mutate(power = anterior_power(z)) |>
    dplyr::rename(measurement = z) |>
    dplyr::cross_join(exam_record) |> #View()
    candidate_contours(exam_peaks)

  # find contour that contains k_max (based on the absolute scale)
  k_max_interval <- findInterval(k_values$k_max, absolute_scale)
  # determine number of contours between every contour and k_max contour
  k_max_distance <- k_max_interval - seq_along(absolute_scale)
  # extract distances for just the candidates
  k_max_distance <- k_max_distance[which(absolute_scale %in% candidates$contour)]
  # add to the data frame
  candidates$k_max_distance <- k_max_distance

  # repeat for k_next (if it exists)
  if ('k_next' %in% colnames(k_values) && !is.na(k_values$k_next)) {
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

