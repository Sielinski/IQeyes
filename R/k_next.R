# find the adjacent points (on the same and adjacent diameters)
# pt is a row from a data frame that contains a ring_diam, angle, and measurement
adjacent_points <- function(source_dat, pt) {

  # find next nearest diameters (note that the vector will include R_min diameter)
  diameters <- source_dat$ring_diam |>
    unique() |>
    sort()
  diam_pt_idx <- which(diameters == pt$ring_diam)
  diam_lower_idx <- ifelse(diam_pt_idx == 1, 1, diam_pt_idx - 1)
  diam_upper_idx <- ifelse(diam_pt_idx == length(diameters), length(diameters), diam_pt_idx + 1)

  # append rows for each diameter
  adjacent_dat <- adjacent_points_diameter(source_dat, pt, diameters[diam_lower_idx])

  adjacent_dat <- adjacent_points_diameter(source_dat, pt, diameters[diam_upper_idx]) |>
    dplyr::bind_rows(adjacent_dat)

  adjacent_dat <- adjacent_points_diameter(source_dat, pt, diameters[diam_pt_idx]) |>
    dplyr::bind_rows(adjacent_dat)

  return(adjacent_dat)

}


# identify points that qualify as adjacent on the target diameter
# i.e., points with next higher or lower (or equal) angle relative
# to the source point
adjacent_points_diameter <- function(source_dat, pt, target_diam) {

  # available angles on the target diameter
  angles <- source_dat |>
    dplyr::filter(ring_diam == target_diam) |>
    dplyr::filter(!is.na(measurement)) |>
    dplyr::select(angle) |>
    unlist() |>
    sort()

  # if there's only one angle, either apex or k_max
  # if there are no angles, there are no available measurements
  if (length(angles) <= 1) {
    if (target_diam == 0) {
      # apex
      adjacent_angles <- 0
    } else {
      # k_max
      adjacent_angles <- NA
    }
  } else {
    if (pt$ring_diam == 0) {
      # if the starting point is the apex, need to consider all angles
      adjacent_angles <- angles
    } else {
      # pick 0° (if exists) and the next two smallest angles, which will often be the same
      angles_delta <- smallest_angle(pt$angle, angles)

      if (min(angles_delta) == 0) {
        angles_idx <- order(angles_delta)[1:3]
      } else {
        angles_idx <- order(angles_delta)[1:2]
      }

      adjacent_angles <- angles[angles_idx]
    }
  }

  adjacent_pts <- source_dat |>
    dplyr::filter(ring_diam == target_diam) |>
    # exclude the point itself
    dplyr::filter(!(ring_diam == pt$ring_diam & angle == pt$angle)) |>
    dplyr::filter(angle %in% adjacent_angles)

  return(adjacent_pts)

}


############
## k_next
############

#' Finds k-next
#' @description
#' Find the secondary peak on a curvature map.
#' @param exam_curvature
#' A data frame containing one row for each curvature \code{measurement} and the
#' same columns as \code{sample_curvature}.
#' @param just_points
#' A Boolean. \code{TRUE} to return the points that comprise the primary peak.
#' @return
#' If \code{just_points} = \code{F}, \code{k_next()} returns a one-row data
#' frame containing the \code{join_fields} and a variety of metadata related to
#' k-next.
#'
#' If \code{just_points} = \code{T}, \code{k_next()} returns a copy of
#' \code{exam_curvature} with an additional column \code{peak_id} identifying
#' the peak that each measured point belongs to.
#' @examples
#' k_next(sample_curvature)
#'
k_next <- function(exam_curvature, just_points = F) {

  # add an identifier flag to the cleansed curvature data
  source_dat <- exam_curvature |>
    #dplyr::inner_join(status_ok, by = join_fields) |>
    dplyr::filter(surface == 'FRONT') |>
    dplyr::filter(!is.na(measurement)) |>
    dplyr::mutate(peak_id = 0)

  # identify k_max (i.e., R_min)
  k_max_dat <- source_dat |>
    dplyr::filter(column_name == 'R_min')

  # if r_min is at the same location as another point, remove the latter
  match_rmin <- source_dat |>
    dplyr::semi_join(k_max_dat, by = c('ring_diam', 'angle'))

  if (nrow(match_rmin) > 1) {
    source_dat <- source_dat |>
      dplyr::anti_join(k_max_dat, by = c('ring_diam', 'angle')) |>
      dplyr::bind_rows(k_max_dat)
  }

  # add curvature measurements to the adjacent points data frame
  peak_tracking <- adjacent_points_dat |>
    dplyr::inner_join(source_dat, by = c('ring_diam_src' = 'ring_diam', 'angle_src' = 'angle')) |> #View()
    dplyr::rename(measurement_src = measurement) |>
    dplyr::select(ring_diam_src, angle_src, measurement_src, ring_diam, angle) |>
    dplyr::inner_join(source_dat, by = c('ring_diam', 'angle')) |>
    dplyr::select(ring_diam_src, angle_src, measurement_src, ring_diam, angle, measurement)

  if (nrow(match_rmin) == 1) {
    # add to k_max to tracking data frame: need both source and adjacent points
    k_max_source <- k_max_dat |>
      dplyr::select(ring_diam, angle, measurement) |>
      dplyr::rename(ring_diam_src = ring_diam,
             angle_src = angle,
             measurement_src = measurement)

    # points that are adjacent to k_max will be flagged as members of that peak
    k_max_adjascent <- adjacent_points(source_dat, k_max_dat) |>
      dplyr::select(ring_diam, angle, measurement)

    # combine k_max with all other data points
    peak_tracking <- dplyr::bind_cols(k_max_source, k_max_adjascent) |>
      dplyr::bind_rows(peak_tracking)
  }

  # create a flag to track the adjacent points that have been processed
  peak_tracking$processed <- F

  # keep track of the current peak
  current_peak <- 1

  # can we incorporate into the loop?
  pts_to_evaluate <- k_max_dat |>
    dplyr::select(ring_diam, angle)

  source_dat$peak_id[which(source_dat$column_name == 'R_min')] <- current_peak

  repeat {

    repeat {

      # identify unprocessed adjacencies to evaluate
      evaluate_dat <- peak_tracking |>
        dplyr::inner_join(pts_to_evaluate, by = c('ring_diam_src' = 'ring_diam', 'angle_src' = 'angle')) |>
        dplyr::filter(!processed)

      if (nrow(evaluate_dat) == 0) break

      # flag the points as evaluated
      evaluate_dat$processed <- T

      peak_tracking <-
        dplyr::rows_update(
          peak_tracking,
          evaluate_dat,
          by = c('ring_diam_src', 'angle_src', 'ring_diam', 'angle')
        )

      # see if curvature is greater (i.e., power is lower)
      evaluate_dat <- evaluate_dat |>
        dplyr::filter(measurement >= measurement_src)

      # identify points/rows in the source data frame that are members of the
      # current peak and record peak_id
      update_source <- source_dat |>
        dplyr::semi_join(evaluate_dat, by = c('ring_diam', 'angle', 'measurement')) |>
        dplyr::filter(peak_id == 0) |>
        dplyr::mutate(peak_id = current_peak)

      # update source data frame
      source_dat <-
        dplyr::rows_update(source_dat,
                    update_source,
                    by = c('ring_diam', 'angle', 'measurement'))

      pts_to_evaluate <- update_source |>
        dplyr::select(ring_diam, angle)
    }

    # find the next r_min
    pts_to_evaluate <- source_dat |>
      dplyr::filter(peak_id == 0) |>
      dplyr::slice_min(measurement, n = 1)

    if (nrow(pts_to_evaluate) == 0) break

    # increment peak number
    current_peak <- current_peak + 1

    # flag next r_min as member of new peak
    pts_to_evaluate$peak_id <- current_peak

    # update source data frame
    source_dat <-
      dplyr::rows_update(source_dat,
                  pts_to_evaluate,
                  by = c('ring_diam', 'angle', 'measurement'))

    pts_to_evaluate <- pts_to_evaluate |>
      dplyr::select(ring_diam, angle)

    # reset peak_tracking
    # no need to process any point that has already belongs to a peak
    peak_tracking <- peak_tracking |>
      dplyr::inner_join(source_dat, by = c('ring_diam' , 'angle', 'measurement')) |>
      dplyr::mutate(processed = ifelse(peak_id == 0, F, T)) |>
      dplyr::select(all_of(colnames(peak_tracking)))
  }

  #if (F) {
  #  # plot the points that were identified
  #  curvature_plot(plot_from_join(join_dat), labels = F) +
  #    ggplot2::geom_point(data = source_dat, aes(x = x, y = y, color = as.factor(peak_id), size = anterior_power(measurement)))
  #}

  if (just_points)
    return(source_dat)

  # check to see if there's another peak
  if(max(source_dat$peak_id) == 1) {
    # no other points on the cornea, so there is no k_next
    k_next <- data.frame(
      k_max = anterior_power(k_max_dat$measurement),
      k_next = NA,
      ring_diam = NA,
      angle = NA,
      k_max_coverage = 100,
      max_next_alignment_deg = NA,
      max_next_distance_mm = NA,
      max_next_delta = NA,
      peak_count = current_peak
    )
  } else {
    # identify the secondary peak, ordered by the number of points in the peak
    # include min_R in summary
    next_peak <- source_dat |>
      dplyr::filter(peak_id > 1) |>
      dplyr::summarize(cnt = dplyr::n(), min_r = min(measurement), .by = peak_id) |>
      # convert radius (R) to power (K) so we can take max
      dplyr::mutate(max_k = anterior_power(min_r)) |>
      dplyr::slice_max(data.frame(cnt, max_k), n = 1)

    # identify how well aligned and how far apart k_max and k_next are
    k_next_dat <- source_dat |>
      dplyr::semi_join(next_peak, by = c('peak_id', 'measurement' = 'min_r')) |>
      dplyr::sample_n(1)

    k_max_points <- source_dat |>
      dplyr::filter(peak_id == 1) |>
      nrow()

    k_next <- data.frame(
      k_max = anterior_power(k_max_dat$measurement),
      k_next = anterior_power(k_next_dat$measurement),
      ring_diam = k_next_dat$ring_diam,
      angle = k_next_dat$angle,
      k_max_coverage = round(k_max_points / nrow(source_dat) * 100, 2),
      max_next_alignment_deg = 180 - smallest_angle(k_max_dat$angle, k_next_dat$angle),
      #max_next_alignment_deg = smallest_angle(k_max_dat$angle, k_next_dat$angle),
      max_next_distance_mm = sqrt((k_max_dat$x - k_next_dat$x) ^ 2 + (k_max_dat$y - k_next_dat$y) ^ 2),
      peak_count = current_peak
    )

    k_next <- k_next |>
      dplyr::mutate(max_next_delta = k_max - k_next)

    # anytime there's more than one point with the same R_min, need to
    # disambiguate the point chosen as k_next
    k_next <- k_next |>
      dplyr::slice_min(order_by = data.frame(max_next_distance_mm, max_next_alignment_deg))

  }

  # recombine and return the join fields (to identify the exam) and k_next info
  exam_curvature |>
    dplyr::select(tidyselect::all_of(join_fields)) |>
    dplyr::slice(1) |>
    dplyr::bind_cols(k_next)
}
