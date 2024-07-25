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
    pull(angle) |>
    #unlist() |>
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


################
## k_next_old
################

#' Find k-next
#' @description
#' Finds the secondary peak on a curvature map.
#'
#' @param exam_curvature
#' A data frame with the same structure as [IQeyes::sample_curvature],
#' containing one row for each curvature (i.e., radius of curvature)
#' \code{measurement}.
#' @param just_points
#' A Boolean. \code{TRUE} to return the points that comprise the primary peak.
#'
#' @return
#' If \code{just_points} = \code{T}, \code{k_next()} returns a copy of
#' \code{exam_curvature} with an additional column \code{peak_id} identifying
#' the peak that each measured point belongs to.
#'
#' If \code{just_points} = \code{F}, \code{k_next()} returns a one-row data
#' frame containing the \code{join_fields} and a variety of metadata related to
#' k-next:
#'
#' \describe{
#'  \item{k_max}{The dioptric power of K-max}
#'  \item{k_next}{The dioptric power of K-next}
#'  \item{ring_diam}{The radial distance of K-next from the vertex.}
#'  \item{angle}{The radial angle of K-next (in degrees).}
#'  \item{k_max_coverage}{The percentage of points on the cornea that descend
#'  from K-max.}
#'  \item{max_next_alignment_deg}{The angle (in degrees) between K-next and the
#'  projection of K-max in its opposite direction. See details below.}
#'  \item{max_next_distance_mm}{The Euclidean distance (in mm) between K-max and
#'  K-next.}
#'  \item{peak_count}{Number of peaks found on the cornea.}
#'  \item{max_next_delta}{The difference in power (in diopters) between K-max
#'  and K-next.}
#' }
#'
#' @details
#' The algorithm to identify K-next is based on measured (i.e., not interpreted)
#' points only, so K-next itself is one of the measured points.
#'
#' \code{max_next_alignment_deg} quantifies the departure from symmetry. When
#' K-max and K-next are symmetrically aligned, their radial axes are 180° apart.
#' \code{max_next_alignment_deg} is calculated as 180° \emph{minus} the smallest
#' angle, so when K-max and K-next are 180° apart, the result is 0°. Generally,
#' the smallest angle between the K-max and K-next will be less than 180°, and
#' the departure from symmetry will be larger. If there is no K-next,
#' \code{max_next_alignment_deg} is \code{NA}.
#'
#' @examples
#' k_next_old(sample_curvature) |>
#'   dplyr::select(-all_of(join_fields))
#'
#' @family K-next
#'
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom dplyr mutate
#' @importFrom dplyr rename
#' @importFrom dplyr rows_update
#' @importFrom dplyr inner_join
#' @importFrom dplyr semi_join
#' @importFrom dplyr anti_join
#' @importFrom dplyr bind_rows
#' @importFrom dplyr bind_cols
#' @importFrom dplyr slice
#' @importFrom dplyr slice_min
#' @importFrom dplyr slice_max
#' @importFrom dplyr sample_n
#' @importFrom dplyr summarize
#' @importFrom dplyr n
#' @importFrom tidyselect all_of
#'
#' @export
k_next_old <- function(exam_curvature, just_points = F) {

  # add an identifier flag to the cleansed curvature data
  source_dat <- exam_curvature |>
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
    # add k_max to tracking data frame: need both source and adjacent points
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


# Function to make a matrix square
make_square <- function(mat) {
  # Get the number of rows and columns
  n_rows <- nrow(mat)
  n_cols <- ncol(mat)

  # Check if the matrix is already square
  if (n_rows == n_cols) return(mat)

  # Determine the size of the new square matrix
  max_dim <- max(n_rows, n_cols)

  # Create a new square matrix with dimensions max_dim x max_dim, filled with NA
  square_mat <- matrix(NA, nrow = max_dim, ncol = max_dim)

  # Copy the original matrix into the new square matrix
  square_mat[1:n_rows, 1:n_cols] <- mat

  # Transfer row and column names?

  return(square_mat)
}


################
## Find Peaks
################

#' Find the peak points of a square numeric matrix.
#'
#' @description
#' Finds the peak points (i.e., max values) and the points that descend from
#' them in an \emph{n}x\emph{n} matrix.
#'
#' @param input_matrix
#' A square numeric matrix. If input_matrix is not square, this function will
#' introduce rows or columns of \code{NA}, as necessary.
#'
#' @return
#' A numeric matrix with the same dimensions as \code{input_matrix}. The values
#' in the matrix identify the peaks and the points that descend from them. The
#' highest peak is identified by the point(s) equal to 1, and the points that
#' descend from that peak are identified by -1. If there is a next highest peak,
#' it is identified by the point(s) equal to 2, and the points that descend from
#' it are identified by -2. And so on until all of the points in
#' \code{input_matrix} have been assigned to a peak.
#'
#' Any points that was \code{NA} in the \code{input_matrix} will be assigned 0.
#'
#' @examples
#' test_m <- matrix(c(12, 12, 11, 12, 12,
#'                    14, 14, 12, 13, 13,
#'                    14, 15, 13, 14, 13,
#'                    13, 14, 13, 14, 12,
#'                    12, 12, 12, 12, 12),
#'                  nrow = 5, byrow = T)
#' find_peaks(test_m)
#'
#' @family K-next
#'
#' @export
find_peaks <- function(input_matrix) {

  input_matrix <- make_square(input_matrix)

  n <- nrow(input_matrix)
  result <- matrix(0, n, n)
  peak_counter <- 1

  # Function to get neighbors
  get_neighbors <- function(x, y, n) {
    neighbors <- list()
    if (x > 1) neighbors <- append(neighbors, list(c(x-1, y)))
    if (x < n) neighbors <- append(neighbors, list(c(x+1, y)))
    if (y > 1) neighbors <- append(neighbors, list(c(x, y-1)))
    if (y < n) neighbors <- append(neighbors, list(c(x, y+1)))
    return(neighbors)
  }

  # Function to perform flood fill
  flood_fill <- function(input_matrix, result, x, y, peak_value, peak_label) {
    stack <- list(c(x, y))
    n <- nrow(input_matrix)

    while (length(stack) > 0) {
      point <- stack[[1]]
      stack <- stack[-1]
      px <- point[1]
      py <- point[2]

      if (result[px, py] == 0 && !is.na(input_matrix[px, py])) {
        result[px, py] <- ifelse(input_matrix[px, py] == peak_value, peak_label, -peak_label)

        neighbors <- get_neighbors(px, py, n)
        for (neighbor in neighbors) {
          nx <- neighbor[1]
          ny <- neighbor[2]
          if (result[nx, ny] == 0 && !is.na(input_matrix[nx, ny]) && input_matrix[nx, ny] <= input_matrix[px, py]) {
            stack <- append(stack, list(c(nx, ny)))
          }
        }
      }
    }
    return(result)
  }

  while (TRUE) {
    # Find the highest unprocessed peak
    max_value <- max(input_matrix * (result == 0), na.rm = TRUE)
    if (max_value == -Inf) break

    peak_pos <- which(input_matrix == max_value & result == 0, arr.ind = TRUE)
    if (nrow(peak_pos) == 0) break
    px <- peak_pos[1, 1]
    py <- peak_pos[1, 2]

    result <- flood_fill(input_matrix, result, px, py, max_value, peak_counter)

    peak_counter <- peak_counter + 1
  }

  return(result)
}


############
## K-next
############

#' Find the peaks on a curvature map
#'
#' @description
#' Finds the points of maximum dioptric power on a curvature map, including
#' K-max, K-next, and any other peaks.

#' @param curvature_m
#' A square numeric matrix. The rows must be labeled with their respective
#' \emph{y}-axis positions, and the columns with their respective \emph{x}-axis
#' positions. The values of the matrix must be the dioptric power at each
#' (\emph{x}, \emph{y}) position.
#' @param ignore_singletons
#' A Boolean. \code{TRUE} to ignore peaks that contain just one point.
#'
#' @return
#' Returns a data frame that contains one row for each peak on the curvature
#' map. The first few columns of the data frame will contain the
#' \code{join_fields}. The remaining columns describe the characteristics of
#' each peak:
#'
#' \describe{
#'  \item{peak_id}{The ordinal of the peak, ranked by the dioptric power of the
#'  peak's apex. The row with \code{peak_id == 1} will be K-max and the row with
#'  \code{peak_id == 2} will be K-next.}
#'  \item{x}{The peak's \emph{x}-axis position. See details.}
#'  \item{y}{The peak's \emph{y}-axis position.}
#'  \item{power}{The dioptric power of the peak.}
#'  \item{pct_coverage}{The percentage of points on the cornea that descend
#'  from the peak.}
#'  \item{r}{The radial distance of the peak from the vertex (in mm).}
#'  \item{theta}{The angle of the peak from the 3 o'clock position (in
#'  radians).}
#'  \item{angle}{The angle of the peak from the 3 o'clock position (in
#'  degrees).}
#'  \item{alignment_deg}{The angle (in degrees) between a peak and the
#'  projection of K-max in its opposite direction. See details.}
#'  \item{distance_to_max}{The Euclidean distance (in mm) between the peak and
#'  K-max.}
#'  \item{diff_power}{The difference in power (in diopters) between the peak
#'  and K-max peak.}
#' }
#'
#' @details
#' Oftentimes, a peak is comprised of several points with the same dioptric
#' power. When that happens, the (\emph{x}, \emph{y}) position of the peak is
#' calculated as the centroid of those points.
#'
#' \code{alignment_deg} quantifies departure from symmetry. When K-max and
#' another peak (particularly K-next) are symmetrically aligned, their radial
#' axes are 180° apart. As such, \code{alignment_deg} is calculated as 180°
#' \emph{minus} the smallest angle, so when K-max and another peak are 180°
#' apart, the result is 0°. Generally, the smallest angle between the K-max and
#' another peak will be less than 180°, and the departure from symmetry will be
#' larger.
#'
#' @examples
#' interpolate_measurements(sample_curvature) |>
#'   dplyr::mutate(power = anterior_power(z)) |>
#'   dplyr::select(x, y, power) |>
#'   reshape2::acast(y ~ x, value.var = 'power') |>
#'   k_next()
#'
#' @family K-next
#'
#' @importFrom sf st_centroid
#' @importFrom dplyr bind_rows
#' @importFrom dplyr mutate
#'
#' @export
k_next <- function(curvature_m, ignore_singletons = T) {

  if (is.null(curvature_m)) {
    return(NA)

  } else {
    # find the peaks
    curvature_peaks <- find_peaks(curvature_m)

    # curvature_peaks may have gotten resized, so fix that
    curvature_peaks <- curvature_peaks[1:nrow(curvature_m), 1:ncol(curvature_m)]

    # count the number of non-zero points, used as the denominator for
    # pct_coverage
    cnt_all_points <- sum(curvature_peaks != 0)

    # create a data frame to collate summary data for each peak
    peaks_summary <- data.frame()

    # count of peaks
    peak_id <- 0

    # get the centroid location and percentage coverage for each peak
    for (i in 1:max(curvature_peaks)) {
      # find the point(s) that comprise the peak
      # peak_location will contain the points that comprise the peak
      current_peak_apex <- arrayInd(which(curvature_peaks == i), dim(curvature_peaks))

      # peak_location contains the peak and all of the points that descend from it
      current_peak_all <- arrayInd(which(abs(curvature_peaks) == i), dim(curvature_peaks))

      if (!ignore_singletons | nrow(current_peak_all) > 1) {
        # increment the identifier
        peak_id <- peak_id + 1

        # turn row and col location(s) into x- and y-axis location(s)
        colnames(curvature_peaks) <- colnames(curvature_m)
        rownames(curvature_peaks) <- rownames(curvature_m)

        y <- rownames(curvature_m)[current_peak_apex[, 1]] |> as.numeric()
        x <- colnames(curvature_m)[current_peak_apex[, 2]] |> as.numeric()

        # all of the peak points will have the same power
        peak_power <- curvature_m[current_peak_apex[1, 1], current_peak_apex[1, 2]]

        # create a data frame of just the peak's points to find their centroid
        current_peak_apex <- data.frame(x,
                                  y,
                                  segment_id = i,
                                  contour = peak_power)

        # determine whether the peak has enough points to use {sf}
        if (nrow(current_peak_apex) > 4) {
          peak_polygon <- contour_to_sf_polygon(current_peak_apex)

          peak_centroid <- sf::st_centroid(peak_polygon) |>
            unlist(use.names = F)
        } else {
          # otherwise, just take the average of x and y as the location
          peak_centroid <- c(mean(x), mean(y))
        }

        # calculate the coverage as the percentage of all points
        peak_coverage <- nrow(current_peak_all) / cnt_all_points

        # create a one-row data frame to describe the peak
        peak_dat <- data.frame(
          peak_id,
          x = peak_centroid[1],
          y = peak_centroid[2],
          power = peak_power,
          pct_coverage = peak_coverage |> round(4)
        )

        # add the peak to the results data frame
        peaks_summary <- peaks_summary |>
          dplyr::bind_rows(peak_dat)
        }
    }

    # calculate radial position (r, theta) and Cartesian angle (angle)
    peaks_summary <- peaks_summary |>
      dplyr::mutate(r = sqrt(x ^ 2 + y ^ 2),
             theta = atan2(y, x) |> within_2pi(),
             angle = rad_to_deg(theta) |> within_360() |> round(1)
      )

    # calculate alignment with K-max
    peaks_summary$alignment_deg <- 180 - smallest_angle(peaks_summary$angle[1], peaks_summary$angle)
    peaks_summary$alignment_deg[1] <- NA

    # calculate distance to K-max
    peaks_summary$distance_to_max <- sqrt((peaks_summary$x[1] - peaks_summary$x) ^ 2 + (peaks_summary$y[1] - peaks_summary$y) ^ 2)
    peaks_summary$distance_to_max[1] <- NA

    # calculate difference in power from K-max
    peaks_summary$diff_power <- peaks_summary$power[1] - peaks_summary$power
    peaks_summary$diff_power[1] <- NA


    return(peaks_summary)
  }
}


#####################
## antipodal_power
#####################

#' Dioptric power at the antipode of K-max
#'
#' @description
#' Returns the dioptric power at the antipode of K-max.
#'
#' @param exam_curvature
#' A data frame containing, at a minimum, the \code{x}- and \code{y}-axis
#' position of each curvature (i.e., radius of curvature) \code{measurement}.
#' One row for each \code{measurement}.
#' @param exam_k_max
#' A data frame containing the [IQeyes::join_fields] and the \code{x}, \code{y},
#' and \code{measurement} of K-max.
#' @param interp
#' A Boolean. \code{TRUE} to interpolate measurements. See details.
#'
#' @return
#' A one-row data frame containing the [IQeyes::join_fields] and the following
#' columns:
#'
#' \describe{
#'  \item{k_max}{Dioptric power of K-max.}
#'  \item{opposite_power}{Dioptric power of the antipode.}
#'  \item{x}{The \emph{x}-axis position of the measurement. See details.}
#'  \item{y}{The \emph{y}-axis position of the measurement.}
#'  \item{opposite_x}{The \emph{x}-axis position of the antipode. See details.}
#'  \item{opposite_y}{The \emph{y}-axis position of the antipode.}
#'  \item{euch_dist}{Euchlidean distance of the selected \code{measurement} from
#'  the location of "exact" antipode of K-max. See details.}
#' }
#'
#' @details
#' The antipode of K-max is the point directly opposite of K-max on the
#' curvature map. Since measurement data for that \emph{exact} position is
#' rarely available, this function finds the point that's nearest to the
#' antipode and returns its position and power.
#'
#' The [IQeyes::join_fields] in the returned data frame come from the
#' \code{exam_k_max}.
#'
#' @family K-next
#'
#' @importFrom dplyr select
#' @importFrom dplyr mutate
#' @importFrom dplyr cross_join
#' @importFrom dplyr slice_min
#' @importFrom tidyselect all_of
#' @importFrom tidyselect starts_with
#'
#' @export
antipodal_power <- function(exam_curvature, exam_k_max, interp = F) {

  if (interp) {
    z_dat <- exam_curvature |>
      interpolate_measurements()
  } else {
    z_dat <- exam_curvature |>
      select(x, y, measurement) |>
      rename(z = measurement)
  }

  # establish the antipode of K-max
  # K-max is at (x, y), the antiode is at (-x, -y)
  exam_k_max <- exam_k_max |>
    dplyr::select(tidyselect::all_of(join_fields), x, y, measurement) |>
    dplyr::mutate(opposite_x = -x,
                  opposite_y = -y,
                  k_max = anterior_power(measurement)) |>
    dplyr::select(tidyselect::all_of(join_fields), tidyselect::starts_with('opposite'), k_max)

  if (nrow(exam_k_max) == 0) {
    warning('K-max not found.')
    return(NA)
  }

  # find the measurement that's closest to the antipodal position of K-max
  antipodal_power <- exam_k_max |>
    dplyr::cross_join(z_dat) |>
    dplyr::mutate(euch_dist = sqrt((x - opposite_x) ^ 2 + (y - opposite_y) ^ 2)) |>
    dplyr::slice_min(order_by = euch_dist, n = 1, with_ties = F, by = tidyselect::all_of(join_fields)) |>
    dplyr::mutate(opposite_power = anterior_power(z)) |>
    dplyr::select(all_of(join_fields), k_max, opposite_power, x, y, opposite_x, opposite_y, euch_dist)

  return(antipodal_power)

}
