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

#' Find the peak points of a numeric matrix.
#'
#' @description
#' Finds the peak points (i.e., max values) and the points that descend from
#' them in a numeric matrix.
#'
#' @param input_matrix
#' A numeric matrix.
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

  n_rows <- nrow(input_matrix)
  n_cols <- ncol(input_matrix)
  result <- matrix(0, n_rows, n_cols)
  peak_counter <- 1

  # Function to get neighbors
  get_neighbors <- function(x, y, n_rows, n_cols) {
    neighbors <- list()
    if (x > 1) neighbors <- append(neighbors, list(c(x-1, y)))
    if (x < n_rows) neighbors <- append(neighbors, list(c(x+1, y)))
    if (y > 1) neighbors <- append(neighbors, list(c(x, y-1)))
    if (y < n_cols) neighbors <- append(neighbors, list(c(x, y+1)))
    return(neighbors)
  }

  # Function to perform flood fill
  flood_fill <- function(input_matrix, result, x, y, peak_value, peak_label) {
    stack <- list(c(x, y))
    n_rows <- nrow(input_matrix)
    n_cols <- ncol(input_matrix)

    while (length(stack) > 0) {
      point <- stack[[1]]
      stack <- stack[-1]
      px <- point[1]
      py <- point[2]

      if (result[px, py] == 0 && !is.na(input_matrix[px, py])) {
        result[px, py] <- ifelse(input_matrix[px, py] == peak_value, peak_label, -peak_label)

        neighbors <- get_neighbors(px, py, n_rows, n_cols)
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
