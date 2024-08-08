## Color Scale ##

# the absolute power scale that is used to determine the contour levels for
# curvature maps. The scale is not linear
absolute_scale <- c(seq(10, 30, by = 2.5),
                    seq(30, 46, by = 0.5),
                    seq(46, 50, by = 1),
                    seq(50, 90, by = 2.5)) |> unique()


# base version of the color scale used for curvature maps
# the colors are taken from a screen shot of Pentacam's American Style scale,
# using PowerPoint's eye dropper
# the radius numbers are based on the same constant used by anterior_power()
absolute_color_scale <- data.frame(
  value = c(10, 20, 30, 32, 34, 36, 38, 40, 42, 44, 46, 50, 60, 70, 80, 90),
  #radius = c(34.5, 22.5,  10.5,  10.0,  9.6,  9.2,  8.8,  8.4,  8.0,  7.6,  7.2,  6.8,  6.2,  5.4,  4.6,  3.8),
  radius = c(33.75, 16.88, 11.25, 10.55, 9.93, 9.38, 8.88, 8.44, 8.04, 7.67, 7.34, 6.75, 5.62, 4.82, 4.22, 3.75),
  color = c(
    '#B428DC',  # 10
    '#8C28B4',  # 20
    '#642896',  # 30
    '#28008C',  # 32
    '#283CDC',  # 34
    '#1478F0',  # 36
    '#00B9FC',  # 38
    '#00F0FC',  # 40
    '#5AE600',  # 42
    '#FFFF00',  # 44
    '#DCC800',  # 46
    '#DC9600',  # 50
    '#DC5000',  # 60
    '#B40000',  # 70
    '#A03232',  # 80
    '#C86464'   # 90
  )
)

# augment the color scale with power values
absolute_color_scale$power <- anterior_power(absolute_color_scale$radius)

# extract the range of n colors from the absolute color scale that spans the min
# and max of a given z-axis (presumed to be in diopters)
# n is the number of colors to include in the scale. If n is not provided,
# provide one color for each member of z
plot_color_scale <- function(z, n = NA) {
  # filter the color scale to values between the min and max measurements
  min_max_color_scale <- absolute_color_scale |>
    dplyr::filter(dplyr::between(value, min(z, na.rm = T), max(z, na.rm = T)))

  if (is.na(n))
    n <- length(z)

  colorRampPalette(min_max_color_scale$color)(n)
}

# get the color ramp for each of the linear segments of the full scale
full_color_scale <- c(
  plot_color_scale(seq(10, 30, by = 2.5)),
  plot_color_scale(seq(30, 46, by = 0.5)),
  plot_color_scale(seq(46, 50, by = 1)),
  plot_color_scale(seq(50, 90, by = 2.5))
) |> unique()

# redefine the absolute color scale based on the full fidelity of the
# absolute scale. calls to plot_color_scale() will come closer to matching the
# Pentacam's output
absolute_color_scale <- data.frame(value = absolute_scale,
                                   color = full_color_scale)


################
## plot_scale
################

#' Contour breaks needed to plot an anterior curvature map
#' @description
#'
#' Returns the subset of contour breaks from \code{absolute_scale} required to
#' plot a curvature map of the values in \code{z}.
#'
#' @param z
#' A numeric vector containing power (in diopters). Currently
#' works for \emph{anterior} surfaces only.
#' @return
#' A numeric vector containing the breaks on the \emph{z}-axis scale (in
#' diopters).
#'
#' @details
#' The contour levels on a curvature map are defined by the breaks in
#' \code{absolute_scale}, which is based on the Pentacam's American Style scale
#' for the dioptric power of anterior curvature maps.
#'
#' @examples
#' sample_curvature$measurement |>
#'   anterior_power() |>
#'   plot_scale()
#'
#' @family plots
#'
#' @importFrom dplyr between
#'
#' @export
plot_scale <- function(z) {
  scale_indices <-
    which(dplyr::between(absolute_scale, min(z, na.rm = T), max(z, na.rm = T)))

  if (min(scale_indices) > 1) {
    if (min(z, na.rm = T) < absolute_scale[scale_indices[1]]) {
      scale_indices <- c(min(scale_indices) - 1, scale_indices)
    }
  }

  if (max(scale_indices) < length(absolute_scale)) {
    if (max(z, na.rm = T) > absolute_scale[scale_indices[length(scale_indices)]]) {
      scale_indices <- c(scale_indices, max(scale_indices) + 1)
    }
  }

  return(absolute_scale[scale_indices])
}


####################
## curvature_plot
####################

#' Generate the curvature plot of the anterior surface of a cornea
#' @description
#' Returns (but does not render) the curvature plot of the anterior surface of a
#' cornea by interpolating measurements from \code{COR-PWR.CSV}.
#'
#' @param exam_curvature
#' A data frame with the same structure as [IQeyes::sample_curvature],
#' containing one row for each curvature (i.e., radius of curvature)
#' \code{measurement}.
#' @param interp
#' A Boolean. \code{TRUE} to interpolate measurements.
#' @param grid
#' A Boolean. \code{TRUE} to include a background grid.
#' @param labels
#' A Boolean. \code{TRUE} to label the \emph{measured} curvature points. Only
#' works when data are interpreted.
#' @param titles
#' A Boolean. \code{TRUE} to include plot title(s).
#' @param axes
#' A Boolean. \code{TRUE} to include \emph{x} and \emph{y} axes.
#' @param legend
#' A Boolean. \code{TRUE} to include a legend/key.
#' @param truncate_circle
#' A Boolean. \code{TRUE} to limit the curvature map to a diameter of
#' \code{circle_dia}.
#' @param circle_dia
#' A numeric value. The diameter of the plotted region of the matrix. Used only
#' if \code{truncate_circle} is \code{TRUE}.
#' @param points
#' A Boolean. \code{TRUE} to show points on the plot at the \emph{x}- and
#' \emph{x}-coordinates of the \emph{measured} points. Only works when data are
#' interpreted.
#' @param greyscale
#' A Boolean. \code{TRUE} to render the contours in greyscale.
#' @param theme_minimal
#' A Boolean. \code{TRUE} to use ggplot's minimal theme (e.g., no background
#' color).
#' @param coord_fixed
#' A Boolean. \code{TRUE} to ensure the \emph{x} and \emph{y} axes remain fixed
#' (i.e., equal), regardless of differences in the available plot area.
#'
#' @return
#' A ggplot object.
#'
#' @details
#' This function currently works for anterior surfaces only because
#' \code{plot_scale()} presumes the curvature maps uses the Pentacam's American
#' scale, which is designed for the dioptric power of anterior surfaces.
#'
#' @examples
#' curvature_plot(sample_curvature)
#'
#' @family plots
#'
#' @importFrom dplyr select
#' @importFrom dplyr mutate
#' @importFrom dplyr filter
#' @importFrom tidyselect all_of
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_contour_filled
#' @importFrom ggplot2 guides
#' @importFrom ggplot2 guide_legend
#' @importFrom ggplot2 geom_spoke
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 geom_label
#' @importFrom ggplot2 theme_void
#' @importFrom ggplot2 scale_fill_grey
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 coord_fixed
#' @importFrom ggplot2 theme_minimal
#' @importFrom ggforce geom_circle
#'
#' @export
curvature_plot <-
  function(exam_curvature,
           interp = T,
           grid = F,
           labels = F,
           titles = T,
           axes = T,
           legend = T,
           truncate_circle = T,
           circle_dia = 9,
           points = F,
           greyscale = F,
           theme_minimal = T,
           coord_fixed = T) {

    if (nrow(exam_curvature) == 0) warning('No data in exam_curvature.')

    # establish join_fields (and surface) from exam_curvature, which will be
    # used for the chart's labels
    exam_record <- exam_curvature[1, ] |>
      dplyr::select(tidyselect::all_of(join_fields), surface)

    if (nrow(exam_record) > 1) warning('More than one exam record contained in exam_curvature.')

    # interpolate data for the entire surface (this returns x, y, z)
    if (interp) {
      plot_dat <- exam_curvature |>
        interpolate_measurements()
    } else {
      plot_dat <- exam_curvature |>
        select(x, y, measurement) |>
        rename(z = measurement)
    }

    # convert curvature radius (R) to power in diopters
    if (exam_record$surface == 'FRONT') {
      plot_dat <- plot_dat |>
        dplyr::mutate(power = round(anterior_power(z), 1))
      # calculate power
      exam_curvature <- exam_curvature |>
        dplyr::mutate(power = round(anterior_power(measurement), 1))
    } else if (exam_record$surface == 'BACK') {
      plot_dat <- plot_dat |>
        dplyr::mutate(power = round(posterior_power(z), 1))
      exam_curvature <- exam_curvature |>
        dplyr::mutate(power = round(posterior_power(measurement), 1))
    }

    if (truncate_circle) {
      radius <- circle_dia / 2
      plot_dat <- plot_dat |>
        dplyr::filter(sqrt((x ^ 2) + (y ^ 2)) <= radius)
      exam_curvature <- exam_curvature |>
        dplyr::filter(sqrt((x ^ 2) + (y ^ 2)) <= radius)
    }

    breaks_vector <- plot_scale(plot_dat$power)

    p <- ggplot2::ggplot() +
        ggplot2::geom_contour_filled(data = plot_dat,
                                     ggplot2::aes(x = x, y = y, z = power),
          # need breaks to match with color scale
          breaks = absolute_scale,
          show.legend = legend,
          na.rm = T
        ) +
        # this shows the curvature scale from high to low (top to bottom)
        ggplot2::guides(fill = ggplot2::guide_legend(reverse = T, title = 'Curvature'))

    if (grid) {
      # a simple radial background
      p <- p +
        ggplot2::geom_spoke(data = radii_degrees, ggplot2::aes(x = 0, y = 0, angle = deg_to_rad(angle), radius = ring_diam), color = 'grey') +
        ggforce::geom_circle(data = radii_degrees, ggplot2::aes(x0 = 0, y0 = 0, r = ring_diam), color = 'grey')
    }

    if (points) {
      if (!interp) {
        warning('Can plot points only when curvature data has been interpolated.')
      } else {
        p <- p +
          ggplot2::geom_point(data = exam_curvature, ggplot2::aes(x = x, y = y))
      }
    }

    if (titles) {
      p <- p +
        ggplot2::labs(
          title = 'Curvature (Power)',
          subtitle = paste0('Patient: ', exam_record$pat_id, ', Eye: ', exam_record$exam_eye, ', Surface: ', exam_record$surface
          )
        )
    }

    if (labels) {
      if (!interp) {
        warning('Can plot labels only when curvature data has been interpolated.')
      } else {
        p <- p +
          ggplot2::geom_label(data = exam_curvature, ggplot2::aes(x = x, y = y, label = power), size = 2)
      }
    }

    if (!axes) {
      p <- p +
        ggplot2::theme_void()
    }

    if (greyscale) {
      p <- p +
        ggplot2::scale_fill_grey()
    } else {
      p <- p +
        ggplot2::scale_fill_manual(values = plot_color_scale(plot_dat$power, length(breaks_vector)))
    }

    if (coord_fixed) {
      p <- p +
        ggplot2::coord_fixed()
    }

    if (theme_minimal) {
      p <- p +
        ggplot2::theme_minimal()
    }

    return(p)
  }


#################
## plot_matrix
#################

#' Plot the detailed matrix data
#'
#' @description
#' Returns (but does not render) a filled contour plot of a matrix dataset.
#'
#' @param exam_matrix
#' A numeric matrix containing the values to plot. The row and column names of
#' the matrix must contain y- and x-axis positions, respectively, of the z-axis
#' values.
#' @param color_scale
#' A data frame containing the \code{breaks} and \code{hex_color} values for the
#' contours of the plot.
#' @param scale_units
#' A character string containing the units of measure for the z-axis values.
#' Used to label the legend.
#' @param invert_order
#' A Boolean. \code{FALSE} (the default) to retain the order of the
#' \code{color_scale}. \code{TRUE} to invert the order.
#' @param truncate_circle
#' A Boolean. \code{TRUE} to limit the plotted region of the matrix to a circle
#' of diameter \code{circle_dia}.
#' @param circle_dia
#' A numeric value. When \code{truncate_circle} is \code{TRUE}, the diameter of
#' the plotted region of the matrix.
#' @param axes
#' A Boolean. \code{TRUE} to display the \emph{x} and \emph{y} axes on the plot.
#' @param legend
#' A Boolean. \code{TRUE} to display a legend/key.
#' @param show_arrows
#' A Boolean. \code{TRUE} to display the point-to-point change in z-axis values
#' as arrows.
#' @param arrows_legend
#' A Boolean. \code{TRUE} to display a legend/key for arrows.
#' @param theme_minimal
#' A Boolean. \code{TRUE} to use ggplot's minimal theme (e.g., no background
#' color).
#' @param coord_fixed
#' A Boolean. \code{TRUE} to ensure the \emph{x} and \emph{y} axes remain fixed
#' (i.e., equal), regardless of differences in the extents of the \emph{x} and
#' \emph{y} axes.
#'
#' @return
#' A ggplot object.
#'
#' @details
#' To plot anterior axial/sagittal curvature, the z-axis values must be in
#' diopters. Use [IQeyes::anterior_power] to convert radii values to diopters.
#'
#' Posterior curvature has not yet been implemented: It requires only the
#' definition of a color scale.
#'
#' Tangential curvature has not yet been implemented.
#'
#' @family plots
#'
#' @importFrom reshape2 melt
#' @importFrom dplyr filter
#' @importFrom dplyr mutate
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_contour_filled
#' @importFrom ggplot2 guides
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 guide_legend
#' @importFrom ggplot2 coord_fixed
#' @importFrom ggplot2 theme_minimal
#' @importFrom ggplot2 theme_void
#'
#' @importFrom pracma gradient
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr arrange
#' @importFrom metR geom_arrow
#' @importFrom metR scale_mag
#' @importFrom ggplot2 labs
#'
#' @export
plot_matrix <- function(exam_matrix,
                        color_scale = c("elevation_colors", "pachymetry_colors", "curvature_colors"),
                        scale_units = "Measurement",
                        invert_order = F,
                        truncate_circle = T,
                        circle_dia = 9,
                        axes = T,
                        legend = T,
                        show_arrows = F,
                        arrows_legend = T,
                        theme_minimal = T,
                        coord_fixed = T) {

  # exam_matrix <- CUR_front[[cur_idx]] |> anterior_power()
  # color_scale <- curvature_colors

  # Convert the matrix to a data frame in long format
  exam_df <- reshape2::melt(exam_matrix, varnames = c('y', 'x'), value.name = "z") |>
    dplyr::filter(!is.na(z))

  if (truncate_circle) {
    # filter the data frame to the Pentacam's viewable range
    exam_df <- exam_df |>
      dplyr::mutate(r = sqrt(x ^ 2 + y ^ 2)) |>
      dplyr::filter(r <= circle_dia / 2)
  }

  # fit the range of the color scale to the min and max of the z-axis values
  hist_dat <- hist(exam_df$z, breaks = color_scale$breaks, plot = F)

  #hist_dat$breaks[hist_dat$counts > 0]
  #min(exam_df$z, na.rm = T)
  #max(exam_df$z, na.rm = T)

  scale_indices <- which(hist_dat$counts > 0)

  if (min(scale_indices) > 1) {
    if (min(exam_df$z, na.rm = T) < hist_dat$breaks[scale_indices[1]]) {
      scale_indices <- c(min(scale_indices) - 1, scale_indices)
    }
  }

  if (max(scale_indices) < length(hist_dat$breaks))
    if (max(exam_df$z, na.rm = T) > hist_dat$breaks[max(scale_indices)]) {
      scale_indices <- c(scale_indices, max(scale_indices) + 1)
    }

  color_scale <- color_scale |>
    dplyr::filter(breaks >= hist_dat$breaks[scale_indices] |> min(),
                  breaks <= hist_dat$breaks[scale_indices] |> max()
    )

  # generate the base plot
  p <- ggplot2::ggplot() +
    ggplot2::geom_contour_filled(data = exam_df,
                        ggplot2::aes(x = x, y = y, z = z),
                        # need breaks to match with color scale
                        breaks = color_scale$breaks,
                        show.legend = legend,
                        na.rm = T
    ) +
    ggplot2::scale_fill_manual(values = color_scale$hex_color) +
    ggplot2::guides(fill = ggplot2::guide_legend(reverse = invert_order, title = scale_units))

  # include x and y axes
  if (!axes) {
    p <- p +
      ggplot2::theme_void()
  }

  # used fixed coordinates
  if (coord_fixed) {
    p <- p +
      ggplot2::coord_fixed()
  }

  # use a minimal theme
  if (theme_minimal) {
    p <- p +
      ggplot2::theme_minimal()
  }

  # show peak arrows
  if (show_arrows) {
    # find the detailed peak data
    exam_peaks_m <- exam_matrix |>
      find_peaks()

    # Compute slope gradients (using finite differences to estimate slopes)
    gradients <- pracma::gradient(exam_matrix)

    # Extract x and y gradients
    grad_x <- gradients[[1]]
    grad_y <- gradients[[2]]

    # Compute gradient magnitude and direction (in radians)
    grad_mag <- sqrt(grad_x ^ 2 + grad_y ^ 2)
    grad_dir <- atan2(grad_y, -grad_x)

    # convert one of the matrices to a data frame
    grad_df <- grad_mag |>
      as.data.frame() |>
      tibble::rownames_to_column('y') |>
      tidyr::pivot_longer(cols = -'y', names_to = 'x', values_to = 'magnitude') |>
      dplyr::mutate(x = as.numeric(x), y = as.numeric(y)) |>
      # put rows of the data frame in the same order that as.vector() function
      # uses when it converts matrices
      dplyr::arrange(x, desc(y))

    # add the other matrices as vector columns in the data frame
    grad_df$matrix_measure = as.vector(exam_matrix)
    grad_df$direction = as.vector(grad_dir)
    grad_df$peak = as.vector(exam_peaks_m) |> abs()

    # after the data frame has been collated, remove NAs
    grad_df <- grad_df |>
      dplyr::filter(!is.na(magnitude))

    if (truncate_circle) {
      # filter the data frame to the Pentacam's viewable range
      grad_df <- grad_df |>
        dplyr::mutate(r = sqrt(x ^ 2 + y ^ 2)) |>
        dplyr::filter(r <= circle_dia / 2)
    }

    p <- p +
      metR::geom_arrow(data = grad_df, aes(x, y,
                                     dx = magnitude * cos(direction),
                                     dy = magnitude * sin(direction),
                                     color = as.factor(peak)),
                 skip = 2,
                 start = 0, direction = 'ccw',
                 show.legend = arrows_legend,
                 arrow.length = 0.75) +
      metR::scale_mag(name = 'Scale') +
      ggplot2::labs(color = 'Peak')

    # if no legend, remove legend/labels for color and scale
    if (!arrows_legend) p <- p + ggplot2::labs(color = NULL, mag = NULL)

  }

  return(p)

}
