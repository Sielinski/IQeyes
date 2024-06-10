## Color Scale ##

# base version of the color scale used for curvature maps
absolute_color_scale <- data.frame(
  value = c(10, 20, 30, 32, 34, 36, 38, 40, 42, 44, 46, 50, 60, 70, 80, 90),
  radius = c(34.5, 22.5, 10.5, 10.0, 9.6, 9.2, 8.8, 8.4, 8.0, 7.6, 7.2, 6.8, 6.2, 5.4, 4.6, 3.8),
  color = c(
    '#B428E1',
    '#8C28B4',
    '#6E28A0',
    '#28008C',
    '#233278',
    '#1978F0',
    '#00BEFF',
    '#00F0FF',
    '#5AE600',
    '#FFFF00',
    '#DCC800',
    '#DC9600',
    '#DC5000',
    '#B40000',
    '#A03232',
    '#C86464'
  )
)

# augment the scale with power values
absolute_color_scale$power <- anterior_power(absolute_color_scale$radius)

# the absolute power scale that is used for curvature maps is not linear
absolute_scale <- c(seq(10, 30, by = 2.5),
                    seq(30, 46, by = 0.5),
                    seq(46, 50, by = 1),
                    seq(50, 90, by = 2.5)) |>
  unique()


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
# pentacam's output
absolute_color_scale <- data.frame(value = absolute_scale,
                                   color = full_color_scale)


# the absolute scale when shown in terms of curvature radius
absolute_curvature_scale <- c(seq(34.5, 10.5, by = -3),
                              seq(10.5, 10, by = -0.125),
                              seq(10, 6.8, by = -0.1),
                              seq(6.8, 6.2, by = -0.15),
                              seq(6.2, 3.8, by = -0.2)) |>
  unique()


# extract the breaks from the absolute scale that span the min and max of the
# relevant z-axis (presumed to be in diopters)
plot_scale <- function(z) {
  scale_indices <-
    which(dplyr::between(absolute_scale, min(z, na.rm = T), max(z, na.rm = T)))

  if (min(z, na.rm = T) < absolute_scale[scale_indices[1]]) {
    if (min(scale_indices) > 1)
      scale_indices <- c(min(scale_indices) - 1, scale_indices)
  }

  if (max(z, na.rm = T) > absolute_scale[scale_indices[length(scale_indices)]]) {
    if (max(scale_indices) < length(absolute_scale))
      scale_indices <- c(scale_indices, max(scale_indices) + 1)
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
#' @param exam_curvature
#' A data frame containing one row for each curvature \code{measurement} and the
#' same columns as [IQeyes::sample_curvature].
#' @param basic
#' A Boolean. \code{TRUE} to forgo contours and colors.
#' @param grid
#' A Boolean. \code{TRUE} to include a background grid.
#' @param labels
#' A Boolean. \code{TRUE} to label the \emph{measured} curvature points.
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
#' A number. \code{TRUE} to include a legend/key.
#' @param points
#' A Boolean. \code{TRUE} to points on the plot at the \emph{x}- and \emph{x}-
#' coordinates of the \emph{measured} points.
#' @param greyscale
#' A Boolean. \code{TRUE} to render the contours in greyscale.
#' @return
#' A ggplot object.
#' @details
#' This function doesn't currently work for posterior surfaces because the
#' \code{plot_scale} function presumes the color scale is on the "absolute
#' scale", which is designed the dioptric power of anterior surfaces.
#' @seealso
#' [IQeyes::plot_scale()]
#' @examples
#' curvature_plot(sample_curvature)
#'
curvature_plot <-
  function(exam_curvature,
           basic = F,
           grid = F,
           labels = T,
           titles = T,
           axes = T,
           legend = T,
           truncate_circle = T,
           circle_dia = 9,
           points = F,
           greyscale = F) {

    if (nrow(exam_curvature) == 0) warning('No data in exam_curvature.')

    # establish join_fields (and surface) from exam_curvature, which will be
    # used for the chart's labels
    exam_record <- exam_curvature |>
      dplyr::select(tidyselect::all_of(join_fields), surface) |>
      unique()

    if (nrow(exam_record) > 1) warning('More than one exam record contained in exam_curvature.')

    # interpolate data for the entire surface (this returns x, y, z)
    interp_df <- exam_curvature |>
      dplyr::select(x, y, measurement) |>
      dplyr::filter(!is.na(measurement)) |>
      interpolate_measurements()

        # convert curvature radius (R) to power in diopters
    if (exam_record$surface == 'FRONT') {
      interp_df <- interp_df |>
        dplyr::mutate(power = round(anterior_power(z), 1))
      curvature_dat <- exam_curvature |>
        dplyr::mutate(power = round(anterior_power(measurement), 1))
    } else if (exam_record$surface == 'BACK') {
      interp_df <- interp_df |>
        dplyr::mutate(power = round(posterior_power(z), 1))
      curvature_dat <- exam_curvature |>
        dplyr::mutate(power = round(posterior_power(measurement), 1))
    }

    if (truncate_circle) {
      radius <- circle_dia / 2
      interp_df <- interp_df |>
        dplyr::filter(sqrt((x ^ 2) + (y ^ 2)) <= radius)
      curvature_dat <- curvature_dat |>
        dplyr::filter(sqrt((x ^ 2) + (y ^ 2)) <= radius)
    }

    breaks_vector <- plot_scale(interp_df$power)

    p <- ggplot2::ggplot(interp_df, ggplot2::aes(x = x, y = y))

    if (!basic) {
      p <- p +
        ggplot2::geom_contour_filled(
          ggplot2::aes(z = power),
          # need breaks to match with color scale
          breaks = breaks_vector,
          show.legend = legend,
          na.rm = T
        ) +
        ggplot2::guides(fill = ggplot2::guide_legend(reverse = T, title = 'Curvature'))
    }

    if (grid) {
      # a simple radial background
      p <- p +
        ggplot2::geom_spoke(data = radii_degrees, ggplot2::aes(x = 0, y = 0, angle = deg_to_rad(angle), radius = ring_diam), color = 'grey') +
        ggforce::geom_circle(data = radii_degrees, ggplot2::aes(x0 = 0, y0 = 0, r = ring_diam), color = 'grey')
    }

    if (points) {
      p <- p +
        ggplot2::geom_point(data = curvature_dat, ggplot2::aes(x = x, y = y))
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
      p <- p +
        ggplot2::geom_label(data = curvature_dat, ggplot2::aes(x = x, y = y, label = power), size = 2)
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
        ggplot2::scale_fill_manual(values = plot_color_scale(interp_df$power, length(breaks_vector)))
    }

    return(p)
  }
