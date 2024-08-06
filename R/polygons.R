###########################
## contour_to_sf_polygon
##########################

#' Convert a contour to a polygon
#'
#' @description
#' Converts a contour to an sf polygon object.
#'
#' @param contour
#' A data frame containing one row for each point of the contour and the same
#' columns as [IQeyes::sample_contour]. If the contour needs to be scaled and/or
#' rotated, the necessary transforms must be performed prior to calling
#' \code{contour_to_sf_polygon()}.
#' @param required_points
#' A number indicating the minimum number of points required to convert a
#' segment of the contour to a feature of the polygon. Must be four or more.
#'
#' @return
#' An sf object. See details.
#'
#' @details
#' Each segment of the contour will become a feature of the polygon. There are
#' two exceptions, however: Segments that are fully contained within another
#' segment will \emph{subtracted} from the resulting polygon, and segments with
#' fewer than \code{required_points} will be disregarded.
#'
#' @examples
#' contour_to_sf_polygon(
#'   contour = get_contour(sample_curvature, interp = T, contour_power = 44.5)
#' )
#'
#' reference_contours |>
#'   dplyr::filter(cluster == 7) |>
#'   contour_to_sf_polygon()
#'
#' @family polygons
#'
#' @importFrom sf st_sf
#' @importFrom sf st_sfc
#' @importFrom sf st_as_sf
#' @importFrom sf st_cast
#' @importFrom sf st_contains
#' @importFrom sf st_union
#' @importFrom sf st_contains
#' @importFrom sf st_difference
#' @importFrom sf st_is_valid
#' @importFrom sf st_make_valid
#' @importFrom dplyr ungroup
#' @importFrom dplyr summarize
#' @importFrom dplyr filter
#' @importFrom dplyr pull
#' @importFrom dplyr select
#' @importFrom dplyr group_by
#'
#' @export
contour_to_sf_polygon <- function(contour, required_points = 4) {

  # contour <- get_contour(sample_curvature, interp = T, contour_power = 46)[1:4, ]

  if (required_points < 4) {
    warning('required_points must be greater than 4: resetting to 4.')
    required_points <- 4
  }

  # ensure that each segment has enough points to become a valid polygon
  valid_segments <- contour |>
    dplyr::ungroup() |>
    dplyr::summarize(cnt = dplyr::n(), .by = segment_id) |>
    dplyr::filter(cnt > required_points) |>
    dplyr::pull(segment_id)

  # extract just the valid segments and necessary columns
  df <- contour |>
    dplyr::filter(segment_id %in% valid_segments) |>
    dplyr::select(x, y, segment_id)

  if(length(valid_segments) == 0) {
    warning('None of the segments are large enough to convert to a polygon')

    empty_sf <- sf::st_sf(
      geometry = sf::st_sfc()
    )
    return(empty_sf)

  } else {
    # Convert the data frame to an sf object with polygons
    sf_segments <- df |>
      dplyr::select(x, y, segment_id) |>
      sf::st_as_sf(coords = c("x", "y"), crs = NA) |>
      dplyr::group_by(segment_id) |>
      dplyr::summarize(do_union = FALSE) |>
      sf::st_cast("POLYGON")

    # if necessary, repair the polygon
    for (i in 1:nrow(sf_segments)) {
      if (!sf::st_is_valid(sf_segments[i, ]))
        sf_segments[i, ] <- sf::st_make_valid(sf_segments[i, ])
    }

    # Initialize a list to hold the final polygons
    final_polygons <- list()

    # go through every polygon
    for (i in 1:nrow(sf_segments)) {
      current_polygon <- sf_segments[i, ]

      is_outer <- TRUE

      # check if the current polygon is enclosed by any other polygon
      for (j in 1:nrow(sf_segments)) {
        if (i != j) {
          test_polygon <- sf_segments[j, ]
          if (sf::st_contains(test_polygon, current_polygon, sparse = FALSE)) {
            is_outer <- FALSE
            break
          }
        }
      }

      # if the polygon is not enclosed by any other polygon, add it to the
      # list of outer polygons
      if (is_outer) {
        final_polygons[[length(final_polygons) + 1]] <- current_polygon
      }
    }

    # Combine all outer polygons
    combined_polygon <- final_polygons |>
      dplyr::bind_rows()

    # check if all segments have been combined
    if (nrow(combined_polygon) < nrow(sf_segments)) {
      # if not, identify the inner segment(s)
      inner_segments <- sf_segments$segment_id[!sf_segments$segment_id %in% combined_polygon$segment_id]

      # remove any identified inner segments
      for (i in seq_along(inner_segments)) {
        combined_polygon <- suppressWarnings(sf::st_difference(combined_polygon, sf_segments[inner_segments[i], ]))
      }
    }

    # get rid of any unnecessary columns
    combined_polygon <- combined_polygon |>
      as.data.frame() |>
      dplyr::select(geometry) |>
      sf::st_as_sf()

    # if necessary, repair the polygon
    # if (!sf::st_is_valid(combined_polygon)) combined_polygon <- sf::st_make_valid(combined_polygon)

    return(combined_polygon)
  }

}


########################
## silhouette_overlap
########################

#' Average overlap between the silhouettes of two contours
#'
#' @description
#' Identifies the overlapping area of two silhouettes, calculates the percentage
#' of each silhouette's area contained by the overlap, and returns the average
#' percentage overlap.
#'
#' Both \code{contour_A} and \code{contour_B} should be data frames containing
#' one row for each point of their respective contours and the same columns as
#' [IQeyes::sample_contour]. If the contours need to be scaled and/or
#' rotated, the necessary transforms must be performed prior to calling
#' \code{silhouette_overlap()}. See [IQeyes:scale_rotate].
#'
#' @param contour_A
#' A data frame containing contour A.
#' @param contour_B
#' A data frame containing contour B.
#' @param show_plot
#' A Boolean. \code{TRUE} to render a plot that illustrates the comparison
#' being made. To access the plot as an object, see [IQeyes:silhouette_plot].
#'
#' @return
#' The average percentage overlap between the two silhouettes.
#'
#' @examples
#' silhouette_overlap(
#'   contour_A = get_contour(sample_curvature, interp = T, contour_power = 44),
#'   contour_B = get_contour(sample_curvature, interp = T, contour_power = 44.5),
#'   show_plot = T
#' )
#'
#' @family polygons
#'
#' @importFrom sf st_intersection
#' @importFrom sf st_area
#' @importFrom ggplot2 labs
#'
#' @export
silhouette_overlap <- function(contour_A, contour_B, show_plot = F) {

  # convert contours to polygons
  poly_A <- contour_to_sf_polygon(contour_A)
  poly_B <- contour_to_sf_polygon(contour_B)

  # identify overlap/intersection
  intersection <- suppressWarnings(sf::st_intersection(poly_A, poly_B))

  # Calculate areas
  area_A <- sf::st_area(poly_A) |> sum()
  area_B <- sf::st_area(poly_B) |> sum()
  area_intersection <- sf::st_area(intersection) |> sum()

  if (show_plot) {
    p <- silhouette_plot(poly_A, poly_B) +
      ggplot2::labs(title = 'Shape-to-shape comparison')

    print(p)
  }

  # return the average percentage overlap
  if (length(area_intersection) == 0) {
  #if (any(length(area_intersection) == 0,
  #        length(area_A) == 0,
  #        length(area_A) == 0)) {
      return(0)
  } else {
    mean_overlap <- mean(c(area_intersection / area_A, area_intersection / area_B), na.rm = T)
    return(mean_overlap)
  }

}


##############################
## silhouette_overlap_group
##############################

#' Average overlap between the silhouettes of an individual contour and a group
#' of contours
#'
#' @description
#' Calculates the average percentage overlap between the silhouette of an
#' individual contour and every member of a group of contours.
#'
#' The individual exam (\code{contour_exam}) should be a data frame with the
#' same structure as [IQeyes::sample_contour], containing one row for each point
#' in the contour. The group (\code{poly_group}) should be a data frame
#' containing the polygon geometries of the group.
#'
#' If the contours need to be scaled and/or
#' rotated, the necessary transforms must be performed prior to calling
#' \code{silhouette_overlap_group()}. See [IQeyes:scale_rotate].
#'
#' @param contour_exam
#' A data frame containing the individual contour.
#' @param poly_group
#' A data frame containing a group of polygons. The data frame
#' \emph{must} include a \code{cluster} column that uniquely identifies
#' each member of the group.
#' @param show_plot
#' A Boolean. \code{TRUE} to render a plot that illustrates the comparisons
#' being made. To access the plot as an object, see [IQeyes:silhouette_plot].
#' @param return_detail
#' A Boolean. \code{FALSE} (the default) to return the index of the group member
#' that has the highest percentage of overlap with the reference silhouette.
#' \code{TRUE} to return a list object with two elements: \code{closest_fit} and
#' \code{overlap}, a vector containing the average percentage of
#' overlap between the reference silhouette and every member of the group.
#'
#' @return
#' If \code{return_detail = F}, the function returns just the index of the
#' group member that has the highest percentage of overlap with
#' the reference silhouette.
#'
#' If \code{return_detail = T}, the data frame will be a list object object with
#' two elements: \code{closest_fit}, the index of the group member with the
#' highest percentage of overlap, and \code{overlap}, a vector containing the
#' percentage overlap between the reference silhouette and every member of the
#' group.
#'
#' @examples
#' silhouette_overlap_group(
#'   contour_exam = get_contour(
#'     sample_curvature,
#'     interp = T,
#'     contour_power = 45.5
#'   ) |> scale_rotate(axs = 33.6),
#'   show_plot = T,
#'   return_detail = T
#' )
#'
#' @family polygons
#'
#' @importFrom dplyr group_by
#' @importFrom dplyr group_split
#' @importFrom sf st_intersection
#' @importFrom sf st_as_sf
#' @importFrom sf st_area
#' @importFrom ggplot2 labs
#'
#' @export
silhouette_overlap_group <- function(contour_exam,
                                     poly_group = reference_polygons,
                                     show_plot = F,
                                     return_detail = F) {

  # convert contours to polygons
  #poly_exam <- contour_to_sf_polygon(contour_exam)
  poly_exam <- st_as_sf(contour_exam)

  # calculate the average overlap of the exam with each member of the group
  mean_overlap_group <- poly_group |>
    dplyr::group_by(cluster) |>
    dplyr::group_split() |>
    lapply(function(x) {
      # identify overlap/intersection
      intersection <- suppressWarnings(sf::st_intersection(x, poly_exam))

      # Calculate areas
      area_A <- sf::st_area(x) |> sum()
      area_B <- sf::st_area(poly_exam) |> sum()
      area_intersection <- sf::st_area(intersection) |> sum()

      # return the average percentage overlap
      if (length(area_intersection) == 0) {
        return(0)
      } else {
        mean_overlap <- mean(c(area_intersection / area_A, area_intersection / area_B),
                             na.rm = T)
        return(mean_overlap)
      }
    }) |>
    unlist()

  closest_fit <- which.max(mean_overlap_group)

  # plot
  if (show_plot) {
    p <- silhouette_plot(poly_exam, poly_group) +
      ggplot2::labs(title = 'Shape-to-shape comparisons')

    print(p)
  }

  if (return_detail) {
    return(list(
      closest_fit = closest_fit,
      overlap = mean_overlap_group
    ))
  } else {
    return(closest_fit)
  }

}


#####################
## silhouette_plot
#####################

#' Plot the overlap of one silhouette with one or more silhouettes
#'
#' @description
#' Plots the overlap of one contour's silhouette with either another silhouette
#' or a group of silhouettes.
#'
#' Plotting requires sf polygon objects. These objects can be created by
#' calling [IQeyes::contour_to_sf_polygon] before calling
#' \code{silhouette_plot()}.
#'
#' @param ref_polygon
#' A data frame containing the polygon object for the exam.
#' @param compare_polygons
#' A data frame containing one or more polygons. If the data frame contains a
#' \code{cluster} column, the plot will be faceted by cluster.
#'
#' @details
#' For convenience, the polygons for the reference shapes
#' ([IQeyes:reference_polygons]) are included in this package.
#'
#' @return
#' A ggplot2 object that that shows the overlapping silhouettes.
#'
#' @examples
#' sample_polygon <- get_contour(
#'   sample_curvature,
#'   interp = T,
#'   add_edge = T,
#'   contour_power = 45.5
#' ) |>
#'   scale_rotate(axs = 33.6) |>
#'   contour_to_sf_polygon()
#'
#' silhouette_plot(
#'   ref_polygon = sample_polygon,
#'   compare_polygons = reference_polygons,
#'   highlight_cluster = 1
#' )
#'
#' @family polygons
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_sf
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_blank
#' @importFrom dplyr mutate
#' @importFrom sf st_as_sf
#'
#' @export
silhouette_plot <- function(ref_polygon, compare_polygons, highlight_cluster = NULL) {

  # if the compare_polygons data frame contains a cluster column, facet the plot
  if ('cluster' %in% colnames(compare_polygons)) {
    p <- ggplot2::ggplot() +
      ggplot2::geom_sf(data = ref_polygon, fill = 'deepskyblue', alpha = 0.3, color = 'black') +
      ggplot2::geom_sf(data = compare_polygons, fill = 'darkslategrey', alpha = 0.3, color = 'black') +
      ggplot2::facet_wrap(~ cluster)

  } else {
    p <- ggplot2::ggplot() +
      ggplot2::geom_sf(data = ref_polygon, fill = 'deepskyblue', alpha = 0.3, color = 'black') +
      ggplot2::geom_sf(data = compare_polygons, fill = 'darkslategrey', alpha = 0.3, color = 'black')
  }

  if (!is.null(highlight_cluster)) {
    if(highlight_cluster %in% compare_polygons$cluster) {

      # add a cluster column to the reference sf object and replot it
      dat <- ref_polygon |>
        as.data.frame() |>
        dplyr::mutate(cluster = highlight_cluster) |>
        sf::st_as_sf()

      # blue and pink add up to purple
      p <- p +
        ggplot2::geom_sf(data = dat,
                         fill = 'yellow', alpha = 0.3, color = 'black')
    }
  }

  # add the label and remove the x- and y-axis elements
  p <- p +
    ggplot2::theme(
        axis.title.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        axis.line.x = ggplot2::element_blank(),
        axis.line.y = ggplot2::element_blank()
      )

  return(p)
}
