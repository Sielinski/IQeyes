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
#' \code{contour_polygon}.
#'
#' @return
#' An sf object.
#'
#' @details
#' If the contour contains segments that are fully enclosed in another segment,
#' this function will \emph{subtract} those segments from the resulting polygon.
#'
#' @examples
#' contour_to_sf_polygon(
#'   contour = get_contour(sample_curvature, contour_power = 44.5)
#' )
#'
#' @family Polygons
#'
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
#' @importFrom purrr reduce
#'
#' @export
contour_to_sf_polygon <- function(contour, required_points = 4) {

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

  # Convert the data frame to an sf object with polygons
  sf_segments <- df |>
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
  combined_polygon <- purrr::reduce(final_polygons, sf::st_union)

  # repeat the st_contains test, but now subtract any inner polygons from
  # the combined polygon
  for (i in 1:nrow(sf_segments)) {
    current_polygon <- sf_segments[i, ]

    for (j in 1:nrow(sf_segments)) {
      if (i != j) {
        test_polygon <- sf_segments[j, ]
        if (sf::st_contains(current_polygon, test_polygon, sparse = FALSE)) {
          combined_polygon <- suppressWarnings(sf::st_difference(combined_polygon, test_polygon))
        }
      }
    }
  }

  # if necessary, repair the polygon
  if (!sf::st_is_valid(combined_polygon))
    combined_polygon <- sf::st_make_valid(combined_polygon)

  return(combined_polygon)
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
#' @param contour_A
#' A data frame containing the contour A.
#' @param contour_B
#' A data frame containing the contour B.
#' @param show_plot
#' A Boolean. \code{TRUE} to render a plot that illustrates the comparison
#' being made. To access the plot as an object, see [IQeyes:silhouette_plot].
#'
#' @return
#' The average percentage overlap between the two silhouettes.
#'
#' @examples
#' silhouette_overlap(
#'   contour_A = get_contour(sample_curvature, contour_power = 44),
#'   contour_B = get_contour(sample_curvature, contour_power = 44.5),
#'   show_plot = T
#' )
#'
#' @family Polygons
#'
#' @importFrom sf st_intersection
#' @importFrom sf st_area
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_sf
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
  area_A <- sf::st_area(poly_A)
  area_B <- sf::st_area(poly_B)
  area_intersection <- sf::st_area(intersection)

  if (show_plot) {
    #p <- ggplot2::ggplot() +
    #  ggplot2::geom_sf(data = poly_A, fill = 'coral', alpha = 0.3, color = 'black') +
    #  ggplot2::geom_sf(data = poly_B, fill = 'deepskyblue', alpha = 0.3, color = 'black') +
    #  ggplot2::labs(title = 'Overlapping silhouettes')

    p <- silhouette_plot(poly_A, poly_B)

    print(p)
  }

  # return the average percentage overlap
  if (length(area_intersection) == 0) {
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
#' @param contour_exam
#' A data frame containing the individual contour.
#' @param contour_group
#' A data frame containing the group of contours. The data frame
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
#'   contour_exam = get_contour(sample_curvature, contour_power = 45.5),
#'   contour_group = canonical_contours,
#'   show_plot = T,
#'   return_detail = T
#' )
#'
#' @family Polygons
#'
#' @export
silhouette_overlap_group <- function(contour_exam, contour_group, show_plot = F, return_detail = F) {

  # convert contours to polygons
  poly_exam <- contour_to_sf_polygon(contour_exam)

  poly_group <- contour_group |>
    dplyr::group_by(dplyr::across(tidyselect::all_of(join_fields))) |>
    dplyr::group_split() |>
    lapply(contour_to_sf_polygon)

  # calculate the average overlap of the exam with each member of the group
  mean_overlap_group <- lapply(poly_group, function(x) {
    # identify overlap/intersection
    intersection <- suppressWarnings(sf::st_intersection(x, poly_exam))

    # Calculate areas
    area_A <- sf::st_area(x)
    area_B <- sf::st_area(poly_exam)
    area_intersection <- sf::st_area(intersection)

    # return the average percentage overlap
    if (length(area_intersection) == 0) {
      return(0)
    } else {
      mean_overlap <- mean(c(area_intersection / area_A, area_intersection / area_B), na.rm = T)
      return(mean_overlap)
    }
  })

  # plot
  if (show_plot) {
    plot_group <- dplyr::bind_rows(poly_group, .id = 'cluster') |>
      as.data.frame() |>
      dplyr::select(cluster, geometry) |>
      sf::st_as_sf()

    p <- silhouette_plot(poly_exam, plot_group)

    #p <- ggplot2::ggplot(data = plot_group) +
    #  ggplot2::geom_sf(fill = 'deepskyblue', alpha = 0.3, color = 'black') +
    #  ggplot2::geom_sf(data = poly_exam, fill = 'coral', alpha = 0.3, color = 'black') +
    #  ggplot2::facet_wrap(~ cluster) +
    #  ggplot2::labs(title = 'Overlapping silhouettes')

    print(p)
  }

  mean_overlap_group <- unlist(mean_overlap_group)

  closest_fit <- which.max(mean_overlap_group)

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
#' For convenience, the polygons for the canonical shapes
#' ([IQeyes:canonical_polygons]) are included in this package.
#'
#' @return
#' A ggplot2 object that that shows the overlapping silhouettes.
#'
#' @examples
#' sample_polygon <- get_contour(sample_curvature, add_edge = T, contour_power = 45.5) |>
#'   scale_rotate(axs = 33.6) |>
#'   contour_to_sf_polygon()
#'
#' silhouette_plot(
#'   ref_polygon = sample_polygon,
#'   compare_polygons = canonical_polygons,
#'   highlight_cluster = 1
#' )
#'
#' @family Polygons
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_polygon
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_blank
#' @importFrom dplyr mutate
#' @importFrom sf st_as_sf
#'
#' @export
silhouette_plot <- function(ref_polygon, compare_polygons, highlight_cluster = NULL) {

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
    ggplot2::labs(title = 'Overlapping silhouettes') +
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
