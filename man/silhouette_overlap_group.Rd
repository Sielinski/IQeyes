% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/polygons.R
\name{silhouette_overlap_group}
\alias{silhouette_overlap_group}
\title{Average overlap between the silhouettes of an individual contour and a group
of contours}
\usage{
silhouette_overlap_group(
  contour_exam,
  poly_group = reference_polygons,
  show_plot = F,
  return_detail = F
)
}
\arguments{
\item{contour_exam}{A data frame containing the individual contour.}

\item{poly_group}{A data frame containing a group of polygons. The data frame
\emph{must} include a \code{cluster} column that uniquely identifies
each member of the group.}

\item{show_plot}{A Boolean. \code{TRUE} to render a plot that illustrates the comparisons
being made. To access the plot as an object, see \link{IQeyes:silhouette_plot}.}

\item{return_detail}{A Boolean. \code{FALSE} (the default) to return the index of the group member
that has the highest percentage of overlap with the reference silhouette.
\code{TRUE} to return a list object with two elements: \code{closest_fit} and
\code{overlap}, a vector containing the average percentage of
overlap between the reference silhouette and every member of the group.}
}
\value{
If \code{return_detail = F}, the function returns just the index of the
group member that has the highest percentage of overlap with
the reference silhouette.

If \code{return_detail = T}, the data frame will be a list object object with
two elements: \code{closest_fit}, the index of the group member with the
highest percentage of overlap, and \code{overlap}, a vector containing the
percentage overlap between the reference silhouette and every member of the
group.
}
\description{
Calculates the average percentage overlap between the silhouette of an
individual contour and every member of a group of contours.

The individual exam (\code{contour_exam}) should be a data frame with the
same structure as \link{sample_contour}, containing one row for each point
in the contour. The group (\code{poly_group}) should be a data frame
containing the polygon geometries of the group.

If the contours need to be scaled and/or
rotated, the necessary transforms must be performed prior to calling
\code{silhouette_overlap_group()}. See \link{IQeyes:scale_rotate}.
}
\examples{
silhouette_overlap_group(
  contour_exam = get_contour(
    sample_curvature,
    interp = T,
    contour_power = 45.5
  ) |> scale_rotate(axs = 33.6),
  show_plot = T,
  return_detail = T
)

}
\seealso{
Other polygons: 
\code{\link{contour_to_sf_polygon}()},
\code{\link{silhouette_overlap}()},
\code{\link{silhouette_plot}()}
}
\concept{polygons}
