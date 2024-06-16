#' canonical_polygons
#'
#' Silhouettes (filled shapes) of the canonical contours.
#'
#' @name canonical_polygons
#' @format A data frame with the following columns:
#' \describe{
#' \item{last_name}{Last name of the patient}
#' \item{first_name}{First name of the patient}
#' \item{pat_id}{Patient ID}
#' \item{exam_date}{Date of the exam}
#' \item{exam_time}{Time of the exam}
#' \item{exam_eye}{'RIGHT' or 'LEFT'}
#' \item{surface}{'FRONT' or 'BACK' surface of the cornea}
#' \item{x}{x-axis coordinate of the measurement}
#' \item{y}{y-axis coordinate of the measurement}
#' \item{contour_id}{Ordinal identifier of the segment(s) that comprise a contour}
#' \item{contour}{Dioptric power of the contour}
#' \item{cluster}{Oridinal identifier of the cluster}
#' \item{shape}{Name of the cluster}
#' }
#'
"canonical_polygons"
