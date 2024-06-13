#' canonical_curvature
#'
#' Curvature data for the canonical shapes (from COR-PWR).
#'
#' @name canonical_curvature
#' @format A data frame with the following columns:
#' \describe{
#' \item{last_name}{Last name of the patient}
#' \item{first_name}{First name of the patient}
#' \item{pat_id}{Patient ID}
#' \item{exam_date}{Date of the exam}
#' \item{exam_time}{Time of the exam}
#' \item{exam_eye}{'RIGHT' or 'LEFT'}
#' \item{surface}{'FRONT' or 'BACK' surface of the cornea}
#' \item{column_name}{Name of the column in COR-PWR that contained the measurement}
#' \item{measurement}{Measured curvature radius}
#' \item{ring_diam}{Ring where the measurement was taken}
#' \item{angle}{Angle where the measurement was taken}
#' \item{x}{x-axis coordinate of the measurement}
#' \item{y}{y-axis coordinate of the measurement}
#' }
#'
"canonical_curvature"
