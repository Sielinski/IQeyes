#' join_fields
#'
#' A character vector that names of the columns that uniquely identify an exam
#' record. Convenient for dplyr commands (e.g., \code{select(all_of(join_fields)))}).
#'
#' @name join_fields
#' @format A character vector of six (6) elements:
#' \describe{
#' \item{last_name}{Last name of the patient}
#' \item{first_name}{First name of the patient}
#' \item{pat_id}{Patient ID}
#' \item{exam_date}{Date of the exam}
#' \item{exam_time}{Time of the exam}
#' \item{exam_eye}{'RIGHT' or 'LEFT'}
#' }
#'
"join_fields"
