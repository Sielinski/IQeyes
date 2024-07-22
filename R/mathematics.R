#####################
## signed_pct_diff
#####################

#' Signed percentage difference
#' @description
#' Calculates the signed percentage difference between two numbers.
#' @param x
#' A number. If \code{x} is 0, \code{y} needs to be non-zero.
#' @param y
#' A number. If \code{y} is 0, \code{x} needs to be non-zero.
#' @return
#' Returns \code{(x - y) / (max(abs(x), abs(y)))}. Note
#' that the result is directional: If \code{x} > \code{y}, then
#' \code{signed_pct_diff} returns a positive value; If \code{x} < \code{y}, then
#' \code{signed_pct_diff} returns a negative value.
#' @examples
#' signed_pct_diff(x = 5, y = 10)
#' signed_pct_diff(5, 10)
#'
#' @export
signed_pct_diff <- function(x, y) {
  if (x == 0 & y == 0) {
    warning('Either x or y needs to be non-zero.')
  } else {
    (x - y) / (max(abs(x), abs(y)))
  }
}

signed_pct_diff <- Vectorize(signed_pct_diff)


###########
## logit
###########

#' Compute the logit transformation of proportions or percentages
#' @description
#' Maps a probability value from \[0, 1\] to a real number in \[−Inf, +Inf\].
#' Useful for mapping a probability to a real number.
#' @param p
#' A numeric probability ranging from 0 to 1.
#' @return
#' Returns \code{log(p / (1 - p))}.
#' @seealso
#' [IQeyes::logistic()]
#' @details
#' Inverse of the logistic function.
#' @examples
#' logit(p = 0.5)
#' logit(1/2)
#'
#' @family mathematics
#'
#' @export
logit <- function(p) {
  if (p < 0 | p > 1) {
    warning('p needs to be between 0 and 1.')
  } else {
    log(p / (1 - p))
  }
}

logit <- Vectorize(logit)


##############
## logistic
##############

#' Compute the logistic transformation of number
#' @description
#' Maps a numeric value from \[−Inf, +Inf\] to \[0, 1\]. Useful for mapping real
#' numbers to a probability.
#' @param x
#' A numeric value.
#' @return
#' Returns \code{1 / (1 + exp(-x))}.
#' @seealso
#' [IQeyes::logit()]
#' @details
#' Inverse of the logit function.
#' @examples
#' logistic(x = 64)
#' logistic(Inf)
#'
#' @family mathematics
#'
#' @export
logistic <- function(x) {
  1 / (1 + exp(-x))
}


#############
## nearest
#############

#' Round a number at the specified grain
#' @description
#' Rounds \code{x} to the nearest multiple of \code{y}.
#' @param x
#' A numeric vector of values to be rounded.
#' @param y
#' A number. The grain for rounding.
#' @return
#' The rounded number(s).
#' @examples
#' nearest(0.75, 0.2)
#'
#' @family mathematics
#'
#' @export
nearest <- function(x, y) round(x / y) * y


####################
## array_position
####################

#' Returns the index/position of a 2d matrix for a given row and column
#' @description
#' This function does the inverse of \code{arrayInd(index, dim(array))}.
#'
#' @param matrix
#' A two-dimensional array.
#' @param row_index
#' The row number.
#' @param col_index
#' The column number.
#' @param row_label
#' A character value matching one of the row names of \code{matrix}.
#' @param col_label
#' A character value matching one of the column names of \code{matrix}.
#'
#' @return
#' Returns a numeric index (of the matrix).
#'
#' @details
#' Requires both a row and column, which can be provided as a \code{row_index}
#' and \code{col_index}, as a \code{row_label} and \code{col_label}, or some
#' combination of an index and a label.
#'
#' @examples
#' # find the index of row 2, column 1
#' matrix(seq(1:9), nrow = 3, ncol = 3, byrow = T) |>
#'   array_position(row_index = 2, col_index = 1)
#'
#' # find the row and column of index 2
#' arrayInd(2, dim(matrix(seq(1:9), nrow = 3, ncol = 3, byrow = T)))
#'
#' @family mathematics
#'
#' @export
array_position <- function(matrix, row_index = NULL, col_index = NULL, row_label = NULL, col_label = NULL) {

  # Get the row and column indices
  if (is.null(row_index) & !is.null(row_label)) row_index <- which(rownames(matrix) == row_label)
  if (is.null(col_index) & !is.null(col_label)) col_index <- which(colnames(matrix) == col_label)

  if (is.null(row_index) | is.null(col_index)) warning('Must provide both a row and a column.')

  # Calculate the linear index
  position <- (col_index - 1) * nrow(matrix) + row_index
  return(position)
}
