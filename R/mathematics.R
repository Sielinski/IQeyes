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
signed_pct_diff <- function(x, y) {
  if (x == 0 & y == 0) {
    warning('Either x or y needs to be non-zero.')
  } else {
    (x - y) / (max(abs(x), abs(y)))
  }
}


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
logit <- function(p) {
  if (p < 0 | p > 1) {
    warning('p needs to be between 0 and 1.')
  } else {
    log(p / (1 - p))
  }
}


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
nearest <- function(x, y) round(x / y) * y
