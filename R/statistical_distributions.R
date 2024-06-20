####################
## sd_half_normal
####################

#' Standard deviation of a half-normal distribution
#' @description
#' Calculates the standard deviation of a sample from a half-normal
#' distribution.
#' @param x
#' A numeric vector of values from a half-normal distribution.
#' @return
#' Returns the standard deviation of \code{x}.
#' @details
#' \code{sd_half_normal} adjusts \code{sd(x)} to account for the characteristics
#' of the half-normal distribution
#' @examples
#' sd_half_normal(abs(rnorm(10)))
#'
#' @family Statistical Functions
#'
#' @export
sd_half_normal <- function(x, ...) sd(x, ...) * sqrt(1 - 2 / pi)


######################
## mean_half_normal
######################

#' Mean of a half-normal distribution
#' @description
#' Calculates the mean (expected value) of a sample from a half-normal
#' distribution.
#' @param x
#' A numeric vector of values from a half-normal distribution.
#' @return
#' Returns the mean of \code{x}.
#' @details
#' \code{mean_half_normal()} uses \code{sd(x)} to account for the characteristics
#' of the half-normal distribution
#' @examples
#' mean_half_normal(abs(rnorm(10)))
#'
#' @family Statistical Functions
#'
#' @export
mean_half_normal <- function(x, ...) sd(x, ...) * sqrt(2 / pi)


########################
## mean_folded_normal
########################

#' Mean of a folded-normal distribution
#' @description
#' Calculates the mean (expected value) of a folded-normal distribution based on
#' the parameters of a normal distribution. Presumes that the folded-normal
#' distribution is the absolute value of the normal distribution.
#' @param mu
#' Mean of the normal distribution.
#' @param sigma
#' Standard deviation of the normal distribution.
#' @return
#' Returns the mean.
#' @examples
#' mean_folded_normal(mu = 3, sigma = 1)
#'
#' @family Statistical Functions
#'
#' @export
mean_folded_normal <- function(mu, sigma) {
  term1 <- sigma * sqrt(2 / pi) * exp(-mu ^ 2 / (2 * sigma ^ 2))
  term2 <- mu * (1 - 2 * pnorm(-mu / sigma))
  return(term1 + term2)
}


######################
## sd_folded_normal
######################

#' Standard deviation of a folded-normal distribution
#' @description
#' Calculates the standard deviation of a folded-normal distribution based on
#' the parameters of a normal distribution. Presumes that the folded-normal
#' distribution is the absolute value of the normal distribution.
#' @param mu
#' Mean of the normal distribution.
#' @param sigma
#' Standard deviation of the normal distribution.
#' @return
#' Returns the standard deviation.
#' @examples
#' sd_folded_normal(mu = 3, sigma = 1)
#'
#' @family Statistical Functions
#'
#' @export
sd_folded_normal <- function(mu, sigma) {
  sqrt(mu ^ 2 + sigma ^ 2 - mean_folded_normal(mu, sigma) ^ 2)
}


##############
## z_to_pct
##############

#' Z-score to a percentage (of a normal distribution)
#' @description
#' Converts a z-score to the percentage of a normal distribution that
#' lie within z standard deviations of the mean.
#' @param z
#' Mean of a normal distribution.
#' @param two_tail
#' A Boolean indicating whether the z-score is based on two (\code{T}) or
#' one (\code{F}) tails. Basically, two_tail = \code{T} cuts off values outside
#' \code{±z * σ}.
#' @return
#' Returns the percentage of a normal distribution.
#' @examples
#' z_to_pct(z = 3, two_tail = T)
#' z_to_pct(c(1, 2, 3), two_tail = F)
#'
#' @family Statistical Functions
#'
#' @export
z_to_pct <- function(z, two_tail = T) {
  # use the two_tail parameter to calculate the percentage of values that lie
  # within z standard deviations from the mean of a normal distribution (in both
  # directions).
  if (two_tail) {
    pnorm(z) - pnorm(-z)
  } else {
    pnorm(z)
  }
}


##############
## pct_to_z
##############

#' Z-score from a percentage (of a normal distribution)
#' @description
#' Calculates a z-score based on the percentage of values from a normal
#' distribution.
#' @param pct
#' Percentage of values.
#' @param two_tail
#' A Boolean indicating whether the z-score should be based on two (\code{T}) or
#' one (\code{F}) tails. Basically, two_tail = \code{T} cuts off values outside
#' \code{±z * σ}.
#' @return
#' Returns the z-score for the percentage of values.
#' @examples
#' pct_to_z(pct = 0.5, two_tail = F)
#' pct_to_z(c(-0.5, 0.5), two_tail = T)
#'
#' @family Statistical Functions
#'
#' @export
pct_to_z <- function(pct, two_tail = T) {
  if (two_tail) {
    upper_tail <- pct + ((1 - pct) / 2)
    qnorm(upper_tail)
  } else {
    qnorm(pct)
  }
}


##############
## z_score
##############

#' Z-score
#' @description
#' Calculates the z-score for a value \code{x} taken from a normal distribution
#' defined by it mean and standard deviation.
#' @param x
#' A value from a normal distribution.
#' @param mean
#' The mean of the normal distribution
#' @param sd
#' The standard deviation of the normal distribution
#' @return
#' Returns the z-score of \code{x}.
#' @examples
#' z_score(0.5, mean = 0, sd = 1)
#'
#' @family Statistical Functions
#'
#' @export
z_score <- function(x, mean, sd) {
  (x - mean) / sd
}


############
## d_mode
############

#' Dominant mode of a sample
#' @description
#' Calculates the mode for a vector of numeric values \code{x}.
#' @param x
#' A vector of numeric values.
#' @return
#' Returns the mode of \code{x}.
#' @examples
#' d_mode(rnorm(1000, mean = 2))
#'
#' @family Statistical Functions
#'
#' @export
d_mode <- function(x) {
  if (length(x) == 1) {
    return(x)
  } else {
    den <- density(x, kernel = 'gaussian')
    return(den$x[which.max(den$y)])
  }
}


##########
## knee
##########

#' Find the knee in a scree plot
#' @description
#' Finds the knee in the curve of a two-dimensional scree plot.
#' @param df
#' A data frame that contains the \emph{x} and \emph{y} coordinates of the points
#' that comprise a scree plot, ordered from small to large values of \emph{x}.
#' The columns must be labeled \code{x} and \code{y}.
#' @return
#' Returns the row of the data frame that contains the knee of the curve.
#' @examples
#' data.frame(x = 1:15, y = sort(rpois(15, 42), decreasing = T)) |>
#'   knee()
#'
#' @family Statistical Functions
#'
#' @export
knee <- function(df) {
  data_pts <- nrow(df)

  if (data_pts <= 3) {
    knee <- df[max(1, round(data_pts / 2, 0)), 1]

  } else {
    # identify the start and end points of the line
    pt_start <- df[1, ]
    pt_end <- df[data_pts, ]

    # calculate slope and intercept
    slope <- (pt_end$y - pt_start$y) / (pt_end$x - pt_start$x)
    intercept <- pt_start$y - slope * pt_start$x

    # plot the function and the sloping line
    #ggplot(df, aes(x = x, y = y)) +
    #  geom_point() +
    #  geom_line() +
    #  geom_abline(slope = slope, intercept = intercept, color = 'red')

    # create a vector to hold the lengths of the perpendicular lines
    perp_length <- rep(NA, data_pts)

    # for the j for-loop: The knee will depend upon the number of points
    # along the x-axis
    x_points <- rep(NA, data_pts)

    perp_slope <- -1 / slope

    for (i in 2:data_pts) {
      # calculate slope and intercept of perpendicular line
      perp_inter <- df[i, 2] - (perp_slope * df[i, 1])  # b = y - mx

      # at the intersection of the perpendicular and sloping lines,
      # both lines have the same y value, so solve (mx + b)[1] = (mx + b)[2] to find x
      x_intersection <-
        (perp_inter - intercept) / (slope - perp_slope)
      # then find y
      y_intersection <- slope * x_intersection + intercept

      # length will equal sqrt(a^2 + b^2)
      perp_length[i] <-
        sqrt((x_intersection - df[i, 1]) ^ 2 + (y_intersection - df[i, 2]) ^ 2)
    }

    #ggplot(df, aes(x = x, y = y)) +
    #  geom_point() +
    #  geom_line() +
    #  geom_abline(slope = slope, intercept = intercept, color = 'red') +
    #  geom_abline(slope = slope, intercept = max_intercept, color = 'blue') +
    #  geom_point(data = df[which.max(perp_length), ], aes(x = x, y = y), color = 'yellow')

    # which is the max?
    knee <- df[which.max(perp_length), 1]
  }

  return(knee)
}
