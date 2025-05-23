% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mathematics.R
\name{nearest}
\alias{nearest}
\title{Round a number at the specified grain}
\usage{
nearest(x, y)
}
\arguments{
\item{x}{A numeric vector of values to be rounded.}

\item{y}{A number. The grain for rounding.}
}
\value{
The rounded number(s).
}
\description{
Rounds \code{x} to the nearest multiple of \code{y}.
}
\examples{
nearest(0.75, 0.2)

}
\seealso{
Other mathematics: 
\code{\link{array_position}()},
\code{\link{logistic}()},
\code{\link{logit}()}
}
\concept{mathematics}
