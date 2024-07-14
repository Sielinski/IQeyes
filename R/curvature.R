####################
## true_net_power
####################

#' True net power
#' @description
#' Calculates the true net power of a cornea based on the front (anterior)
#' and rear (posterior) refractive indices.
#' @param R_front
#' Front refractive index.
#' @param R_rear
#' Front refractive index.
#' @return
#' Returns the optical power of the cornea.
#' @references
#' \emph{Pentacam Interpretation Guide}, p 9.
#' @examples
#' true_net_power(7, 5)
#'
#' @family Curvature
#'
#' @export
true_net_power <- function(R_front, R_back){
  tnp <- ((376 / R_front) - (040 / R_back)) |>
    round(1)
  return(tnp)
}


####################
## anterior_power
####################

#' Anterior power
#' @description
#' Calculates the dioptric power of the anterior (front) surface of a cornea
#' based on the radius of curvature.
#' @param R
#' Front refractive index.
#' @return
#' Returns the optical power of the anterior surface of a cornea.
#' @references
#' \emph{Pentacam Interpretation Guide}, pp 6-7.
#' @examples
#' anterior_power(7)
#'
#' @family Curvature
#'
#' @export
anterior_power <- function(R) 337.5 / R


#####################
## posterior_power
#####################

#' Posterior power
#' @description
#' Calculates the dioptric power of the posterior (rear) surface of a cornea
#' based on the radius of curvature.
#' @param R
#' Rear refractive index.
#' @return
#' Returns the optical power of the posterior surface of a cornea.
#' @references
#' \emph{Pentacam Interpretation Guide}, p 7.
#' @examples
#' posterior_power(5)
#'
#' @family Curvature
#'
#' @export
posterior_power <- function(R) 30 / R  # 1.376/1.336 ≈ 0.30


#####################
## kappa_from_radius
#####################

#' Kappa from the radius of curvature
#' @description
#' Calculate curvature (kappa, K) from the radius of curvature (R).
#' @param R
#' Radius of curvature.
#' @return
#' Returns the curvature.
#' @references
#' \emph{Pentacam Interpretation Guide}, p 7.
#' @details
#' Curvature (K) is different from the conic constant (Q)
#' @examples
#' kappa_from_radius(7)
#'
#' @family Curvature
#'
#' @export
kappa_from_radius <- function(R) 1/R


