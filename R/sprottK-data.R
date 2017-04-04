#' Data included in package GPoM
#'
#' @name sprottK time series
#' @docType data
#' 
#' @description
#' The Sprott-K system is the 3-dimensional chaotic system \cr
#' \eqn{dx/dt = x z - a y} \cr
#' \eqn{dy/dt = x + b y} \cr
#' \eqn{dz/dt = x - z}, \cr
#' discovered by Julien C. Sprott in 1994 [1].
#' The following parameters and initial conditions
#' were used in the present data set:\cr
#' a = 0.5, b = 0.2 \cr
#' and (x0, y0, z0) = (-0.01170987, 1.010009, -4.522642e-05).\cr
#' The following four columns are provided:\cr
#' (1) time t, (2) x(t), (3) y(t) and (4) z(t).\cr
#' For this parameterization, the Sprott-K system
#' produces a periodic behavior of period 4.
#' 
#' @author Sylvain Mangiarotti, Flavie Le Jean.
#' 
#' @references
#' [1] J.C. Sprott, Some simple flows, 1994.
#' Physical Review E, 50(2), 647-650.
#' 
#' @keywords data
"sprottK"
