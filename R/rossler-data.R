#' This is data to be included in my package
#'
#' @name Rossler-1976 data set
#' @docType data
#' 
#' @description
#' The Rossler system is the 3-dimensional chaotic system \cr
#' \eqn{dx/dt = - y - z} \cr
#' \eqn{dy/dt = x + a y} \cr
#' \eqn{dz/dt = b + z (x - c)}, \cr
#' discovered by Otto Rossler in 1976 [1].
#' The following parameters an initial conditions
#' were used in the present data set:\cr
#' a = 0.520, b = 2, c = 4 \cr
#' and (x0, y0, z0) = (-0.04298734, 1.025536, 0.09057987).\cr
#' The following four columns are provided:\cr
#' (1) time t, (2) x(t), (3) y(t) and (4) z(t).\cr
#' For this parameterization, the Rossler system produces
#' a chaotic behavior non-coherent in phase.
#' 
#' @author Sylvain Mangiarotti, Flavie Le Jean, 
#' Malika Chassan, Laurent Drapeau, Mireille Huc.
#' 
#' @references
#' [1] P. Rossler, 1976.
#' An Equation for Continuous Chaos, 
#' Physics Letters, 71A, 2-3, 155-157. \cr
#' 
#' @keywords data
"rossler"
