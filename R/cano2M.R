#' @title cano2M : Converts a model in canonical form into a matrix form
#'
#' @description Converts the vectorial formulation of canonical models
#' into a matrix formulation (that is, including explicitely all the
#' equations). For both input, the list of terms follows the convention
#' defined by \code{poLabs}.
#'
#' @inheritParams poLabs
#' @inheritParams derivODE2
#' @inheritParams gloMoId
#'
#' @param poly A vector of coefficients corresponding to the regressor
#' of the canonical function
#'
#' @author Sylvain Mangiarotti, Mireille Huc
#'
#' @seealso \code{\link{drvSucc}}, \code{\link{gPoMo}}, \code{\link{poLabs}}
#'
#' @examples
#' # A vector of polynomial terms corresponding to a canonical form:
#' polyTerms <- c(0.2,0,-1,0.5,0,0,0,0,0,0)
#' # Convert this vector into a matrix formulation with all the equations: 
#' K <- cano2M(3,2,polyTerms)
#' # Visualize the equations:
#' visuEq(K,3,2)
#' 
#' @export
#' 
#' @return \code{Kmod} A matrix with nVar columns of the complete description of the equations.
#' The first columns relates to the canonical part dX1/dt = X2, dX2/dt = X3 etc. and
#' the column is the polynomial term itself 
cano2M <- function(nVar, dMax, poly, dMin = 0) {
  # check
  pMax <- length(poly)
  if (pMax != d2pMax(nVar, dMax, dMin = dMin))
    stop('Polynomial function size incompatible with formulation poLabs(nVar, dMax).')
  # initiate the output
  Kmod <- matrix(0, ncol = nVar, nrow = pMax)
  # (1) put the polynomial function in the last equation (i = nVar)
  Kmod[, nVar] <- poly
  # (2) search of the regressor X[i+1] such as dX[i]/dt = X[i+1]
  # and put it in the ith equations
  for (i in 1:(nVar - 1)) {
    ireg <- d2pMax(nVar - i - 1, dMax) + 1
    Kmod[ireg, i] <- 1
  }
  # model output
  Kmod
}
