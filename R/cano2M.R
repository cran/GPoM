#' @title canoToM : Converts canonical formulation into matricial formulation 
#'
#' @description Converts the vectorial formulation of canonical models
#' into a matricial formulation. For both input, the list of terms 
#' follows the convention given in `poLabs`.
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
#' @examples
#' polyTerms <- c(0.2,0,-1,0.5,0,0,0,0,0,0)
#' K <- cano2M(3,2,polyTerms)
#' visuEq(3,2,K)
#' 
#' @export
cano2M <- function(nVar, dMax, poly) {
  # check
  pMax <- length(poly)
  if (pMax != d2pMax(nVar, dMax))
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
