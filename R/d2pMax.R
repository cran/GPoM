#' @title Provides the number of polynomial terms \code{pMax}
#' given \code{dMax} and \code{nVar}
#'
#' @description Computes the number of polynomial terms \code{pMax}
#'              used to formulate an equation given
#'              the maximal polynomial degree \code{dMax}
#'              and the number of variables \code{nVar}
#'              following the conventions as defined by fuction \code{poLabs}.
#'
#' @inheritParams gloMoId
#' @inheritParams regOrd
#' @param dMaxKnown The maximum polynomial degree \code{dMax}
#'
#' @return The number \code{pMax} of polynomial terms used to code
#' a polynomial equation
#'
#' @author Sylvain Mangiarotti
#'
#' @seealso \code{\link{gloMoId}}, \code{\link{gPoMo}}, \code{\link{poLabs}}
#'
#' @examples
#' #############
#' # Example 1 #
#' #############
#' # Maximum polynomial degree ?
#' # number of variables:
#' nVar <- 3
#' # polynomial degree:
#' dMax <- 3
#' # The maximal polynomial degree used for coding the polynomial is:
#' d2pMax(nVar,dMax)
#'
#' @export
#'
d2pMax <- function(nVar, dMaxKnown, dMin = 0) {

    if (dMin == 0) {
      pM <- choose(nVar + dMaxKnown, nVar)}
    else {
      pM = dim(regOrd(nVar, dMaxKnown, dMin=dMin))[2]
    }
  
    pM
}
