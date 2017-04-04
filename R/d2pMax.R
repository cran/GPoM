#' @title d2pMax : provides pMax given dMax and nVar
#'
#' @description Gives the maximum polynomial degree pMax
#'              given the number of polynomial terms dMax
#'              and the system dimension nVar
#'
#' @inheritParams gloMoId
#' @param dMaxKnown The number of polynoms terms to retrieve
#'
#' @return The maximum polynomial degree \code{dMax} used to code the polynomial
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
d2pMax <- function(nVar, dMaxKnown) {
  
  pM <- choose(nVar + dMaxKnown, nVar)
  
  pM
}
