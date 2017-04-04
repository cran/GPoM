#' @title p2dMax : provides dMax given pMax and nVar
#'
#' @description Search the maximum polynomial degree dMax
#'              given the number of polynomial terms pMax
#'              and the system dimension nVar
#'
#'
#' @inheritParams gloMoId
#' @param pMaxKnown  The number of polynoms terms to retrieve
#'
#' @return The maximum polynomial degree \code{dMax} used to code the polynomial
#'
#' @author Sylvain Mangiarotti, Laurent Drapeau
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
#' # size of the polynomial vector:
#' pMax <- 10
#' # The maximal polynomial degree used for coding the polynomial is:
#' p2dMax(nVar,pMax)
#'
#' #############
#' # Example 2 #
#' #############
#' # for pMax = 462 and nVar = 6, then dMax is:
#' p2dMax(6,462)
#' # indeed:
#' length(poLabs(nVar=6, dMax=5))
#'
#' @export
#'
p2dMax <- function(nVar, pMaxKnown) {

  # initialisation
  pM <- dM <- 0
  nIter = 500

  # iterative searching loop
  while (pM != pMaxKnown & dM != nIter) {
    dM <- dM + 1
    pM <- d2pMax(nVar, dM)
  }

  # stop if pM was not found
  if (dM == nIter) {
    block <- paste("Impossible to find 'dMax' from 'nVar' and 'pMax'
                   with ", nIter, "iterations.", sep="")
    stop(block)
  }

  # retreived maximum polynomial degree:
  dM
}
