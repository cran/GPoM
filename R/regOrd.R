#' @title regOrd : Defines the conventional order for regressors
#' in the polynomial formulation
#'
#' @description Defines the order convention for regressors in
#' polynomial. It is formulated as a matrix of explonents. Each
#' column of the matrix (a,b,c, ...) corresponds to one regressor 
#' defined as the product of X1^a * X2^b * X3^c, etc.
#'
#' @param nVar The number of variables
#' @param dMax The maximum degree allowed for the polynomial
#'
#' @return A matrix of exponent. Each column correspond to one
#' polynomial term. Each line correspond to one variable.
#' For example, a column of three exponents (0,2,1) correspond
#' to regressor X1^0 * X2^2 * X3^1, that is X2^2 * X3.
#'
#' @seealso \code{\link{poLabs}}
#'
#' @author Sylvain Mangiarotti
#'
#' @export
regOrd <- function(nVar,dMax) {
  # total number of regressors
  pMax <- d2pMax(nVar, dMax)
  # prepare the matrix of exponents
  #
  pExpo <- as.matrix(0:dMax)
  for (i in 1:(nVar-1)) {
    pExpotmp <- pExpo
    tmpEx <- pExpotmp[,1] * 0
    pExpo <- tmpEx2 <- cbind(pExpotmp, tmpEx)
    for (j in 1:dMax) {
      tmpEx <- pExpotmp[,1] * 0 + j
      tmpEx2 <- cbind(pExpotmp, tmpEx)
      pExpo <- rbind(pExpo, tmpEx2)
    }
  }
  # select regressors such as degree <= dMax
  ltdMax <- colSums(t(pExpo)) <= dMax
  pExpo <- t(pExpo[ltdMax,])[nVar:1,]
  # rownames
  noms <- c("X1")
  for (i in 2:nVar) noms <- c(noms, paste("X", i, sep=""))
  rownames(pExpo) <- noms
  # return
  pExpo
}
