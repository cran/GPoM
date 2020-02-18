#' @title Generate the conventional order for polynomial terms
#' in a the polynomial formulation
#'
#' @description Generate the conventional order of the polynomial
#' terms for the polynomial description.
#' It is formulated as a matrix of exponents: Each
#' column of the matrix (a,b,c, ...) corresponds to a product
#' of the \code{nVar} available variables X1, X2, X3, etc.,
#' that is, \eqn{X1^a X2^b X3^c}, etc.
#'
#' @param nVar The number of variables
#' @param dMax The maximum degree allowed in the formulation
#'
#' @return A matrix of exponents. Each column corresponds to one
#' polynomial term. Each line correspond to the exponent of one
#' variable.
#' For example, a column of three exponents \code{(0,2,1)} corresponds
#' to the monomial \code{X1^0 * X2^2 * X3^1}, that is \eqn{X2^2 X3}.
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
