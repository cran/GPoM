#' @title Gram-Schmidt procedure
#'
#' @description Computes regressors coefficients
#' using the Gram-Schmidt procedure.
#'
#' @param polyK One list including \code{$Y} and \code{$phy} with:
#' \code{$Y} a matrix for which the ith column will be used
#' to add one orthogonal vector to the (i-1)th vectors of the
#' current orthogonal base;
#' and \code{$phy} such as the current orthogonal base is
#' given by the (i-1)th first columns of matrix \code{polyK$phy}.
#' @param ivec Defines i, the current vector of \code{polyK$Y} and
#' the current orthogonal base of \code{pParam$phy}.
#' @param weight The weighing vector.
#'
#' @return \code{uNew} The model parameterization, that is:
#' The residual orthogonal vector that can be included into
#' the current orthogonal base. If the current base is empty,
#' \code{uNew} is equal to the input vector of \code{$Y};
#' if the base is complete, \code{uNew} equals 0.
#'
#' @author Sylvain Mangiarotti
#'
GSproc <- function(polyK, ivec, weight = NULL) {
  # initiate uNew and keep a memory in u
  uNew <- u <- polyK$Y[, ivec]
  if (ivec > 1) {
    # Gram-Schmidt iterations
    for (i in 1:(ivec - 1)) {
      v <- polyK$phi[,i]
      # compute the norm
      normV <- wInProd(v, v, weight)
      # compute the projection
      proj <- wInProd(u, v, weight) / normV * v
      # compute the new vector
      uNew <- uNew - proj
    }
  }
  # output vector
  uNew
}
