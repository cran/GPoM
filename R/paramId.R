#' @title For parameter Identification
#'
#' @description Estimate the polynomial coefficients.
#'
#' @param allForK The list of input parameters
#' @param drv The derivative (on the equation left hand)
#' @param weight The weighting series
#'
#' @importFrom stats lsfit
#' @return \code{allForK} The initial list completed with the model parameters.
#'
#' @author Sylvain Mangiarotti
#'
paramId <- function(allForK, drv, weight) {

  nRegModif <- dim(allForK$A)[1]
  for (i in 1:nRegModif) {
    allForK$phi[,i] <- GSproc(allForK, i, weight)
    allForK$A[i,1:i] <- lsfit(allForK$Y[,1:i],
                               allForK$phi[,i],
                               intercept = F)$coefficients
  }
  #
  cc <- vector("numeric", length = nRegModif)
  for (i in 1:nRegModif) {
    cc[i] <- wInProd(drv,
                    allForK$phi[,i],
                    weight) /
            wInProd(allForK$phi[,i],
                    allForK$phi[,i],
                    weight)
  }
  K <- cc
  for (i in 1:nRegModif) {
    K[i] <- t(cc[i:nRegModif]) %*% as.matrix(allForK$A[i:nRegModif,i])
  }
  resTot <- sum((weight * drv - weight * allForK$Y %*% as.matrix(K))^2)
  # prepare output
  allForK$K <- K
  allForK$resTot <- resTot
  # return
  invisible(allForK)
}
