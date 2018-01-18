#' @title Estimates the monomial time series
#'
#' @description Creates time series by multiplying given time series among them.
#'
#' @inheritParams autoGPoMoSearch
#' @param series A matrix containing the original time series from which
#' the monomials are built. Each column corresponds to one given variable.
#' @param pReg A matrix filled, for each column, with powers of time series
#' used to create.
#'
#' @return \code{rpFull} A matrix of time series. Each column corresponds to one
#' regressor such as \eqn{X_1^2 X_3 X_4}
#'
#' @author Sylvain Mangiarotti
#'
#' @examples
#' data(TSallMod_nVar3_dMax2)
#' sprottK <- as.matrix(TSallMod_nVar3_dMax2$SprK$reconstr)[,2:4]
#' dMax <- 2
#' nVar <- dim(sprottK)[2]
#'
#' #Example 1
#' polySeries2 <- regSeries(nVar, dMax, sprottK)
#'
#' #Example 2
#' p <- c(1,3,1)
#' polySeries2 <- regSeries(nVar, dMax, sprottK, pReg=p)
#'
#' @export
regSeries <- function(nVar, dMax, series, pReg = NULL) {

  #
  if (is.vector(series)) {
    series <- t(series)
  }
  if (is.data.frame(series)) {
    series <- as.matrix(series)
  }

  #Determination des puissances auxquelles elever les series
  if (is.null(pReg)) {
    if (is.null(dMax)) {
      stop("'dMax' or 'pReg' is required.")
    }
    pReg <- regOrd(nVar,dMax)
  } else {
    if (is.vector(pReg)) {
      pReg <- as.matrix(pReg)
    }
  }

  # Compute the regressors time series
  RpFull <- NULL
  for (i in 1:nVar) {
    # computation is performed variable by variable
    Rp <- c()
    for (k in 1:dim(pReg)[2]) {
      # compute Xi^a corresponding to regressor k
      # exponent is in pExpo
      R1 <- series[,i]^pReg[i,k]
      Rp <- cbind(Rp, R1)
    }
    # take into account the product of the ith variable
    # at each new iteration such as: . * Xi^n
    if (is.null(RpFull)) {
      RpFull <- Rp
    }
    else {
      RpFull <- RpFull * Rp
    }
  }
  # return
  RpFull
}
