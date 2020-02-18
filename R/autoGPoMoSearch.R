#' @title Automatic search of polynomial Equations
#'
#' @description This algorithm aims to get an ensemble of possible
#' models which integrability will be tested later with function
#' \code{autoGPoMoTest}. By default, all the terms are considered
#' available (Some of the terms can be excluded intentionally
#' using the option \code{filterReg}).
#' The maximum size of the equation depends on the model dimension
#' \code{nVar}, and on the maximum polynomial degree \code{dMax}.
#' The algorithm removes polynomial terms one by one using a
#' leave-one-out method.
#'
#' @inheritParams  gloMoId
#' @inheritParams  drvSucc
#' @inheritParams  poLabs
#'
#' @param data Input Time series: Each column is one time series
#' that corresponds to one variable.
#' @param dt Time sampling of the input series.
#' @param underSamp Number of points used for undersampling the data.
#' For \code{undersamp = 1} the complete time series is used.
#' For \code{undersamp = 2}, only one data out of two is kept, etc.
#'
#' @return A list of two matrices:
#' @return \code{$filtMemo} describes the selected terms
#' (1 if the term is used, 0 if not)
#' @return \code{$KMemo} provides the corresponding coefficients
#'
#' @author Sylvain Mangiarotti, Flavie Le Jean
#'
#' @seealso \code{\link{autoGPoMoTest}}, \code{\link{gPoMo}},
#'          \code{\link{findAllSets}}, \code{\link{poLabs}}
#' 
#' @export
#' @examples
#' # Load data
#' data('RosYco')
#' # Search for potential models
#' filt = autoGPoMoSearch(RosYco[,2], nVar = 3, dMax = 2,
#'                        dt = 1/125, show = 1)
#' # As an example, the equations of the fourth line has the following terms:
#' poLabs(nVar = 3, dMax = 2)[filt$filtMemo[5,] != 0]
#' # which coefficients correspond to
#' cbind(filt$KMemo[5,], poLabs(nVar = 3, dMax = 2))[filt$filtMemo[5,] != 0,]
#'
autoGPoMoSearch <- function (data, dt, nVar = nVar, dMax = dMax, weight = NULL,
                             show = 0, underSamp = NULL, filterReg = NULL)
{
  if (!is.null(underSamp)) {
    sechdata = data[seq(1, dim(data)[1], by = underSamp), ]
  }
  else {
    sechdata = data
  }
  if (min(dim(as.matrix(sechdata))) == 1) {
    reg <- gloMoId(sechdata, nVar= nVar, dt = dt, dMax = dMax, weight = weight,
                 show = 0, filterReg = filterReg == 1)
  }
  else {
    reg <- gloMoId(sechdata, dt = dt, dMax = dMax, weight = weight,
                 show = 0, filterReg = filterReg == 1)
  }
  testFin = 0
  Memo <<- list()
  filtMemo <- matrix(1, nrow = 1, ncol = length(reg$filterReg))
  filtMemo[1,] <- reg$filterReg
  reg <- gloMoId(reg$init, dt = dt, dMax = dMax, weight = weight,
               show = 0, filterReg = reg$filterReg == 1)
  KMemo <- reg$K
  while (!testFin) {
    getRidOf <- (reg$resSsMod/reg$resTot - 1) == min((reg$resSsMod/reg$resTot - 1))
    if (sum(getRidOf) > 0 & sum(reg$filterReg[reg$filterReg == 1]) > 2) {
      reg$filterReg[reg$filterReg == 1][getRidOf == 1] = 0
      #            poLabs(nVar, dMax, reg$filterReg == 0)
      reg <- gloMoId(reg$init, dt = dt, dMax = dMax, weight = weight, show = 0,
                   filterReg = reg$filterReg == 1)
    }
    else if (sum(getRidOf) > 0 & sum(reg$filterReg[reg$filterReg == 1]) == 2) {
      reg$filterReg[reg$filterReg == 1][getRidOf == 1] <- 0
    }
    else {
      testFin = 1
    }
    if (sum(reg$filterReg) != sum(filtMemo[dim(filtMemo)[1],])) {
      filtMemo <- rbind(filtMemo, reg$filterReg)
      Ktmp <- reg$filterReg
      Ktmp[reg$filterReg==1] <- reg$K[reg$filterReg==1]
      KMemo <- rbind(KMemo, Ktmp)
    }
  }
  #
  tab <- regSeries(nVar, dMax, reg$init[,1:nVar])[, reg$filterReg ==  1]
  coeff <- lsfit(tab, reg$init[, nVar + 1], intercept = F)$coefficients
  KMemo[dim(KMemo)[1],][reg$filterReg ==  1] <- coeff
  #
  # ajout 17/10/2018 (to have all the monomial equal to 0)
  # modified for test 21/03/2019
  KMemo <- rbind(KMemo, Ktmp * 0)
  filtMemo <- rbind(filtMemo, Ktmp * 0)
  
  # preparing data to return back
  Memo$filtMemo <- filtMemo
  Memo$KMemo <- KMemo
  Memo
}



