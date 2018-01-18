#' @title Detection of limit cycles of period-1
#'
#' @description This algorithm aim to detect period-1 limit cycles
#' from trajectories in the phase sapce considered in a
#' bidimensional projection.
#'
#' @inheritParams autoGPoMoSearch
#' @inheritParams autoGPoMoTest
#'
#' @param data A matrix of the trajectory in a 2D space (if more than
#' two columns are provided, only the two first columns are considered)
#' @param LimCyclThreshold The detection threshold
#' @param show Indicates the deepness of the feedback (from 0 to 2)
#'
#' @importFrom graphics plot lines
#' @importFrom stats spline
#'
#' @author Sylvain Mangiarotti
#'
#' @return Indicates if a limit cycle is detected (1) or not (0)
#'
#' @seealso \code{\link{autoGPoMoTest}}
#'
# ###########
# # Example #
# ###########
#
# # Two time series are prepared
# tin <- seq(0, 20, by = 0.01)
# # the first one is periodic of period-1
# x1 <- sin(2*tin)
# # the second one is periodic of period-2
# x2 <- sin(2*tin) + 0.5*cos(4*tin)
# #
# # Their plots
# dev.new()
# plot(tin, x1, type = 'l', col = 'green')
# lines(tin, x2, type = 'l', col = 'red')
# #
# # Their derivatives are computed
# drv1 <- drvSucc(tin = tin, serie = x1, nDeriv = 3)
# drv2 <- drvSucc(tin = tin, serie = x2, nDeriv = 3)
# # their phase portraits are plotted
# plot(drv2$seriesDeriv[,1], drv2$seriesDeriv[,2], type = 'l', col = 'red')
# lines(drv1$seriesDeriv[,1], drv1$seriesDeriv[,2], type = 'l', col = 'green')
# #
# # The detection of period-1 limit cycle is then tested
# # The period-1 time series is detected as period-1
# detectP1limCycl(drv1$seriesDeriv, LimCyclThreshold = 0.01, show = 0)
# # The period-2 time series is not detected as period-1
# detectP1limCycl(drv2$seriesDeriv, LimCyclThreshold = 0.01, show = 0)
#
detectP1limCycl <- function (data, LimCyclThreshold = 0.01, show=2)
{
  izP1limCycl <- 0
  bary <- c(mean(data[, 1]), mean(data[, 2]))
  ro <- (data[, 1] - bary[1])^2 + (data[, 2] - bary[2])^2
  phi <- atan2(data[, 2] - bary[2], data[, 1] - bary[1])
  if (show == 2) {
    plot(phi, ro, col='blue')
  }
  ndata <- dim(data)[1]
  what <- (abs(phi - c(phi[1], phi[-ndata])) > 2)
  if (sum(what) > 3) {
    num <- rep(1:ndata)[rep(1:ndata) * what != 0]
    nChg <- length(num)
    roLast <- ro[(num[nChg - 1] + 1):(num[nChg] - 2)]
    phiLast <- phi[(num[nChg - 1] + 1):(num[nChg] - 2)]
    roAnte <- ro[(num[nChg - 2] + 1):(num[nChg - 1] - 2)]
    phiAnte <- phi[(num[nChg - 2] + 1):(num[nChg - 1] - 2)]
    if (length(phiLast) > 5) {
      phiInt <- phiLast[2:(length(phiLast) - 2)]
      if (sum(phiAnte > min(phiInt)) != 0 & sum(phiAnte <
                                                max(phiInt)) != 0) {
        roLastInt <- spline(phiLast, roLast, xout = phiInt)
        roAnteInt <- spline(phiAnte, roAnte, xout = phiInt)
        diff <- roLastInt$y - roAnteInt$y
        #
        if (show == 2) {
          plot(phiInt,roLastInt$y,main="ro=f(phi)", col='green')
          lines(phiInt,roAnteInt$y,main="ro=f(phi)", col='red')
          plot(phiInt,diff,main="ro=f(phi)", col='blue')
        }
        #
        if (sqrt(t(diff) %*% diff) <= LimCyclThreshold) {
          izP1limCycl <- 1
        }
      }
    }
  }
  izP1limCycl
}


