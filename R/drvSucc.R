#' @title drvSucc : Computes the successive derivatives of a time series
#'
#' @description Computes the successive derivatives from one single time series,
#' using the Savitzky-Golay algorithm (1964).
#'
#' @inheritParams gloMoId
#'
#' @param serie A single time series provided as a single vector.
#' @param nDeriv The number of derivatives to be computed from
#' the input \code{series}. The resulting number of time series
#' obtained as an output will thus be nVar = nDeriv + 1.
#' @param tstep Time sampling of the input \code{series}. Used
#' only if time vector \code{tin} is not provided.
#' @param winL Number (odd) of points of the local window used
#' for computing the derivatives along the input time series
#' (\code{series}). The Savitzky-Golay filter is used for this
#' purpose [1,2].
#'
#' @return A list containing:
#' $serie The original time serie
#' $tin The time vector corresponding to the original time series
#' $tstep The time step (assumed to be invariant)
#' $tout The time vector of the output series
#' ^seriesDeriv A matrix containing the original variable (as smoothed
#' by the filtering process) and its \code{nDeriv} + 1 first
#' derivatives (note that winL values of the original time series
#' will be lost both (winL - 1)/2 at the begining and at the end
#' of the time series due to boundary effect).
#'
#' @references
#' [1] Savitzky, A.; Golay, M.J.E.,
#' Smoothing and Differentiation of Data by Simplified Least Squares Procedures.
#' Analytical Chemistry 36 (8), 1627-1639, 1964.\cr
#' [2] Steinier J., Termonia Y., Deltour, J.
#' Comments on smoothing and differentiation of data by simplified least square procedure.
#' Analytical Chemistry 44 (11): 1906-1909, 1972. \cr
#'
#' @author Sylvain Mangiarotti, Mireille Huc
#'
#' @seealso \code{\link{gloMoId}}, \code{\link{gPoMo}}, \code{\link{poLabs}}
#'
#' @examples
#' #############
#' # Example 1 #
#' #############
#' tin <- seq(0, 5, by = 0.01)
#' data <- 2 * sin(5*tin)
#' dev.new()
#' par(mfrow = c(3, 1))
#' # numerical derivation
#' drv <- drvSucc(tin = tin, nDeriv = 2, serie = data, winL = 5)
#' #
#' # plot original and filtered series
#' plot(tin, data, type='l', col = 'black', xlab = 't', ylab = 'x(t)')
#' lines(drv$tout, drv$seriesDeriv[,1], lty = 3, lwd = 3, col = 'green')
#' #
#' # analytic 1st derivative
#' firstD <- 10 * cos(5 * tin)
#' # plot both
#' plot(tin, firstD, type = 'l', col = 'black', xlab = 't', ylab = 'dx/dt')
#' lines(drv$tout, drv$seriesDeriv[,2], lty = 3, lwd = 3, col = 'green')
#' #
#' # analytic 2nd derivative
#' scdD <- -50 * sin(5 * tin)
#' # plot both
#' plot(tin, scdD, type = 'l', col = 'black', xlab = 't', ylab = 'd2x/dt2')
#' lines(drv$tout, drv$seriesDeriv[,3], lty=3, lwd = 3, col = 'green')
#'
#' #############
#' # Example 2 #
#' #############
#' # load data
#' data("Ross76")
#' tin <- Ross76[,1]
#' data <- Ross76[,2]
#'
#' # Apply derivatives
#' drvOut <- drvSucc(tin, data, nDeriv=4)
#' dev.new()
#' par(mfrow = c(3, 1))
#' # original and smoothed variable:
#' plot(drvOut$tin, drvOut$serie,
#'      type='p', cex = 1, xlab = 'time', ylab = 'x(t)')
#' lines(drvOut$tout, drvOut$seriesDeriv[,1], type='p', col='red')
#' lines(drvOut$tout, drvOut$seriesDeriv[,1], type='l', col='red')
#' # 1st derivative:
#' plot(drvOut$tout, drvOut$seriesDeriv[,2],
#'      type='p', col='red', xlab = 'time', ylab = 'dx(t)/dt')
#' lines(drvOut$tout, drvOut$seriesDeriv[,2], type='l', col='red')
#' # 2nd derivative:
#' plot(drvOut$tout, drvOut$seriesDeriv[,3],
#'      type='p', col='red', xlab = 'time', ylab = 'd2x(t)/dt2')
#' lines(drvOut$tout, drvOut$seriesDeriv[,3], type='l', col='red')
#'
#'
#' @export
drvSucc <- function(tin = NULL, serie, nDeriv, weight = NULL, tstep = NULL, winL=9) {
  drvOut <- list()
  if (!is.vector(serie)) stop("serie should be a vector")
  if (is.null(tin)) {
    if (is.null(tstep)) {
      warning("tin and tstep are chosen by default such as:
                tstep = 0.01 and tin = seq(1:length(serie), by=tstep)")
      tstep = 0.01
    }
    tin = seq(from = 0, to = length(serie), by = tstep)
  }

  if (winL < (2 * nDeriv + 1)) {
    warning("winL too small compared to nDeriv, winL is reset to (2*nDeriv+1)")
    winL = (2 * nDeriv + 1)
  }

  drvOut$tin <- tin
  drvOut$serie <- serie
  if (is.null(tstep)) {
    drvOut$tstep <- tstep <- tin[3] - tin[2]
  }
  else {
    if(tin[3] - tin[2] != tstep) stop("tin and tstep incompatible")
  }
  #

  drvOut$tout <- tin[((winL+1)/2):(length(tin) - (winL+1)/2)]
  if (!is.null(weight)) {
    finalWeight <- weight
    for (i in ((winL+1)/2):(length(tin) - (winL+1)/2) ) {
      finalWeight[i] = min(weight[(i-(winL+1)/2): (i+(winL+1)/2)])
    }
    drvOut$Wout <- finalWeight[((winL+1)/2):(length(tin) - (winL+1)/2)]
  }
  else {
    drvOut$Wout <- NULL
  }
  drvOut$seriesDeriv <- compDeriv(serie, nDeriv, tstep, winL = winL)

  drvOut
}
