#' @title drvSucc : Computes the successive derivatives of a time series
#'
#' @description Computes the successive derivatives from one single time series,
#' using the Savitzky-Golay algorithm (1964).
#'
#' @inheritParams gloMoId
#'
#' @param serie A single time series provided as a single vector.
#' @param nDeriv The number of derivatives to be computed from
#' the input time series. The resulting number of
#' time series obtained in output will be \code{nDeriv + 1}.
#' @param tstep Sampling time of the input time series. Used
#' only if time vector \code{tin} is not provided.
#' @param winL Number (exclusively odd number) of points of
#' the local window used for computing the derivatives along
#' the input time series. The Savitzky-Golay filter is used for
#' this purpose [1,2].
#'
#' @return A list containing:
#' @return $serie The original time serie
#' @return $tin The time vector containing the dates corresponding
#' to the original time series
#' @return $tstep The time step (assumed to be regular)
#' @return $tout The time vector of the output series
#' @return seriesDeriv A matrix containing the original time series
#' (smoothed by the filtering process) in the first column
#' and its \code{nDeriv + 1} successive derivatives in the next ones.
#' Note that \code{winL} values of the original time series will be lost,
#' that is \code{(winL - 1)/2} at the begining and \code{(winL - 1)/2}
#' at the end of the time series due to a computation boundary effect).
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
#' @seealso \code{\link{gloMoId}}, \code{\link{gPoMo}}, \code{\link{poLabs}}, \code{\link{compDeriv}}
#'
#' @examples
#' #############
#' # Example 1 #
#' #############
#' # Generate a time series:
#' tin <- seq(0, 5, by = 0.01)
#' data <- 2 * sin(5*tin)
#' dev.new()
#' oldpar <- par(no.readonly = TRUE)    
#' on.exit(par(oldpar))  
#' par(mfrow = c(3, 1))
#' # Compute its derivatives:
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
#' # load data:
#' data("Ross76")
#' tin <- Ross76[,1]
#' data <- Ross76[,2]
#'
#' # Compute the derivatives
#' drvOut <- drvSucc(tin, data, nDeriv=4)
#' dev.new()
#' oldpar <- par(no.readonly = TRUE)    
#' on.exit(par(oldpar))  
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
