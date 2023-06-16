#' @title Computes the successive derivatives of a time series
#'
#' @description Computes the successive derivatives from one single time series,
#' with the Savitzky-Golay approach (1964).
#'
#' @param TS A single time series provided as a single vector.
#' @param nDrv The number of derivatives to be computed from the input
#' \code{series}. The resulting number of outpout time series
#' will thus be \code{nVar = nDrv + 1}.
#' @param tstep Sampling Time of the input time series \code{TS}.
#' @param winL The local window length used for computing
#' the derivatives [1-2].
#'
#' @return drv A matrix containing the original variable (smoothed by the
#' filtering process) in the first comlumn and its \code{nDrv}+1 first derivatives
#' in the next columns (note that \code{winL} values of the original time series
#' will be lost both at the begining and the end of the time series due to boundary
#' effect).
#'
#' @references
#' [1] Savitzky, A.; Golay, M.J.E.,
#' Smoothing and Differentiation of Data by Simplified Least Squares Procedures.
#' Analytical Chemistry 36 (8), 1627-1639, 1964.\cr
#' [2] Steinier J., Termonia Y., Deltour, J.
#' Comments on smoothing and differentiation of data by simplified least square procedure.
#' Analytical Chemistry 44 (11): 1906-1909, 1972. \cr
#'
#' @author Sylvain Mangiarotti
#'
#' @seealso \code{\link{gloMoId}}, \code{\link{gPoMo}}, \code{\link{poLabs}}
#'
#'
#' @examples
#' # load data:
#' data(NDVI)
#'
#' # Compute the derivatives:
#' drv <- compDeriv(NDVI[,1], nDrv = 3, tstep = 1/125)
#'
#' @export
compDeriv <- function(TS, nDrv, tstep, winL = 9) {

  TS <- as.matrix(TS)

  # Build derivative filter
  dFlt <- bDrvFilt(nDrv = nDrv, tstep = tstep, winL = winL)

  # Apply filter along the original time series
  drv <- NULL
  for (i in 1:(length(TS) - winL)) {
    drv <- cbind(drv,
                 factorial(0:nDrv) * (dFlt %*% TS[i:(i + winL - 1), 1])
                 )
  }
  drv <- t(drv)
  drv
}
