#' @title bDrvFilt : Build Derivative Filter
#'
#' @description Build the Savitzky-Golay derivative filter (Savitzky-Golay, 1964).
#'
#' @param nDrv The number of derivatives to be computed.
#' @param tstep Time sampling.
#' @param winL The local window length to be used for computing the derivatives [1].
#'
#' @references
#' [1] Savitzky, A.; Golay, M.J.E.,
#' Smoothing and Differentiation of Data by Simplified Least Squares Procedures.
#' Analytical Chemistry 36 (8), 1627-1639, 1964.\cr
#'
#' @author Sylvain Mangiarotti
#'
# #############
# # Example 1 #
# #############
# # Apply derivatives
# dFlt <- bDrvFilt(nDrv = 3, tstep = 0.1, winL = 5)
#'
#' @export
bDrvFilt <- function(nDrv, tstep, winL = 9) {
  
  if (winL %% 2 == 0) {
    stop('filter length parameter winL should be an odd number')
  }
  # Build derivative filter
  X <- NULL
  halfW <- (winL - 1) / 2
  for (i in 0:nDrv) {
    X <- cbind(X, ((-halfW:halfW) * tstep)^i)
  }
  dFlt <- solve(t(X) %*% X) %*% t(X)
  # return
  dFlt
}
