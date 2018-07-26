#' @title Concat Concatenates separated time series
#'
#' @description The aim of this code is to provide, from a set of multiple
#' time series, a single concatenated time series for applying the global
#' modeling technique to all the time time series in association.
#'
#' @inheritParams  gloMoId
#'
#' @param svrlTS All separated time series.
#'
#' @return \code{concaTS} The concatenated time series.
#'
#' @importFrom grDevices dev.new
#'
#' @author Sylvain Mangiarotti, Mireille Huc
#'
#' @examples
#' # load data
#' data("svrlTS")
#' # Concatenate the data set into a single time series
#' winL = 55
#' concaTS <- concat(svrlTS, winL = winL)
#' # Plot the concatenated time series
#' plot(concaTS$sglTS$TS[,1], concaTS$sglTS$TS[,2],
#'      main = 'Concatenated time series',
#'      xlab = 'Time (concatenated)', ylab = 'y(t)',
#'      type = 'l', col = 'gray')
#' lines(concaTS$sglTS$TS[concaTS$sglTS$W == 1,1],
#'       concaTS$sglTS$TS[concaTS$sglTS$W == 1,2], type = 'p', col = 'green', cex = 0.5)
#' lines(concaTS$sglTS$TS[concaTS$sglTS$W == 0,1],
#'       concaTS$sglTS$TS[concaTS$sglTS$W == 0,2], type = 'p', col = 'red', cex = 0.5)
#' lines(concaTS$sglTS$TS[,1], concaTS$sglTS$W, type = 'l')
#' \dontrun{
#' # The concatenated data set can be used for global modelling:
#' GPout1 <- gPoMo(data = concaTS$sglTS$TS[,2], tin = concaTS$sglTS$TS[,1],
#'                 dMax = 2, nS = 3, winL = winL, weight = concaTS$sglTS$W, show = 1,
#'                 IstepMin = 10, IstepMax = 6000, nPmin = 11, nPmax = 11, method = 'rk4')
#' }
#'
#' @references S. Mangiarotti, F. Le Jean, M. Huc & C. Letellier, 2016.
#' Global modeling of aggregated and associated chaotic dynamics,
#' Chaos, Solitons & Fractals, 83, 82-96.
#'
#' @export
#'
concat <- function (svrlTS, winL = 9){
  
  # How many time series in input
  nTS <- length(svrlTS)
  
  for (i in 1:nTS) {
    if (is.null(svrlTS[[i]]$W)) {
      svrlTS[[i]][[2]] <- svrlTS[[i]][[1]][,1] * 0 + 1
    }
  }
  
  # number of data lost at each boundary
  lost <- (winL - 1) / 2
  
  tps0 <- 0
  W0 <- c()
  sglTS <- list()
  sglTS$W <- c()
  for (i in 1:nTS) {
    # extract time vector
    tps <- svrlTS[[i]][[1]][,1] - svrlTS[[i]][[1]][1,1]
    # extract data time series
    TS  <- svrlTS[[i]][[1]][,2]
    
    # if weight function is not provided:
    if (is.null(svrlTS[[i]]$W)) {
      W   <- tps * 0 + 1
    }
    # and if it is:
    else {
      # extract weight vector (0 or 1)
      W   <- svrlTS[[i]][[2]]
    }
    
    # reject boudary values
    W[1:lost] <- 0
    W[(length(W)-lost+1):length(W)] <- 0
    # concatenate the last time series at the end of the single time series
    sglTS$TS <- rbind(sglTS$TS, cbind(tps + tps0, TS))
    W0 <- c(W0, W)
    # new tps0
    tps0 <- tail(sglTS$TS[,1],1) + tps[2] - tps[1]
  }
  sglTS$W <- W0

  
  # check the interspace
  N <- dim(sglTS$TS)[1]
  maxdif <- max(sglTS$TS[2:N,1] - sglTS$TS[1:(N-1),1])
  mindif <- min(sglTS$TS[2:N,1] - sglTS$TS[1:(N-1),1])
  if (  (maxdif - mindif) / maxdif * 100 > 0.01 ) {
    stop('Stop: irregular time vector.')
  }
  
  concaTS <- list()
  # input data is also kept for memory
  concaTS$svrlTS <- svrlTS
  concaTS$sglTS <- sglTS
  return(concaTS)
}

