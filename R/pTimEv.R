#' @title Model stationnary testing
#'
#' @description Estimate the parameters variations 
#' of a model of canonical form considering a sliding
#' window on an external dataset.
#'
#' @importFrom stats quantile

#' @inheritParams  gPoMo
#' @inheritParams  gloMoId
#' @inheritParams  drvSucc
#' @inheritParams  poLabs
#' @inheritParams  autoGPoMoSearch
#'
#' @param TS The time series to be tested
#' @param TSdate The time vector
#' @param wlength The window length
#' @param onestep Step length between two estimations
#' @param whatTerms The terms to be considered in the analysis. Note 
#' that these are organised following the convention defined by
#' poLabs(nVar,dMax). Since only the structure is required, if
#' coefficients are provided, these are transformed to 1.
#' @param removeExtr Ratio of estimated values to be removed (if
#' chosen equal to 0.1, only 90% of the estimates of smaller
#' disersion will be kept)
#'
#' @author Sylvain Mangiarotti
#'
#' @return A list containing:
#' @return \code{$slidingoutGM} An n*(pMax+1) matrix presenting
#' the pMax estimated parameters p1(t), p2(t) etc. column by column.
#' The residual signal epsilon(t) is provided in the last (i.e. pMax + 1)
#' column. Each line correspond to one date provided in \code{$TSdate}
#' @return \code{$TSdate} A time vector relating to the estimates
#' presented in \code{$slidingoutGM}
#' @return \code{$W} A vector providing the output values that
#' can kept (=1) or must be removed (=0)
#' @return \code{$whatTerms} A vector recalling the terms
#' taken into account in the analysis (their order refers
#' to \code{poLabs(nVar,dMax)} function)
#' @return \code{$param} A vector with the parameter values
#' used to apply the function: nVar, dMax, wlength, onestep,
#' removeExtr
#'
#' @seealso \code{\link{autoGPoMoSearch}}, \code{\link{gPoMo}}, \code{\link{poLabs}}
#'
#' @examples
#' #Example
#' data(TS)
#' plot(TS[,1], TS[,2], type='l')
#' nVar <- 3
#' dMax <- 2
#' pMax <- choose(nVar+dMax,dMax)
#' whatTerms <- c(1,1,0,1,1,1,1,1,1,1)
#'  
#' # apply pTimEv
#' statio <- pTimEv(TS[,2], nVar, dMax, TS[,1], whatTerms = whatTerms, 
#'                  wlength = 1000, onestep = 20, removeExtr = 0.15)
#' # Plot the results
#' dev.new()
#'   layout(matrix(1:12, nrow=4, ncol=3, byrow = TRUE))
#'   what <- which(statio$whatTerms!=0)
#'   for (i in what) {
#'       plot(statio$TSdate[statio$W==1], statio$slidingoutGM[statio$W==1,i],
#'            xlab='TSdate', ylab='coeff', main=poLabs(nVar,dMax)[i])
#'      }
#'   plot(statio$TSdate[statio$W==1], statio$slidingoutGM[statio$W==1,pMax+1],
#'        xlab='date', ylab='Epsilon', main='Resid', log = 'y')
#'        
#' @export
pTimEv <- function(TS, nVar, dMax, TSdate, whatTerms = NULL, weight = NULL,
                       wlength = 1000, onestep = 100, removeExtr = 1) {
    #
    pMax <- choose(nVar+dMax,dMax)
    if (is.null(whatTerms) ) {whatTerms <- rep(1,pMax)}
    if (is.null(weight) ) weight <- TSdate * 0 + 1
    #
    # Calculate the number of time steps
    nw <- floor((length(TS) - wlength) / onestep) + 1
    #
    # prepare the analysis matrix
    slidingoutGM <- matrix(0, ncol=2+length(whatTerms), nrow=nw)
    #
    # loop on the sliding window
    for (i in 1:nw) {
      #
      # Extracts the local time series
      series <- TS[(((i-1)*onestep + 1) : ((i-1)*onestep + wlength))]
      W  <- weight[(((i-1)*onestep + 1) : ((i-1)*onestep + wlength))]
      #
      # Estimates the local coefficients
      outGM <- gloMoId(series, dt = TSdate[2]-TSdate[1], nVar = nVar, dMax = dMax, filterReg = whatTerms, weight = W)
      #
      # Keep the results in the analysis matrix
      slidingoutGM[i,] <- cbind(TSdate[(i-1)*onestep + 1] + (TSdate[2]-TSdate[1])*wlength/2,
                                t(c(outGM$K,outGM$resTot)))
    }
    
    W <- rep(1, dim(slidingoutGM)[1])
    what <- (outGM$K != 0) * (1:pMax)
    what <- what[what!=0]
    for (i in what) {
      lesQ <- quantile(slidingoutGM[,i+1], probs = c(0,removeExtr,1-removeExtr,1))
      W[slidingoutGM[,i+1] <= min(lesQ[2], lesQ[3])] <- 0
      W[slidingoutGM[,i+1] >= max(lesQ[2], lesQ[3])] <- 0
    }
    
    outlist <- list()
    outlist$slidingoutGM <- slidingoutGM[,2:dim(slidingoutGM)[2]]
    outlist$TSdate <- slidingoutGM[,1]
    outlist$W <- W
    outlist$whatTerms <- whatTerms
    outlist$param <- c(nVar, dMax, wlength, onestep, removeExtr)
    
    outlist
  }
  