#' @title Numerical Integration polynomial ODEs with Multiple eXternal forcing
#'
#' @description Function for the numerical integration
#' of Ordinary Differential Equations of polynomial form
#' including single or Multiple external forcing
#'
#' @inheritParams autoGPoMoSearch
#' @param Istep The number of integration time steps. By default,
#' Istep = 1000
#' @param onestep The time step to be used for numerical integration
#' @param KDf The nonautonomous model in its matrix formulation,
#' NA (i.e. not available) values should be provided for forcing
#' variables provided as an external signal
#' @param v0 The initial conditions. Its length should be in agreement
#' with the dynamical system dimension. Therefore, 0 or NA can be
#' provided for external forcing
#' @param method integration method. By default 'rk4' is used
#' @param extF A matrix providing the time vector in the first column,
#' and time series of each forcing in the next ones
#'
#' @return A list of two variables: \cr
#' @return \code{$KDf} The nonautonomous model in its matrix formulation \cr
#' @return \code{$reconstr} The integrated trajectory (first column is the time,
#' next columns are the model variables)
#'
#' @author Sylvain Mangiarotti
#'
#' @examples
#' #############
#' # Example 1 #
#' #############
#' # build a non autonomous model
#' nVar = 4
#' dMax = 3
#' gamma = 0.05
#' KDf=matrix(0, nrow = d2pMax(nVar = nVar, dMax = dMax), ncol = nVar)
#' KDf[11,1]  = 1
#' KDf[2,2]  = 1
#' KDf[5,2]  = 1
#' KDf[11,2]  = -gamma
#' KDf[35,2] = -1
#' KDf[2,3]  = NA
#' KDf[2,4]  = NA
#' visuEq(K = KDf, substit = c('x', 'y', 'u', 'v'))
#' 
#' # build an external forcing
#' # number of integration time step
#' Istep <- 500
#' 
#' # time step
#' smpl <- 1 / 20
#' 
#' # output time vector
#' tvec <- (0:(Istep-1)) * smpl
#' 
#' # angular frequency (for periodic forcing)
#' omega = 0.2
#' 
#' # half step time vector (for Runge-Kutta integration)
#' tvecX <- (0:(Istep*2-2)) * smpl / 2
#' # generate the forcing (here variables u and v)
#' extF = cbind(tvecX, -0.1 * cos(tvecX * omega), 0.05 * cos(tvecX * 16/3*omega))
#' 
#' # decimate the data
#' extFrs <- extF[seq(1,dim(extF)[1],by=50),]
#' extFrs <- rbind(extFrs,extF[dim(extF)[1],])
#' 
#' # Initial conditions to be used (external variables can be set to 0)
#' etatInit <- c(-0.616109362 , -0.126882584 , NA, NA)
#' # model integration
#' out <- numiMultiX(nVar, dMax, Istep=Istep, onestep=smpl, KDf=KDf,
#'                        extF,
#'                        v0=etatInit, method="rk4")
#' outrs <- numiMultiX(nVar, dMax, Istep=Istep, onestep=smpl, KDf=KDf,
#'                        extFrs,
#'                        v0=etatInit, method="rk4")
#' dev.new()
#' par(mfrow = c(2, 2), # 2 x 2 pictures on one plot
#'     pty = "s")
#' plot(out$reconstr[,2],out$reconstr[,3],
#'     xlab = 'x(t)', ylab = 'y(t)', type = 'l', col = 'red')
#' lines(outrs$reconstr[,2],outrs$reconstr[,3],
#'     xlab = 'x(t)', ylab = 'y(t)', type = 'l', col = 'green')
#' plot(out$reconstr[,2],out$reconstr[,4],
#'     xlab = 'x(t)', ylab = 'u(t)', type = 'l', col = 'red')
#' plot(out$reconstr[,4],out$reconstr[,5],
#'     xlab = 'u(t)', ylab = 'v(t)', type = 'l', col = 'red')
#' 
#' @seealso \code{\link{derivODE2}}, \code{\link{numicano}}, \code{\link{numinoisy}}
#'
#' @export
numiMultiX = function(nVar, dMax,
                      Istep=1000, onestep=1/125,
                      KDf,
                      extF = extF,
                      v0=NULL,
                      method="rk4") {
  
  pMax <- d2pMax(nVar, dMax)
  # check dimensions
  if (dim(KDf)[2] != nVar) {
    stop("nVar (=",nVar,") does not match with the model dimension (=",dim(KDf)[2],")")
  }
  if (length(v0) != nVar) {
    stop("v0 length (=",length(v0),") does not match with the model dimension (=",dim(KDf)[2],")")
  }
  # CHECK (and adapt when possible) the external forcing by spline interpolation
  # If the external forcing extF does not have the exact size
  # extF is resampled at the required double time sampling
  # for integration:
  extFmemo <- extF # keep it for memo (return)
  if (dim(extF)[1] != (Istep - 1) * 2 + 1) {
    rsplextF <- matrix(0, ncol = dim(extF)[2], nrow = (Istep - 1) * 2 + 1)
    for (ivar in 1:(dim(extF)[2])) {
      tmp <- spline(x = extF[,1], y = extF[,ivar], xout = (0:(Istep*2-2))*onestep/2)
      rsplextF[,ivar] <- tmp$y
    }
    extF <- rsplextF
  }
  # then check the size consistency
  if (dim(extF)[1] != (Istep - 1) * 2 + 1 ) {
    stop("The forcing time length (",dim(extF)[1]," lignes) does not match with the required double time sampling (=",(Istep - 1) * 2 + 1,")")
  }
  # check integer
  if (method!="rk4") {
    warning("Only 4th order Runge-Kutta method 'rk4' is supported by the numiMultiX function")
  }

  
  # Prepare the external forcing
  # output time vector
  tvec <- (0:(Istep-1)) * onestep
  # hald step time vector (for Runge-Kutta integration)
#  tvecX <- (0:(Istep*2 - 2)) * onestep / 2
  # generate the forcing (here variables u and v)
#  extF = cbind(tvecX, -0.1 * cos(tvecX * omega), 0.05 * cos(tvecX * 16/3*omega))
  #
  # Numerical integration
  reconstr <- reconstr2 <- ode(v0, tvec, derivODEwMultiX,
                  KDf, extF = extF, method = 'rk4')
  # Reconstruction of the output
  iToUpdate <- which(is.na(colSums(KDf)))
  reconstr[, iToUpdate + 1] <- extF[(1:Istep) * 2 - 1,(1:length(iToUpdate)+1)]
  
  outNumiMultiX <- list()
  outNumiMultiX$input$extF <- extFmemo
  outNumiMultiX$input$v0 <- v0
  outNumiMultiX$extF <- extF
  outNumiMultiX$KDf <- KDf
  outNumiMultiX$reconstr <- reconstr
  
  return(outNumiMultiX)
  
}
