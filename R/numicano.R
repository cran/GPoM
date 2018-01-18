#' @title Numerical Integration of models in ODE of polynomial form
#'
#' @description Function for the numerical integration
#' of Ordinary Differential Equations of polynomial form.
#'
#' @inheritParams autoGPoMoSearch
#' @param Istep The number of integration time steps
#' @param onestep Time step length
#' @param KL Matrix formulation of the model to integrate numerically
#' @param PolyTerms Vectorial formulation of the model (only for models
#' of canonical form)
#' @param v0 The initial conditions (a vector which length should correspond
#' to the model dimension \code{nVar})
#' @param method The integration method (See package \code{deSolve}),
#' by default \code{method = 'rk4'}.
#'
#' @return A list of two variables: \cr
#' @return \code{$KL} The model in its matrix formulation \cr
#' @return \code{$reconstr} The integrated trajectory (first column is the time,
#' next columns are the model variables)
#'
#' @author Sylvain Mangiarotti
#'
#' @examples
#' #############
#' # Example 1 #
#' #############
#' # For a model of general form (here the rossler model)
#' # model dimension:
#' nVar = 3
#' # maximal polynomial degree
#' dMax = 2
#' # Number of parameter number (by default)
#' pMax <- d2pMax(nVar, dMax)
#' # convention used for the model formulation
#' poLabs(nVar, dMax)
#' # Definition of the Model Function
#' a = 0.520
#' b = 2
#' c = 4
#' Eq1 <- c(0,-1, 0,-1, 0, 0, 0, 0, 0, 0)
#' Eq2 <- c(0, 0, 0, a, 0, 0, 1, 0, 0, 0)
#' Eq3 <- c(b,-c, 0, 0, 0, 0, 0, 1, 0, 0)
#' K <- cbind(Eq1, Eq2, Eq3)
#' # Edition of the equations
#' visuEq(nVar, dMax, K)
#' # initial conditions
#' v0 <- c(-0.6, 0.6, 0.4)
#' # model integration
#' reconstr <- numicano(nVar, dMax, Istep=1000, onestep=1/50, KL=K,
#'                       v0=v0, method="ode45")
#' # Plot of the simulated time series obtained
#' dev.new()
#' plot(reconstr$reconstr[,2], reconstr$reconstr[,3], type='l',
#'       main='phase portrait', xlab='x(t)', ylab = 'y(t)')
#'
#' \dontrun{
#' #############
#' # Example 2 #
#' #############
#' # For a model of canonical form
#' # model dimension:
#' nVar = 4
#' # maximal polynomial degree
#' dMax = 3
#' # Number of parameter number (by default)
#' pMax <- d2pMax(nVar, dMax)
#' # Definition of the Model Function
#' PolyTerms <- c(281000, 0, 0, 0, -2275, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#'                861, 0, 0, 0, -878300, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
#' # terms used in the model
#' poLabs(nVar, dMax, PolyTerms!=0)
#' # initial conditions
#' v0 <- c(0.54, 3.76, -90, -5200)
#' # model integration
#' reconstr <- numicano(nVar, dMax, Istep=500, onestep=1/250, PolyTerms=PolyTerms,
#'                      v0=v0, method="ode45")
#' # Plot of the simulated time series obtained
#' plot(reconstr$reconstr[,2], reconstr$reconstr[,3], type='l',
#'      main='phase portrait', xlab='x', ylab = 'dx/dt')
#' # Edition of the equations
#' visuEq(nVar, dMax, reconstr$KL)
#'}
#'
#' @seealso \code{\link{derivODE2}}, \code{\link{numinoisy}}
#'
#' @export
numicano = function(nVar, dMax, Istep=1000, onestep=1/125, KL=NULL, PolyTerms=NULL,
                    v0=NULL, method="rk4") {

  pMax <- d2pMax(nVar, dMax)
  # check integer
  if (is.null(KL) & is.null(PolyTerms)) {
    stop("more information is required (either KL or PolyTerms)")
  }
  if ((!is.null(KL)) & (!is.null(PolyTerms))) {
    stop("too much information is given (chose either KL or PolyTerms)")
  }
  if (is.null(KL) & (!is.null(PolyTerms))) {
    # Definition of the KL Model structure from PolyTerms
    KL <- matrix(0, nrow=pMax, ncol=nVar)
    KL[,nVar] <- PolyTerms
    for (i in 2:nVar) {
      nx <- sum((poLabs(nVar,dMax) == paste("X", i, " ", sep="")) * rep(1:pMax))
      KL[nx,i-1] = 1
    }
  }

  # numerical integration:
  reconstr <- ode(v0, (0:(Istep-1))*onestep, derivODE2, KL, method=method)

  outNumicano <- list()
  outNumicano$KL <- KL
  outNumicano$reconstr <- reconstr

  return(outNumicano)

}
