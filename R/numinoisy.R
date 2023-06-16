#' @title Generates time series of deterministic-behavior
#' with stochatic perturbations (measurement and/or dynamical noise)
#'
#' @description Generates time series from Ordinary Differential Equations
#' perturbed by dynamical and/or measurement noises
#'
#' @importFrom stats sd rnorm
#'
#' @param x0 The initial conditions. Should be a vector which size must be equal
#' to the model dimension \code{dim(K)[2]} (the number of variables of the
#' model defined by matrix \code{K}).
#' @param t A vector providing all the dates for which the output are expected.
#' @param K The Ordinary Differential Equations used to model the dynamics.
#' The number of column should correspond to the number of variables, the
#' number of lines to the number of parameters following the convention
#' defined by \code{poLabs(nVar,dMax)}.
#' @param varData A vector of size \code{nVar} providing the caracteristic
#' variances of each variable of the dynamical systems in ODE defined
#' by matrix \code{K}.
#' If not provided, this variance is automatically estimated.
#' @param txVarBruitA A vector defining the ratio of ADDITIVE noise
#' for each variable of the dynamical system in ODE. The additive noise is
#' added at the end of the numerical integration process. The ratio is
#' defined relatively to the signal variance of each variable.
#' @param txVarBruitM A vector defining the ratio of DYNAMICAL
#' noise for each variable of the dynamical system in ODE. This
#' noise is a perturbation added at each numerical integration step. The
#' ratio is defined relatively to the signal variance of each variable.
#' @param varBruitA A vector defining the variance of ADDITIVE noise
#' for each variable of the dynamical system in ODE. The additive noise is
#' added at the end of the numerical integration process.
#' @param varBruitM A vector defining the variance of DYNAMICAL
#' noise for each variable of the dynamical system in ODE. This
#' noise is a perturbation added at each numerical integration step.
#' @param taux Generates random gaps in time series. Parameter \code{taux}
#' defines the ratio of data to be kept (e.g. for \eqn{taux=0.75}, 75 percents
#' of the data are kept).
#' @param freq Subsamples the time series. Parameter \code{freq} defines the
#' periodicity of data kept (e.g. for \eqn{freq=3}, 1 data out of 3 is kept).
#' @param variables Defines which variables must be generated.
#' @param method Defines the numerical integration method to be used.
#' The fourth-order Runge-Kutta method is used by default
#' (\code{method = 'rk4'}). Other method may be used (such as \code{'ode45'}
#' or \code{'lsoda'}), see function \code{ode} from package \code{deSolve}
#' for details.
#'
#' @return A list of two variables: \cr
#' @return \code{$donnees} The integrated trajectory (first column is the time,
#' next columns are the model variables) \cr
#' @return \code{$bruitM} The level of dynamical noise \cr
#' @return \code{$bruitA} The level of additive noise \cr
#' @return \code{$vectBruitM} The vector of the dynamical noise used to produce
#' the time series \cr
#' @return \code{$vectBruitA} The vector of the additive noise used to produce
#' the time series \cr
#' @return \code{$ecart_type} The level standard deviation \cr
#'
#' @author Sylvain Mangiarotti, Malika Chassan
#'
#' @examples
#' #############
#' # Example 1 #
#' #############
#' # Rossler Model formulation
#' # The model dimension
#' nVar = 3
#'  # maximal polynomial degree
#' dMax = 2
#' a = 0.520
#' b = 2
#' c = 4
#' Eq1 <- c(0,-1, 0,-1, 0, 0, 0, 0, 0, 0)
#' Eq2 <- c(0, 0, 0, a, 0, 0, 1, 0, 0, 0)
#' Eq3 <- c(b,-c, 0, 0, 0, 0, 0, 1, 0, 0)
#' K <- cbind(Eq1, Eq2, Eq3)
#' # Edit the equations
#' visuEq(K, nVar, dMax)
#' # initial conditions
#' v0 <- c(-0.6, 0.6, 0.4)
#' # output time required
#' timeOut = (0:800)/50
#' # variance of additive noise
#' varBruitA = c(0,0,0)^2
#' # variance of multiplitive noise
#' varBruitM = c(2E-2, 0, 2E-2)^2
#' # numerical integration with noise
#' intgr <- numinoisy(v0, timeOut, K, varBruitA = varBruitA, varBruitM = varBruitM, freq = 1)
#' # Plot of the simulated time series obtained
#' dev.new()
#' plot(intgr$donnees[,2], intgr$donnees[,3], type='l',
#'       main='phase portrait', xlab='x(t)', ylab = 'y(t)')
#' dev.new()
#' oldpar <- par(no.readonly = TRUE)    
#' on.exit(par(oldpar))  
#' par(mfrow = c(3, 1))
#' plot(intgr$donnees[,1], intgr$donnees[,2], type='l',
#'       main='phase portrait', xlab='x(t)', ylab = 'y(t)')
#' lines(intgr$donnees[,1], intgr$vectBruitM[,2]*10, type='l',
#'       main='phase portrait', xlab='x(t)', ylab = 'e(t)*10', col='red')
#' plot(intgr$donnees[,1], intgr$donnees[,3], type='l',
#'       main='phase portrait', xlab='x(t)', ylab = 'y(t)')
#' lines(intgr$donnees[,1], intgr$vectBruitM[,3]*10, type='l',
#'       main='phase portrait', xlab='x(t)', ylab = 'e(t)*10', col='red')
#' plot(intgr$donnees[,1], intgr$donnees[,4], type='l',
#'       main='phase portrait', xlab='x(t)', ylab = 'y(t)')
#' lines(intgr$donnees[,1], intgr$vectBruitM[,4]*10, type='l',
#'       main='phase portrait', xlab='x(t)', ylab = 'e(t)*10', col='red')
#'
#' @export
numinoisy <- function(x0, t, K,
                      varData = NULL, txVarBruitA = NULL, txVarBruitM = NULL,
                      varBruitA = NULL, varBruitM = NULL, taux=NULL, freq=NULL,
                      variables=NULL, method=NULL){
  
  # check dimensions
  if (length(x0) != dim(K)[2]) {
    stop("v0 length (=",length(x0),") does not match with the model dimension (=",dim(K)[2],")")
  }
  
  # slvn 24/09/2015
  # EITHER txVarBruitM OR sigBruitM should be given


  if(is.null(taux)& is.null(freq))  stop("saisir un param d'echantillonage")
  if(!is.null(taux) & !is.null(freq)) stop("choisir un seul param d'echantillonage")
  # ajout slvn 24/09/2015
  # check if EITHER txVarBruitM OR varBruitM is given
  if ( is.null(varBruitM) & is.null(txVarBruitM)) stop("either txVarBruitM or varBruitM should be provided")
  # if none is given then NO multiplicative noise
  if ( !is.null(varBruitM) & !is.null(txVarBruitM)) txVarBruitM <- BruitM <- x0*0
  # check if EITHER txVarBruitA OR varBruitA is given
  if ( is.null(varBruitA) & is.null(txVarBruitA)) stop("either txVarBruitA or varBruitA should be provided")
  # if none is given then NO additive noise
  if ( !is.null(varBruitA) & !is.null(txVarBruitA)) txVarBruitA <- BruitA <- x0*0


  ## Ajout de bruit multiplicatif ##
  data <- odeBruitMult2(x0, t, K, varData = varData, txVarBruitM = txVarBruitM, varBruitM = varBruitM, method)
  donnees = data$reconstr
  ecTy = data$ecart_type
  sigma = sqrt(data$bruitM)
  vectBruitM = data$vectBruitM


  ## Ajout de bruit additif ##
  L = length(t)
  nVar = length(x0)
  bruitADD = matrix(data = rnorm(L*nVar,0,1), nrow = L, ncol = nVar)

  # add conditional slvn 24/09/2015
  # if varBruitM is NOT given (it is NULL), it is computed based on txVarBruitM
  rho = varBruitA
  if ( is.null(varBruitA) ) {
    # if ecTy is necessary and NOT given it is recomputed
    if ( is.null(ecTy) ) {
      rec = deSolve::ode(x0, t, derivODE2, K, method)
      ecTy = apply(rec[,2:(nVar+1)],2,sd)
    }
    rho = txVarBruitA * ecTy^2
  }
  RHO = matrix(data = rho, nrow = L, ncol = nVar, byrow=TRUE)

  vectBruitA <- matrix(0,nrow=length(t),ncol=nVar+1)
  vectBruitA[,1] = t
  vectBruitA[,2:(nVar+1)] <- bruitADD*sqrt(RHO)

  donnees[,2:(nVar+1)] = donnees[,2:(nVar+1)] + vectBruitA[,2:(nVar+1)]

  ## Sous echantillonnage ##
  # (pour chq variable separement)


  # Aleatoirement
  if(!is.null(taux)){
    tx = floor(taux*length(t))
    vect = (1:length(t))
    indices = matrix(nrow = tx, ncol = nVar)
    for ( i in 1:nVar){
      indices[,i] = sort(sample(vect)[1:tx])
      temp = donnees[indices[,i],i+1]
      donnees[,i+1] = NA
      donnees[indices[,i],i+1] = temp
    }
  }

  # A frequence fixe
  if(!is.null(freq)){
    nb = floor(length(t)/freq)
    indices = matrix(data = (1:nb)*freq, nrow = nb, ncol = nVar)
    temp = donnees[indices[,1],2:(nVar+1)]
    donnees[,2:(nVar+1)] = NA
    donnees[indices[,1], 2:(nVar+1)] = temp
  }


  ## Selection de variables ##
  if(!is.null(variables)){
    donnees[ , (2:(nVar+1))[-variables] ]= NA
  }

  ## Erreur sur chaque variable
  DONNEES = matrix(nrow = length(t), ncol = 2*nVar+1)
  DONNEES[,1:(nVar+1)] <- donnees
  for (i in 1:nVar){
    DONNEES[indices[,i], nVar+1+i] = rho[i] + sigma[i]
  }


  # Output
  list(donnees = DONNEES, indices = indices, bruitM = sigma, bruitA = rho, vectBruitM = vectBruitM,vectBruitA = vectBruitA, ecart_type=ecTy)
}

