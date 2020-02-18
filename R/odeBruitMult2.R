#' @title For the numerical integration of ordinary
#' differential equations with dynamical noise.
#'
#' @description A subfunction for the numerical integration of Ordinary
#' Differential Equations provided in a generic polynomial form.
#' Model formulation follows the convention defined
#' by function \code{poLabs}.
#'
#' @inheritParams derivODE2
#' @inheritParams numinoisy
#'
#' @importFrom stats sd rnorm
#'
#' @param x0 Initial conditions
#' @param method Numerical method used in the integration process.
#' (see \code{ode} function in \code{deSolve} package for details).
#'
#' @seealso \code{\link{numinoisy}}
#'
#' @author Sylvain Mangiarotti, Malika Chassan
#'
odeBruitMult2 <- function(x0, t, K, varData = NULL, txVarBruitM = NULL, varBruitM = NULL, method = NULL){

	nVar = length(x0)

	# ajout slvn 24/09/2015
	# check if EITHER txBruitM OR sigBruitM is given
	if ( is.null(varBruitM) & is.null(txVarBruitM)) stop("either txVarBruitM or varBruitM should be provided")
	if ( !is.null(varBruitM) & !is.null(txVarBruitM)) stop("either txVarBruitM or varBruitM should be provided, not both")

	rec = ode (x0, t, derivODE2, K, method)
	vectBruitM <- matrix(0,nrow=length(t),ncol=nVar+1)
	vectBruitM[,1] = t
	if ( is.null(varData) ) {
	  ecTy = apply(rec[floor(length(t)/2):length(t),2:(nVar+1)],2,sd)
	}
	else {
	  ecTy = sqrt(varData)
	}

	# add conditional slvn 24/09/2015
	# if varBruitM is NOT given (it is NULL), it is computed based on txBruitM
	if ( is.null(varBruitM) ) {
	  varBruitM = txVarBruitM * ecTy^2
	}


	nul = which(varBruitM == 0)
	if ( length(nul) == nVar ){
	  vectBruitM <- rec
	  vectBruitM[,2:(nVar+1)] <- 0
	  }

	if ( length(nul) < nVar ){

	  x = x0
	  rec = c(t[1], x0)

	  for (i in 1:(length(t)-1) ){

	    temp = ode (x, t[i:(i+1)], derivODE2, K, method)
	    rec = rbind (rec, temp[2:dim(temp)[1],])

	    x = rec[dim(rec)[1], 2:(nVar+1)]
	    bruitM = matrix(data = rnorm(nVar,0,1), nrow = 1, ncol = nVar)

	    x = x + bruitM * sqrt(varBruitM)

	    vectBruitM[i+1,2:(nVar+1)] <- bruitM * sqrt(varBruitM)

	  }
	}

list(reconstr = rec, ecart_type = ecTy, bruitM = varBruitM, vectBruitM= vectBruitM )

}

