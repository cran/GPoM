#' @title deriveODE2 : A Subfonction for the numerical integration
#' of polynomial equations in the generic form defined by function
#' \code{poLabs}
#'
#' @description This function provides the one step integration of
#' polynomial Ordinary Differential Equations (ODE). This function
#' requires the function \code{ode} ("deSolve" package).
#'
#' @inheritParams gloMoId
#'
#' @param t All the dates for which the result of the numerical integration
#' of the model will have to be provided
#' @param x Current state vector (input from which the next state will
#' be estimated)
#' @param K is the model: each column corresponds to one equation
#' which organisation is following the convention given by function
#' \code{poLabs} which requires the definition of the model dimension
#' \code{nVar} (i.e. the number of variables) and the maximum polynomial
#' degree \code{dMax} allowed.
#' @param regS Current states of each polynomial terms used
#' in \code{poLabs}. These states can be deduced from the current
#' state vector x (using function \code{regSeries}). When available,
#' it can be provided as an input to avoid unecessary computation.
#'
#' @author Sylvain Mangiarotti
#'
#' #@export
derivODE2 <- function(t, x, K, regS = NULL) {

    if (is.null(regS)) {
      # if not provided regressors are recomputed
      # from state vector x:
      # nVar and dMax are required use the appropriate convention
      dMax <- p2dMax(length(x), dim(K)[1])
      # compute the regressors values
      regS <- regSeries(length(x), dMax, x)
      temp <- NULL
    }
    else {
      temp <- regS[2:dim(regS)[1],]
      # only the first state is used
      regS <- regS[1,]
    }
  # compute the increment
  derives <- regS %*% K

  # prepare output
  list(derives, regS = temp)
}
