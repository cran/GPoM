#' @title A subfonction for the numerical integration
#' of polynomial equations provided in a generic form
#' following the convetion defined by function \code{poLabs}.
#'
#' @description This function provides the one step integration of
#' polynomial Ordinary Differential Equations (ODE). This function
#' requires the function \code{ode} (\code{deSolve} package).
#'
#' @inheritParams regOrd
#'
#' @param t All the dates for which the result of the numerical
#' integration of the model must be provided
#' @param x Current state vector (input from which the next state will
#' be estimated)
#' @param K A matrix providing the model description:
#' each column corresponds to one equation which polynomial organisation
#' is following the convention defined by function \code{poLabs}.
#' @param regS Current states of each polynomial terms used
#' in \code{poLabs}. These states can be deduced from the current
#' state vector \code{x} (using the function \code{regSeries}).
#' When available, it can be provided as an input to avoid
#' unecessary computation.
#'
#' @author Sylvain Mangiarotti
#'
#' @seealso \code{\link{numicano}}, \code{\link{numinoisy}}
#'
derivODE2 <- function(t, x, K, dMin = 0, regS = NULL) {

    if (is.null(regS)) {
      # if not provided regressors are recomputed
      # from state vector x:
      # nVar and dMax are required use the appropriate convention
      dMax <- p2dMax(length(x), dim(K)[1], dMin=dMin)
      ###dMax <- 2
      # compute the regressors values
      regS <- regSeries(length(x), dMax, x, dMin = dMin)
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
