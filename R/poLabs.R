#' @title Polynomial labels order
#'
#' @description Defines the order of the polynomial labels given
#' the number of variables \code{nVar} and the maximum polynomial
#' degree \code{dMax}.
#' @seealso \code{\link{visuEq}}
#'
#' @inheritParams regOrd
#'
#' @param Xnote Enables to defines the notation used for the variable,
#' by default \code{Xnote = 'X'}.
#' @param findIt A vector of selected terms.
#'
#' @return \code{lbls} A vector of characters. Each element is the
#' expression of one polynomial term, such as \eqn{X_1^2 X_3 X_4}
#'
#' @author Sylvain Mangiarotti
#'
#' @examples
#' #Regressor order for three variables \eqn{(X1,X2,X3)} (nVar = 3) for a maximum
#' #polynomial degree equal to 2 (dMax = 2): poLabs(3,2)
#' #and for two variables only : poLabs(2,2)
#'
#' # For a quadratic equation of two variables,
#' # the polynomial \deqn{P(X1,X2) = 0.5 + 0.3 X1 -0.25 X1 X2}
#' # could thus be written as a vector Pvec such as:
#' Pvec = c(0.5, 0, 0, 0.3, -0.25, 0)
#' # considering the convention corresponding to
#' poLabs(2,2)
#' # Indeed:
#' poLabs(2, 2, findIt = Pvec!=0)
#' # An alternative notation can be used with parameter Xnote
#' poLabs(2, 2, findIt = Pvec!=0, Xnote = 'w')
#' # or also
#' poLabs(2, 2, findIt = Pvec!=0, Xnote = c('x','y'))
#'
#' @export
poLabs <- function(nVar, dMax, dMin = 0, findIt = NULL, Xnote = 'X') {

  # notation
  if (length(Xnote) == 1) {
    Xnote <- rep(Xnote, nVar)
    for (j in 1:nVar) {
      Xnote[j] <- paste(Xnote[j],j, sep="")
    }
  }
  # test Xnote length
  if (length(Xnote) != nVar) {
    stop('Xnote should be either one single character or a  nVar length vector of character')
  }
  # number of regressors
  ##pMax <- d2pMax(nVar, dMax)
  # regressors Order :
  pExpo <- regOrd(nVar, dMax, dMin = dMin)
  # number of regressors
  pMax <- dim(pExpo)[2]
  #
  lbls <- vector("character",pMax)
  #
  # by definition, first regressor is the constant
  lbls[1] <- "ct"
  #
  # other regressors are reconstructed one by one
  for (i in 2:pMax) {
    # each factor of the regressor
    for (j in 1:nVar) {
      if (pExpo[j,i] > 0) {
        # exponent equal to 1 are eluded
        if (pExpo[j,i] == 1)
            lbls[i] <- paste(lbls[i],Xnote[j]," ",sep="")
        # exponent greater than 1 are noted
        if (pExpo[j,i] > 1)
           lbls[i] <- paste(lbls[i],Xnote[j],"^",pExpo[j,i]," ",sep="")
      }
    }
  }
  if (dMin < 0) {
    for (i in 2:pMax) {
      if (pExpo[1,i] == -1) {
        if (lbls[i] == "") {lbls[i] = "1"}
        lbls[i] <- paste(lbls[i], "/", Xnote[1])
      }
    }
  }
  #
  if (!is.null(findIt)) {
    lbls <- lbls[findIt]
  }
  # return
  lbls
}
