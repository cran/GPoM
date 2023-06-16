#' @title Displays the models Equations
#'
#' @description Displays the model equations for a polynomial model
#' which description is provided as a matrix \code{K}, each column
#' corresponding to one equation. The coefficients of the polynomial
#' terms are given following the order defined by function \code{poLabs}.
#' The matrix can also be provided in a list \code{K},
#' in this case, the matrix should be located in \code{K$model[[selecmod]]}
#' where selecmod should be provided as input parameter.
#'
#' @inheritParams poLabs
#' @inheritParams derivODE2
#' @inheritParams gloMoId
#'
#' @param substit Applies subtitutions to the default values:
#' for \code{substit = 0} (default value), variables are chosen
#' as \code{X1}, \code{X2}, ...
#' for \code{substit = 1}, variable \code{X1}, \code{X2}, ...
#' will be replaced by \code{x}, \code{y}, \code{z}, ...
#' for \code{substit = 2}, the codes provides a LaTex-like
#' formulation of the model.
#' The variables name can also be defined explicitely as follows:
#' for \code{substit = c('x', 'H', 'T1')}, variables \code{X1},
#' \code{X2}, \code{X3} ... will be replaced by \code{x}, \code{H}
#' and \code{T1}.
#' @param approx The number of extra digits to be used:
#' for \code{approx = FALSE} (default value) digits are
#' edited with double precision;
#' for \code{approx = TRUE}, only the minimum number of digits is
#' edited (in order to have all the terms different from 0)
#' for \code{approx} = 1, 2, etc. then respectively 1, 2, etc.
#' digits are added to the minimum number of digits corresponding
#' to \code{approx = TRUE}.
#' @param selecmod An integer providing the number in the sublist
#' when the model matrix is provided in a list. Should not be
#' provided (or NULL) if the model matrix is provided directly.
#'
#' @author Sylvain Mangiarotti
#'
#' @examples
#' #EQUATIONS VISUALISATION
#' # number of variables:
#' nVar <- 3
#' # maximum polynomial degree:
#' dMax <- 2
#' # polynomial organization:
#' poLabs(nVar,dMax)
#' # model construction
#' KL = matrix(0, ncol = 3, nrow = 10)
#' KL[1,1] <- KL[2,2] <- 1
#' KL[4,1] <- -1
#' KL[5,3] <- -0.123456789
#' # Equations visualisation:
#' # (a) by default, variables names X1, X2, X3 are used
#' visuEq(KL, nVar, dMax)
#' # (b) for susbstit=1, variables names x, y, y are used instead
#' visuEq(KL, nVar, dMax, approx = TRUE, substit=1)
#' # (c) the name of the variables can also be chosen manualy
#' visuEq(KL, nVar, dMax, approx = 3, substit=c('U', 'V', 'W'))
#'
#' # A canonical model can be provided as a single vector
#' polyTerms <- c(0.2,0,-1,0.5,0,0,0,0,0,0)
#' visuEq(polyTerms, 3,2)
#'
#' @export
#' @return N A list of Nvar elements presenting a set of equtions,
#' each equation being an element of this list and written as
#' a character strings

visuEq <- function (K, nVar=NULL, dMax=NULL, dMin = 0, substit = 0, approx = FALSE, selecmod = NULL) {
  if (is.list(K)) {
    if (is.null(selecmod)) {
      stop('The model number should be provided using option: selecmod = ')
      }
    if (length(selecmod)!=1) {
      stop('selecmod should be an integer. Only a single model can be displayed')
    }
    K <- K$models[[selecmod]]
#    nVar <- dim(K)[2]
#    dMax <- p2dMax(nVar = nVar, pMaxKnown = dim(K)[1])
  }
  if (is.vector(K)) {
    K <- cano2M(nVar, dMax, K, dMin = dMin)
  }
  if (is.matrix(K)) {
    if (is.null(nVar)) nVar <- dim(K)[2]
    if (is.null(dMax)) dMax <- p2dMax(nVar = nVar, pMaxKnown = dim(K)[1], dMin=dMin)
  }
  # Check compatibility
  if (dim(K)[2] != nVar) {
    stop('nVar = ', nVar,
         ' is not compatible with the model dimension dim(K)[2] = ',
         dim(K)[2])
  }
  if (d2pMax(nVar,dMax, dMin=dMin) != dim(K)[1]) {
    stop('dMax = ', dMax,
         ' is not compatible with the model size dim(K)[1] = ',
         dim(K)[1])
  }
  # for approximation with all terms
    if (approx == TRUE) {
      lePlusPetit <- min(abs(K)[abs(K) != 0])
      digits <- abs(floor(log10(lePlusPetit)))
      K <- round(K, digits)
    }
  else if ((approx %% 1) == 0 & approx != FALSE) {
    lePlusPetit <- min(abs(K)[abs(K) != 0])
    digits <- abs(floor(log10(lePlusPetit)))
    K <- round(K, digits + approx)
  }
  else if ((approx %% 1) != 0 & !is.logical(approx)) {
    stop("parameter 'approx' should be either a logical or an integer value")
  }
  #
  N <- list()
  for (i in 1:nVar) {
    if (!is.numeric(substit[1])) {
      varsList = substit
    }
    else {
      if (substit == 0) {
        varsList = c("X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9")
      }
      else if (substit == 1) {
        varsList = c("x", "y", "z", "w", "u", "v", "r", "s", "q")
      }
      else if (substit == 2) {
        varsList = c("X_1", "X_2", "X_3", "X_4", "X_5", "X_6", "X_7", "X_8", "X_9")
        
      }
      else {
        stop("substit must be either a list of nVar values or equal to 0 or 1 or 2")
      }
    }
    dN <- paste("d", varsList[i],"/dt =", sep="")

    M <- cbind(K[K[,i]!=0,i], poLabs(nVar, dMax, dMin = dMin, findIt=K[,i]!=0, Xnote=varsList[1:nVar]), "+")
 
    M <- matrix(t(M), nrow=1, ncol=length(M))
    M <- M[,1:(dim(M)[2]-1)]
    M <- cbind(dN, t(M))
    M <- paste(M,collapse=" ")
    for (i in 1:nVar) {    
       M <- gsub(paste("1 ", varsList[i]), varsList[i], M, fixed=TRUE)
    }
    M <- gsub(" + -", " -", M, fixed=TRUE)
    M <- gsub("+  -", " -", M, fixed=TRUE)
    M <- gsub("+ -", " -", M, fixed=TRUE)
    M <- gsub(" +-", " -", M, fixed=TRUE)
    M <- gsub("ct", "", M)
    M <- gsub("Ct", "", M)
    if (substit[1] == 2) {
      # substitutions needed to agree with latex notation
      M <- gsub("e-0", "e-", M)
      M <- gsub("e+0", "e+", M)
      M <- gsub("e", ".10^{", M)
      M <- gsub("+-", "-", M)
      M <- gsub("+ -", "-", M)
    }

    # plot expression ...
    message(M, "\n\n")
    N <- append(N, M)
  }
  invisible(N)
}






