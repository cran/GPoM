#' @title visuEq : Displays the model Equations
#'
#' @description Displays the model equations for a polynomial model
#' which description is provided as a matrix K, each column corresponding
#' to an equation, each equation being corresponding to a list of terms
#' following a convention given in `poLabs`.
#'
#' @inheritParams poLabs
#' @inheritParams derivODE2
#' @inheritParams gloMoId
#'
#' @param substit applies subtitutions:
#' for substit = 0 (default value), variables are chosen as X1, X2, ...
#' for substit = 1, variables X1, X2, ... are replaces by x,y,z, ...
#' for substit = 2, the codes provides a LaTex-like formulation of the model
#' @param approx number of digits to be used:
#' for approx = FALSE (default value) digits are edited with double precision 
#' for approx = TRUE, only the minimum number of digits is edited (in order to
#' have all the terms different from 0)
#' for approx = 1, 2, etc. then respectively 1, 2, etc. digits are added to
#' the minimum number of digits corresponding to approx = TRUE.
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
#' visuEq(nVar, dMax, KL)
#' # (b) for susbstit=1, variables names x, y, y are used instead
#' visuEq(nVar, dMax, KL, approx = TRUE, substit=1)
#' # (c) the name of the variables can also be chosen manualy
#' visuEq(nVar, dMax, KL, approx = 3, substit=c('U', 'V', 'W'))
#'
#' # A canonical model can be provided as a single vector
#' polyTerms <- c(0.2,0,-1,0.5,0,0,0,0,0,0)
#' visuEq(3,2,KL)
#'
#' @export
visuEq <- function (nVar, dMax, K, substit = 0, approx = FALSE) {
  if (is.vector(K)) {
    K <- cano2M(nVar, dMax, K)
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
    dN <- paste("dX", i,"/dt =", sep="")
    M <- cbind(K[K[,i]!=0,i], poLabs(nVar, dMax, K[,i]!=0), "+")
    M <- matrix(t(M), nrow=1, ncol=length(M))
    M <- M[,1:(dim(M)[2]-1)]
    M <- cbind(dN, t(M))
    M <- paste(M,collapse=" ")
    M <- gsub(" 1 X", " X", M, fixed=TRUE)
    M <- gsub("-1 X", "-X", M, fixed=TRUE)
    M <- gsub(" + -", " -", M, fixed=TRUE)
    M <- gsub("+  -", " -", M, fixed=TRUE)
    M <- gsub("+ -", " -", M, fixed=TRUE)
    M <- gsub(" +-", " -", M, fixed=TRUE)
    M <- gsub("ct", "", M)
    M <- gsub("Ct", "", M)
    if (substit[1] == 1) {
      M <- gsub("X1", "x", M)
      M <- gsub("X2", "y", M)
      M <- gsub("X3", "z", M)
      M <- gsub("X4", "w", M)
      M <- gsub("X5", "u", M)
      M <- gsub("X6", "v", M)
      M <- gsub("X7", "r", M)
      M <- gsub("X8", "s", M)
      M <- gsub("X9", "q", M)
    }
    if (substit[1] == 2) {
      M <- gsub("X1", "X_1", M)
      M <- gsub("X2", "X_2", M)
      M <- gsub("X3", "X_3", M)
      M <- gsub("X4", "X_4", M)
      M <- gsub("X5", "X_5", M)
      M <- gsub("X6", "X_6", M)
      M <- gsub("X7", "X_7", M)
      M <- gsub("X8", "X_8", M)
      M <- gsub("X9", "X_9", M)
      M <- gsub("e-0", "e-", M)
      M <- gsub("e+0", "e+", M)
      M <- gsub("e", ".10^{", M)
      M <- gsub("+-", "-", M)
      M <- gsub("+ -", "-", M)
    }
    if (!is.numeric(substit[1])) {
      M <- gsub("X1", substit[1], M)
      M <- gsub("X2", substit[2], M)
      M <- gsub("X3", substit[3], M)
      M <- gsub("X4", substit[4], M)
      M <- gsub("X5", substit[5], M)
      M <- gsub("X6", substit[6], M)
      M <- gsub("X7", substit[7], M)
      M <- gsub("X8", substit[8], M)
      M <- gsub("X9", substit[9], M)
    }
    # plot expression ...
    cat(M, "\n\n")
    N[i] <- M
  }

  invisible(N)
}






