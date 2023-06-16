#' @title extractEq : Extracts Equations from one system
#'
#' @description Combines equations of different sources
#' into a single system. During this combination, the polynomial
#' maximal degree can be either imposed or optimized to reduce
#' the model size. All the input have to follow
#' the convention defined by \code{poLabs}.
#'
#' @param KL A model, provided as a matrix.
#' @param eqVect A vector of integers, providing the equations
#' numbers to be kept in the output equation system. 
#' 
#' @author Mireille Huc
#' 


extractEq <- function(KL, eqVect) {
  
  # dimensions of the original system
  pMaxIn= dim(KL)[1]
  nVarIn = dim(KL)[2]  
  
  # number of polynomial terms of the original system
  #pMaxIn = d2pMax(nVarIn, dMax)
  dMax = p2dMax(nVarIn, pMaxIn)
  
  # number of variables of the extracted system
  nVarOut = length(eqVect)

  poLabsOrder = colSums(regOrd(nVarIn, dMax))
  dMaxOut = max(ifelse(abs(KL) > 0, poLabsOrder, 0))

  # the dimensions of the final system are then known, the final matrix can be built
  pMaxOut = d2pMax(nVarOut, dMaxOut)
  KL_new= matrix(0L, nrow = pMaxOut, ncol = nVarOut)
  
  rows = match(poLabs(nVarOut, dMaxOut, Xnote = letters[eqVect]),
               poLabs(nVarIn, dMax, Xnote  = letters[1:nVarIn]))

  for (i in 1:pMaxOut) {
     for (j in 1:nVarOut) {
       KL_new[i, j] = KL[rows[i], eqVect[j]]
     }
  }

  KL_new
}
