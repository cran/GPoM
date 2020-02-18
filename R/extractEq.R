extractEq <- function(KL, eqExtract) {
  
  # dimensions of the original system
  pMaxIn= dim(KL)[1]
  nVarIn = dim(KL)[2]  
  # number of polynomial terms of the original system
  #pMaxIn = d2pMax(nVarIn, dMax)
  dMax = p2dMax(nVarIn, pMaxIn)
  
  # number of variables of the extracted system
  nVarOut = length(eqExtract)

  poLabsOrder = colSums(regOrd(nVarIn, dMax))
  dMaxOut = max(ifelse(abs(KL) > 0, poLabsOrder, 0))

  # the dimensions of the final system are then known, the final matrix can be built
  pMaxOut = d2pMax(nVarOut, dMaxOut)
  KL_new= matrix(0L, nrow = pMaxOut, ncol = nVarOut)
  
  rows = match(poLabs(nVarOut, dMaxOut, Xnote = letters[eqExtract]),
               poLabs(nVarIn, dMax, Xnote  = letters[1:nVarIn]))

  for (i in 1:pMaxOut) {
     for (j in 1:nVarOut) {
       KL_new[i, j] = KL[rows[i], eqExtract[j]]
     }
  }

  KL_new
}