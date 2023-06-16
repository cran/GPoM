#' @title combiEq : Combine Equations from different sources
#'
#' @description Combines equations of different sources
#' into a single system. During this combination, the polynomial
#' maximal degree can be either imposed or optimized to reduce
#' the model size. All the input have to follow
#' the convention defined by \code{poLabs}.
#'
#' @param allKL A list of models, each provided as a matrix.
#' A single matrix can also be provided, it will be transformed
#' into a list containing a single matrix.
#' @param eqOrder A list of vector, providing each the equations
#' number (relating to the input models) to be kept in the
#' output equation system. If not provided, all the equations
#' are kept.
#' A single matrix can also be provided, it will be transformed
#' into a list containing a single matrix.
#' @param dMaxOut The maximal polynomial degree of the output
#' equation system (if not provided, this degree is deduced
#' from the input models)
#'
#' @author Sylvain Mangiarotti
#'
#' @seealso \code{\link{gPoMo}}, \code{\link{poLabs}}
#'
#' @examples
#' # Load models
#' data("allMod_nVar3_dMax2")
#' # Display equations of system 1
#' visuEq(nVar = 3, dMax = 2, K = allMod_nVar3_dMax2$NH86, substit = 1)
#' # Display equations of system 2
#' visuEq(nVar = 3, dMax = 2, K = allMod_nVar3_dMax2$R76, substit = 1)
#' # put the two systems in a list
#' allK <- list()
#' allK[[1]] <- allMod_nVar3_dMax2$NH86
#' allK[[2]] <- allMod_nVar3_dMax2$R76
#' 
#' # Example 1: reformulate two autonomous system in a single matrix 
#' visuEq(K = allK[[1]], substit = c('u', 'v', 'w'))
#' visuEq(K = allK[[2]], substit = c('X', 'Y', 'Z'))
#' Knew <- combiEq(allK)
#' visuEq(K = Knew, substit = c('u', 'v', 'w', 'X', 'Y', 'Z'))
#' 
#' # Example 2
#' inXnote = list()
#' inXnote[[1]] <- c('u', 'v', 'w')
#' inXnote[[2]] <- c('X', 'Y', 'Z')
#' visuEq(K = allK[[1]], substit = inXnote[[1]])
#' visuEq(K = allK[[2]], substit = inXnote[[2]])
#' XnoteOut = c('X', 'Y', 'Z', 'u', 'v', 'w')
#' Knew2 <- combiEq(allK, eqOrder=c(4,5,6,1,2,3))
#' visuEq(K = Knew2, substit = XnoteOut)
#' 
#' # Example 3
#' inXnote = list()
#' inXnote[[1]] <- c('u', 'v', 'w')
#' inXnote[[2]] <- c('X', 'Y', 'Z')
#' visuEq(K = allK[[1]], substit = inXnote[[1]])
#' visuEq(K = allK[[2]], substit = inXnote[[2]])
#' XnoteOut = c('u', 'X', 'v', 'Y', 'w', 'Z')
#' Knew3 <- combiEq(allK, eqOrder=c(1,4,2,5,3,6), dMaxOut = 3)
#' visuEq(K = Knew3, substit = XnoteOut)
#' @export
#' 
#' @return KLout A matrix of the combined model
combiEq <- function(allKL, eqOrder = NULL, dMaxOut = NULL) {
  
    if (is.matrix(allKL)) {
      tmp <- list()
      tmp$allKL <- allKL
      allKL <- tmp
    }
  
    nbSystems = length(allKL)

    # sets the final system actual value of dMax (either provided as dMaxOut or retrieved as the
    # maximum of dMax of all systems) 
    if (is.null(dMaxOut)) {
      dMax = 0
      dMaxOut = 0
      for (i in 1:nbSystems) {
        KL = allKL[[i]]
        pMax_i = dim(KL)[1]
        nVar_i = dim(KL)[2]
        dMax_i = p2dMax(nVar_i, pMax_i)
        poLabsOrder = colSums(regOrd(nVar_i, dMax_i))
        dMaxOut_i = max(ifelse(abs(KL) > 0, poLabsOrder, 0))
        #dMax = max(dMax, dMax_i)
        dMaxOut = max(dMaxOut, dMaxOut_i)
      }
    }
    # else {
    #   dMax = dMaxOut
    # }    
    # sets the final system number of variables nVar
    nVar = 0
    for (i in 1:nbSystems) {
      nVar = nVar + dim(allKL[[i]])[2]
    }    

    # sets the total number of polynomial terms
    pMax = d2pMax(nVar, dMaxOut)
    
    # sets the actual order of variables
    if (!is.null(eqOrder)) { 
        XnoteOut = letters[eqOrder]
    }
    else {
        XnoteOut = letters[1:nVar]
    }
    KL_tmp = NULL
    firstIndex = 1
    for (i in 1:nbSystems) {
       nVar_i = dim(allKL[[i]])[2]
       dMax_i = p2dMax(dim(allKL[[i]])[2], dim(allKL[[i]])[1])
       
       # find positions of system 1 rows in final system

       rows = match(poLabs(nVar_i, dMax_i, Xnote  = letters[(firstIndex):(firstIndex+nVar_i-1)]), 
                    poLabs(nVar, dMaxOut, Xnote = XnoteOut))
       
       firstIndex = firstIndex + nVar_i
       # creates a new matrix with the system re written with pMax lines
       KL_i_new = matrix(0L, nrow = pMax, ncol = nVar_i)

      for (j in 1:length(rows)) {
          KL_i_new[rows[[j]],] = allKL[[i]][j,]
      }
      if (is.null(KL_tmp)) {
        KL_tmp = KL_i_new
      }
       else {
         KL_tmp = cbind(KL_tmp, KL_i_new)
       }
    }

    KLout <- KL_tmp

    
    if (!is.null(eqOrder)) {
      KLout <- list()
      KL_reordered = matrix(0L, nrow = pMax, ncol = length(eqOrder))
      for (j in 1:length(eqOrder)) {
          KL_reordered[, j] = KL_tmp[, eqOrder[j]]
      }
      KLout= KL_reordered
    }
    

    
    KLout
}
