#' @title combiEq : Combine Equations from different sources
#'
#' @description Combines equations of different sources
#' into a single system. During this combination, the polynomial
#' maximal degree can be either imposed or optimized to reduce
#' the model size. All the input have to follow
#' the convention defined by \code{poLabs}.
#'
#' @param inK A list of models, each provided as a matrix.
#' A single matrix can also be provided, it will be transformed
#' into a list containing a single matrix.
#' @param inXnote A list of vectors with the names of the input
#' variables for each model. If not provided, default notation
#' is used: "X1", "X2", etc.
#' A single matrix can also be provided, it will be transformed
#' into a list containing a single matrix.
#' @param eqNum A list of vector, providing each the equations
#' number (relating to the input models) to be kept in the
#' output equation system. If not provided, all the equations
#' are kept.
#' A single matrix can also be provided, it will be transformed
#' into a list containing a single matrix.
#' @param XnoteOut A vector with the names of the output variables.
#' If not provided, default notation is used considering that
#' the variables of the input models are all different
#' @param nVarOut The dimension of the output equation system
#' (if not provided, this degree is deduced from the input models)
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
#' Knew2 <- combiEq(allK, inXnote = inXnote, XnoteOut = XnoteOut)
#' visuEq(K = Knew2, substit = XnoteOut)
#' 
#' # Example 3
#' inXnote = list()
#' inXnote[[1]] <- c('u', 'v', 'w')
#' inXnote[[2]] <- c('X', 'Y', 'Z')
#' visuEq(K = allK[[1]], substit = inXnote[[1]])
#' visuEq(K = allK[[2]], substit = inXnote[[2]])
#' XnoteOut = c('u', 'X', 'v', 'Y', 'w', 'Z')
#' Knew3 <- combiEq(allK, inXnote = inXnote, XnoteOut = XnoteOut, dMaxOut = 3)
#' visuEq(K = Knew3, substit = XnoteOut)
#' 
#' # Example 4
#' dim(Knew3)
#' inXnote = c('x', 'X', 'y', 'Y', 'z', 'Z')
#' visuEq(K = Knew3, substit = inXnote)
#' XnoteOut = c('X', 'Y', 'Z')
#' Knew4 <- combiEq(Knew3, inXnote = inXnote, XnoteOut = XnoteOut)
#' dim(Knew4)
#' visuEq(K = Knew4, substit = XnoteOut)
#' 
#' @export
combiEq <- function(inK, inXnote = NULL, eqNum = NULL,
                    XnoteOut = NULL, nVarOut = NULL, dMaxOut = NULL) {

#############################################################
#############################################################
### TESTS
  
### inK, inXnote and eqNum should be lists
  if (is.list(inK) == 0) {
    tmp <- inK
    inK <- list()
    inK[[1]] <- tmp
  }
  if (is.null(inXnote) == 0 & is.list(inXnote) == 0) {
    tmp <- inXnote
    inXnote <- list()
    inXnote[[1]] <- tmp
  }
  if (is.null(eqNum) == 0 & is.list(eqNum) == 0) {
    tmp <- eqNum
    eqNum <- list()
    eqNum[[1]] <- tmp
  }
  
### Calculate all the input models:
  # (1) number of variables nVar,
  # (2) number of formulation terms pMax
  # and (3) polynomial degree formulation dMax
  inKnVar <- c()
  inKpMax <- c()
  inKdMax <- c()
  for (i in 1:length(inK)) {
    inKnVar[i] <- dim(inK[[i]])[2]
    inKpMax[i] <- dim(inK[[i]])[1]
    inKdMax[i] <- p2dMax(nVar = inKnVar[i], pMaxKnown = inKpMax[i])
  }
  
### Prepare the variables names (if not provided)
  # or test their coherence 
  if (is.null(inXnote)) {
    # prepare variables labels
    inXnote <- list()
    varLab <- poLabs(sum(inKnVar),1)[-1][sum(inKnVar):1]
    tmp1 <- c(0,inKnVar[1:(length(inKnVar)-1)])
    tmp2 <- inKnVar
    for (i in 1:length(inK)) {
      inXnote[[i]] <- varLab[(sum(tmp1[1:i])+1):sum(tmp2[1:i])]
    }
  }
  else {
    # test coherence
    for (i in 1:length(inK)) {
      if (sum(inKnVar[i]) != length(inXnote[[i]])) {
        stop('Input models dimensions is incompatible with variables labels number')
      }
    }
  }
  
### Prepare the variables selection (if not provided, all are kept)
  # or test their coherence
  if (is.null(eqNum) & is.null(XnoteOut)) {
    # all the variables are kept
    XnoteOut <- unlist(inXnote)
    #
    eqNum <- list()
    tmp1 <- c(0,inKnVar[1:(length(inKnVar)-1)])
    tmp2 <- inKnVar
    for (i in 1:length(inK)) {
      eqNum[[i]] <- (sum(tmp1[1:i])+1):sum(tmp2[1:i])
    }
    nVarNew <- length(XnoteOut)
    icol = 1:nVarNew
    tokeep <- rep(0,sum(inKnVar))
    tokeep[icol] <- order(icol)
  }
  else if (is.null(eqNum) & is.null(XnoteOut) == 0) {
    # XnoteOut is used to prepare EqNum
    icol <- c()
    for (i in 1:length(XnoteOut)) {
      icol <- c(icol,which(unlist(inXnote) == XnoteOut[i]))
    }
    tokeep <- rep(0,sum(inKnVar))
    tokeep[icol] <- order(icol)
    tmp <- tokeep
    for (imod in 1:length(inK)) {
      eqNum[[imod]] <- tmp[1:inKnVar[imod]]
      tmp <- tmp[(inKnVar[imod]+1):length(tokeep)]
    }
    nVarNew <- length(XnoteOut)
  }
  else if (is.null(eqNum) == 0 & is.null(XnoteOut)) {
    # eqNum is used to prepare XnoteOut
    XnoteOut <- c()
    tmp <- unlist(inXnote)
    for (i in 1:length(tmp)) {
      XnoteOut <- c(XnoteOut, tmp[unlist(eqNum)])
    }
    nVarNew <- length(XnoteOut)
    icol <- c()
    for (i in 1:length(XnoteOut)) {
      icol <- c(icol,which(unlist(inXnote) == XnoteOut[i]))
    }
    tokeep <- rep(0,sum(inKnVar))
    tokeep[icol] <- order(icol)
  }
  else {
    # test coherence
    for (i in 1:length(inK)) {
      if (length(eqNum[[i]]) > inKnVar[i]) {
        stop('Input models dimensions is incompatible with output variables number')
      }
    }
    icol <- c()
    for (i in 1:length(XnoteOut)) {
      icol <- c(icol,which(unlist(inXnote) == XnoteOut[i]))
    }
    tokeep <- rep(0,sum(inKnVar))
    tokeep[icol] <- order(icol)
  }
  
### output model dimension
  if (is.null(nVarOut)) {
    # prepare the output model dimension
    nVarOut <- nVarNew
  }
  else {
    # test coherence
    if (nVarNew != nVarOut) {
      stop('Output model dimension is incompatible with nVarOut')
    }
  }
  
### test for coherence of dMaxOut to be added
  dmx <- 0
  if (is.null(dMaxOut)) {
    # calculate the minimum polynomial degree required
    for (imod in 1:length(inK)) {
      for (ieq in 1:sum(eqNum[[imod]]!=0)) {
        # test if one equation is concerned
        if (sum(eqNum[[imod]] != 0) != 0) {
          # take one equation
          oneEq <- inK[[imod]][,eqNum[[imod]] != 0][,ieq]
          # what coeff not zero
          notz <- which(oneEq != 0)
          # max polynomial degree?
          dmx <- max(dmx, colSums(regOrd(inKnVar[imod], inKdMax[imod]))[notz])
        }
      }
    }
    dMaxOut <- dmx
  }
  else {
    # test coherence
    if (dMaxOut < max(inKdMax)) {
      warning('Output polynomial degree is potentially incompatible with dMaxOut')
    }
  }

#############################################################
#############################################################
### MAIN PROGRAM
  
  # New model struture
  pMaxOut <- d2pMax(nVar = nVarOut, dMaxKnown = dMaxOut)
  
  KLnew <- matrix(data = 0, nrow = pMaxOut, ncol = nVarOut)
  ieqNew = 1
  for (imod in 1:length(inK)) {
    for (ieq in 1:length(inXnote[[imod]])) {
      for (term in (which(inK[[imod]][,which(eqNum[[imod]]!=0)[ieq]] != 0))) {
        inew <- which(poLabs(nVar = nVarOut,
                             dMax = dMaxOut,
                             Xnote = XnoteOut)
                          ==
                      poLabs(nVar = inKnVar[imod],
                             dMax = inKdMax[imod],
                             Xnote = inXnote[[imod]])[term]
                      )
        KLnew[inew,eqNum[[imod]][which(eqNum[[imod]]!=0)][ieq]] <- inK[[imod]][,which(eqNum[[imod]]!=0)[ieq]][term]
      }
    ieqNew <- ieqNew + 1
    }
  }

  KLnew <- KLnew[,tokeep[icol]]
  #
  KLnew
}
