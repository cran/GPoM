#' @title subSysD : Sub-systems Disentangling
#'
#' @description Detect, disentangle and reformulate Sub-systems
#' from an ensemble of equations.
#'
#'
#' @param inK A list of models, each provided as a matrix.
#' A single matrix can also be provided, it will be transformed
#' into a list containing a single matrix.
#' @param inXnote A vector with the names of the input
#' variables. If not provided, default notation
#' is used: "X1", "X2", etc.
#'
#' @author Sylvain Mangiarotti
#'
#' @seealso \code{\link{gPoMo}}, \code{\link{poLabs}}, \code{\link{combiEq}}
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
#' # Example 1 (two independant subsystems)
#' # take two separate systems and mix them
#' inXnote = list()
#' inXnote[[1]] <- c('u', 'v', 'w')
#' inXnote[[2]] <- c('X', 'Y', 'Z')
#' visuEq(K = allK[[1]], substit = inXnote[[1]])
#' visuEq(K = allK[[2]], substit = inXnote[[2]])
#' XnoteOut = c('u', 'X', 'v', 'Y', 'w', 'Z')
#' Knew3 <- combiEq(allK, inXnote = inXnote, XnoteOut = XnoteOut, dMaxOut = 3)
#' visuEq(K = Knew3, substit = XnoteOut)
#' # Disentangle the subsystems from the mixed equations
#' dstgl <- subSysD(Knew3, inXnote = XnoteOut)
#' ## Optional
#' # library(igraph)
#' # g1<-graph.adjacency(dstgl$FM);
#' # l <- layout_with_fr(g1)
#' # plot(g1, edge.arrow.siez = .4, edge.curved=.4, vertex.label=XnoteOut, layout = l)
#' 
#' # Example 2 (one subsystem included in the other)
#' Kduff <- matrix(0, ncol = 4, nrow = 35)
#' Kduff[11,1] <- Kduff[5,2] <- Kduff[2,3] <- 1
#' Kduff[35,2] <- -1
#' Kduff[11,2] <- -0.05
#' Kduff[5,4] <- 2 * acos(-1) / 6.2
#' Xnote <- c("x", "y", "u", "v")
#' visuEq(Kduff, substit = Xnote)
#' dstgl2 <- subSysD(Kduff, inXnote = Xnote)
#' 
#' @export
subSysD <- function(inK, inXnote = NULL) {
  
  #############################################################
  #############################################################
  ### TESTS
  
  nVar <- dim(inK)[2]
  pMax <- dim(inK)[1]
  dMax <- p2dMax(nVar = nVar, pMaxKnown = pMax)
  
  ### Prepare the variables names (if not provided)
  # or test their coherence (if provided)
  if (is.null(inXnote)) {
    # prepare variables labels
    varLab <- rev(poLabs(nVar,1)[-1])
  }
  else {
    varLab <- inXnote
    if (nVar != length(varLab)) {
      stop('Input models dimensions is incompatible with variables labels number')
    }
  }

  #############################################################
  #############################################################
  ### MAIN PROGRAM
  
  # prepare Fluence matrix
  FM <- matrix(0, nrow = nVar, ncol = nVar)
  for (i in 1:nVar) {
    terms <- which(inK[,i] != 0)
    terms <- terms[terms != 1]
    for (j in 1:length(terms)) {
      whatVar <- which(regOrd(nVar = nVar, dMax = dMax)[,terms[j]] >= 1)
      FM[whatVar,i] <- 1
    }
  }
  
  # Detec autonomous
  FMtmp <- FM
  for (i in 1:nVar) {
    FMtmp[i,i] <- 1
  }
  # Individual Variable Dependance IVD
  IVD <- list()
  for (i in 1:nVar) {
    IVD[[i]] <- FMtmp[,i]
    for (j in 1:nVar) {
      what <- which(IVD[[i]] == 1 & FMtmp[,j] == 1)
      what <- what[what != 1]
      if (length(what)!=0) {
        for (k in 1:length(what)) {
          IVD[[i]][which(FMtmp[,what[k]] == 1)] <- 1
        }
      }
    }
  }
  
  # SubSystems disentangling
  allSubS <- list()
  allSubS[[1]] <- IVD[[1]]
  iplus = 1
  for (i in 1:length(IVD)) {
    keep = 1
    for (j in 1:length(allSubS)) {
      if (length(IVD[[i]]) == length(allSubS[[j]])) {
        if(sum(IVD[[i]] != allSubS[[j]]) == 0) {
          keep = 0
        }
      }
    }
    if (keep == 1) {
      iplus = iplus + 1
      allSubS[[iplus]] <- IVD[[i]]
    }
  }
  
  # Display the subsystems
  cat(length(allSubS), 'sub-systems can be distinguished:', "\n")
  for (isys in 1:length(allSubS)) {
    #vlab <- rev(poLabs(nVar,1)[-1])
    Eq <- combiEq(inK = inK,
                  inXnote = varLab,
                  XnoteOut = varLab[which(allSubS[[isys]]!=0)])
    # Display
    cat('### Sub-system', isys, ':', "\n")
    visuEq(K = Eq, substit = varLab[which(allSubS[[isys]]!=0)])
  }
  
  subS <- list()
  subS$FM <- FM
  subS$FMtmp <- FMtmp
  subS$IVD <- IVD
  subS$allSubS <- allSubS
  
  subS
}
