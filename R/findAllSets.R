#' @title findAllSets : find all potential systems (a set of
#' equations) from a given ensemble of potential equations
#'
#' @description For each equation to be retrieved, an ensemble
#' of potential formulation is given. For instance, three formulations
#' are provided for equation (1), one for equation (2) and two for
#' equation (3). In this case, six possible sets can be obtained
#' from it.
#' The aim of this program is to formulate all the potential
#' systems from the provided formulation.
#'
#' @inheritParams gPoMo
#'
#' @param allFilt A list with:
#' - a matrix allFilt$Xi of possible formulations for each variable Xi
#' - a vector allFilt$Npi providing the number of parameters of each formulation Xi
#'
#' @author Sylvain Mangiarotti
#'
#' @examples
#' #############
#' # Example 1 #
#' #############
#' allFilt <- list()
#' allFilt$Np1 <- 1         # only one formulation with one single parameter
#' allFilt$Np2 <- c(3, 4)   # two potential formulations, one with respectively
#' # three and four parameters
#' allFilt$Np3 <- c(2, 4)   # two potential formulations, one with respectively
#' # two and four parameters
#' # formulations for variables Xi:
#' allFilt$X1 <- t(as.matrix(c(0,0,0,1,0,0,0,0,0,0)))
#' allFilt$X2 <- t(matrix(c(0,-0.85,0,-0.27,0,0,0,0.46,0,0,
#'                          0,-0.64,0,0,0,0,0,0.43,0,0),
#'                        ncol=2, nrow=10))
#' allFilt$X3 <- t(matrix(c(0, 0.52,  0, -1.22e-05,  0, 0, 0.99, 5.38e-05, 0, 0,
#'                          0, 0.52, 0, 0, 0, 0, 0.99, 0, 0, 0),
#'                        ncol=2, nrow=10))
#' findAllSets(allFilt, nS=c(3), nPmin=1, nPmax=14)
#'
#' @export
findAllSets <- function (allFilt, nS=c(3), nPmin=1, nPmax=14)
{
  nVar <- sum(nS)

  # All sets of combined equations from an ensemble of possible equations
  #
  ncolMax <- 10000
  maxNp <- 1
  propMod <- matrix(NA, nrow = nVar, ncol=ncolMax)
  propMod[1:nVar, 1] <- 1:nVar
  Np <- NULL
  for (i in 1:nVar) {
    block <- paste("propMod[i,2:3] <- dim(allFilt$X", i, ")",sep="")
    eval((parse(text = block)))
    block <- paste("Np <- allFilt$Np", i, "",sep="")
    eval((parse(text = block)))
    propMod[i,4:(length(Np)+3)] <- Np
    maxNp <- max(maxNp,length(Np))
  }
  possEq <- propMod[,1:3]
  propMod <- propMod[,4:(maxNp+3)]

  # when one / when more than one equation?
  mth1 <- possEq[possEq[,2] != 1,1]
  onl1 <-  possEq[possEq[,2] == 1,1]
  #
  # equations combinations
  allSets <- matrix(rep(1:possEq[mth1[1],2]), nrow = possEq[mth1,2], ncol = 1)
  allNp <- propMod[mth1[1],allSets]
  allSets <- allSets[allNp <= nPmax]
  #allSets <- allSets[allNp >= nPmin & allNp <= nPmax] # slvn 16/01/2017
  allNp <- propMod[mth1[1],allSets]
  if (2<=length(mth1)) {
    for (i in 2:length(mth1)) {
      allSetsMemo <- allSets
      allNpMemo <- allNp
      allSets <- cbind(allSets,rep(1,dim(as.matrix(allSets))[1]))
      allNp <- cbind(allNp,rep(propMod[mth1[i],1],dim(as.matrix(allNp))[1]))
      for (j in 2:possEq[mth1[i],2]) {
        if (propMod[mth1[i],j] <= nPmax) {
        #if (propMod[mth1[i],j] >= nPmin & propMod[mth1[i],j] <= nPmax) { # slvn 16/01/2017
          tmp <- cbind(allSetsMemo,rep(j,dim(as.matrix(allSetsMemo))[1]))
          tmpNp <- cbind(allNpMemo,rep(propMod[mth1[i],j],dim(as.matrix(allNpMemo))[1]))
          allSets <- rbind(allSets, tmp)
          allNp <- rbind(allNp, tmpNp)
        }
      }
      lesbons <- (colSums(t(allNp)) <= nPmax)
      #lesbons <- (colSums(t(allNp)) >= nPmin & colSums(t(allNp)) <= nPmax) # slvn 16/01/2017
      allSets <- allSets[lesbons,]
      allNp <- allNp[lesbons,]
    }
  }

  # add equations with single models
  SetsNp <- list()
  SetsNp$Np <- SetsNp$Sets <- matrix(NA, nrow=dim(as.matrix(allSets))[1], ncol=nVar)
  SetsNp$Sets[,mth1] <- allSets
  SetsNp$Np[,mth1] <- allNp
  SetsNp$Sets[,onl1] <- 1
  SetsNp$Np[,onl1] <- possEq[onl1,2]
  # check size
  lesbons <- (colSums(t(SetsNp$Np)) >= nPmin & colSums(t(SetsNp$Np)) <= nPmax)
  SetsNp$Sets <- SetsNp$Sets[lesbons, ]
  SetsNp$Np <- SetsNp$Np[lesbons, ]

  if (is.vector(SetsNp$Sets)) {
    SetsNp$Sets = t(SetsNp$Sets)
    SetsNp$Np = t(SetsNp$Np)
  }
  SetsNp
}



