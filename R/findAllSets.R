#' @title Find all possible sets of equation combinations
#' considering an ensemble of possible equation.
#'
#' @description For each equation to be retrieved, an ensemble
#' of potential formulation is given. For instance, if three possible
#' formulations are provided for equation (1), one for equation (2)
#' and two for equation (3). In this case, six (i.e. 3*1*2) possible
#' sets of equations can be obtained from these potential formulations.
#' The aim of this program is to formulate all the potential
#' systems from the individual formulations provided of the
#' individual equations.
#'
#' @inheritParams gPoMo
#'
#' @param allFilt A list with:
#' (1) A matrix \code{allFilt$Xi} of possible formulations
#' for each equation (corresponding to variable \code{Xi});
#' And (2) a vector \code{allFilt$Npi} providing the number of
#' polynomial terms contained in each formulation.
#'
#' @author Sylvain Mangiarotti
#'
#' @examples
#' #############
#' # Example 1 #
#' #############
#' # We build an example
#' allFilt <- list()
#' # For equation 1 (variable X1)
#' allFilt$Np1 <- 1         # only one formulation with one single parameter
#' # For equation 2 (variable X2)
#' allFilt$Np2 <- c(3, 2)   # two potential formulations, with respectively three and four parameters
#' # For equation 3 (variable X3)
#' allFilt$Np3 <- c(4, 2)   # two potential formulations, with respectively two and four parameters
#' # Formulations for variables Xi:
#' # For X1:
#' allFilt$X1 <- t(as.matrix(c(0,0,0,1,0,0,0,0,0,0)))
#' # For X2:
#' allFilt$X2 <- t(matrix(c(0,-0.85,0,-0.27,0,0,0,0.46,0,0,
#'                          0,-0.64,0,0,0,0,0,0.43,0,0),
#'                        ncol=2, nrow=10))
#' # For X3:
#' allFilt$X3 <- t(matrix(c(0, 0.52,  0, -1.22e-05,  0, 0, 0.99, 5.38e-05, 0, 0,
#'                          0, 0.52, 0, 0, 0, 0, 0.99, 0, 0, 0),
#'                        ncol=2, nrow=10))
#' # From these individual we can retrieve all possible formulations
#' findAllSets(allFilt, nS=c(3), nPmin=1, nPmax=14)
#' # if only formulations with seven maximum number of terms are expected:
#' findAllSets(allFilt, nS=c(3), nPmin=1, nPmax=7)
#'
#' @seealso \code{\link{autoGPoMoSearch}}
#'
#' @export
#' 
#' @return SetsNp A list of two matrices
#' $Sets A matrix defining all the sets the equation combination
#' (each line provides a combination, for instance, a line with 1,2,2
#' means the first equation of allFilt$X1, the second one of allFilt$X2
#' and the second one of allFilt$X3)
#' $Np  A matrix providing the number of parameters of all equation
#' combination (each line provides the number of parameter of the selected
#' equations)
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
    block <- paste("propMod[i,2:3] <- dim(as.matrix(allFilt$X", i, "))",sep="")
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



