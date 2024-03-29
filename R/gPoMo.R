#' @title Generalized Polynomial Modeling
#'
#' @seealso \code{\link{gloMoId}}, \code{\link{autoGPoMoSearch}},
#'          \code{\link{autoGPoMoTest}}
#'
#' @description Algorithm for a Generalized Polynomial formulation
#' of multivariate Global Modeling.
#' Global modeling aims to obtain multidimensional models
#' from single time series [1-2].
#' In the generalized (polynomial) formulation provided in this
#' function, it can also be applied to multivariate time series [3-4].
#'
#' Example:\cr
#' Note that \code{nS} provides the number of dimensions used from each variable
#'
#' case I \cr
#' For \code{nS=c(2,3)} means that 2 dimensions are reconstructed from variable 1:
#' the original variable \code{X1} and its first derivative \code{X2}), and
#' 3 dimensions are reconstructed from variable 2: the original
#' variable \code{X3} and its first and second derivatives \code{X4} and \code{X5}.
#' The generalized model will thus be such as: \cr
#' \eqn{dX1/dt = X2}\cr
#' \eqn{dX2/dt = P1(X1,X2,X3,X4,X5)}\cr
#' \eqn{dX3/dt = X4}\cr
#' \eqn{dX4/dt = X5}\cr
#' \eqn{dX5/dt = P2(X1,X2,X3,X4,X5).}\cr
#'
#' case II\cr
#' For \code{nS=c(1,1,1,1)} means that only the original variables
#' \code{X1}, \code{X2}, \code{X3} and \code{X4} will be used.
#' The generalized model will thus be such as:\cr
#' \eqn{dX1/dt = P1(X1,X2,X3,X4)}\cr
#' \eqn{dX2/dt = P2(X1,X2,X3,X4)}\cr
#' \eqn{dX3/dt = P3(X1,X2,X3,X4)}\cr
#' \eqn{dX4/dt = P4(X1,X2,X3,X4).}\cr
#'
#' @importFrom utils tail
#'
#'@examples
#' #Example 1
#' data("Ross76")
#' tin <- Ross76[,1]
#' data <- Ross76[,3]
#' dev.new()
#' out1 <- gPoMo(data, tin = tin, dMax = 2, nS=c(3), show = 1,
#'               IstepMax = 1000, nPmin = 9, nPmax = 11)
#' visuEq(out1$models$model1, approx = 4)
#' 
#' 
#'\donttest{
#' #Example 2
#' data("Ross76")
#' tin <- Ross76[,1]
#' data <- Ross76[,3]
#' # if some data are not valid (vector 'weight' with zeros)
#' W <- tin * 0 + 1
#' W[1:100] <- 0
#' W[700:1500] <- 0
#' W[2000:2800] <- 0
#' W[3000:3500] <- 0
#' dev.new()
#' out2 <- gPoMo(data, tin = tin, weight = W,
#'                  dMax = 2, nS=c(3), show = 1,
#'                  IstepMax = 6000, nPmin = 9, nPmax = 11)
#' visuEq(out2$models$model3, approx = 4)
#'}
#'
#'
#'\donttest{
#' #Example 3
#' data("Ross76")
#' tin <- Ross76[,1]
#' data <- Ross76[,2:4]
#' dev.new()
#' out3 <- gPoMo(data, tin=tin, dMax = 2, nS=c(1,1,1), show = 1,
#'               IstepMin = 10, IstepMax = 3000, nPmin = 7, nPmax = 8)
#' # the simplest model able to reproduce the observed dynamics is model #5
#' visuEq(out3$models$model5, approx = 3, substit = 1) # the original Rossler system is thus retrieved
#'}
#'
#'\donttest{
#' #Example 4
#' data("Ross76")
#' tin <- Ross76[,1]
#' data <- Ross76[,2:3]
#' # model template:
#' EqS <- matrix(1, ncol = 3, nrow = 10)
#' EqS[,1] <- c(0,0,0,1,0,0,0,0,0,0)
#' EqS[,2] <- c(1,1,0,1,0,1,1,1,1,1)
#' EqS[,3] <- c(0,1,0,0,0,0,1,1,0,0)
#' visuEq(EqS, substit = c('X','Y','Z'))
#' dev.new()
#' out4 <- gPoMo(data, tin=tin, dMax = 2, nS=c(2,1), show = 1,
#'       EqS = EqS, IstepMin = 10, IstepMax = 2000,
#'       nPmin = 9, nPmax = 11)
#' visuEq(out4$models$model2, approx = 2, substit = c("Y","Y2","Z"))
#'}
#'
#'\donttest{
#' #Example 5
#' # load data
#' data("TSallMod_nVar3_dMax2")
#' #multiple (six) time series
#' tin <- TSallMod_nVar3_dMax2$SprK$reconstr[1:400,1]
#' TSRo76 <- TSallMod_nVar3_dMax2$R76$reconstr[,2:4]
#' TSSprK <- TSallMod_nVar3_dMax2$SprK$reconstr[,2:4]
#' data <- cbind(TSRo76,TSSprK)[1:400,]
#' dev.new()
#' # generalized Polynomial modelling
#' out5 <- gPoMo(data, tin = tin, dMax = 2, nS = c(1,1,1,1,1,1),
#'               show = 0, method = 'rk4',
#'               IstepMin = 2, IstepMax = 3,
#'               nPmin = 13, nPmax = 13)
#'
#' # the original Rossler (variables x, y and z) and Sprott (variables u, v and w)
#' # systems are retrieved:
#' visuEq(out5$models$model347, approx = 4,
#'        substit = c('x', 'y', 'z', 'u', 'v', 'w'))
#' # to check the robustness of the model, the integration duration
#' # should be chosen longer (at least IstepMax = 4000)
#'}
#'
#' @author Sylvain Mangiarotti, Flavie Le Jean, Mireille Huc
#'
#' @references
#' [1] Gouesbet G. & Letellier C., 1994. Global vector-field reconstruction by using a
#' multivariate polynomial L2 approximation on nets, Physical Review E, 49 (6),
#' 4955-4972. \cr
#' [2] Mangiarotti S., Coudret R., Drapeau L. & Jarlan L., Polynomial search and
#' Global modelling: two algorithms for modeling chaos. Physical Review E, 86(4),
#' 046205. \cr
#' [3] Mangiarotti S., Le Jean F., Huc M. & Letellier C., Global Modeling of aggregated
#' and associated chaotic dynamics. Chaos, Solitons and Fractals, 83, 82-96. \cr
#' [4] S. Mangiarotti, M. Peyre & M. Huc, 2016.
#' A chaotic model for the epidemic of Ebola virus disease
#' in West Africa (2013-2016). Chaos, 26, 113112. \cr
#'
#' @inheritParams  gloMoId
#' @inheritParams  drvSucc
#' @inheritParams  poLabs
#' @inheritParams  autoGPoMoSearch
#' @inheritParams  autoGPoMoTest
#'
#' @param dtFixe Time step used for the analysis. It should correspond
#' to the sampling time of the input data.
#' Note that for very large and very small time steps, alternative units
#' may be used in order to stabilize the numerical computation.
#' @param nS A vector providing the number of dimensions used for each
#' input variables (see Examples 1 and 2). The dimension of the resulting
#' model will be \code{nVar = sum(nS)}.
#' @param EqS Model template including all allowed regressors.
#' Each column corresponds to one equation. Each line corresponds to one
#' polynomial term as defined by function \code{poLabs}.
#' @param AndManda AND-mandatory terms in the equations (all the provided
#' terms should be in the equations).
#' @param OrMandaPerEq OR-mandatory terms per equations (at least one of
#' the provided terms should be in each equation).
#' @param nPmin Corresponds to the minimum number of parameters (and thus
#' of polynomial term) allowed.
#' @param nPmax Corresponds to the maximum number of parameters (and thus
#' of polynomial) allowed.
#' @param nPminPerEq Corresponds to the minimum number of parameters (and thus
#' of polynomial term) allowed per equation.
#' @param nPmaxPerEq Corresponds to the maximum number of parameters (and thus
#' of polynomial) allowed per equation.
#' @param verbose Gives information (if set to 1) about the algorithm
#' progress and keeps silent if set to 0.
#'
#' @return A list containing:
#' @return \code{$tin}        The time vector of the input time series
#' @return \code{$inputdata}  The input time series
#' @return \code{$tfiltdata}  The time vector of the filtered time series (boudary removed)
#' @return \code{$filtdata}   A matrix of the filtered time series with its derivatives
#' @return \code{$okMod}      A vector classifying the models: diverging models (0), periodic
#' models of period-1 (-1), unclassified models (1).
#' @return \code{$coeff}      A matrix with the coefficients of one selected model
#' @return \code{$models}     A list of all the models to be tested \code{$mToTest1},
#' \code{$mToTest2}, etc. and all selected models \code{$model1}, \code{$model2}, etc.
#' @return \code{$tout}       The time vector of the output time series (vector length
#' corresponding to the longest numerical integration duration)
#' @return \code{$stockoutreg} A list of matrices with the integrated trajectories
#' (variable \code{X1} in column 1, \code{X2} in 2, etc.) of all the models \code{$model1},
#' \code{$model2}, etc.
#'
#' @import deSolve
#' @export
#'
#' @seealso \code{\link{autoGPoMoSearch}}, \code{\link{autoGPoMoTest}}, \code{\link{visuOutGP}},
#'          \code{\link{poLabs}}, \code{\link{predictab}}, \code{\link{drvSucc}}
#'
gPoMo <- function (data, tin = NULL, dtFixe = NULL, dMax = 2, dMin = 0,
                   nS=c(3), winL = 9,
                   weight = NULL, show = 1, verbose = 1,
                   underSamp = NULL,
                   EqS = NULL,
                   AndManda = NULL, OrMandaPerEq = NULL,
                   IstepMin = 2, IstepMax = 2000,
                   nPmin = 1, nPmax = 14,
                   tooFarThr = 4, FxPtThr = 1E-8, LimCyclThr = 1E-6,
                   nPminPerEq = 1, nPmaxPerEq = NULL,
                   method = 'rk4')
{

  nVar = sum(nS)
  pMax <- d2pMax(nVar, dMax, dMin=dMin)

  if (is.vector(data)) data <- as.matrix(data)

  if (IstepMax < IstepMin) {
    stop("Integration steps are inconsistent: IstepMax < IstepMin")
  }
  
  if (is.null(tin) & is.null(dtFixe)) {
    message("when neither input time vector 'tin' nor time step 'dtFixe' are given")
    message("'dtFixe' is taken such as dtFixe=0.01")
    dtFixe = 0.01
    tin = 0:(dim(data)[1]-1)*dtFixe
  }
  else if (is.null(tin) & !is.null(dtFixe)) {
    tin = 0:(dim(data)[1]-1)*dtFixe
  }
  else if (!is.null(tin) & is.null(dtFixe)) {
    if (is.matrix(tin)) {
      if(dim(tin)[1] == 1 | dim(tin)[2] == 1) {
        tin <- as.vector(tin)
      }
      else {
        stop("Input time vector dimension 'tin': (", dim(tin)[1], " ", dim(tin)[2],") incompatible")
      }
    }
    # if time step is regular
    if (max(abs(diff(tin))) != 0) dtFixe = tin[3] - tin[2]
  }
  else if (!is.null(tin) & !is.null(dtFixe)) {
    message("Input time vector 'tin' and  time step 'dtFixe' are inconsistent")
    stop
  }

  if (is.null(weight)) {
    weight <- tin * 0 + 1
  }
  
  if (is.vector(data)) data <- as.matrix(data)
  
  # if nPmaxPerEq is not provided it is chosen equal to nPmax
  # in order to keep all the possible equations
  if (is.null(nPmaxPerEq)) nPmaxPerEq <- nPmax
  # convert nPminPerEq and nPmaxPerEq to vectors if required
  # (the minimum number of term will be the same for all equations)
  if (length(nPminPerEq) == 1) nPminPerEq <- rep(nPminPerEq,nVar)
  if (length(nPmaxPerEq) == 1) nPmaxPerEq <- rep(nPmaxPerEq,nVar)
  
  if(sum((nPmaxPerEq - nPminPerEq) < 0) > 0) {
    stop("stop: nPminPerEq and nPmaxPerEq are inconsistent")
  }

  # Compute the first nS[i] derivatives dXi/dt for each time series Xi
    # output size
    if (max(nS)+1 <= 4) {
      toutref <- tin[((winL+1)/2):(length(tin) - (winL+1)/2)]
    }
    if (max(nS)+1 > 4 & max(nS)+1 <= 8) {
      toutref <- tin[((2*winL+1)/2):(length(tin) - (2*winL+1)/2)]
    }
    # initiate output
    derivedleft <- derivedright <- matrix(0, ncol = 0, nrow = length(toutref))
    # compute derivatives right and left for each time series
    for (i in 1:length(nS)) {
      objse <- data[,i]
      derivserie <- drvSucc(tin, objse, nDeriv=nS[i]+1, weight = weight,
                            tstep = dtFixe, winL = winL)
      derivdata <- derivserie$seriesDeriv
      nfirst <- which(derivserie$tout == toutref[1])
      nlast <- which(derivserie$tout == tail(toutref,1))
      tout <- derivserie$tout[nfirst:nlast] # modified 11/07/2022
      Wout <- derivserie$Wout[nfirst:nlast] # modified 11/07/2022
      derivedright <- cbind(derivedright, derivdata[nfirst:nlast,1:nS[i]])
      derivedleft <- cbind(derivedleft, derivdata[nfirst:nlast,(nS[i]+1)])

    }

  # build model's structure
  #
  # equations will be all defined in allFilt
  allFilt <- list()
  #
  quelco <- vector(length = 0)
  for (i in 1:length(nS)) {
    quelco <- c(quelco,(length(quelco)+1):(length(quelco)+nS[i])+1)
  }
  #
  # equations resultating from canonical form: dXi/dt = Xi+1
  for (i in 1:nVar) {
    filt <- t(as.matrix(rep(0,pMax)))
    nx <- sum((poLabs(nVar,dMax,dMin = dMin) == paste("X", quelco[i], " ", sep="")) * rep(1:pMax))
    
    #nx <- sum((poLabs(nVar,dMax) == paste("X", i, " ", sep="")) * rep(1:pMax))
    #reg = regOrd(nVar, dMax)
    #nx = which((reg[quelco[i],] == 1) &  (colSums(reg) == 1)) 
    
    filt[nx] = 1
    block <- paste("allFilt$X", i, " <- filt", sep="")
    eval((parse(text = block)))
    #
    # equations sizes
    block <- paste("allFilt$Np", i, " <- as.vector(sum(filt))", sep="")
    eval((parse(text = block)))
  }
  #
  # equations for which models should be found:
  notreadymod <- nS[1]
  for (i in 2:length(nS)) notreadymod <- c(notreadymod, sum(nS[1:i]))
  #
  # Search of possible equation dXi/dt = f(X1,X2,..,Xn) for each variable Xi
  for (i in 1:length(nS)) {
    coherentdata <- cbind(derivedright[,], as.matrix(derivedleft)[,i])

    # model prestructure
    if (is.null(EqS)) {
      filterReg <- rep(1,pMax)
    }
    else {
      filterReg <- as.vector(EqS[,sum(nS[1:i])])
    }

    Memo = autoGPoMoSearch(coherentdata, dt=dtFixe, nVar = nVar, dMax = dMax, dMin = dMin, weight = Wout,
                           show = show, underSamp = underSamp, filterReg = filterReg)
    filt <- Memo$filtMemo
    K <- Memo$KMemo
    block <- paste("allFilt$X", sum(nS[1:i]), " <- K", sep="")
    eval((parse(text = block)))
    #
    # equations sizes
    block <- paste("allFilt$Np", sum(nS[1:i]), " <- as.vector(colSums(t(filt)))", sep="")
    eval((parse(text = block)))
  }
  
  # Remove equations individually too small or too big
  for (i in 1:nVar) {
    # Take the information from the possible models provided in input
    block <- paste("PotModParam <- allFilt$X", i, "",sep="")
    eval((parse(text = block)))
    whatlines <- (colSums(t(PotModParam)!=0) >= nPminPerEq[i]
                  & colSums(t(PotModParam)!=0) <= nPmaxPerEq[i])
    PotModParam <- PotModParam[whatlines,]
    if (is.vector(PotModParam)) PotModParam <- t(as.matrix(PotModParam))
    block <- paste("allFilt$X", i, " <- as.matrix(PotModParam)",sep="")
    eval((parse(text = block)))
    block <- paste("allFilt$Np", i, " <- as.matrix(colSums(t(PotModParam)!=0))",sep="")
    eval((parse(text = block)))
    }
  
  # Remove equations that do not include the And-mandatory terms
  if (is.null(AndManda) == 'FALSE') {
    for (i in 1:nVar) {
      # Take the information from the possible models provided in input
      block <- paste("PotModParam <- allFilt$X", i, "",sep="")
      eval((parse(text = block)))
      tmpM <- rep(AndManda[,i],dim(PotModParam)[1])
      dim(tmpM) <- c(dim(PotModParam)[2],dim(PotModParam)[1])
      whatlines <- which(colSums(t(PotModParam!=0) * tmpM) == colSums(tmpM))
      PotModParam <- PotModParam[whatlines,]
      if (is.vector(PotModParam)) PotModParam <- t(as.matrix(PotModParam))
      block <- paste("allFilt$X", i, " <- as.matrix(PotModParam)",sep="")
      eval((parse(text = block)))
      block <- paste("allFilt$Np", i, " <- as.matrix(colSums(t(PotModParam)!=0))",sep="")
      eval((parse(text = block)))
    }
  }
  
  # Remove equations that do not include the OR-mandatory terms
  if (is.null(OrMandaPerEq) == 'FALSE') {
    for (i in 1:nVar) {
      # Take the information from the possible models provided in input
      block <- paste("PotModParam <- allFilt$X", i, "",sep="")
      eval((parse(text = block)))
      tmpM <- rep(OrMandaPerEq[,i],dim(PotModParam)[1])
      dim(tmpM) <- c(dim(PotModParam)[2],dim(PotModParam)[1])
      whatlines <- which(colSums(t(PotModParam!=0) * tmpM) == 1)
      PotModParam <- PotModParam[whatlines,]
      if (is.vector(PotModParam)) PotModParam <- t(as.matrix(PotModParam))
      block <- paste("allFilt$X", i, " <- as.matrix(PotModParam)",sep="")
      eval((parse(text = block)))
      block <- paste("allFilt$Np", i, " <- as.matrix(colSums(t(PotModParam)!=0))",sep="")
      eval((parse(text = block)))
    }
  }
  
  # Sets of combined equations
  #
  SetsNp <- findAllSets(allFilt, nS=nS,
                        nPmin=nPmin, nPmax=nPmax)

  # reorder
  nParam <- colSums(t(SetsNp$Np))
  neworder <- order(nParam)

  # prepare the complete models formulation
  allToTest <- list()
  for (i in 1:dim(SetsNp$Sets)[1]) {
    nm <- SetsNp$Sets[neworder[i],]
    # extract equations
    mToTest <- allFilt$X1[nm[1],]
    for (j in 2:nVar) {
      block <- paste("mToTest <- cbind(mToTest, allFilt$X", j,"[nm[", j,"],])", sep="")
      eval((parse(text = block)))
      #      }
    }
    block <- paste("allToTest$mToTest", i," <- mToTest", sep="")
    eval((parse(text = block)))
  }


  # Select valid initial conditions, i.e. such as W != 0
  numValidIC <- min(which(Wout == 1))
  # test integration
  out <- autoGPoMoTest(as.matrix(coherentdata),
                       tin = NULL, dt=dtFixe, nVar = nVar, dMax = dMax, dMin = dMin,
                       show = show, verbose = verbose,
                       allKL = allToTest,
                       numValidIC = numValidIC,
                       weight = Wout,
                       IstepMin = IstepMin,
                       IstepMax = IstepMax, tooFarThr = tooFarThr,
                       LimCyclThr = LimCyclThr,
                       FxPtThr = FxPtThr, method = method)


  # ouput
  out$inputdata <- data
  out$tin <- tin
  out$weight <- weight
  out$filtdata <- as.matrix(coherentdata)
  out$tfiltdata <- tout
  out$Wfiltdata <- Wout
  out
}

