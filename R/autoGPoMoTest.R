#' @title Tests the numerical integrability of models and
#' classify their dynamical regime
#'
#' @description Tests the numerical integrability of
#' provided models (these may have been obtained with
#' function \code{autoGPoMoSearch}),
#' and classify these models as Divergent, Fixed Points,
#' Periodic or not Unclassified (potentially chaotic).
#'
#' @inheritParams  gPoMo
#' @inheritParams  gloMoId
#' @inheritParams  drvSucc
#' @inheritParams  poLabs
#' @inheritParams  autoGPoMoSearch
#'
#' @importFrom graphics par plot text lines
#' @importFrom stats lsfit
#'
#' @param allKL A list of all the models \code{$mToTest1},
#' \code{$mToTest2}, etc. to be tested. Each model is provided
#' as a matrix.
#' @param numValidIC Line number of the first valid initial
#' conditions, that is, such as weight is not equal to zero.
#' @param IstepMin The minimum number of integration step to start
#' of the analysis (by default \code{IstepMin = 10}).
#' @param IstepMax The maximum number of integration steps for
#' stopping the analysis (by default \code{IstepMax = 10000}).
#' @param tooFarThr Divergence threshold, maximum value
#' of the model trajectory compared to the data standard
#' deviation. By default a trjactory is too far if
#' the distance to the center is larger than four times the variance
#' of the input data.
#' @param LimCyclThr Threshold used to detect the limit cycle.
#' @param FxPtThr Threshold used to detect fixed points.
#' @param method The integration technique used for the numerical
#' integration. By default, the fourth-order Runge-Kutta method
#' (\code{method = 'rk4'}) is used. Other methods such as 'ode45'
#' or 'lsoda' may also be chosen. See package \code{deSolve}
#' for details.
#'
#' @author Sylvain Mangiarotti, Flavie Le Jean
#'
#' @return A list containing:
#' @return \code{$okMod}      A vector classifying the models: diverging models (0),
#' periodic models of period-1 (-1), unclassified models (1).
#' @return \code{$okMod}      A matrix classifying the model variables: diverging variable (0),
#' period-1 variable (-1), period-2 variable (-2), fixed point variable (2), unclassified models (1).
#' @return \code{$coeff}      A matrix with the coefficients of one selected model
#' @return \code{$models}     A list of all the models to be tested \code{$mToTest1},
#' \code{$mToTest2}, etc. and of all selected models \code{$model1}, \code{$model2}, etc.
#' @return \code{$tout}       The time vector of the output time series (vector length
#' corresponding to the longest numerical integration duration)
#' @return \code{$stockoutreg} A list of matrices with the integrated trajectories
#' (variable \code{X1} in column 1, \code{X2} in 2, etc.) for all the models
#' \code{$model1}, \code{$model2}, etc.
#'
#' @seealso \code{\link{autoGPoMoSearch}}, \code{\link{gPoMo}}, \code{\link{poLabs}}
#'
#' @examples
#' #Example
#' # Load data:
#' data('RosYco')
#' # Structure choice
#' data('allToTest')
#' # Test the models
#' outGPT <- autoGPoMoTest(RosYco, nVar= 3, dMax = 2, dt = 1/125, show=1,
#'                         allKL = allToTest, IstepMax = 60)
#'
#' @export
autoGPoMoTest <- function (data, tin = NULL, dt=NULL,
                           nVar = nVar, dMax = dMax,
                           show = 1, verbose = 1,
                           allKL = allKL, numValidIC = 1, weight = NULL, IstepMin = 10,
                           IstepMax = 10000,
                           tooFarThr = 4, FxPtThr = 1E-8, LimCyclThr = 1E-6,
                           method = 'rk4')
{
  # PREPARE AND TEST THE INPUTS
  data <- as.matrix(data)
  if (IstepMax < IstepMin) {
    stop("Integration steps are inconsistent: IstepMax < IstepMin")
  }
  if (is.null(tin) & is.null(dt)) {
    cat("when neither input time vector 'tin' nor time step 'dt' are given")
    cat("'dtFixe' is taken such as dt=0.01")
    dt = 0.01
    tin = 0:(dim(data)[1]-1)*dt
  }
  else if (is.null(tin) & !is.null(dt)) {
    tin = 0:(dim(data)[1]-1)*dt
  }
  else if (!is.null(tin) & is.null(dt)) {
    # if time step is regular
    if (max(abs(diff(tin))) != 0) dt = tin[2]-tin[1]
  }
  else if (!is.null(tin) & !is.null(dt)) {
    cat("input time vector 'tin' and  time step 'dt' are inconsistent")
    stop
  }
  if (is.null(weight)) {
    weight <- tin * 0 + 1
  }
  tout <- NULL

  # BASIC STATISTICS OF THE INPUT DATA
  datamin <- apply(data[, 1:nVar], 2, min)
  datamax <- apply(data[, 1:nVar], 2, max)
  datamoy <- apply(data[, 1:nVar], 2, mean)
  datasd  <- apply(data[, 1:nVar], 2, sd)
  #
  #    derivdata <- reg$init
  #    if (is.null(tooFarThr)) {
  #        tooFarThr <- 4 * datasd
  #    }
  #    if (is.null(LimCyclThreshold)) {
  #        LimCyclThreshold <- datasd[1] * 0.1
  #    }
  #    if (is.null(FxPtThr)) {
  #        FxPtThr <- datasd[1] * 0.01
  #    }
  #
  # Compute pMx
  pMax <- d2pMax(nVar, dMax)
  # Not sure this is really necessary:
  listMod <- 1:length(allKL)
  # by default all the input models are unclassified (=1)
  okMod <- rep(1,length(allKL))
  # by default all the variables are unclassified (=1)
  okVar <- matrix(1, ncol = length(allKL), nrow = nVar)
  # Open a new window for plotting the figures
  if (show == 1) dev.new()
  #op <- par(mfrow = c(4, 6), pty = "s", mar=c(5,3,2,2)+0.1)
  op <- par(mfrow = c(4, 6), mar=c(5,3,2,2)+0.1)
  #
  # Initiate the integration time step number
  Istep <- IstepMin
  # Initiate the memory of the models initial conditions
  StockInitState <- list()
  # Get the initial conditions from the original data matrix
  InitStates <- matrix(data=1, nrow = length(listMod), ncol = 1) %*%
    as.vector(data[numValidIC, 1:nVar])
  # reference time
  ptm <- NULL
  #kmod <- matrix(NA, nrow = pMax*0, ncol = nVar)
  # model matrix creation
  KL <- NULL
  
  
  #### LOOP FOR NUMERICAL INGRATION ON INCREASING DURATION LENGTHS (BEGINING)
  while (sum(okMod == 1) != 0 & Istep <= IstepMax) {
    
    # Display information preparation (begining)
    if (verbose == 1) {
      if (is.null(ptm)) {
        # Number of model under test
        block <- paste("### For Istep = ", Istep,
                       " (max: ", IstepMax,
                       "), models to test: ", sum(okMod == 1),
                       " / ", length(listMod), sep="")
      }
      else {
        # Number of model under test
        Tnext <- (proc.time()[3] - ptm) / sum(oldOkMod == 1) * sum(okMod == 1) * 2
        hr <- Tnext %/% 3600
        min <- floor((Tnext - hr * 3600) %/% 60)
        sec <- round(Tnext  - hr * 3600 - min * 60, digits = 2)
        block <- paste("### For Istep = ", Istep,
                       " (max: ", IstepMax,
                       "), models to test: ", sum(okMod == 1),
                       " / ", length(listMod),
                       ". Runtime: ~ ",
                       hr, "h ",
                       min, "min ",
                       sec, "s ",
                       sep="")
      }
      cat(block, "\n")
      ptm <- proc.time()[3]
      oldOkMod <- okMod
    }
    # Display information preparation (end)
    
    ### PREPARE THE LOOP ON ALL THE MODELS
    #
    # initiate what models have to be considered in the loop: i.e. such as okMod == 1
    whatModel <- (okMod == 1) * rep(1:length(okMod == 1))
    whatModel <- whatModel[whatModel != 0]
    # initiate the model matrix and set it to NA
    kmod <- matrix(NA, nrow = pMax*0, ncol = nVar)
    #
    
    ### LOOP ON ALL THE MODELS (BEGINING)
    for (i in 1:length(whatModel)) {
      
      # Current model number
      iMod <- whatModel[i]
      #
      # Load the model
      block <- paste("KL <- allKL$mToTest", iMod, sep="")
      eval((parse(text = block)))
      #
      # Try integration and integrate
      outreg <- try(deSolve::ode(InitStates[iMod, ], (0:Istep) *
                          dt, derivODE2, KL, verbose = 0,
                          method = method), silent = TRUE)
      options(warn = 0)
      #
      # Detect unclassified variables
      whatVar <- which(okVar[,iMod] == 1)
      #
      # How many unclassified variables
      nb <- length(whatVar)
      
      # If at least 2 unclassified variables then plot
      # (one part of the model remains unclassified)
      if (show >= 1 & nb >= 2) {
        # Model phase portrait
        plot(outreg[, 1+whatVar[1]], outreg[, 1+whatVar[2]], type = "l",
             main = paste("#", iMod, "(", sum(KL!=0
             ), "p)"), col = "red",
             xlab=paste("I=",Istep," (X", whatVar[1],",X",whatVar[2],")", sep = ""),
             ylab="")
        # Add the model initial conditions on the plot
        lines(outreg[1, 1+whatVar[1]], outreg[1, 1+whatVar[2]], type = "p",
              col = "red")
        # Add the original phase portrait for comparison
        data0 <- data
        data0[weight == 0,2] <- NaN
        lines(data0[, whatVar[1]], data0[, whatVar[2]], type = "l", col = 'gray')
        
        # If more than 3 variables, then plot another projection (begining)
        if (nVar > 3) {
          # Model phase portrait
          plot(outreg[, 1+whatVar[nb-1]], outreg[, 1+whatVar[nb]], type = "l",
               main = paste("#", iMod, "(", sum(KL!=0
               ), "p)"), col = "green",
               xlab=paste("I=",Istep," (X", whatVar[nb-1],",X",whatVar[nb],")", sep = ""),
               ylab="")
          # Put the model initial conditions
          lines(outreg[1, 1+whatVar[nb-1]], outreg[1, 1+whatVar[nb]], type = "p",
                col = "green")
          # Plot the original phase portrait
          data0 <- data
          data0[weight==0,2] <- NaN
          lines(data0[, whatVar[nb-1]], data0[, whatVar[nb]], type = "l", col = 'gray')
        }
        # If more than 3 variables, then plot another projection (end)
        
      }
      # If at least 2 unclassified variables then plot (end)
      
      
      # INFORMATION ABOUT THE RESULTS
      # (Initial part: prepare the place on the figure)
      if (show >= 1) {
        plot(0, 0, type = "p", axes = 0,
             xlab = "", ylab = "",
             main = paste("#", iMod, "(", sum(KL!=0), "p)"),
             col = "white")
      }
      
      # DETECT if integration failed (begining)
      if (inherits(outreg, "try-error")) {
        # the model is rejected
        okMod[iMod] <- 0
        # all the variables are rejected
        okVar[,iMod] <- 0
        # Failing is indicated
        # INFORMATION ABOUT THE RESULTS
        if (show >= 1) {
          #plot(0, 0, type = "p",
          #     main = paste("#", iMod, "(", sum(KL!=0
          #     ), "p)"), col = "white", xlab=paste("I=",Istep), ylab="")
          text(0, 0.9, "integration failed")
        }
      }
      # DETECT if integration failed (end)
      
      # ALL OTHER CASES:
      # Too far or Divergent (0)
      # Stable Fixed Point (2)
      # Period-1 limit cycles (-1)
      # Period-2 limit cycle (-2) etc.
      #
      else {
        # How many time steps
        nlast <- dim(outreg)[1]
        # Distance relatively to the standard deviation (for each variable)
        ratio <- (outreg[nlast, 2:(nVar + 1)] - datamoy) / datasd
        # Test for variable rejection
        iznogood <- ratio[abs(ratio) > tooFarThr]

        # 0: DIVERGENT or TOO FAR
        if (is.na(sum(outreg)) | sum(iznogood) != 0) {
          # Variable identification:
          okVar[abs(ratio) > tooFarThr,iMod] <- 0
          okVar[is.na(sum(outreg)),iMod] <- 0
        }
        # Model identification (if all the variables are too far):
        if (sum(okVar[,iMod] == 0) == nVar) okMod[iMod] <- 0
        
        # 2: STABLE FIXED POINT
        isSFP <- abs(outreg[nlast, 2:(nVar+1)] - outreg[1,2:(nVar+1)]) / datasd
        if (length(which(isSFP <= FxPtThr) >= 1)) {
          # Variable identification:
          whatFPt <- which(abs(outreg[nlast,2:(nVar+1)]-outreg[1,2:(nVar+1)])/datasd <= FxPtThr)
          okVar[whatFPt,iMod] <- 2
        }
        # Model identification (if all the variables Stable Fixed Point):
        if (sum(okVar[,iMod] == 2) == nVar) okMod[iMod] <- 2
        
        # -1: LIMIT CYCLES
          izP1 <- list()
          izP1$status <- 0
          izP1err <- rep(100,nVar)
          #izP1 <- detectP1limCycl(outreg[,2:3], LimCyclThreshold = LimCyclThreshold, show = show)
          # The behavior is considered through a two-by-two variable study
          for (k in 1:(nVar-1)) {
            for (l in (k+1):nVar) {
              # test the periodicity:
              outP <- testP(data = outreg[,c(k,l)+1] / datasd[c(k,l)],
                            wthresh = LimCyclThr,
                            fxPtThresh = FxPtThr,
                            show = 0)
              # memo
              izP1err[k] <- min(izP1err[k], outP$errP1)
              izP1err[l] <- min(izP1err[l], outP$errP1)
              # variable identification:
              if (outP$status == -1) {
                okVar[c(k,l),iMod] <- -1
              }
            }
          }
          # If all the variables are P1 then the model is classified as P1
          if (sum(okVar[,iMod] == -1 & okMod[iMod] != 0) == nVar) izP1$status <- -1
          
          # If only 2 or less variables are unclassified
          if (sum(okVar[,iMod] == 1) <= 2 & okMod[iMod] != 0) {
            
            # If at least two are classified as P2, the model is classified P2
            if (sum(okVar[,iMod] == -2) >= 2) {
              okMod[iMod] <- -2
            }
            # If at least two are classified as P1, the model is classified P1
            else if (sum(okVar[,iMod] == -1) >= 2) {
              okMod[iMod] <- -1
            }
            # If at least two are classified as FixPt, the model is classified FixPt
            else if (sum(okVar[,iMod] == 2) >= 2) {
              okMod[iMod] <- 2
            }
            # Else model remains unclassified
            else if (okMod[iMod] != 0) {
              okMod[iMod] <- 1
            }
          }

          # If the model is unclassified
          if (okMod[iMod] == 1) {
            kmod <- rbind(kmod, KL)
            #
            eval(parse(text = paste("allKL$model",
                                    iMod, "<- KL", sep = "")))
            #
            tout <- outreg[, 1]
          }

        # INFORMATION ABOUT THE RESULTS (begining)
        if (show >= 1) {
          # Header
          text(-1.0, 0.97, paste("Dist:"), col = 'blue', adj = c(0,0))
          text(-0.2, 0.97, paste("FixPt:"), col = 'blue', adj = c(0.5,0))
          text(0.3, 0.97, paste("P1:"), col = 'blue', adj = c(0.5,0))
          text(1.05, 0.97, paste("what:"), col = 'blue', adj = c(1,0))
          # Conditions
          text(-1.0, 0.87, paste("<",tooFarThr), col = 'green', adj = c(0,0))
          text(-0.2, 0.87, paste(">",FxPtThr), col = 'green', adj = c(0.5,0))
          text(0.3, 0.87, paste(">",LimCyclThr), col = 'green', adj = c(0.5,0))
          #
          lines(c(-1,1), c(0.85,0.85), col = 'blue')

          # TOO FAR?
          for (i in 1:nVar) {
            # divergent
            if (is.nan(ratio[i])) {
              block <- paste("text(-1.0,", 0.9 - 0.2 * i,
                             ", paste(NaN), col = 'red', adj = c(0,0))", sep="")
            }
            # not too far
            else if (abs(ratio[i]) <= 0.8 * tooFarThr) {
              block <- paste("text(-1.0,", 0.9 - 0.2 * i,
                             ", paste(signif(",abs(ratio[i]),
                             ",digits=1)), col = 'green', adj = c(0,0))", sep="")
            }
            # almost too far
            else if (abs(ratio[i]) > 0.8 * tooFarThr
                     & abs(ratio[i]) <= tooFarThr) {
              block <- paste("text(-1.0,", 0.9 - 0.2 * i,
                             ", paste(signif(",ratio[i],
                             ",digits=1)), col = 'orange', adj = c(0,0))", sep="")
            }
            # too far
            else {
              block <- paste("text(-1.0,", 0.9 - 0.2 * i,
                             ", paste(signif(",ratio[i],
                             ",digits=1)), col = 'red', adj = c(0,0))", sep="")
            }
            # display the text
            eval((parse(text = block)))
            
            # Stable Fixed Point (SFP) ?
              # SFP not detected
              if (isSFP[i] >= 2 * FxPtThr
                  & is.nan(ratio[i]) == 'FALSE') {
                block <- paste("text(-0.2,", 0.9 - 0.2 * i,
                               ", paste(signif(",isSFP[i],
                               ",digits=2)), col = 'green', adj = c(0.5,0))", sep="")
              }
              # SFP almost detected
              else if (isSFP[i] < 2 * FxPtThr
                       & isSFP[i] > FxPtThr
                       & is.nan(ratio[i]) == 'FALSE') {
                block <- paste("text(-0.2,", 0.9 - 0.2 * i,
                               ", paste(signif(",isSFP[i],
                               ",digits=2)), col = 'orange', adj = c(0.5,0))", sep="")
              }
              # SFP detected
              else if (isSFP[i] <= FxPtThr
                       & is.nan(ratio[i]) == 'FALSE') {
                block <- paste("text(-0.2,", 0.9 - 0.2 * i,
                               ", paste(signif(",isSFP[i],
                               ",digits=2)), col = 'red', adj = c(0.5,0))", sep="")
              }
              # display the text
              eval((parse(text = block)))
              
              # Limit Cycle ?
              # Periodic cycle Not detected
              if ((izP1err[i] >= LimCyclThr) == 'TRUE') {
                izP1err[i]
                # Too short data set
                if ((izP1err[i] == 100) == 'TRUE') {
                  block <- paste("text(0.4,", 0.9 - 0.2 * i,
                                 ", '-', col = 'red', adj = c(0.5,0))", sep="")
                }
                else {
                  block <- paste("text(0.4,", 0.9 - 0.2 * i,
                               ", paste(signif(",izP1err[i],
                               ",digits=2)), col = 'green', adj = c(0.5,0))", sep="")
                }
              }
              # Periodic cycle Detected
              else if ((izP1err[i] < LimCyclThr) == 'TRUE') {
                izP1err[i]
                  block <- paste("text(0.4,", 0.9 - 0.2 * i,
                                 ", paste(signif(",izP1err[i],
                                 ",digits=2)), col = 'red', adj = c(0.5,0))", sep="")
              }
              eval((parse(text = block)))
              
              # BILAN
              if (okVar[i,iMod] == 1) {
                block <- paste("text(1.05,", 0.9 - 0.2 * i,
                               ",", okVar[i,iMod],
                               ", col = 'green', adj = c(1,0))", sep="")
              }
              else {
                block <- paste("text(1.05,", 0.9 - 0.2 * i,
                               ",", okVar[i,iMod],
                               ", col = 'red', adj = c(1,0))", sep="")
              }
              eval((parse(text = block)))
              
          }
          
          ############################
          # OkMod
          if (okMod[iMod] == 1) {
            block <- paste("text(0.9,", 0.7 - 0.2 * i,
                           ",", okMod[iMod],
                           ", col = 'green', adj = c(1,0))", sep="")
            lines(0.81, 0.77 - 0.2 * i,
                  type = 'p', pch = 22, col = 'green', cex = 4)
          }
          else {
            block <- paste("text(0.9,", 0.7 - 0.2 * i,
                           ",", okMod[iMod],
                           ", col = 'red', adj = c(1,0))", sep="")
            lines(0.81, 0.77 - 0.2 * i,
                  type = 'p', pch = 25, col = 'red', cex = 5)
          }
          eval((parse(text = block)))
          
        }
        # INFORMATION ABOUT THE RESULTS (end)
        
        
        if (sum(okVar[,iMod] == 1) == 0) okMod[iMod] <- 0
        InitStates[iMod, ] <- outreg[nlast, 2:(nVar +1)]
        eval(parse(text = paste("StockInitState$model",
                                iMod, "<- outreg[, 1:(nVar + 1)]", sep = "")))
      }
    }
    ### LOOP ON ALL THE MODELS (END)
    
    Istep <- 2 * Istep
    if (show == 1) {
      if (nVar <= 3) nsubplot <- 2*sum(okMod == 1)
      else nsubplot <- 3* sum(okMod == 1)
      if (nsubplot > 16) {
        #op <- par(mfrow = c(4, 6), pty = "s")
        op <- par(mfrow = c(4, 6))
      }
      else {
        if (nsubplot > 0) {
          #op <- par(mfrow = c(ceiling(sqrt(nsubplot)),
          #                    ceiling(sqrt(nsubplot))), pty = "s")
          op <- par(mfrow = c(ceiling(sqrt(nsubplot)),
                              ceiling(sqrt(nsubplot))))
        }
      }
    }
  }
  #### LOOP FOR NUMERICAL INGRATION ON INCREASING DURATION LENGTHS (END)
  
  
  if (verbose == 1) {
    # Number of unclassified models
    block <- paste("### Number of unclassified models: ", sum(okMod == 1),
                   " / ", length(listMod), sep="")
    cat(block, "\n")
  }
  # return
  list(okMod = okMod, okVar = okVar, tout = tout,
       stockoutreg = StockInitState, coeff = kmod,
       models = allKL)
}
