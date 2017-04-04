#' @title autoGPoMoTest: Test models numerical integrability & classify the attractor
#' reached at the convergence (from chosen initial conditions)
#'
#' @description Test the numerical integrability of models (of polynomial structure)
#' possibly obtained with function autoGPoMoSearch, and classify these models as
#' divergent, Fixed Points, Periodic or not classificable (potentially chaotic).
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
#' @param allKL A list of the all the models \code{$mToTest1},
#' \code{$mToTest2}, etc. to be tested. Each model is provided
#' as a matrix.
#' @param IstepMin Minimum step of integration at the beginning
#' of the analysis (by default IstepMin=10).
#' @param IstepMax Maximum step of integration before stopping the
#' analysis, with the interface this value can be changed during
#' the analysis
#' @param tooFarThreshold Divergence threshold, maximum times number
#' that the model can be greater than the data standard deviation,
#' without being removed from the analysis
#' @param LimCyclThreshold Limit cycle threshold. Minimum neighbors
#' distance between two integration steps, without being removed
#' from the analysis
#' @param fixedPtThreshold Limit cycle threshold. Minimum neighbors
#' distance between two integration steps, without being removed
#' from the analysis
#' @param method The integration technique used for the numerical
#' integration. Default is 'lsoda'. Others such as 'rk4' or 'ode45'
#' may also be used. See package deSolve for details.
#'
#' @author Sylvain Mangiarotti, Flavie Le Jean
#'
#' @return A list containing:
#' @return \code{$okMod}      A vector classifying the models: diverging models (0), periodic
#' models of period-1 (-1), unclassified models (1).
#' @return \code{$coeff}      A matrix with the coefficients of one selected model
#' @return \code{$models}     A list of all the models to be tested \code{$mToTest1},
#' \code{$mToTest2}, etc. and all selected models \code{$model1}, \code{$model2}, etc.
#' @return \code{$tout}       The time vector of the output time series (vector length
#' corresponding to the longest numerical integration duration)
#' @return \code{$stockoutreg} A list of matrices with the integrated trajectories
#' (variable \code{X1} in column 1, \code{X2} in 2, etc.) of all the models \code{$model1}, \code{$model2}, etc.
#'
#' @examples
#' #Examples
#' data('RosYco')
#' # Structure choice
#' data('allToTest')
#' outGPT <- autoGPoMoTest(RosYco, nVar= 3, dMax = 2, dt = 1/125, show=1,
#'                         allKL = allToTest, IstepMax = 60)
#'
#' @export
autoGPoMoTest <- function (data, tin = NULL, dt=NULL,
                           nVar = nVar, dMax = dMax,
                           show = 1, verbose = 1,
                           allKL = allKL, IstepMin = 10,
                           IstepMax = 10000, tooFarThreshold = 4, LimCyclThreshold = 0.0,
                           fixedPtThreshold = 1E-8, method = 'lsoda')
{

  data <- as.matrix(data)
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
  tout <- NULL

  datamin <- apply(data[, 1:nVar], 2, min)
  datamax <- apply(data[, 1:nVar], 2, max)
  datamoy <- apply(data[, 1:nVar], 2, mean)
  datasd <- apply(data[, 1:nVar], 2, sd)
  #    derivdata <- reg$init
  #    if (is.null(tooFarThreshold)) {
  #        tooFarThreshold <- 4 * datasd
  #    }
  #    if (is.null(LimCyclThreshold)) {
  #        LimCyclThreshold <- datasd[1] * 0.1
  #    }
  #    if (is.null(fixedPtThreshold)) {
  #        fixedPtThreshold <- datasd[1] * 0.01
  #    }
  pMax <- d2pMax(nVar, dMax)
  listMod <- 1:length(allKL)
  okMod <- rep(1,length(allKL))
  op <- par(mfrow = c(4, 6), pty = "s", mar=c(5,3,2,2)+0.1)
  Istep <- IstepMin
  StockInitState <- list()
  InitStates <- matrix(data=1, nrow = length(listMod), ncol = 1) %*%
    as.vector(data[1, 1:nVar])
  # reference time
  ptm <- NULL
  #kmod <- matrix(NA, nrow = pMax*0, ncol = nVar)
  KL <- NULL
  while (sum(okMod) != 0 & Istep < IstepMax) {
    if (verbose == 1) {
      if (is.null(ptm)) {
        # Number of model under test
        block <- paste("### For Istep = ", Istep,
                       " (max: ", IstepMax,
                       "), models to test: ", sum(okMod),
                       " / ", length(listMod), sep="")
      }
      else {
        # Number of model under test
        Tnext <- (proc.time()[3] - ptm) / sum(oldOkMod) * sum(okMod) * 2
        hr <- Tnext %/% 3600
        min <- floor((Tnext - hr * 3600) %/% 60)
        sec <- round(Tnext  - hr * 3600 - min * 60, digits = 2)
        block <- paste("### For Istep = ", Istep,
                       " (max: ", IstepMax,
                       "), models to test: ", sum(okMod),
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
    #
    whatModel <- (okMod == 1) * rep(1:length(okMod))
    whatModel <- whatModel[whatModel != 0]
    kmod <- matrix(NA, nrow = pMax*0, ncol = nVar)
    for (i in 1:length(whatModel)) {
      iMod <- whatModel[i]
      block <- paste("KL <- allKL$mToTest", iMod, sep="")
      eval((parse(text = block)))
      outreg <- try(deSolve::ode(InitStates[iMod, ], (0:Istep) *
                          dt, derivODE2, KL, verbose = 0,
                          method = method), silent = TRUE)
      options(warn = 0)
      if (inherits(outreg, "try-error")) {
        okMod[iMod] <- 0
        if (show >= 1) {
          plot(0, 0, type = "p",
               main = paste("#", iMod, "(", sum(KL!=0
               ), "p)"), col = "white", xlab=paste("I=",Istep), ylab="")
          text(0, 0.9, "integration failed")
        }
      }
      else {
        nlast <- dim(outreg)[1]

        ratio <- (outreg[nlast, 2:(nVar + 1)]-datamoy)/datasd
        iznogood <- ratio[abs(ratio)>4]
        if (is.na(sum(outreg)) | sum(iznogood) != 0 |
            abs(outreg[nlast, 2] - outreg[1, 2]) <= fixedPtThreshold) {
          okMod[iMod] <- 0
          if (show >= 1) {
            plot(0, 0, type = "p",
                 main = paste("#", iMod, "(", sum(KL!=0
                 ), "p)"), col = "white", xlab=paste("I=",Istep), ylab="")
            if (is.na(sum(outreg))) {
              text(0, 0.9, "DIVERGENT", col = 'orange')
            }
            else if (sum(iznogood) != 0) {
              text(0., 0.9, paste("TOO FAR: |sig/ref| > 4"), col = 'orange')
              for (i in 1:nVar) {
                if (abs(ratio[i]) <= 4) {
                  block <- paste("text(0,", 0.9 - 0.2 * i,
                                 ", paste('variable ", i,
                                 " :',signif(",ratio[i],
                                 ",digits=2)))", sep="")
                }
                else {
                  block <- paste("text(0,", 0.9 - 0.2 * i,
                                 ", paste('variable ", i,
                                 " :',signif(",ratio[i],
                                 ",digits=2)), col = 'red')", sep="")
                }
                eval((parse(text = block)))
              }
            }
            else if (abs(outreg[nlast, 2] - outreg[1,
                                                   2]) <= fixedPtThreshold) {
              text(0, 0.7, "SFpt")
              text(0, 0.5, paste(abs(outreg[nlast, 2] - outreg[1,2]),
                                 " <= ", fixedPtThreshold))
              okMod[iMod] <- 2
            }
          }
        }
        else {
          izP1 <- detectP1limCycl(outreg[,2:3], LimCyclThreshold = LimCyclThreshold, show = show)
          if (izP1 == 1) {
            okMod[iMod] <- -1
            if (show >= 1) {
              plot(0, 0, type = "p",
                   main = paste("#", iMod, "(", sum(KL!=0
                   ), "p)"), col = "white", xlab=paste("I=",Istep), ylab="")
              text(0, 0.9, "P1cycl:")
              text(0, 0.7, paste("distance < LimCyclThreshold"))
              text(0, 0.5, LimCyclThreshold)
            }
          }
          else {
            kmod <- rbind(kmod, KL)
            #
            eval(parse(text = paste("allKL$model",
                                    iMod, "<- KL", sep = "")))
            #
            tout <- outreg[, 1]
            if (show >= 1) {
              plot(outreg[, 2], outreg[, 3], type = "l",
                   main = paste("#", iMod, "(", sum(KL!=0
                   ), "p)"), col = "red", xlab=paste("I=",Istep), ylab="")
              lines(data[, 1], data[, 2], type = "l")
              if (nVar > 3) {
                plot(outreg[, nVar-1], outreg[, nVar], type = "l",
                     main = paste("#", iMod, "(", sum(KL!=0
                     ), "p)"), col = "green", xlab=paste("I=",Istep), ylab="")
                lines(data[, nVar-2], data[, nVar-1], type = "l")
              }
            }
          }
        }
        InitStates[iMod, ] <- outreg[nlast, 2:(nVar +
                                                 1)]
        eval(parse(text = paste("StockInitState$model",
                                iMod, "<- outreg[, 2:(nVar + 1)]", sep = "")))
      }
    }
    Istep <- 2 * Istep
    if (show == 1) {
      if (nVar <= 3) nsubplot <- sum(okMod)
      else nsubplot <- 2* sum(okMod)
      if (nsubplot > 16) {
        op <- par(mfrow = c(4, 6), pty = "s")
      }
      else {
        if (nsubplot > 0) {
          op <- par(mfrow = c(ceiling(sqrt(nsubplot)),
                              ceiling(sqrt(nsubplot))), pty = "s")
        }
      }
    }
  }
  if (verbose == 1) {
    # Number of unclassified models
    block <- paste("### Number of unclassified models: ", sum(okMod),
                   " / ", length(listMod), sep="")
    cat(block, "\n")
  }
  # return
  list(okMod = okMod, tout = tout, stockoutreg = StockInitState, coeff = kmod,
       models = allKL)
}
