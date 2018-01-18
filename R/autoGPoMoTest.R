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
#' @param IstepMin The minimum number of integration step to start
#' of the analysis (by default \code{IstepMin = 10}).
#' @param IstepMax The maximum number of integration steps for
#' stopping the analysis (by default \code{IstepMax = 10000}).
#' @param tooFarThreshold Divergence threshold, maximum value
#' of the model trajectory compared to the data standard
#' deviation. By default a trjactory is too far if
#' the distance to the center is larger than four times the variance
#' of the input data.
#' @param LimCyclThreshold Threshold used to detect the limit cycle.
#' @param fixedPtThreshold Threshold used to detect fixed points.
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
                           allKL = allKL, IstepMin = 10,
                           IstepMax = 10000, tooFarThreshold = 4, LimCyclThreshold = 0.0,
                           fixedPtThreshold = 1E-8, method = 'rk4')
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
  while (sum(okMod == 1) != 0 & Istep < IstepMax) {
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
    #
    whatModel <- (okMod == 1) * rep(1:length(okMod == 1))
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

        ratio <- (outreg[nlast, 2:(nVar + 1)] - datamoy) / datasd
        iznogood <- ratio[abs(ratio) > tooFarThreshold]
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
                if (abs(ratio[i]) <= tooFarThreshold) {
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
      if (nVar <= 3) nsubplot <- sum(okMod == 1)
      else nsubplot <- 2* sum(okMod == 1)
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
    block <- paste("### Number of unclassified models: ", sum(okMod == 1),
                   " / ", length(listMod), sep="")
    cat(block, "\n")
  }
  # return
  list(okMod = okMod, tout = tout, stockoutreg = StockInitState, coeff = kmod,
       models = allKL)
}
