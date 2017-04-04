#' @title predictab : estimate the models predictability of gPoMo output
#'
#' @description The algorithm aims to estimate automatically the predictability
#' of the models obtained with `gPoMo`.
#'
#' @inheritParams  gloMoId
#' @inheritParams  drvSucc
#' @inheritParams  poLabs
#'
#' @param ogp The output list obtained from gPoMo.
#' @param selecmod a vector of the model selected.
#' @param id the type of model to identify. id = 1 correspond to the unidentified
#' models, that is, potentialy chaotic models).
#' @param fullt Time vector of the data set for which predictability will be tested
#' @param fulldata Data set for which predictability will be tested
#' @param hp Time vector of the horizon of prediction
#' @param Nech Number of simulations
#' @return ErrmodAll A list of matrix $Errmod1, $Errmod2, etc.
#' providing the forecasting error of each model 1, 2, etc.
#' Each column corresponds to one simulation starting from
#' specific initial condition. Each line corresponds to
#' the the horizon of prediction.
#' Vectors corresponding to the initial condition time tE
#' and the horizon of prediction hpE are also provided
#' in $tE and $hpE, respectively.
#'
#'
#' @importFrom grDevices dev.new
#' @importFrom graphics lines par plot
#' @importFrom stats qt
#' 
#' @author Sylvain Mangiarotti
#'
#' @examples
#' # load data
#' data("Ross76")
#' # time vector
#' tin <- Ross76[seq(1, 3000, by = 8), 1]
#' # single time series
#' data <- Ross76[seq(1, 3000, by = 8), 3]
#' # dev.new()
#' # plot(tin, data, xlab = 'time', ylab = 'y(t)')
#'
#' # global modelling
#' # results are put in list outputGPoM
#' outputGPoM <- gPoMo(data[1:300], tin = tin[1:300], dMax = 2, nS=c(3),
#'                     show = 0, method = 'rk4',
#'                     nPmin = 10, nPmax = 12,
#'                     IstepMin = 150, IstepMax = 151)
#' #
#' visuOutGP(outputGPoM)
#'
#' ###########################
#' # and test predictability #
#' ###########################
#' outpred <- predictab(outputGPoM, hp = 15, Nech = 30)
#'
#' # manual visualisation of the outputs (e.g. for model 1):
#' dev.new()
#' image(outpred$tE, outpred$hpE, t(outpred$Errmod1),
#' xlab = 't', ylab = 'hp', main = 'Errmod1')
#'
#' @export
predictab <- function (ogp, fullt = NULL, fulldata = NULL,
                       hp = NULL, Nech = 50, show = 1,
                       selecmod = NULL, id = 1)
{
  if (is.null(selecmod)) {
    # How many models could not be identified as non-chaotic
    nmod <- sum(ogp$okMod == id)
    cat(nmod, 'models identified as id = ', id, "\n")
    #
    # what is the reference number of these models?
    which(ogp$okMod == id)
    seemod <- which(ogp$okMod == id)
  }
  else {
    nmod <- length(selecmod)
    seemod <- selecmod
  }
  #
  nVar <- dim(ogp$models$mToTest1)[2]
  pMax <- dim(ogp$models$mToTest1)[1]
  dMax <- p2dMax(nVar, pMax)
  # maximum integration time steps of each tested models
  modInfo <- matrix(0, ncol = 3, nrow = length(ogp$ok))
  # initiate output
  ErrmodAll <- list()
  IstepMax <- NULL
  Np <- NULL
  for (imod in 1:length(ogp$ok) ) {
    block <- paste("IstepMax <- dim(ogp$stockoutreg$model", imod, ")[1]", sep="")
    eval((parse(text = block)))
    block <- paste("Np <- sum(ogp$models$mToTest", imod, "!=0)", sep="")
    eval((parse(text = block)))
    modInfo[imod, 1] = ogp$ok[imod]
    modInfo[imod, 2] = Np
    if (!is.null(IstepMax)) modInfo[imod, 3] = IstepMax
  }

  if (nmod <= 24) {
    if (show == 1) {
      # plot data and models time series
      dev.new()
      op <- par(mfrow = c(4, 6), pty = "m")
      if (nmod <= 20) op <- par(mfrow = c(4, 5), pty = "m")
      if (nmod <= 18) op <- par(mfrow = c(3, 6), pty = "m")
      if (nmod <= 15) op <- par(mfrow = c(3, 5), pty = "m")
      if (nmod <= 12) op <- par(mfrow = c(3, 4), pty = "m")
      if (nmod <= 9)  op <- par(mfrow = c(3, 3), pty = "m")
      if (nmod <= 4)  op <- par(mfrow = c(2, 2), pty = "m")
      if (nmod <= 1)  op <- par(mfrow = c(1, 1), pty = "m")
    }

    Kmod <- NULL
    for (imod in seemod)  {
      # models are tested one by one
      #
      # prepare the model:
      # equations
      block <- paste("Kmod <- ogp$models$model", imod, sep="")
      eval((parse(text = block)))
      # data
      if (is.null(fulldata) | is.null(fullt)) {
        fullt <- ogp$tfiltdata
        fulldata <- ogp$filtdata
      }
      else {
        if (min(dim(as.matrix(fulldata))) == 1) {
          stop('variable fulldata should be provided with its derivatives')
        }
      }
      # maximum prediction steps
      if (is.null(hp)) hp <- round(dim(fulldata)[1] / (Nech + 1))
      # nombre of sample to consider
      if (Nech > (dim(fulldata)[1] - hp - 1)) {
        warning('Nech and hp values is uncoherent: Nech is modified')
        Nech <- dim(fulldata)[1] - hp - 1
      }
      # prepare matrix output

      # start predictions
      if (show == 1) {
        #dev.new()
        block0 <- paste("Model ", imod, " (Np = ", modInfo[imod, 2], ")", sep="")
        plot(fullt, fulldata[,1], type='l', xlab = 't', ylab = 'f(t)', main = block0)
        lines(fullt, fulldata[,1], type='p', cex = 0.5)
      }
      iinit <- seq(1, length(fullt) - hp, by = floor((length(fullt) - hp - 1) / Nech))
      Errmod = c()
      for (i in 1:Nech) {
        outNumi <- numicano(nVar = nVar, dMax = dMax, KL = Kmod,
                            Istep = hp, v0 = fulldata[iinit[i],1:nVar],
                            method = 'rk4',
                            onestep = ogp$tfiltdata[2] - ogp$tfiltdata[1])
        Err <- outNumi$reconstr[,2] - fulldata[iinit[i]:(iinit[i] + hp - 1),1]
        if (show == 1) {
          lines(fullt[iinit[i]] + outNumi$reconstr[,1], outNumi$reconstr[,2],
                type='l', col='red')
          lines(fullt[iinit[i]] + outNumi$reconstr[,1], outNumi$reconstr[,2],
                type='p', col='red')
          lines(fullt[iinit[i]] + outNumi$reconstr[,1], Err, col= 'green')
        }
        #
        # stockage
        Errmod <- cbind(Errmod, Err)
      }
      block <- paste("ErrmodAll$Errmod", imod," <- Errmod", sep="")
      eval((parse(text = block)))
    }

    if (show == 1) {
      dev.new()
      op <- par(mfrow = c(4, 6), pty = "m")
      if (nmod <= 20) op <- par(mfrow = c(4, 5), pty = "m")
      if (nmod <= 18) op <- par(mfrow = c(3, 6), pty = "m")
      if (nmod <= 15) op <- par(mfrow = c(3, 5), pty = "m")
      if (nmod <= 12) op <- par(mfrow = c(3, 4), pty = "m")
      if (nmod <= 9)  op <- par(mfrow = c(3, 3), pty = "m")
      if (nmod <= 4)  op <- par(mfrow = c(2, 2), pty = "m")
      if (nmod <= 1)  op <- par(mfrow = c(1, 1), pty = "m")
    }
    for (imod in seemod) {
      block0 <- paste("Model ", imod, " (Np = ", modInfo[imod, 2], ")", sep="")
      hpE <- fullt[iinit[1]:(iinit[1] + hp - 1)]
      tE <- fullt[iinit]
      if (show == 1) {
        block <- paste("image(tE, hpE, t(ErrmodAll$Errmod", imod,
                       "), xlab = 't', ylab = 'hp', main = block0)", sep="")
        eval((parse(text = block)))
      }
    }
    ErrmodAll$tE <- tE
    ErrmodAll$hpE <- hpE

    if (show == 1) {
      dev.new()
      op <- par(mfrow = c(4, 6), pty = "m")
      if (nmod <= 20) op <- par(mfrow = c(4, 5), pty = "m")
      if (nmod <= 18) op <- par(mfrow = c(3, 6), pty = "m")
      if (nmod <= 15) op <- par(mfrow = c(3, 5), pty = "m")
      if (nmod <= 12) op <- par(mfrow = c(3, 4), pty = "m")
      if (nmod <= 9)  op <- par(mfrow = c(3, 3), pty = "m")
    }
    for (imod in seemod) {
      block0 <- paste("Model ", imod, " (Np = ", modInfo[imod, 2], ")", sep="")
      hpE <- fullt[iinit[1]:(iinit[1] + hp - 1)]
      if (show == 1) {
        block <- paste("plot(hpE, ErrmodAll$Errmod", imod,
                       "[,i], type = 'l', xlab = 'hp', ylab = 'err',
                     ylim = c(min(ErrmodAll$Errmod", imod,
                       "), max(ErrmodAll$Errmod", imod,
                       ")), col = 'green', main = block0)",
                       sep="")
        eval((parse(text = block)))
        for (i in 1:Nech) {
          block <- paste("lines(hpE, ErrmodAll$Errmod", imod,
                         "[,i], xlab = 'hp', ylab = 'err', col = 'green', main = block0)",
                         sep="")
          eval((parse(text = block)))
        }
      }
      block <- paste("qt <- apply(ErrmodAll$Errmod", imod,
                     ", 1, stats::quantile)",
                     sep="")
      eval((parse(text = block)))
      if (show == 1) {
        lines(hpE, qt[2,], lty = 2)
        lines(hpE, qt[3,], lty = 1)
        lines(hpE, qt[4,], lty = 2)
      }
    }

    if (show == 1) {
      dev.new()
      op <- par(mfrow = c(4, 6), pty = "m")
      if (nmod <= 20) op <- par(mfrow = c(4, 5), pty = "m")
      if (nmod <= 18) op <- par(mfrow = c(3, 6), pty = "m")
      if (nmod <= 15) op <- par(mfrow = c(3, 5), pty = "m")
      if (nmod <= 12) op <- par(mfrow = c(3, 4), pty = "m")
      if (nmod <= 9)  op <- par(mfrow = c(3, 3), pty = "m")
    }
    for (imod in seemod) {
      block0 <- paste("Model ", imod, " (Np = ", modInfo[imod, 2], ")", sep="")
      hpE <- fullt[iinit[1]:(iinit[1] + hp - 1)]
      if (show == 1) {
        block <- paste("plot(hpE, ErrmodAll$Errmod", imod,
                       "[,i]^2, type = 'l', xlab = 'hp', ylab = 'err^2',
                     ylim = c(min(ErrmodAll$Errmod", imod,
                       "^2), max(ErrmodAll$Errmod", imod,
                       "^2)), col = 'green', main = block0)",
                       sep="")
        eval((parse(text = block)))
        for (i in 1:Nech) {
          block <- paste("lines(hpE, ErrmodAll$Errmod", imod,
                         "[,i]^2, xlab = 'hp', ylab = 'err', col = 'green', main = block0)",
                         sep="")
          eval((parse(text = block)))
        }
      }
      block <- paste("qt <- apply(ErrmodAll$Errmod", imod,
                     "^2, 1, stats::quantile)",
                     sep="")
      if (show == 1) {
        eval((parse(text = block)))
        lines(hpE, qt[2,], lty = 2)
        lines(hpE, qt[3,], lty = 1)
        lines(hpE, qt[4,], lty = 2)
      }
    }

  }

  return(ErrmodAll)
}
