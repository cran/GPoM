#' @title Estimate the models performance obtained with \code{GPoMo}
#' in term of predictability
#'
#' @description The algorithm aims to estimate automatically the forecasting
#' performances of the models obtained with \code{gPoMo}.
#'
#' @inheritParams  gloMoId
#' @inheritParams  drvSucc
#' @inheritParams  poLabs
#'
#' @param ogp The output list obtained from function \code{gPoMo}.
#' @param selecmod A vector of the model selected.
#' @param id The type of model to identify. \code{id = 1} corresponds
#' to unidentified models, that is, potentialy chaotic.
#' @param fullt Time vector of the data set for which predictability
#' will be tested
#' @param fulldata Data set for which predictability will be tested
#' @param hp Time vector of the horizon of prediction
#' @param Nech Number of simulations
#' @param intSimStep Internal number of simulation steps
#' @param selV Selected variable for the analysis
#' @param na.rm Indicates if the \code{NA} should be removed
#' (\code{na.rm = TRUE}) or not (\code{na.rm = FALSE}).
#'
#' @return \code{ErrmodAll} A list of matrix \code{$Predmod1},
#' \code{$Predmod2}, etc. and \code{$Errmod1}, \code{$Errmod2}, etc.
#' providing respectively the forecasting and the forecasting error
#' of models 1, 2, etc.
#' Each column corresponds to one simulation starting from
#' a specific initial condition. Each line corresponds to
#' one horizon of prediction.
#' Vectors corresponding to the initial condition time \code{tE}
#' and the horizon of prediction \code{hpE} are also provided
#' in \code{$tE} and \code{$hpE}, respectively.
#' The percentiles of the distributions of error growth
#' are provided in \code{qt} (0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)
#' and of absolute error growth in \code{qt2} (0.5, 0.75, 0.9, 0.95, 0.98, 0.99).
#'
#'
#' @importFrom grDevices dev.new
#' @importFrom graphics lines par plot
#' @importFrom stats qt
#'
#' @author Sylvain Mangiarotti, Mireille Huc
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
#'
predictab <- function (ogp, fullt = NULL, fulldata = NULL,
                        hp = NULL, Nech = 50, intSimStep = NULL,
                        show = 1, selecmod = NULL, id = 1,
                        selV = 1,na.rm = FALSE){

  if (is.null(selecmod)) {
    # How many models could not be identified as non-chaotic
    nmod <- sum(ogp$okMod == id)
    message(nmod, 'models identified as id = ', id, "\n")
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
  if (is.null(Nech) & is.null(hp) & is.null(intSimStep)) {
    Nech <- 10
    hp <- floor(dim(fulldata)[1] / Nech)
    intSimStep <- floor((dim(fulldata)[1] - hp) / Nech)
  }
  if (is.null(hp) & is.null(intSimStep)) {
    hp <- floor(dim(fulldata)[1] / Nech)
    intSimStep <- floor((dim(fulldata)[1] - hp) / Nech)
  }
  #
  nVar <- dim(ogp$models$mToTest1)[2]
  pMax <- dim(ogp$models$mToTest1)[1]
  dMax <- p2dMax(nVar, pMax)
  # maximum integration time steps of each tested models
  modInfo <- matrix(0, ncol = 3, nrow = length(ogp$okMod))
  # initiate output
  ErrmodAll <- list()
  IstepMax <- NULL
  Np <- NULL
  for (imod in 1:length(ogp$okMod) ) {
    block <- paste("IstepMax <- dim(ogp$stockoutreg$model", imod, ")[1]", sep="")
    eval((parse(text = block)))
    block <- paste("Np <- sum(ogp$models$mToTest", imod, "!=0)", sep="")
    eval((parse(text = block)))
    modInfo[imod, 1] = ogp$okMod[imod]
    modInfo[imod, 2] = Np
    if (!is.null(IstepMax)) modInfo[imod, 3] = IstepMax
  }

  if (nmod <= 24) {
    if (show == 1) {
      # plot data and models time series
      dev.new()
      oldpar <- par(no.readonly = TRUE)    
      on.exit(par(oldpar))            
      op <- par(mfrow = c(4, 6), pty = "m")
      if (nmod <= 20) op <- par(mfrow = c(4, 5), pty = "m")
      if (nmod <= 18) op <- par(mfrow = c(3, 6), pty = "m")
      if (nmod <= 15) op <- par(mfrow = c(3, 5), pty = "m")
      if (nmod <= 12) op <- par(mfrow = c(3, 4), pty = "m")
      if (nmod <= 9)  op <- par(mfrow = c(3, 3), pty = "m")
      if (nmod <= 6)  op <- par(mfrow = c(3, 2), pty = "m")
      if (nmod <= 4)  op <- par(mfrow = c(2, 2), pty = "m")
      if (nmod == 2)  op <- par(mfrow = c(2, 1), pty = "m")
      if (nmod == 1)  op <- par(mfrow = c(1, 1), pty = "m")
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
#      if (is.null(hp)) hp <- round(dim(fulldata)[1] / (Nech + 1))
      if (is.null(hp)) hp <- floor( dim(fulldata)[1] - Nech*intSimStep + 1 )
      if (is.null(intSimStep)) intSimStep <- floor( (dim(fulldata)[1] - hp) / Nech )
      if (is.null(Nech)) Nech <- floor( (dim(fulldata)[1] - hp + 1) / intSimStep )
      # nombre of sample to consider
      if (Nech > (dim(fulldata)[1] - hp - 1)) {
        warning('Nech and hp values is uncoherent: Nech is modified in Nech = ', Nech)
        Nech <- dim(fulldata)[1] - hp - 1
      }
      # prepare matrix output

      # start predictions
      if (show == 1) {
        block0 <- paste("Model ", imod, " (Np = ", modInfo[imod, 2], ")", sep="")
        plot(fullt, fulldata[, selV], type = "l", xlab ="t",
        ylab = "f(t)", main = block0)
        lines(fullt, fulldata[, selV], type = "p", cex = 0.5)
      }
      iinit <- seq(1, length(fullt) - hp,
                   by = intSimStep)[1:Nech] #changed line
      Errmod = c()
      rErrmod = c()
      Predmod = c()
      for (i in 1:Nech) {
        outNumi <- numicano(nVar = nVar, dMax = dMax, KL = Kmod,
                            Istep = hp, v0 = fulldata[iinit[i],1:nVar],
                            method = 'rk4',
                            onestep = ogp$tfiltdata[2] - ogp$tfiltdata[1])
        selV1 <- selV +1
        originalDat <- fulldata[iinit[i]:(iinit[i] + hp - 1), selV]
        Err1 <- outNumi$reconstr[,selV1] - originalDat
#        Err <- outNumi$reconstr[,2] - fulldata[iinit[i]:(iinit[i] + hp - 1),1]
        rErr1 <- (outNumi$reconstr[,selV1] - originalDat) / originalDat * 100
        Pred <- outNumi$reconstr[,2:(nVar+1)]
        if (show == 1) {
          lines(fullt[iinit[i]] + outNumi$reconstr[,1], outNumi$reconstr[, selV1],
                type = "l", col = "red")
          lines(fullt[iinit[i]] + outNumi$reconstr[,1], outNumi$reconstr[, selV1],
                type = "p", col = "red", cex = 0.5)
          lines(fullt[iinit[i]] + outNumi$reconstr[,1], Err1,
                col = "green")
        }
        #
        # stockage
        Errmod <- cbind(Errmod, Err1)
        rErrmod <- cbind(rErrmod, rErr1)
        Predmod <- cbind(Predmod, Pred)
      }
      block <- paste("ErrmodAll$Errmod", imod," <- Errmod", sep="")
      eval((parse(text = block)))
      block <- paste("ErrmodAll$rErrmod", imod," <- rErrmod", sep="")
      eval((parse(text = block)))
      block <- paste("ErrmodAll$Predmod", imod," <- Predmod", sep="")
      eval((parse(text = block)))
    }

    if (show == 1) {
      dev.new()
      oldpar <- par(no.readonly = TRUE)    
      on.exit(par(oldpar))
      op <- par(mfrow = c(4, 6), pty = "m")
      if (nmod <= 20) op <- par(mfrow = c(4, 5), pty = "m")
      if (nmod <= 18) op <- par(mfrow = c(3, 6), pty = "m")
      if (nmod <= 15) op <- par(mfrow = c(3, 5), pty = "m")
      if (nmod <= 12) op <- par(mfrow = c(3, 4), pty = "m")
      if (nmod <= 9)  op <- par(mfrow = c(3, 3), pty = "m")
      if (nmod <= 6)  op <- par(mfrow = c(3, 2), pty = "m")
      if (nmod <= 4)  op <- par(mfrow = c(2, 2), pty = "m")
      if (nmod == 2)  op <- par(mfrow = c(2, 1), pty = "m")
      if (nmod == 1)  op <- par(mfrow = c(1, 1), pty = "m")
    }
    for (imod in seemod) {
      onestep <- ogp$tfiltdata[2] - ogp$tfiltdata[1]
      block0 <- paste("Model ", imod, " (Np = ", modInfo[imod, 2], ")", sep="")
      hpE <- fullt[iinit[1]:(iinit[1] + hp - 1)] - fullt[1]
      tE <- fullt[iinit] - fullt[1]
      ntE <- NULL
      block <- paste("ntE <- dim(ErrmodAll$Errmod", imod,")[2]", sep="")
      eval((parse(text = block)))
      tE <- tE[1:ntE]
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
      oldpar <- par(no.readonly = TRUE)    
      on.exit(par(oldpar))
      op <- par(mfrow = c(4, 6), pty = "m")
      if (nmod <= 20) op <- par(mfrow = c(4, 5), pty = "m")
      if (nmod <= 18) op <- par(mfrow = c(3, 6), pty = "m")
      if (nmod <= 15) op <- par(mfrow = c(3, 5), pty = "m")
      if (nmod <= 12) op <- par(mfrow = c(3, 4), pty = "m")
      if (nmod <= 9)  op <- par(mfrow = c(3, 3), pty = "m")
      if (nmod <= 6)  op <- par(mfrow = c(3, 2), pty = "m")
      if (nmod <= 4)  op <- par(mfrow = c(2, 2), pty = "m")
      if (nmod == 2)  op <- par(mfrow = c(2, 1), pty = "m")
      if (nmod == 1)  op <- par(mfrow = c(1, 1), pty = "m")
    }
    for (imod in seemod) {
      block0 <- paste("Model ", imod, " (Np = ", modInfo[imod, 2], ")", sep="")
      hpE <- fullt[iinit[1]:(iinit[1] + hp - 1)] - fullt[1]
      if (show == 1) {
        block <- paste("plot(hpE, ErrmodAll$Errmod", imod,
                       "[,i], type = 'l', xlab = 'hp', ylab = expression(err(t)),
                     ylim = c(min(ErrmodAll$Errmod", imod,", na.rm = na.rm),
                              max(ErrmodAll$Errmod", imod,", na.rm = na.rm)),
                       col = 'green', main = block0)",
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
                     ", 1, stats::quantile, probs = c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95), na.rm=na.rm)",
                     sep="")
      eval((parse(text = block)))
      block <- paste("ErrmodAll$QTmod", imod,"$qt <- qt", sep="")
      eval((parse(text = block)))
      #
      block <- paste("rqt <- apply(ErrmodAll$rErrmod", imod,
                     ", 1, stats::quantile, probs = c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95), na.rm=na.rm)",
                     sep="")
      eval((parse(text = block)))
      block <- paste("ErrmodAll$QTmod", imod,"$rqt <- rqt", sep="")
      eval((parse(text = block)))
      #
      if (show == 1) {
        lines(hpE, qt[1,], lty = 2, lwd = 2)
        lines(hpE, qt[2,], lty = 4, lwd = 2)
        lines(hpE, qt[3,], lty = 3, lwd = 2)
        lines(hpE, qt[4,], lty = 1, lwd = 3)
        lines(hpE, qt[5,], lty = 3, lwd = 2)
        lines(hpE, qt[6,], lty = 4, lwd = 2)
        lines(hpE, qt[7,], lty = 2, lwd = 2)
      }
    }

    if (show == 1) {
      dev.new()
      oldpar <- par(no.readonly = TRUE)    
      on.exit(par(oldpar))
      op <- par(mfrow = c(4, 6), pty = "m")
      if (nmod <= 20) op <- par(mfrow = c(4, 5), pty = "m")
      if (nmod <= 18) op <- par(mfrow = c(3, 6), pty = "m")
      if (nmod <= 15) op <- par(mfrow = c(3, 5), pty = "m")
      if (nmod <= 12) op <- par(mfrow = c(3, 4), pty = "m")
      if (nmod <= 9)  op <- par(mfrow = c(3, 3), pty = "m")
      if (nmod <= 6)  op <- par(mfrow = c(3, 2), pty = "m")
      if (nmod <= 4)  op <- par(mfrow = c(2, 2), pty = "m")
      if (nmod == 2)  op <- par(mfrow = c(2, 1), pty = "m")
      if (nmod == 1)  op <- par(mfrow = c(1, 1), pty = "m")
    }
    for (imod in seemod) {
      block0 <- paste("Model ", imod, " (Np = ", modInfo[imod, 2], ")", sep="")
      hpE <- fullt[iinit[1]:(iinit[1] + hp - 1)] - fullt[1]
      if (show == 1) {

        block <- paste("plot(hpE, abs(ErrmodAll$Errmod", imod,
                       "[,i]), type = 'l', xlab = 'hp', ylab = expression(italic(abs(err))),
                        ylim = c(min(abs(ErrmodAll$Errmod", imod, "), na.rm = na.rm),
                              max(abs(ErrmodAll$Errmod", imod, "), na.rm = na.rm)),
                       col = 'green', main = block0)",
                       sep="")
        eval((parse(text = block)))
        for (i in 1:Nech) {
          block <- paste("lines(hpE, abs(ErrmodAll$Errmod", imod,
                         "[,i]), col = 'green', main = block0)",
                         sep="")
          eval((parse(text = block)))
        }
      }
      qt2 <- NULL
      block <- paste("qt2 <- apply(abs(ErrmodAll$Errmod", imod,
                     "), 1, stats::quantile, probs = c(0.5, 0.75, 0.9, 0.95, 0.98, 0.99), na.rm = na.rm)",
                     sep="")
      eval((parse(text = block)))
      block <- paste("ErrmodAll$QTmod", imod,"$qt2 <- qt2", sep="")
      eval((parse(text = block)))
      #
      block <- paste("rqt2 <- apply(abs(ErrmodAll$rErrmod", imod,
                     "), 1, stats::quantile, probs = c(0.5, 0.75, 0.9, 0.95, 0.98, 0.99), na.rm = na.rm)",
                     sep="")
      eval((parse(text = block)))
      block <- paste("ErrmodAll$QTmod", imod,"$rqt2 <- rqt2", sep="")
      eval((parse(text = block)))
      #
      if (show == 1) {
        eval((parse(text = block)))
        lines(hpE, qt2[2,], lty = 2, lwd = 2)
        lines(hpE, qt2[3,], lty = 1, lwd = 2)
        lines(hpE, qt2[4,], lty = 2, lwd = 2)
      }
    }

  }
  # input data is also kept for memory
  ErrmodAll$inputtime <- fullt
  ErrmodAll$inputdata <- fulldata
  ErrmodAll$hp <- hp
  ErrmodAll$Nech <- Nech
  return(ErrmodAll)
}
