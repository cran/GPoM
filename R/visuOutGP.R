#' @title visuOutGP : get a quick information of gPoMo output
#'
#' @description The algorithm aims to get a quick information about
#' the outputs obtained with gPoMo.
#'
#' @param ogp The output list obtained from gPoMo.
#' @param selecmod A vector of the selected model. Maximum 24 models can be
#' presented at the same time.
#' @param id The type of model to identify. \code{id = 1} corresponds to
#' the unidentified models, that is, potentialy chaotic models).
#' @param prioMinMax Gives the priority for the plots among: \code{"data"},
#' \code{"model"}, \code{"dataonly"} and \code{"modelonly"}.
#' @param opt3D Provides a 3D plot (x,y,z) when \code{opt = 'TRUE'}
#' (the \code{rgl} library is required).
#' @param maxPages The maximum of pages to be displayed (4 by default,
#' but this may be insufficient when too many models remain)
#' @param seeEq Indicates if equations should be displayed (seeEq = 1,
#' by default) or not (seeEq = 0).
#'
#' @return A Matrix describing the terms composing each model by row. The first
#' row corresponds to the model detection (1 unclarified, 2 diverging, 0 is fixed
#' point, -n with n an integer, is period-n cycle' )
#'
#' @author Sylvain Mangiarotti
#'
#' @importFrom grDevices dev.new
#' @importFrom graphics lines par plot
#' @import rgl
#' @export
#' @examples
#' # load data
#' data("Ross76")
#' # # time vector
#' tin <- Ross76[seq(1, 3000, by = 8), 1]
#' # single time series
#' data <- Ross76[seq(1, 3000, by = 8), 3]
#' dev.new()
#' plot(tin, data, type = 'l', main = 'Observed time series')
#' # global modelling
#' # results are put in list outputGPoM
#' outputGPoM <- gPoMo(data, tin=tin, dMax = 2, nS=c(3), show = 0,
#'                     nPmin = 9, nPmax = 12, method = 'rk4',
#'                     IstepMin = 200, IstepMax = 201)
#' visuOutGP(outputGPoM)
#'
visuOutGP <- function (ogp, selecmod = NULL, id = 1,
                       prioMinMax = 'data', opt3D = 'TRUE',
                       maxPages = NULL, seeEq = 1)
{
  if (is.null(selecmod)) {
    # How many models could not be identified as non-chaotic
    nmod <- sum(ogp$okMod == id)
    cat(nmod, 'models identified as id = ', id, ":", "\n")
    cat(which(ogp$okMod == 1), "\n")
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
  if (is.null(selecmod)) {
    maxPages = 4
  }
  #
  nVar <- dim(ogp$models[[1]])[2]
  pMax <- dim(ogp$models[[1]])[1]
  dMax <- p2dMax(nVar, pMax)
  # maximum integration time steps of each tested models
  modInfo <- matrix(0, ncol = 3, nrow = 0)
  IstepMax <- NULL
  Np <- NULL
  for (imod in 1:length(ogp$okMod) ) {
    block <- paste("IstepMax <- dim(ogp$stockoutreg$model", imod, ")[1]", sep="")
    eval((parse(text = block)))
    block <- paste("Np <- sum(ogp$models$mToTest", imod, "!=0)", sep="")
    eval((parse(text = block)))
#    modInfo[imod, 1] = ogp$okMod[imod]
#    modInfo[imod, 2] = Np
    if (!is.null(IstepMax)) {
      modInfo = rbind(modInfo, c(ogp$okMod[imod], Np, IstepMax))
    }
    else {
      modInfo = rbind(modInfo, c(ogp$okMod[imod], Np, 0))
    }
  }
  
  # rewrite the filtered data with NaN where non valid values
  fdat <- ogp$filtdata
  for (i in 1:dim(ogp$filtdata)[2]) {
    fdat[ogp$Wfiltdata==0,i] <- NaN
  }
  # find first valid value
  firstOk <- which(is.nan(ogp$Wfiltdata)==0)[1]

  


  if (nmod <= 24) {
    # plot data and models phase portraits
    dev.new()
    oldpar <- par(no.readonly = TRUE)    
    on.exit(par(oldpar))  
    op <- par(mfrow = c(4, 6), pty = "s")
    if (nmod <= 20) op <- par(mfrow = c(4, 5), pty = "s")
    if (nmod <= 18) op <- par(mfrow = c(3, 6), pty = "s")
    if (nmod <= 15) op <- par(mfrow = c(3, 5), pty = "s")
    if (nmod <= 12) op <- par(mfrow = c(3, 4), pty = "s")
    if (nmod <= 9)  op <- par(mfrow = c(3, 3), pty = "s")
    if (nmod <= 6)  op <- par(mfrow = c(3, 2), pty = "m")
    if (nmod <= 4)  op <- par(mfrow = c(2, 2), pty = "m")
    if (nmod == 2)  op <- par(mfrow = c(2, 1), pty = "m")
    if (nmod == 1)  op <- par(mfrow = c(1, 1), pty = "m")
    for (imod in seemod)  {
      block0 <- paste("Model ", imod, " (Np = ", modInfo[imod, 2], ")", sep="")
      if (prioMinMax == 'data') {
        plot(fdat[,1], fdat[,2], type='l', col='green',
           xlab='X1', ylab='X2', main = block0)
        lines(fdat[firstOk,1], fdat[firstOk,2],
              type='p', col='green',
              xlab='X1', ylab='X2', main = block0, cex = 1.2)
        block <- paste(
            "lines(ogp$stockoutreg$model", imod,
            "[,2], ogp$stockoutreg$model", imod,
            "[,3], type='l', col='red')", sep="")
        eval((parse(text = block)))
        block <- paste(
          "lines(ogp$stockoutreg$model", imod,
          "[1,2], ogp$stockoutreg$model", imod,
          "[1,3], type='p', cex = 1.2, col='red')", sep="")
        eval((parse(text = block)))
        if (opt3D) {
          plot3d(fdat[,1], fdat[,2], fdat[,3], type='l', col='green',
               xlab='X1', ylab='X2', zlab='X3', main = block0)
          plot3d(fdat[firstOk,1], fdat[firstOk,2], fdat[firstOk,3],
                 type='p', size = 5, col='green',
                 xlab='X1', ylab='X2', zlab='X3', main = block0, add = TRUE)
          block <- paste(
            "plot3d(ogp$stockoutreg$model", imod,
            "[,2], ogp$stockoutreg$model", imod,
            "[,3], ogp$stockoutreg$model", imod,
            "[,4], type='l', col='red', add = TRUE)", sep="")
          eval((parse(text = block)))
          block <- paste(
            "plot3d(ogp$stockoutreg$model", imod,
            "[1,2], ogp$stockoutreg$model", imod,
            "[1,3], ogp$stockoutreg$model", imod,
            "[1,4], type='p', size = 5, col='red', add = TRUE)", sep="")
          eval((parse(text = block)))
        }
      }
      else if (prioMinMax == 'dataonly') {
        plot(fdat[,1], fdat[,2], type='l', col='green',
             xlab='X1', ylab='X2', main = block0)
        lines(fdat[firstOk,1], fdat[firstOk,2],
              type='p', col='green',
              xlab='X1', ylab='X2', main = block0, cex = 1.2)
        if (opt3D) {
            plot3d(fdat[,1], fdat[,2], fdat[,3],
                   type='l', col='green',
                   xlab='X1', ylab='X2', zlab='X3', main = block0)
            plot3d(fdat[firstOk,1], fdat[firstOk,2], fdat[firstOk,3],
                   type='p', size = 5, col='green',
                   xlab='X1', ylab='X2', zlab='X3', main = block0)
        }
      }
      else if (prioMinMax == 'model') {
        block <- paste(
          "plot(ogp$stockoutreg$model", imod,
          "[,2], ogp$stockoutreg$model", imod,
          "[,3], type='l', col='red', xlab='X1', ylab='X2', main = block0)", sep="")
        eval((parse(text = block)))
        block <- paste(
          "lines(ogp$stockoutreg$model", imod,
          "[1,2], ogp$stockoutreg$model", imod,
          "[1,3], type='p', cex = 1.2, col='red', xlab='X1', ylab='X2', main = block0)", sep="")
        eval((parse(text = block)))
        lines(fdat[,1], fdat[,2], type='l', col='green')
        lines(fdat[firstOk,1], fdat[firstOk,2],
              type='p', col='green', cex = 1.2)
        if (opt3D) {
          block <- paste(
            "plot3d(ogp$stockoutreg$model", imod,
            "[,2], ogp$stockoutreg$model", imod,
            "[,3], ogp$stockoutreg$model", imod,
            "[,4], type='l', col='red', xlab='X1', ylab='X2', zlab='X3', main = block0)", sep="")
          eval((parse(text = block)))
          block <- paste(
            "plot3d(ogp$stockoutreg$model", imod,
            "[1,2], ogp$stockoutreg$model", imod,
            "[1,3], ogp$stockoutreg$model", imod,
            "[1,4], type='p', size = 5, col='red', xlab='X1', ylab='X2', zlab='X3', main = block0, add = TRUE)", sep="")
          eval((parse(text = block)))
          plot3d(fdat[,1], fdat[,2], fdat[,3],
                 type='l', col='green', add = TRUE)
          plot3d(fdat[firstOk,1], fdat[firstOk,2], fdat[firstOk,3],
                 type='p', size = 5, col='green', add = TRUE)
        }
      }
      else if (prioMinMax == 'modelonly') {
        block <- paste(
          "plot(ogp$stockoutreg$model", imod,
          "[,2], ogp$stockoutreg$model", imod,
          "[,3], type='l', col='red', xlab='X1', ylab='X2', main = block0)", sep="")
        eval((parse(text = block)))
        block <- paste(
          "plot(ogp$stockoutreg$model", imod,
          "[1,2], ogp$stockoutreg$model", imod,
          "[1,3], type='p', cex = 1.2, col='red', xlab='X1', ylab='X2', main = block0)", sep="")
        eval((parse(text = block)))
        if (opt3D) {
            block <- paste(
              "plot3d(ogp$stockoutreg$model", imod,
              "[,2], ogp$stockoutreg$model", imod,
              "[,3], ogp$stockoutreg$model", imod,
              "[,4], type='l', col='red', xlab='X1', ylab='X2', zlab='X3', main = block0)", sep="")
            eval((parse(text = block)))
            block <- paste(
              "plot3d(ogp$stockoutreg$model", imod,
              "[1,2], ogp$stockoutreg$model", imod,
              "[1,3], ogp$stockoutreg$model", imod,
              "[1,4], type='p', size = 5, col='red', xlab='X1', ylab='X2', zlab='X3', main = block0, add = TRUE)", sep="")
            eval((parse(text = block)))
        }
      }
    }
  }
  else {
    if (ceiling(nmod / 24) > maxPages) {
      # All the pages are can be displayed
      for (iPage in 1:maxPages) {
          visuOutGP(ogp,
                    selecmod = which(ogp$okMod == 1)[((iPage-1)*24+1):((iPage)*24)],
                    seeEq = seeEq)
        }
      # TOO MANY MODELS:
      # Maximum 'maxPages' pages are displayed by default
      warning('Too many models nmod = ', nmod,'. By default the fonction ',
           'can display four pages 24 models each (that is nmod <= 96).',
           'Too have more page, please use the option maxPages = ',
           ceiling(nmod / 24))
    }
    else {
      maxPages <- ceiling(nmod / 24)
      # All the pages are can be displayed
      for (iPage in 1:maxPages) {
        if (iPage < maxPages) {
          visuOutGP(ogp,
                    selecmod = which(ogp$okMod == 1)[((iPage-1)*24+1):((iPage)*24)],
                    seeEq = seeEq)
        }
        else if (iPage == maxPages) {
          visuOutGP(ogp,
                    selecmod = which(ogp$okMod == 1)[((iPage-1)*24+1):nmod],
                    seeEq = seeEq)
        }
      }
    }
  }

  if (seeEq == 1) {
    # See equations
    cat('Equations of', "\n")
    #
    ic <- NULL
    for (imod in seemod)  {
      block <- paste("### Model ", imod,
                     " (Np = ", modInfo[imod,2],
                     "):", sep="")
      cat(block, "\n")
      block <- paste("ic <- ogp$stockoutreg$model", imod,
                     "[1,]", sep="")
      eval((parse(text = block)))
      cat("### Initial conditions: ", "\n")
      block <- paste(ic, sep="")
      cat(block, "\n")
      cat("### Equations: ", "\n")
      block <- paste(
        "visuEq(ogp$models$model", imod,
        ")", sep="")
      eval((parse(text = block)))
    }
  }

#  cat('Models classification, size and number of integration steps:', "\n")
  invisible(t(modInfo))
}
