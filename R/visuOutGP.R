#' @title visuOutGP : get a quick information of gPoMo output
#'
#' @description The algorithm aims to get a quick information about
#' the outputs obtained with gPoMo.
#'
#' @inheritParams  gloMoId
#' @inheritParams  drvSucc
#' @inheritParams  poLabs
#'
#' @param ogp The output list obtained from gPoMo.
#' @param selecmod a vector of the model selected.
#' @param id the type of model to identify. id = 1 correspond to the unidentified
#' models, that is, potentialy chaotic models).
#' @param prioMinMax gives the priority for the plots among: "data", "model",
#' "dataonly" and "modelonly".
#' @param opt3D provides a 3D plot (x,y,z) when 'TRUE' (rgl library required).
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
#' plot(tin, data, type = 'l')
#' # global modelling
#' # results are put in list outputGPoM
#' outputGPoM <- gPoMo(data, tin=tin, dMax = 2, nS=c(3), show = 0,
#'                     nPmin = 9, nPmax = 12, method = 'rk4',
#'                     IstepMin = 200, IstepMax = 201)
#' visuOutGP(outputGPoM)
#'
visuOutGP <- function (ogp, selecmod = NULL, id = 1,
                       prioMinMax = 'data', opt3D = 'TRUE')
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
    for (imod in seemod)  {
      block0 <- paste("Model ", imod, " (Np = ", modInfo[imod, 2], ")", sep="")
      #
      if (is.matrix(ogp$inputdata)) inputseries <- ogp$inputdata[,1]
      if (prioMinMax == 'data') {
        plot(ogp$tin, inputseries, cex = 0.7, col='black',
             main = block0, xlab='t', ylab='y(t)')
        lines(ogp$tfiltdata, ogp$filtdata[,1], type='l', col='green')
        #
        block <- paste(
          "lines(ogp$tout[1:modInfo[", imod,
          ",3]], ogp$stockoutreg$model", imod,
          "[,1], type='l', col='red')", sep="")
        eval((parse(text = block)))
      }
      else if (prioMinMax == 'dataonly') {
        plot(ogp$tin, inputseries, cex = 0.7, col='black',
             main = block0, xlab='t', ylab='y(t)')
        lines(ogp$tfiltdata, ogp$filtdata[,1], type='l', col='green')
      }
      else if (prioMinMax == 'model') {
        #
        block <- paste(
          "plot(ogp$tout[1:modInfo[", imod,
          ",3]], ogp$stockoutreg$model", imod,
          "[,1], type='l', col='red', main = block0, xlab='t', ylab='y(t)')", sep="")
        eval((parse(text = block)))
        lines(ogp$tin, inputseries, cex = 0.7, col='black', type='o')
        lines(ogp$tfiltdata, ogp$filtdata[,1], type='l', col='green')
      }
      else if (prioMinMax == 'modelonly') {
        #
        block <- paste(
          "plot(ogp$tout[1:modInfo[", imod,
          ",3]], ogp$stockoutreg$model", imod,
          "[,1], type='l', col='red', main = block0, xlab='t', ylab='y(t)')", sep="")
        eval((parse(text = block)))
      }
    }
  }

  if (nmod <= 24) {
    # plot data and models phase portraits
    dev.new()
    op <- par(mfrow = c(4, 6), pty = "s")
    if (nmod <= 20) op <- par(mfrow = c(4, 5), pty = "s")
    if (nmod <= 18) op <- par(mfrow = c(3, 6), pty = "s")
    if (nmod <= 15) op <- par(mfrow = c(3, 5), pty = "s")
    if (nmod <= 12) op <- par(mfrow = c(3, 4), pty = "s")
    if (nmod <= 9)  op <- par(mfrow = c(3, 3), pty = "s")
    if (nmod <= 4)  op <- par(mfrow = c(2, 2), pty = "s")
    if (nmod <= 1)  op <- par(mfrow = c(1, 1), pty = "s")
    for (imod in seemod)  {
      block0 <- paste("Model ", imod, " (Np = ", modInfo[imod, 2], ")", sep="")
      if (prioMinMax == 'data') {
        plot(ogp$filtdata[,1], ogp$filtdata[,2], type='l', col='green',
           xlab='y(t)', ylab='dy(t)/dt', main = block0)
        block <- paste(
            "lines(ogp$stockoutreg$model", imod,
            "[,1], ogp$stockoutreg$model", imod,
            "[,2], type='l', col='red')", sep="")
        eval((parse(text = block)))
        if (opt3D) {
          plot3d(ogp$filtdata[,1], ogp$filtdata[,2], ogp$filtdata[,3], type='l', col='green',
               xlab='y(t)', ylab='dy(t)/dt', zlab='d2y(t)/dt2', main = block0)
          block <- paste(
            "plot3d(ogp$stockoutreg$model", imod,
            "[,1], ogp$stockoutreg$model", imod,
            "[,2], ogp$stockoutreg$model", imod,
            "[,3], type='l', col='red', add = TRUE)", sep="")
          eval((parse(text = block)))
        }
      }
      else if (prioMinMax == 'dataonly') {
        plot(ogp$filtdata[,1], ogp$filtdata[,2], type='l', col='green',
             xlab='y(t)', ylab='dy(t)/dt', main = block0)
        if (opt3D) {
            plot3d(ogp$filtdata[,1], ogp$filtdata[,2], ogp$filtdata[,3], type='l', col='green',
                   xlab='y(t)', ylab='dy(t)/dt', zlab='d2y(t)/dt2', main = block0)
        }
      }
      else if (prioMinMax == 'model') {
        block <- paste(
          "plot(ogp$stockoutreg$model", imod,
          "[,1], ogp$stockoutreg$model", imod,
          "[,2], type='l', col='red', xlab='y(t)', ylab='dy(t)/dt', main = block0)", sep="")
        eval((parse(text = block)))
        lines(ogp$filtdata[,1], ogp$filtdata[,2], type='l', col='green')
        if (opt3D) {
          block <- paste(
            "plot3d(ogp$stockoutreg$model", imod,
            "[,1], ogp$stockoutreg$model", imod,
            "[,2], ogp$stockoutreg$model", imod,
            "[,3], type='l', col='red', xlab='y(t)', ylab='dy(t)/dt', zlab='d2y(t)/dt2', main = block0)", sep="")
          eval((parse(text = block)))
          plot3d(ogp$filtdata[,1], ogp$filtdata[,2], ogp$filtdata[,3],
                 type='l', col='green', add = TRUE)
        }
      }
      else if (prioMinMax == 'modelonly') {
        block <- paste(
          "plot(ogp$stockoutreg$model", imod,
          "[,1], ogp$stockoutreg$model", imod,
          "[,2], type='l', col='red', xlab='y(t)', ylab='dy(t)/dt', main = block0)", sep="")
        eval((parse(text = block)))
        if (opt3D) {
            block <- paste(
              "plot3d(ogp$stockoutreg$model", imod,
              "[,1], ogp$stockoutreg$model", imod,
              "[,2], ogp$stockoutreg$model", imod,
              "[,3], type='l', col='red', xlab='y(t)', ylab='dy(t)/dt', zlab='d2y(t)/dt2', main = block0)", sep="")
            eval((parse(text = block)))
        }
      }
    }
  }

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
      "visuEq(", nVar,
      ",",dMax,
      ", ogp$models$model", imod,
      ")", sep="")
    eval((parse(text = block)))
  }

  cat('Models classification, size and number of integration steps:', "\n")
  invisible(t(modInfo))
}
