#
# TestVisuOutput
context("TestVisuOutput")
# This file aims to show obtain a
# quick view of the results obtained
# with gPoMo

#devtools::wd(path='./tests/testthat/')

test_that("First example", {  # load data
  data("Ross76")
  # time vector
  tin <- Ross76[seq(1, 3000, by = 8), 1]
  # single time series
  data <- Ross76[seq(1, 3000, by = 8), 3]
  dev.new()
  plot(tin, data)
  # global modelling
  # results are put in list outputGPoM
  outputGPoM <- gPoMo(data, tin=tin, dMax = 2, nS=c(3), show = 0,
                nPmax = 12, IstepMin = 400, IstepMax = 401)
  #
  ######################
  # automatic ##########
  ######################
  visuOutGP(outputGPoM)

  ######################
  # manual    ##########
  ######################
  # How many models could not be identified as non-chaotic
  sum(outputGPoM$okMod)
  #
  # what is the reference number of these models?
  which(outputGPoM$okMod == 1)
  #
  # maximum integration time steps check
  dim(outputGPoM$stockoutreg$model7)[1]
  #
  # plot data and models time series
  dev.new()
  op <- par(mfrow = c(2, 2), pty = "s")
  plot(outputGPoM$tin, outputGPoM$inputdata,
       xlab='t', ylab='y(t)', cex = 0.7, main = 'Model 7')
  lines(outputGPoM$tfiltdata, outputGPoM$filtdata[,1], type='l', col='green')
  lines(outputGPoM$tout, outputGPoM$stockoutreg$model7[,1], type='l', col='red')
  #
  plot(outputGPoM$tin, outputGPoM$inputdata,
       xlab='t', ylab='y(t)', cex = 0.7, main = 'Model 9')
  lines(outputGPoM$tfiltdata, outputGPoM$filtdata[,1], type='l', col='green')
  lines(outputGPoM$tout, outputGPoM$stockoutreg$model9[,1], type='l', col='red')
  #
  plot(outputGPoM$tin, outputGPoM$inputdata,
       xlab='t', ylab='y(t)', cex = 0.7, main = 'Model 10')
  lines(outputGPoM$tfiltdata, outputGPoM$filtdata[,1], type='l', col='green')
  lines(outputGPoM$tout, outputGPoM$stockoutreg$model10[,1], type='l', col='red')
  #
  # plot data and models phase portraits
  dev.new()
  op <- par(mfrow = c(2, 2), pty = "s")
  plot(outputGPoM$filtdata[,1], outputGPoM$filtdata[,2], type='l', col='green',
       xlab='y(t)', ylab='dy(t)/dt', main = 'Model 7')
  lines(outputGPoM$stockoutreg$model7[,1], outputGPoM$stockoutreg$model7[,2],
        type='l', col='red')
  #
  plot(outputGPoM$filtdata[,1], outputGPoM$filtdata[,2], type='l', col='green',
       xlab='y(t)', ylab='dy(t)/dt', main = 'Model 9')
  lines(outputGPoM$stockoutreg$model9[,1], outputGPoM$stockoutreg$model9[,2],
        type='l', col='red')
  #
  plot(outputGPoM$filtdata[,1], outputGPoM$filtdata[,2], type='l', col='green',
       xlab='y(t)', ylab='dy(t)/dt', main = 'Model 10')
  lines(outputGPoM$stockoutreg$model10[,1], outputGPoM$stockoutreg$model10[,2],
        type='l', col='red')
  #
  # See equations
  visuEq(3,2, outputGPoM$models$model7, approx = 3)
  #
  visuEq(3,2, outputGPoM$models$model9, approx = 3)
  #
  visuEq(3,2, outputGPoM$models$model10, approx = 3)
})
