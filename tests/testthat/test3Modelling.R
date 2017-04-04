#
# TestModelling
context("TestModelling")
# This file aims to show how to obtain
# models of Ordinary Differential Equations
# of polynomial form from either single
# or multiple time series

#devtools::wd(path='./tests/testthat/')

test_that("First example", {
  # load data
  data("Ross76")
  # time vector
  tin <- Ross76[,1]
  # single time series
  data <- Ross76[,3]
  # generalized Polynomial modelling
  dev.new()
  out1 <- gPoMo(data, tin=tin, dMax = 2, nS=c(3), show = 1,
                nPmin = 9, nPmax = 11)
  expect_equal_to_reference(out1, "./Meta/TestModelling_out1.rds")

})

test_that("Second example", {

  # load data
  data("Ross76")
  # time vector
  tin <- Ross76[,1]
  # multiple (three) time series
  data <- Ross76[,2:4]
  # generalized Polynomial modelling
  dev.new()
  out2 <- gPoMo(data, tin=tin, dMax = 2, nS=c(1,1,1),  show = 1,
                         IstepMin = 10, IstepMax = 3000, nPmin = 7, nPmax = 8)
  # the simplest model able to reproduce the observed dynamics is model #5
  visuEq(3, 2, out2$models$model5, approx = 4) # the original Rossler system is thus retrieved
  # comparison with machine 2 (goeland)
  #expect_equal_to_reference(out2, "./Meta/TestModelling_out2goeland_2mars2017.rds")
  # comparison with reference results 1 (machine 1)
#  expect_equal_to_reference(out2, "./Meta/TestModelling_out2.rds")
})


test_that("Third example", {

  # load data
  data("Ross76")
  # time vector
  tin <- Ross76[,1]
  # multiple (two) time series
  data <- Ross76[,2:3]
  # model template:
  EqS <- matrix(1, ncol = 3, nrow = 10)
  EqS[,1] <- c(0,0,0,1,0,0,0,0,0,0)
  EqS[,2] <- c(1,1,0,1,0,1,1,1,1,1)
  EqS[,3] <- c(0,1,0,0,0,0,1,1,0,0)
  visuEq(3, 2, EqS, substit = c('X','Y','Z'))
  dev.new()
  # generalized Polynomial modelling
  out3 <- gPoMo(data, tin=tin, dMax = 2, nS=c(2,1),
                show = 1, verbose = 1,
                EqS = EqS,
                IstepMin = 10, IstepMax = 2000, nPmin = 9, nPmax = 11)
  expect_equal_to_reference(out3, "./Meta/TestModelling_out3.rds")
})

test_that("Fourth example", {

  # load data
  data(sprottK)
  data(rossler)
  # multiple (six) time series
  data <- cbind(rossler,sprottK)[1:400,]
  dev.new()
  # generalized Polynomial modelling
  out4 <- gPoMo(data, dtFixe = 1/20, dMax = 2, nS=c(1,1,1,1,1,1),
                show = 1, method = 'rk4',
                IstepMin = 10, IstepMax = 11,
                nPmin = 13, nPmax = 13)
  # Only one system can reproduce the phase portraits
  # of the two original systems (Rössler-1976 and Sprott-K):
  visuEq(6, 2, out4$models$model347, approx = 2, substit = 1)
  # to check the robustness of the model, the integration duration
  # should be chosen longer (at least IstepMax = 4000)
  expect_equal_to_reference(out4, "./Meta/TestModelling_out4.rds")
})
