context("visuEq")

test_that("First example", {
  nVar <- 3
  dMax <- 2
  KL = matrix(0, ncol = 3, nrow = 10)
  KL[1,1] <- KL[2,2] <- 1
  KL[4,1] <- -1
  KL[5,3] <- -0.123456789
  # by default, variables names X1, X2, X3 are used
  M <- visuEq(nVar, dMax, KL, approx = 2)
  expect_equal(M[[1]], "dX1/dt = 1  -X2 ")
  expect_equal(M[[2]], "dX2/dt = X3 ")
  expect_equal(M[[3]], "dX3/dt = -0.123 X2 X3 ")
  })

test_that("Second example", {
  nVar <- 3
  dMax <- 2
  KL = matrix(0, ncol = 3, nrow = 10)
  KL[1,1] <- KL[2,2] <- 1
  KL[4,1] <- -1
  KL[5,3] <- -0.123456789
  # by default, variables names X1, X2, X3 are used
  M <- visuEq(nVar, dMax, KL, substit = 1, approx = 2)
  expect_equal(M[[1]], "dx/dt = 1  -y ")
  expect_equal(M[[2]], "dy/dt = z ")
  expect_equal(M[[3]], "dz/dt = -0.123 y z ")
})

test_that("Third example", {
  nVar <- 3
  dMax <- 2
  KL = matrix(0, ncol = 3, nrow = 10)
  KL[1,1] <- KL[2,2] <- 1
  KL[4,1] <- -1
  KL[5,3] <- -0.123456789
  # by default, variables names X1, X2, X3 are used
  M <- visuEq(nVar, dMax, KL, substit = 2, approx = 2)
  expect_equal(M[[1]], "dX_1/dt = 1 -X_2 ")
  expect_equal(M[[2]], "dX_2/dt = X_3 ")
  expect_equal(M[[3]], "dX_3/dt =-0.123 X_2 X_3 ")
})
