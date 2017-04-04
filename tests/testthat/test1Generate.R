#
# TestGenerate
context("TestGenerate")
# This file aims to generate time series
# from Ordinary Differential Equations
# of polynomial form

# The convention used to formulate the model
# numerically is given by function poLabs(nVar, dMax)
# with nVar the model dimension and dMax the maximum
# polynomial value

test_that("Generate example", {
# The model dimension
nVar = 3
# maximal polynomial degree
dMax = 2

# polynomial structuration given by
poLabs(nVar, dMax)

# Number of parameter number pMax
length(poLabs(nVar, dMax))
# is also given by
pMax <- d2pMax(nVar, dMax)

# Based on this convention, one single ordinary
# differential equation (ODE) of polynomial form
# with N < or = nVar variables and of degree q < or = dMax
# can be formulated using one single vector
# For example, equation
# dX1/dt = 1 + 2 X1 + 3 X1X3 + 4 X2^2
# will be formulated as
# c(1, 0, 0, 0, 0, 4, 2, 3, 0, 0)
cbind(c(1, 0, 0, 0, 0, 4, 2, 3, 0, 0), poLabs(nVar, dMax))

######################
# Example 1 ##########
######################

# a set of N equations can thus be represented as
# a matrix of pMax lines by nVar columns
#
# For example the Rossler system will be defined as follows:
# parameters
a = 0.520
b = 2
c = 4
# equations
Eq1 <- c(0,-1, 0,-1, 0, 0, 0, 0, 0, 0)
Eq2 <- c(0, 0, 0, a, 0, 0, 1, 0, 0, 0)
Eq3 <- c(b,-c, 0, 0, 0, 0, 0, 1, 0, 0)
# model formulation
K <- cbind(Eq1, Eq2, Eq3)
# Equations can be edited with function visuEq()
# substit = 1 is used to have 'x', 'y', 'z' instead of 'X_1', 'X_2', 'X_3'
visuEq(nVar, dMax, K, substit = 1)

# initial conditions
v0 <- c(-0.6, 0.6, 0.4)
# model integration
reconstr <- numicano(nVar, dMax, Istep=5000, onestep=1/50, KL=K,
                     v0=v0, method="ode45")
# Plot of the simulated time series obtained
dev.new()
plot(reconstr$reconstr[,2], reconstr$reconstr[,3], type='l',
     main='phase portrait', xlab='x(t)', ylab = 'y(t)')

######################
# Example 2 ##########
######################

# For a model of canonical form, model can be defined
# by one single polynomial
#
# model dimension:
nVar = 4
# maximal polynomial degree
dMax = 3
# Number of parameter number (by default)
pMax <- d2pMax(nVar, dMax)
# Definition of the Model Function
PolyTerms <- c(281000, 0, 0, 0, -2275, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               861, 0, 0, 0, -878300, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
# terms used in the model
poLabs(nVar, dMax, PolyTerms!=0)
# initial conditions
v0 <- c(0.54, 3.76, -90, -5200)
# model integration
reconstr <- numicano(nVar, dMax, Istep=500, onestep=1/250, PolyTerms=PolyTerms,
                     v0=v0, method="ode45")
# Plot of the simulated time series obtained
plot(reconstr$reconstr[,2], reconstr$reconstr[,3], type='l',
     main='phase portrait', xlab='x', ylab = 'dx/dt')
# Edition of the equations
visuEq(nVar, dMax, reconstr$KL)

######################
# Example 3 ##########
######################

# for generating a model perturbed by a muliplicative noise
# The model dimension
nVar = 3
# maximal polynomial degree
dMax = 2
# Rossler Model formulation
a = 0.520
b = 2
c = 4
Eq1 <- c(0,-1, 0,-1, 0, 0, 0, 0, 0, 0)
Eq2 <- c(0, 0, 0, a, 0, 0, 1, 0, 0, 0)
Eq3 <- c(b,-c, 0, 0, 0, 0, 0, 1, 0, 0)
K <- cbind(Eq1, Eq2, Eq3)

# Edit the equations
visuEq(nVar, dMax, K)

# initial conditions
v0 <- c(-0.6, 0.6, 0.4)

# output time required
timeOut = (0:5000)/50

# variance of additive noise
varBruitA = c(0,0,0)^2

# variance of multiplitive noise
varBruitM = c(2E-2, 0, 2E-2)^2

# numerical integration with noise
intgr <- numinoisy(v0, timeOut, K, varBruitA = varBruitA, varBruitM = varBruitM, freq = 1)

# Plot of the simulated time series obtained
dev.new()
plot(intgr$donnees[,2], intgr$donnees[,3], type='l',
     main='phase portrait', xlab='x(t)', ylab = 'y(t)')
dev.new()
par(mfrow = c(3, 1))
plot(intgr$donnees[,1], intgr$donnees[,2], type='l',
     main='phase portrait', xlab='x(t)', ylab = 'y(t)')
lines(intgr$donnees[,1], intgr$vectBruitM[,2]*10, type='l',
      main='phase portrait', xlab='x(t)', ylab = 'e(t)*10', col='red')
plot(intgr$donnees[,1], intgr$donnees[,3], type='l',
     main='phase portrait', xlab='x(t)', ylab = 'y(t)')
lines(intgr$donnees[,1], intgr$vectBruitM[,3]*10, type='l',
      main='phase portrait', xlab='x(t)', ylab = 'e(t)*10', col='red')
plot(intgr$donnees[,1], intgr$donnees[,4], type='l',
     main='phase portrait', xlab='x(t)', ylab = 'y(t)')
lines(intgr$donnees[,1], intgr$vectBruitM[,4]*10, type='l',
      main='phase portrait', xlab='x(t)', ylab = 'e(t)*10', col='red')
})
