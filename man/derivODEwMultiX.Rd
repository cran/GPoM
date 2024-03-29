% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/derivODEwMultiX.R
\name{derivODEwMultiX}
\alias{derivODEwMultiX}
\title{deriveODEwMultiX : A Subfonction for the numerical integration
of polynomial equations in the generic form defined by function
\code{poLabs} and with External Forcing F(t)}
\usage{
derivODEwMultiX(t, x, K, extF, regS = NULL)
}
\arguments{
\item{t}{All the dates for which the result of the numerical integration
of the model will have to be provided}

\item{x}{Current state vector (input from which the next state will
be estimated)}

\item{K}{is the model: each column corresponds to one equation
which organisation is following the convention given by function
\code{poLabs} which requires the definition of the model dimension
\code{nVar} (i.e. the number of variables) and the maximum polynomial
degree \code{dMax} allowed. The last Equation correspond to the
forcing variable that is artificially set to 0.}

\item{extF}{is the external forcing. It is defined by two columns.
The first colomn correspond to time t. The second column to F(t) the
forcing at time t. Note that when launching the integration function
ode, the forcing F(t) should be provided with a sampling time twice
the sampling time used in t (because rk4 method will always use an
intermediate time step).}

\item{regS}{Current states of each polynomial terms used
in \code{poLabs}. These states can be deduced from the current
state vector x (using function \code{regSeries}). When available,
it can be provided as an input to avoid unecessary computation.}
}
\value{
xxx
}
\description{
This function provides the one step integration of
polynomial Ordinary Differential Equations (ODE). This function
requires the function \code{ode} ("deSolve" package).
This function has to be run with the Runge-Kutta method (method = 'rk4')
}
\examples{
# build a non autonomous model
nVar = 4
dMax = 3
omega = 0.2
gamma = 0.05
KDf=matrix(0, nrow = d2pMax(nVar = nVar, dMax = dMax), ncol = nVar)
KDf[11,1]  = 1
KDf[2,2]  = 1
KDf[5,2]  = 1
KDf[11,2]  = -gamma
KDf[35,2] = -1
KDf[2,3]  = NA
KDf[2,4]  = NA
visuEq(K = KDf, substit = c('x', 'y', 'u', 'v'))
#
# Prepare the external forcing
# number of integration time step
Istep <- 500
# time step
smpl <- 1 / 20
# output time vector
dater <- (0:Istep) * smpl
# hald step time vector (for Runge-Kutta integration)
daterdbl <- (0:(Istep*2 + 1)) * smpl / 2
# generate the forcing (here variables u and v)
extF = cbind(daterdbl, -0.1 * cos(daterdbl * omega), 0.05 * cos(daterdbl * 16/3*omega))
#
# Initial conditions to be used (external variables can be set to 0)
etatInit <- c(-0.616109362 , -0.126882584 , 0, 0)
#
# Numerical integration
reconstr2 <- ode(etatInit, dater, derivODEwMultiX,
                 KDf, extF = extF, method = 'rk4')
# Reconstruction of the output
nVarExt <- dim(extF)[2] - 1
reconstr2[,(nVar - nVarExt + 2):(nVar + 1)] <- extF[(0:Istep+1)*2, 2:(nVarExt+1)]

}
\author{
Sylvain Mangiarotti
}
