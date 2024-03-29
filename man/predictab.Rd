% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predictab.R
\name{predictab}
\alias{predictab}
\title{Estimate the models performance obtained with \code{GPoMo}
in term of predictability}
\usage{
predictab(
  ogp,
  fullt = NULL,
  fulldata = NULL,
  hp = NULL,
  Nech = 50,
  intSimStep = NULL,
  show = 1,
  selecmod = NULL,
  id = 1,
  selV = 1,
  na.rm = FALSE
)
}
\arguments{
\item{ogp}{The output list obtained from function \code{gPoMo}.}

\item{fullt}{Time vector of the data set for which predictability
will be tested}

\item{fulldata}{Data set for which predictability will be tested}

\item{hp}{Time vector of the horizon of prediction}

\item{Nech}{Number of simulations}

\item{intSimStep}{Internal number of simulation steps}

\item{show}{Provide (2) or not (0-1) visual output during
the running process.}

\item{selecmod}{A vector of the model selected.}

\item{id}{The type of model to identify. \code{id = 1} corresponds
to unidentified models, that is, potentialy chaotic.}

\item{selV}{Selected variable for the analysis}

\item{na.rm}{Indicates if the \code{NA} should be removed
(\code{na.rm = TRUE}) or not (\code{na.rm = FALSE}).}
}
\value{
\code{ErrmodAll} A list of matrix \code{$Predmod1},
\code{$Predmod2}, etc. and \code{$Errmod1}, \code{$Errmod2}, etc.
providing respectively the forecasting and the forecasting error
of models 1, 2, etc.
Each column corresponds to one simulation starting from
a specific initial condition. Each line corresponds to
one horizon of prediction.
Vectors corresponding to the initial condition time \code{tE}
and the horizon of prediction \code{hpE} are also provided
in \code{$tE} and \code{$hpE}, respectively.
The percentiles of the distributions of error growth
are provided in \code{qt} (0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)
and of absolute error growth in \code{qt2} (0.5, 0.75, 0.9, 0.95, 0.98, 0.99).
}
\description{
The algorithm aims to estimate automatically the forecasting
performances of the models obtained with \code{gPoMo}.
}
\examples{
# load data
data("Ross76")
# time vector
tin <- Ross76[seq(1, 3000, by = 8), 1]
# single time series
data <- Ross76[seq(1, 3000, by = 8), 3]
# dev.new()
# plot(tin, data, xlab = 'time', ylab = 'y(t)')

# global modelling
# results are put in list outputGPoM
outputGPoM <- gPoMo(data[1:300], tin = tin[1:300], dMax = 2, nS=c(3),
                    show = 0, method = 'rk4',
                    nPmin = 10, nPmax = 12,
                    IstepMin = 150, IstepMax = 151)
#
visuOutGP(outputGPoM)

###########################
# and test predictability #
###########################
outpred <- predictab(outputGPoM, hp = 15, Nech = 30)

# manual visualisation of the outputs (e.g. for model 1):
dev.new()
image(outpred$tE, outpred$hpE, t(outpred$Errmod1),
xlab = 't', ylab = 'hp', main = 'Errmod1')

}
\author{
Sylvain Mangiarotti, Mireille Huc
}
