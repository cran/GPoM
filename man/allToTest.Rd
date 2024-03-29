% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/allToTest-data.R
\docType{data}
\name{allToTest}
\alias{allToTest}
\title{A list providing the description of six models tested by
the function \code{autoGPoMoTest}.}
\format{
An object of class \code{list} of length 6.
}
\usage{
allToTest
}
\description{
List of 6 models available for tests (by \code{autoGPoMoTest}).
Each model (\code{$mToTest1}, \code{$mToTest2}, etc.) is provided as a matrix
of dimension 10 * 3. Each column corresponds to one equation.
The order of the coefficients follows the conventions defined
by \code{poLabs(nVar = 3, dMax = 2)}.
}
\examples{
###########
# example #
###########
data("allToTest")
# 6 models are available in this list:
names(allToTest)
# The parameter of their formulation (nVar and dMax)
# can be retrieved:
nVar <- dim(allToTest$mToTest6)[2]
dMax <- p2dMax(nVar = 3, pMaxKnown = dim(allToTest$mToTest6)[1])
# Their equation can be edited as follows:
visuEq(allToTest$mToTest6, nVar, dMax, approx = 2)

}
\author{
Sylvain Mangiarotti, Mireille Huc
}
\keyword{data}
