% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/combiEq.R
\name{combiEq}
\alias{combiEq}
\title{combiEq : Combine Equations from different sources}
\usage{
combiEq(allKL, eqOrder = NULL, dMaxOut = NULL)
}
\arguments{
\item{allKL}{A list of models, each provided as a matrix.
A single matrix can also be provided, it will be transformed
into a list containing a single matrix.}

\item{eqOrder}{A list of vector, providing each the equations
number (relating to the input models) to be kept in the
output equation system. If not provided, all the equations
are kept.
A single matrix can also be provided, it will be transformed
into a list containing a single matrix.}

\item{dMaxOut}{The maximal polynomial degree of the output
equation system (if not provided, this degree is deduced
from the input models)}
}
\value{
KLout A matrix of the combined model
}
\description{
Combines equations of different sources
into a single system. During this combination, the polynomial
maximal degree can be either imposed or optimized to reduce
the model size. All the input have to follow
the convention defined by \code{poLabs}.
}
\examples{
# Load models
data("allMod_nVar3_dMax2")
# Display equations of system 1
visuEq(nVar = 3, dMax = 2, K = allMod_nVar3_dMax2$NH86, substit = 1)
# Display equations of system 2
visuEq(nVar = 3, dMax = 2, K = allMod_nVar3_dMax2$R76, substit = 1)
# put the two systems in a list
allK <- list()
allK[[1]] <- allMod_nVar3_dMax2$NH86
allK[[2]] <- allMod_nVar3_dMax2$R76

# Example 1: reformulate two autonomous system in a single matrix 
visuEq(K = allK[[1]], substit = c('u', 'v', 'w'))
visuEq(K = allK[[2]], substit = c('X', 'Y', 'Z'))
Knew <- combiEq(allK)
visuEq(K = Knew, substit = c('u', 'v', 'w', 'X', 'Y', 'Z'))

# Example 2
inXnote = list()
inXnote[[1]] <- c('u', 'v', 'w')
inXnote[[2]] <- c('X', 'Y', 'Z')
visuEq(K = allK[[1]], substit = inXnote[[1]])
visuEq(K = allK[[2]], substit = inXnote[[2]])
XnoteOut = c('X', 'Y', 'Z', 'u', 'v', 'w')
Knew2 <- combiEq(allK, eqOrder=c(4,5,6,1,2,3))
visuEq(K = Knew2, substit = XnoteOut)

# Example 3
inXnote = list()
inXnote[[1]] <- c('u', 'v', 'w')
inXnote[[2]] <- c('X', 'Y', 'Z')
visuEq(K = allK[[1]], substit = inXnote[[1]])
visuEq(K = allK[[2]], substit = inXnote[[2]])
XnoteOut = c('u', 'X', 'v', 'Y', 'w', 'Z')
Knew3 <- combiEq(allK, eqOrder=c(1,4,2,5,3,6), dMaxOut = 3)
visuEq(K = Knew3, substit = XnoteOut)
}
\seealso{
\code{\link{gPoMo}}, \code{\link{poLabs}}
}
\author{
Sylvain Mangiarotti
}
