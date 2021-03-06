% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cano2M.R
\name{cano2M}
\alias{cano2M}
\title{cano2M : Converts a model in canonical form into a matrix form}
\usage{
cano2M(nVar, dMax, poly)
}
\arguments{
\item{nVar}{The number of variables}

\item{dMax}{The maximum degree allowed in the formulation}

\item{poly}{A vector of coefficients corresponding to the regressor
of the canonical function}
}
\description{
Converts the vectorial formulation of canonical models
into a matrix formulation (that is, including explicitely all the
equations). For both input, the list of terms follows the convention
defined by \code{poLabs}.
}
\examples{
# A vector of polynomial terms corresponding to a canonical form:
polyTerms <- c(0.2,0,-1,0.5,0,0,0,0,0,0)
# Convert this vector into a matrix formulation with all the equations: 
K <- cano2M(3,2,polyTerms)
# Visualize the equations:
visuEq(K,3,2)

}
\seealso{
\code{\link{drvSucc}}, \code{\link{gPoMo}}, \code{\link{poLabs}}
}
\author{
Sylvain Mangiarotti, Mireille Huc
}
