#' Numerical description of a list of eighteen three-dimensional chaotic
#' sytems (see vignette \code{VII_Retro-Modelling})
#'
#' @name allMod_nVar3_dMax2 data set
#' @docType data
#'
#' @encoding UTF-8
#'
#' @description A list named \code{allMod_nVar3_dMax2} of matrix
#' providing the numerical description of eighteen three-dimensional
#' chaotic systems: \cr
#' Lorenz-1963 (\code{$L63}), Rössler-1976 (\code{$R76}), Burke & shaw 1981 (\code{$BS81}),
#' Lorenz-1984 (\code{$L84}), Nosé & Hooer 1986 (\code{$NH86}), Genesio & Tosi 1992 (\code{$GT92}),
#' Spott systems 1994 (\code{$SprF}, \code{$SprH}, \code{$SprK}, \code{$SprO}, \code{$SprP}, \code{$SprG},
#' \code{$SprM}, \code{$SprQ}, \code{$SprS}),
#' Chlouverakis & Sprott 2004 (\code{$CS2004}), Li 2007 (\code{$Li2007})
#' and the Cord system by Aguirre & Letellier 2012 (\code{$Cord2012}).
#' Each dynamical system is provided as a matrix:
#' each column corresponds to one equation,
#' each lines to the polynomial coefficients
#' which order is following the convetion defined
#' by function \code{poLabs(nVar = 3, dMax = 2)}.
#'
#' @author Sylvain Mangiarotti, Mireille Huc.
#'
#' @references All the references are provided in
#' vignette \code{VII_retro-modelling}.
#'
#' @keywords data
"allMod_nVar3_dMax2"
