#' This is data to be included in my package
#'
#' @name allToTest
#' @docType data
#' @author Sylvain Mangiarotti, Mireille Huc
#' @description List of 6 models to be tested by autoGPoMoTest
#' Each model $mToTest1, $mToTest2, etc. is provided as a matrix
#' of coefficients of dimension 10 * 3 that corresponds to a
#' formulation of three equations (nVar = 3) of maximal polynomial
#' degree dMax = 2. Each column column corresponds to one equation.
#' The order of the terms is given by \code{poLabs(nVar = 3, dMax = 2)}.
#' @examples
#' ###########
#' # example #
#' ###########
#' data("allToTest")
#' # 6 models are available in this list:
#' names(allToTest)
#' # The parameter of their formulation (nVar and dMax)
#' # can be retrieved:
#' nVar <- dim(allToTest$mToTest6)[2]
#' dMax <- p2dMax(nVar = 3, pMaxKnown = dim(allToTest$mToTest6)[1])
#' # Their equation can be edited as follows:
#' visuEq(nVar, dMax, allToTest$mToTest6, approx = 2)
#'
#' @keywords data
"allToTest"
