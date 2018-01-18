#' @title Weighted inner product
#'
#' @description Computes weighted inner products.
#'
#' @param A1 The input matrix 1.
#' @param A2 The input matrix 2.
#' @param weight The weighting vector.
#'
#' @return \code{inP} The weighted inner product.
#'
#' @examples
#' ############
#' #Example 1 #
#' ############
#' A1 = c(0,1,2,0,1,3)
#' A2 = c(1,2,0,0,4,1)
#' wInProd(A1, A2)
#'
#' ############
#' #Example 2 #
#' ############
#' A1 = c(0,1,2,0,1,3)
#' A2 = c(1,2,0,0,4,1)
#' w = c(0,0,0,1,1,1)
#' wInProd(A1, A2, weight = w)
#'
#' @export
wInProd <- function(A1, A2, weight=NULL) {
  # by default, that is if weight is unspecified (weight = NULL)
  # then the weight is taken equal to 1 everywhere
  if (is.null(weight)) weight <- matrix(1, ncol = 1, nrow = length(A2))

  # weighting
  inP <- (as.matrix(A2) * weight)
  # inner product
  inP <- t(as.matrix(A1)) %*% inP
  # return
  inP
}
