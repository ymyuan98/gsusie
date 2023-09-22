#' @title Auxilary computation functions
#' 
#' @description This file defines the auxilary computation functions. 
#' 
#' @details \code{compute_colstats} sets three attributes for matrix X: 
#' \code{attr(X,'scaled:center')} is a p-vector of column means of X if 
#' \code{center=TRUE}, a p-dim vector of zeros otherwise; \code{'attr(X,'scaled:scale')} 
#' is a p-dim vector of column standard deviations of X if \code{scale=TRUE}, 
#' a p-dim vector of ones otherwise; \code{'attr(X,'d')} is a p-dim
#' vector of column sums of \code{X.standardized^2}, where \code{X.standardized} 
#' is the matrix X centered by \code{attr(X,'scaled:center')} and scaled 
#' by \code{attr(X,'scaled:scale')}.
#' 
#' \code{compute_Xb} offers the matrix multiplication of 
#' an (n x p) matrix \code{X} and a (p x 1) array/matrix {b}, \ie \eqn{Xb}.  
#' The output is a (n x 1) matrix. 
#' 
#' \code{compute_XtWY} offers the matrix multiplication of \eqn{X'WY} 

require(Matrix)

#' @keywords internal
#' 
compute_colstats <- function(X, center = FALSE, scale = TRUE) {
    out <- list()

    if (!center) {
        out$cm <- colMeans(X)
    } else {
        out$cm <- rep(0, times = ncol(X))
    }

    if (scale) {
        out$csd <- apply(X, 2, sd)
        # X.standardized <- scale(X)
    } else {
        out$csd <- rep(1, times = ncol(X))
        # X.standardized <- X
    }

    # out$d <- colSums(X.standardized^2)
    return(out)
}

#' @keywords internal
#'
compute_Xb <- function(X, b) {
    if (!is.matrix(X))
        stop("Input X should be a matrix")
    if (!is.matrix(b)) {
        b <- as.matrix(as.array(b))
        if (dim(b)[2] != 1)
            stop("Input b should be a (p x 1) matrix")
    }
    if (dim(X)[2] != dim(b)[1])
        stop("2nd dim of X and 1st dim of b do not match")

    return(X %*% b)
}

#' #' @keywords internal
#' #'
#' compute_XtWX <- function(X, W) { 
#'     # Simplify calculation since W is a diagnoal matrix
#'     if (!is.matrix(X))
#'         stop("Input X should be a matrix")
#'     if (!is.matrix(W))
#'         stop("Input W should be a matrix")
#'     if (!isSymmetric(W))
#'         stop("Input W should be a symmetric square matrix")
#'     if (dim(X)[1] != dim(W)[1])
#'         stop("1st dim of X and dim of W do not match")
#' 
#'     return(t(X) %*% W %*% X)
#' }

#' #' @keywords internal
#' #'
#' compute_XtWY <- function(X, W, Y) {
#'     if (!is.matrix(X))
#'         stop("Input X should be a matrix")
#'     if (!is.matrix(W))
#'         stop("Input W should be a matrix")
#'     if (!isSymmetric(W))
#'         stop("Input W should be a symmetric square matrix")
#'     if (!is.matrix(Y)) {
#'         Y <- as.matrix(as.array(Y))
#'         if (dim(Y)[2] != 1)
#'         stop("Input Y should be a (n x 1) matrix")
#'     }
#' 
#'     if (dim(X)[1] != dim(W)[1])
#'         stop("1st dim of X and dim of W do not match")
#'     if (dim(W)[1] != dim(Y)[1])
#'         stop("Dim of W and 1st dim of Y do not match")
#'     
#'     return(t(X) %*% W %*% Y)
#' }
