#' @title Auxilary computation functions
#'
#' @description This file defines the auxiliary computation functions.
#'
#' @details
#' \code{compute_colstats} returns a list with two attributes:
#' \code{out$cm} is a p-vector of column-wise means of X if
#' \code{center=TRUE}, a p-dim vector of zeros otherwise;
#' \code{out$csd} is a p-dim vector of column-wise standard
#' deviations of X if \code{scale=TRUE}, a p-dim vector of ones otherwise;
#'
#' \code{compute_Xb} offers the matrix multiplication of
#' an (n x p) matrix \code{X} and a (p x 1) array/matrix {b}.
#' The output is a (n x 1) matrix.
#'
#' @importFrom stats sd
#'
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
    } else {
        out$csd <- rep(1, times = ncol(X))
    }

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
