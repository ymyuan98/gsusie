#' @title Remove abnormal points from the current matrix/array
#'
#' @param abn_index an array of indices of abnormal points.
#'
#' @param object a matrix or a vector that may contain abnormal points.
#' If \code{object} is a matrix, then \code{abn_index} indicates
#' the rows to be removed.
#'
#' @keywords internal
remove_abnormal <- function(abn_index, object) {
  if (!(is.matrix(object) | is.vector(object)))
    stop("Input object should be a matrix or a vector")

  abn_index <- unique(abn_index)

  if (is.matrix(object)) {
    len <- nrow(object)
  } else { # is.vector(object)
    len <- length(object)
  }

  if (!is.null(abn_index)) {
    sub_index <- (!(1:len) %in% abn_index)
  } else {
    sub_index <- rep(TRUE, times = len)
  }

  if (is.matrix(object)) {
    return(object[sub_index, , drop = F])
  } else {  # is.vector(object)
    return(object[sub_index])
  }

}
