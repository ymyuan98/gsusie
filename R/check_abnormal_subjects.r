#' @title Sensitivity check
#'
#' @description
#'
#' In iterative weighted least square methods,
#' each data point is transformed into a corresponding pseudo-response with
#' a corresponding pseudo-variance.
#' These transformations would inevitably point to some abnormal subjects
#' whose pseudo-response or pseudo-variance is an abnormal value,
#' such as Inf/-Inf and NaN(?).
#'
#' We take these subjects with abnormal pseudo-response/pseudo-variance
#' as ``abnormal subjects''.
#' \code{check_abnormal_subjects()} return the indices of abnormal
#' subjects (if any).
#'

check_abnormal_subjects <- function(values) {
    # check if any infinite
    abn_idx <- which(is.infinite(values) | is.nan(values))
    if (length(abn_idx) == 0) abn_idx <- NULL

    return(abn_idx = abn_idx)
}
