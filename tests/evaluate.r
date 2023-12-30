#' @title Evaluation of simulation results
#' @description This Rscript contains functions for evaluating the simulation
#' results, including false discovery rate, discovery rate, sensitivity,
#' and specificity.

#' evaluate true positive rate (TPR) / sensitivity / recall
evaluate_tpr <- function(test_result, ground_truth) {
  if (is.logical(test_result))  test_result <- which(test_result)
  if (is.logical(ground_truth)) ground_truth <- which(ground_truth)

  tp <- test_result %in% ground_truth
  tpr <- sum(tp) / length(ground_truth)
  return(tpr)
}

#' evaluate precision / positive predictive value (PPV)
evaluate_ppv <- function(test_result, ground_truth) {
  if (is.logical(test_result))  test_result <- which(test_result)
  if (is.logical(ground_truth)) ground_truth <- which(ground_truth)

  tp <- test_result %in% ground_truth
  ppv <- sum(tp) / length(test_result)
  return(ppv)
}

#' evaluate false discovery rate (FDR)
evaluate_fdr <- function(test_result, ground_truth) {
  fdr <- 1 - evaluate_ppv(test_result, ground_truth)
  return(fdr)
}

#' evaluate accuracy
evaluate_acc <- function(test_result, ground_truth) {
  if (!(is.logical(test_result) && is.logical(ground_truth))) {
    stop("Either the input test result or the ground truth is not logical")
  }
  if (length(test_result) != length(ground_truth)) {
    stop("The lengths of input test result and ground truth are not the same")
  }

  numerator <- sum(test_result == ground_truth)
  denominator <- length(ground_truth)
  return(numerator / denominator)
}
