#' @title Aggreagate fitted gsusie results
#' @description aggregate fitted gsusie results especially
#' when choosing the grid optimization (?) for the coefficient prior variance
#' to refine the fitted results.
#' 
#' If \code{length(gs) == 1}, indicating that we fix one coefficient prior
#' variance, the output is \code{gs} itself (removing the outer list number). 
#' Otherwise, we have \code{length(gs)} sets of gsusie results
#' fitted under various \code{grid_prior_variance_value}.
#' If \code{length(gs) > 1}, indicating that we have multiple gsusie results
#' yielded from multiple coefficient prior variance, 
#' the output is the best gsusie fit with the highest ELBO. 
#' 
#' NoNoNo...
#' If \code{method = "best"}, the output is the gsusie fit with the 
#' highest ELBO among \code{length(gs)} sets of gsusie fits; 
#' if \code{method = "weighted_sum"}, the output is the sum 
#' of \code{length(gs)} sets of gsusie fits weighted by their ELBOs. 
#'

aggregate_gsusie_results <- function(gs) {
# aggregate_gsusie_results <- function(gs, method = c("best", "weighted_sum")) {

    # return the result if length(gs) == 1
    if (length(gs) == 1) {
        out <- gs[[1]]
        return(out)
    }

    # method <- match.arg(method, c("best", "weighted_sum"))
    elbos <- sapply(gs, function(x) {x$elbo[x$niter]})
    # if (method == "best") {
        best_case <- which.max(elbos)
        out <- gs[[best_case]]
        # out$aggregate_method <- "best"
    # } else {  # if method == "weighted_sum"
    #     ws <- elbos
    # }

    # class(out) <- "gsusie"  ##?

    return(out)
}


    # out <- list()
    # # extract the shared elements
    # for (name in c("pie", "family", "null_index")) {
    #     out[[name]] <- gs[[1]][[name]]
    # }

    # # aggregate alpha
    # out$alpha <- matrix(0, nrow = nrow(gs[[1]]$alpha),
    #                        ncol = ncol(gs[[1]]$alpha))
    # out$sigma02 <- rep(NA, times = length(gs))
    # denom <- sum(sapply(gs, function(x){x$converged == T}))
    # for (ls_idx in 1 : length(gs)) {
    #     if (gs[[ls_idx]]$converged) {
    #         out$alpha <- out$alpha + gs[[ls_idx]]$alpha
    #         out$mu    <- out$mu    + gs[[ls_idx]]$mu
    #         out$mu2   <- out$mu2   + gs[[ls_idx]]$mu2
    #     } else {
    #         warning(paste("Drop the fitted alpha under coefficient",
    #                       "prior variance of", gs[[ls_idx]]$sigma02))
    #     }
    #     out$sigma02[ls_idx] <- gs[[ls_idx]]$sigma02
    # }
    # out$alpha <- out$alpha / denom
    # out$mu    <- out$mu    / denom
    # out$mu2   <- out$mu2   / denom

    # # for the rest of the elements...
    # for (ls_idx in 1 : length(gs)) {
    #     gs_parts <- list()
    #     for (name in c("sigma02", "alpha", "mu", "mu2", "betahat",
    #                    "lbf", "lbf_variable", "converged", "abn_subjects",
    #                    "elbo", "loglik_exact", "loglik_apprx", "niter",
    #                    "null_index")) {
    #         gs_parts[[name]] <- gs[[ls_idx]][[name]]
    #         class(gs_parts) <- "gsusie"
    #     }
    # out <- append(out, list(gs_parts))
    # rm(gs_parts)
    # }
    
    # names(out)[(length(out)-length(gs)+1) : length(out)] <- 
    #     paste0("case", 1 : length(gs))
