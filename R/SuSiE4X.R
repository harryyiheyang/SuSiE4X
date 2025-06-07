#' SuSiE4X: Iterative SuSiE fitting for main and interaction effects (dense matrix only)
#'
#' This function performs iterative estimation of additive (Z) and main (X) effects, as well as
#' interaction effects between selected X components and Z, using SuSiE on sufficient statistics.
#'
#' For large sample size, it uses a blockwise strategy (like in `bam`) to compute `crossprod(X)` efficiently
#' without requiring the full X'X in memory at once.
#'
#' The algorithm alternates between updating the effects from Z (via linear projection),
#' from X (via susie_suff_stat), and interaction terms (via susie_suff_stat after constructing selected interactions).
#' After convergence, a post-hoc linear model of y ~ Z is fit with offset(etaX + etaW) to estimate
#' the significance of Z effects adjusted for X and interactions.
#'
#' @param X An n × p matrix of covariates.
#' @param Z An n × q matrix of covariates. If Z = NULL, only GxG interaction will be consider.
#' @param y An n-vector of responses.
#' @param n_threads Integer. Number of threads used in parallel computations, such as blockwise matrix crossproduct or matrix–vector multiplication. Default is 1 (no parallelism). Increasing this can speed up computation on multicore machines when \code{is.large = TRUE}.
#' @param Lmain Number of SuSiE effects for main effect X (default = 5).
#' @param Linteraction Number of SuSiE effects for interactions (default = 5).
#' @param max.iter Maximum number of outer iterations (default = 10).
#' @param max.eps Convergence threshold on η (default = 1e-5).
#' @param min.iter Minimum number of outer iterations (default = 3).
#' @param susie.iter Maximum number of iterations for each SuSiE fit (default = 300).
#' @param verbose A logical indicator of whether to display iteration diagnostics (default = TRUE).
#' @param ... Additional arguments passed to \code{susie_suff_stat}, such as \code{prior_variance}, \code{null_weight}, \code{standardize}, etc. Note that these tuning parameters will only be used in the fine-mapping of marginal effects.
#'
#' @return A list containing results from the iterative SuSiE fitting procedure:
#' \describe{
#'   \item{iter}{Total number of outer iterations completed before convergence or reaching maximum.}
#'   \item{error}{Numeric vector of max change in fitted values across iterations.}
#'   \item{fitX}{SuSiE result for main effects.}
#'   \item{fitW}{SuSiE result for interaction terms.}
#'   \item{fitZ}{Final linear model of y ~ Z with offset.}
#' }
#'
#' @importFrom Matrix colMeans crossprod
#' @importFrom stats var lm
#' @importFrom susieR susie_suff_stat coef.susie
#' @importFrom CppMatrix matrixMultiply matrixVectorMultiply
#' @importFrom graphics text
#' @importFrom stats lm coef
#'
#' @export
SuSiE4X <- function(X, Z=NULL, y, crossprodX=NULL,
                    n_threads = 1, Lmain = 5, Linteraction = 5,
                    max.iter = 15, max.eps = 1e-5, min.iter = 3,
                    susie.iter = 300, verbose = TRUE,...) {

if (is.null(colnames(X))){
colnames(X) <- paste0("X", seq_len(ncol(X)))
}

if(is.null(Z)==F){
nameZ=colnames(Z)
if (is.null(colnames(Z))){
nameZ <- paste0("Z", seq_len(ncol(Z)))
}
Z=demean(Z)

Run_GGE(X = X, Z = Z, y = y, crossprodX=crossprodX,
                Lmain = Lmain, Linteraction = Linteraction,
                max.iter = max.iter, max.eps = max.eps, min.iter = min.iter,
                susie.iter = susie.iter, verbose = verbose,n_threads=n_threads,...)
}else{
Run_GG(X = X, y = y,
        Lmain = Lmain, Linteraction = Linteraction,
        max.iter = max.iter, max.eps = max.eps, min.iter = min.iter,
        susie.iter = susie.iter, verbose = verbose,n_threads=n_threads,...)
}

}
