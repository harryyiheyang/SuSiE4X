#' Blockwise crossproduct for large dense matrices
#'
#' Computes X^T X by row blocks to reduce memory usage.
#'
#' @param X A numeric matrix.
#' @param n_threads Number of threads for parallel computation (default = 1).
#' @param block_size Number of rows per block (default = 10000).
#'
#' @return A p × p matrix equal to X^T X.
#' @export
blockwise_crossprod <- function(X, n_threads = 1L, block_size = 10000L) {
  .Call('_SuSiE4X_blockwise_crossprod', PACKAGE = 'SuSiE4X', X, n_threads, block_size)
}

