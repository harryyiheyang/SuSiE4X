#' Crossprod with direct single-thread and blockwise parallel paths
#'
#' Computes `crossprod(X)` or `crossprod(X, Z)`. When `n_threads <= 1`,
#' it uses direct dense multiplication from CppMatrix; otherwise it uses
#' the row-blocked parallel C++ implementation.
#'
#' @param X Numeric matrix.
#' @param Z Optional second matrix. If provided, computes `crossprod(X, Z)`.
#' @param n_threads Number of threads. Defaults to 4.
#' @param block_size Row block size for the parallel path.
#' @export
blockwise_crossprod <- function(X, Z = NULL, n_threads = 4L, block_size = 10000L) {
  if (!is.matrix(X) || !is.numeric(X)) stop("X must be a numeric matrix.")
  if (!is.null(Z) && (!is.matrix(Z) || !is.numeric(Z))) {
    stop("Z must be a numeric matrix.")
  }
  if (!is.numeric(n_threads) || length(n_threads) != 1L ||
      !is.finite(n_threads) || n_threads < 1) {
    stop("n_threads must be a positive numeric scalar.")
  }
  if (!is.numeric(block_size) || length(block_size) != 1L ||
      !is.finite(block_size) || block_size < 1) {
    stop("block_size must be a positive numeric scalar.")
  }
  n_threads <- as.integer(n_threads)
  block_size <- as.integer(block_size)

  if (n_threads <= 1L) {
    if (is.null(Z)) {
      return(CppMatrix::matrixMultiply(X, X, transA = TRUE))
    }
    return(CppMatrix::matrixMultiply(X, Z, transA = TRUE))
  }

  if (is.null(Z)) {
    return(.Call(`_SuSiE4I_blockwise_crossprod`, X, n_threads, block_size))
  }
  .Call(`_SuSiE4I_blockwise_crossprod2`, X, Z, n_threads, block_size)
}
