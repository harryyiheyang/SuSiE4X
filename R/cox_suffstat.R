#' Cox sufficient statistics for one design block, projected onto the nuisance complement
#'
#' Internal helper shared by the Cox runners. Given a design block Xblk to fine-map,
#' the current linear predictor eta, and nuisance covariates Znui (an intercept
#' is always appended for the projection), it builds the Breslow observed information
#' and score statistic for the Xblk block after projecting out the nuisance columns,
#' then returns the projected XtX and Xty with n_eff equal to the number of events.
#'
#' @keywords internal
#' @noRd
cox_suffstat_block <- function(Xblk, eta, Znui, surv_time, surv_status,
                               n_threads = 1, ridge = 1e-6,
                               block_size = 10000L) {
  Xblk <- as.matrix(Xblk)
  n <- nrow(Xblk)
  p <- ncol(Xblk)

  if (is.null(Znui)) {
    ZI <- matrix(1, n, 1)
    colnames(ZI) <- "Intercept"
  } else {
    Znui <- as.matrix(Znui)
    ZI <- cbind(Intercept = 1, Znui)
  }
  q <- ncol(ZI)

  XZE <- cbind(Xblk, eta, ZI)

  ss <- cox_suffstat(X = XZE, eta = eta, time = surv_time,
                     status = as.integer(surv_status), n_threads = n_threads)
  a <- as.numeric(ss$a)
  B <- as.matrix(ss$B)
  XZEty <- as.numeric(ss$Xty)
  n_eff <- ss$d

  XZEa <- XZE * sqrt(a)
  A <- blockwise_crossprod(XZEa, n_threads = n_threads,
                           block_size = block_size)
  BtB <- blockwise_crossprod(B, n_threads = n_threads,
                             block_size = block_size)
  XZEtXZE <- A - BtB
  XZEtXZE <- (XZEtXZE + t(XZEtXZE)) / 2

  idxX <- seq_len(p)
  idxE <- p + 1L
  idxZ <- p + 1L + seq_len(q)

  XtX <- XZEtXZE[idxX, idxX, drop = FALSE]
  XtE <- XZEtXZE[idxX, idxE, drop = FALSE]
  XtZ <- XZEtXZE[idxX, idxZ, drop = FALSE]
  ZtZ <- XZEtXZE[idxZ, idxZ, drop = FALSE]
  ZtX <- XZEtXZE[idxZ, idxX, drop = FALSE]
  ZtE <- XZEtXZE[idxZ, idxE, drop = FALSE]
  XtM <- XZEty[idxX]

  Zinv_ZtX <- solve_with_ridge(ZtZ, ZtX, ridge = ridge)
  Zinv_ZtE <- solve_with_ridge(ZtZ, ZtE, ridge = ridge)

  XtX_proj <- XtX - matrixMultiply(XtZ, Zinv_ZtX)
  XtE_proj <- as.vector(XtE - matrixVectorMultiply(XtZ, Zinv_ZtE))

  Xty <- XtE_proj + XtM
  XtX <- (XtX_proj + t(XtX_proj)) / 2
  XtX_pre_ridge <- XtX
  diag(XtX) <- diag(XtX) + ridge

  dXtX <- diag(XtX)
  zhat <- Xty / sqrt(dXtX)
  R <- cov2cor(XtX)
  R <- (R + t(R)) / 2

  list(
    z = zhat, R = R, n_eff = n_eff, XtX = XtX,
    XtX_pre_ridge = XtX_pre_ridge, Xty = Xty
  )
}
