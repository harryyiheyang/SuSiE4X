test_that("Gaussian fit stops on the third consecutive no-CS iteration", {
  set.seed(20260719)
  n <- 100L
  p <- 6L
  Q <- qr.Q(qr(cbind(Intercept = 1, matrix(rnorm(n * (p + 1L)), nrow = n))))
  X <- Q[, 2:(p + 1L), drop = FALSE] * sqrt(n)
  y <- Q[, p + 2L] * sqrt(n)
  colnames(X) <- paste0("X", seq_len(p))

  susie_para <- list(coverage = 0.99, min_abs_corr = 1, max_iter = 20L)
  fit <- SuSiE4I(
    X = X, y = y, family = gaussian(), scale_data = FALSE,
    L_main = 2L, L_int = 1L, include_x_squared = TRUE,
    susie_para_main = susie_para, susie_para_int = susie_para,
    max_iter = 8L, min_iter = 6L, max_eps = 0,
    n_threads = 1L, verbose = FALSE
  )

  expect_equal(fit$diagnostics$iterations, 3L)
  expect_null(fit$main_discoveries)
  expect_s3_class(fit$fitX, "susie")
  expect_false(is.null(fit$fitJoint))
})
