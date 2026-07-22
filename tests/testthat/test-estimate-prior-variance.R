test_that("estimate_prior_variance returns the two prior scales", {
  expect_identical(formals(estimate_prior_variance)$max_iter, 100L)
  expect_identical(formals(estimate_prior_variance)$tol, 1e-6)
  expect_identical(formals(estimate_prior_variance)$r_thres, 0)

  set.seed(11)
  n <- 2000
  p <- 20
  X <- scale(matrix(rnorm(n * p), n, p))
  colnames(X) <- paste0("rs", seq_len(p))
  Z <- scale(matrix(rnorm(n * 2), n, 2))
  PRS <- matrix(rnorm(n * 22), n, 22)
  colnames(PRS) <- as.character(seq_len(22))
  SNPInfo <- data.frame(
    SNP = colnames(X),
    CHR = rep(c("1", "2"), each = p / 2)
  )
  b <- rnorm(p, sd = sqrt(0.02 / p))
  y <- as.numeric(Z %*% c(0.4, -0.2) + rowSums(PRS) * 0.02 + X %*% b +
                    rnorm(n))

  out <- estimate_prior_variance(
    X = X, y = y, Z = Z, PRS = PRS, SNPInfo = SNPInfo,
    family = gaussian(), n_threads = 1
  )

  expect_identical(
    names(out), c("prior_variance", "scaled_prior_variance", "diagnostics")
  )
  expect_true(is.finite(out$prior_variance) && out$prior_variance > 0)
  expect_equal(
    out$scaled_prior_variance,
    out$prior_variance / out$diagnostics$var_y_ss
  )
  expect_s3_class(out$diagnostics, "data.frame")
  expect_equal(nrow(out$diagnostics), 1)
  expect_equal(out$diagnostics$r_thres, 0)

  out_correlated <- estimate_prior_variance(
    X = X, y = y, Z = Z, PRS = PRS, SNPInfo = SNPInfo,
    family = gaussian(), n_threads = 1, r_thres = 0.01
  )
  expect_true(is.finite(out_correlated$prior_variance))
  expect_gt(out_correlated$prior_variance, 0)
  expect_equal(out_correlated$diagnostics$r_thres, 0.01)
  expect_lte(out_correlated$diagnostics$xtx_density, 0.5)

  out_no_prs <- estimate_prior_variance(
    X = X, y = y, Z = Z, PRS = NULL, SNPInfo = SNPInfo,
    family = gaussian(), n_threads = 1
  )
  out_no_prs_correlated <- estimate_prior_variance(
    X = X, y = y, Z = Z, PRS = NULL, SNPInfo = SNPInfo,
    family = gaussian(), n_threads = 1, r_thres = 0.01
  )
  expect_true(is.finite(out_no_prs$prior_variance))
  expect_true(is.finite(out_no_prs_correlated$prior_variance))
  expect_equal(out_no_prs_correlated$diagnostics$r_thres, 0.01)
})

test_that("estimate_prior_variance requires exact SNP and chromosome matches", {
  X <- matrix(rnorm(400), 100, 4)
  colnames(X) <- paste0("rs", 1:4)
  PRS <- matrix(rnorm(2200), 100, 22)
  colnames(PRS) <- as.character(seq_len(22))
  SNPInfo <- data.frame(SNP = colnames(X), CHR = c("1", "1", "2", "X"))
  Z <- matrix(rnorm(200), 100, 2)

  expect_error(
    estimate_prior_variance(X, rnorm(100), Z, PRS, SNPInfo, gaussian()),
    "exactly match a PRS column"
  )

  SNPInfo$CHR[4] <- "2"
  expect_error(
    estimate_prior_variance(
      X, rnorm(100), Z, PRS, SNPInfo, gaussian(), r_thres = 1
    ),
    "r_thres"
  )
  expect_error(
    estimate_prior_variance(
      X, rnorm(100), NULL, PRS, SNPInfo, gaussian()
    ),
    "Z must be a non-empty"
  )
})
