test_that("Gaussian paraPen uses phi/V and matches the ridge posterior", {
  set.seed(20260717)
  n <- 90L
  X <- cbind(x1 = rnorm(n), x2 = rnorm(n))
  phi <- 3.2
  V <- c(Main_CS1 = 0.7, Main_CS2 = 0.4)
  y <- as.numeric(X %*% c(0.8, -0.4) + rnorm(n, sd = sqrt(phi)))
  dat <- data.frame(y = y, Main_CS1 = X[, 1L], Main_CS2 = X[, 2L])

  fit <- SuSiE4I:::mgcv_fit_fixed_ridge(
    "y", names(V), dat, gaussian(), V,
    dispersion = phi, mgcv_model = "gam"
  )
  A <- cbind(Intercept = 1, X)
  closed <- solve(
    crossprod(A) / phi + diag(c(0, 1 / V), nrow = 3L),
    crossprod(A, y) / phi
  )
  pen <- attr(fit, "refit_penalty")

  expect_equal(unname(coef(fit)[names(V)]), as.numeric(closed[-1L]),
               tolerance = 1e-10)
  expect_equal(pen$working_precision, phi / V)
  expect_equal(unname(fit$full.sp), 1)

  wrong <- solve(
    crossprod(A) / phi + diag(c(0, 1 / (phi * V)), nrow = 3L),
    crossprod(A, y) / phi
  )
  expect_gt(max(abs(as.numeric(closed[-1L]) - as.numeric(wrong[-1L]))), 1e-4)
})

test_that("exact Gaussian ridge matches bam coefficients edf and dispersion", {
  src <- paste(deparse(body(SuSiE4I:::gaussian_ridge_refit)), collapse = "\n")
  expect_true(grepl("blockwise_crossprod", src, fixed = TRUE))
  expect_true(grepl("CppMatrix::matrixMultiply", src, fixed = TRUE))
  expect_true(grepl("CppMatrix::matrixSolve", src, fixed = TRUE))
  expect_true(grepl("CppMatrix::matrixVectorMultiply", src, fixed = TRUE))

  set.seed(20260720)
  n <- 500L
  Z <- matrix(rnorm(n * 2L), n, 2L)
  X <- matrix(rnorm(n * 3L), n, 3L)
  pred <- data.frame(Z1 = Z[, 1L], Z2 = Z[, 2L],
                     Main_CS1 = X[, 1L], Main_CS2 = X[, 2L],
                     Int_CS1 = X[, 3L])
  V <- c(Main_CS1 = 0.3, Main_CS2 = 0.5, Int_CS1 = 0.4)
  phi <- 2.1
  y <- 0.2 + as.matrix(pred) %*% c(0.1, -0.1, 0.5, -0.3, 0.2) +
    rnorm(n, sd = sqrt(phi))

  exact <- SuSiE4I:::gaussian_ridge_refit(
    y, pred, V, phi, n_threads = 2L, block_size = 100L
  )
  dat <- cbind(data.frame(y = as.numeric(y)), pred)
  bam_fit <- SuSiE4I:::mgcv_fit_fixed_ridge(
    "y", names(pred), dat, gaussian(), V,
    dispersion = phi, mgcv_model = "bam"
  )

  expect_equal(exact$coefficients, coef(bam_fit), tolerance = 1e-9)
  expect_equal(exact$edf, sum(bam_fit$edf), tolerance = 1e-9)
  expect_equal(exact$dispersion, bam_fit$sig2, tolerance = 1e-9)
  expect_null(exact$hat_matrix)
})

test_that("Gaussian ridge fixed point agrees with the final bam scale", {
  set.seed(20260721)
  n <- 600L
  pred <- data.frame(
    Z1 = rnorm(n), Main_CS1 = rnorm(n), Int_CS1 = rnorm(n)
  )
  V <- c(Main_CS1 = 0.08, Int_CS1 = 0.12)
  y <- 0.3 + as.matrix(pred) %*% c(0.2, 0.6, -0.4) + rnorm(n, sd = 1.5)
  exact <- SuSiE4I:::gaussian_ridge_refit(
    y, pred, V, dispersion = 1, fixed_point = TRUE
  )
  dat <- cbind(data.frame(y = as.numeric(y)), pred)
  bam_fit <- SuSiE4I:::mgcv_fit_fixed_ridge(
    "y", names(pred), dat, gaussian(), V,
    dispersion = exact$dispersion, mgcv_model = "bam"
  )

  expect_equal(bam_fit$sig2, exact$dispersion, tolerance = 1e-9)
  expect_lt(abs(exact$dispersion - exact$penalty_dispersion), 1e-8)
})

test_that("binomial phi one reduces the ridge scale to 1/V", {
  set.seed(20260718)
  n <- 140L
  z <- rnorm(n)
  x <- rnorm(n)
  dat <- data.frame(
    y = rbinom(n, 1, plogis(0.3 * z + 0.8 * x)),
    Z1 = z, Main_CS1 = x
  )
  fit <- SuSiE4I:::mgcv_fit_fixed_ridge(
    "y", c("Z1", "Main_CS1"), dat, binomial(),
    c(Main_CS1 = 0.4), dispersion = 1, mgcv_model = "gam"
  )
  pen <- attr(fit, "refit_penalty")

  expect_equal(pen$working_precision, c(Main_CS1 = 2.5))
  expect_false("Z1" %in% names(pen$V))
  expect_true(all(c("Z1", "Main_CS1") %in% names(coef(fit))))
})

test_that("main and interaction terms map to their component V", {
  fitX <- list(V = c(0.2, 0.4))
  fitW <- list(V = c(0.3, 0.6, 0.9))
  out <- SuSiE4I:::refit_penalty_variance(
    fitX, fitW, c("Main_CS2", "Int_CS3")
  )
  expect_equal(out, c(Main_CS2 = 0.4, Int_CS3 = 0.9))

  shared <- SuSiE4I:::refit_penalty_variance(
    list(V = c(0.5, 0.5)), list(V = c(0.8, 0.8)),
    c("Main_noncs_res", "Int_noncs_res")
  )
  expect_equal(shared, c(Main_noncs_res = 0.5, Int_noncs_res = 0.8))
  component_noncs <- SuSiE4I:::refit_penalty_variance(
    fitX, fitW, c("Main_noncs_res", "Int_noncs_res")
  )
  expect_equal(
    component_noncs,
    c(Main_noncs_res = 0.4, Int_noncs_res = 0.9)
  )
})

test_that("environment CS terms map to fitZ while raw Z stays unpenalized", {
  fitX <- list(V = c(0.2, 0.4))
  fitW <- list(V = c(0.3, 0.6))
  fitZ <- list(V = c(0.7, 0.9))
  terms <- SuSiE4I:::refit_penalty_terms(c(
    "Z1", "Env_CS2", "Z_noncs_res", "Main_CS1", "Int_CS2"
  ))
  out <- SuSiE4I:::refit_penalty_variance(
    fitX, fitW, terms, fitZ = fitZ
  )

  expect_equal(
    out, c(Env_CS2 = 0.9, Z_noncs_res = 0.9, Main_CS1 = 0.2, Int_CS2 = 0.6)
  )
  expect_false("Z1" %in% names(out))
})

test_that("ordinal logit refit uses mgcv ocat with fixed V", {
  set.seed(20260720)
  n <- 180L
  z <- rnorm(n)
  x <- rnorm(n)
  latent <- 0.3 * z + 0.8 * x + stats::rlogis(n)
  y <- ordered(cut(latent, c(-Inf, -0.5, 0.5, Inf), labels = FALSE))
  pred <- data.frame(Z1 = z, Main_CS1 = x)
  dat <- cbind(data.frame(y = as.integer(y)), pred)
  fit <- SuSiE4I:::mgcv_fit_fixed_ridge(
    "y", colnames(pred), dat, mgcv::ocat(R = nlevels(y)),
    penalty_V = c(Main_CS1 = 0.5), dispersion = 1,
    mgcv_model = "gam"
  )
  pen <- attr(fit, "refit_penalty")

  expect_s3_class(fit, "gam")
  expect_equal(pen$V, c(Main_CS1 = 0.5))
  expect_equal(pen$working_precision, c(Main_CS1 = 2))
  expect_false("Z1" %in% names(pen$V))
})

test_that("Cox ridge uses theta 1/V and leaves nuisance terms unpenalized", {
  set.seed(20260719)
  n <- 160L
  z <- rnorm(n)
  x <- rnorm(n)
  te <- rexp(n, exp(0.25 * z + 0.7 * x))
  tc <- rexp(n, 0.5)
  dat <- data.frame(Z1 = z, Main_CS1 = x)
  fit <- SuSiE4I:::cox_fit_fixed_ridge(
    pmin(te, tc), te <= tc, dat, c(Main_CS1 = 0.5)
  )
  pen <- attr(fit, "refit_penalty")

  expect_s3_class(fit, "coxph.penal")
  expect_equal(pen$theta, c(Main_CS1 = 2))
  expect_false(pen$scale)
  expect_true(all(c("Z1", "Main_CS1") %in% names(coef(fit))))
  expect_false("Z1" %in% names(pen$V))
})
