test_that("the public interface is ordered and minimal", {
  nm <- names(formals(SuSiE4I))
  expect_lt(match("L_main", nm), match("L_int", nm))
  expect_lt(match("L_int", nm), match("L_env", nm))
  expect_lt(match("susie_para_main", nm), match("susie_para_int", nm))
  expect_lt(match("susie_para_int", nm), match("susie_para_env", nm))
  expect_equal(formals(SuSiE4I)$max_iter, 10)
  expect_false(any(c(
    "init_cor_method", "prior_variance", "susie_iter", "ridge",
    "eta_clip_range", "min_etaX_var", "min_etaW_var", "continue_no_cs",
    "noncs_w", "..."
  ) %in% nm))
})

test_that("every family uses separate GG, GGE, and select-env runners", {
  runners <- c(
    "Run_GG", "Run_GGE", "Run_GGE_Select",
    "Run_GG_GLM", "Run_GGE_GLM", "Run_GGE_Select_GLM",
    "Run_GG_OCAT", "Run_GGE_OCAT", "Run_GGE_Select_OCAT",
    "Run_GG_CLM", "Run_GGE_CLM", "Run_GGE_Select_CLM",
    "Run_GG_Cox", "Run_GGE_Cox", "Run_GGE_Select_Cox",
    "Run_GG_ZIP", "Run_GGE_ZIP", "Run_GGE_Select_ZIP"
  )
  expect_true(all(vapply(runners, exists, logical(1), mode = "function")))

  gg <- runners[grepl("^Run_GG(_|$)", runners)]
  gge <- runners[grepl("^Run_GGE_", runners) &
                   !grepl("^Run_GGE_Select", runners)]
  select <- runners[grepl("^Run_GGE_Select", runners)]
  expect_true(all(vapply(gg, function(nm) {
    !"Z" %in% names(formals(get(nm)))
  }, logical(1))))
  expect_true(all(vapply(gge, function(nm) {
    "Z" %in% names(formals(get(nm)))
  }, logical(1))))
  expect_true(all(vapply(select, function(nm) {
    "Z" %in% names(formals(get(nm)))
  }, logical(1))))
  expect_false(any(vapply(runners, function(nm) {
    "select_env" %in% names(formals(get(nm)))
  }, logical(1))))
})

test_that("the public family dispatch reaches the 18 explicit runners", {
  runners <- c(
    "Run_GG", "Run_GGE", "Run_GGE_Select",
    "Run_GG_GLM", "Run_GGE_GLM", "Run_GGE_Select_GLM",
    "Run_GG_OCAT", "Run_GGE_OCAT", "Run_GGE_Select_OCAT",
    "Run_GG_CLM", "Run_GGE_CLM", "Run_GGE_Select_CLM",
    "Run_GG_Cox", "Run_GGE_Cox", "Run_GGE_Select_Cox",
    "Run_GG_ZIP", "Run_GGE_ZIP", "Run_GGE_Select_ZIP"
  )
  mocks <- stats::setNames(lapply(runners, function(nm) {
    force(nm)
    function(...) nm
  }), runners)
  do.call(testthat::local_mocked_bindings,
          c(mocks, list(.package = "SuSiE4I")))

  X <- matrix(seq_len(24), nrow = 8L)
  Z <- matrix(seq_len(16), nrow = 8L)
  run_case <- function(family, expected, y, status = NULL,
                       use_z = FALSE, select_env = FALSE) {
    args <- list(
      X = X, y = y, family = family,
      L_main = 1, L_int = 1, L_env = 1,
      max_iter = 1, min_iter = 0, verbose = FALSE
    )
    if (use_z) args$Z <- Z
    if (select_env) args$select_env <- TRUE
    if (!is.null(status)) args$status <- status
    expect_identical(do.call(SuSiE4I, args), expected)
  }

  structures <- list(
    list(FALSE, FALSE, "GG"),
    list(TRUE, FALSE, "GGE"),
    list(TRUE, TRUE, "GGE_Select")
  )
  for (s in structures) {
    run_case("gaussian", paste0("Run_", s[[3L]]), seq_len(8),
             use_z = s[[1L]], select_env = s[[2L]])
    run_case("binomial", paste0("Run_", s[[3L]], "_GLM"),
             rep(c(0, 1), 4), use_z = s[[1L]], select_env = s[[2L]])
    run_case("clm_logit", paste0("Run_", s[[3L]], "_OCAT"),
             ordered(rep(1:4, each = 2)),
             use_z = s[[1L]], select_env = s[[2L]])
    run_case("clm_probit", paste0("Run_", s[[3L]], "_CLM"),
             ordered(rep(1:4, each = 2)),
             use_z = s[[1L]], select_env = s[[2L]])
    run_case("cox", paste0("Run_", s[[3L]], "_Cox"), seq_len(8),
             status = rep(c(0, 1), 4),
             use_z = s[[1L]], select_env = s[[2L]])
    run_case("zip", paste0("Run_", s[[3L]], "_ZIP"),
             rep(c(0, 1, 0, 2), 2),
             use_z = s[[1L]], select_env = s[[2L]])
  }
})

test_that("Gaussian public dispatch reaches all three dedicated runners", {
  set.seed(20260719)
  n <- 120L
  X <- matrix(rnorm(n * 4L), n)
  Z <- matrix(rnorm(n * 2L), n)
  y <- 1.5 * X[, 1L] + 3 * Z[, 1L] + rnorm(n, sd = 0.5)

  fit_x <- SuSiE4I(
    X, y = y, family = gaussian(), include_x_squared = TRUE,
    L_main = 1, L_int = 1, max_iter = 3, min_iter = 1, verbose = FALSE
  )
  fit_xz <- SuSiE4I(
    X, Z, y = y, family = gaussian(), include_x_squared = TRUE,
    L_main = 1, L_int = 1, max_iter = 3, min_iter = 1, verbose = FALSE
  )
  fit_select <- SuSiE4I(
    X, Z, y = y, family = gaussian(), select_env = TRUE,
    include_x_squared = TRUE, L_main = 1, L_int = 1, L_env = 1,
    max_iter = 3, min_iter = 1, verbose = FALSE
  )

  expect_s3_class(fit_x$fitJoint, "gam")
  expect_s3_class(fit_xz$fitJoint, "gam")
  expect_s3_class(fit_select$fitJoint, "gam")
})

test_that("Gaussian public dispatch forwards crossprodX and L.init", {
  set.seed(20260720)
  n <- 120L
  X <- matrix(rnorm(n * 4L), n)
  Z <- matrix(rnorm(n * 2L), n)
  y <- 1.5 * X[, 1L] + 3 * Z[, 1L] + rnorm(n, sd = 0.5)

  expect_error(
    SuSiE4I(
      X, y = y, family = gaussian(), crossprodX = diag(3),
      include_x_squared = TRUE, L_main = 1, L_int = 1,
      max_iter = 3, min_iter = 1, verbose = FALSE
    ),
    "crossprodX must be"
  )
  expect_error(
    SuSiE4I(
      X, Z, y = y, family = gaussian(), select_env = TRUE,
      crossprodX = crossprod(X), L.init = 0, include_x_squared = TRUE,
      L_main = 1, L_int = 1, L_env = 1,
      max_iter = 3, min_iter = 1, verbose = FALSE
    ),
    "L.init must be"
  )
})

test_that("Gaussian runners remain independent from GLM runner bodies", {
  gaussian_files <- c("R/Run_GG.R", "R/Run_GGE.R", "R/Run_GGE_Select.R")
  glm_files <- c("R/Run_GG_GLM.R", "R/Run_GGE_GLM.R", "R/Run_GGE_Select_GLM.R")
  gaussian_src <- vapply(
    gaussian_files, function(path) {
      paste(readLines(file.path("..", "..", path)), collapse = "\n")
    },
    character(1)
  )
  glm_src <- vapply(
    glm_files, function(path) {
      paste(readLines(file.path("..", "..", path)), collapse = "\n")
    },
    character(1)
  )

  expect_false(any(grepl("[.]run_(gg|gge|gge_select)_current", gaussian_src)))
  expect_false(any(grepl("is_gaussian|gaussian_ridge_refit", glm_src)))
})

test_that("diagnostics use a one-row data frame", {
  d <- SuSiE4I:::make_diagnostics(3, c(0.2, 0.1), proc.time()[["elapsed"]])

  expect_s3_class(d, "data.frame")
  expect_equal(nrow(d), 1)
  expect_identical(names(d), c("iterations", "eps", "runtime_seconds"))
  expect_equal(d$iterations, 3L)
  expect_equal(d$eps, 0.1)
})

test_that("stage defaults differ only where intended", {
  main <- SuSiE4I:::.susie_default_para("main")
  int <- SuSiE4I:::.susie_default_para("int")
  env <- SuSiE4I:::.susie_default_para("env")

  expect_equal(main$coverage, 0.95)
  expect_equal(int$coverage, 0.95)
  expect_equal(env$coverage, 0.95)
  expect_equal(main$min_abs_corr, 0.5)
  expect_false(main$standardize)
  expect_false(int$standardize)
  expect_false(env$standardize)
  expect_equal(main$scaled_prior_variance, 2)
  expect_true(main$estimate_prior_variance)
  expect_equal(main$estimate_prior_method, "optim")
  expect_equal(main$residual_variance_upperbound, 1.01)

  gaussian_main <- SuSiE4I:::.susie_default_para(
    "main", gaussian = TRUE, residual_variance = 0.7
  )
  expect_false(gaussian_main$estimate_residual_variance)
  expect_equal(gaussian_main$residual_variance, 0.7)
  expect_false("residual_variance_lowerbound" %in% names(gaussian_main))
  expect_false("residual_variance_upperbound" %in% names(gaussian_main))
})

test_that("warm-up fixes V while retaining other SuSiE controls", {
  expect_warning(
    user <- SuSiE4I:::.resolve_susie_para(
      list(
        scaled_prior_variance = 1.5,
        estimate_prior_variance = TRUE,
        estimate_residual_variance = FALSE,
        max_iter = 7,
        coverage = 0.8
      ),
      "susie_para_main"
    ),
    "interpreted as an absolute"
  )
  structural <- list(
    XtX = diag(2), Xty = c(1, 0), yty = 1, n = 10, L = 2
  )

  warm <- SuSiE4I:::.susie_iteration_args(
    user, structural, "main", iter = 2, min.iter = 2
  )
  post <- SuSiE4I:::.susie_iteration_args(
    user, structural, "main", iter = 3, min.iter = 2
  )

  expect_equal(warm$scaled_prior_variance, 2)
  expect_false(warm$standardize)
  expect_false(warm$estimate_prior_variance)
  expect_false(warm$estimate_residual_variance)
  expect_equal(warm$max_iter, 7)
  expect_equal(warm$coverage, 0.8)
  expect_equal(post$scaled_prior_variance, 1.5 / (1 / 9))
  expect_false(post$standardize)
  expect_true(post$estimate_prior_variance)
  expect_false(post$estimate_residual_variance)
  expect_equal(post$max_iter, 7)
  expect_equal(post$coverage, 0.8)
})

test_that("absolute prior variance is applied only after warm-up", {
  structural <- list(
    XtX = diag(2), Xty = c(2, 0.5), yty = 120, n = 21, L = 2
  )
  fixed <- SuSiE4I:::.resolve_susie_para(list(
    prior_variance = 0.7, estimate_prior_variance = FALSE, max_iter = 50
  ), "susie_para_main")
  warm <- SuSiE4I:::.susie_iteration_args(
    fixed, structural, "main", iter = 1, min.iter = 2
  )
  post <- SuSiE4I:::.susie_iteration_args(
    fixed, structural, "main", iter = 3, min.iter = 2
  )

  expect_false("prior_variance" %in% names(warm))
  expect_equal(warm$scaled_prior_variance, 2)
  expect_false(warm$estimate_prior_variance)
  expect_equal(warm$max_iter, 50)
  expect_equal(post$scaled_prior_variance, 0.7 / (120 / 20))
  expect_false(post$estimate_prior_variance)
  fit <- do.call(susieR::susie_ss, post)
  expect_equal(as.numeric(fit$V), rep(0.7, 2), tolerance = 1e-8)

  adaptive <- SuSiE4I:::.resolve_susie_para(list(
    prior_variance = 0.7, estimate_prior_variance = TRUE
  ), "susie_para_main")
  initial <- SuSiE4I:::.susie_iteration_args(
    adaptive, structural, "main", iter = 1, min.iter = 2
  )
  adaptive_post <- SuSiE4I:::.susie_iteration_args(
    adaptive, structural, "main", iter = 3, min.iter = 2
  )
  expect_false(initial$estimate_prior_variance)
  expect_equal(initial$scaled_prior_variance, 2)
  expect_true(adaptive_post$estimate_prior_variance)
  expect_equal(adaptive_post$scaled_prior_variance, 0.7 / 6)
})

test_that("prior scale validation is explicit and scaled warning occurs once", {
  expect_error(
    SuSiE4I:::.resolve_susie_para(list(
      prior_variance = 1, scaled_prior_variance = 1
    ), "susie_para_int"),
    "cannot contain both"
  )
  expect_warning(
    scaled <- SuSiE4I:::.resolve_susie_para(
      list(scaled_prior_variance = 1, estimate_prior_variance = FALSE),
      "susie_para_main"
    ),
    "interpreted as an absolute"
  )
  expect_equal(scaled$prior_variance, 1)
  expect_false("scaled_prior_variance" %in% names(scaled))
  structural <- list(
    XtX = diag(2), Xty = c(1, 0), yty = 90, n = 10, L = 2
  )
  args <- SuSiE4I:::.susie_iteration_args(
    scaled, structural, "int", iter = 3, min.iter = 2
  )
  expect_equal(args$scaled_prior_variance, 0.1)
  fit <- do.call(susieR::susie_ss, args)
  expect_equal(as.numeric(fit$V), rep(1, 2), tolerance = 1e-8)
})

test_that("prior variance controls are independent across stages", {
  main <- SuSiE4I:::.resolve_susie_para(
    list(prior_variance = 0.7, estimate_prior_variance = FALSE),
    "susie_para_main"
  )
  int <- SuSiE4I:::.resolve_susie_para(NULL, "susie_para_int")
  env <- SuSiE4I:::.resolve_susie_para(
    list(prior_variance = 0.3, estimate_prior_variance = TRUE),
    "susie_para_env"
  )

  expect_equal(main$prior_variance, 0.7)
  expect_false("prior_variance" %in% names(int))
  expect_equal(env$prior_variance, 0.3)
  expect_true(env$estimate_prior_variance)
})

test_that("mixed mode may continue without orphan interaction candidates", {
  expect_null(get_pairwise_interactions(NULL))
  expect_false(interaction_design_available(
    NULL, iter = 3, min_iter = 2, allow_empty = TRUE
  ))
  expect_error(
    interaction_design_available(NULL, iter = 3, min_iter = 2),
    "No interaction candidates"
  )
})

test_that("NULL CS overrides resolve to current susieR native defaults", {
  structural <- list(
    XtX = diag(2), Xty = c(1, 0), yty = 1, n = 10, L = 2
  )
  args <- SuSiE4I:::.susie_iteration_args(
    list(coverage = NULL, min_abs_corr = NULL),
    structural, "main", iter = 3, min.iter = 2
  )

  expect_equal(args$coverage, 0.95)
  expect_equal(args$min_abs_corr, 0.5)
})

test_that("fit stores the effective CS configuration", {
  structural <- list(
    XtX = diag(2), Xty = c(1, 0), yty = 1, n = 10, L = 2
  )
  native <- SuSiE4I:::.fit_susie_stage(
    structural = structural,
    susie_para = list(coverage = NULL, min_abs_corr = NULL, max_iter = 20),
    stage = "main", iter = 3, min.iter = 2
  )
  expect_equal(native$sets$requested_coverage, 0.95)
  expect_identical(
    native$cs_config,
    list(coverage = 0.95, min_abs_corr = 0.5)
  )

  fit <- SuSiE4I:::.fit_susie_stage(
    structural = structural,
    susie_para = list(coverage = 0.8, min_abs_corr = 0.6, max_iter = 20),
    stage = "main", iter = 3, min.iter = 2
  )

  expect_equal(fit$sets$requested_coverage, 0.8)
  expect_identical(fit$cs_config, list(coverage = 0.8, min_abs_corr = 0.6))
})

test_that("Cox and OCAT diagnostics retain the pre-ridge cross-product", {
  set.seed(42)
  X <- matrix(rnorm(180), 60, 3)
  Z <- matrix(rnorm(120), 60, 2)
  ridge <- 1e-6

  cox_ss <- cox_suffstat_block(
    X, rep(0, 60), Z, rexp(60), rbinom(60, 1, 0.7),
    nuisance_precision = numeric(0),
    n_threads = 1, ridge = ridge
  )
  expect_equal(diag(cox_ss$XtX) - diag(cox_ss$XtX_pre_ridge),
               rep(ridge, 3))

  ocat_ss <- ocat_suffstat_block(
    X, rep(1:3, each = 20), rep(0, 60), Z, c(-1, 1),
    n_threads = 1, ridge = ridge
  )
  expect_equal(diag(ocat_ss$XtX) - diag(ocat_ss$XtX_pre_ridge),
               rep(ridge, 3))
})

test_that("invalid and structural SuSiE names stop", {
  expect_error(
    SuSiE4I:::.resolve_susie_para(list(L = 2), "susie_para_main"),
    "cannot set structural"
  )
  expect_error(
    SuSiE4I:::.resolve_susie_para(
      list(not_a_susie_argument = 1), "susie_para_int"
    ),
    "Unknown susieR::susie_ss parameter"
  )
})

test_that("Cox environment initialization uses 1se, min, then correlation", {
  Z <- cbind(
    z1 = c(-2, -1, 0, 1, 2),
    z2 = c(1, -1, 1, -1, 1),
    z3 = c(2, 0, -2, 0, 2)
  )
  y <- c(1, 2, 3, 4, 5)
  expect_identical(
    select_initial_cox_env(Z, y, c(0, 2, 0), c(3, 0, 0)), 2L
  )
  expect_identical(
    select_initial_cox_env(Z, y, numeric(3), c(3, 0, 0)), 1L
  )
  expect_identical(
    select_initial_cox_env(Z, y, numeric(3), numeric(3)), 1L
  )
})

test_that("continuous environment initialization returns a credible set", {
  set.seed(10)
  Z <- scale(matrix(rnorm(300 * 5), 300, 5))
  y <- Z[, 1] + rnorm(300, sd = 0.4)
  fit <- initial_continuous_env(Z, y, Lenv = 2, n_threads = 1)

  expect_s3_class(fit$fitZ, "susie")
  expect_true(is.matrix(fit$ZCS))
  expect_equal(nrow(fit$ZCS), nrow(Z))
})
