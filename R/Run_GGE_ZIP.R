Run_GGE_ZIP <- function(X, Z, y,
                        family = mgcv::ziP(),
                        mgcv_model = NULL,
                        Lmain, Lint,
                        max.iter, min.iter, max.eps,
                        susie_para_main, susie_para_int,
                        noint_env = NULL,
                        verbose = TRUE, n_threads = 1,
                        L.init = 1,
                        x_noncs_var = 0.1,
                        w_noncs_var = 0.1,
                        noncs_max_abs_cor = 0.9,
                        include_x_squared = FALSE,
                        suff_block_size = 10000L,
                        returnModel = FALSE) {
  run_start <- proc.time()[["elapsed"]]
  n <- NROW(y)
  p <- ncol(X)
  y <- as.numeric(y)
  if (length(y) != nrow(X)) stop("Length(y) must equal nrow(X).")
  if (any(!is.finite(y)) || any(y < 0) ||
      any(abs(y - round(y)) > sqrt(.Machine$double.eps))) {
    stop("ZIP y must be finite non-negative integer counts.")
  }
  family <- zip_prepare_family(family)
  suff_block_size <- validate_suff_block_size(suff_block_size)
  ridge <- 1e-8
  eta_clip_range <- c(-50, 50)
  min_etaX_var <- 1e-7
  min_etaW_var <- 1e-7

  if (is.null(dim(Z))) Z <- matrix(Z, ncol = 1)
  if (nrow(Z) != n) stop("nrow(Z) must equal nrow(X).")
  if (is.null(colnames(Z))) colnames(Z) <- paste0("Z", seq_len(ncol(Z)))
  nZ <- ncol(Z)
  Zmodel <- Z

  response_info <- mgcv_prepare_response(y, family)
  warm <- mgcv_greedy_warm_start(
    X = X, response_info = response_info, Z = Zmodel,
    family = family, L.init = L.init,
    mgcv_model = mgcv_model
  )
  fit_final <- warm$fit

  g <- c()
  err <- Inf
  beta <- rep(0, p)
  alpha <- rep(0, nZ)
  theta_raw <- zip_theta_raw(fit_final$family)
  XCS <- NULL
  XCS_refit <- NULL
  WCS <- NULL
  WCS_refit <- NULL
  W <- NULL
  fitX <- NULL
  fitW <- NULL
  work <- NULL

  fitX_no_cs_streak <- 0L
  for (iter in seq_len(max.iter)) {
    beta_prev <- beta
    alpha_prev <- alpha
    theta_prev <- theta_raw

    work <- zip_extract_working(fit_final, y, eta_clip_range = eta_clip_range)


    ZI_main <- zip_nuisance_design(n, Zmodel, WCS_refit)
    ssX <- weighted_residual_suffstats(
      X = X, z = work$z, ZI = ZI_main, weights = work$weights,
      n_threads = n_threads, ridge = ridge, block_size = suff_block_size
    )
    fitX <- .fit_susie_stage(
      structural = list(
        XtX = ssX$XtX, Xty = ssX$Xty, yty = ssX$yty,
        n = max(0.95 * n, work$n_eff), L = Lmain
      ),
      susie_para = susie_para_main, stage = "main",
      iter = iter, min.iter = min.iter
    )
    beta <- coef.susie(fitX)[-1]

    CSdt <- summary(fitX)$vars
    x_component <- build_component_design_from_fit(
      X, fitX, "Main_CS"
    )
    cs_indices <- x_component$cs_indices
    fitX_no_cs_streak <- if (length(cs_indices)) 0L else fitX_no_cs_streak + 1L
    main_no_cs <- is.null(x_component$design)
    if (is.null(x_component$design)) {
        noncs_main <- build_full_noncs_refit_term(X, fitX)
        XCS <- NULL
        if (is.null(noncs_main)) {
          XCS_refit <- NULL
        } else {
          XCS_refit <- matrix(noncs_main, ncol = 1)
          colnames(XCS_refit) <- "Main_noncs_res"
        }
        W <- NULL
        WCS <- NULL
        WCS_refit <- NULL
        fitW <- NULL
    } else {
      XCS <- x_component$design
      XCS_refit <- XCS
      {
        main_noncs <- build_noncs_refit_term(
          X = X, fit = fitX, CSdt = CSdt, cs_indices = cs_indices,
          XCS = XCS, noncs_var = x_noncs_var,
          min_eta_var = min_etaX_var,
          noncs_max_abs_cor = noncs_max_abs_cor, cor_design = Zmodel
        )
        XCS_refit <- append_noncs_refit_term(
          XCS_refit, main_noncs, "Main_noncs_res",
          corr_threshold = noncs_max_abs_cor
        )
      }
    }

      XCS_W <- XCS_refit
      W <- get_pairwise_interactions(
        XCS_W, Z = Zmodel, noint_env = noint_env,
        include_x_squared = if (main_no_cs) FALSE else include_x_squared
      )
      if (!interaction_design_available(W, iter, min.iter, allow_empty = main_no_cs)) {
        W <- NULL
        fitW <- NULL
        WCS <- NULL
        WCS_refit <- NULL
      } else {
      ZI_int <- zip_nuisance_design(n, Zmodel, XCS_refit)
      ssW <- weighted_residual_suffstats(
        X = W, z = work$z, ZI = ZI_int, weights = work$weights,
        n_threads = n_threads, ridge = ridge, block_size = suff_block_size
      )
      fitW <- .fit_susie_stage(
        structural = list(
          XtX = ssW$XtX, Xty = ssW$Xty, yty = ssW$yty,
          n = max(0.95 * n, work$n_eff), L = Lint
        ),
        susie_para = susie_para_int, stage = "int",
        iter = iter, min.iter = min.iter
      )

      CSdt_w <- summary(fitW)$vars
      cs_indices_w <- sort(unique(CSdt_w$cs[CSdt_w$cs > 0]))
      w_component <- build_component_design_from_fit(
        W, fitW, "Int_CS"
      )
      WCS <- w_component$design
      WCS_refit <- WCS

      {
      w_noncs <- build_w_noncs_refit_term(
        W = W, fitW = fitW, WCS = WCS,
        etaX = matrixVectorMultiply(X, beta),
        XCS = XCS, Z = Zmodel,
        w_noncs_var = w_noncs_var, min_etaW_var = min_etaW_var,
        noncs_max_abs_cor = noncs_max_abs_cor
      )
      if (!is.null(w_noncs)) {
        if (is.null(WCS_refit)) {
          WCS_refit <- matrix(w_noncs, ncol = 1)
          colnames(WCS_refit) <- "Int_noncs_res"
        } else {
          WCS_refit <- append_noncs_refit_term(
            WCS_refit, w_noncs, "Int_noncs_res",
            corr_threshold = noncs_max_abs_cor
          )
        }
      }
      }
      }

    pred <- mgcv_predictor_data(
      Z = Zmodel, Xextra = cbind(XCS_refit, WCS_refit), n = n
    )
    Data <- cbind(response_info$data, pred)
    refit_dispersion <- mgcv_refit_dispersion(fit_final)
    penalty_names <- refit_penalty_terms(colnames(pred))
    penalty_V <- refit_penalty_variance(
      fitX, fitW, penalty_names
    )
    fit_final <- mgcv_fit_fixed_ridge(
      response_info$response, colnames(pred), Data, family, penalty_V,
      dispersion = refit_dispersion, mgcv_model = mgcv_model
    )

    coefs <- clean_coef(coef(fit_final))
    if (length(coefs) >= 1L + nZ) {
      alpha <- coefs[2:(1 + nZ)]
    }
    theta_raw <- zip_theta_raw(fit_final$family)

    errX <- sqrt(mean((beta - beta_prev)^2))
    errZ <- sqrt(mean((alpha - alpha_prev)^2))
    errT <- if (!is.null(theta_prev) && !is.null(theta_raw)) {
      sqrt(mean((theta_raw - theta_prev)^2))
    } else {
      0
    }
    err <- max(errX, errZ, errT)
    g[iter] <- err
    if (verbose) {
      cat(sprintf(
        "Iteration %d: err = %.3e, n_eff = %.1f, med_W = %.3g\n",
        iter, err, work$n_eff, work$info$med_weight
      ))
    }
    if (fitX_no_cs_streak >= 3L) {
      if (verbose) cat("No main credible set detected in 3 consecutive iterations; stopping.\n")
      break
    }
    if (err < max.eps && iter > min.iter) {
      if (verbose) cat("Converged!\n")
      break
    }
  }

  pred <- mgcv_predictor_data(
    Z = Zmodel, Xextra = cbind(XCS_refit, WCS_refit), n = n
  )
  Dat <- cbind(response_info$data, pred)
  refit_dispersion <- mgcv_refit_dispersion(fit_final)
  penalty_names <- refit_penalty_terms(colnames(pred))
  penalty_V <- refit_penalty_variance(
    fitX, fitW, penalty_names
  )
  fit_final <- mgcv_fit_fixed_ridge(
    response_info$response, colnames(pred), Dat, family, penalty_V,
    dispersion = refit_dispersion, mgcv_model = mgcv_model
  )
  fit_final$n_eff <- work$n_eff
  G <- summary(fit_final)$p.table
  MainIndex <- Identifying_MainEffect(fitX, colnames(X))
  MainIndex <- safe_add_p(MainIndex, G)
  IntIndex <- Identifying_IntEffect(fitW, colnames(W))
  IntIndex <- filter_noncs_interactions(IntIndex)
  IntIndex <- safe_add_p(IntIndex, G)
  if (verbose) {
    plot(g, type = "o", col = "black", pch = 16,
         xlab = "Iteration", ylab = "Max Parameter Change",
         main = "Convergence Trace (ZIP)")
    for (i in seq_along(g)) {
      text(x = i, y = g[i], labels = formatC(g[i], format = "e", digits = 1),
           pos = 3, cex = 0.7, col = "red")
    }
  }

  AA <- list(
    diagnostics = make_diagnostics(iter, g, run_start),
    fitX = fitX,
    fitW = fitW,
    fitJoint = fit_final,
    main_discoveries = MainIndex,
    interaction_discoveries = IntIndex
  )
  if (returnModel) AA$FinalModel <- Dat
  AA$report <- extract_direction_table(AA, G)
  return(AA)
}
