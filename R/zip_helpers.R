zip_is_family <- function(family) {
  inherits(family, "family") &&
    grepl("zero inflated poisson", tolower(paste(family$family, collapse = " ")),
          fixed = TRUE)
}

zip_prepare_family <- function(family = NULL) {
  if (is.null(family)) return(mgcv::ziP())
  if (is.character(family)) {
    if (length(family) != 1L) stop("family must be a scalar string or a ZIP family object.")
    family_string <- tolower(family)
    if (family_string %in% c("zip", "ziP", "zero-inflated-poisson",
                             "zero_inflated_poisson", "zero inflated poisson")) {
      return(mgcv::ziP())
    }
  }
  if (!zip_is_family(family)) {
    stop("family must be mgcv::ziP() or family = 'zip'.")
  }
  family
}

zip_theta_raw <- function(family) {
  if (!is.function(family$getTheta)) return(NULL)
  theta_raw <- family$getTheta(FALSE)
  if (is.null(theta_raw) || length(theta_raw) != 2L ||
      any(!is.finite(theta_raw))) {
    return(NULL)
  }
  as.numeric(theta_raw)
}

zip_theta <- function(family) {
  if (!is.function(family$getTheta)) return(NULL)
  theta <- family$getTheta(TRUE)
  if (is.null(theta) || length(theta) != 2L || any(!is.finite(theta))) {
    return(NULL)
  }
  as.numeric(theta)
}

zip_family_b <- function(family) {
  theta_raw <- zip_theta_raw(family)
  theta <- zip_theta(family)
  if (is.null(theta_raw) || is.null(theta)) return(0)
  b <- theta[2] - exp(theta_raw[2])
  if (!is.finite(b) || b < 0) b <- 0
  as.numeric(b)
}

zip_nuisance_design <- function(n, ...) {
  out <- matrix(1, nrow = n, ncol = 1)
  colnames(out) <- "Intercept"
  parts <- list(...)
  for (part in parts) {
    if (is.null(part)) next
    part <- as.matrix(part)
    if (nrow(part) != n) stop("Nuisance design has incompatible row count.")
    if (ncol(part) == 0L) next
    out <- cbind(out, part)
  }
  out
}

weighted_residual_suffstats <- function(X, z, ZI, weights,
                                        n_threads = 1,
                                        ridge = 1e-8,
                                        block_size = 10000L) {
  weighted_projected_suffstats(
    X = X, y = z, ZI = ZI, weights = weights,
    n_threads = n_threads, ridge = ridge,
    block_size = block_size
  )
}

zip_extract_working <- function(fit, y, eta_clip_range = c(-50, 50)) {
  gamma <- as.numeric(fit$linear.predictors)
  gamma <- pmin(pmax(gamma, eta_clip_range[1]), eta_clip_range[2])
  theta_raw <- zip_theta_raw(fit$family)
  theta <- zip_theta(fit$family)
  if (is.null(theta_raw)) stop("Could not extract raw theta from mgcv::ziP() family.")

  wt <- fit$prior.weights
  if (is.null(wt)) wt <- rep(1, length(y))
  dd <- fit$family$Dd(as.numeric(y), gamma, theta_raw, wt = wt, level = 0)
  score <- -0.5 * as.numeric(dd$Dmu)
  W <- 0.5 * as.numeric(dd$EDmu2)

  bad <- !is.finite(gamma) | !is.finite(score) | !is.finite(W) | W <= 0
  if (mean(bad) > 0.9) stop("Too many invalid ZIP working observations.")
  z <- gamma
  ok <- !bad
  z[ok] <- gamma[ok] + score[ok] / W[ok]
  z[!is.finite(z)] <- 0
  W[bad] <- 0

  denom <- sum(W^2)
  if (!is.finite(denom) || denom <= 0) stop("All ZIP working weights are zero.")
  n_eff <- (sum(W)^2) / denom

  list(
    z = z,
    weights = W,
    gamma = gamma,
    score = score,
    n_eff = n_eff,
    theta = theta,
    theta_raw = theta_raw,
    info = list(
      n_eff = n_eff,
      min_weight = min(W[W > 0]),
      med_weight = stats::median(W[W > 0]),
      max_weight = max(W),
      pos_rate = mean(as.numeric(y) > 0),
      b = zip_family_b(fit$family)
    )
  )
}
