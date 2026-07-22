library(CppMatrix)
library(susieR)
library(survival)
library(devtools)
library(MASS)
document()

ARcov <- function(p, rho) {
  s <- rho^(0:(p - 1))
  return(toeplitz(s))
}

CScov <- function(p, rho) {
  diag(p) * (1 - rho) + matrix(rho, p, p)
}

source("example/otherfunction.R")

## -----------------------------------------------------------------------------
## Simulate right-censored survival data from a Cox proportional-hazards model.
##   hazard(t | X, Z) = lambda0 * exp(eta),  eta = Z alpha + X beta + interactions
## Event times: exponential baseline, T = -log(U) / (lambda0 * exp(eta)).
## Censoring:   independent exponential, tuned to a target censoring fraction.
## -----------------------------------------------------------------------------
sim_cox <- function(n, p, lambda0 = 0.1, cens_rate = 0.05) {
  R <- kronecker(CScov(5, 0.3), ARcov(2, 0.5))
  X <- mvrnorm(n = n, mu = runif(p, 0, 1), R)
  colnames(X) <- paste0("X", seq_len(p))
  Z <- mvrnorm(n = n, mu = runif(15, 0, 1), diag(15))
  colnames(Z) <- paste0("UKBB", 1:15)

  alpha0 <- rep(0, 15)
  alpha0[1:5] <- rnorm(5, 0, 1 / sqrt(5))
  beta0 <- rnorm(p, 0, sqrt(5e-5))
  ind <- sort(sample(p, 3))
  beta0[ind] <- 0.25

  ## two interactions: X_ind1 * X_ind2 (GxG) and X_ind1 * UKBB1 (GxE)
  eta <- matrixVectorMultiply(Z, alpha0) +
    matrixVectorMultiply(X, beta0) +
    (X[, ind[1]] * X[, ind[2]] + X[, ind[1]] * Z[, 1]) / 2

  U <- runif(n)
  Ttrue <- -log(U) / (lambda0 * exp(eta))
  Cens  <- rexp(n, rate = cens_rate)
  time   <- pmin(Ttrue, Cens)
  status <- as.integer(Ttrue <= Cens)

  true_int_gxg <- paste0("X", ind[1], "*X", ind[2])     # GxG truth
  true_int_gei <- paste0("UKBB1", "*X", ind[1])         # GxE truth

  list(X = X, Z = Z, time = time, status = status,
       ind = ind, true_int_variable = c(true_int_gxg, true_int_gei))
}

set.seed(1)
n <- 1000
p <- 10

## =============================================================================
## Single run: inspect that fine-mapping recovers the truth
## =============================================================================
dat <- sim_cox(n = n, p = p)
cat(sprintf("Events: %d / %d (%.1f%% censored)\n",
            sum(dat$status), n, 100 * mean(dat$status == 0)))
cat("True main effects (columns of X):", dat$ind, "\n")
cat("True interactions:", paste(dat$true_int_variable, collapse = ", "), "\n\n")

## --- 1. Fixed Z, GxG + GxE; UKBB11..15 excluded from GEI search ---
fit1 <- SuSiE4I(
  X = dat$X, Z = dat$Z, y = dat$time, status = dat$status,
  family = "cox",
  select_env = FALSE, noint_env = 11:15,
  L_main = 5, L_int = 5, n_threads = 4,
  max_iter = 15, susie_iter = 500,
  coverage_main = 0.9, coverage_int = 0.9,
   verbose = TRUE
)

cat("\n--- Run_GGE_Cox: main_index ---\n");        print(fit1$main_discoveries)
cat("\n--- Run_GGE_Cox: interaction_index ---\n"); print(fit1$interaction_discoveries)
cat("\n--- Run_GGE_Cox: direction report ---\n");  print(fit1$report)

## --- 2. Z also fine-mapped ---
fit2 <- SuSiE4I(
  X = dat$X, Z = dat$Z, y = dat$time, status = dat$status,
  family = "cox",
  select_env = TRUE,
  L_main = 5, L_int = 5, L_env = 10, n_threads = 4,
  max_iter = 15, susie_iter = 500,
  coverage_main = 0.9, coverage_int = 0.9,
   verbose = FALSE
)
cat("\n--- Run_GGE_Select_Cox: main_index ---\n");        print(fit2$main_discoveries)
cat("\n--- Run_GGE_Select_Cox: env_index ---\n");         print(fit2$environment_discoveries)
cat("\n--- Run_GGE_Select_Cox: interaction_index ---\n"); print(fit2$interaction_discoveries)

## --- 3. No Z, GxG only ---
fit3 <- SuSiE4I(
  X = dat$X, Z = NULL, y = dat$time, status = dat$status,
  family = "cox",
  L_main = 5, L_int = 5, n_threads = 4,
  max_iter = 15, susie_iter = 500,
  coverage_main = 0.9, coverage_int = 0.9,
  estimate_residual_variance = FALSE, verbose = FALSE
)
cat("\n--- Run_GG_Cox: main_index ---\n");        print(fit3$main_discoveries)
cat("\n--- Run_GG_Cox: interaction_index ---\n"); print(fit3$interaction_discoveries)

## =============================================================================
## Repeated simulation: TP / TN rates for main and interaction recovery
##   (GG_Cox only sees GxG, so it is scored against the GxG interaction alone)
## =============================================================================
nrep <- 20
Main_TP <- Main_TN <- Int_TP <- Int_TN <- matrix(0, nrep, 3)

for (iter in 1:nrep) {
  d <- sim_cox(n = n, p = p)
  gxg_only <- d$true_int_variable[1]  # the single GxG truth for the GG model

  f1 <- SuSiE4I(X = d$X, Z = d$Z, y = d$time, status = d$status,
                      family = "cox",
                      select_env = FALSE, noint_env = 11:15,
                      L_main = 5, L_int = 5, n_threads = 4, max_iter = 15,
                      susie_iter = 500, coverage_main = 0.9, coverage_int = 0.9,
                       verbose = FALSE)
  f2 <- SuSiE4I(X = d$X, Z = d$Z, y = d$time, status = d$status,
                      family = "cox",
                      select_env = TRUE,
                      L_main = 5, L_int = 5, L_env = 10, n_threads = 4, max_iter = 15,
                      susie_iter = 500, coverage_main = 0.9, coverage_int = 0.9,
                       verbose = FALSE)
  f3 <- SuSiE4I(X = d$X, Z = NULL, y = d$time, status = d$status,
                      family = "cox",
                      L_main = 5, L_int = 5, n_threads = 4, max_iter = 15,
                      susie_iter = 500, coverage_main = 0.9, coverage_int = 0.9,
                       verbose = FALSE)

  e1 <- SuSiE4I_evaluate(true_main_index = d$ind, true_int_variable = d$true_int_variable, f1)
  e2 <- SuSiE4I_evaluate(true_main_index = d$ind, true_int_variable = d$true_int_variable, f2)
  e3 <- SuSiE4I_evaluate(true_main_index = d$ind, true_int_variable = gxg_only,            f3)

  Main_TP[iter, ] <- c(e1$tp_main, e2$tp_main, e3$tp_main)
  Main_TN[iter, ] <- c(e1$tn_main, e2$tn_main, e3$tn_main)
  Int_TP[iter, ]  <- c(e1$tp_int,  e2$tp_int,  e3$tp_int)
  Int_TN[iter, ]  <- c(e1$tn_int,  e2$tn_int,  e3$tn_int)

  if (iter %% 5 == 0) cat(sprintf("rep %d done\n", iter))
}

colnames(Main_TP) <- colnames(Main_TN) <- colnames(Int_TP) <- colnames(Int_TN) <-
  c("GGE_Cox", "Select_Cox", "GG_Cox")

cat("\n===== Final Summary (mean over reps) =====\n")
cat("Main TP (power):       \n"); print(colMeans(Main_TP))
cat("Main TN (1 - false pos):\n"); print(colMeans(Main_TN))
cat("Int  TP (power):       \n"); print(colMeans(Int_TP))
cat("Int  TN (1 - false pos):\n"); print(colMeans(Int_TN))
