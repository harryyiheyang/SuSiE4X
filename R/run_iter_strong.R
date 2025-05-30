run_iter_strong <- function(X, Z, y, crossprodX=NULL, Lmain, Linteraction, max.iter, min.iter, max.eps, susie.iter,verbose=T,n_threads=1,...) {

n <- length(y)
Xmean <- colMeans(X)

if(is.null(crossprodX)){
XtX <- blockwise_crossprod(X,X,n_threads) - n * tcrossprod(Xmean)
}else{
XtX = crossprodX - n * tcrossprod(Xmean)
}
ZtZ <- blockwise_crossprod(Z,Z,n_threads)
etaX <- etaZ <- etaW <- 0
fitX=NULL
g=c()
beta=as.vector(X[1,])*0
alpha=as.vector(Z[1,])*0
gamma=0

for (iter in 1:max.iter) {
gamma_prev <- gamma
beta_prev <- beta
alpha_prev <- alpha

## --- update etaZ ---
rZ <- y - etaX - etaW
Zty <- crossprod(Z, rZ)
alpha <- as.vector(solve(ZtZ, Zty))
etaZ <- matrixVectorMultiply(Z, alpha)

## --- update etaX ---
rX <- y - etaZ - etaW
Xty <- as.vector(crossprod(X, rX))
yty4X <- var(rX) * n
fitX <- susie_suff_stat(XtX = XtX, Xty = Xty, yty = yty4X, n = n, L = Lmain, max_iter = susie.iter, estimate_prior_method = "EM")
beta <- coef.susie(fitX)[-1]
etaX <- matrixVectorMultiply(X, beta)
etaX <- etaX - mean(etaX)

## --- update etaW ---
rW <- y- etaZ - etaX
XCS <- matrixMultiply(X, t(as.matrix(fitX$alpha)))
colnames(XCS) <- paste0("Main_CS", seq_len(ncol(XCS)))
cs <- get_active_indices(fitX)
if(length(cs)==0){
  stop("No credible set of marginal effects detected. SuSiE4X stops.")
}
XCS <- as.matrix(XCS[, cs, drop = FALSE])
G = cbind(Z,XCS)
W <- get_pairwise_interactions(XCS, Z[, -1])
GtG = blockwise_crossprod(G,G,n_threads)
GtW = blockwise_crossprod(G,W,n_threads)
ProjPart = matrixMultiply(G,(solve(GtG)%*%(GtW)))
W= W - ProjPart
WtW <- blockwise_crossprod(W,W,n_threads)
Wty <- as.vector(crossprod(W, rW))
yty4W <- var(rW) * n
fitW <- susie_suff_stat(XtX = WtW, Xty = Wty, yty = yty4W, n = n, L = Linteraction, max_iter = susie.iter, estimate_prior_method = "EM")
gamma <- coef.susie(fitW)[-1]
etaW <- matrixVectorMultiply(W, gamma)

## --- check convergence ---
errX <- sqrt(mean((beta - beta_prev)^2))
errZ <- sqrt(mean((alpha - alpha_prev)^2))
err = c(errX, errZ)
g[iter]=max(err)
if (verbose == TRUE){
cat(paste0("iteration ",iter,"\n"))
}
if (max(err) < max.eps&iter>min.iter) {
  break
}
}

## --- post-hoc: linear model of Z with offset ---
fitZ <- lm(y ~ Z-1, offset = etaX + etaW)

MainIndex=Identifying_MainEffect(fitX,colnames(X))
IntIndex=Identifying_IntEffect(fitW,colnames(W))


if (verbose == TRUE) {
plot(g,
     type = "o",
     col = "black",
     pch = 16,
     xlab = "Iteration",
     ylab = "Max Parameter Change",
     main = "Convergence Trace of SuSiE4X (Max |Δ| in alpha and beta)")

for (i in seq_along(g)) {
text(x = i,
     y = g[i],
     labels = formatC(g[i], format = "e", digits = 1),
     pos = 3,        #
     cex = 0.7,      #
     col = "red")
}
}
list(iter=iter,
   error=g,
   fitX = fitX,
   fitW = fitW,
   fitZ = fitZ,
   main_index=MainIndex,
   interaction_index=IntIndex)
}
