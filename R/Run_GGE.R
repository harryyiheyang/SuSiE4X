Run_GGE <- function(X, Z, y,crossprodX=NULL, Lmain, Linteraction, max.iter, min.iter, max.eps, susie.iter,verbose=T,n_threads=1,...) {

n <- length(y)

Z_Res=X_Res=etaX=etaZ=etaW=0
ZtZ = blockwise_crossprod(X=Z,n_threads=n_threads)
if(is.null(crossprodX)){
XtX = blockwise_crossprod(X=X,n_threads=n_threads)
}else{
XtX=crossprodX
}

fitX=NULL
g=c()
beta=0
alpha=0
gamma=0
Ymean=mean(y)

for (iter in 1:max.iter) {
gamma_prev <- gamma
beta_prev <- beta
alpha_prev <- alpha

## --- update etaZ ---
rZ <- y - etaW - etaX
Zty <- crossprod(Z, rZ)
alpha <- as.vector(solve(ZtZ, Zty))

## --- update etaX ---
rX <- y - etaW - etaZ
Xty <- as.vector(crossprod(X, rX))
yty4X <- sum(rX^2)
fitX <- susie_suff_stat(XtX = XtX, Xty = Xty, yty = yty4X, n = n, L = Lmain, max_iter = susie.iter, estimate_prior_method = "EM")
beta <- coef.susie(fitX)[-1]

## --- update etaW ---
rW <- y - etaX - etaZ
XCS <- matrixMultiply(X, t(as.matrix(fitX$alpha)))
colnames(XCS) <- paste0("Main_CS", seq_len(ncol(XCS)))
cs <- sort(get_active_indices(fitX))
if(length(cs)==0){
  stop("No credible set of marginal effects detected. SuSiE4X stops.")
}
XCS <- as.matrix(XCS[, cs, drop = FALSE])
G = cbind(1,Z,XCS)
W <- get_pairwise_interactions(XCS, Z)
GtG = blockwise_crossprod(X=G,n_threads=n_threads)
GtW = blockwise_crossprod(X=G,Z=W,n_threads=n_threads)
ProjPart = matrixMultiply(G,(solve(GtG)%*%(GtW)))
W_Res = W - ProjPart
WtW <- blockwise_crossprod(X=W_Res,n_threads=n_threads)
Wty <- as.vector(crossprod(W_Res, rW))
yty4W <- sum(rW^2)
fitW <- susie_suff_stat(XtX = WtW, Xty = Wty, yty = yty4W, n = n, L = Linteraction, max_iter = susie.iter, estimate_prior_method = "EM")

## Update WCS
cs_W <- sort(get_active_indices(fitW))
if(length(cs)>0){
WCS <- matrixMultiply(W, t(as.matrix(fitW$alpha)))
WCS <- as.matrix(WCS[, cs_W, drop = FALSE])
}else{
WCS=NULL
}

if(is.null(WCS)){
fit_final=lm(y~Z+XCS+WCS)
coefs=coef(fit_final)
alpha=coefs[1:(ncol(Z)+1)]
XCSbeta=coefs[(ncol(Z)+2):(ncol(Z)+ncol(XCS)+1)]
WCSbeta=coefs[-c(1:(ncol(Z)+ncol(XCS)+1))]
etaZ=matrixVectorMultiply(cbind(1,Z),alpha)
etaX=matrixVectorMultiply(XCS,XCSbeta)
etaW=matrixVectorMultiply(WCS,WCSbeta)
}else{
fit_final=lm(y~Z+XCS)
coefs=coef(fit_final)
alpha=coefs[1:(ncol(Z)+1)]
XCSbeta=coefs[(ncol(Z)+2):(ncol(Z)+ncol(XCS)+1)]
etaZ=matrixVectorMultiply(cbind(1,Z),alpha)
etaX=matrixVectorMultiply(XCS,XCSbeta)
etaW=0
}

## --- check convergence ---
errX <- sqrt(mean((beta - beta_prev)^2))
errZ <- sqrt(mean((alpha - alpha_prev)^2))
err = c(errX, errZ)
g[iter]=max(err)
if (verbose == TRUE){
cat(paste0("Iteration ",iter,"\n"))
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
   interaction_index=IntIndex,
   nameW=colnames(W))
}
