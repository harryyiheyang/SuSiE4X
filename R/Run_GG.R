Run_GG <- function(X, y,crossprodX=NULL, Lmain, Linteraction, max.iter, min.iter, max.eps, susie.iter,verbose=T,n_threads=1,...) {

  n <- length(y)

  Xmean=colMeans(X)
  etaX=etaZ=etaW=0
  meanY=mean(y)
  if(is.null(crossprodX)){
    XtX = blockwise_crossprod(X=X,n_threads=n_threads)
  }else{
    XtX=crossprodX
  }
  fitX=NULL
  g=c()
  beta=0
  gamma=0

  for (iter in 1:max.iter) {
    gamma_prev <- gamma
    beta_prev <- beta

    ## --- update etaX ---
    rX <- y - etaW - meanY
    Xty <- as.vector(crossprod(X, rX))
    yty4X <- sum(rX^2)
    fitX <- susie_suff_stat(XtX = XtX, Xty = Xty, yty = yty4X, n = n, L = Lmain,
                            max_iter = susie.iter, estimate_prior_method = "EM",
                            X_colmeans=Xmean,y_mean=mean(rX))
    beta <- coef.susie(fitX)[-1]
    etaX = matrixVectorMultiply(X,beta)+coef.susie(fitX)[1]

    ## --- update etaW ---
    rW <- y - etaX - meanY
    XCS <- matrixMultiply(X, t(as.matrix(fitX$alpha)))
    colnames(XCS) <- paste0("Main_CS", seq_len(ncol(XCS)))
    cs <- sort(get_active_indices(fitX))
    if(length(cs)==0){
      stop("No credible set of marginal effects detected. SuSiE4X stops.")
    }
    XCS <- as.matrix(XCS[, cs, drop = FALSE])
    G = cbind(1,XCS)
    W <- get_pairwise_interactions(XCS)
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

    if(is.null(WCS)==F){
      fit_final=lm(y~XCS+WCS)
      coefs=coef(fit_final)
      meanY=coefs[1]
      XCSbeta=coefs[2:(ncol(XCS)+1)]
      WCSbeta=coefs[-c(1:(ncol(XCS)+1))]
      etaX=matrixVectorMultiply(XCS,XCSbeta)
      etaW=matrixVectorMultiply(WCS,WCSbeta)
    }else{
      fit_final=lm(y~XCS)
      coefs=coef(fit_final)
      meanY=coefs[1]
      XCSbeta=coefs[-1]
      etaX=matrixVectorMultiply(XCS,XCSbeta)
      etaW=0
    }

    ## --- check convergence ---
    err <- sqrt(mean((beta - beta_prev)^2))
    g[iter]=max(err)
    if (verbose == TRUE){
      cat(paste0("Iteration ",iter,"\n"))
    }
    if (max(err) < max.eps&iter>min.iter) {
      break
    }
  }

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
       fitJoint = fit_final,
       main_index=MainIndex,
       interaction_index=IntIndex,
       nameW=colnames(W))
}
