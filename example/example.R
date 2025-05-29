library(CppMatrix)
library(susieR)
library(devtools)
library(MASS)
library(RAMP)
load_all()
data("EURblock")
#Rini=matrixMultiply(EURblock$U,t(EURblock$U)*EURblock$lambda)
ind=get_consecutive(dim(Rini)[1],1000)
R=Rini[ind,ind]
n=1000
X=mvrnorm(n=n,mu=rep(0,1000),R)
Xmean=colMeans(X)
Z=mvrnorm(n=n,mu=rep(0,15),diag(15))
#colnames(X) <- paste0("X", seq_len(ncol(X)))
#colnames(Z) = paste0("Z", seq_len(ncol(Z)))

alpha0=rnorm(15,0,1/sqrt(15))
beta0=rep(0,1000)
ind=sample(1000,3)
beta0[ind]=0.5
eta=matrixVectorMultiply(Z,alpha0)+matrixVectorMultiply(X,beta0)+X[,ind[1]]*X[,ind[2]]
y=eta+rnorm(n,0,1)
Lmain=5;Linteraction=5;max.iter = 50;max.eps = 1e-5;susie.iter = 300;verbose=T

fit_SuSiE4X=SuSiE4X(X=X,Z=Z,y=y,Lmain=5,Linteraction=5,n_threads=4)
fit_RAMP=RAMP(X=cbind(X,Z),y,penalty="MCP",n.lambda=50)
