library(CppMatrix)
library(susieR)
library(devtools)
library(MASS)
library(RAMP)
library(hierNet)
library(glinternet)
load_all()
ARcov=function(p,rho){
s=c(1:p)
for(i in 1:p){
s[i]=rho^(i-1)
}
return(toeplitz(s))
}
CScov=function(p,rho){
diag(p)*(1-rho)+matrix(rho,p,p)
}
n=1000
p=10
R=kronecker(CScov(5,0.3),ARcov(2,0.5))
X=mvrnorm(n=n,mu=runif(p,0,1),R)
Xmean=colMeans(X)
Z=mvrnorm(n=n,mu=runif(15,0,1),diag(15))
alpha0=rnorm(15,0,1/sqrt(15))
beta0=rep(0,p)
ind=sample(p,3)
beta0[ind]=0.5
eta=matrixVectorMultiply(Z,alpha0)+matrixVectorMultiply(X,beta0)+X[,ind[1]]*X[,ind[2]]+X[,ind[1]]*X[,ind[3]]
y=eta+rnorm(n,0,2)
Lmain=5;Linteraction=5;max.iter = 50;min.iter = 3;max.eps = 1e-5;susie.iter = 300;verbose=T;crossprodX=NULL;n_threads=10

fit_SuSiE4X=SuSiE4X(X=X,Z=Z,y=y,Lmain=5,Linteraction=5,n_threads=10,strong=T)
fit_RAMP=RAMP(X=cbind(X,Z),y,penalty="MCP",n.lambda=50)
fit_glinternet=glinternet.cv(X=cbind(X,Z),Y=y,numLevels=rep(1,ncol(X)+ncol(Z)),nLambda=50,numCores=10,verbose=F)
