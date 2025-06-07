get_consecutive <- function(max, number) {
start <- sample(1:(max - number + 1), 1)
return(start:(start + number - 1))
}
###############################################################################
demean=function(X){
Xname=colnames(X)
p=ncol(X)
if(is.null(p)){
X=X-mean(X)
}else{
Xmean=colMeans(X)
X=sweep(X,2,Xmean,"-")
}
colnames(X)=Xname
return(X)
}
###############################################################################
get_pairwise_interactions <- function(W, Z=NULL) {
n <- nrow(W)
p <- ncol(W)
q <- ncol(Z)
#
colnames_W <- colnames(W)
colnames_Z <- colnames(Z)
# --- W × W interactions ---
ww_cols <- p * (p + 1) / 2
WW <- matrix(NA, n, ww_cols)
colnames_WW <- character(ww_cols)
#
col_idx <- 1
for (i in 1:p) {
for (j in i:p) {
WW[, col_idx] <- W[, i] * W[, j]
colnames_WW[col_idx] <- paste0(colnames_W[i], "*", colnames_W[j])
col_idx <- col_idx + 1
}
}
# --- Z × W interactions ---
if(is.null(q)==F){
zw_cols <- q * p
ZW <- matrix(NA, n, zw_cols)
colnames_ZW <- character(zw_cols)
#
col_idx <- 1
for (i in 1:q) {
for (j in 1:p) {
ZW[, col_idx] <- Z[, i] * W[, j]
colnames_ZW[col_idx] <- paste0(colnames_Z[i], "*", colnames_W[j])
col_idx <- col_idx + 1
}
}
#
out <- cbind(WW, ZW)
colnames(out) <- c(colnames_WW, colnames_ZW)
}else{
out = WW
colnames(out)=colnames_WW
}
return(out)
}
###############################################################################
get_active_indices <- function(fit) {
cs = tryCatch(summary(fit), error = function(e) NULL)
if (!is.null(cs) && length(cs$cs) > 0) {
active_idx = unique(unlist(cs$cs$cs))
}else{
active_idx=NULL
}
return(active_idx)
}
###############################################################################
Identifying_MainEffect=function(fit,nam){
summ=summary(fit)$vars
g=unique(summ$cs[which(summ$cs>0)])
S=list()
for(i in g){
indi=which(summ$cs==i)
a=summ$variable[indi]
b=data.frame(Index=a,Variable=nam[summ$variable[indi]],CS=paste0("Main_CS",i))
S[[i]]=b
}
return(do.call(rbind,S))
}
###############################################################################
Identifying_IntEffect=function(fitW,namW){
summ=summary(fitW)$vars
if(length(which(summ$cs>0))>0){
g=unique(summ$cs[which(summ$cs>0)])
S=list()
for(i in g){
indi=which(summ$cs==i)
a=summ$variable[indi]
b=data.frame(Index=a,Variable=namW[summ$variable[indi]],CS=paste0("Int_CS",i))
S[[i]]=b
}
return(do.call(rbind,S))
}else{
return(NULL)
}
}

ProjectRes=function(A,B,inercept=F,n_threads){
if(inercept==T){
B=cbind(1,B)
}
BtB = blockwise_crossprod(X=B,n_threads=n_threads)
BtA = blockwise_crossprod(B,A,n_threads)
ProjPart = matrixMultiply(B,(solve(BtB)%*%(BtA)))
return(ProjPart)
}
