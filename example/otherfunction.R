glinternet_evaulate=function(true_main_index,true_int_variable,fit_glinternet){
i_1Std <- which(fit_glinternet$lambdaHat1Std == fit_glinternet$lambda)
coefs <- coef(fit_glinternet$glinternetFit)[[i_1Std]]
tp_main=ifelse(length(setdiff(true_main_index,coefs$mainEffects$cont))==0,1,0)
tn_main=ifelse(length(setdiff(coefs$mainEffects$cont,true_main_index))==0,1,0)
if(is.null(coefs$interactions$contcont)){
tp_int=ifelse(is.null(true_int_variable),1,0)
tn_int=1
}else{
p=nrow(coefs$interactions$contcont)
int_variable=c()
for(i in 1:p){
int_variable[i]=paste0("X",coefs$interactions$contcont[i,1],"X",coefs$interactions$contcont[i,2])
}
tp_int=ifelse(length(setdiff(true_int_variable,int_variable))==0,1,0)
tn_int=ifelse(length(setdiff(int_variable,true_int_variable))==0,1,0)
}
g=data.frame(tp_main=tp_main,tn_main=tn_main,tp_int=tp_int,tn_int=tn_int)
return(g)
}


RAMP_evaulate=function(true_main_index,true_int_variable,fit_RAMP){
tp_main=ifelse(length(setdiff(true_main_index,fit_RAMP$mainInd))==0,1,0)
tn_main=ifelse(length(setdiff(fit_RAMP$mainInd,true_main_index))==0,1,0)
tp_int=ifelse(length(setdiff(true_int_variable,fit_RAMP$interInd))==0,1,0)
tn_int=ifelse(length(setdiff(true_int_variable,fit_RAMP$interInd))==0,1,0)
g=data.frame(tp_main=tp_main,tn_main=tn_main,tp_int=tp_int,tn_int=tn_int)
return(g)
}

SuSiE4X_evaulate=function(true_main_index,true_int_variable,fit_SuSiE4X){
int=int_tptn(true_int_variable, fit_SuSiE4X$main_index, fit_SuSiE4X$interaction_index)
main=main_tptn(true_main_index,fit_SuSiE4X$main_index)
g=data.frame(tp_main=main$tp,tn_main=main$tn,tp_int=int$tp,tn_int=int$tn)
return(g)
}

################################################################################
int_tptn <- function(true_int_variable, main_index, int_index) {

if(is.null(int_index)){
tn=1
tp=ifelse(is.null(true_int_variable),1,0)
return(list(tp=tp, tn=tn))
}else{

# Step 1: Create mappings from CS to variable names, and from variable names to CS.
cs_to_var <- split(main_index$Variable, main_index$CS)
var_to_cs <- setNames(main_index$CS, main_index$Variable)

# Step 2: Function to expand interactions from CS-level into individual variable pairs.
expand_cs_pair <- function(part1, part2) {
vars1 <- if(grepl("^Main_CS", part1)) cs_to_var[[part1]] else part1
vars2 <- if(grepl("^Main_CS", part2)) cs_to_var[[part2]] else part2

if(is.null(vars1) || is.null(vars2)) return(character(0))

unique(c(
paste0(rep(vars1, each=length(vars2)), "*", rep(vars2, times=length(vars1))),
paste0(rep(vars2, each=length(vars1)), "*", rep(vars1, times=length(vars2)))
))
}

# Step 3: Identify CS pairs supported by true interactions.
supported_cs_pairs <- c()
for (pair in true_int_variable) {
vars <- unlist(strsplit(pair, "\\*"))

# Determine corresponding CS; if environment variable (e.g., PM25), use directly.
cs1 <- if(vars[1] %in% names(var_to_cs)) var_to_cs[vars[1]] else vars[1]
cs2 <- if(vars[2] %in% names(var_to_cs)) var_to_cs[vars[2]] else vars[2]

# Add both directions for symmetry.
supported_cs_pairs <- c(supported_cs_pairs,
                      paste(cs1, cs2, sep="*"),
                      paste(cs2, cs1, sep="*"))
}

# Step 4: Identify unsupported CS pairs (possible sources of false positives).
unsupported_cs_pairs <- setdiff(int_index$Variable, supported_cs_pairs)

# Step 5: Expand unsupported CS pairs into individual interactions (potential false positives).
risky_interactions <- unlist(lapply(strsplit(unsupported_cs_pairs, "\\*"),
                                function(cs) expand_cs_pair(cs[1], cs[2])))

# Step 6: Expand all detected interactions for TP calculation.
detected <- unlist(lapply(strsplit(int_index$Variable, "\\*"),
                      function(cs) expand_cs_pair(cs[1], cs[2])))

# Step 7: Evaluate TP (true positive) and TN (true negative).
# TP: all true interactions must be detected.
tp <- as.integer(all(true_int_variable %in% detected))

# TN: no risky (unsupported) interactions should be detected.
tn <- as.integer(length(setdiff(risky_interactions, true_int_variable)) == 0)

return(list(tp=tp, tn=tn))
}
}
main_tptn=function(true_main_index,main_index){

# if there is no true main effect, tp = tn = 1 only if main_index is also null
if(is.null(true_main_index)){
tp=tn=ifelse(is.null(main_index)==1,1,0)
}

tp=ifelse(length(setdiff(true_main_index,main_index$Index))==0,1,0)
true_main_cs=unique(main_index$CS[match(true_main_index,main_index$Index)])
tn = ifelse(length(setdiff(main_index$Index, true_main_index)) == 0, 1, 0)

return(list(tp=tp,tn=tn))
}
