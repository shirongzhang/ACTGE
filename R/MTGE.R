#'  Multiple Testing for Gene-Environment Interactions (MTGE)
#'
#'  This function returns P values of the interaction effect between individual SNPs 
#'  and exposure status. P values include unadjusted P values, P values adjusted by MTGE and Sidak. Users need to input a set of genetic variants G, outcome Y, 
#'  exposure E, and a few covariates to asjust W and outcome type, "B" for binary and ="C" for continuous.
#'
#'
#'@param G  n by p matrix of p SNPs for n subjects
#'@param Y  length n vector of outcome for n subjects
#'@param E  length n vector of exposure for n subjects
#'@param W  n by q data frame of q covariates to ajust for n subjects
#'@param outcomeType  "B" for binary outcome and "C" for continuous outcome
#'@param maf_cut  SNPs with Minor Allele Frequency lower than this threshold 
#'would be removed from the set, by default 0, i.e no SNP removed
#'@param missing_cut  SNPs with missing rate higher than this cut-off value
#'would be removed from the set, by default 0.15. Missing values would be imputed by average
#'@return unadjusted and adjusted P values for each individual SNP
#'\item{p_unadjsted}{unadjusted P values}
#'\item{p_MTGE}{adjusted P values using MTGE}
#'\item{p_Sidak}{adjusted P values using Sidak}
#'@author Shirong Zhang  <shirongz@usc.edu>
#'@author Juan Pablo Lewinger  <lewinger@usc.edu>
#'@examples 
#'## Generate data for examples
#'set.seed(100)
#'G1=rbinom(2000,1,0.3)
#'G2=rbinom(2000,1,0.4)
#'G3=rbinom(2000,1,0.25)
#'G4=rbinom(2000,1,0.5)
#'G=cbind(G1,G2,G3,G4)
#'E=rbinom(2000,1,0.3)
#'Y=rbinom(2000,1,0.2)
#'
#'## Example without covariates to adjust
#'GE_multiple(G,Y,E,W=NULL,outcomeType="B")
#'
#'## Example with covariates to adjust
#'age=rnorm(2000,45,2.6)
#'bmi=rnorm(2000,11.2,2.2)
#'GE_multiple(G,Y,E,W=cbind.data.frame(age,bmi),outcomeType="B")
#'
#'## Example if users want to change default settings
#'Y=rnorm(2000,2,4)
#'GE_multiple(G,Y,E,W=cbind.data.frame(age,bmi),outcomeType="C",maf_cut=0.05,missing_cut=0.1)
#'@details The input G, Y, E, W need to be numerical, with same number of observations.Note that Y and E cannot
#'contain missing values, and please remove those observations with missing in Y or E. 
#'We suggest users impute G before using this function, though SNPs with missing rate lower than 
#'missing_cut would be imputed in the function with mean value for that SNP. 
#'@references Shirong Zhang, Juan Pablo Lewinger Multiple Testing for Gene-Environment Interactions
#'@export
GE_multiple=function(G,Y,E,W,outcomeType,maf_cut=0,missing_cut=0.15){
  
  library("MASS")
  library("mvtnorm")
  
  # stop if there is no G,Y or E input
  if(is.null(G)){stop("There is no G input")}
  if(is.null(Y)){stop("There is no Y input")}
  if(is.null(E)){stop("There is no Y input")}
  
  G=as.matrix(G); if(!is.null(W)){W=as.data.frame(W)}
  
  # stop if there is NA in Y or E. Users need to remove those observations
  if(any(is.na(Y))){stop("There is NA in Y, please get rid of those observations!")}
  if(any(is.na(E))){stop("There is NA in E, please get rid of those observations!")}
  
  # stop if there are not same number of observations in G, Y, E, W
  if(length(unique(c(length(Y),length(E),nrow(G),nrow(W))))>1){stop("There is not equal number of obervations in G,Y,E,W")}
  
  # remove high missing rate SNP
  miss.rate=colMeans(is.na(G))
  miss.SNP=which(miss.rate>missing_cut)
  n.miss=length(miss.SNP)
  if(n.miss>0) {G=G[,-miss.SNP];warning(paste(n.miss,"SNP removed due to high missing rate"))}
  if(ncol(G)==0){stop("There is no SNP left in the set")}
  
  # remove low frequency allele SNPs
  allele.freq=colMeans(G,na.rm=TRUE)/2
  maf=ifelse(allele.freq>0.5,1-allele.freq,allele.freq)
  n.low.freq=length(which(maf<maf_cut))
  if(n.low.freq>0){G=G[,-which(maf<maf_cut)];warning(paste(n.low.freq,"SNP removed due to low allele frequency"))}
  if(ncol(G)==0){stop("There is no SNP left in the set")}
  
  # impute missing SNP with mean
  impute=function(x){x[is.na(x)]=mean(x,na.rm=TRUE);return(x)}
  G=apply(G,2,impute)
  
  # remove monomorphic G or G by E
  is.monomorphic=function(x){ length(unique(x))==1}
  mono.SNP=which(apply(G,2,is.monomorphic)==TRUE )
  if(length(mono.SNP)!=0){G=G[,-mono.SNP,drop=FALSE]}
  mono.GE=which(apply(G*E,2,is.monomorphic)==TRUE )
  if(length(mono.GE)!=0){G=G[,-mono.GE,drop=FALSE]}
  
  n.mono=length(mono.SNP)+length(mono.GE)
  
  if(n.mono>0){warning(paste(n.mono,"SNP removed due to mornomorphic: no variation in main effect or G*E"))}
  # can further remove G with almost mornomorphic, need to see how to set cutoff
  if(ncol(G)==0){stop("There is no SNP left in the set")}
  if(ncol(G)==1){stop("There is 1 SNP left in the set")}
  
  nSNPs = ncol(G)
  Z_GE_MeanCov=GE_Cov(G,Y,E,W,outcomeType)
  effect_size1=as.matrix(Z_GE_MeanCov$Zbetas)
  if(length(effect_size1)==1){
    Pv=round(pnorm(abs(effect_size1),lower.tail=F)*2,digits=8)
    stop("There is only 1 SNP in the set. The P value for SNP*E is ",Pv)}
  
  
  Pvals <- pchisq(effect_size1^2, df=1, lower=FALSE)     
  p_Sidak=1-(1-Pvals)^nSNPs
  p=numeric(nSNPs)
  for(i in 1:nSNPs){
    p[i]=1-pmvnorm(lower=rep(qnorm(Pvals[i]/2),nSNPs) ,upper=rep(qnorm(1-Pvals[i]/2),nSNPs),sigma=Z_GE_MeanCov$betaCor)
  }
  
  bindP=cbind.data.frame(Pvals,p,p_Sidak)
  colnames(bindP)=c("p_unadjsted","p_MTGE","p_Sidak")
  if(!is.null(colnames(G))){rownames(bindP)=colnames(G)}
  
  return(bindP)
  
}

# new version solve the issue of error in tcrossprod.. [ge.index,ge.index] out of bound script. by adding 1e-12 to 
# omitted GE terms due to low frequency in glm
GE_Cov=function(G,Y,E,W,outcomeType=c("B","C")){ 
  SNPs=G # if want to remove below 3 lines and use G,Y,E directly, need to change "SNPs" in fits also
  y = Y
  exposure=E
  N = length(y)
  nSNPs = ncol(SNPs)
  #nConfound=ifelse(is.null(W),0,ncol(W))
  
  fits = vector('list', nSNPs)
  Covinf=vector("list",nSNPs)
  covBetas = matrix(0, ncol= nSNPs, nrow=nSNPs)
  varBetas = numeric(nSNPs) 
  
  #   if(!is.null(W)){
  #     for (i in 1:ncol(SNPs)) { 
  #       if(outcomeType=="B"){ fits[[i]]=glm(y~SNPs[,i]*exposure+as.matrix(W),family=binomial)}
  #       else if(outcomeType=="C"){fits[[i]]=glm(y~SNPs[,i]*exposure+as.matrix(W),family=gaussian)}
  #     }
  if(!is.null(W)){
    for (i in 1:ncol(SNPs)) { 
      if(outcomeType=="B"){ fits[[i]]=glm(y~SNPs[,i]*exposure+.,family=binomial,data=W)} # . contains all variables in W
      else if(outcomeType=="C"){fits[[i]]=glm(y~SNPs[,i]*exposure+.,family=gaussian,data=W)}
    }
  }else{
    for (i in 1:ncol(SNPs)) { 
      if(outcomeType=="B"){fits[[i]]=glm(y~SNPs[,i]*exposure,family=binomial)}
      else if(outcomeType=="C"){fits[[i]]=glm(y~SNPs[,i]*exposure,family=gaussian)}
    }
  }
  
  Zbetas= sapply(fits, FUN=function(fit)ifelse("SNPs[, i]:exposure" %in% rownames(summary(fit)$coefficients), summary(fit)$coefficients["SNPs[, i]:exposure",3],0))
  auxInfo <- lapply(fits, covInfo)
  
  #   # below can have problem when W contains factor variables with level>2
  #   for(k in 1:nSNPs){
  #     for(l in 1:nSNPs){ covBetas[k, l]=tcrossprod(auxInfo[[k]], auxInfo[[l]])[4+nConfound,4+nConfound]} # for main effects, use 2 rather than 4
  #   } 
  
  # sometimes W contains factor variables, nconfound can be more than 1, say, race=1,2,3,4, nconf=3, not 1
  nvar=numeric(nSNPs)
  for(g in 1:nSNPs){nvar[g]=nrow(summary(fits[[g]])$coef)}
  ge.index=max(nvar)
  
  geOmit=which(nvar<ge.index)
  if(length(geOmit)>0){
    for(h in 1:length(geOmit)){
      auxInfo[[geOmit[h]]]=rbind(auxInfo[[geOmit[h]]],1e-16)
      warning(paste("SNP", geOmit[h]," omitted for GE term in glm model due to low frequency!"))
    }
  }
  #ge.index=nrow(summary(fits[[1]])$coef)
  for(k in 1:nSNPs){
    for(l in 1:nSNPs){ covBetas[k, l]=tcrossprod(auxInfo[[k]], auxInfo[[l]])[ge.index,ge.index]} # for main effects, use 2 rather than 4
  } 
  
  covBetas = (covBetas + t(covBetas))/2 # symmetrize 
  corBetas = cov2cor(covBetas)
  
  list(betaCor=corBetas,Zbetas=Zbetas)
}


covInfo=function(fit){
  X <- model.matrix(fit)
  if (any(is.na(coef(fit)))) {X=X[,-which(is.na(coef(fit)))]}
  U <- residuals(fit, type="response")*X ### why multiply residuals by design matrix, to get score??
  # 1 by n                       n by 4, 4 parameters (Intercept) snp smk snp:smk
  # but the dimension is 2000 by 4
  iV <- vcov(fit) #4by4 #Calculate Variance-Covariance Matrix for a Fitted Model Object
  tcrossprod(iV, U)
}
