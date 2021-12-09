#'  Adaptive Combination Tests for Set-Based Gene-Environment Interactions
#'
#'  This function returns P values of the interaction effect between a set of SNPs 
#'  and exposure status. P values include ACT, ACT weighted by screening statistics in both Z scale 
#'  and Chi-Square scale. Users need to input a set of genetic variants G, outcome Y, 
#'  exposure E, and a few covariates to asjust W and outcome type, "B" for binary and ="C" for continuous.
#'
#'
#'@param G  n by p matrix of p SNPs for n subjects
#'@param Y  length n vector of outcome for n subjects
#'@param E  length n vector of exposure for n subjects
#'@param W  n by q data frame of q covariates to ajust for n subjects
#'@param outcomeType  "B" for binary outcome and "C" for continuous outcome
#'@param GeneID  length p vector indicating GeneID for each SNP, i.e which Gene does each SNP belong to? If you don't want a two-level analysis, you can ignore it
#'@param maf_cut  SNPs with Minor Allele Frequency lower than this threshold 
#'would be removed from the set, by default 0, i.e no SNP removed
#'@param missing_cut  SNPs with missing rate higher than this cut-off value
#'would be removed from the set, by default 0.15. Missing values would be imputed by average
#'@param nReps  number of resampling to calculate ACT P Value, by default 10000. 
#'We recommend to use a small value, say 1000, to see the significance. If significant, then try larger value to get more accurate P Values.
#'@return nGene: number of Genes
#'@return P.Values: A list of ACT P values
#'\item{p.ACT.Z.weight.pos}{ACT P value based on positive part of marginal Z statistics weighted by screening statistics}
#'\item{p.ACT.Z}{ACT P value based on marginal Z statistics weighted by screening statistics}
#'\item{p.ACT.Z.weight}{ACT P value based on marginal Z statistics}
#'\item{p.ACT.X}{ACT P value based on marginal Chi-Square statistics}
#'\item{p.ACT.X.weight}{ACT P value based on marginal Chi-Square statistics weighted by screening statistics}
#'@return Z.Statistics: Z Statistics used to calculate ACT P Values
#'\item{GE.interaction}{Z statistics of marginal GE interaction}
#'\item{GY.association}{Z statistics of marginal GY association}
#'\item{GE.correlation}{Z statistics of marginal GE correlation}
#'\item{Weight}{screening Z statistics used as weight,stronger of GY association/GE correlation}
#'@return GE.interaction,GY.association,GE.correlation summary 
#'\item{Estimate}{Coefficient Estimate}
#'\item{Std.Error}{Standard Error}
#'\item{z value}{Z statistics}
#'\item{Pr(>|z|)}{P Value}

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
#'GeneID=c("Gene1","Gene1","Gene2","Gene3")
#'## Example without covariates to adjust
#'GE_ACT(G,Y,E,W=NULL,outcomeType="B",GeneID)
#'
#'## Example with covariates to adjust
#'age=rnorm(2000,45,2.6)
#'bmi=rnorm(2000,11.2,2.2)
#'GE_ACT(G,Y,E,W=cbind.data.frame(age,bmi),outcomeType="B")
#'
#'## Example if users want to change default settings
#'Y=rnorm(2000,2,4)
#'GE_ACT(G,Y,E,W=cbind.data.frame(age,bmi),outcomeType="C",maf_cut=0.05,missing_cut=0.1)
#'@details The input G, Y, E, W need to be numerical, with same number of observations.Note that Y and E cannot
#'contain missing values, and please remove those observations with missing in Y or E. 
#'We suggest users impute G before using this function, though SNPs with missing rate lower than 
#'missing_cut would be imputed in the function with mean value for that SNP. 
#'@references Shirong Zhang, Juan Pablo Lewinger Adaptive set-based methods for Gene-Environment Interactions
#'@export
GE_ACT=function(G,Y,E,W,outcomeType,GeneID=NULL,maf_cut=0,missing_cut=0.15,nReps=10000){
  
  library("MASS")
  library("mvtnorm")
  
  if(is.null(GeneID)){GeneID=rep("Gene",ncol(G))}
  if(length(GeneID)!=ncol(G)){stop("length of GeneID is not equal to number of SNPs")}
  GeneName=unique(GeneID);nGene=length(GeneName)
  
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
  if(n.miss>0) {G=G[,-miss.SNP];GeneID=GeneID[-miss.SNP];warning(paste(n.miss,"SNP removed due to high missing rate"))}
  if(ncol(G)==0){stop("There is no SNP left in the set")}
  
  # remove low frequency allele SNPs
  allele.freq=colMeans(G,na.rm=TRUE)/2
  maf=ifelse(allele.freq>0.5,1-allele.freq,allele.freq)
  low.SNP=which(maf<maf_cut)
  n.low.freq=length(low.SNP)
  if(n.low.freq>0){G=G[,-low.SNP];GeneID=GeneID[-low.SNP];warning(paste(n.low.freq,"SNP removed due to low allele frequency"))}
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
  
  simnull=SimNull1(nSNPs, nReps=nReps,corr=Z_GE_MeanCov$betaCor)
  
  screen=screening(Y,G,E,W)
  #sqrtR=Z_GE_MeanCov$sqrtR
  
  weight=screen$Bz
  
  ## Z version
  ZrealAndNull=cbind(effect_size1,simnull$Z)
  
  
  ## function body to return 2 level p values for pathway
  # empty matrix to store real and null P values from each gene
  ACT.wZ.pos=ACT.Z=weight.ACT.Z=ACT.X=weight.ACT.X=matrix(0,nGene,nReps+1)
  # matrix to return all 10 ACT p values for each gene
  p.Gene.ACT=matrix(0,5,nGene+2)
  colnames(p.Gene.ACT)=c("Pathway.1level","Pathway.2level",paste0("Gene_",unique(GeneID)))
  rownames(p.Gene.ACT)=c("p.ACT.Z.weight.pos", "p.ACT.Z", "p.ACT.Z.weight",  "p.ACT.X", "p.ACT.X.weight")

  
  for(i in 1:nGene){
    weight.gene=weight[which(GeneID==GeneName[i])]
    Z.gene=ZrealAndNull[which(GeneID==GeneName[i]),,drop=FALSE]
    outp=getP2(Z.gene,weight.gene)
    p.Gene.ACT[,2+i]=outp[,1]
    
    ACT.wZ.pos[i,]=outp["p.ACT.wZ.pos",]
    ACT.Z[i,]=outp["p.ACT.Z",]
    weight.ACT.Z[i,]=outp["p.weight.ACT.Z",]
    
    ACT.X[i,]=outp["p.ACT.X",]
    weight.ACT.X[i,]=outp["p.weight.ACT.X",]
  }
  
  # calculate 2 level p values
  
  p.ACT.wZ.pos.2=ACT(-log(ACT.wZ.pos))$p_best_p[1]
  
  p.ACT.Z.2=ACT(-log(ACT.Z))$p_best_p[1]
  p.weight.ACT.Z.2=ACT(-log(weight.ACT.Z))$p_best_p[1]
  
  p.ACT.X.2=ACT(-log(ACT.X))$p_best_p[1]
  p.weight.ACT.X.2=ACT(-log(weight.ACT.X))$p_best_p[1]
  
  path.2=c(p.ACT.wZ.pos.2,p.ACT.Z.2,p.weight.ACT.Z.2,p.ACT.X.2,p.weight.ACT.X.2)
  
  p.Gene.ACT[,2]=path.2
  
  # calculate 1 level pathway p value
  path.1=getP2(ZrealAndNull,weight)
  p.Gene.ACT[,1]=path.1[,1]
  
  Z.Statistics=cbind(effect_size1,screen$Mz,screen$Cz,screen$Bz)
  colnames(Z.Statistics)=c("GE.interaction","GY.association","GE.correlation","Weight")
  if(!is.null(colnames(G))){rownames(Z.Statistics)=colnames(G)}
  
  list (nGene=nGene,P.Values=p.Gene.ACT,Z.Statistics=Z.Statistics,GE.interaction=Z_GE_MeanCov$GE.inter,GY.association=screen$M,GE.correlation=screen$C)
}

getP2=function(ZrealAndNull,weight){
  # screenZ=screen$Bz;XrealAndNull;weight=screen$Bz^2
  
  # Z part
  p.ACT.Z=ACT(ZrealAndNull)$p_best_p
  p.weight.ACT.Z=ACT(ZrealAndNull*weight)$p_best_p
  wZ=ZrealAndNull*weight
  wZ.pos=wZ.neg=wZ; wZ.pos[which(wZ<0,arr.ind=TRUE)]=0
  p.ACT.wZ.pos=ACT(wZ.pos)$p_best_p
  
  # X part
  XrealAndNull=ZrealAndNull^2
  p.ACT.X=ACT(XrealAndNull)$p_best_p
  p.weight.ACT.X=ACT(XrealAndNull*(weight^2))$p_best_p
  return(rbind(p.ACT.wZ.pos,p.ACT.Z,p.weight.ACT.Z,p.ACT.X,p.weight.ACT.X))
}



# below GE_ACT does not have GeneID to do 2 level analysis, used for the package before 9/13/2016
# GE_ACT=function(G,Y,E,W,outcomeType,maf_cut=0,missing_cut=0.15,nReps=10000){
#   
#   library("MASS")
#   library("mvtnorm")
#   
#   # stop if there is no G,Y or E input
#   if(is.null(G)){stop("There is no G input")}
#   if(is.null(Y)){stop("There is no Y input")}
#   if(is.null(E)){stop("There is no Y input")}
#   
#   G=as.matrix(G); if(!is.null(W)){W=as.data.frame(W)}
#   
#   # stop if there is NA in Y or E. Users need to remove those observations
#   if(any(is.na(Y))){stop("There is NA in Y, please get rid of those observations!")}
#   if(any(is.na(E))){stop("There is NA in E, please get rid of those observations!")}
#   
#   # stop if there are not same number of observations in G, Y, E, W
#   if(length(unique(c(length(Y),length(E),nrow(G),nrow(W))))>1){stop("There is not equal number of obervations in G,Y,E,W")}
#   
#   # remove high missing rate SNP
#   miss.rate=colMeans(is.na(G))
#   miss.SNP=which(miss.rate>missing_cut)
#   n.miss=length(miss.SNP)
#   if(n.miss>0) {G=G[,-miss.SNP];warning(paste(n.miss,"SNP removed due to high missing rate"))}
#   if(ncol(G)==0){stop("There is no SNP left in the set")}
#   
#   # remove low frequency allele SNPs
#   allele.freq=colMeans(G,na.rm=TRUE)/2
#   maf=ifelse(allele.freq>0.5,1-allele.freq,allele.freq)
#   n.low.freq=length(which(maf<maf_cut))
#   if(n.low.freq>0){G=G[,-which(maf<maf_cut)];warning(paste(n.low.freq,"SNP removed due to low allele frequency"))}
#   if(ncol(G)==0){stop("There is no SNP left in the set")}
#   
#   # impute missing SNP with mean
#   impute=function(x){x[is.na(x)]=mean(x,na.rm=TRUE);return(x)}
#   G=apply(G,2,impute)
#   
#   # remove monomorphic G or G by E
#   is.monomorphic=function(x){ length(unique(x))==1}
#   mono.SNP=which(apply(G,2,is.monomorphic)==TRUE )
#   if(length(mono.SNP)!=0){G=G[,-mono.SNP,drop=FALSE]}
#   mono.GE=which(apply(G*E,2,is.monomorphic)==TRUE )
#   if(length(mono.GE)!=0){G=G[,-mono.GE,drop=FALSE]}
#   
#   n.mono=length(mono.SNP)+length(mono.GE)
#   
#   if(n.mono>0){warning(paste(n.mono,"SNP removed due to mornomorphic: no variation in main effect or G*E"))}
#   # can further remove G with almost mornomorphic, need to see how to set cutoff
#   if(ncol(G)==0){stop("There is no SNP left in the set")}
#   if(ncol(G)==1){stop("There is 1 SNP left in the set")}
#   
#   nSNPs = ncol(G)
#   Z_GE_MeanCov=GE_Cov(G,Y,E,W,outcomeType)
#   effect_size1=as.matrix(Z_GE_MeanCov$Zbetas)
#   if(length(effect_size1)==1){
#     Pv=round(pnorm(abs(effect_size1),lower.tail=F)*2,digits=8)
#     stop("There is only 1 SNP in the set. The P value for SNP*E is ",Pv)}
#   
#   simnull=SimNull1(nSNPs, nReps=nReps,corr=Z_GE_MeanCov$betaCor)
#   
#   screen=screening(Y,G,E,W)
# 
#   screenZ=screen$Bz
#   weight=screen$Bz
#   
#   ## Z version
#   ZrealAndNull=cbind(effect_size1,simnull$Z)
#   ACT.w.s.Z=getP(screenZ,ZrealAndNull,weight)
#   names(ACT.w.s.Z)=paste0(names(ACT.w.s.Z),".Z")
#   
#   wZ=ZrealAndNull*weight
#   wZ.pos=wZ.neg=wZ; wZ.pos[which(wZ<0,arr.ind=TRUE)]=0
#   p.ACT.wZ.pos=ACT(wZ.pos)$p_best_p[1]
#   
#   
#   ## X version
#   XrealAndNull=ZrealAndNull^2
#   ACT.w.s.X=getP(screenZ=screen$Bz,XrealAndNull,weight=screen$Bz^2)
#   names(ACT.w.s.X)=paste0(names(ACT.w.s.X),".X")
#   
#   P.Values=c(p.ACT.Z.weight.pos=p.ACT.wZ.pos,ACT.w.s.Z,ACT.w.s.X)
#   #return(P.Values)
#   Z.Statistics=cbind(effect_size1,screen$Mz,screen$Cz,screen$Bz)
#   colnames(Z.Statistics)=c("GE interaction","GY association","GE correlation","weight")
#   if(!is.null(colnames(G))){rownames(Z.Statistics)=colnames(G)}
#   list(P.Values=unlist(P.Values),Z.Statistics=Z.Statistics,GE.interaction=Z_GE_MeanCov$GE.inter,GY.association=screen$M,GE.correlation=screen$C)
# }

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
  
  #Zbetas= sapply(fits, FUN=function(fit)ifelse("SNPs[, i]:exposure" %in% rownames(summary(fit)$coefficients), summary(fit)$coefficients["SNPs[, i]:exposure",3],0))
  GE.inter=t(sapply(fits, FUN=function(fit){if("SNPs[, i]:exposure" %in% rownames(summary(fit)$coefficients)) summary(fit)$coefficients["SNPs[, i]:exposure",] else 0}))
  if(!is.null(colnames(G))){rownames(GE.inter)=colnames(G)}
  Zbetas=GE.inter[,3]
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
  
  list(betaCor=corBetas,Zbetas=Zbetas,GE.inter=GE.inter)
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

# screening function can be added a W matrix to adjust confounding factor
screening=function(Y,G,E,W=NULL){
  
  nSNPs=ncol(G)
  if(setequal(unique(Y),c(0,1))){outcomeType="B"}else{outcomeType="C"}
  if(setequal(unique(E),c(0,1))){exposureType="B"}else{exposureType="C"}
  
  # Marginal Y-G association Screening  
  M=matrix(0,nSNPs,4)
  fits = vector('list', nSNPs)
  
  if(!is.null(W)){
    for (i in 1:nSNPs) { 
      if(outcomeType=="B"){ fits[[i]]=glm(Y~G[,i]+.,family=binomial,data=W)}
      else if(outcomeType=="C"){fits[[i]]=glm(Y~G[,i]+.,family=gaussian,data=W)}
      if("G[, i]" %in% rownames(summary(fits[[i]])$coefficients)){
        #Mz[i]=summary(fits[[i]])$coefficients["G[, i]",3]
        M[i,]=summary(fits[[i]])$coefficients["G[, i]",]
      }
    }
  }else{
    for (i in 1:nSNPs) { 
      if(outcomeType=="B"){fits[[i]]=glm(Y~G[,i],family=binomial)}
      else if(outcomeType=="C"){fits[[i]]=glm(Y~G[,i],family=gaussian)}
      if("G[, i]" %in% rownames(summary(fits[[i]])$coefficients)){
        M[i,]=summary(fits[[i]])$coefficients["G[, i]",]
      }
    }
  }  
  
  # Marginal E-G association Screening  
  C=matrix(0,nSNPs,4)
  fits1 = vector('list', nSNPs)
  
  if(!is.null(W)){
    for (i in 1:nSNPs) { 
      if(exposureType=="B"){ fits1[[i]]=glm(E~G[,i]+.,family=binomial,data=W)}
      else if(exposureType=="C"){fits1[[i]]=glm(E~G[,i]+.,family=gaussian,data=W)}
      if("G[, i]" %in% rownames(summary(fits1[[i]])$coefficients)){
        C[i,]=summary(fits1[[i]])$coefficients["G[, i]",]
      }
    }
  }else{
    for (i in 1:nSNPs) { 
      if(exposureType=="B"){fits1[[i]]=glm(E~G[,i],family=binomial)}
      else if(exposureType=="C"){fits1[[i]]=glm(E~G[,i],family=gaussian)}
      if("G[, i]" %in% rownames(summary(fits1[[i]])$coefficients)){
        C[i,]=summary(fits1[[i]])$coefficients["G[, i]",]
      }
    }
  }  
  
  # get combined Zs for weights, which ever larger in abs value is selected
  colnames(M)=colnames(C)=c("Estimate","Std.Error","z value", "Pr(>|z|)")
  if(!is.null(colnames(G))){rownames(M)=rownames(C)=colnames(G)}
  Mz=M[,3];Cz=C[,3]
  Bz=Cz;Bz[which(abs(Cz)<abs(Mz))]=Mz[which(abs(Cz)<abs(Mz))]
  list(M=M,C=C,Mz=Mz,Cz=Cz,Bz=Bz)
}



# getP short has only weighted and non weighted Z X
getP=function(screenZ,XrealAndNull,weight){
  p.ACT=ACT(XrealAndNull)$p_best_p[1]
  p.weight.ACT=ACT(XrealAndNull*weight)$p_best_p[1]
  list(p.ACT=p.ACT,p.weight.ACT=p.weight.ACT)
}


# ACT.pos takes only positive X, usually is chisquared version XrealAndNull, 
# a p by n matrix, or length n vector if only 1 snp in set. returns p_best_p[1] 
# is p value of ACT.X of the real data, which is the first column of the X matrix
ACT.pos=function(X){
  # if X is a vector, means only 1 row, 1 snp
  if(!is.matrix(X) | nrow(X)==1){
    ecdf_X <- ecdf(X[2:length(X)]); p=1-ecdf_X(X);bestp=p
  }else{
    X=CumSum(apply(X,2,sort,decreasing=TRUE))
    ecdf_X <- Ecdf(X[,2:ncol(X)]) 
    p=1-t(sapply(1:nrow(X), function(i) ecdf_X[[i]](X[i,])))
    bestp =apply(p,2,min)
  }  
  ecdf_best_p<-Ecdf(bestp[2:length(bestp)])
  p_best_p=ecdf_best_p(bestp) 
  
  if(is.matrix(p)) {p_RTP=p[,1]} else {p_RTP=p[1]}
  list(p_best_p=p_best_p,p_RTP=p_RTP) 
}

# ACT can also return p.pos, p.neg
ACT=function(X){
  if(all(X>=0)) {ACT.pos(X)}
  else if(all(X<0)) {ACT.pos(-X)}
  else{    
    X.pos=X.neg=X; X.pos[which(X<0,arr.ind=TRUE)]=0; X.neg[which(X>=0,arr.ind=TRUE)]=0
    pos=ACT.pos(X.pos)$p_best_p; neg=ACT(-X.neg)$p_best_p; mix=ACT.pos(-log(rbind(pos,neg)))
    list(p_best_p=mix$p_best_p,p_RTP=mix$p_RTP,p_pos=pos[1],p_neg=neg[1])
  }
}  


SimNull1=function(nSNPs, nReps, corr, corr.type){
  dat0 <- GenData(nSNPs, nReps, rep(0, nSNPs), corr, corr.type)
  sumX0 <- CumSum(dat0$X) 
  ecdf_sumX0 <- Ecdf(sumX0)
  ####################################
  p=1-t(sapply(1:nrow(sumX0), function(i) ecdf_sumX0[[i]](sumX0[i,])))
  bestp = apply(p,2,min)
  ecdf_best_p<-Ecdf(bestp)
  list(Z=dat0$Z,Pvals=dat0$P,X0=dat0$X,sumX0=sumX0,sumX=ecdf_sumX0,ecdf_best_p0=ecdf_best_p)  # give ecdf0 list, can be used in simAlternative
}
GenData=function(nSNPs, nReps, effect_size, corr, corr.type){   
  if (is.matrix(corr)) corrMat = corr else corrMat = Covar(nSNPs, corr, corr.type) 
  # if corr is scalar construct correlation matrix from model
  zStats <- t(mvrnorm(nReps, effect_size, Sigma=corrMat))   
  # each replicates is an nReps column vector 
  chi2Stats <-apply(zStats^2, 2, sort, decreasing=TRUE)    
  # sorts chi-square stats within replicate
  Pvals <- pchisq(chi2Stats, df=1, lower=FALSE)     
  # pvals for the chi-squares (automatically sorted after sorting chi2's stats)
  list(Z=zStats, covZ=corrMat, X=chi2Stats, P=Pvals)
}

CumSum=function(stats) apply(stats, 2, cumsum)
Ecdf=function(stats) {if (is.matrix(stats))  apply(stats, 1, ecdf) else ecdf(stats)}
covInfo=function(fit){
  X <- model.matrix(fit)
  if (any(is.na(coef(fit)))) {X=X[,-which(is.na(coef(fit)))]}
  U <- residuals(fit, type="response")*X ### why multiply residuals by design matrix, to get score??
  # 1 by n                       n by 4, 4 parameters (Intercept) snp smk snp:smk
  # but the dimension is 2000 by 4
  iV <- vcov(fit) #4by4 #Calculate Variance-Covariance Matrix for a Fitted Model Object
  tcrossprod(iV, U)
}
