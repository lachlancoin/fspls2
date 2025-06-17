##SET APPROPRIATE LIB PATHS 

##SHOULD BE RUN FROM WHERE GIT CLONED TO
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.1/")
options(bigmemory.allow.dimnames=TRUE)
library(jsonlite)
library(R6)
library(Matrix)
library(glmnet)
library(tidyr)
library(MASS);
library(ggplot2)
library(pROC); 
library(nnet)  ## for multinomial
library("binom") ## for plotting
library(confintr)
#library(wCorr)
##LOAD CODE
##SHOULD RUN FROM FSPLS2 directory
if(rev(strsplit(getwd(),"/")[[1]])[1]!="fspls2")stop("not in right directory")
src1=grep(".R$",dir("./R",rec=T),v=T)
invisible(try(lapply(paste("./R",src1,sep="/"), function(x) {print(x);source(x)})))

if(TRUE){
rdats = grep("rdat",dir("data",full=T),v=T)
examples = lapply(rdats, function(rdat){
  load(rdat)
  .convert(data, max=nrow(data$data),factor=length(grep("multi", rdat))>0)
})
names(examples)=lapply(rdats, function(str)strsplit(str,"/")[[1]][2])
}else{
## LOAD EXAMPLE DATASETS
example_files= grep("_data.rds",dir("data", full=T),v=T)
names(example_files) = lapply(example_files, function(x)strsplit(x,"/")[[1]][2])
examples = lapply(example_files, readRDS)
}
datasets =examples[4]
options("fspls.types"=
          fromJSON('{"gaussian": ["correlation","var"],"binomial":"AUC","multinomial":"AUC","ordinal" : "AUC_all"}'))
runAll<-function(datasets, randomise=F,pthresh = 0.001){
  y = datasets[[1]]$y
  #y1=y[sample.int(nrow(y), nrow(y)),1]
  #datasets[[1]]$y = as.matrix(data.frame(list(y = y[,1], y1 = y1,y2=y[,1])))
  flags = list(pthresh = pthresh, nrep=10,batch=0, train=names(datasets)[1],topn=20,beam=1)
  ## MAKE THE FSPLS DATA OBJECT
  datas =datasEnv$new(datasets,flags=flags) 
  if(randomise) datas$randomise()
  phens=datas$pheno()
  
  ## FIND VARIABLES
  variables = datas$select(phens, flags,verbose=F)
  #names(variables)[1]="x"
  print(variables)
  ## FIT MODELS
  all_models = datas$makeAllModels(variables, phens, flags)
  eval = datas$evaluateAllModels(all_models, phens, flags)
##PLOT
  if(is.null(eval)) return (NULL)
  eval2 = .calcEval1(eval)
  ggps = .plotEval1(eval2, rename=F,len=1)
#for multinomial or ordinal
#ggps = .plotEval1(eval, rename=F,grid="subpheno~cv")
ggps
if(TRUE) return(ggps)


##VISUALISE PREDICTIONS
predictions =datas$extractPredictions(all_models,phens, flags, CV = F);
#aa=roc(predictions[[2]]$y, predictions[[2]]$X0)
.plotArea(predictions, rename=F,datas$datas[[1]]$family[1])


varnames = variables[[length(variables)]]$var
projOut=datas$getProjectedData(varnames)
variance_before = datas$getVariance()
variance_after =   lapply(projOut, function(p1){
  lapply(p1, function(p2){
    sparse_variance(p2)
  })
})

##PRINT VARIANCE BEFORE AND AFTER
for(j in 1:length(datas$datas)){
  for(k in 1:length(datas$datas[[j]]$data)){
    inds = match(unlist(lapply(varnames, function(v)v[2])),colnames(datas$datas[[j]]$data[[k]]))
    print("before")
    df = data.frame(cbind(variance_before[[j]][[k]][inds],variance_after[[j]][[k]][inds]))
    names(df) = c("before","after")
    print(df)
  }
}

}

ggp_all = lapply(1:length(examples), function(i) runAll(examples[i], randomise=F, pthresh = 0.001))
names(ggp_all) = names(examples)
#runAll(datasets[1])

