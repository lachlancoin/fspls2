##SET APPROPRIATE LIB PATHS 

##SHOULD BE RUN FROM WHERE GIT CLONED TO
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.1/")
options(bigmemory.allow.dimnames=TRUE)

library(jsonlite)
library(R6)
library(Matrix)
library(glmnet)
library(tidyr)
library(pROC); 
library(wCorr)
library(MASS);
library(ggplot2)
library(nnet)  ## for multinomial
library("binom") ## for plotting
#optional packages
#library(bigmemory)
#library(bigalgebra)
#library(biganalytics)

#library(readr);
#library(dplyr)
#library(readxl)
#library(base64enc)
#library(arrow)
#library(binom); 
#library(data.table)
#library(R.utils)
#library(httr);

##LOAD CODE
if(rev(strsplit(getwd(),"/")[[1]])[1]!="fspls2")stop("not in right directory")
src1=grep(".R$",dir("./R",rec=T),v=T)
invisible(try(lapply(paste("./R",src1,sep="/"), function(x) {print(x);source(x)})))

#rawl = .readRawlinsonData('coin_multiclass_data');
path="~/github/FSPLS-publication-repo/input"
print(dir(path,full=T,rec=T))
rawl = .readRawlinsonData(filenames=list(golub = 'coin_data/coin_multiclass_data.prepd.Rds'), path= path)
rawl = .readRawlinsonData(filenames=list(golub = 'ng_data/ng_counts.prepd.Rds'), path= path)


flags = list(pthresh = 1e-5, nrep=10,batch=0, train=names(rawl$datasets)[1],topn=100,beam=1)

## MAKE THE FSPLS DATA OBJECT
datas =datasEnv$new(rawl,flags=flags) 


phens=datas$pheno()[1]
## FIND VARIABLES
variables = datas$select(phens, flags)
print(variables)
## FIT MODELS
all_models = datas$makeAllModels(variables, phens, flags)
eval = datas$evaluateAllModels(all_models, phens, flags)

##PLOT
ggps = .plotEval1(eval, rename=F, len=1)
#for multinomial or ordinal
#ggps = .plotEval1(eval, rename=F,grid="subpheno~cv", sep="pheno")
ggps


##VISUALISE PREDICTIONS
predictions =datas$extractPredictions(all_models,phens, flags, CV = F);
#aa=roc(predictions[[2]]$y, predictions[[2]]$X0)
.plotArea(predictions, rename=F)


### get projection
varnames = variables[[length(variables)]]$var
projOut=datas$getProjectedData(varnames)

indices = match(lapply(varnames, function(v)v[2]), dimnames(datas$datas$golub$data$rna)[[2]])
variance_before = datas$getVariance()
variance_after =   lapply(projOut, function(p1){
  lapply(p1, function(p2){
    apply(p2,2,var,na.rm=T)
  })
})
head(sort(unlist(variance_before),decr=T))
head(sort(unlist(variance_after),decr=T))



