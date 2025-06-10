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

## LOAD EXAMPLE DATASETS
example_files= grep("_data.rds",dir("data", full=T),v=T)
names(example_files) = lapply(example_files, function(x)strsplit(x,"/")[[1]][2])
examples = lapply(example_files, readRDS)

datasets =examples[2]
#datasets[[1]]$y =Matrix(datasets[[1]]$y, sparse=T)

flags = list(pthresh = 5e-2, nrep=5,batch=0, train=names(datasets)[1],topn=20,beam=1)

## MAKE THE FSPLS DATA OBJECT
datas =datasEnv$new(datasets,flags=flags) 
phens=datas$pheno()[1]
## FIND VARIABLES
variables = datas$select(phens, flags)
#names(variables)[1]="x"
print(variables)
## FIT MODELS
all_models = datas$makeAllModels(variables, phens, flags)

options("fspls.types"=
          fromJSON('{"gaussian": ["correlation","var"],"binomial":"AUC","multinomial":"AUC_all","ordinal" : "AUC_all"}'))

eval = datas$evaluateAllModels(all_models, phens, flags)

##PLOT
eval2 = .calcEval1(eval,sep="subpheno")
ggps = .plotEval1(eval2, rename=F,len=1)
#for multinomial or ordinal
#ggps = .plotEval1(eval, rename=F,grid="subpheno~cv")
ggps


##VISUALISE PREDICTIONS
predictions =datas$extractPredictions(all_models,phens, flags, CV = F);
#aa=roc(predictions[[2]]$y, predictions[[2]]$X0)
.plotArea(predictions, rename=F,datas$datas[[1]]$family[1])




