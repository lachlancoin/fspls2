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
library(bigmemory)
library(bigalgebra)
library(biganalytics)

#library(readr);
#library(dplyr)
#library(readxl)
#library(base64enc)
#library(arrow)
#library(binom); 
#library(data.table)
#library(R.utils)
#library(httr);
#maxsize=1.0 * 1e9
##LOAD CODE
if(rev(strsplit(getwd(),"/")[[1]])[1]!="fspls2")stop("not in right directory")
src1=grep(".R$",dir("./R",rec=T),v=T)
invisible(try(lapply(paste("./R",src1,sep="/"), function(x) {print(x);source(x)})))

# ## loads data
##FUNCTION TO CONVERT TO RIGHT FORMAT FOR THIS VERSION
.convert<-function(data,mode="rna", expand=F){
  names(data) =gsub("x_data","data",tolower(names(data)))
  if(is.list(data$y)){
    data$y = as.matrix(data$y)
    if(nrow(data$y)!=nrow(data$data)) stop("!!")
    dimnames(data$y)[[1]] = dimnames(data$data)[[1]]
  }
  if(!is.matrix(data$y)) {
    data$y =  data.frame(matrix(data$y, dimnames = list(dimnames(data$data)[[1]],"y")))
        if(is.character(data$y[1,1])){
          if(expand){
            data$y = .makeMultiClass(data$y)
            
          }else{
            data$y[[1]] = factor(data$y[[1]])
          }
        }
  }
  dataset = list(data$data);
  names(dataset)=mode
  y = data$y
 
  list(dataset=dataset,y=y)
}

######from rawlinson paper
## get data from this repo https://github.com/dn-ra/FSPLS-publication-repo via git clone
  path="~/github/FSPLS-publication-repo/input"
  files = grep(".Rds",dir(path,f=T,rec=T),v=T)
  names(files) = lapply(files, function(f) gsub(".prepd.Rds","",rev(strsplit(f,"/")[[1]])[1]))
  data1 = .convert(data=readRDS(files[['coin_multiclass_data']]), mode="rna", expand=F)
  data2  =.convert( data=readRDS(files[['coin_validaiton_data']]), mode="rna", expand=F)
  datasets = list(cohort1=data1$dataset,cohort2 = data2$dataset)
  ys = list(cohort1=data1$y, cohort2 = data2$y)
  ##GOLUB
  data1 =  .convert(readRDS(files[['golub']]), mode="rna", expand=F)
  datasets = list(cohort=data1$dataset) 
  ys = list(cohort=data1$y) 
  
flags = list(pthresh = 1e-3, nrep=1,batch=0, train=names(datasets)[1],topn=20,beam=1)

## MAKE THE FSPLS DATA OBJECT
datas =datasEnv$new(datasets, ys,flags=flags) 


phens=datas$pheno()[1]
## FIND VARIABLES
variables = datas$select(phens, flags)

## FIT MODELS
all_models = datas$makeAllModels(variables, phens, flags)
eval = datas$evaluateAllModels(all_models, phens, flags)



#print(attr(eval0,"translate"))
#ggps = .plotEval1(eval, rename=F, len=1)
#for multinomial or ordinal
ggps = .plotEval1(eval, rename=F,grid="subpheno~cv", sep="pheno")
ggps



predictions =datas$extractPredictions(all_models,phens, flags, CV = F);
#aa=roc(predictions[[2]]$y, predictions[[2]]$X0)
.plotArea(predictions, rename=F)
ggps


### get projection
varnames = variables[[length(variables)]]$var
datas$getProjectedData(varnames)
