##SET APPROPRIATE LIB PATHS 

##SHOULD BE RUN FROM WHERE GIT CLONED TO
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.1/")
if(rev(strsplit(getwd(),"/")[[1]])[1]!="fspls2")stop("not in right directory")
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

library(SeuratObject)
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
src1=grep(".R$",dir("./R",rec=T),v=T)
invisible(try(lapply(paste("./R",src1,sep="/"), function(x) {print(x);source(x)})))

# ## loads data
##FUNCTION TO CONVERT 



pbmc = readRDS("data/pbmc_20K.rds")
counts = t(pbmc@assays$RNA$counts)
sums = colSums(counts)
counts1 = getSparseSubMatrix(counts, which(sums>10000))
variance = sparse_variance(counts1)

meta=pbmc@meta.data[,match(c( "predicted_labels_broad", "predicted_labels_fine"),names(pbmc@meta.data))]
for(k in 1:ncol(meta))meta[[k]] = factor(meta[[k]])
ys=list(pbmc=meta)
datasets = list(pbmc=list(counts=counts1))
mats = lapply(datasets, function(d) lapply(d, function(d1).getSparseMatrices(d1, hasNA=F)))
flags = list(pthresh = 1e-5, max=50, nrep=5,batch=0, train=names(datasets)[1],topn=20,beam=1,verbose=T)
datas =datasEnv$new(NULL, ys,mats=mats,flags=flags) 
phens=datas$pheno()
## FIND VARIABLES
variables = datas$select(phens, flags)
variables1 = variables[unlist(lapply(variables, function(var) "full" %in% names(var$inds)))]
print(unlist(lapply(variables1[[length(variables1)]]$var, function(v) v[2])))

## FIT MODELS
all_models = datas$makeAllModels(variables, phens, flags)
full_models = lapply(all_models, function(all_models1) all_models1[['full']])
full_models = full_models[unlist(lapply(full_models, length))>0]
model_size= unlist(lapply(names(full_models), function(x) length(strsplit(x,";")[[1]])))
final_model = full_models[[which.max(model_size)]]

model_weights = unlist(lapply(final_model, function(mod1){
  phen1= names(mod1$betas)
  names(phen1) = phen1
  lapply(phen1, function(p1){
    bet=mod1$betas[[p1]]
    varn = data.frame(t(data.frame(lapply(names(mod1$var_names), function(str)strsplit(str,"\\.")[1]))))
    names(varn)=c("type","var")
    df=cbind(varn,bet)
    df
  })
}),rec=F)

library(writexl)
write_xlsx(model_weights,"weights.xlsx")

eval = datas$evaluateAllModels(all_models, phens, flags)



##PLOT
#ggps = .plotEval1(eval, rename=F, len=1)
#for multinomial or ordinal
ggps = .plotEval1(eval, rename=T,grid="cv~subpheno", sep="pheno", txtsize=1,len=0)
pdf("plots.pdf",width=30,height=30)
for(ggp in ggps) print(ggp)
dev.off()


##VISUALISE PREDICTIONS
#predictions =datas$extractPredictions(all_models,phens, flags, CV = F);
#.plotArea(predictions, rename=F)


### get projection
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



