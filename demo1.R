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

library(writexl)  ## to save weights
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

.expandAllvAll<-function(meta, colnmes){
  rown = rownames(meta)
  res_all = data.frame(unlist(lapply(colnmes, function(colnme){
      y = meta[[colnme]]
      levs = levels(y)
      inds = 1:length(levs)
      names(inds) = levs[inds]
      unlist(lapply(inds[-1], function(ind){
        inds2 = 1:(ind-1)
        names(inds2) = levs[inds2]
        lapply(inds2, function(ind2){
           y1 = rep(NA,  length(y))
           y1[y==levs[[ind]]]=0
           y1[y==levs[[ind2]]]=1
           y1
        })
      }),rec=F)
  }),rec=F))
  rownames(res_all) = rown
  res_all
}
meta = .expandAllvAll(meta,"predicted_labels_fine")
ys=list(pbmc=meta)
datasets = list(pbmc=list(counts=counts1))
mats = lapply(datasets, function(d) lapply(d, function(d1).getSparseMatrices(d1, hasNA=F)))
flags = list(pthresh = 1e-5, max=100, nrep=1,batch=0, train=names(datasets)[1],topn=20,beam=1,verbose=T)
datas =datasEnv$new(NULL, ys,mats=mats,flags=flags) 
phens=datas$pheno()
## FIND VARIABLES
variables = datas$select(phens, flags)
variables1 = variables[unlist(lapply(variables, function(var) "full" %in% names(var$inds)))]
print(unlist(lapply(variables1[[length(variables1)]]$var, function(v) v[2])))

## FIT MODELS
all_models = list()
all_models = 
  datas$makeAllModels(variables[seq(1,length(variables),by=20)], phens, flags, all_models)



eval = datas$evaluateAllModels(all_models, phens, flags)

## GET WEIGHTS FROM FULL MODEL
final_model = .getFinalModel(all_models, target_size = "max")
model_weights = .extractWeights(final_model)
names(model_weights) = gsub("pbmc.","",names(model_weights))
combined = list(combined=cbind(model_weights[[1]][,1:2], data.frame(lapply(model_weights, function(mw)mw[,3]))))
outw = paste0("weights",max(eval$numvars),".xlsx");
write_xlsx(combined,outw)

##PLOT
#ggps = .plotEval1(eval, rename=F, len=1)
 ggps=.plotEval1(eval, rename=F, len=1, grid="pheno")

#for multinomial or ordinal
#ggps = .plotEval1(eval, rename=T,grid="cv~subpheno", sep="pheno", txtsize=1,len=0)
out = paste0("plot",max(eval$numvars),".pdf");
pdf(out,width=30,height=30)
for(ggp in ggps) print(ggp)
dev.off()


##VISUALISE PREDICTIONS
predictions =datas$extractPredictions(all_models[c(10,20,30,40,50)],phens, flags, CV = F);
.plotArea(predictions, rename=F, grid="pheno", p_incl = "Memory.B.cells")


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



