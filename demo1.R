#library(reticulate)

#use_python("/apps/easybuild-2022/easybuild/software/Compiler/GCCcore/11.3.0/Python/3.11.3/bin/python")
#library(anndata)
#setwd("/data/scratch/projects/punim1068/Zhitao")
#data = read_h5ad("PanCancer_igt_s9_fine_counts.h5ad", backed='r') 

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
library(readxl)
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

weights_old=read_xlsx("weights42.xlsx")

if(!is.null(weights_old)){
  counts1 = counts1[,dimnames(counts1)[[2]] %in% weights_old$var]
}
datasets = list(pbmc=list(counts=counts1))
mats = lapply(datasets, function(d) lapply(d, function(d1).getSparseMatrices(d1, hasNA=F)))
flags = list(pthresh = 1e-10, max=100, nrep=10,batch=0, train=names(datasets)[1],topn=20,beam=1,verbose=T,all_v_all=F)



datasAll =datasEnv$new(NULL, ys,mats=mats,flags=flags) 


phens=datasAll$pheno()
phens$all = phens$all[2]
if(!is.null(weights_old)){
  #quicker to use existing variables if provided
  vars_all = datasAll$convert(weights_old$var,phens)
}else{
  vars_all = datasAll$select(phens, flags, verbose=T)
}

## FIT MODELS
options("fspls.verbose1"=F); options("fspls.check"=F)
phens_ = phens
#phens_[[1]][[1]] = phens[[1]][[1]][10]
#all_models_ = datas$makeAllModels(vars_all, phens_, flags)
all_models = datasAll$makeAllModels(vars_all, phens, flags,verbose=T)

#betas_save = all_models_[[1]]$`counts.CD74;counts.FTH1;counts.CCL5;counts.HLA-DRB1;counts.RPS12;counts.LYZ;counts.GNLY;counts.NIBAN1;counts.IGHM;counts.BANK1`$all$full$pbmc$betas$binomial
#print(betas_save)
#compare=lapply(1:length(all_models[[1]]), function(i){
#  cbind(   all_models_[[1]][[i]]$all$full$pbmc$betas$binomial[,1],   all_models[[1]][[i]]$all$full$pbmc$betas$binomial[,10])
#})
#compare1=lapply(1:length(all_models[[1]]), function(i){
#  cbind(   all_models_[[1]][[i]]$all$full$pbmc$constants_proj$binomial[1],   all_models[[1]][[i]]$all$full$pbmc$constants_proj$binomial[10])
#})

#betas_save1 = all_models[[1]]$`counts.CD74;counts.FTH1;counts.CCL5;counts.HLA-DRB1;counts.RPS12;counts.LYZ;counts.GNLY;counts.NIBAN1;counts.IGHM;counts.BANK1`$all$full$pbmc$betas$binomial[,10,drop=F]
#print(betas_save1)
#cbind(betas_save, betas_save1)

options("diff_thresh"=0.1)

eval = datasAll$evaluateAllModels(all_models, phens, flags)
aa=subset(eval, pheno=="Non.classical.monocytes.cd8_Tm"  )
grep("Naive.B", unique(eval$pheno),v=T)
eval1 = subset(eval, numvars==100)

## GET WEIGHTS FROM FULL MODEL
final_models = .getFinalModel(all_models$y, target_size = "max")
model_weights = .extractWeights(final_models)
outw = paste0("weights",max(eval$numvars),".xlsx");
write_xlsx(model_weights,outw)

##PLOT
#ggps = .plotEval1(eval, rename=F, len=1)
eval1 = .calcEval1(eval)
 ggps2=.plotEval1(eval1,  grid0="pheno", showranges=T, scales="fixed", sep="cv_full", grid1="")
ggps2$`CV=avg`
#for multinomial or ordinal
#ggps = .plotEval1(eval, rename=T,grid="cv~subpheno", sep="pheno", txtsize=1,len=0)
out = paste0("plot_multinom",max(eval$numvars),".pdf");
pdf(out,width=30,height=30)
print(ggps2$`CV=avg`)
#for(ggp in ggps) print(ggps)
dev.off()


##VISUALISE PREDICTIONS
predictions =datasAll$extractPredictions(all_models,phens, flags,liab=F,CV = F);
len = length(predictions[[1]][[1]]$pbmc)

ggps2 = .plotArea1(predictions, rename=F,subset = len,
                   grid="pheno", p_incl = "Memory.B.cells")


### get projection
variables = vars_all$y$variables
varnames = variables[[length(variables)]]
projOut=datasAll$getProjectedData(varnames)



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



