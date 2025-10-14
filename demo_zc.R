setOptions<-function(config){
  if(!file.exists(config))stop("need config file specifying USER,PASS and KEY")
  config = t(read.csv(config,sep="=",header=F))
  config1 = as.list(config[2,])
  names(config1) = config[1,]
  config1 = rev(config1)
  config1 = config1[!duplicated(names(config1))]
  options(config1)
  config1
}
##SET APPROPRIATE LIB PATHS 
opts = setOptions("~/.sfx")
rlib = opts$rlib
dbDir = opts$dbDir
github=opts$github
.libPaths(rlib)


options(bigmemory.allow.dimnames=TRUE)
##SHOULD BE RUN FROM WHERE GIT CLONED TO

library(jsonlite)
library(R6)
library(Matrix)
library(glmnet)
library(tidyr)
library(pROC); 
#library(wCorr)
library(MASS);
library(ggplot2)
library("binom") ## for plotting


#library(SeuratObject)

#library(writexl)  ## to save weights
#library(readxl)
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
#src1=grep(".R$",dir("./R",rec=T),v=T)
#invisible(try(lapply(paste("./R",src1,sep="/"), function(x) {print(x);source(x)})))

options("fspls.types"=
          fromJSON('{"gaussian": ["correlation","rms"],"binomial":["AUC","correlation"],"multinomial":"AUC","ordinal" : "AUC_all"}'))
setwd(paste0(github,"/streamformatics/R"))
source("libs.R")

### FROM ZHITAO
datasets_ = readRDS("/data/scratch/projects/punim1068/Zhitao/covid_pbmcl_day_list_fspls2.rds")
datasets1 = lapply(datasets_, function(d){
 d1 = d[[1]] 
 nme1 = unlist(lapply(dimnames(d1$y)[[1]], function(nme)strsplit(nme,"_D")[[1]][1]))
 nme2 = unlist(lapply(dimnames(d1$dataset[[1]])[[1]], function(nme)strsplit(nme,"_D")[[1]][1]))
 dimnames(d1$y)[[1]] = nme1
 dimnames(d1$dataset[[1]])[[1]] = nme1
 d1
})
vals = unique(sort(unlist(lapply(datasets1, function(d){
  dimnames(d$y)[[1]]
}))))
y = datasets1[[1]]$y
datas = lapply(datasets1, function(d) d$dataset[[1]])
datasets = list(comb=list(y=y, dataset=datas))

#datasets=datasets_$input_D_neg_1_pbmc_fspls2
transform_y=getYTransform(pow = c(1), offset=0.1, norm=1000, n_random=10)

flags = list(pthresh = 1.0, nrep=0,batch=1,
#             incl = "input_D_neg_1_pbmc_fspls2",
             topn=30,beam=1,all_v_all=T,x_transform=T, transform_y = toJSON(transform_y))
## MAKE THE FSPLS DATA OBJECT
datasAll =datasEnv$new(datasets,flags=flags)
phens=datasAll$pheno()
## FIND VARIABLES
variables = datasAll$select(phens$all, flags, verbose=T)
vars_all_full=.extractFullVars(variables); print(vars_all_full$variables)

print(variables$variables)
## FIT MODELS
all_models = datasAll$makeAllModels(variables, phens$all, flags)
eval1 = datasAll$evaluateAllModels(all_models)
subset(eval1, cv_full == "CV= TRUE FULL= TRUE" & measure=="AUC" & numvars==1)
ggps1=.plotEval2(eval1,legend=T, grid1=c("subpheno","pheno"), grid0="measure",
                 shape_color=c("data"),sep_by=c("cv_full"), showranges=T,
                 scales="free",title =names(phens)[1], title1="pheno" )
ggps1$`CV= FALSE FULL= TRUE`
ggps1$`CV= TRUE FULL= FALSE`
ggps1$`CV= TRUE FULL= TRUE`
object = list(datas=datas, vars_all = variables,all_models = all_models)
#saveRDS(object,"/data/scratch/projects/punim1068/Zhitao/object.rds")

####


# ## loads data
##FUNCTION TO CONVERT 

rds=readRDS("/data/scratch/projects/punim1068/Zhitao/covid_pbmcl_day_list_fspls2.rds")
datasets = list(pbmc=list(counts=counts1))
mats = lapply(datasets, function(d) lapply(d, function(d1).getSparseMatrices(d1, hasNA=F)))

datasets = rds$input_D_neg_1_pbmc_fspls2$input_D_neg_1_pbmc$dataset
ys = rds$input_D_neg_1_pbmc_fspls2$input_D_neg_1_pbmc$y
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



