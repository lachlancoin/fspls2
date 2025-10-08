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


flags = list(pthresh = 0.2, nrep=10,batch=0, max=50,topn=100,beam=2,all_v_all=F, one_v_rest=F,x_transform=T)
#flags[['transform']] = '{"x" :"function(x) x","exp" :"function(x) exp(x)", "x3":"function(x) x^3","1x":"function(x) 1/x"}'
#flags[['transform']] = '{"x" :"function(x) x","log" :"function(x) log1p(x)"}'
#flags[['transform']] = '{"x" :"function(x) x","exp" :"function(x) exp(x)"}'
flags$transform_y = toJSON(getYTransform(pows = seq(0.2,2.0,by=0.4) ,offset=0.1, n_random=10,norm=10))

#flags[['transform']] = '{"x" :"function(x) x"}'

#flags[['transform']] =toJSON(getXTransform(c(seq(-1,-0.2,by=0.2),seq(0.1,0.9,by=.1), seq(1,2,by=.5))))#  '{"x" :"function(x) x", "log":"function(x) log1p(x)"}'


rawl = .readRawlinsonData(filenames=list(golub = 'coin_data/coin_multiclass_data.prepd.Rds'), path= path)
rawl = .readRawlinsonData(filenames=list(golub = 'ng_data/ng_counts.prepd.Rds'), path= path)
rawl = .readRawlinsonData(filenames=list(golub = 'golub_data/golub.prepd.Rds'), path= path)
rawl = .readRawlinsonData(filenames=list(golub = 'alvez_data/alvez_data.prepd.Rds'), path= path)
dir("/home/unimelb.edu.au/lcoin/github/FSPLS-publication-repo/output/")
rds = readRDS("/home/unimelb.edu.au/lcoin/github/FSPLS-publication-repo/output/alvez_kfold_results_unweighted.Rds")


## MAKE THE FSPLS DATA OBJECT
vars = apply(rawl$golub$dataset$rna,2,var)
rawl$golub$dataset$rna = rawl$golub$dataset$rna[,vars>quantile(vars)[1]] ## remove low variance cols
datasAll =datasEnv$new(rawl,flags=flags) 
datasAll$updateTransforms( toJSON(getYTransform(pows = seq(1) ,expX=T,n_random=10,norm=10)))

phens=datasAll$pheno(sep=F)
#apply(datasAll$datas$golub$y$multinomial.y,2,sum)
phens = phens[[1]]
#phens$all$binomial.multiway = phens$all$binomial.multiway[2]
## FIND VARIABLES
#transform_y=getYTransform(n_random=10)

vars_all = datasAll$select(phens, flags,verbose=T)
vars_all1=.extractFullVars(vars_all)

## FIT MODELS
options("fspls.verbose1"=T); options("fspls.CHECK"=F)
all_models = datasAll$makeAllModels(vars_all, phens, flags)
eval = datasAll$evaluateAllModels(all_models, phens, flags)
##PLOT
ggps1=.plotEval2(eval,legend=T, grid1="subpheno", grid0="measure",
                 shape_color=c("data","transf"),sep_by=c("cv_full"), showranges=T,
                 scales="free",title =names(phens)[1], title1="pheno" ) 

ggps1
##VISUALISE PREDICTIONS
predictions0 =datasAll$extractPredictions(all_models,phens, flags, CV = F, liab=F);

predictions1 =datasAll$extractPredictions(all_models,phens, flags, CV = T, liab=F);

#aa=roc(predictions[[2]]$y, predictions[[2]]$X0)
ggps = .plotArea1(predictions1, rename=F, max=10)
ggps

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



