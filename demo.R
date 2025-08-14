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

##to store signatures
library(DBI); 
library(RSQLite);

library(cowplot)
#library(wCorr)
##LOAD CODE
##SHOULD RUN FROM FSPLS2 directory
if(rev(strsplit(getwd(),"/")[[1]])[1]!="fspls2")stop("not in right directory")
src1=grep(".R$",dir("./R",rec=T),v=T)
invisible(try(lapply(paste("./R",src1,sep="/"), function(x) {print(x);source(x)})))

if(FALSE){
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
datasets =examples[4]; pthresh = 0.001 ; randomise=F; duplicate=F
options("fspls.types"=
          fromJSON('{"gaussian": ["correlation","rms"],"binomial":"AUC","multinomial":"AUC","ordinal" : "AUC_all"}'))

func_y = list(y="function(y) y")  #'{"y":"function(y) y","y1":"function(y) y^2"} '

runAll<-function(datasets, randomise=F,pthresh = 0.001, duplicate=F,
                 func_y = list(y="function(y) y")){
  y = datasets[[1]]$y
  flags = list(pthresh = pthresh, max=10,nrep=10,batch=0, train=names(datasets)[1],topn=20,beam=1,all_v_all=F,  project=T,
               useoffset=T,useglm=T,#quantiles = toJSON(c(0,0.1)),
               transform_y =toJSON(func_y) )
  ## MAKE THE FSPLS DATA OBJECT
  if(duplicate){
  y2 = do.call(cbind,replicate(2,datasets[[1]]$y,simplify=F))
  y2[,2] = y2[,2] + rnorm(nrow(y2))
  dimnames(y2)[[2]] = paste("y",1:ncol(y2),sep=".")
  datasets[[1]]$y = y2
  }
  datasAll =datasEnv$new(datasets,flags=flags) 
  sigs = datasAll$getSigDB(reload=T)
  
  sigs$datas(); 
  #sigs$clear_all()
  if(randomise) datasAll$randomise()
  phens=datasAll$pheno()
  
  ## FIND VARIABLES
  vars_all = datasAll$select(phens, flags,verbose=T)
  
  sigs$saveVars(vars_all, phens, flags)
  sigs$experiments()
  
  vars_all1 = sigs$loadVars(flags, phens)
  print(vars_all1)
  options("fspls.verbose1"=T)
  
  all_models = datasAll$makeAllModels(vars_all, phens, flags)
  sigs$saveModels(all_models, flags,phens,"")
  all_models1 = sigs$loadModels(flags, phens ,"")
  eval0 = datasAll$evaluateAllModels(all_models1, phens, flags)
  
  sigs$saveEval(eval0, flags, phens, "")
  eval1 = sigs$loadEval(flags,phens,"")
  #eval = subset(eval, model!="avg")
  
  if(FALSE){ ## show cum plots
      predictions0 =datasAll$extractPredictions(all_models,phens, flags, CV = F, liab=F, data_nme = names(datasAll$datas)[[1]]);
      ggp_pred0=.plotArea1(predictions0, rename=T,max_vars=44)
      ggp_pred0
  }
  if(is.null(eval)) return (NULL)
  eval3 = .calcEval1(eval0)
  ggps = .plotEval1(eval3,  legend=T)
#for multinomial or ordinal
#ggps = .plotEval1(eval, rename=F,grid="subpheno~cv")
ggps
if(TRUE) return(ggps)


##VISUALISE PREDICTIONS
#predictions =datas$extractPredictions(all_models,phens, flags, CV = F);
#aa=roc(predictions[[2]]$y, predictions[[2]]$X0)
#.plotArea(predictions, rename=F,datas$datas[[1]]$family[1])


varnames = variables$variables[[length(variables$variables)]]
projOut=datasAll$getProjectedData(varnames)
variance_before = datasAll$getVariance()
variance_after =   lapply(projOut, function(p1){
  lapply(p1, function(p2){
    sparse_variance(p2)
  })
})

##PRINT VARIANCE BEFORE AND AFTER
varnames1 = unlist(lapply(varnames, function(v)match(v[[1]], names(datasAll$datas[[1]]$data))))
for(j in 1:length(datasAll$datas)){
  for(k in 1:length(datasAll$datas[[j]]$data)){
   varnames2 = varnames[varnames1==k]
    inds = match(unlist(lapply(varnames2, function(v)v[2])),colnames(datasAll$datas[[j]]$data[[k]]))
    print("before")
    df = data.frame(cbind(variance_before[[j]][[k]][inds],variance_after[[j]][[k]][inds]))
    names(df) = c("before","after")
    print(df)
  }
}

}

ggp_all =list()
#examples = examples[grep("multinomial", names(examples), inv=T)] # broken for multinomail
for(i in 1:length(examples)){  
  print(i)
  ggp_all[[i]] = runAll(examples[i], randomise=F, pthresh = 0.0001)
}
plot_grid(ggp_all[[1]][[1]], ggp_all[[2]][[1]], ggp_all[[3]][[1]],ggp_all[[4]][[1]])
names(ggp_all) = names(examples)
#runAll(datasets[1])

