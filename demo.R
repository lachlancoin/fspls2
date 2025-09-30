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

library(bigmemory)
library(bigalgebra)
options(bigmemory.allow.dimnames=TRUE)
#library(wCorr)
##LOAD CODE
##SHOULD RUN FROM FSPLS2 directory
if(rev(strsplit(getwd(),"/")[[1]])[1]!="fspls2")stop("not in right directory")
src1=grep(".R$",dir("./R",rec=T),v=T)

invisible(try(lapply(paste("./R",src1,sep="/"), function(x) {print(x);source(x)})))

.runAll1<-function(datasAll, flags){
  options("x_transform"=T)  ## apply transform_y to x instead of y
  phens=datasAll$pheno()$all
  vars_all = datasAll$select(phens, flags, verbose=T)
  vars_all1=.extractFullVars(vars_all)
  all_models = datasAll$makeAllModels(vars_all,flags=flags,db=NULL)
  eval1 = datasAll$evaluateAllModels(all_models,db=NULL)
  
  ggps1=.plotEval2(eval1,legend=T, grid1="pheno", grid0="measure",
                   shape_color=c("data","transf"),sep_by=c("cv_full"), showranges=T,
                   scales="free",title =names(phens)[1], title1="pheno" ) #, grid="pheno~cv_full",showranges = F)
  ggps1
  list(vars_all = vars_all, all_models = all_models, ggps1 = ggps1, vars_all1 = vars_all1, eval = eval1)
}

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
datasets =examples[3]; pthresh = 0.05 ; randomise=F; duplicate=F
options("fspls.types"=
          fromJSON('{"gaussian": ["correlation","rms"],"binomial":"AUC","multinomial":"AUC","ordinal" : "AUC_all"}'))
#transform_y=list(x=c("function(y) y","function(y) y"),exp= c("function(y) sign(y)*y^2", "function(y) sign(y)*abs(y)^(1/2)"))

runAll<-function(datasets, randomise=F,pthresh = 0.001, duplicate=F,transform_y =getYTransform(n_random=1)){
                 
  y = datasets[[1]]$y
  flags = list(pthresh = pthresh, max=10,nrep=1,batch=0, train=names(datasets)[1],topn=20,beam=1,all_v_all=F,  project=T,
               stop_y="rand",bigmatrix=T,
               useoffset=T,useglmnet=T#quantiles = toJSON(c(0,0.1))
              )
  #flags[['transform']] =toJSON(getXTransform(c(seq(-1,-0.2,by=0.2),seq(0.1,0.9,by=.1), seq(1,2,by=.5))))#  '{"x" :"function(x) x", "log":"function(x) log1p(x)"}'
  #flags$transform=toJSON(getXTransform(seq(.1,1,by=.1)))
 #flags$transform=toJSON(getXTransform(1))
  pows=1
  pows = seq(-1,1.8,by=0.4)
  pows = pows[pows>0]
 flags$transform_y = toJSON(getYTransform(pows = pows ,n_random=10))
 flags$transform=NULL
  ## MAKE THE FSPLS DATA OBJECT
  if(duplicate){
      y2 = do.call(cbind,replicate(2,datasets[[1]]$y,simplify=F))
      y2[,2] = y2[,2] + rnorm(nrow(y2))
      dimnames(y2)[[2]] = paste("y",1:ncol(y2),sep=".")
      datasets[[1]]$y = y2
  }
  datasAll =datasEnv$new(datasets,flags=flags,convertToBigMatrix = T) 
  #datasAll =datasEnv$new(datasets,flags=flags,convertToBigMatrix = F) 
  
  r1 = .runAll1(datasAll_bm, flags)
  r2 = .runAll1(datasAll, flags )
  nme1 = names(r1$all_models$models)
  nme2 = names(r2$all_models$models)
  match(nme1,nme2)
  datasAll = datasAll_bm
  
 
  if(TRUE) return(ggps1)
  ggps1$`CV=avg`
  ggps1$`CV= FALSE FULL= TRUE`
  if(FALSE){ ## show cum plots
    predictions0 =datasAll$extractPredictions(all_models, CV = F, liab=F, data_nme = names(datasAll$datas)[[1]]);
    fam=datasAll$datas[[1]]$family[[1]]
    area_p = .getAreaPlot1(predictions0,families=fam)
    ###area_p0$value= funcstr1(area_p0$value) only for transformed values, should put this in extractPredictions
    #aa=roc(predictions[[2]]$y, predictions[[2]]$X0)
    ggp_pred0=.plotArea(area_p, rename=F)
    ggp_pred0
     # predictions0 =datasAll$extractPredictions(all_models,phens, flags, CV = F, liab=F, data_nme = names(datasAll$datas)[[1]]);
    #  ggp_pred0=.plotArea1(predictions0, rename=T,max_vars=44)
    #  ggp_pred0
  }


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

