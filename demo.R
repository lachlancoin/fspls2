
##SET APPROPRIATE LIB PATHS 
#.libPaths("~/R/x86_64-pc-linux-gnu-library/4.1/")
options(bigmemory.allow.dimnames=TRUE)
library(jsonlite)
library(R6)
library(Matrix)
library(glmnet)
library(tidyr)
library(MASS);
library(ggplot2)
library(pROC); 
library(binom) ## for plotting
library(confintr)

##to store signatures
library(DBI); 
library(RSQLite);
library(cowplot)
library(bigmemory)


##SHOULD RUN FROM FSPLS2 directory
if(rev(strsplit(getwd(),"/")[[1]])[1]!="fspls2")stop("not in right directory")
src1=grep(".R$",dir("./R",rec=T),v=T)
invisible(try(lapply(paste("./R",src1,sep="/"), function(x) {print(x);source(x)})))

## LOAD EXAMPLE DATASETS - USE RDS
if(FALSE){
  ##LOAD RDAT DATA, BUT THIS HAS NOT BEEN INCLUDED IN REPO
  rdats = grep("rdat",dir("data",full=T),v=T)
  examples = lapply(rdats, function(rdat){
    load(rdat)
    .convert(data, nme = rdat, max=nrow(data$data),factor=length(grep("multi", rdat))>0)
  })
  names(examples)=lapply(rdats, function(str)strsplit(str,"/")[[1]][2])
}else{
  example_files= grep("_data.rds",dir("data", full=T),v=T)
  names(example_files) = lapply(example_files, function(x)strsplit(x,"/")[[1]][2])
  examples = lapply(example_files, readRDS)
  for(k2 in 1:length(examples)) examples[[k2]]$nme = names(examples)[[k2]]
}



##HOW TO EVALUATE DIFFERENT DATASETS
options("fspls.types"= fromJSON('{"gaussian": ["correlation","rms"],"binomial":["AUC","area","max_diff","max_diff_x"],"multinomial":["AUC","area","max_diff","max_diff_x"],"ordinal" : "AUC_all"}'))

#SET UP TRANSFORMATIONS
pows = c(0.5,1,1.5)
transform_y=getYTransform(pow = pows, offset=0.1, norm=1000, n_random=10)

#SET UP FLAGS
flags = list(pthresh = 0.5, max=10,nrep=1,batch=0,topn=10,beam=2,all_v_all=F,  project=T,  stop_y="rand",x_transform=T,
             transform_y = toJSON(transform_y), useoffset=T,useglmnet=T,loadPV=T
)

datasets =examples[2];
#runAll<-function(datasets){
  datasH = lapply(datasets, function(d)dataH$new(d,nme=d$nme, flags=flags))
  nmesH = names(datasH); names(nmesH) = nmesH
  datasAll =datasEnv$new(datasH,flags=flags) 
  phens=datasAll$pheno()$all
  datasAll$update(phens, flags)
  
  
  vars_all0=datasH[[1]]$select(datasAll, phens , flags)
  #datasAll$clear_db(drop=T, exclude=c(), recursive=T)# this clears the attached dbs
  vars_all = datasAll$select(phens, flags, verbose=T,useDB=T)
  ##extract the variables for the full model only
  vars_all1=.extractFullVars(vars_all)
  all_models = lapply(nmesH, function(nmeh){
    dh=datasAll$datasH[[nmeh]]
   am=dh$makeAllModels(vars_all,useDB=F)
  })
  eval1= .merge1_new(lapply(nmesH, function(nmeh){
    all_modelsh = all_models[[nmeh]]
    dh = datasAll$datasH[[nmeh]]
    dh$evaluateAllModels(all_modelsh, useDB=F)
  }), addName="data")
  
  ggps1=.plotEval2(eval1,legend=T, grid1=c("subpheno","pheno"), grid0="measure",linetype="beam", ##"full_model"
                   shape_color=c("data","transf"),sep_by=c("cv_full"), showranges=T,
                   scales="free",title =names(phens)[1], title1="pheno" ) #, grid="pheno~cv_full",showranges = F)
  pdf("output.pdf")
  lapply(ggps1, print)
  dev.off()
  
  if(TRUE){ ## show cumulative plots for calibration
    ##extract all predictions
    predictions0 = lapply(nmesH, function(nmeh){
      dh=datasAll$datasH[[nmeh]]
      am=dh$extractPredictions(all_models[[nmeh]], CV = F, liab=F);
    })
      
    families=datasAll$datasH[[1]]$data$family[[1]]
    area_p = .getAreaPlot1(predictions0,families=families)
    ###area_p0$value= funcstr1(area_p0$value) only for transformed values, should put this in extractPredictions
    #aa=roc(predictions[[2]]$y, predictions[[2]]$X0)
    ggp_pred0=.plotArea(area_p, rename=F)
    ggp_pred0
  }
#}

#examples = examples[grep("multinomial", names(examples), inv=T)] # broken for multinomail
#for(i in 1:length(examples)){  
#  print(i)
#  ggp_all[[i]] = runAll(examples[i], randomise=F, pthresh = 0.0001)
#}
#plot_grid(ggp_all[[1]][[1]], ggp_all[[2]][[1]], ggp_all[[3]][[1]],ggp_all[[4]][[1]])
#names(ggp_all) = names(examples)
#runAll(datasets[1])

