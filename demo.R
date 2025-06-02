
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.1/")
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
#library(GenomicRanges)
#library(dplyr)
#library(readxl)
#library(base64enc)
#library(arrow)
#library(binom); 
#
#

#library(data.table)
#library(R.utils)
#library(httr);
#maxsize=1.0 * 1e9
##LOAD CODE
src1=grep(".R$",dir("./R",rec=T),v=T)
invisible(try(lapply(paste("./R",src1,sep="/"), function(x) {print(x);source(x)})))



family="binomial"
#family="gaussian"
## -- load demo dataset -- ##
if(family=="gaussian"){
    ## gaussian dataset
    datf = "data/z.tbod.gaussian.11.dat.rdat"
    pv_thresh = 0.05
}else if(family=="binomial"){
    ## binomial dataset
    datf = "data/w.tbod.binomial.11.dat.rdat"
    pv_thresh = 0.01
}else if(family=="multinomial"){
    ## 5 categories multinomial dataset
    datf = "data/r.tbod.multinomial.totcat5.11.dat.rdat"
    if(categories==3)
    ## 3 categories multinomial dataset
    datf = "data/r.tbod.multinomial.totcat3.11.dat.rdat"
    pv_thresh = 0.001
    refit = TRUE
}
# ## loads data
##FUNCTION TO CONVERT TO RIGHT FORMAT FOR THIS VERSION
.convert<-function(data,mode="rna", expand=F){
  names(data) =gsub("x_data","data",tolower(names(data)))
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

###DEMO DATA 
load(datf)
data1 =.convert(data)
datasets = list(cohort1=data1$dataset) #,cohort2 = data2$dataset)
ys = list(cohort1=data1$y) #, cohort2 = data2$y)



######from rawlinson paper
if(FALSE){
  path="/home/unimelb.edu.au/lcoin/github/FSPLS-publication-repo/input"
  files = grep(".Rds",dir(path,f=T,rec=T),v=T)
  names(files) = lapply(files, function(f) gsub(".prepd.Rds","",rev(strsplit(f,"/")[[1]])[1]))
  data1 = .convert(data=readRDS(files[['coin_multiclass_data']]), mode="rna", expand=F)
  data2  =.convert( data=readRDS(files[['coin_validaiton_data']]), mode="rna", expand=F)
  datasets = list(cohort1=data1$dataset,cohort2 = data2$dataset)
  ys = list(cohort1=data1$y, cohort2 = data2$y)
}

flags = list(pthresh = 1e-2, nrep=5,batch=0, train="cohort1",topn=20,beam=1)


datas =datasEnv$new(datasets, ys,flags=flags) 


phens=names(ys$cohort1)

variables = datas$select(phens, flags)

all_models = datas$makeAllModels(variables, phens, flags)

eval = datas$evaluateAllModels(all_models, phens, flags)



#print(attr(eval0,"translate"))
ggps = .plotEval1(eval, rename=F, len=1)
#for multinomial or ordinal
#ggps = .plotEval1(eval0, "subpheno~cv", sep="pheno")
ggps



##GET ROC
full_models = getModels(all_models, "full")

predictions =datas$extractPredictions(all_models,phens, flags, CV = T);
#aa=roc(predictions[[2]]$y, predictions[[2]]$X0)
.plotArea(predictions, rename=F)

#
#ggplot(area_plot2, aes(x=knots, y=value))+geom_line()#+geom_point(aes(size=counts))

#+ggtitle(rdsfile)
