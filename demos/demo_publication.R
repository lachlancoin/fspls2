##SET APPROPRIATE LIB PATHS 

##SHOULD BE RUN FROM WHERE GIT CLONED TO
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.6/")
options(bigmemory.allow.dimnames=TRUE)

{
  library(RColorBrewer)
  library(jsonlite)
library(R6)
library(Matrix)
library(glmnet)
library(tidyr)
library(pROC); 
#library(wCorr)
library(MASS);
library(ggplot2)
library(nnet)  ## for multinomial
  library(DBI); 
  library(RSQLite);
  library(cowplot)
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
}
##LOAD CODE
if(rev(strsplit(getwd(),"/")[[1]])[1]!="fspls2")stop("not in right directory")
src1=grep(".R$",dir("./R",rec=T),v=T)
invisible(try(lapply(paste("./R",src1,sep="/"), function(x) {print(x);source(x)})))

#rawl = .readRawlinsonData('coin_multiclass_data');
path="~/github/FSPLS-publication-repo/input"
print(dir(path,full=T,rec=T))


transform_y=getYTransform(pow = 1,  n_random=1, perm=F)
flags = list(pthresh = 0.05, max=600,nrep=10,batch=0,topn=50,beam=1,all_v_all=T,  project=T,  stop_y="rand",x_transform=T,
             checkRMSV=F,  ## for checking RMSV when building model to ensure its improving
             transform_y = toJSON(transform_y), useoffset=T,useglmnet=T,loadPV=T, angles_only=F, get_plots=F)
options("x_transform"="NA")

names = c('coin_data/coin_multiclass_data.prepd.Rds', 'ng_data/ng_counts.prepd.Rds','golub_data/golub.prepd.Rds', 'alvez_data/alvez_data.prepd.Rds')
nme = names[1]
rawl = .readRawlinsonData(filenames=list(golub =nme), path= path)
dir("/home/unimelb.edu.au/lcoin/github/FSPLS-publication-repo/output/")
#rds = readRDS("/home/unimelb.edu.au/lcoin/github/FSPLS-publication-repo/output/alvez_kfold_results_unweighted.Rds")


## MAKE THE FSPLS DATA OBJECT
vars = apply(rawl$golub$dataset$rna,2,var)
rawl$golub$dataset$rna = rawl$golub$dataset$rna[,vars>quantile(vars)[1]] ## remove low variance col


dbDir1="./"; dbDir="./"
datasH = lapply(names(rawl), function(n){
  dataH$new(rawl[[n]],     dbDir = dbDir1, flags=flags, useDB=T)
  })
names(datasH) = names(rawl)

analysis =analysisEnv$new(dbDir=dbDir, flags=flags) ;

analysis$clear_db(drop=T, exclude=c(), datasH=datasH)# this clears the attached dbs
phens=datasH[[1]]$pheno()$all
#phens[[1]] = phens[[1]][3]
#phens2 = list(phens[[1]][3]); names(phens2) = names(phens)
flags[['data_types']] =toJSON(names(datasH[[1]]$data$data))
dh = datasH[[1]] 
dh$update(phens, flags, force=T);
vars_all=dh$select(analysis, phens , flags, verbose=F, useDB=F)
comb_plot = attr(vars_all,"plots")
if(!is.null(comb_plot))traj_plots=plot_traj_all(comb_plot, y="cumulative",keep_best = 5,txtsize=8,step=10)
vars_all1=.extractFullVars(vars_all)
all_modelsh= dh$makeAllModels(vars_all,phens=phens,useDB=F, verbose=T)

eval1=  dh$evaluateAllModels(all_modelsh,phens = phens, useDB=F, verbose=T)
#}), addName="data")
#eval1w=eval1 %>% pivot_wider(names_from = pheno, values_from=mid)

ggps2=.plotEval2(eval1,legend=T, grid1=c("subpheno","pheno"), grid0=c("measure","beam"),linetype="beam", ##"full_model"
                 shape_color=c("data","transf"),sep_by=c("cv_full"), showranges=T,
                 scales="free",title =names(phens)[1], title1="pheno" ) #, grid="pheno~cv_full",showranges = F)




ggps2[[1]]


