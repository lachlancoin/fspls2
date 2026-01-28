.libPaths("~/R/x86_64-pc-linux-gnu-library/4.1/")


{
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
  ##library(bigmemory)
  library(data.table)
  #library(glmnet);
}
options(bigmemory.allow.dimnames=TRUE)
##SHOULD RUN FROM FSPLS2 directory
## this is just to source the scripts
setwd("~/github/fspls2")  ## to source code
src1=grep(".R$",dir("./R",rec=T),v=T)
invisible(try(lapply(paste("./R",src1,sep="/"), function(x) {print(x);source(x)})))
## set directory with data;
setwd("~/Data/Jingni")  ## where data is located

## function to help match strings
.slug<-function(str){
  gsub("-","_",str)
}

meta_data = data.table::fread("NITRIC_Metadata_09112023.csv.gz")
comp_index = match("comp_outcome",names(meta_data) )

## this may be different for different input
counts0= data.table::fread( "NITRIC_corrected_counts.csv.gz")
gene_names = counts0$gene_id
sample_names = .slug(names(counts0)[-1])
counts = t(counts0[,-1]) ## so that transcripts are in columns.  Removes first column since its gene ids
dimnames(counts)= list(sample_names, gene_names)
meta_data$Sequencing_Sample_ID = .slug(meta_data$Sequencing_Sample_ID)
mi1 = match(sample_names,meta_data$Sequencing_Sample_ID)
meta_data = meta_data[mi1[!is.na(mi1)],,drop=F];
rownames(meta_data) = meta_data$Sequencing_Sample_ID
counts = counts[!is.na(mi1),,drop=F]

subind = which(dimnames(meta_data)[[2]] %in% c("comp_outcome"))[1]
y=as.matrix(meta_data[,..subind,drop=F])
rownames(y) = rownames(meta_data)

nitric_data = list(nme="NITRIC", dataset = list(rna = counts), y=y)





## NOW RUN FSPLS
pows = c(0.5,1,1.5)
transform_y=getYTransform(pow = pows, offset=0.1, norm=1000, n_random=10)

#SET UP FLAGS
flags = list(pthresh = 0.05, max=10,nrep=10,batch=0,topn=10,beam=1,all_v_all=F,  project=T,  stop_y="rand",x_transform=T,
             pheno_balance = T,
             transform_y = toJSON(transform_y), useoffset=T,useglmnet=T,loadPV=T
)


datasets = list(NITRIC=nitric_data)
datasH = lapply(datasets, function(d)dataH$new(d,nme=d$nme, hasNA=F, convertToBigMatrix=F, flags=flags))
#lapply(datasH, function(dh)dh$clear_db(T))
nmesH = names(datasH); names(nmesH) = nmesH
datasAll =analysisEnv$new(flags=flags) 

phens=datasH[[1]]$pheno()$all
flags[['data_types']] =toJSON(names(datasH[[1]]$data$data))
lapply(datasH, function(dh) dh$update(phens, flags))
dh = datasH[[1]] 
vars_all=dh$select(datasAll, phens , flags, verbose=T)

all_models = lapply(datasH, function(dh){
  dh$makeAllModels(vars_all,useDB=F)
})
eval1= .merge1_new(lapply(nmesH, function(nmeh){
  all_modelsh = all_models[[nmeh]]
  dh = datasH[[nmeh]]
  dh$evaluateAllModels(all_modelsh, useDB=F)
}), addName="data")

ggps1=.plotEval2(eval1,legend=T, grid1=c("subpheno","pheno"), grid0="measure",linetype="beam", ##"full_model"
                 shape_color=c("data","transf"),sep_by=c("cv_full"), showranges=T,
                 scales="free",title =names(phens)[1], title1="pheno" ) #, grid="pheno~cv_full",showranges = F)
