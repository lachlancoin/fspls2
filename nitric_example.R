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
setwd("~/github/fspls2_main")  ## to source code
src1=grep(".R$",dir("./R",rec=T),v=T)
invisible(try(lapply(paste("./R",src1,sep="/"), function(x) {print(x);source(x)})))
## set directory with data;
setwd("~/Data/Jingni")  ## where data is located

## function to help match strings
.slug<-function(str){
  gsub("-","_",str)
}
ids = data.table::fread("ids_discovery_and_validation.csv.gz")
ids$Sequencing_Sample_ID = .slug(ids$Sequencing_Sample_ID)
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
mi2 = match(sample_names, ids$Sequencing_Sample_ID)

ids = ids[mi2[!is.na(mi2)],,drop=F]
meta_data = meta_data[mi1[!is.na(mi1)],,drop=F];
rownames(meta_data) = meta_data$Sequencing_Sample_ID
counts = counts[!is.na(mi1),,drop=F]

subind = which(dimnames(meta_data)[[2]] %in% c("comp_outcome"))[1]
y=as.matrix(meta_data[,..subind,drop=F])
rownames(y) = rownames(meta_data)


discovery = which(ids$validation==0)
validation = which(ids$validation==1)





###optional normalisation
#counts = apply(counts, c(1,2), function(x) log(x+0.1))

nitric_data = list(nme="NITRIC", dataset = list(rna = counts[discovery,,drop=F]), y=y[discovery,,drop=F])
nitric_validation= list(nme="NITRIC", dataset = list(rna = counts[validation,,drop=F]), y=y[validation,,drop=F])





## NOW RUN FSPLS
#pows = c(0.5,1,1.5,2.0)
pows = c(1)
transform_y=getYTransform(pow = pows, offset=0.1, norm=1000, n_random=10)

#SET UP FLAGS
flags = list(pthresh = 0.05, max=10,nrep=1,batch=0,topn=20,beam=1,all_v_all=F,  project=T,  stop_y="rand",x_transform=T,
             pheno_balance = T,transform_y = toJSON(transform_y), useoffset=T,useglmnet=T,loadPV=T
             
)

flags$quantiles = '[0.5]' 

datasets = list(NITRIC=nitric_data)
datasets_val = list(validation = nitric_validation);


#  flags[['transform']] =toJSON(getXTransform(c(seq(-1,-0.2,by=0.2),seq(0.1,0.9,by=.1), seq(1,2,by=.5))))#  '{"x" :"function(x) x", "log":"function(x) log1p(x)"}'


datasAll =datasEnv$new(datasets,flags=flags, hasNA=F) 
datasAll_validation =datasEnv$new(datasets_val,flags=flags, hasNA=F) 

datasAll$updateTransforms(toJSON(transform_y))
phens=datasAll$pheno()$all
options("x_transform"="NA")
## FIND VARIABLES
vars_all = datasAll$select(phens, flags, verbose=T)
#vars_all1=.extractFullVars(vars_all)


all_models = datasAll$makeAllModels(vars_all,flags=flags)


#  sigs$saveModels(all_models)
#all_models1 = sigs$loadModels(flags, phens ,"")
eval1 = datasAll$evaluateAllModels(all_models)

eval_validation = datasAll_validation$evaluateAllModels(all_models)


# sigs$saveEval(eval0, flags, phens, "")
#eval1 = sigs$loadEval(flags,phens,"")
#eval = subset(eval, model!="avg")
ggps1=.plotEval2(eval_validation,legend=T, grid1="subpheno", grid0="measure",
                 shape_color=c("data","transf"),sep_by=c("cv_full"), showranges=T,
                 scales="free",title =names(phens)[1], title1="pheno" ) #, grid="pheno~cv_full",showranges = F)
ggps1

