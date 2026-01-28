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
.sqrt<-function(x){
  sqrt(abs(x))*sign(x)
}

ids = data.table::fread("ids_discovery_and_validation.csv.gz")
ids$Sequencing_Sample_ID = .slug(ids$Sequencing_Sample_ID)
meta_data = data.table::fread("NITRIC_Metadata_09112023.csv.gz")
p1 = c("death_28day","comp_outcome","vfd","los_picu","picu_lcos_any")

comp_index = match(p1,names(meta_data) )

## this may be different for different input
counts0= data.table::fread( "NITRIC_corrected_counts.csv.gz")
gene_names = counts0$gene_id
sample_names = .slug(names(counts0)[-1])
counts = t(counts0[,-1]) ## so that transcripts are in columns.  Removes first column since its gene ids
dimnames(counts)= list(sample_names, gene_names)
meta_data$Sequencing_Sample_ID = .slug(meta_data$Sequencing_Sample_ID)
mi1 = match(sample_names,meta_data$Sequencing_Sample_ID)
mi2 = match(sample_names, ids$Sequencing_Sample_ID)

ids1 = ids[mi2[!is.na(mi2)],,drop=F]
meta_data1 = meta_data[mi1[!is.na(mi1)],,drop=F];
rownames(meta_data1) = meta_data1$Sequencing_Sample_ID

levs = unique(meta_data1$Sampling_time)
pre_inds = which(meta_data1$Sampling_time==levs[1])
post_inds = which(meta_data1$Sampling_time==levs[2])

counts = counts[!is.na(mi1),,drop=F]
#counts = counts ^(0.5);  ## variance stabilising transformation of count data


meta_pre = meta_data1[pre_inds,]
meta_post = meta_data1[post_inds,]
mi3 = match(meta_pre$record_id, meta_post$record_id)

meta_pre1 = meta_pre[!is.na(mi3),]
meta_post1 = meta_post[mi3[!is.na(mi3)],]
rownames(meta_pre1) = meta_pre1$record_id
rownames(meta_post1) = meta_post1$record_id
counts_pre1 = counts[pre_inds,][!is.na(mi3),]
counts_post1 = counts[post_inds,][mi3[!is.na(mi3)],]
rownames(counts_pre1) = rownames(meta_pre1)
rownames(counts_post1) = rownames(meta_post1)

counts_diff = counts_post1 - counts_pre1

ids1_pre = ids1[pre_inds,][!is.na(mi3),]

subind = sort(comp_index)
y_pre=as.matrix(meta_pre1[,..subind,drop=F])
y_post =as.matrix(meta_post1[,..subind,drop=F])
y = y_post
rownames(y) = meta_pre1$record_id


discovery = which(ids1_pre$validation==0)
validation = which(ids1_pre$validation==1)
print(length(discovery))
print(length(validation))


train_d =  list(pre = counts_pre1[discovery,], post = counts_post1[discovery,], diff = counts_diff[discovery,])
#train_d1 = lapply(train_d, .sqrt)
#names(train_d1) = paste0("sqrt.",names(train_d))

test_d =  list(pre = counts_pre1[validation,], post = counts_post1[validation,], diff = counts_diff[validation,])
#test_d1 = lapply(test_d, .sqrt)
#names(test_d1) = paste0("sqrt.",names(test_d))

nitric_data = list(nme="NITRIC", 
            dataset=train_d,#c(train_d, train_d1),
              y=y[discovery,,drop=F])
nitric_validation= list(nme="NITRIC", 
                        dataset=test_d, #c(test_d, test_d1),
                        y=y[validation,,drop=F])



dim(nitric_data$y)

## NOW RUN FSPLS
pows = c(0.5,1,1.5,2.0)
#pows = 1
transform_y=getYTransform(pow = pows, offset=0.1, norm=1000, n_random=10)
nmes_ = names(nitric_data$dataset)
#SET UP FLAGS
flags = list(pthresh = 0.05, max=10,nrep=1,batch=0,topn=20,beam=1,all_v_all=F,  project=T,  stop_y="rand",x_transform=T,
             pheno_balance = T,transform_y = toJSON(transform_y), useoffset=T,useglmnet=T,loadPV=T,
            data_types = toJSON(list("pre"=grep("pre",nmes_,v=T),"all"=nmes_))
             
)
#optional - not set datatypes
flags$data_types=toJSON(list(all = grep("sqrt",nmes_, inv=T,v=T)))
flags$quantiles = '[0.001]' 

datasets = list(NITRIC=nitric_data)
datasets_val = list(validation = nitric_validation);

options("max_ordinal"=10)
#  flags[['transform']] =toJSON(getXTransform(c(seq(-1,-0.2,by=0.2),seq(0.1,0.9,by=.1), seq(1,2,by=.5))))#  '{"x" :"function(x) x", "log":"function(x) log1p(x)"}'


datasAll =datasEnv$new(datasets,flags=flags, hasNA=F) 
datasAll_validation =datasEnv$new(datasets_val,flags=flags, hasNA=F) 

datasAll$updateTransforms(toJSON(transform_y))

phens=datasAll$pheno()$all
#phens = phens[1];
phens
options("x_transform"="NA")
## FIND VARIABLES
vars_all = datasAll$select(phens, flags, verbose=T)
  #vars_all1=.extractFullVars(vars_all)


all_models = datasAll$makeAllModels(vars_all,flags=flags)


#  sigs$saveModels(all_models)
#all_models1 = sigs$loadModels(flags, phens ,"")
eval1 = datasAll$evaluateAllModels(all_models)

eval_validation = datasAll_validation$evaluateAllModels(all_models)
eval2 = rbind(eval1, eval_validation)



# sigs$saveEval(eval0, flags, phens, "")
#eval1 = sigs$loadEval(flags,phens,"")
#eval = subset(eval, model!="avg")
ggps1=.plotEval2(eval2,legend=T, grid1="pheno", grid0="measure",
                 shape_color=c("data","transf"),sep_by=c("cv_full","measure"), showranges=T,
                 scales="free",title =names(phens)[1], title1="pheno" ) #, grid="pheno~cv_full",showranges = F)
ggps1

