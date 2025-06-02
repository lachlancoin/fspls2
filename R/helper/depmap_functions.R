
.getParams<-function(data_len){
  dt = NULL #rep(getOption("family","gaussian"), data_len) #unlist(  lapply(datas,.getDataType))
  fspls.max = as.numeric(getOption("fspls.max",20))
  fspls.beam = unlist(getOption("fspls.beam",c(2,4)))  # expand by 2 but only keep4
  fspls.beam = unlist(getOption("fspls.beam",c(1,1)))  # expand by 2 but only keep4
  
  fspls.batch =  as.numeric(getOption("fspls.batch",0))
  fspls.nrep = as.numeric(getOption("fspls.nrep",0))
  fspls.threshv= as.numeric(getOption("fspls.threshv",0.999))
  if(length(fspls.nrep)==1) fspls.nrep = rep(fspls.nrep[[1]],2)
  if(length(fspls.batch)==1) fspls.batch = rep(fspls.batch[[1]],2)
  
  list(maxvars=fspls.max, beam=fspls.beam, CI=T, nrep=fspls.nrep[[1]],
       family=dt,
       
       gaussians=rep("gaussians",data_len),
       mult = rep(1, data_len),
       reweight=getOption("reweight",F),
       max_null = getOption("fspls.max_null",0),
       type = getOption("fspls.method","slow"),
       types_=getOption("fspls.type",list(gaussian="correlation",binomial="AUC", 
                                          "ordinal"="rank_correlation", "multinomial"="AUC_all")),
       mean_thresh = getOption("fspls.mean_thresh",NULL),printAll=T,
       batch=fspls.batch[[1]], misclass=F, doplot=T,threshv=fspls.threshv,
       search_projection=getOption("search_projection",T),
       betas_projection= getOption("beta_projection",T))
}

##gens1 records NA

.getNormalisingFactors=function(gens,gens2,incl_var=F, incl_non_zero=F, set_NA_to_zero=F){
  ## 
  print("get norm")
  var=NULL
  non_zero=NULL
# na_len = NULL
  if(typeof(gens)=="S4"){
    norm = -1 * apply(gens, 2, function(g) sqrt(sum((g - mean(g, na.rm = T))^2, na.rm = T)))
    mean_x = apply(gens,2, function(g) mean(g, na.rm=T))
 #   na_len = apply(genes, 2, function(g) length(which(is.na(g))))
    if (incl_var) var = apply(gens, 2, function(g) var(g, na.rm = T))
    if (incl_non_zero) non_zero = apply(gens, 2, function(g) sum(g != 0, na.rm = T) / length(g))
    if(set_NA_to_zero){  ##only do if type double
     for(kk in 1:nrow(gens)){
       isna = is.na(gens[kk,])
       gens[kk,isna]=mean_x[isna]
       gens2[kk,isna]=1
     }
    }
  }else{
    norm = -1 *biganalytics::apply(gens,2, function(g) sqrt(sum((g-mean(g, na.rm=T))^2),na.rm=T))
   mean_x = -1 *biganalytics::apply(gens,2, function(g) mean(g, na.rm=T))
    if(incl_var) var = biganalytics::apply(gens,2, function(g) var(g,na.rm=T))
    if(incl_non_zero)nonzero = biganalytics::apply(gens,2, function(g) length(which(g!=0))/length(g))
   
  
   
  }
  
  list( norm=norm,var=var,non_zero=non_zero, mean_x = mean_x)
}


##FOLLOWING FUNCTIONS FOR RECURSIVELY FITTING RESIDUAL




###THIS READS A CSV FILE AS BIG MATRIX
##subsets are different subsets to read as matrices
##this for reading depmap data
.readCSVAsBigMatrix<-function(infile,subsets, head=T, 
                              transformations =fromJSON( '{
  "ans" : "function(x) 2*sqrt(x+3/8)",
                               "log" : "function(x) log1p(x)",
                               "bin" : "function(x) sqrt(1e6)*asin(sqrt(x/1e6))"
  }'),
                              sep=",",redo=F){
  print(paste("reading",infile))
  indir=paste(infile,".dir",sep="")
  if(redo && file.exists(indir)){
    removed=unlist(lapply(dir(indir,full=T), file.remove) )
    print("removed")
    print(removed)
    dir(dir(indir))
  }else if(!file.exists(indir)){
    dir.create(indir)
  }
  backingpath=indir
  rds=paste(indir,"file.rds",sep="/")
  if(redo)file.remove(rds)
  resu1=NULL
  trans_to_do = names(transformations)
  subsets_to_do = names(subsets)
  if(file.exists(rds)){
    resu1=(readRDS(rds))
     #subsets_to_do = subsets_to_do[!(subsets_to_do %in% names(resu1$zdesc))]
  #   bms = unlist(lapply(resu1$zdesc, function(zd){
  #       bm = try(attach.big.matrix(zd[[1]][['vals']]))
  #       inherits(bm,"try-error")
  #   }))
     subsets_to_do = subsets_to_do[which(!(subsets_to_do %in% names(resu1$zdesc)))]
   if(length(subsets_to_do)>0){
   #  resu1$zdesc= resu1$zdesc[which(!bms)]
  #   resu1$norm_factors= resu1$norm_factors[which(!bms)]
   }else{
     trans_to_do =   trans_to_do[!(trans_to_do %in%  names(resu1$zdesc[[1]]))]
   }
    # print("subsets_to_do")
    # print(subsets_to_do)
    if(length(trans_to_do)==0 && length(subsets_to_do)==0) return(resu1)
  }
  connection = gzfile(infile,"rb")
  line = readLines(connection,1)
  header=strsplit(line,split=sep)[[1]]
  ncols = length(header)-1
 
  nmes_s = if(length(subsets_to_do)>0) subsets_to_do else names(subsets)
  names(nmes_s)= nmes_s
  nmes_t = names(transformations)
  names(nmes_t) = nmes_t
  types = list(vals=options()$bigmemory.default.type, isNA="integer")
  matrices=lapply(nmes_s, function(subset_nme){
    subset = subsets[[subset_nme]]
    matrices1=lapply(nmes_t, function(trans_nme){
      if(!is.null(resu1) && !is.null(resu1$zdesc[[subset_nme]]) ){
       zd= resu1$zdesc[[subset_nme]][[trans_nme]]
        if(!is.null(zd)){
          print(paste("done", subset_nme, trans_nme))
         return(attach.big.matrix(zd))
        }
      }
      lapply(types, function(typ){
      backingfile = paste(trans_nme,subset_nme,typ,"bf",sep=".")
      descriptorfile = paste(trans_nme,subset_nme,typ,"df",sep=".")
      unlist(backingfile); unlist(descriptorfile)
      big.matrix(
        nrow = length(subset),
        ncol = ncols,
        dimnames=list(subset,header[-1]),
        #  ncol=length(mut_indices_list)+(end-num_to_rem),
        type = typ,
        init = if(typ=="integer") 0 else NA,
        separated = FALSE,
        backingfile =  backingfile,
        backingpath =  backingpath,
        descriptorfile = descriptorfile,
        binarydescriptor = FALSE,
        shared = options()$bigmemory.default.shared
      )
      })
    })
    matrices1
  })
  names(nmes_s)=NULL
  samp_all=  unlist(lapply(nmes_s, function(nme1){
    v1=subsets[[nme1]]
    v = rep(nme1, length(v1))
    names(v)= v1
    v
  }))
  
  i=1
  # nmes = rep("", nlines)
  
  ##fill matrix
  trans1 = lapply(transformations,  function(t) eval(str2lang(t)))
                    
  while(length(line)>0){
    if(i%%100==1)print(i)
    line = readLines(connection,1)
    if(length(line)==0){
      break;
    }
    row = strsplit(line,sep)[[1]]
    if(row[[1]] %in% names(samp_all)){
      type_nme=samp_all[[row[1]]]
    
        mut_matrix = matrices[[type_nme]]
      if(!is.null(mut_matrix)){
        i1=which(subsets[[type_nme]]==row[1])
        if(length(row)!=ncols+1) warning(paste("problem with row",i, length(row)-1, ncols))
        rowv = as.numeric(row[-1])
        for(nme1 in trans_to_do){  #names(mut_matrix)){
          mut_matrix[[nme1]][['vals']][i1,1:length(row)-1]=trans1[[nme1]](rowv)
          #mut_matrix[[nme1]][['isNA']][i1,1:length(row)-1][is.na(rowv)] =1
        }
      }
    }
    i=i+1
  }
  close(connection)
 # matrices2 = unlist(matrices, recursive=F)
  zdesc=lapply(matrices,function(m) lapply(m,function(m1) lapply(m1,describe)))
  if(!is.null(resu1) &&  length(subsets_to_do)>0){
    zdesc = c(resu1$zdesc, zdesc)
  }
  resu1 = list(zdesc = zdesc, rdsfile = rds)
  saveRDS(resu1, rds)
  return(resu1)
}
.addNormFactors<-function(resu1){
  zdesc = resu1$zdesc
  resu1$zdesc = lapply(zdesc, function(zd1){
    lapply(zd1, function(zd2){
      if(is.null(attr(zd2,"norm"))){
        mat = attach.big.matrix( zd2[['vals']])
        mat1 = attach.big.matrix( zd2[['isNA']])
        #list( norm=norm,var=var,non_zero=non_zero, mean_x = mean_x)
        norm_factors=.getNormalisingFactors(mat,mat1,set_NA_to_zero=T)
        attr(zd2,"norm")= norm_factors$norm
        attr(zd2,"mean_x")= norm_factors$mean_x
        attr(zd2,"var")= norm_factors$var
        attr(zd2,"non_zero")= norm_factors$non_zero
      }
      zd2
    })
  })
 
  saveRDS(resu1, resu1$rdsfile)
  return(resu1)
}
.getDepmapPhenotypes<-function(data_directory, inpfile_, sampfile, treat_file){
  .convert=function(x) gsub("\\.","-",strsplit(x,"\\.\\.")[[1]][1])
  inpfile = paste(data_directory,unlist(inpfile_) ,sep="/")
  names(inpfile) = names(inpfile_)
  sample_file = paste(data_directory,sampfile,sep="/")
  treatment_file = paste(data_directory,treat_file,sep="/")
  samples = read.csv(sample_file,head=T, row.names=1)  
  treatment=read.csv(treatment_file,head=T,row.names=1)
  phen=read.csv(inpfile[['pheno']] ,head=T, row.names=1)
  nme2Treatment=unlist(lapply(dimnames(phen)[[2]],.convert ))
  na_inds =is.na(treatment$name) & !is.na(treatment$broad_id)
  treatment$name[na_inds]=treatment$broad_id[na_inds]
  names(nme2Treatment) =treatment$name[ match(nme2Treatment,treatment$broad_id)]
  incl_samps = dimnames(phen)[[1]]  %in% dimnames(samples)[[1]]
  phen = phen[incl_samps,]
  incl_samps1 = dimnames(samples)[[1]] %in% dimnames(phen)[[1]]
  samples = samples[incl_samps1,]
  ##samples in same order as phens
  samples = samples[match(dimnames(phen)[[1]], dimnames(samples)[[1]]),]
  dimnames(phen)[[2]] = names(nme2Treatment)
  could_not_match = list()
  na_inds=which(is.na(names(nme2Treatment)))
  sampleID = dimnames(phen)[[1]]
  list(drugs = cbind(sampleID, phen), mat2 = cbind(sampleID, samples))
}
##if redo is true then it recalculates everything from scratch
.readAllDepmap<-function(data_directory, inpfile_, sampfile, treat_file,
                         incl_all = F,
                         transformations =fromJSON( '{
                                "x" : "function(x) x",
                              "ans" : "function(x) 2*sqrt(x+3/8)",
                               "log" : "function(x) log1p(x)",
                               "bin" : "function(x) sqrt(1e6)*asin(sqrt(x/1e6))"
  }'),
                         transformations_pheno =fromJSON( '{
                           "x" : "function(x) x"
                         }'),
                         min_thresh=10,
                         redo=F,types=NULL,
                         .convert=function(x) gsub("\\.","-",strsplit(x,"\\.\\.")[[1]][1])
                         
){
  rdsfile = paste(data_directory, "all.rds",sep="/")
  
  if(file.exists(rdsfile)){
    if(redo)file.remove(rdsfile)
    else return(readRDS(rdsfile))
  }
  inpfile = paste(data_directory,unlist(inpfile_) ,sep="/")
  names(inpfile) = names(inpfile_)
  sample_file = paste(data_directory,sampfile,sep="/")
  
  treatment_file = paste(data_directory,treat_file,sep="/")
  samples = read.csv(sample_file,head=T, row.names=1)  
  treatment=read.csv(treatment_file,head=T,row.names=1)
  phen=read.csv(inpfile[['pheno']] ,head=T, row.names=1)
  nme2Treatment=unlist(lapply(dimnames(phen)[[2]],.convert ))
  names(nme2Treatment) =treatment$name[ match(nme2Treatment,treatment$broad_id)]
  incl_samps = dimnames(phen)[[1]]  %in% dimnames(samples)[[1]]
  phen = phen[incl_samps,]
  incl_samps1 = dimnames(samples)[[1]] %in% dimnames(phen)[[1]]
  samples = samples[incl_samps1,]
  ##samples in same order as phens
  samples = samples[match(dimnames(phen)[[1]], dimnames(samples)[[1]]),]
 all_list= list(all=dimnames(samples)[[1]])
  if(!is.null(types)) samples = samples[samples$primary_disease %in% types,]
  samples$primary_disease = factor( samples$primary_disease)
  
  levs = levels(samples$primary_disease)
  
  names(levs) = levs
  subsets = lapply(levs,function(x){
    dimnames(samples)[[1]][which(samples$primary_disease==x)]
  })
  if(incl_all){
    subsets = c(subsets,all_list)
  }
  lens=unlist(lapply(subsets,length))
  ##minimum of 10 samples per subset
  subsets = subsets[lens>min_thresh]
  #dat_file="merged.dat"
  
  ## READ IN DATA FILES
  names(subsets)= gsub(" ","_",gsub("/",",", names(subsets)))
  ##following will be slow to read first time as it is creating the big.matrix, but fast after that
 ## assume first is pheno file
   pheno_files=lapply(inpfile[1], .readCSVAsBigMatrix,subsets, head=T,sep=",",redo=redo, transformations=transformations_pheno)  
   input_files=lapply(inpfile[-1], .readCSVAsBigMatrix,subsets, head=T,sep=",",redo=redo, transformations=transformations)  
   input_files = lapply(input_files, .addNormFactors) 
  resu = list(input_files = input_files, pheno_files = pheno_files,
              transformations = transformations,
              subsets = subsets,nme2Treatment = nme2Treatment, samples = samples, treatment=treatment)
    saveRDS(resu, file=rdsfile)
   resu
}


#data_directory="../../examples/FSPLS_external_validation/RAPIDS_sepsis"
.readAllRapids<-function(data_directory,subsets, 
                         inpfile_= list( expression="RAPIDS_Sepsis_htseq_stranded_count.csv"),
                         sep=",",
                         transformations =fromJSON( '{
  "x" :  "function(x) x",
  "ans" : "function(x) 2*sqrt(x+3/8)",
                               "log" : "function(x) log1p(x)",
                               "bin" : "function(x) sqrt(1e6)*asin(sqrt(x/1e6))"
  }'),
                         redo=F,
                         rdsname="all.rds",
                         .convert=function(x) gsub("^X","",gsub("_count.txt","",x))  ## this is conversion to apply to csv sample names so they match
                         
                         #.convert=function(x) gsub("\\.","-",strsplit(x,"\\.\\.")[[1]][1])
                         
){
  rdsfile = paste(data_directory, rdsname,sep="/")
  if(file.exists(rdsfile)){
    if(redo)file.remove(rdsfile)
    else return(readRDS(rdsfile))
  }
  inpfile = paste(data_directory,unlist(inpfile_) ,sep="/")
  names(inpfile) = names(inpfile_)
 
  #.readCSVAsBigMatrix1<-function(infile,subsets, head=T, sep=",",redo=F,
                                
  input_files=lapply(inpfile, .readCSVAsBigMatrix1,subsets,transformations=transformations,
                     head=T,sep=sep,redo=redo, .convert=.convert)  
  input_files = lapply(input_files, .addNormFactors) 
  
  resu = list(input_files = input_files, transformations = transformations,
              subsets = subsets,samples=samples, data=data_directory, 
              name= rev(strsplit(data_directory,"/")[[1]])[1]
  )
  saveRDS(resu, file=rdsfile)
  resu
}

##convert converts names in the csv  so to match the subse
.readCSVAsBigMatrix1<-function(infile,subsets, head=T, sep=",",redo=F,
                               transformations =fromJSON( '{
  "ans" : "function(x) 2*sqrt(x+3/8)",
                               "log" : "function(x) log1p(x)",
                               "bin" : "function(x) sqrt(1e6)*asin(sqrt(x/1e6))"
  }'),
                               .convert=function(x) gsub("^X","",gsub("_count.txt","",x))  ## this is conversion to apply to csv sample names so they match
                               
                                 ){
  print(paste("reading",infile))
  indir=paste(infile,".dir",sep="")
  if(!file.exists(indir))dir.create(indir)
  backingpath=indir
  rds=paste(indir,"file.rds",sep="/")
  if(redo )file.remove(rds)
  resu1=NULL
  trans_to_do = names(transformations)
  if(file.exists(rds)){
    resu1=(readRDS(rds))
    trans_to_do =   trans_to_do[!(trans_to_do %in%  names(resu1$zdesc[[1]]))]
    if(length(trans_to_do)==0) return(resu1)
  }
  
  csv=read.delim(infile, sep=sep)
  subrows = grep("__",csv$gene_id,inv=T) ## remove rows starting with __
  rownames=csv$gene_id[subrows]
  csv = t(csv[subrows,-1])
  dimnames(csv)[[2]] =rownames 
 # .conv=  eval(str2lang())
  dimnames(csv)[[1]]=   unlist(lapply(dimnames(csv)[[1]], .convert))
  if(length(which(duplicated(dimnames(csv)[[1]])))>0) stop(" duplicated after conversion")
  sums = apply(csv,1,sum)
  csv1 = 1e6*csv/sums
#  ncols = length(subrows)
  nmes_s = names(subsets)
  names(nmes_s)= nmes_s
  nme_trans=names(transformations)
  names(nme_trans)=nme_trans
  types = list(vals=options()$bigmemory.default.type, isNA="integer")
  matrices_all=lapply(nmes_s, function(subset_nme){
 #   matrices=lapply(nmes_s, function(subset_nme){
    samps = subsets[[subset_nme]]
    mi = match(samps, dimnames(csv1)[[1]])
    if(length(which(is.na(mi)))>0){
        print(samps[is.na(mi)])
        warning(" cannot match all samples")
    }
    samps = samps[!is.na(mi)]
    csv3 = csv1[mi[!is.na(mi)],,drop=F]
    matrices=lapply(nme_trans, function(trans){
      trans1=  eval(str2lang(transformations[[trans]]))
     
      lapply(types, function(typ){
        if(typ=="integer"){
          csv4 = array(0, dim =dim(csv3), dimnames = dimnames(csv3))
        }else{
          csv4 = apply(csv3,c(1,2), trans1)
        }
      if(!is.null(resu1) && !is.null(resu1$zdesc[[subset_nme]][[typ]]) ){
        zd= resu1$zdesc[[subset_nme]][[trans]]
        if(!is.null(zd)){
          print(paste("done", subset_nme, trans))
          return(attach.big.matrix(zd))
        }
      }
      print(paste(subset_nme,trans, typ))
      
      backingfile = paste(trans,subset_nme,typ,"bf",sep=".")
      descriptorfile = paste(trans,subset_nme,typ,"df",sep=".")
      if(redo){
        file.remove(paste(backingpath,backingfile,sep="/"));
        file.remove(paste(backingpath,descriptorfile,sep="/"));
      }
      if(file.exists(paste(backingpath,backingfile,sep="/"))) return(NULL)
      unlist(backingfile); unlist(descriptorfile)
      print(paste(backingfile, descriptorfile))
      as.big.matrix(
        csv4,
       type = typ,
       separated = FALSE,
       backingfile =  backingfile,
       backingpath =  backingpath,
       descriptorfile = descriptorfile,
       binarydescriptor = FALSE,
       shared = options()$bigmemory.default.shared
     )
      })
  })
    matrices
  })
  #matrices_all = unlist(matrices_all1, rec=F)
  zdesc=lapply(matrices_all,function(m) lapply(m,function(m1) lapply(m1,describe)))
  if(!is.null(resu1) &&  length(subsets_to_do)>0){
    zdesc = c(resu1$zdesc, zdesc)
  }
  resu1 = list(zdesc = zdesc, rdsfile = rds)
  saveRDS(resu1, rds)
  
  return(resu1)
}
.extractPheno<-function(depmapData, nme1="all"){
  .conv1=function(x) strsplit(x,"::")[[1]][1]
  phenos =attach.big.matrix(depmapData$pheno_files$pheno$zdesc[[nme1]]$x$vals)
  dimnames(phenos)[[2]] = unlist(lapply(dimnames(phenos)[[2]], .conv1))
  phenos = data.frame(phenos[])
  dn2 = gsub("\\.","-",dimnames(phenos)[[2]])
  dn1 = depmapData$nme2Treatment
  mi2=match(dn2,dn1)
  if(length(which(is.na(mi2)))==0){
    dimnames(phenos)[[2]] = names(dn1)[mi2]
  }
  phenos
}
.findCorrelatedDrugs<-function(drug, depmapData){
  #this finds correlated vars
  phenos=.extractPheno(depmapData)
  ind1 = match(drug,dimnames(phenos)[[2]])
  cors_all = sapply(1:ncol(phenos), function(i1) cor(phenos[,c(i1, ind1)], use="complete.obs")[1,2] )
  names(cors_all) = dimnames(phenos)[[2]]
  cors_all =sort(cors_all, decr=T)
  return(cors_all)
  dimnames(phenos)[[2]][which(abs(cors_all)>thresh)]
  unlist(lapply(dimnames(phenos)[[2]], function(n) names(depmapData$nme2Treatment)[n]))
  
  cor1 = cor(phenos, use="pairwise.complete")
  pheno_covar=unlist(lapply(phenos,function(v)cor(cbind(v, v1),use="pairwise.complete")[1,2]))
  pheno_covar=sort(pheno_covar,decr=T)[1:2]
  unlist(lapply(match(names(pheno_covar),dimnames(phenos)[[2]]), function(n) names(depmapData$nme2Treatment)[n]))
}



##if redo is true then it recalculates everything from scratch
.readAllSanger<-function(data_directory, inpfile_, sampfile,
                         incl_all = F,
                         transformations =fromJSON( '{
                                "x" : "function(x) x",
                              "ans" : "function(x) 2*sqrt(x+3/8)",
                               "log" : "function(x) log1p(x)",
                               "bin" : "function(x) sqrt(1e6)*asin(sqrt(x/1e6))"
  }'),
                         transformations_pheno =fromJSON( '{
                           "x" : "function(x) x"
                         }'),
                         min_thresh=10, ##original is 10, change to 4,contain at least 5 samples
                         redo=F,types=NULL,
                         .convert=function(x) gsub("\\.","-",strsplit(x,"\\.\\.")[[1]][1])
                         
){
  rdsfile = paste(data_directory, "all.rds",sep="/")
  
  if(file.exists(rdsfile)){
    if(redo)file.remove(rdsfile)
    else return(readRDS(rdsfile))
  }
  inpfile = paste(data_directory,unlist(inpfile_) ,sep="/")
  names(inpfile) = names(inpfile_)
  sample_file = paste(data_directory,sampfile,sep="/")
  
  samples = read.csv(sample_file,head=T, row.names=1)
  phen=read.csv(inpfile[['pheno']] ,head=T, row.names=1,check.names = FALSE)
  colnames(phen)=unlist(lapply(colnames(phen), function(x) gsub("\\.","-",strsplit(x,"\\.\\.")[[1]][[1]])))
  
  incl_samps = dimnames(phen)[[1]]  %in% dimnames(samples)[[1]]
  phen = phen[incl_samps,]
  incl_samps1 = dimnames(samples)[[1]] %in% dimnames(phen)[[1]]
  samples = samples[incl_samps1,]
  ##samples in same order as phens
  samples = samples[match(dimnames(phen)[[1]], dimnames(samples)[[1]]),]
  all_list= list(all=dimnames(samples)[[1]])
  if(!is.null(types)) samples = samples[samples$cancer_type %in% types,]
  samples$cancer_type = factor( samples$cancer_type )
  levs = levels(samples$cancer_type)
  names(levs) = levs
  subsets = lapply(levs,function(x){
    dimnames(samples)[[1]][which(samples$cancer_type==x)]
  })
  if(incl_all){
    subsets = c(subsets,all_list)
  }
  lens=unlist(lapply(subsets,length))
  ##minimum of 10 samples per subset
  subsets = subsets[lens>min_thresh]
  #dat_file="merged.dat"
  
  ## READ IN DATA FILES
  names(subsets)= gsub(" ","_",gsub("/",",", names(subsets)))
  ##following will be slow to read first time as it is creating the big.matrix, but fast after that
  ## assume first is pheno file
  pheno_files=lapply(inpfile[1], .readCSVAsBigMatrix,subsets, head=T,sep=",",redo=redo, transformations=transformations_pheno)
  input_files=lapply(inpfile[-1], .readCSVAsBigMatrix,subsets, head=T,sep=",",redo=redo, transformations=transformations)
  input_files = lapply(input_files, .addNormFactors)
  
  resu = list(input_files = input_files, pheno_files = pheno_files,
              transformations = transformations,
              subsets = subsets,nme2Treatment = NULL, samples = samples, treatment=NULL)
  saveRDS(resu, file=rdsfile)
  resu
}


