depmapEnv<-R6Class("depmapEnv", public = list(
  drugs1="list",
  drugs2="character",
  depmap="character",
  
  initialize=function(depmap){
    self$depmap = depmap
    drugs = paste(depmap,"primary-screen-replicate-treatment-info.csv.gz",sep="/")
    drugs1 = read.delim(drugs,sep=",")
    drugs1 = drugs1[!duplicated(drugs1$broad_id),]
    self$drugs2 = unlist(lapply(drugs1$broad_id, function(str){
      paste(tolower(strsplit(str,"-")[[1]][1:5]),collapse=".")
    }))
    self$drugs1 = drugs1
  },
  getPhens=function(datasAll, topphens,sep=T){
    code = self$convertToCode(topphens)
    dn1= dimnames(datasAll$datas[[1]]$y$gaussian)[[2]]
    code1=unlist(lapply(code, function(c) grep(c, dn1,v=T)))
    
    phen1 = datasAll$pheno(sep=sep, sep_group=F, code=code1)
    phen1
  },
  convertFromCode=function(topphen){
    topphen1 = unlist(lapply(topphen, function(str) paste(strsplit(str,"\\.")[[1]][1:5],collapse=".")))
    self$drugs1$name [match(topphen1, self$drugs2)]
  },
  convertToCode=function(topphen){
    dr=self$drugs2[match(topphen, self$drugs1$name)]
    names(dr) =topphen
    if(length(is.na(dr))>0){
      print( topphen[is.na(dr)])
    }
    dr = dr[!is.na(dr)]
    code=unlist(lapply(dr,function(dr1)grep(dr1,self$drugs2,v=T)))
    code
  },
  importData=function(dist1, filenames=list(rna="CCLE_expression.csv.gz",cn="CCLE_gene_cn.csv.gz"  )){
    flags = list()
    depmap = self$depmap
    lapply(names(filenames), function(nme){
      print(nme)
      mat = data.frame(data.table::fread(paste(depmap,filenames[[nme]],sep="/"),sep=","))
      rownames(mat) = mat[,1]
      colnames(mat) = unlist(lapply(colnames(mat), function(str) strsplit(str,"\\..")[[1]][1]))
      #dimnames(mat)[[1]] = mat[,1]
      print("importing")
      dist1$importMat(mat[,-1], filenames[[nme]], opts$USER, flags, type=nme)
      print("done")
    })
  },
  importPheno=function(dist1, 
                       phenos1 = dist1$getPheno("drugs"),
                       nme1 = "primary-screen-replicate-collapsed-logfold-change.csv.gz"){
    depmap = self$depmap
    filenme1 = paste(depmap,nme1,sep="/")
    flags = list()
    flags[['sep']]=','; flags[['id']]="V2"
    phenos1$upload_pheno(filenme1,nme1,opts$USER,flags, samples=dist1$sampleID())
    phenos2 = dist1$getPheno("cells")
    nme2 = "sample_info.csv.gz"
    filenme2 = paste(depmap,nme2,sep="/")
    flags$slug_sample = '["function(x) x"]'
    flags$id = "DepMap_ID";flags$sep=","
    phenos2$upload_pheno(filenme2, nme2, opts$USER, flags, samples = dist1$sampleID())
  }
  
))
