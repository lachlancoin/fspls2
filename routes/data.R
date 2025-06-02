envir = environment()
#dbDir = getOption("dbDir","/app/db")
keys = envir$keys
dbDir = envir$dbDir
endpoint=envir$endpoint
endpoint="fspls/data.R"
all = allEnv$new(dbDir, keys, endpoint)
#all_dist = allEnv$new(dbDir, keys,"dist.R")


#* Load the data
#* @param bigmatrix  use big matrix if available
#* @post /load/<org>/<project>
function( req,org,project) {
  user=keys$validate(req,org,project,endpoint)
  if(is.null(user)) return(list(msg="unauthorised access"))
  multipart <- mime::parse_multipart(req)
  flags1 = .readFlags(multipart$flags,'{}')
 
  datas=all$get1(org,project,typesn=NULL, flags= flags1)
  if(is.null(datas)) return(list(msg="unauthorised access"))
  datas$dims()
}


#* Return phenotypes
#* @get /pheno/<org>/<project>
function( req,org,project) {
  user=keys$validate(req,org,db,endpoint)
  if(is.null(user)) return(list(msg="unauthorised access"))
  datas =all$get1( org,project)
  if(is.null(datas))return(list(msg="unauthorised access"))
  datas$pheno()
}
#* Return phenotypes
#* @get /categories/<org>/<project>
function( req,org,project) {
  user=keys$validate(req,org,project,endpoint)
  if(is.null(user)) return(list(msg="unauthorised access"))
  dbs = keys$dbs(org,project)
  keys$cats(org, project,dbs[[1]])
}


#* Return the score
#* @post /angles/<org>/<project>
function( req,org,project) {
  user=keys$validate(req,org,project,endpoint)
  if(is.null(user)) return(list(msg="unauthorised access"))
  datas =all$get1( org,project)
  if(is.null(datas))return(list(msg="unauthorised access"))
  multipart <- mime::parse_multipart(req)
  flags = fromJSON(multipart$flags)
  phens=fromJSON(multipart$phens)
  vars = fromJSON(multipart$vars)
  datas$angles(vars,phens,flags)  #  getAngles=function(train,vars,phens,incl,k, type,topn){
}





#* Return pvalue of last variable in context of preceding
#* @param vars  vars=fromJSON ({"rna_star.IFI27.rna_star.RTN2":{"rna_star.IFI27":["rna_star","IFI27"],"rna_star.RTN2":["rna_star","RTN2"]}} )
#* @param k fold replicate
#* @post /makeModels/<org>/<project>
function( req,org,project) {
  user=keys$validate(req,org,project,endpoint)
  if(is.null(user)) return(list(msg="unauthorised access"))
  datas =all$get1( org,project)
  if(is.null(datas))return(list(msg="unauthorised access"))
#  if(is.null(dim(datas[[1]]$y))) return (list(msg="need to define y matrix"))
  
  multipart <- mime::parse_multipart(req)
  phens=fromJSON(multipart$phens)
  flags = fromJSON(multipart$flags)
  if(!is.null(multipart$vars)){
    vars = fromJSON(multipart$vars)
  }else{
      f=multipart$upload$datapath
      vars=read_json(f)
  }
  
  datas$makeModels(vars,phens, flags)
}

##makes the model and does evaluation
#* @post /evaluateVars/<org>/<project>
#* @serializer tsv list(type="text/plain; charset=UTF-8")
function( req,org,project) {
  user=keys$validate(req,org,project,endpoint)
  if(is.null(user)) return(list(msg="unauthorised access"))
  
  datas =all$get1( org,project)
  if(is.null(datas))return(list(msg="unauthorised access"))
  multipart <- mime::parse_multipart(req)
  phens=fromJSON(multipart$phens)
  flags = fromJSON(multipart$flags)
  #f=multipart$upload$datapath
  if(!is.null(multipart$vars)){
    vars = fromJSON(multipart$vars)
  }else{
    f=multipart$upload$datapath
    vars=read_json(f)
  }
  all_models = datas$makeModels(vars,phens, flags)
  datas$evaluateModels(all_models,phens,flags)
}





#* @post /select/<org>/<project>
function( req,org,project) {
  user=keys$validate(req,org,project,endpoint)
  if(is.null(user)) return(list(msg="unauthorised access"))
  datas =all$get1( org,project)
  if(is.null(datas))return(list(msg="unauthorised access"))
  multipart <- mime::parse_multipart(req)
  phens=fromJSON(multipart$phens)
  flags = fromJSON(multipart$flags)
  datas$select(phens, flags)
}


#* @post /evaluate/<org>/<project>
#* @serializer tsv list(type="text/plain; charset=UTF-8")
function( req,org,project) {
  user=keys$validate(req,org,project,endpoint)
  if(is.null(user)) return(list(msg="unauthorised access"))
  
  datas =all$get1( org,project)
  if(is.null(datas))return(list(msg="unauthorised access"))
  multipart <- mime::parse_multipart(req)
  phens=fromJSON(multipart$phens)
  flags = fromJSON(multipart$flags)
  
  f=multipart$upload$datapath
  
  if(!is.null(f)){
    all_models=read_json(f,simplifyVector = T)
  }else{
    all_models = fromJSON(multipart$models)
  }
  datas$evaluateModels(all_models,phens,flags)
}



#* @post /plot/<org>
#* @serializer png
function( req,org) {
  user=keys$validate(req,org)
  if(is.null(user)) return(list(msg="unauthorised access"))
  multipart <- mime::parse_multipart(req)
  # return(names(multipart$upload))
  f=multipart$upload$datapath
  eval=read.delim(f,header=T,sep="\t")
  print(eval)
  ggp=.plotEval(eval)
  print(ggp)
}


#* Return pvalue of last variable in context of preceding
#* @param vars   which variables included, vars = {"rna_star.IFI27.rna_star.RTN2":{"rna_star.IFI27":["rna_star","IFI27"],"rna_star.RTN2":["rna_star","RTN2"]}} 
#* @param k fold replicate
#* @post /pvalue/<org>/<project>
function( req,org,project) {
  user=keys$validate(req,org,db,project)
  if(is.null(user)) return(list(msg="unauthorised access"))
  multipart <- mime::parse_multipart(req)
  phens=fromJSON(multipart$phens)
  flags = fromJSON(multipart$flags)
  vars = fromJSON(multipart$vars)
  datas =all$get1( org,project)
  if(is.null(datas))return(list(msg="unauthorised access"))
  datas$pvalues(vars,phens,flags)
}