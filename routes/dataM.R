envir = environment()
#dbDir = getOption("dbDir","/app/db")
keys = envir$keys
dbDir = envir$dbDir
endpoint=envir$endpoint
#all = allEnv$new(dbDir, keys, endpoint)
#all_dist = allEnv$new(dbDir, keys,"dist.R")


#* Load the data /<org>/<project>
#* @param bigmatrix  use big matrix if available
#* @post /load/<org>
function( req,org) {
  user=keys$validate(req,org,NULL,endpoint)
  if(is.null(user)) return(list(msg="unauthorised access"))
  multipart <- mime::parse_multipart(req)
  urls1 = .readURLS(multipart$urls, self$endpoint1,"load")
  #  msim=all$get1(org,project,typesn=dbs,useBigMatrix=bigmatrix, reload=reload)
  # if(is.null(msim)) return(list(msg="unauthorised access"))
  lapply(urls1$urls, function(url){
    .POST1(req,url, body = list(flags = multipart$flags))
  })
  
}


#* Return phenotypes <org>/<project>
#* @get /pheno/<org>
function( req,org) {
  user=keys$validate(req,org,NULL,endpoint)
  if(is.null(user)) return(list(msg="unauthorised access"))
  multipart <- mime::parse_multipart(req)
  urls1 = .readURLS(multipart$urls, self$endpoint1,"pheno")
  #  msim=all$get1(org,project,typesn=dbs,useBigMatrix=bigmatrix, reload=reload)
  # if(is.null(msim)) return(list(msg="unauthorised access"))
  lapply(urls1$urls, function(url){
    .GET1(req,url)
  })
}
#* Return phenotypes  /<org>/<project>
#* @get /categories/<org>
function( req,org) {
  user=keys$validate(req,org,NULL,endpoint)
  if(is.null(user)) return(list(msg="unauthorised access"))
  multipart <- mime::parse_multipart(req)
  urls1 = .readURLS(multipart$urls, self$endpoint1,"categories")
  lapply(urls1$urls, function(url){
    .GET1(req,url)
  })
}


#* Return the score /<org>/<project>
#* @post /angles/<org>
function( req,org) {
  user=keys$validate(req,org,NULL,endpoint)
  if(is.null(user)) return(list(msg="unauthorised access"))
  multipart <- mime::parse_multipart(req)
  urls1 = .readURLS(multipart$urls, self$endpoint1,"angles")
  lapply(urls1$urls, function(url){
    .POST1(req,url, body =list(vars = multipart$vars, phens=multipart$phens, flags=multipart$flags))
  })
}





#* Return pvalue of last variable in context of preceding  /<org>/<project>
#* @post /makeModels/<org>
function( req,org) {
  user=keys$validate(req,org,NULL,endpoint)
  if(is.null(user)) return(list(msg="unauthorised access"))
  multipart <- mime::parse_multipart(req)
  if(!is.null(multipart$vars)){
    vars = multipart$vars
  }else{
      f=multipart$upload$datapath
      vars=toJSON(read_json(f))
  }
  urls1 = .readURLS(multipart$urls, self$endpoint1,"makeModels")
  lapply(urls1$urls, function(url){
    .POST1(req,url, body =list(vars = vars, phens=multipart$phens, flags=multipart$flags))
  })
}

##makes the model and does evaluation
#* @post /evaluateVars/<org>
#* @serializer tsv list(type="text/plain; charset=UTF-8")
function( req,org) {
  user=keys$validate(req,org,NULL,endpoint)
  if(is.null(user)) return(list(msg="unauthorised access"))
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
  urls1 = .readURLS(multipart$urls, self$endpoint1,"makeModels")
  urls2 = .readURLS(multipart$urls, self$endpoint1,"evaluate")
  nmes = names(urls1$urls)
  names(nmes)=nmes
  lapply(nmes, function(nme){
    url1 = urls1[[nme]]; url2 = urls2[[nme]]
    all_models = .POST1(req,url1, body =list(vars = vars, phens=multipart$phens, flags=multipart$flags))
    .POST1(req,url2, body =list(models = all_models, phens=multipart$phens, flags=multipart$flags))
  })
}




#* Learn variables /<org>/<project>
#* @param phens   phens=names(datas[[1]]$y)[1]  this trains separately
#* @post /select/<org>
function( req,org) {
  user=keys$validate(req,org,NULL,endpoint)
  if(is.null(user)) return(list(msg="unauthorised access"))
  multipart <- mime::parse_multipart(req)
  urls_angles = .readURLS(multipart$urls, self$endpoint1,"angles")
  urls_pvalues=.readURLS(multipart$urls, self$endpoint1,"pvalue")
  flags = .readFlags(multipart$flags)
  topn = .readFlag(flags,'topn', 20)
  return_type = .readFlag(flags,'return','model') ##model,eval,plot
  train_nme = .readFlag(flags,'train', names(datas1))
  test_nme = .readFlag(flags,'test', names(datas1))
  maxsize=.readFlag(flags,'max',50)
  num_pvals = min(topn, 10)
  #type='slow1'  #.readFlag(flags, 'type', 'slow1') 
  incl = .readFlag(flags,'data_types',names(datas1[[1]]$data))
  names(incl )= incl
  logpvthresh = log(.readFlag(flags,"pthresh",1e-5))
  vars_l = list(mStateObj$new(c(),c(), NULL))
  logpv=-100
  k = .readFlag(flags, 'rep',length(datas1[[1]]$train))
  beam = .readFlag(flags,"beam",2)
  while(logpv<logpvthresh && length(vars_l[[1]]$var)<maxsize){
    angles_all = lapply(vars_l, function(vars){
      
      angles=lapply(urls_angles, function(url){
        .POST1(req,url, body =list(vars = vars$var, phens=phens, incl=incl,k=k, type="slow1"))
      })
#      angles=lapply(train_nme, function(data_nme) datas1[[data_nme]]$getAngles1(phens,vars$var,incl=incl,k=k, type=self$type))
      comb=.summariseAngles(.combineAngles(angles,topn=topn),topn)
      
      lapply(urls_pvalues, function(url){
        .POST1(req,url, body =list(vars = multipart$vars, phens=phens, flags=toJSON(flags)))
      })
      nxt_vars = lapply(1:num_pvals, function(ik){
       
        pv =  .getPvsAll(phens,datas1[names(datas1) %in% train_nme], c(vars$var,comb[ik]),k)
        mStateObj$new(comb[ik],  .sumChisq(pv) , prev_i=vars)
      })
      names(nxt_vars) = names(comb)[1:length(nxt_vars)]  
      nxt_vars[order(unlist(lapply(nxt_vars, function(nv)nv$cumpv)))]          
    })
    angles_all1 = unlist(angles_all, rec=F)
    ord = order(unlist(lapply(angles_all1, function(nv)nv$cumpv)))
    logpv = angles_all1[[ord[1]]]$logpv
    print(logpv)
    if(logpv<=logpvthresh){
      vars_l = angles_all1[ord][1:beam]
    }
  }
  #lapply(datas1, function(d)d$saveParquet())
  lapply(vars_l, function(v) v$var)
  
  
#  lapply(urls1$urls, function(url){
#    .POST1(req,url, body =list(phens=multipart$phens, flags=multipart$flags))
#  })
}

## evaluate /<org>/<project>
#* @post /evaluate/<org>
#* @serializer tsv list(type="text/plain; charset=UTF-8")
function( req, org) {
  user=keys$validate(req,org,NULL,endpoint)
  if(is.null(user)) return(list(msg="unauthorised access"))
  multipart <- mime::parse_multipart(req)
  f=multipart$upload$datapath
  if(!is.null(f)){
    all_models=toJSON(read_json(f,simplifyVector = T))
  }else{
    all_models = multipart$models
  }
  urls1 = .readURLS(multipart$urls, self$endpoint1,"evaluate")
  lapply(urls1$urls, function(url){
    .POST1(req,url, body =list(models = all_models, phens=multipart$phens, flags=multipart$flags))
  })
}




#* Return pvalue of last variable in context of preceding /<org>/<project>
#* @post /pvalue/<org>
function( req,org) {
  user=keys$validate(req,org, NULL, endpoint)
  if(is.null(user)) return(list(msg="unauthorised access"))
  multipart <- mime::parse_multipart(req)
  urls1 = .readURLS(multipart$urls, self$endpoint1,"pvalue")
  lapply(urls1$urls, function(url){
    .POST1(req,url, body =list(vars = multipart$vars, phens=multipart$phens, flags=multipart$flags))
  })
}