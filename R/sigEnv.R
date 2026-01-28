##used to filter 
data_keys = fromJSON(getOption("data_keys",'["bigmatrix", "all_v_all","merge","duplicate_ordinal","genes_incls"]'))

expt_keys = fromJSON(getOption("expt_keys",'["project", "useoffset","useglmnet","pthresh","genes_incls","transform_y","batch","beam","nrep","stop_y","max","min","topn"]'))


.orderFlags<-function(flags, keys =c()){
  if(length(keys)>0) flags = flags[names(flags) %in% keys]
  if(length(flags)==0) return(flags)
  flags[order(names(flags))]
}
toJSON1<-function(flags, keys = c()){
  if(length(keys)>0) flags = flags[names(flags) %in% keys]
  if(typeof(flags)!="list") return(toJSON(sort(flags)))
  toJSON(.orderFlags(flags))
}
.compareVars<-function(vars_all, vars_all2){
  
  vt_all = .convertVarsToTable(vars_all)
  vt_all2 = .convertVarsToTable(vars_all2)
  mi1 = match(vt_all$model, vt_all2$model)
  mi2 = match(vt_all2$model, vt_all$model)
  print(mi1); print(mi2)
}



.convertVarsToTable<-function(vars_all1, expt_id=0){
  .merge1_new(lapply(names(vars_all1$variables), function(nme1){
    vars1 = data.frame(list(experiment_id = expt_id, 
                            model = nme1,
                            variables=toJSON(vars_all1$variables[[nme1]]), 
                            inds = toJSON(vars_all1$inds[[nme1]]), 
                            beam = vars_all1$beam,
                      #      transf=toJSON(vars_all1$transf[[nme1]]),
                            cumpv = toJSON(vars_all1$cumpv[[nme1]])))
  }))
}

.compareModels<-function(all_models, all_models2){
  tbl = .convertModelsToTable(all_models$models, 0)
  tbl2 = .convertModelsToTable(all_models2$models, 0)
  match(tbl$var_names, tbl2$var_names)
}
.modelToRow<-function(all_models5){
 # print(all_models5)
  varnames = lapply(all_models5$var_names,paste, collapse=".")
  if(length(varnames)>0)  names(varnames)=1:length(varnames)
 var_names = all_models5$var_names
 names(var_names) = varnames
  vn5=data.frame(list(
                     var_names = toJSON(var_names),
                     varnames = toJSON(varnames),
                     constants_proj = toJSONM(all_models5$constants_proj),
                     Wall = toJSONM(all_models5$Wall),
                     nvar = length(all_models5$var_names), 
                     mean_x = toJSON(all_models5$mean_x),
                     pvs = toJSON(all_models5$pvs),
                     pvs_all = toJSON(all_models5$pvs_all),
                     # transf= toJSON(all_models5$transf),
                     betas_proj = toJSONM(all_models5$betas_proj)
  ))
  vn5
}
.modelFromRow<-function(vn5){
 if(nrow(vn5)>1) {
   print(vn5);
   stop()
   warning("only expecting one row")
 }
  res=list(
    betas_proj = fromJSONM(vn5$betas_proj[[1]]),
    var_names =fromJSON(vn5$var_names[[1]]),
    varnames = unlist(fromJSON(vn5$varnames[[1]])),
    Wall = fromJSONM(vn5$Wall[[1]]),
    mean_x = fromJSON(vn5$mean_x[[1]]),
    pvs = fromJSON(vn5$pvs[[1]]),
    pvs_all = fromJSON(vn5$pvs_all[[1]]),
    #transf = fromJSON(vn5$transf[[1]]),
    constants_proj=fromJSONM(vn5$constants_proj[[1]])
  )
  #names(res$var_names) = res$varnames
  #names(res$varnames)=NULL
  res[unlist(lapply(res, length))>0]
}




#aa3=.splitAll(vn1, nmes,  .modelFromRow)




##converting to and from matrices.  can be either a matrix or an list of matrices
toJSONM<-function(matr){
  if(is.null(matr)) return(toJSON(c()))
  if(typeof(matr)=="list"){
    attr = lapply(matr, function(m1) attributes(m1))
    return(toJSON(list(attr=attr,m=matr)))
  }else if(typeof(matr)=="double"){
    attr =  attributes(matr)
    json=toJSON(list(attr=attr,  m = matr ), digits=NA)
    return(json)
  }else{
    stop("could not transform")
    return(toJSON(matr))
  }
}
setAttr<-function(mat1, attr1){
  if(typeof(mat1)=="list" && typeof(attr1)=="list"  && is.null(attr1$names)){
   res= lapply(1:length(attr1), function(k){
      setAttr(mat1[[k]], attr1[[k]])
    })
   names(res) = names(attr1)
   
   return(res)
  }else{
    nme = attr1$names
    dimn = attr1$dimnames
    if(!is.null(dimn) && !is.list(dimn)){
      attr1$dimnames=as.list(dimn)
    }
    if(!is.null(nme)  && length(nme)==length(mat1)){
      
    }else if(is.list(mat1) && length(mat1)==1 ){
      mat1 = mat1[[1]]
    }
    attributes(mat1) = attr1
  }
  mat1
}

fromJSONM<-function(json){
  ab1 = fromJSON(json, simplifyMatrix = T)
  attr1 = ab1$attr
  mat1 = ab1$m
  if(is.null(ab1$attr)) return(ab1)
  mat1=setAttr(mat1,attr1)
  mat1
}
.correctFlags<-function(fl1,remove=c()){
  if(!is.null(fl1[['pheno_balance']])){
    fl1[['pheno_balance']]=TRUE
  }
  fl1 = fl1[!(names(fl1) %in% remove)]
  fl1
  
}


.diffFlags<-function(flags0, flags1){
  if(toJSON(flags0)==toJSON(flags1)) return (NULL)
  nmes = c(names(flags0), names(flags1))
  nmes = nmes[!duplicated(nmes)]
  names(nmes)=nmes
  comp=lapply(nmes, function(nme){
    f1 = flags0[[nme]]
    f2 = flags1[[nme]]
    if(toJSON(f1)==toJSON(f2)) return(NULL)
    list(a=f1,b=f2)
  })
  comp[unlist(lapply(comp, length))>0]
}


sigEnv<-R6Class("sigEnv", public = list(
  mydb="S4",
  dir="character",
  dbfile="character",
   sigsdir = "character",
 subnme="character",
 data_flags = "list",
 dims = "list",
 #data_id="character",
 user="", ## default user
  initialize=function(dbDir,subnme,flags, dims, user="",
                      clear=FALSE){ #there is duplication in phenosdir and dbDir .. fix later
    if(!file.exists(dbDir)) dir.create(dbDir,recursive=T)
    self$subnme = subnme
    self$user=user
    self$sigsdir=paste(dbDir,subnme, sep="/")
    dir.create(self$sigsdir,recursive=T)
    self$dbfile=paste(self$sigsdir,paste("signatures",subnme,"sqlite",sep="."),sep="/")
    self$mydb=  dbConnect(RSQLite::SQLite(),self$dbfile,flags=SQLITE_RWC )
    if(clear) self$drop_all();
    self$data_flags =.orderFlags(flags[names(flags) %in% data_keys])
    self$dims = dims
  #  self$data_id = self$getDataID(flags,dims, add_new=T)
  },
 data_id=function(){
   self$getDataID(self$data_flags, self$dims)   
 },
  updateFlags=function(flags, flags1,phens = NULL, user=self$user){
    expt_id = self$getExpt(flags=flags1, phens = phens,  user=user,add_new=F)
    if(length(expt_id)>0){
      for(expt in expt_id){
        dbExecute(self$mydb, 'UPDATE experiment set flags=:flags where experiment_id =:expt_id',list(expt_id=expt, flags=toJSON(flags1)))
      }
    }
  },
 #updateTransforms(transform_y){
#   stop(" this not implemented until we split out the transform_y")
# },
 saveEval=function(eval2,flags,phens, user=self$user,replace=T){
   expt_id = self$getExpt(flags, phens,user,add_new=T)
   hasEval="eval" %in% self$tbls()
   if(replace & hasEval){
     dbExecute(self$mydb, 'DELETE FROM eval where experiment_id =:expt_id',list(expt_id=expt_id))
   }
   eval3 = .calcEval1(eval2, rename=F)
   eval31 = eval3 %>% tibble::add_column(experiment_id=expt_id)
 
   try(dbWriteTable(self$mydb, "eval", eval31,overwrite=!hasEval,append=hasEval))
   return(list(msg="success"))
 },
 loadEval=function(flags=NULL, phens=NULL,   user=self$user){
   hasEval="eval" %in% self$tbls()
   if(!hasEval) return(NULL)
   expt_id = self$getExpt(flags, phens, user, add_new=F)
   if(is.null(expt_id)) return(NULL)
   eval1 =  dbGetQuery(self$mydb, 'SELECT * from eval where experiment_id=:exptid',list(exptid=expt_id))
   #vn[,names(vn)!="experiment_id"]
   eval1
 },
 loadModels = function(flags, phens, user=self$user){
   if(!("models" %in% self$tbls())) return(NULL)
   expt_id = self$getExpt(flags, phens, user,add_new=F)
   if(is.null(expt_id)) return(NULL)
   combined =  dbGetQuery(self$mydb, 'SELECT * from models where experiment_id=:exptid',list(exptid=expt_id))
   vn1=combined
   #c("model_name","rep")

   all_models1 =  .splitAll(combined, c("beam","model_name","rep"),.modelFromRow)
   #})
   list(models=all_models1, flags = flags, phens = phens, trainedOn=self$subnme)
 },

 saveAngles=function(flags,phens, data_nme, comb_angs1, varnames,k){
   expt_id = self$getExpt(flags, phens)
   
   tbls = dbListTables(self$mydb)
   combined= list(experiment_id=expt_id, data=data_nme, angles = toJSONM(comb_angs1), varnames = toJSON(varnames),k=k)
   #combined1 = combined[names(combined) %in% c("experiment_id","data","varnames")]
   hasModel = "angles" %in% tbls
   if(hasModel){
    dbExecute(self$mydb, 'DELETE FROM angles where experiment_id =:expt_id',list(expt_id=expt_id))
   }
   
   #    }))
   try(dbWriteTable(self$mydb, "angles", data.frame(combined),overwrite=!hasModel,append=hasModel))
 },
 loadAngles=function(expt_id, varnames,k){
   tbls = dbListTables(self$mydb)
   r1= list(experiment_id=expt_id,  k=k,varnames = toJSON(varnames))
   #combined1 = combined[names(combined) %in% c("experiment_id","data","varnames")]
   hasModel = "angles" %in% tbls
   if(!hasModel) return (NULL)
   combined =  dbGetQuery(self$mydb, 'SELECT * from angles where experiment_id=:experiment_id and varnames=:varnames and k=:k',r1)
   inds = 1:nrow(combined); names(inds) = combined$data
   if(length(duplicated(combined$data))>1) stop("duplicates")
   lapply(inds, function(i){
     fromJSONM(combined$angles[[i]])
   })
 },
#    self$sigs$savePvals(flags, phens, data_nme, ri, varnames,k1, useCurrVarnames = F)

savePvals=function(flags,phens, data_nme, ri, varnames,k,useCurrVarnames=F){
  expt_id = self$getExpt(flags, phens,add_new=T)
  #print(paste("saving pv",expt_id, k, data_nme, toJSON(varnames), useCurrVarnames))
  varn1 = toJSON(varnames)
  tbls = self$tbls()
  combined=.merge_all(ri,c("transf","param","var") , .modelToRow)%>% tibble::add_column(prev_var=varn1, data=data_nme, experiment_id = expt_id,k=k)
  #if(getOption("verbose",F)) print(combined);
  if(useCurrVarnames) combined$prev_var= combined$var_names
  #apply(combined, 2, function(x) length(table(x))/length(x))
  ##follow to check
 # aa= .splitAll(combined, c("data","transf","param","var"),.modelFromRow)
  
  #combined= .convertModelsToTable1(ri, expt_id=expt_id,debug=F) %>% tibble::add_column(prev_var=varn1, data=data_nme)
  hasModel = "pvals" %in% tbls
  if( hasModel){
    li1 = list(expt_id=expt_id, k=k,data=data_nme, prev_var=toJSON(varnames))
    if(useCurrVarnames) li1$prev_var==unlist(lapply(combined$var_names, toJSON))
    dbExecute(self$mydb, 'DELETE FROM pvals where experiment_id =:expt_id AND data=:data AND prev_var=:prev_var AND k=:k' ,
           li1)
  }
  
  #aa = .convertTableToModels(combined)
  try(dbWriteTable(self$mydb, "pvals", combined,overwrite=!hasModel,append=hasModel))
  invisible(list(msg="success", expt_id = expt_id))
},
pvals=function(expt_id){
  combined =  dbGetQuery(self$mydb, 'SELECT * from pvals where experiment_id=:experiment_id ',list(experiment_id=expt_id))
  combined
},
loadPrev=function(expt_id, prev_i, k, data_nme = self$subnme){
  tbls = dbListTables(self$mydb)
  #combined1 = combined[names(combined) %in% c("experiment_id","data","varnames")]
  hasModel = "pvals" %in% tbls
  if(!hasModel) return(prev_i)
  r1= list(experiment_id=expt_id,  varnames = toJSON(prev_i$var_names),k=k, data = data_nme)
  combined =  dbGetQuery(self$mydb, 'SELECT * from pvals where experiment_id=:experiment_id and prev_var=:varnames and k=:k and data =:data',r1)
  if(nrow(combined)==0) return(prev_i)
  prev_i2=.modelFromRow(combined[1,])
  prev_i2
},
loadPvals=function(expt_id, varnames,k){
  if(!is.null(varnames) && is.null(names(varnames))) names(varnames) = lapply(varnames, paste, collapse=".")
  tbls = dbListTables(self$mydb)
  r1= list(experiment_id=expt_id,  varnames = toJSON(varnames),k=k)
  #combined1 = combined[names(combined) %in% c("experiment_id","data","varnames")]
  hasModel = "pvals" %in% tbls
  if(!hasModel) return (NULL)
  combined =  dbGetQuery(self$mydb, 'SELECT * from pvals where experiment_id=:experiment_id and prev_var=:varnames and k=:k',r1)
  #aa=.convertTableToModels1(combined)
  
  aa= .splitAll(combined, c("data","transf","param","var"),.modelFromRow)
aa 
},

 saveModels=function(all_models_, trainedOn=all_models_$trainedOn,
                     flags=all_models_$flags, phens=all_models_$phens,  user=self$user,replace=T){
   debug=getOption("fspls.debug",FALSE)
   tbls = dbListTables(self$mydb)
   hasModel = "models" %in% tbls
   expt_id = self$getExpt(flags, phens,  user,add_new=T)
   if(replace && hasModel){
     dbExecute(self$mydb, 'DELETE FROM models where experiment_id =:expt_id',list(expt_id=expt_id))
   }
   #transform_y=""
     all_models1 = all_models_$models
     combined=.merge_all(all_models1,c("beam","model_name","rep") , .modelToRow)%>% tibble::add_column(experiment_id = expt_id)
      try(dbWriteTable(self$mydb, "models", combined,overwrite=!hasModel,append=hasModel))
   return(list(msg="success"))
 },
#      vn =  dbGetQuery(self$mydb, 'SELECT * from data where user=:user AND flags=:flags AND names=:names AND types=:types AND dims=:dims',expt1)
# expt= data.frame(list(user=user,  flags=toJSON1(data_flags), transform_y = toJSON1(transform_y),
# names =toJSON1(data_names), types=toJSON1(data_types), dims = toJSON1(dims)))
getDataID=function(flags,dims,  add_new=F, select="data_id", user=self$user){
  flags = .orderFlags(flags, data_keys)
  expt_new = data.frame(list(data_id = 0,user=user, date=date(), 
                             flags=toJSON(flags),dims=toJSON(dims)
  ))
    if(!("data" %in% self$tbls()) || add_new  ){
     
      if(!("data" %in% self$tbls()) ||length(which(is.na(match(names(expt_new),names(self$datas()))))!=0) ){
        try(dbWriteTable(self$mydb, "data", expt_new[c(),,drop=F],overwrite=T,append=F))
      }
    }
  str1 = paste('SELECT',select,'from data')
  str = paste(str1,'where user=:user')
  expt = data.frame(list(user=user ))

  
  if(!is.null(flags)) { str = paste(str, "AND flags=:flags"); expt$flags=toJSON1(flags) }
  if(!is.null(dims)){str = paste(str, "AND dims=:dims"); expt$dims=toJSON1(dims) }
  
  vn =  dbGetQuery(self$mydb, str,expt)
  expt_id = NULL
  if(nrow(vn)==0 ){
    vn =  dbGetQuery(self$mydb,str1)
    expt_id = max(0, vn[,1]+1)
    expt_new$data_id = expt_id
    try(dbWriteTable(self$mydb, "data", expt_new,overwrite=F,append=T))
    vn =  dbGetQuery(self$mydb, str,expt)
  }
  if(nrow(vn)==0) return(NULL)
  vn[1,1]
},
 getExpt=function(flags=NULL,phens=NULL, user=self$user, select="experiment_id",
                  add_new =F){  #check c
   #transform_y = ""
   if(!is.null(flags)){
     flags = .orderFlags(flags, expt_keys)
   }
   data_id = self$data_id();
   expt_new = NULL; 
   if(!("experiment" %in% self$tbls()) || add_new){
       expt_new = data.frame(list(experiment_id=0,data_id =data_id,user=user, date=date(), 
                                  flags=toJSON(flags),phens=toJSON(phens)
                                    ))
       if(!("experiment" %in% self$tbls()) ){
            try(dbWriteTable(self$mydb, "experiment", expt_new[c(),,drop=F],overwrite=T,append=F))
       }
   }
   str1 = paste('SELECT',select,'from experiment')
   str = paste(str1,'where data_id=:data_id AND user=:user')
   expt = data.frame(list(data_id = data_id,user=user ))
   if(!is.null(flags)) { str = paste(str, "AND flags=:flags"); expt$flags=toJSON(flags) }
   if(!is.null(phens)){str = paste(str, "AND phens=:phens"); expt$phens=toJSON(phens) }
   vn =  dbGetQuery(self$mydb, str,expt)
   expt_id = NULL
    if(nrow(vn)==0 && add_new){
       vn =  dbGetQuery(self$mydb,str1)
       expt_id = max(0, vn[,1]+1)
       expt_new$experiment_id = expt_id
       try(dbWriteTable(self$mydb, "experiment", expt_new,overwrite=F,append=T))
       vn =  dbGetQuery(self$mydb, str,expt)
    }
  if(nrow(vn)==0) return(NULL)
   vn[[1]]
 },
  loadVars = function(flags,phens,user=self$user){ ##extracts variables
    expt_id = self$getExpt(flags, phens, user=user,add_new=F)
    if(!("vars" %in% self$tbls())) return (NULL)
    if(is.null(expt_id)) return(NULL)
    vn0 =  dbGetQuery(self$mydb, 'SELECT * from vars where experiment_id=:exptid',list(exptid=expt_id))
    if(nrow(vn0)==0) return(NULL)
    beams = vn0$beam; names(beams)=beams
    lapply(beams, function(bm){
      vn1 = subset(vn0, beam==bm)
      vn1 = vn1[!duplicated(vn1$model),,drop=F]
      if(max(table(vn1$model))>1) warning("not unique")
      models = unique(vn1$model); names(models) = models
      
      ## add beam
      res2 = list(
      variables = lapply(models, function(mod)  fromJSON(vn1$variables[vn1$model==mod])) ,# lapply(vn1$variables[vn1$model==mod], function(x) fromJSON(x))),
      inds = lapply(models, function(mod) fromJSON(vn1$inds[vn1$model==mod])), #lapply(vn1$inds[vn1$model==mod], function(x) fromJSON(x))),
#      transf=lapply(models, function(mod) fromJSON(vn1$transf[vn1$model==mod])),
      cumpv = lapply(models, function(mod) fromJSON(vn1$cumpv[vn1$model==mod])#,lapply(vn1$cumpv[vn1$model==mod], function(x) fromJSON(x))) function(x) fromJSON(x)))
      ),
      flags=flags,
      phens =phens, 
      beam = bm,
      db=self$subnme
      )
      res2
    })
  },
 phens=function( user=self$user, flags =NULL){
   transform_y1  = self$getExpt(flags=NULL, phens=phens,  user=user,add_new=F, select="phens")
   lapply(transform_y1, fromJSON)
 },
 
 flags=function( user=self$user, flags1 =NULL, phens = NULL){
   flags_all1 = self$getExpt(flags=NULL, phens=phens, user=user,add_new=F, select="flags")
   if(is.null(flags_all1)) return(NULL)
   flags_all=lapply(flags_all1, fromJSON)
   if(is.null(flags1)) return(flags_all)
     comp1=lapply(flags_all, function(flags0){
       .diffFlags(flags0, flags1)
     })
     o=order(unlist(lapply(comp1, function(comp){ nchar(toJSON(comp))})))
     flags2=flags_all[[o[[1]]]]
     flags2
 },
get_data_flags=function( user=self$user, nmes = self$data_names, types = self$data_types, dims = self$dims){
  query="SELECT data.*, experiment_id from data inner join experiment on experiment.data_id = data.data_id  where data.user=:user AND data.dims =:dims AND data.types=:types AND data.names =:names";
  vn =  dbGetQuery(self$mydb, query, list(user=user, dims = toJSON1(dims),names = toJSON(nmes), types = toJSON1(types)))
  tbl =table(vn$data_id)
  data_ids = names(tbl); names(data_ids) = data_ids
  vn_all = .merge1_new(lapply(data_ids, function(did){
    vn2 = subset(vn, data_id ==did )
    vn3 = vn2[1,]
    vn3[['count']]=nrow(vn2)
    vn3[names(vn3) %in% c("flags","count")]
  }))
   vn_all[order(vn_all$count, decreasing=T),]
},
 clear_results=function(flags,phens,  user=self$user){
   expt_id = self$getExpt(flags=flags, phens=phens,  user=user,add_new=F)
   tbls = self$tbls()
   if(length(expt_id)==0) return(NULL)
   if(length(expt_id)>1) stop("should be just one experiment")
   if ("vars" %in% tbls) dbExecute(self$mydb, 'DELETE FROM vars where experiment_id =:expt_id',list(expt_id=expt_id)) 
   if ("models" %in% tbls) dbExecute(self$mydb, 'DELETE FROM models where experiment_id =:expt_id',list(expt_id=expt_id))   
   if ("eval" %in% tbls) dbExecute(self$mydb, 'DELETE FROM eval where experiment_id =:expt_id',list(expt_id=expt_id))   
   if ("experiment" %in% tbls) dbExecute(self$mydb, 'DELETE FROM experiment where experiment_id =:expt_id',list(expt_id=expt_id))   
 },
 saveVars = function(vars_all,
                       user=self$user, replace=T){
   tbls = self$tbls()
   hasVars = "vars" %in% tbls
   flags=vars_all[[1]]$flags;
   phens=vars_all[[1]]$phens
   expt_id = self$getExpt(flags=flags, phens=phens,  user=user,add_new=T)
   if(replace & hasVars){
     dbExecute(self$mydb, 'DELETE FROM vars where experiment_id =:expt_id',list(expt_id=expt_id))
   }
   for(vars_all1 in vars_all){
   res1 = .convertVarsToTable(vars_all1, expt_id=expt_id)
   
   #}))
   try(dbWriteTable(self$mydb, "vars", res1,overwrite=!hasVars,append=hasVars))
   hasVars=T
   }
   return(list(msg="success"))
 },
  
  drop_all = function(exclude = c()){
    tbls = self$tbls()
    tbls = tbls[!(tbls %in% exclude)]
    lapply(tbls, function(tbl){
      str = paste("DROP table ",tbl)
      try(dbExecute(self$mydb, str))
    })
  },
 
 clear_all_user = function(user,vars = T, models=T, evals=T, experiment =T, data=T){
   tbls = self$tbls();
  expts = self$experiments(user)   
  if(length(expts)==0) return(list(msg="none"))
  exp = list(expt_id=expts)
  v = list(user=user)
   if(experiment && "experiment" %in% tbls)try(dbExecute(self$mydb, 'DELETE  from experiment where user=:user', v))
   if(data && "data" %in% tbls) try(dbExecute(self$mydb, 'DELETE FROM data where user=:user',v))
   if(vars && "vars" %in% tbls) try(dbExecute(self$mydb, 'DELETE FROM vars where experiment_id=:expt_id',exp))
   if(models && "models" %in% tbls) try(dbExecute(self$mydb, 'DELETE FROM models where experiment_id=:expt_id',exp))
   if(evals && "eval" %in% self$tbls()) try(dbExecute(self$mydb, 'DELETE FROM eval where experiment_id=:expt_id',exp))
  list(msg="success")
 },
 tbls = function(){
   return(dbListTables(self$mydb))
 },
 correct2=function(){
   if(TRUE) stop("no longer used")
   models =  dbGetQuery(self$mydb, 'SELECT * from models')
   models2=models
   
   try(dbWriteTable(self$mydb, "models",models,overwrite=T,append=F))
   
 },
 correct1=function(){
   if(TRUE) stop("no longer used")
   vn =  dbGetQuery(self$mydb, 'SELECT * from vars')
   vn2=vn
   for(k in 1:nrow(vn)){
     vn$inds[[k]] = toJSON(fromJSON(vn$inds[[k]])$all)
     vn$cumpv[[k]] = toJSON(fromJSON(vn$cumpv[[k]])$all)
   }
   try(dbWriteTable(self$mydb, "vars", vn,overwrite=T,append=F))
   
 },
  experiments = function(user=NULL,all_cols=F, phens = NULL){
    data_id = self$data_id()
    did1 = list(data_id=data_id)
    if(is.null(user)){
   vn =  dbGetQuery(self$mydb, 'SELECT * from experiment where data_id=:data_id', did1)
    }else if(is.null(phens)){
      vn =  dbGetQuery(self$mydb, 'SELECT * from experiment where user=:user AND data_id=:data_id',list(user=user, data_id = data_id))
    }else{
      vn =  dbGetQuery(self$mydb, 'SELECT * from experiment where user=:user AND data_id=:data_id AND phens =:phens',list(user=user, data_id = data_id, 
                                                                                                                          phens=toJSON(phens) ))
    }
    if(all_cols) return(vn)
   vn$experiment_id
 },
 datas = function(user=NULL){
   if(is.null(user)){
     vn =  dbGetQuery(self$mydb, 'SELECT * from data')
   }else{
     vn =  dbGetQuery(self$mydb, 'SELECT * from experiment where user=:user',list(user=user))
   }
   vn
 },
 models = function(){
     vn =  dbGetQuery(self$mydb, 'SELECT * from models')
   vn
 },
 vars = function(){
   vn =  dbGetQuery(self$mydb, 'SELECT * from vars')
   vn
 },
 evals = function(flags = NULL, phens = NULL){
   hasEval="eval" %in% self$tbls()
   if(!hasEval) return(NULL)
   expt=self$getExpt(flags=flags, phens = phens)
   vn =  dbGetQuery(self$mydb, 'SELECT * from eval')
   evals=subset(vn, experiment_id %in% expt)
   evals
 }
 
))

