
.orderFlags<-function(flags){
  flags[order(names(flags))]
}
toJSON1<-function(flags){
  if(typeof(flags)!="list") return(toJSON(sort(flags)))
  toJSON(.orderFlags(flags))
}


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
    if(is.list(mat1) && length(mat1)==1 ){
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
 data_names = "list",
 data_types = "list",
 phenos = "list",
 data_id="character",
 user="", ## default user
  initialize=function(dbDir,subnme,user="",
                      clear=FALSE){ #there is duplication in phenosdir and dbDir .. fix later
    if(!file.exists(dbDir)) dir.create(dbDir,recursive=T)
    self$subnme = subnme
    self$user=user
    self$sigsdir=paste(dbDir,subnme, sep="/")
    dir.create(self$sigsdir,recursive=T)
    self$dbfile=paste(self$sigsdir,paste("signatures",subnme,"sqlite",sep="."),sep="/")
    self$mydb=  dbConnect(RSQLite::SQLite(),self$dbfile,flags=SQLITE_RWC )
    if(clear) self$drop_all();
    self$data_id=NULL
  },
  updateFlags=function(flags, flags1,phens = NULL, transform_y=NULL,  user=self$user){
    expt_id = self$getExpt(flags=flags1, phens = phens, transform_y = transform_y, user=user,add_new=F)
    if(length(expt_id)>0){
      for(expt in expt_id){
        dbExecute(self$mydb, 'UPDATE experiment set flags=:flags where experiment_id =:expt_id',list(expt_id=expt, flags=toJSON(flags1)))
      }
    }
  },
  updateData=function(  user=self$user,data_flags = list(), data_names = list(), data_types = list(),phenos = list(), dims=list()){
    self$data_flags = data_flags
    self$data_names = data_names
    self$data_types = data_types
    self$phenos = phenos
    tbls = dbListTables(self$mydb)
    expt= data.frame(list(user=user,  flags=toJSON1(data_flags), names =toJSON1(data_names), types=toJSON1(data_types), dims = toJSON1(dims)))
    if(!("data" %in% tbls)){  
      self$data_id=0
      expt$data_id=self$data_id; expt$date=date();
      try(dbWriteTable(self$mydb, "data", expt,overwrite=T,append=F))
    }else{
      vn =  dbGetQuery(self$mydb, 'SELECT data_id from data where user=:user AND flags=:flags AND names=:names AND types=:types AND dims=:dims',expt)
      if(nrow(vn)>0){
        self$data_id = vn$data_id[[1]]
      }else{
        vn =  dbGetQuery(self$mydb, 'SELECT data_id from data')
        self$data_id = max(0,1+max(vn$data_id))
        expt$data_id=self$data_id; expt$date=date();
        try(dbWriteTable(self$mydb, "data", expt,overwrite=F,append=T))
      }
    }
 },
 saveEval=function(eval2,flags,phens,transform_y,  user=self$user,replace=T){
   expt_id = self$getExpt(flags, phens, transform_y, user,add_new=T)
   hasEval="eval" %in% self$tbls()
   if(replace & hasEval){
     dbExecute(self$mydb, 'DELETE FROM eval where experiment_id =:expt_id',list(expt_id=expt_id))
   }
   eval3 = .calcEval1(eval2, rename=F)
   eval31 = eval3 %>% tibble::add_column(experiment_id=expt_id)
 
   try(dbWriteTable(self$mydb, "eval", eval31,overwrite=!hasEval,append=hasEval))
   return(list(msg="success"))
 },
 loadEval=function(flags=NULL, phens=NULL, transform_y = NULL,  user=self$user){
   hasEval="eval" %in% self$tbls()
   if(!hasEval) return(NULL)
   expt_id = self$getExpt(flags, phens, user,transform_y = transform_y, add_new=F)
   if(is.null(expt_id)) return(NULL)
   eval1 =  dbGetQuery(self$mydb, 'SELECT * from eval where experiment_id=:exptid',list(exptid=expt_id))
   #vn[,names(vn)!="experiment_id"]
   eval1
 },
 loadModels = function(flags, phens,transform_y,   user=self$user){
   if(!("models" %in% self$tbls())) return(NULL)
   expt_id = self$getExpt(flags, phens, user,transform_y = transform_y, add_new=F)
   if(is.null(expt_id)) return(NULL)
   combined =  dbGetQuery(self$mydb, 'SELECT * from models where experiment_id=:exptid',list(exptid=expt_id))
   vn1=combined
 #  trans_y = unique(vn$transform_y); names(trans_y)=trans_y
 ##debugging
   #m_nme = vn$model_nme[[1]] ;p_g = vn$pheno_group[[1]]; rep1 = vn$rep[[1]]; ton = vn$trainedOn[[1]];
   #vn5 = subset(vn, model_nme==m_nme & pheno_group==p_g & rep==rep1 & trainedOn==ton)
   #all_models1 = lapply(trans_y, function(ty){
    # vn1 = subset(vn, transform_y ==ty)
     mod_nme = unique(vn1$model_nme); names(mod_nme)=mod_nme
    all_models1 =  lapply(mod_nme, function(m_nme){
      #print(m_nme)
       vn2 = subset(vn1, model_nme==m_nme)
       #pheno_group = unique(vn2$pheno_group); names(pheno_group)=pheno_group
       #lapply(pheno_group, function(p_g){
         vn3 = vn2 #subset(vn2, pheno_group==p_g)
         reps = vn3$rep; names(reps) = reps
         lapply(reps, function(rep1){
           vn4 = subset(vn3,rep==rep1)
           trainedOn = vn4$trainedOn; names(trainedOn)=trainedOn
           lapply(trainedOn, function(ton){
             vn5 = subset(vn4, trainedOn==ton)
            res2 = list(
              betas = fromJSONM(vn5$betas[[1]]),
            var_names = fromJSON(vn5$var_names[[1]]),
            transf = fromJSON(vn5$transf[[1]]),
            constants_proj=fromJSONM(vn5$constants_proj[[1]])
            )
            res2
           })
           
         })
     })
   #})
   list(models=all_models1, flags = flags, phens = phens, transform_y=transform_y, db=self$subnme)
 },
 saveModels=function(all_models, 
                     flags=all_models$flags, phens=all_models$phens,transform_y=all_models$transform_y,  user=self$user,replace=T){
   debug=getOption("fspls.debug",FALSE)
   all_models1 = all_models$models
   tbls = dbListTables(self$mydb)
   hasModel = "models" %in% tbls
   expt_id = self$getExpt(flags, phens, transform_y, user,add_new=T)
   if(replace && hasModel){
     dbExecute(self$mydb, 'DELETE FROM models where experiment_id =:expt_id',list(expt_id=expt_id))
   }
    combined= .merge1_new(lapply(names(all_models1), function(nme1){
        all_models2 = all_models1[[nme1]]
        .merge1_new(lapply(names(all_models2), function(nme2){
          all_models3 = all_models2[[nme2]]
          .merge1_new(lapply(names(all_models3), function(nme3){
              all_models5 = all_models3[[nme3]]
              #.merge1_new(lapply(names(all_models4), function(nme4){
              #  all_models5 = all_models4[[nme4]]
              if(debug) print(paste(nme1,nme2,nme3))
               data.frame(list(experiment_id = expt_id, nvar = length(all_models5$var_names),
                                        var_names = toJSON(all_models5$var_names),
                                        constants_proj = toJSONM(all_models5$constants_proj),
                                      transf= toJSON(all_models5$transf),
                                        betas = toJSONM(all_models5$betas), 
                  model_nme=nme1, rep=nme2,trainedOn=nme3)  ) 
              #}))
              }))
          }))
        }))
   #    }))
  try(dbWriteTable(self$mydb, "models", combined,overwrite=!hasModel,append=hasModel))
   return(list(msg="success"))
 },
 getExpt=function(flags=NULL,phens=NULL,transform_y=NULL,  user=self$user, select="experiment_id",
                  add_new =F){  #check c
   if(is.null(self$data_id)) stop("no data_id")
   expt_new = NULL; 
   if(!("experiment" %in% self$tbls()) || add_new){
       expt_new = data.frame(list(experiment_id=0,data_id = self$data_id,user=user, date=date(), flags=toJSON(flags),phens=toJSON(phens), transform_y=toJSON(transform_y)
                                    ))
       if(!("experiment" %in% self$tbls()) ){
            try(dbWriteTable(self$mydb, "experiment", expt_new[c(),,drop=F],overwrite=T,append=F))
       }
   }
   str1 = paste('SELECT',select,'from experiment')
   str = paste(str1,'where data_id=:data_id AND user=:user')
   expt = data.frame(list(data_id = self$data_id,user=user ))
   if(!is.null(flags)) { str = paste(str, "AND flags=:flags"); expt$flags=toJSON(flags) }
   if(!is.null(phens)){str = paste(str, "AND phens=:phens"); expt$phens=toJSON(phens) }
   if(!is.null(transform_y)){str = paste(str, "AND transform_y=:transform_y"); expt$transform_y=toJSON(transform_y) }
   vn =  dbGetQuery(self$mydb, str,expt)
   expt_id = NULL
    if(nrow(vn)==0 && add_new){
       vn =  dbGetQuery(self$mydb,str1)
       expt_id = max(0, vn[,1]+1)
       expt_new$experiment_id = expt_id
       if(is.null(expt_new$transform_y))   expt_new$transform_y =   toJSON(c("function(y) y","function(y) y"))
       try(dbWriteTable(self$mydb, "experiment", expt_new,overwrite=F,append=T))
       vn =  dbGetQuery(self$mydb, str,expt)
    }
  if(nrow(vn)==0) return(NULL)
   vn[[1]]
 },
  loadVars = function(flags,phens,transform_y,  user=self$user){ ##extracts variables
    expt_id = self$getExpt(flags, phens, transform_y = transform_y, user=user,add_new=F)
    if(is.null(expt_id)) return(NULL)
    vn1 =  dbGetQuery(self$mydb, 'SELECT * from vars where experiment_id=:exptid',list(exptid=expt_id))
    #ty = trans_y[[1]]
      vn1 = vn1[!duplicated(vn1$model),,drop=F]
      if(max(table(vn1$model))>1) warning("not unique")
      models = unique(vn1$model); names(models) = models
      res2 = list(
      variables = lapply(models, function(mod)  fromJSON(vn1$variables[vn1$model==mod])) ,# lapply(vn1$variables[vn1$model==mod], function(x) fromJSON(x))),
      inds = lapply(models, function(mod) fromJSON(vn1$inds[vn1$model==mod])), #lapply(vn1$inds[vn1$model==mod], function(x) fromJSON(x))),
      transf=lapply(models, function(mod) fromJSON(vn1$transf[vn1$model==mod])),
      cumpv = lapply(models, function(mod) fromJSON(vn1$cumpv[vn1$model==mod])#,lapply(vn1$cumpv[vn1$model==mod], function(x) fromJSON(x))) function(x) fromJSON(x)))
      ),
      flags=flags,
      phens =phens, 
      db=self$subnme,
      transform_y = transform_y
      )
      res2
  },
 phens=function( user=self$user, flags =NULL, transform_y = NULL){
   transform_y1  = self$getExpt(flags=NULL, phens=phens, transform_y = NULL, user=user,add_new=F, select="phens")
   lapply(transform_y1, fromJSON)
 },
 transform_y=function( user=self$user, flags =NULL, phens = NULL){
   transform_y1  = self$getExpt(flags=NULL, phens=phens, transform_y = NULL, user=user,add_new=F, select="transform_y")
   lapply(transform_y1, fromJSON)
 },
 flags=function( user=self$user, flags1 =NULL, phens = NULL, transform_y = NULL){
   flags_all1 = self$getExpt(flags=NULL, phens=phens, transform_y = transform_y, user=user,add_new=F, select="flags")
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
 clear_results=function(flags,phens, transform_y, user=self$user){
   expt_id = self$getExpt(flags=flags, phens=phens, transform_y = transform_y, user=user,add_new=F)
   tbls = self$tbls()
   if(length(expt_id)==0) return(NULL)
   if(length(expt_id)>1) stop("should be just one experiment")
   if ("vars" %in% tbls) dbExecute(self$mydb, 'DELETE FROM vars where experiment_id =:expt_id',list(expt_id=expt_id)) 
   if ("models" %in% tbls) dbExecute(self$mydb, 'DELETE FROM models where experiment_id =:expt_id',list(expt_id=expt_id))   
   if ("eval" %in% tbls) dbExecute(self$mydb, 'DELETE FROM eval where experiment_id =:expt_id',list(expt_id=expt_id))   
   if ("experiment" %in% tbls) dbExecute(self$mydb, 'DELETE FROM experiment where experiment_id =:expt_id',list(expt_id=expt_id))   
 },
 saveVars = function(vars_all,flags=vars_all$flags,phens=vars_all$phens,transform_y=vars_all$transform_y,
                       user=self$user, replace=T){
   tbls = self$tbls()
   hasVars = "vars" %in% tbls
   expt_id = self$getExpt(flags=flags, phens=phens, transform_y = transform_y, user=user,add_new=T)
   if(replace & hasVars){
     dbExecute(self$mydb, 'DELETE FROM vars where experiment_id =:expt_id',list(expt_id=expt_id))
   }
   vars_all1 = vars_all
   res1 = #.merge1_new(lapply(names(vars_all), function(t_y){
     .merge1_new(lapply(names(vars_all1$variables), function(nme1){
       vars1 = data.frame(list(experiment_id = expt_id, 
                               model = nme1,
                                             variables=toJSON(vars_all1$variables[[nme1]]), 
                               inds = toJSON(vars_all1$inds[[nme1]]), 
                               transf=toJSON(vars_all1$transf[[nme1]]),
                               cumpv = toJSON(vars_all1$cumpv[[nme1]])))
      }))
   #}))
   try(dbWriteTable(self$mydb, "vars", res1,overwrite=!hasVars,append=hasVars))
   return(list(msg="success"))
 },
  
  drop_all = function(){
    lapply(self$tbls(), function(tbl){
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
   
   try(dbWriteTable(self$mydb, "models",models[,!(names(models)%in% c("transform_y"))],overwrite=T,append=F))
   
 },
 correct1=function(){
   if(TRUE) stop("no longer used")
   vn =  dbGetQuery(self$mydb, 'SELECT * from vars')
   vn2=vn
   for(k in 1:nrow(vn)){
     vn$inds[[k]] = toJSON(fromJSON(vn$inds[[k]])$all)
     vn$cumpv[[k]] = toJSON(fromJSON(vn$cumpv[[k]])$all)
   }
   try(dbWriteTable(self$mydb, "vars", vn[,!(names(vn) %in% "transform_y"),drop=F],overwrite=T,append=F))
   
 },
 correct=function(){
   if(TRUE) stop("no longer used")
   expt1=self$experiments(all_cols=T)
   expt2 = expt1
   if(is.null(expt1[['transform_y']])){
     expt1 = expt1%>% tibble::add_column(transform_y="")
   }
   for(k in 1:nrow(expt1)){
     fl1 = fromJSON(expt1$flags[[k]])
     fp1 = fromJSON(expt1$phens[[k]])
     expt1$flags[[k]] = toJSON(.correctFlags(fl1,remove=c( "transform_y", "transform_y_inverse")))
     expt1$phens[[k]] = toJSON(fp1[[1]])
#     expt1$data_types[[k]] = fl1$data_types
     expt1$transform_y =  toJSON(c("function(y) y","function(y) y"))
       toJSON(list(transform_y = list(y = "function(y) y"), transform_y_inverse=list(y = "function(y) y")))
   }
   try(dbWriteTable(self$mydb, "experiment", expt1,overwrite=T,append=F))
 },
 
  experiments = function(user=NULL,all_cols=F, phens = NULL){
    did1 = list(data_id=self$data_id)
    if(is.null(user)){
   vn =  dbGetQuery(self$mydb, 'SELECT * from experiment where data_id=:data_id', did1)
    }else if(is.null(phens)){
      vn =  dbGetQuery(self$mydb, 'SELECT * from experiment where user=:user AND data_id=:data_id',list(user=user, data_id = self$data_id))
    }else{
      vn =  dbGetQuery(self$mydb, 'SELECT * from experiment where user=:user AND data_id=:data_id AND phens =:phens',list(user=user, data_id = self$data_id, 
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
 evals = function(flags = NULL, transform_y=NULL, phens = NULL){
   hasEval="eval" %in% self$tbls()
   if(!hasEval) return(NULL)
   expt=self$getExpt(flags=flags, transform_y = transform_y, phens = phens)
   vn =  dbGetQuery(self$mydb, 'SELECT * from eval')
   evals=subset(vn, experiment_id %in% expt)
   evals
 }
 
))

