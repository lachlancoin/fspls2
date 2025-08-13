
.orderFlags<-function(flags){
  flags[order(names(flags))]
}
toJSON1<-function(flags){
  if(typeof(flags)!="list") return(toJSON(sort(flags)))
  toJSON(.orderFlags(flags))
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
  initialize=function(dbDir,subnme,
                      clear=FALSE){ #there is duplication in phenosdir and dbDir .. fix later
    if(!file.exists(dbDir)) dir.create(dbDir,recursive=T)
    self$subnme = subnme
    self$sigsdir=paste(dbDir,subnme, sep="/")
    dir.create(self$sigsdir,recursive=T)
    self$dbfile=paste(self$sigsdir,paste("signatures",subnme,"sqlite",sep="."),sep="/")
    self$mydb=  dbConnect(RSQLite::SQLite(),self$dbfile,flags=SQLITE_RWC )
    if(clear) self$drop_all();
    self$data_id=NULL
  },
  updateData=function( user="",data_flags = list(), data_names = list(), data_types = list(),phenos = list(), dims=list()){
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
 saveEval=function(eval2,flags,phens,user="",replace=T){
   expt_id = self$getExpt(flags, phens, user,add_new=T)
   hasEval="eval" %in% self$tbls()
   if(replace & hasEval){
     dbExecute(self$mydb, 'DELETE FROM eval where experiment_id =:expt_id',list(expt_id=expt_id))
   }
   eval3 = .calcEval1(eval2)
   eval31 = eval3 %>% tibble::add_column(experiment_id=expt_id)
 
   try(dbWriteTable(self$mydb, "eval", eval31,overwrite=!hasEval,append=hasEval))
   return(list(msg="success"))
 },
 loadEval=function(flags, phens, user=""){
   hasEval="eval" %in% self$tbls()
   if(!hasEval) return(NULL)
   expt_id = self$getExpt(flags, phens, user,add_new=F)
   if(is.null(expt_id)) return(NULL)
   vn =  dbGetQuery(self$mydb, 'SELECT * from eval where experiment_id=:exptid',list(exptid=expt_id))
   vn[,names(vn)!="experiment_id"]
 },
 loadModels = function(flags, phens, user=""){
   if(!("models" %in% self$tbls())) return(NULL)
   expt_id = self$getExpt(flags, phens, user,add_new=F)
   if(is.null(expt_id)) return(NULL)
   vn =  dbGetQuery(self$mydb, 'SELECT * from models where experiment_id=:exptid',list(exptid=expt_id))
   trans_y = unique(vn$transform_y); names(trans_y)=trans_y
 ##debugging
   #t_y = trans_y[[1]]; m_nme = vn$model_nme[[1]] ;p_g = vn$pheno_group[[1]]; rep1 = vn$rep[[1]]; ton = vn$trainedOn[[1]];
   #vn5 = subset(vn, model_nme==m_nme & pheno_group==p_g & rep==rep1 & trainedOn==ton)
   all_models1 = lapply(trans_y, function(ty){
     vn1 = subset(vn, transform_y ==ty)
     mod_nme = unique(vn1$model_nme); names(mod_nme)=mod_nme
     lapply(mod_nme, function(m_nme){
       vn2 = subset(vn1, model_nme==m_nme)
       pheno_group = unique(vn2$pheno_group); names(pheno_group)=pheno_group
       lapply(pheno_group, function(p_g){
         vn3 = subset(vn2, pheno_group==p_g)
         reps = vn3$rep; names(reps) = reps
         lapply(reps, function(rep1){
           vn4 = subset(vn3,rep==rep1)
           trainedOn = vn4$trainedOn; names(trainedOn)=trainedOn
           lapply(trainedOn, function(ton){
             vn5 = subset(vn4, trainedOn==ton)
            res2 = list(betas = fromJSONM(vn5$betas[[1]]),
            var_names = fromJSON(vn5$var_names[[1]]),
            constants_proj=fromJSONM(vn5$constants_proj[[1]]))
            res2
           })
           
         })
       })
     })
   })
   all_models1
 },
 saveModels=function(all_models, flags, phens, user="",replace=T){
 
   tbls = dbListTables(self$mydb)
   hasModel = "models" %in% tbls
   expt_id = self$getExpt(flags, phens, user,add_new=T)
   if(replace && hasModel){
     dbExecute(self$mydb, 'DELETE FROM models where experiment_id =:expt_id',list(expt_id=expt_id))
   }
  # t_y = names(all_models)[[1]]; nme1 = names(all_models[[1]]); nme2 = names(all_models[[1]][[1]]); nme3 = names(all_models[[1]][[1]][[1]]); nme4 = names(all_models[[1]][[1]][[1]][[1]])
  # all_models5 = all_models[[1]][[1]][[1]][[1]][[1]]
  # all_models4 = all_models[[1]][[1]][[1]][[1]]
  combined=.merge1_new(lapply(names(all_models), function(t_y){
     all_models1 = all_models[[t_y]]
     .merge1_new(lapply(names(all_models1), function(nme1){
        all_models2 = all_models1[[nme1]]
        .merge1_new(lapply(names(all_models2), function(nme2){
          all_models3 = all_models2[[nme2]]
          .merge1_new(lapply(names(all_models3), function(nme3){
              all_models4 = all_models3[[nme3]]
              .merge1_new(lapply(names(all_models4), function(nme4){
                all_models5 = all_models4[[nme4]]
               data.frame(list(transform_y = t_y,experiment_id = expt_id, nvar = length(all_models5$var_names),
                                        var_names = toJSON(all_models5$var_names),
                                        constants_proj = toJSONM(all_models5$constants_proj),
                                        betas = toJSONM(all_models5$betas), 
                  model_nme=nme1, pheno_group=nme2, rep=nme3,trainedOn=nme4)  ) 
              }))
              }))
          }))
        }))
       }))
  try(dbWriteTable(self$mydb, "models", combined,overwrite=!hasModel,append=hasModel))
   return(list(msg="success"))
 },
 getExpt=function(flags,phens, user, add_new =F){  #check c
   if(is.null(self$data_id)) stop("no data_id")
   if(!("experiment" %in% self$tbls())){
       experiment = data.frame(list(experiment_id=0,data_id = 0,user="", date=date(), flags="",phens="" ))
       try(dbWriteTable(self$mydb, "experiment", experiment[c(),,drop=F],overwrite=T,append=F))
   }
   data_id = self$data_id
   expt = data.frame(list(data_id = self$data_id,user=user, phens = toJSON1(phens), flags=toJSON1(flags) ))
   vn =  dbGetQuery(self$mydb, 'SELECT experiment_id from experiment where data_id=:data_id AND user=:user AND flags=:flags AND phens=:phens',expt)
    expt_id = NULL
    if(nrow(vn)==0){
     if(add_new){
       vn =  dbGetQuery(self$mydb, 'SELECT experiment_id from experiment')
       expt_id = max(0, vn$experiment_id+1) 
       expt$experiment_id = expt_id
       expt$date = date()
       try(dbWriteTable(self$mydb, "experiment", expt,overwrite=F,append=T))
     }
   }else{
     expt_id = vn$experiment_id[1]
   }
   expt_id
 },
  loadVars = function(flags,phens,user=""){ ##extracts variables
    expt_id = self$getExpt(flags, phens, user,add_new=F)
    if(is.null(expt_id)) return(NULL)
    vn =  dbGetQuery(self$mydb, 'SELECT * from vars where experiment_id=:exptid',list(exptid=expt_id))
    trans_y = unique(vn$transform_y); names(trans_y)=trans_y
    #ty = trans_y[[1]]
    vars_all1 = lapply(trans_y, function(ty){
      vn1 = subset(vn, transform_y ==ty)
      vn1 = vn1[!duplicated(vn1$model),,drop=F]
      if(max(table(vn1$model))>1) warning("not unique")
      models = unique(vn1$model); names(models) = models
      res2 = list(
      variables = lapply(models, function(mod)  fromJSON(vn1$variables[vn1$model==mod])) ,# lapply(vn1$variables[vn1$model==mod], function(x) fromJSON(x))),
      inds = lapply(models, function(mod) fromJSON(vn1$inds[vn1$model==mod])), #lapply(vn1$inds[vn1$model==mod], function(x) fromJSON(x))),
      cumpv = lapply(models, function(mod) fromJSON(vn1$cumpv[vn1$model==mod])#,lapply(vn1$cumpv[vn1$model==mod], function(x) fromJSON(x))) function(x) fromJSON(x)))
      )
      )
      res2
    })
    vars_all1
  },
 saveVars = function(vars_all,flags,phens,
                      user="", replace=T){
   tbls = self$tbls()
   hasVars = "vars" %in% tbls
   expt_id = self$getExpt(flags, phens, user,add_new=T)
   if(replace & hasVars){
     dbExecute(self$mydb, 'DELETE FROM vars where experiment_id =:expt_id',list(expt_id=expt_id))
   }
   res1 = .merge1_new(lapply(names(vars_all), function(t_y){
     vars_all1 = vars_all[[t_y]]
     .merge1_new(lapply(names(vars_all1$variables), function(nme1){
       vars1 = data.frame(list(transform_y = t_y,experiment_id = expt_id, 
                               model = nme1,
                                             variables=toJSON(vars_all1$variables[[nme1]]), 
                               inds = toJSON(vars_all1$inds[[nme1]]), 
                               cumpv = toJSON(vars_all1$cumpv[[nme1]])))
      }))
   }))
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
 
  experiments = function(user=NULL,all_cols=F){
    if(is.null(user)){
   vn =  dbGetQuery(self$mydb, 'SELECT * from experiment')
    }else{
      vn =  dbGetQuery(self$mydb, 'SELECT * from experiment where user=:user',list(user=user))
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
 evals = function(){
   hasEval="eval" %in% self$tbls()
   if(!hasEval) return(NULL)
   vn =  dbGetQuery(self$mydb, 'SELECT * from eval')
   vn
 }
 
))

