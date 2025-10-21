## this is a class which holds a data object and interacts with the coordination node
.getAllSparseMatrices<-function(data, hasNA=T, convertToBigMatrix=F){
  rn = unlist(lapply(data, function(d1) rownames(d1)))
  rn = rn[!duplicated(rn)]
  lapply(data, function(mat){
    .getSparseMatrices(mat, hasNA=hasNA, convertToBigMatrix = convertToBigMatrix,rn = rn)
  })
  
}
##this function removes NAs
## if no NA matrixNA is just empty matrix
.getSparseMatrices<-function(mat, hasNA=T, convertToBigMatrix=F,rn = rownames(mat)){
  mi1 =  match(rownames(mat), rn)
  mi0 =  match(rn,rownames(mat))
  newNA=T
  if(length(mi1)==length(rn)){
    if(max(abs(apply(cbind(mi1, 1:length(rn)),1,diff)))==0) newNA=F
  }
  newNA = length(which(is.na(mi0))>0)
  
  if(!hasNA& !newNA){
    if(convertToBigMatrx){
      m2=matrix(0, nrow = nrow(mat), ncol = ncol(mat))
      res1 = list(matrix = as.big.matrix(mat),
                  matrixNA = as.big.matrix(m2)
      )
    }else{
      res1 = list(matrix = mat,
                  matrixNA = Matrix(0,nrow(mat) , ncol(mat), sparse = T)
      )
    }
  }else{
    if(newNA){
      m1=apply(mat,2,function(v){
        v1 = v[mi0]
        mv = mean(v, na.rm=T)
        v1[is.na(v1)]=mv
        v1
      })
      rownames(m1) = rn
      m2 = apply(mat,2,function(v){
        v1 = rep(1, length(rn))
        v1[mi1[!is.na(v)]]=0
        v1
      })
      rownames(m2) = rn
      
    }else{
      m1=apply(mat,2,function(v){
        mv = mean(v, na.rm=T)
        v[is.na(v)]=mv
        v
      })
      m2 = apply(mat,2,function(v){
        v1 = rep(0, length(v))
        v1[is.na(v)]=1
        v1
      })
    }
    if(convertToBigMatrix){
      res1 = list(matrix= as.big.matrix(m1), matrixNA = as.big.matrix(m2))
    }else{
      res1 = list( matrix = Matrix(m1, sparse=T),matrixNA = Matrix(m2, sparse=TRUE)) 
    }
  }
  res1
}

getFullModels<-function(all_models){
  .nonZero(lapply(all_models, function(all_model1){
    .nonZero(lapply(all_model1, function(all_model2){
      .nonZero(lapply(all_model2, function(all_model3){
        all_model3[names(all_model3) %in% "full"]
      }))
    }))
  }))
}


.getFamily<-function(y_mat, family1=NULL, max_ordinal=20){
  types = attr(y_mat, "types")
  
  if(!is.null(types)){
    indst = 1:length(types)
    names(indst) = names(y_mat)
    family = (lapply(indst, function(i){
      typ = types[[i]]
      if(typ=="double") return("gaussian")
      if(typ=="boolean") return("binomial")
      if(typ=="integer") {
        tbl=if(is.list(y_mat)) table(y_mat[[i]]) else table(y_mat[,i])
        if(length(tbl)<=2) return("binomial")
        if(length(tbl)>max_ordinal) return("gaussian")
        return("ordinal")
      }
      if(typ=="character"){
        tbl=if(is.list(y_mat)) table(y_mat[[i]]) else table(y_mat[,i])
        if(length(tbl)<=2) return("binomial")
        return("multinomial")
      }
    }))
    subinds = unlist(lapply(family, length))==0
    return(family)
  }
  if(typeof(y_mat)=="list"){
    nmey = names(y_mat); names(nmey)=nmey
    family = lapply(nmey,function(ynme){
      y  = y_mat[[ynme]]
      if(!is.numeric(y)){
        y = as.factor(y);
      }
      if(is.factor(y)){
        return(if(length(levels(y))<=2 )"binomial" else "multinomial")
      }
      if(length(unique(y[!is.na(y)]))<=2) return("binomial")
      vals = unique(y)
      if(sum(abs(vals-round(vals)), na.rm=T)<1e-9) {
        if(length(table(vals))>max_ordinal) return("gaussian")
        return("ordinal")
      }
      return("gaussian")
    })
    return(family)
  }else{
    
    famsy = (apply(y_mat,2,function(y){
      if(!is.numeric(y)){
        y = as.factor(y);
      }
      if(is.factor(y)){
        return(if(length(levels(y))<=2 )"binomial" else "multinomial")
      }
      if(length(unique(y[!is.na(y)]))<=2) return("binomial")
      vals = unique(y)
      if(sum(abs(vals-round(vals)), na.rm=T)<1e-9){
        if(length(vals)>max_ordinal) return("gaussian")
        return("ordinal")
      }
      return("gaussian")
    }))
    famsy
  }
}





dataH<-R6Class("dataH", public = list(
  data ="environment",
  sigsdir="character",
  sigs="environment",
 type="character",  
 flags="list",
 transform_y="character",
 nme="character",
  initialize=function(
    d,
    nme = "none",
    y = d$y,
    flags = list(),
    convertToBigMatrix=F,
    hasNA=T,
    mat = .getAllSparseMatrices(d$data,hasNA=hasNA, convertToBigMatrix=convertToBigMatrix),
    family= .getFamily(d$y),
    dbDir="./",
    memDir=NULL){
    self$flags = flags
    transform_y =.readFlag(flags, "transform_y",toJSON(list(x=list(unvfunc="function(y,param) y",func="function(y,param) y", param=1))))
    ### MAKE SIGNATURE DIRECTORY
    self$nme=nme
    self$sigsdir=paste(dbDir,paste0("fspls_signatures__",nme,sep="/"))
    self$sigs = list()
    #####
    self$data = 
      dataObj$new(mat, nme,dbDir,flags,  transform_y = transform_y,
                  incl_full=T,seed = getOption("seed",42), memDir=if(is.null(memDir)) NULL else paste(memDir, nme,sep="/"))
    self$type="slow1"
    types_all = getOption("types_all",names(self$data$data))
    names(types_all) = types_all
    batch=.readFlag(flags, "batch",0)
    all_v_all = .readFlag(flags,"all_v_all",F)
    one_v_rest = .readFlag(flags,"one_v_rest",F)
    nrep = .readFlag(flags,"nrep",if(batch>0) 0 else 1)
    pheno_balance=.readFlag(flags,"pheno_balance",NULL)
    varn = getOption("varn",c())
    #invisible(lapply(1:length(datas), function(ik) {
      y1 = d$y #ys[[ik]] #dists[[ik]]$updateYdb(cats[['cats']])
     # family = families[[ik]]
      
      ##need to work on all_v_all
      self$data$updateY(y1, family=family, CHECK=T, all_v_all=all_v_all, one_v_rest = one_v_rest)
      self$sigsdir=paste(dbDir,"fspls_signatures1",sep="/")
      dir.create(self$sigsdir, recursive=F, showWarnings=F)
      self$sigs=   sigEnv$new(self$sigsdir,nme, clear=F)
      # dims1 = self$dims()
      transform_y1 = self$sigs$updateData(data_flags = self$flags, 
                                   data_names =nme, 
                                   data_types = names(self$data$data),
                                   dims = self$data$dims(),
                                   transform_y = fromJSON(transform_y)
      )
      transform_y = toJSON1(transform_y1)
      self$data$transforms = .convertToTransform(transform_y)
      self$transform_y = transform_y1
  },
 updateTransform=function(transform_y){
   stop("not working yet, need to consider whats happening in db")
 },
 clear_db=function(drop=F, exclude="vars"){
   if(drop){
     self$sigs$drop_all(exclude=exclude)
   }else{
     warning("need to set drop=T if you are sure, this will delete all saved signatures")
   }
   
 },
 pheno=function(maxpheno=1e9,sep=F, sep_group = F, exclude=NULL, code=NULL, memb=NULL){
   res = self$data$pheno(maxpheno=maxpheno, sep=sep, sep_group = sep_group, code = code,memb=memb);
   if(!is.null(exclude)){
     lapply(res, function(res2){
       lapply(res2, function(res1) res1[-grep(exclude,res1)])
     })
   }
   
   res
 },
 dims=function(){
   self$data$dims()
 },
 nreps=function(){
   ncol(self$data$looc$incl)
 },
 getAngles=function(varnames){
   varnames = vars_l1[[1]]$var_names; 
   type = self$type
     #print(data_nme); 
     angleH=list(angles=
                   self$data$getAngles1(phens,varnames,incl=incl,k=k, type=type),
                 cols_incl = self$datasH[[data_nme]]$data$cols_incl(var_thresh[[data_nme]],incl, g_incl,qq, excl=varnames)) ### fix 
     .combineAngles1(angleH, incl, topn=topn, onlyAll=onlyAll, excl=varnames)
 },
 cats = function(maxpheno = 1e9){
   self$data$cats(maxpheno)
 },
 getPvsAll=function(subphens, prev_i, b_i_name,k, #   prev_i = vars_l1[[nmed]]
                      Wall =lapply(subphens, function(f) matrix(nrow=0,ncol=0)),
                      useglm=F ,inv_transform=getOption("x_transform",T),
                      project=T, useoffset=T){
     d = self$data
     family = strsplit(names(subphens)[[1]],"\\.")[[1]][1]
     if(family=="multinomial") useoffset=F
     prev_i1 = d$makeNextModel(prev_i,b_i_name,subphens,k, Wall,family, ypred=NULL, 
                               inv_transform=inv_transform,
                               project=project, useglm=useglm, logpthresh =0, useoffset=useoffset)
     prev_i1
#   pvs_all = pvs_all[unlist(lapply(pvs_all, length))>0]
  # pvs_all
 },
plotData=function(vars_all1, phens1 = vars_all1$phens, all_types=F, transform_x = NULL, violin=F, assoc=F){
  df4= #.merge1_new( 
   # lapply(self$datas, function(d) 
      self$data$plotData(vars_all1, phens1 = phens1, all_types=all_types, transform_x = transform_x, violin=violin, assoc=assoc)
               #     addName="dataset")
  df4 = df4%>%tibble::add_column(dataset=self$nme);
  facet= if(!is.null(df4[['transform']]) ) "transform~pheno1" else "pheno1"
  df4$y = factor(df4$y)
  gene_levs = unlist(lapply(unlist(vars_all1$variables,rec=F), function(x) x[[2]]))
  gene_levs = gene_levs[!duplicated(gene_levs)]
  df4$gene = factor(df4$gene, levels = gene_levs)
  df5=df4%>% separate(col="pheno",sep="\\.", into=c("family","pheno1"))
  df6=subset(df5, family=="gaussian")
  df7=(subset(df5, family!="gaussian"))
  # color= if(length(unique(df4$dataset))>1) "dataset" else #"pheno"
  ggp = NULL; ggp1 = NULL
  if(nrow(df7)>0){
    ggp<-ggplot(df7, aes(x=y, y=value, color=gene, shape=data, linetype=dataset))+facet_wrap(facet, scales="free");#+ggtitle(unlist(phens1))
    if(violin){
      ggp<-ggp+geom_violin()+geom_point()
    }else{
      ggp<-ggp+geom_boxplot()
    }
  }
  if(nrow(df6)>0){
    prbs=c(0.33,0.5,0.66)
    df6_1 =df6 %>% unite("comb",gene,data,family,pheno1,dataset,sep="__")
    comb1=unique(df6_1$comb); names(comb1)=comb1
    quants=.merge1_new(lapply(comb1,function(c1){
      df6_2 = subset(df6_1, comb==c1)
      df6_2$y = as.numeric(as.character(df6_2$y))
      df6_2 = df6_2[order(df6_2$y),]
      yv = unique(df6_2$y)
      names(yv)=yv
      .merge1_new(lapply(yv, function(yv1){
        df6_3 = subset(df6_2, y<=yv1)
        q1=quantile(df6_3$value,na.rm=T,probs = prbs)
        df_=data.frame(t(data.frame(q1)))
        df_
      }),addName="y")
    }), addName="comb")
    quants1 = quants %>% separate("comb", sep="__", into=c("gene","data","family","pheno1","dataset"))
    names(quants1) = gsub("\\.","",names(quants1))
    quants1$y = as.numeric(quants1$y)
    ggp1<-ggplot(quants1, aes(x=y, y=X50, ymin=X33, ymax = X66,color=gene, fill=gene,shape=data, linetype=dataset))+facet_wrap(facet, scales="free");#+ggtitle(unlist(phens1))
    ggp1<-ggp1+geom_line()+geom_ribbon(alpha=0.1)
  }
  
  list("binomial"=ggp,"gaussian"=ggp1)
},
updateWeights=function(subphens = self$pheno()[[1]][1]){ ## upweights low count values
  for(k in 1:length(self$datas)){
    self$datas[[k]]$updateWeights(subphens)
  } 
},
getProjectedData=function(varnames){
#  lapply(self$datas, function(d){
  d = self$data
    d$getProjectedData(varnames = varnames);      
 # })
},
getVariance=function(varnames){
  d = self$data
    d$getVariance();      
},
updateTrain=function( phens, flags,  verbose=F){
    self$data$updateTrain( phens,flags,verbose=verbose)
},
updateLOOC=function( phens, flags,varn=c(),force=F, verbose=F){
  self$data$updateLOOC( phens,flags,varn=varn,force=force, verbose=verbose); ### update training object - updates all
},
makeModels=function(vars2, inds, phens,flags){
  d = self$data
  logpthresh= log(.readFlag(flags,"pthresh",1e-3))
  project=.readFlag(flags,"project",TRUE)
  useoffset=.readFlag(flags,"useoffset",TRUE)
#  train_nme = .readFlag(flags,'train', names(datas)[1])
  if(length(which(train_nme %in% names(self$datas)))==0)train_nme = names(self$datas)[[1]]
  verbose=.readFlag(flags,"verbose",FALSE)
  if(!is.null(flags[['useglm']])) stop("define useglmnet not useglm")
  useglm=.readFlag(flags,"useglmnet",TRUE)
  inds1 = inds#[[nmes_inds1]]
  phens1 = phens#[[nmes_inds1]]
  #k=inds1[[1]]; d = datas[[1]]
  mods1 = lapply(inds1, function(k){
    # print(k)
   # lapply(datas[names(datas) %in% train_nme], function(d){
    l3=list(
      mods = d$makeModels(phens1, vars2,k,logpthresh = logpthresh,project=project,
                          flags=flags,checkRMSV=F,
                          useglm=useglm, useoffset=useoffset))
    names(l3) =self$nme
      l3
    #})
    #})
  })
  if(length(mods1)==0) stop("!!")
  models=mods1
  #})
  vars = names(models[[1]][[1]])
  names(vars) = vars
  models2 = lapply(vars, function(v){
    # lapply(models, function(models1){
    m3 = lapply(models, function(m){
      m2 = lapply(m, function(m1) m1[[v]])
      m2[unlist(lapply(m2, length))>0]
    })
    m3[unlist(lapply(m3, length))>0]
    #})
  })
  #  print(names(models2))
  models2
},
makeAllModels=function(vars_all, phens=vars_all$phens, flags=vars_all$flags, verbose=F, max = 1e6,
                       user="",useDB=T,
                       db=vars_all$db){
  sigDB = if(useDB) self$sigs else NULL
  if(!is.null(sigDB) ){
    all_models =try( sigDB$loadModels(flags,phens))
    if(inherits(all_models,"try-error")) {
      print(paste("problem reading from DB .. recalculating"))
    }else if(!is.null(all_models) && length(all_models$models)>0){
      return(all_models)
    }
  }
  self$updateLOOC(phens,flags)
  logpthresh= log(.readFlag(flags,"pthresh",1e-3))
  project=.readFlag(flags,"project",TRUE)
  vars = vars_all#[[nme_v_all]]
  all_models = list()
  variables = vars$variables
  var_inds = vars$inds
  rem_inds = 1:self$nreps()
  names(rem_inds) = as.character(rem_inds)
  names(rem_inds)[which(rem_inds==self$nreps())]="full"
  all_models = self$makeModels(list(),rem_inds , phens, flags)
  
  if(length(variables)==0) return(list(models=all_models, flags = flags, phens = phens, db=db))
  ord = order(unlist(lapply(variables, length)),decreasing=T)
  variables = variables[ord]
  var_inds = var_inds[ord]
#   v_nme = names(vars_all$variables)[1]; max=10; verbose=T; k=1;variables =vars_all$variables; 
  
  for(v_nme in names(variables)){
    # print(v_nme)
    if(verbose)print(v_nme)
    vars2 = variables[[v_nme]]
    #     var_transf=strsplit(names(transf[[v_nme]])[[1]],"_")[[1]]
    #    if(verbose)print(var_transf)
    vars2 = vars2[1:min(length(vars2), max)]
    inds =var_inds[[v_nme]]
    nme_ = paste(names(vars2),collapse=";")
    models1 = all_models[[nme_]]
    if(is.null(models1)){
      models1 = self$makeModels( vars2, inds,phens,flags)
      for(k in 1:length(models1)){
        mod1 =   all_models[[names(models1)[[k]]]]
        if(is.null(mod1)){
          all_models[[names(models1)[[k]]]] = models1[[k]]
        }else{
          # for(p_nme in names(models1[[k]])){
          mod2 = all_models[[names(models1)[[k]]]]#[[p_nme]]
          if(is.null(mod2)){
            all_models[[names(models1)[[k]]]]= models1[[k]]#[[p_nme]]#[[p_nme]] 
            
          }else{
            for(r_nme in names(models1[[k]])){ #[[p_nme]])){
              mod3 = all_models[[names(models1)[[k]]]][[r_nme]] #[[p_nme]][[r_nme]]
              if(is.null(mod3)){
                all_models[[names(models1)[[k]]]][[r_nme]]= models1[[k]][[r_nme]] #[[p_nme]][[r_nme]] #[[p_nme]][[r_nme]] 
                
              }else{
                print("already calculated!")
              }
            }
          }
          #}
        }
      }
    }else{ ## fill in gaps
      #for(p_i in 1:length(inds)){
      # p_nme = names(inds)[p_i]
      models2 = all_models[[nme_]]#[[p_nme]]  
      if(is.null(models2)){
        models2 = self$makeModels( vars2, inds#[p_i]
                                   ,phens#[which(names(phens) %in% p_nme)]
                                   ,flags)
        for(k in 1:length(models2)){
          all_models[[names(models2)[[k]]]] = models2[[k]]#[[p_nme]]#[[p_nme]]
          
        }
      }else{
        subinds= inds[which(is.na(match(names(inds), names(all_models[[nme_]]))))]
        if(length(subinds)>0){ ## missing inds
          models3 = self$makeModels( vars2, subinds,phens#[which(names(phens) %in% p_nme)]
                                     ,flags)
          for(nme1_ in names(models3)){
            for(r_i in names(models3[[nme1_]])){
              all_models[[nme1_]][[r_i]] = models3[[nme1_]][[r_i]]
            }
          }
        }
      }
      #}
      
    }
  }
  all_models_=list(models=all_models, flags = flags, phens = phens, db=db)
  #})
  if(!is.null(sigDB) ){
    sigDB$saveModels (all_models_)
  }
  #combined_models
  all_models_
},
evaluateAllModels=function(all_modelsh, phens=all_modelsh$phens,flags=all_modelsh$flags,verbose=F,useDB=T,
                           db=all_modelsh$db, user="", inv_transform_y=!.readFlag(flags,"x_transform",T)){ ## different folds with same variables
  sigDB = if(useDB) self$sigs else NULL
  if(!is.null(sigDB) ){
    eval1 = sigDB$loadEval(flags,phens,)
    if(!is.null(eval1) && nrow(eval1)>0){
      return(eval1)
    }
  }
  self$updateLOOC(phens, flags)
  if(length(all_modelsh$models)==0) return(NULL)
  #nme_d2 = .readFlag(flags,"test",names(self$datas))
  #names(nme_d2) = nme_d2
  all_models_y = all_modelsh$models#[[mod_nme]]
 # eval1 =  .merge1_new(lapply(nme_d2, function(nme1){
    #print(nme1)
    d = self$data
    eval1 = d$evaluateAllModels(all_models_y,phens,flags, verbose=verbose, inv_transform_y=inv_transform_y) %>% tibble::add_column(data=self$nme)
    #if(inherits(resd,"try-error")) {
    #  print(resd)
    #  print(paste("problem", nme1))
    #  stop("!!")
    #  return(NULL)
    #}
#    resd
#  }),addName="data")
  #}),addName="transform_y")
  if(is.null(eval1)) return(NULL)
  #  eval1 = subset(eval1, model!="avg")
  eval2 = eval1%>% pivot_wider(names_from="submeasure") #%>% tibble::add_column(transform_y=strsplit(transform_y[[1]]," ")[[1]][2])
  #  isfull=eval2$model %in% full_model_nmes
  #  eval2%>%tibble::add_column(isfull=isfull)
  if(!is.null(sigDB) ){
    sigDB$saveEval(eval2, flags,phens)
    eval1 = sigDB$loadEval(flags,phens)
    return(eval1)
  }
  #print("HH")
  #  if(T) return(eval2)
  .calcEval1(eval2, rename=F)
}
    
)
)
  