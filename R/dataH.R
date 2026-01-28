
.combineAngles1<-function(angleH, incl, flags,excl=list()){ 
  topn = .readFlag(flags,'topn', 20)
  onlyAll = .readFlag(flags,'only_all',F)
  angles1=angleH$angles;cols_incl1=angleH$cols_incl 
  nme_trans = names(angles1[[1]][[1]]); names(nme_trans) = nme_trans
  #  nmes_angs1 = names(angles1); names(nmes_angs1)=nmes_angs1
  names(incl) = incl
  comb_all2=lapply(nme_trans, function(nme_t1){
    nme_pow = names(angles1[[1]][[1]][[nme_t1]]); names(nme_pow)=nme_pow
    lapply(nme_pow, function(nme_p1){
      comb_all=lapply(incl, function(inc1){
        ang1 = angles1[[inc1]]
        if(is.null(ang1)) return(NULL)
        col_incl = cols_incl1[[inc1]]
        ang2=ang1[[1]][[nme_t1]][[nme_p1]]
        cs = colSums(ang2)
        if(length(ang1)>1){
          for(jk in 1:length(ang1)){
            cs = cs+colSums(ang1[[jk]][[nme_t1]][[nme_p1]])
          }
        }
        excl1 = excl[unlist(lapply(excl, function(ex) ex[3]==nme_t1 && ex[1] == inc1 && ex[4] ==nme_p1))]
        if(length(excl1)>0){
          col_incl[which(names(col_incl) %in% unlist(lapply(excl1, function(ex) ex[2])))]=F
        }
        cs[col_incl]
      })
      top_angles=whichpart1(comb_all, n=topn, return_scores=T)
      t1 = .merge1_new(lapply(top_angles, function(ta){
        data.frame(list(names = names(ta), value=ta))
      }),addName="data_type")
      t1 = t1[order(t1$value),]
      subset(t1, value<999)
    })
  })
}

## this is a class which holds a data object and interacts with the coordination node
.getAllSparseMatrices<-function(data, hasNA=T, convertToBigMatrix=F){
  rn = unlist(lapply(data, function(d1) rownames(d1)))
  rn = rn[!duplicated(rn)]
  print("getting sparse matrices")
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
    if(convertToBigMatrix){
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
 data_id="character",
 #transform_y="character",
 #var_t = "list",
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
    nme=sub("/",".",nme)
    self$flags = flags
   # transform_y =.readFlag(flags, "transform_y",toJSON(list(x=list(unvfunc="function(y,param) y",func="function(y,param) y", param=1))))
    ### MAKE SIGNATURE DIRECTORY
    self$nme=nme
    self$sigsdir=paste(dbDir,paste0("fspls_signatures__",nme,sep="/"))
    self$sigs = list()
    #####
    self$data = 
      dataObj$new(mat, nme,dbDir,flags,  
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
      dims1 = list(self$data$dims()); names(dims1) = self$nme
      self$sigs=   sigEnv$new(self$sigsdir,nme,flags, dims1, clear=F)
      
    
  },
 updateTransform=function(transform_y){
   ##probably no longer relevant
   stop("no longer relevant")
         transform_y1 = self$sigs$updateData(data_flags = self$flags, 
                                             data_names =self$nme, 
                                             data_types = names(self$data$data),
                                             dims = self$data$dims(),
                                             transform_y = fromJSON(transform_y)
         )
         transform_y = toJSON1(transform_y1)
         self$data$transforms = .convertToTransform(transform_y)
         self$transform_y = transform_y1
    },
 var_thresh = function(qq_t){
  lapply(self$data$vars, function(v) quantile(v, qq_t))
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
 update=function(phens, flags, verbose=F){
     self$updateLOOC(phens, flags, verbose=verbose)
     self$updateTrain(phens, flags,verbose=verbose)
   ##updated after updateLOOC
   nreps1 =self$nreps()
   nreps = 1:nreps1
   names(nreps) = nreps
   nreps
 },
 select=function(datasAll, phens,flags,
                 ## expt_id specific to this database .. might be diff for global
                 verbose=F, useDB=T ){#c(y="function(y) y","function(y) y")
   if(is.null(flags[['data_types']])) stop("cannot be NULL")
   nreps = self$update(phens, flags, verbose=verbose);
   if( useDB ){
     vars_all = self$sigs$loadVars(flags, phens)
     if(!is.null(vars_all)) return(vars_all)
   }
   vars_l_todo = datasAll$getTodo(flags, phens)
   expt_id=self$sigs$getExpt(flags, phens, add_new=T)
   variables=lapply(nreps, function(k1){
     if(verbose) print(paste("cv",k1,"of",length(nreps)))
     self$select_k(datasAll, phens,flags, k1, expt_id,
                   vars_l_todo,verbose=verbose)
   })
   vars_all=post_process(variables,flags,phens)
#   if(length(vars_all$variables)==0) return(vars_all)
   if(useDB){
     self$sigs$saveVars(vars_all,replace=T)   #saving local
   }
   vars_all
 },
 select_k=function(datasAll,phens,flags, k1,expt_id,
                   vars_l_todo ,
                   verbose=F){
   stop_y = .readFlag(flags, 'stop_y',"rand")
   logpvthresh = log(.readFlag(flags,"pthresh",0.1))
   beam= log(.readFlag(flags,"beam",1))
  saveAngles=F
  # vars_l = datasAll$nextVars(expt_id, flags)
   while(length(vars_l_todo$todo1)>0){
     
     comb_ = self$multiAnglesAndPv(phens, k1,flags,expt_id, vars_l_todo, saveAngles=saveAngles, verbose=verbose)
     data_nme=self$nme
     vars_l_todo=datasAll$savePvalsAndNextVars(flags,phens,vars_l_todo,comb_,data_nme,  k1,logpvthresh,beam)
     #datasAll$savePvals(expt_id,k1, self$nme, vars_l_todo$vars_l,comb_)
     #vars_l_todo = datasAll$nextVars(vars_l_todo, expt_id, k1,logpvthresh,beam)
   }
  vars_l_todo$vars_l 
 },
 res_inner=function(comb_,prev_i, flags,k, expt_id){
  prev_i2 = self$sigs$loadPrev(expt_id, prev_i, k, data_nme = self$nme)
  if(is.null(prev_i2)) prev_i2 = prev_i
   nme_comb = names(comb_); names(nme_comb) = nme_comb
   #nme_c1 = nme_comb[[1]]; nme_p1 = names(comb_[[nme_c1]])[[1]]; ik=1
   res_inner=lapply(nme_comb, function(nme_c1){
     nmesp1 = names(comb_[[nme_c1]]); names(nmesp1) = nmesp1
     lapply(nmesp1, function(nme_p1){
       comb = comb_[[nme_c1]][[nme_p1]]
       if(nrow(comb)==0) return(NULL)
       num_pvals1 = nrow(comb)
       inds1p = 1:num_pvals1; names(inds1p) = comb$names[1:length(inds1p)]
       nxt_vars = lapply(inds1p, function(ik){
         b_i_name = c(comb$data_type[[ik]], comb$names[[ik]], nme_c1,nme_p1)
          nv= self$getPvsAll(phens,prev_i2, b_i_name,k,  prev_i2$Wall,flags)
         if(inherits(nv,"try-error")) {
           print(paste(nme_c1, "error"))
           return(NULL)
         }
         nv
       })
      
     })
   })
res_inner   
 },
multiAnglesAndPv=function(phens,  k1,flags, expt_id, vars_l_todo,
                          saveAngles=F, verbose=F){
  if(is.null(expt_id)) stop("expt_id is NULL")
  vars_l = vars_l_todo$vars_l
  todo1=vars_l_todo$todo1[[1]]
  incl=todo1$incl
  g_incl = todo1$g_incl
  qq_t = todo1$qq
   invisible( lapply(vars_l, function(prev_i){
      varnames = prev_i$var_names; 
      comb_ = self$anglesAndPv(phens, prev_i, incl, k1, g_incl, qq_t, flags,expt_id, saveAngles=saveAngles, verbose=verbose)
      comb_
    }))
},
anglesAndPv=function(phens, prev_i, incl, k1, g_incl, qq_t, flags, expt_id, saveAngles=F, verbose=F){
  varnames = prev_i$var_names
 
  comb_=self$combinedAngles(phens, varnames, incl, k1,  g_incl, qq_t, flags)
  if(saveAngles) return(comb_)
    #self$sigs$saveAngles(expt_id, data_nme, comb_angs1,varnames ) 
  ri = self$res_inner(comb_,prev_i,flags,k1, expt_id)
  self$sigs$savePvals(flags,phens, self$nme, ri, varnames,k1,useCurrVarnames=T)
  ri_out=lapply(ri, function(ri1){
    lapply(ri1, function(ri2){
      lapply(ri2, function(ri3){
        ri4 = ri3$simplify()
        list(pvs = ri3$pvs, pvs_all = ri3$pvs_all, var_names = ri3$var_names)
      })
    })
  })
  ri_out
},
 combinedAngles=function(phens, varnames, incl, k, g_incl, qq_t, flags){ #phens, varnames, incl=incl, k=k, type=type
  type=self$type
  var_t = self$var_thresh(qq_t)
   angleH=list(angles=
                 self$data$getAngles1(phens,varnames,incl=incl,k=k, type=type),
               cols_incl = self$data$cols_incl(var_t,incl, g_incl,excl=varnames)) ### fix 
   .combineAngles1(angleH, incl, flags, excl=varnames)
 },
 getPvsAll=function(subphens, prev_i, b_i_name,k, #   prev_i = vars_l1[[nmed]]
                      Wall =lapply(subphens, function(f) matrix(nrow=0,ncol=0)),
                    flags){
                      #useglm=F ,inv_transform=getOption("x_transform",T),
                      #project=T, useoffset=T){
   inv_transform=T
   project=.readFlag(flags,"project",T)
   useoffset=.readFlag(flags,"useoffset",T)
   useglm = .readFlag(flags,'useglmnet',T)
   
     d = self$data
     family = strsplit(names(subphens)[[1]],"\\.")[[1]][1]
     if(family=="multinomial") useoffset=F
     prev_i1 = d$makeNextModel(prev_i,b_i_name,subphens,k, family, ypred=NULL, 
                               project=project, useglm=useglm, logpthresh =0, useoffset=useoffset)
     prev_i1
#   pvs_all = pvs_all[unlist(lapply(pvs_all, length))>0]
  # pvs_all
 },
extractPredictions=function(all_modelsh,phens=all_modelsh$phens, flags=all_modelsh$flags, CV = FALSE, liab=T){
  self$updateLOOC(phens, flags)
  all_models_y = all_modelsh$models
  
  res3 = lapply(all_models_y, function(all_models_){
      d = self$data
      if(is.null(d)){
        print(d_nme)
        stop("!!")
      }
       d$extractPredictions(all_models_, phens, flags, CV=CV, liab=liab)
  })
  predictions0=res3[unlist(lapply(res3, function(x) length(x[[1]])))>0]
  predictions0 # 
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
#  if(length(which(train_nme %in% names(self$datas)))==0)train_nme = names(self$datas)[[1]]
  verbose=.readFlag(flags,"verbose",FALSE)
  if(!is.null(flags[['useglm']])) stop("define useglmnet not useglm")
  useglm=.readFlag(flags,"useglmnet",TRUE)
  inds1 = inds#[[nmes_inds1]]
  phens1 = phens#[[nmes_inds1]]
  #k=inds1[[1]]; d = datas[[1]]
  mods1 = lapply(inds1, function(k){
    # print(k)
   # lapply(datas[names(datas) %in% train_nme], function(d){
      mods = d$makeModels(phens1, vars2,k,logpthresh = logpthresh,project=project,
                          flags=flags,checkRMSV=F,
                          useglm=useglm, useoffset=useoffset)
    mods
    #})
    #})
  })
  if(length(mods1)==0) stop("!!")
  models=mods1
  #})
  vars = names(models[[1]])
  names(vars) = vars
  models2 = lapply(vars, function(v){
    # lapply(models, function(models1){
    m3 = lapply(models, function(m){
      m[[v]]
      #m2[unlist(lapply(m2, length))>0]
    })
    m3[unlist(lapply(m3, length))>0]
    #})
  })
  #  print(names(models2))
  models2
},
makeAllModels=function(vars_all0, 
                       phens=vars_all0[[1]]$phens, flags=vars_all0[[1]]$flags, verbose=F, max = 1e6,
                       user="",useDB=T){
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
  beams = names(vars_all0); names(beams)=beams
  all_models_full=lapply(beams, function(beam){
    vars_all = vars_all0[[beam]]
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
  #v_nme = names(vars_all$variables)[1]; max=10; verbose=T; k=1;variables =vars_all$variables; 
  
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
  all_models

  })
  all_models_=list(models=all_models_full, flags = flags, phens = phens, trainedOn=self$nme)
  if(useDB ){
    sigDB$saveModels (all_models_)
  }
  #combined_models
  all_models_
},
evaluateAllModels=function(all_modelsh, phens=all_modelsh$phens,flags=all_modelsh$flags,verbose=F,useDB=T, user=""){ ## different folds with same variables
  sigDB = if(useDB) self$sigs else NULL
  if(!is.null(sigDB) ){
    eval1 = sigDB$loadEval(flags,phens,)
    if(!is.null(eval1) && nrow(eval1)>0){
      return(eval1)
    }
  }
  inv_transform_y=F
  self$updateLOOC(phens, flags)
  if(length(all_modelsh$models)==0) return(NULL)
  #nme_d2 = .readFlag(flags,"test",names(self$datas))
  #names(nme_d2) = nme_d2
  all_models_y = all_modelsh$models#[[mod_nme]]
 # eval1 =  .merge1_new(lapply(nme_d2, function(nme1){
    #print(nme1)
    d = self$data
    eval1 =   .merge1_new(lapply(all_models_y, function(all_models_y1){
   d$evaluateAllModels(all_models_y1,phens,flags, verbose=verbose) %>% tibble::add_column(data=self$nme, trainedOn=all_modelsh$trainedOn)#%>% tibble::add_column(trainedOn=self$nam)
  }), addName="beam")  #if(inherits(resd,"try-error")) {
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
  
