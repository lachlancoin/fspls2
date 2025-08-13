
.getPvsAll<-function(subphens,datas1, vars_l1, b_i_name,k, funcst, project=T, useoffset=T){
  transform_func = eval(str2lang(funcst))
  #var_new = c(prev_i$var,b_i)
  nmesd = names(datas1); names(nmesd)=nmesd
  #nmed = nmesd[[1]]
  pvs_all= lapply(nmesd, function(nmed){
    d = datas1[[nmed]]
    prev_i = vars_l1[[nmed]]
    family = strsplit(names(subphens)[[1]],"\\.")[[1]][1]
    if(family=="multinomial") useoffset=F
    prev_i1 = try(d$makeNextModel(prev_i,b_i_name,subphens,k, transform_func,family, ypred=NULL, project=project, useglm=T, logpthresh =0, useoffset=useoffset))
    if(inherits(prev_i1,"try-error")) {
    #  print(paste("problem", nmed))
      return(NULL)
    }
     prev_i1
  })
  pvs_all = pvs_all[unlist(lapply(pvs_all, length))>0]
  pvs_all
}
.nonZero<-function(am2){
  am2[unlist(lapply(am2, length))>0]
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
    return(lapply(y_mat,function(y){
      if(!is.numeric(y)){
        y = as.factor(y);
      }
      if(is.factor(y)){
        return(if(length(levels(y))<=2 )"binomial" else "multinomial")
      }
      if(length(unique(y[!is.na(y)]))<=2) return("binomial")
      vals = unique(y)
      if(sum(abs(vals-round(vals)))<1e-9) {
        if(length(table(vals))>max_ordinal) return("gaussian")
        return("ordinal")
      }
      return("gaussian")
    }))
  }else{
   
  return(apply(y_mat,2,function(y){
    if(!is.numeric(y)){
      y = as.factor(y);
    }
    if(is.factor(y)){
       return(if(length(levels(y))<=2 )"binomial" else "multinomial")
    }
    if(length(unique(y[!is.na(y)]))<=2) return("binomial")
    vals = unique(y)
    if(sum(abs(vals-round(vals)))<1e-9) return("ordinal")
    return("gaussian")
    }))
  }
}

##this function removes NAs
## if no NA matrixNA is just empty matrix
.getSparseMatrices<-function(mat, hasNA=T){
  if(!hasNA){
   res1 = list(matrix = mat,
               matrixNA = Matrix(0,nrow(mat) , ncol(mat), sparse = T)
                    )
  }else{
  
  res1 = list(
    matrix = Matrix(apply(mat,2,function(v){
      mv = mean(v, na.rm=T)
      v[is.na(v)]=mv
      v
    }), sparse=T),
    matrixNA = Matrix(apply(mat,2,function(v){
      v1 = rep(0, length(v))
      v1[is.na(v)]=1
      v1
    }), sparse=TRUE)
  )
  }
  res1
}

datasEnv<-R6Class("datasEnv", public = list(
  datas = "list",
  type="character",
  sigsdir="character",
  sigs="list",
  flags ="list",
  initialize=function(
              datasets,
              ys = lapply(datasets, function(d) d$y),
              flags = list(),
              mats = lapply(datasets, function(d) lapply(d$data, function(d1).getSparseMatrices(d1))),
              families=lapply(ys, function(d) .getFamily(d)),
              dbDir="./",
                      memDir=NULL){
    self$flags = flags
    ### MAKE SIGNATURE DIRECTORY
    self$sigsdir=paste(dbDir,"fspls_signatures",sep="/")
    
    self$sigs = list()
   
    #####
    
    if(!is.null(flags$transform) & typeof(mats[[1]][[1]][[1]])=="S4"){
      transforms =fromJSON (flags$transform)
      trans_matrs = lapply(mats,function(mats1){
        tm1 = unlist(lapply(mats1, function(mats2){
          aa=lapply(transforms, function(str){
            print(str)
            func=  eval(str2lang(str))
            range=-5:5
            test_v = apply(cbind(range,func(range)),1,function(v) abs(v[2]-v[1]))
            if(max(test_v,na.rm=T)<0.001){
              ##is identity
              return(mats2) 
            }
            mat3 =  func(mats2$matrix)
            na_ind = which(is.na(mat3))
            na_m = mats2$matrixNA
            if(length(na_ind)>0){
              stop("added NAs")
              na_m[na_ind] = 1
            }
            list(matrix = func(mats2$matrix), matrixNA = na_m)
          })
        }),rec=F)
        names(tm1) = gsub("\\.","_",names(tm1))
        tm1
      })
      mats = trans_matrs
    }
    
    
    datas = lapply(names(mats), function(nme){
      dataObj$new(mats[[nme]], incl_full=T,seed = getOption("seed",42), memDir=if(is.null(memDir)) NULL else paste(memDir, nme,sep="/"))
    })
  
    names(datas) = names(mats)
    self$type="slow1"
    types_all = getOption("types_all",names(datas[[1]]$data))
    names(types_all) = types_all
 #   var_threshs=  lapply(types_all, function(v) .readFlag(flags, "var_quantile",0.00))
  #  
    batch=.readFlag(flags, "batch",0)
    all_v_all = .readFlag(flags,"all_v_all",F)
    one_v_rest = .readFlag(flags,"one_v_rest",F)
    nrep = .readFlag(flags,"nrep",if(batch>0) 0 else 1)
    invisible(lapply(datas, function(data) data$init1(nrep=nrep,  batch = batch)))
    varn = getOption("varn",c())
   
    invisible(lapply(1:length(datas), function(ik) {
      y1 = ys[[ik]] #dists[[ik]]$updateYdb(cats[['cats']])
      family = families[[ik]]
      datas[[ik]]$updateY(y1, family=family, CHECK=T, all_v_all=all_v_all, one_v_rest = one_v_rest)
  #    datas[[ik]]$updateYdb(dists[[ik]]$mydb, cats[['cats']])
  #    missing_vals = self$updateY(y1, family=family, CHECK=T)
      
      datas[[ik]]$initTrain(varn=varn)
    #data[[ik]]$initY()
     return(NULL)
    }))
    self$datas = datas
  },
  updateWeights=function(subphens = self$pheno()[[1]][1]){ ## upweights low count values
   for(k in 1:length(self$datas)){
     self$datas[[k]]$updateWeights(subphens)
   } 
  },
   getSigDB=function(nme1,reload= F, clear=F){
     if(is.null(nme1)) return(NULL)
   if(reload || is.null(self$sigs[[nme1]])){
     self$sigs[[nme1]]=   sigEnv$new(self$sigsdir,nme1, clear=clear)
     self$sigs[[nme1]]$updateData(data_flags = self$flags, 
                                  data_names =names(self$datas), 
                                  data_types = names(self$datas[[1]]$data),
                                  dims = self$dims(),
                                  phenos =self$datas[[1]]$pheno(), user="")
     #self$sigs[[nme1]]$data_id
   }
   self$sigs[[nme1]]
 },
  update=function( flags = list()){
    types_all = getOption("types_all",names(self$datas[[1]]$data))
    names(types_all) = types_all
  #  var_threshs=  lapply(types_all, function(v) .readFlag(flags, "var_quantile",0.00))
    batch=.readFlag(flags, "batch",0)
    
    nrep = .readFlag(flags,"nrep",if(batch>0) 0 else 1)
    invisible(lapply(self$datas, function(data) data$init1(nrep=nrep,  batch = batch)))
    varn = getOption("varn",c())
    for(ik in 1:length(self$datas)){
      self$datas[[ik]]$initTrain(varn=varn)
    }
  },
  
  randomise=function(){
    for(i in 1:length(self$datas)){
      self$datas[[i]]$randomise();
    }
  },
  
  
  
  
  
  dims=function(){
    lapply(self$datas, function(data){
      list(
        nonNA=lapply(data$data, function(d){
        list(dim = dim(d),bigmatrix=is.big.matrix(d))
      }),
      na=lapply(data$dataNA, function(d){
        list(dim = dim(d),bigmatrix=is.big.matrix(d))
      }))
    })
  },
 cats = function(maxpheno = 1e9){
    lapply(self$datas, function(d) d$cats(maxpheno))
 },
  pheno=function(maxpheno=1e9,sep=F, sep_group = F, exclude=NULL){
   res = self$datas[[1]]$pheno(maxpheno=maxpheno, sep=sep, sep_group = sep_group);
  if(!is.null(exclude)){
    lapply(res, function(res2){
    lapply(res2, function(res1) res1[-grep(exclude,res1)])
    })
  }
   res
  },

  angles=function(vars,phens,flags){
    datas = self$datas
    topn = .readFlag(flags,'topn', 100) 
    train_nme = .readFlag(flags,'train', names(self$datas))
    names(train_nme) = train_nme
    #  type=.readFlag(flags, 'type', 'slow1') 
    incl = .readFlag(flags,'data_types',names(datas[[1]]$data))
    names(incl )= incl
    k = .readFlag(flags, 'rep',length(datas[[1]]$train)) ## this is the full data, not cv
    angles_all = lapply(vars, function(vars1){
      angles=lapply(train_nme, function(data_nme) datas[[data_nme]]$getAngles1(phens,vars1,incl=incl,k=k, type=self$type))
      cols_incl = lapply(train_nme, function(data_nme)datas1[[data_nme]]$cols_incl)
      comb=.combineAngles(angles,cols_incl,topn=topn)
      as.list(comb)       
    })
    angles_all
  },
  makeAllModels=function(vars_all, phens, flags, verbose=F, max = 1e6,db=NULL){
    
    sigDB = self$getSigDB(db)
    
    if(!is.null(sigDB) ){
      all_models = sigDB$loadModels(flags,phens)
      if(!is.null(all_models) && length(all_models)>0){
        return(all_models)
      }
    }
    
    nmes_vars_all = names(vars_all); names(nmes_vars_all) = nmes_vars_all
    func_strs = fromJSON(.readFlag(flags,"transform_y",'{"y":"function(y) y"}'))
    #nme_v_all = nmes_vars_all[[1]]
    
    logpthresh= log(.readFlag(flags,"pthresh",1e-3))
    project=.readFlag(flags,"project",TRUE)
    #nme_v_all = nmes_vars_all[[1]]; v_nme = names(vars_all[[1]]$variables)[1]; max=10; verbose=T; k=1;variables =vars_all[[nme_v_all]]$variables;  v_nme = names(variables)[[1]]
    combined_models=lapply(nmes_vars_all, function(nme_v_all){
      if(verbose) print(nme_v_all)
      vars = vars_all[[nme_v_all]]
      func_str = func_strs[[nme_v_all]]
      all_models = list()
      variables = vars$variables
      var_inds = vars$inds
    #  print(var_inds);stop("!!")
      if(length(variables)==0) return(list())
      ord = order(unlist(lapply(variables, length)),decreasing=T)
      variables = variables[ord]
      var_inds = var_inds[ord]
     for(v_nme in names(variables)){
     #  print(v_nme)
       if(verbose)print(v_nme)
       vars2 = variables[[v_nme]]
       vars2 = vars2[1:min(length(vars2), max)]
       inds =var_inds[[v_nme]]
       nme_ = paste(names(vars2),collapse=";")
       models1 = all_models[[nme_]]
       if(is.null(models1)){
        models1 = self$makeModels( vars2, inds,phens,func_str,flags)
        for(k in 1:length(models1)){
          mod1 =   all_models[[names(models1)[[k]]]]
          if(is.null(mod1)){
             all_models[[names(models1)[[k]]]] = models1[[k]]
          }else{
            for(p_nme in names(models1[[k]])){
              mod2 = all_models[[names(models1)[[k]]]][[p_nme]]
              if(is.null(mod2)){
                all_models[[names(models1)[[k]]]][[p_nme]] = models1[[k]][[p_nme]]
              }else{
                for(r_nme in names(models1[[k]][[p_nme]])){
                   mod3 = all_models[[names(models1)[[k]]]][[p_nme]][[r_nme]]
                   if(is.null(mod3)){
                     all_models[[names(models1)[[k]]]][[p_nme]][[r_nme]] = models1[[k]][[p_nme]][[r_nme]]
                   }else{
                     print("already calculated!")
                   }
                }
              }
            }
          }
        }
      }else{ ## fill in gaps
        for(p_i in 1:length(inds)){
          p_nme = names(inds)[p_i]
          models2 = all_models[[nme_]][[p_nme]]  
          if(is.null(models2)){
            models2 = self$makeModels( vars2, inds[p_i],phens[which(names(phens) %in% p_nme)],func_str,flags)
            for(k in 1:length(models2)){
              all_models[[names(models2)[[k]]]][[p_nme]] = models2[[k]][[p_nme]]
            }
          }else{
            for(r_i in 1:length(names(inds[[p_i]]))){
              r_nme = names(inds[[p_i]])[r_i]
              models3 = all_models[[nme_]][[p_nme]][[r_nme]]  
              if(is.null(models3)){
                inds2 = lapply(inds[p_i], function(indsk)indsk[r_i])
                models3 = self$makeModels( vars2, inds2,phens[which(names(phens) %in% p_nme)],func_str, flags)
                for(k in 1:length(models3)){
                  all_models[[names(models3)[[k]]]][[p_nme]][[r_nme]] = models3[[k]][[p_nme]][[r_nme]]
                }
              }
            }
          }
        }
       
      }
    }
    all_models
    })
    if(!is.null(sigDB) ){
      sigDB$saveModels (combined_models, flags, phens)
    }
    combined_models
  },
  makeModels=function(vars2, inds, phens,func_str, flags){
    datas=self$datas
#    print(inds)
    logpthresh= log(.readFlag(flags,"pthresh",1e-3))
    project=.readFlag(flags,"project",TRUE)
    useoffset=.readFlag(flags,"useoffset",TRUE)
    train_nme = .readFlag(flags,'train', names(datas)[1])
    if(length(which(train_nme %in% names(self$datas)))==0)train_nme = names(self$datas)[[1]]
    verbose=.readFlag(flags,"verbose",FALSE)
    #k = .readFlag(flags, 'rep',length(datas[[1]]$train))
    nmes_inds = names(inds);names(nmes_inds)=nmes_inds
    #nmes_inds1 = nmes_inds[[1]];d =datas[[1]]
    models = lapply(nmes_inds, function(nmes_inds1){
      inds1 = inds[[nmes_inds1]]
      phens1 = phens[[nmes_inds1]]
      #k=inds1[[1]]
        mods1 = lapply(inds1, function(k){
            lapply(datas[names(datas) %in% train_nme], function(d){
              mods = d$makeModels(phens1, vars2,k,logpthresh = logpthresh,project=project, func_str=func_str, useoffset=useoffset)
            })
          #})
        })
        if(length(mods1)==0) stop("!!")
        mods1
    })
    vars = names(models[[1]][[1]][[1]])
    names(vars) = vars
    models2 = lapply(vars, function(v){
      lapply(models, function(models1){
      m3 = lapply(models1, function(m){
        m2 = lapply(m, function(m1) m1[[v]])
        m2[unlist(lapply(m2, length))>0]
      })
      m3[unlist(lapply(m3, length))>0]
    })
    })
  #  print(names(models2))
    models2
  },
  getProjectedData=function(varnames){
    lapply(self$datas, function(d){
      d$getProjectedData(varnames = varnames);      
    })
  },
  getVariance=function(varnames){
    lapply(self$datas, function(d){
      d$getVariance();      
    })
  },
  post_process=function(variables, flags_out = list()){
    full_index = length(variables)
    #lens = unlist(lapply(variables, function(x) length(x$var)))
    #if(max(lens)==0) return(list())
    vars_all = list()
    vars_all1 = list()
    vars_all2 = list()
    names(variables) = 1:length(variables)
    for(repn in names(variables)){
      full = repn==full_index
      var1 = variables[[repn]]
      for(phenn in names(var1)){
        var2 = var1[[phenn]][[1]]$var_names
        names(var2) = var1[[phenn]][[1]]$varnames
        cumpv = lapply(var1, function(vv) attr(vv,"cumpv"))
        varn = paste(names(var2),collapse=";") #paste(names(var2), collapse=";")
        if(is.null(vars_all[[varn]])){
          vars_all[[varn]] = list()
          vars_all1[[varn]] = var2
          vars_all2[[varn]] = list()
        }
        repn1 = as.list(as.numeric(repn))
        repn2 = as.list(cumpv)
        names(repn1) = if(full) "full" else repn
        names(repn2) = if(full) "full" else repn
        if(is.null(vars_all[[varn]][[phenn]])){
          vars_all[[varn]][[phenn]] =repn1
          vars_all2[[varn]][[phenn]] =repn2
        }else{
          vars_all[[varn]][[phenn]] = c(vars_all[[varn]][[phenn]] , repn1)
          vars_all2[[varn]][[phenn]] = c(vars_all2[[varn]][[phenn]] , repn2)
        }
      }
    }
    vars_combined = list(variables = vars_all1, inds = vars_all,cumpv=vars_all2, flags=flags_out) 
  
    vars_combined
  },
  convert1=function(variables, phens){
    nreps1 =unlist(lapply(self$datas, function(x) ncol(x$looc$incl)))
    inds = lapply(variables, function(v){
      lapply(phens, function(phens1){
        ii = 1:nreps1[[1]]
        names(ii)= ii
        names(ii)[[length(ii)]]="full"
        as.list(ii)
      })
    })
    vars_ = list(variables = variables, inds = inds)
    vars_
  },
  convert=function(genes_incl, phens){
    variables = lapply(names(self$datas[[1]]$data), function(nme){
      d1 = self$datas[[1]]$data[[nme]]
      genes1=genes_incl[which(genes_incl %in% colnames(d1))]
      vars2 = lapply(genes1,function(x){
        c(nme,x)
      })
      names(vars2) = lapply(vars2, function(vv)paste(vv,collapse="."))
      vars2
    })
    names(variables) = lapply(variables, function(l2)paste(names(l2), collapse=";"))
  result=self$convert1(variables, phens)
  func_str = fromJSON(.readFlag(flags,"transform_y",'{"y":"function(y) y"}'))
  lapply(func_str, function(xx) result)
  },
  select=function(phens,flags,verbose=F, db=NULL ){
    sigDB = self$getSigDB(db)
    if(!is.null(sigDB) ){
      vars_all = sigDB$loadVars(flags, phens)
      if(!is.null(vars_all)){
        return(vars_all)
      }
    }
    nreps1 =ncol(self$datas[[1]]$looc$incl)
    datas1 = self$datas
    project=.readFlag(flags,"project",T)
    useoffset=.readFlag(flags,"useoffset",T)
    topn = .readFlag(flags,'topn', 20)
    onlyAll = .readFlag(flags,'only_all',F)
    train_nme = .readFlag(flags,'train', names(datas1))
    train_nme = train_nme[train_nme %in% names(datas1)]
    quantiles = sort(fromJSON(.readFlag(flags, "quantiles","[0]")),decreasing=T)
    genes_incls=fromJSON(.readFlag(flags,"genes_incls",'{"all":["all"]}')) #,getOption("genes_incls",NULL)
    if(length(train_nme)==0) train_nme = names(datas1)[[1]]
    names(train_nme) = train_nme
    maxsize=.readFlag(flags,'max',50)
    num_pvals = min(topn, 10)
    incls = fromJSON(.readFlag(flags,'data_types',"{}")) 
    if(length(incls) == 0 )incls = list("all"=names(datas1[[1]]$data))
    logpvthresh = log(.readFlag(flags,"pthresh",1e-5))
    logpv=-100
    nreps = 1:nreps1
    names(nreps) = nreps
    beam = .readFlag(flags,"beam",1)
    func_str = fromJSON(.readFlag(flags,"transform_y",'{"y":"function(y) y"}'))
    flags_out = list(train=train_nme,  max=maxsize, 
                     nreps = nreps, beam=beam,func_str = func_str,
                     topn = num_pvals, pthresh =exp(logpvthresh), data_types = incls)
    phens_index = 1:length(phens)
    names(phens_index) = names(phens)
    var_thresh = lapply(train_nme, function(data_nme){
      lapply(datas1[[data_nme]]$vars, function(v) quantile(v, quantiles))
    })
   #funcst = func_str[[1]]; k=1;p_index = phens_index[[1]]; g_incl  = genes_incls[[1]]; qq =1; incl = incls[[1]]; data_nme = train_nme[[1]]
    vars_combined=lapply(func_str,function(funcst){
      print(funcst)
     variables=lapply(nreps, function(k){
      print(paste("cv",k,"of",length(nreps)))
      jj1=0
      incls_all = unique(unlist(incls))
      #invisible(lapply(train_nme, function(data_nme) datas1[[data_nme]]$update(k, funcst, phens,incls_all))); ### update training object
      res2=  lapply(phens_index, function(p_index){
        if(verbose) print(paste("pindex",p_index,k,length(nreps)))
        subphens = phens[[p_index]]
        invisible(lapply(train_nme, function(data_nme) datas1[[data_nme]]$update(k, funcst, subphens,incls_all))); ### update training object
        if(FALSE) cat(p_index); cat("\t");
        vars_l =list(lapply(train_nme,function(xx) stateObj$new(subphens, NULL,NULL,NULL,NULL,k, var=c(), varnames=c(), W_all = NULL)))
        #vars_l1 = vars_l[[1]]
        for(incl in incls){
          if(verbose) print(incl)
        for(g_incl in genes_incls){
         for(qq in 1:length(quantiles)){
         while(logpv<logpvthresh && length(vars_l[[1]][[1]]$var)<maxsize){
          angles_all = lapply(vars_l, function(vars_l1){
            nxt_vars1 =  tryCatch({
              varnames = vars_l1[[1]]$var_names; type = self$type
            angles=lapply(train_nme, function(data_nme) {
              #print(data_nme); 
              datas1[[data_nme]]$getAngles1(subphens,varnames,incl=incl,k=k, type=type)
              })
            names(angles) = train_nme
            cols_incl = lapply(train_nme, function(data_nme)datas1[[data_nme]]$cols_incl(var_thresh[[data_nme]],incl, g_incl,qq)) ### fix 
             comb=.combineAngles(angles,cols_incl,incl,topn=topn, onlyAll = onlyAll)
         
            if(nrow(comb)==0) return(NULL)
            num_pvals1 = min(num_pvals, nrow(comb))
            inds1p = 1:num_pvals1; names(inds1p) = comb$names[1:length(inds1p)]
            nxt_vars = lapply(inds1p, function(ik){
            #  print(ik)
               b_i_name = c(comb$data_type[[ik]], comb$names[[ik]])
               nv = .getPvsAll(subphens,datas1[names(datas1) %in% train_nme], vars_l1, b_i_name,k, funcst, project = project, useoffset=useoffset)
              attr(nv,"cumpv")= .sumChisq(unlist(lapply(nv, function(nv1){
                 unlist(nv1$pvs)
               })))
              nv
              #mStateObj$new(comb[ik],  .sumChisq(pv) , prev_i=prev_i)
            })
            pvs_list = unlist(lapply(nxt_vars, function(nv) attr(nv,"cumpv")))
            #print(pvs_list)
            #subinds1 = pvs_list<=logpvthresh
            #if(length(which(subinds1)==0)) return(NULL)
            #nxt_vars = nxt_vars[subinds1]
            #pvs_list = pvs_list[subinds1]
            nxt_vars[order(pvs_list)] 
            },error=function(w){
              print(w)
              print("error")
              return(NULL)
            })
            nxt_vars1
            
          })
          angles_all = angles_all[unlist(lapply(angles_all,length))>0]
          angles_all1 = unlist(angles_all, rec=F)
          ord = order(unlist(lapply(angles_all1, function(nv)attr(nv,"cumpv"))))
          logpv = attr(angles_all1[[ord[1]]],"cumpv")
          #print(logpv)
          if(logpv<=logpvthresh){
            vars_l = angles_all1[ord][1:beam]
          }
          if(verbose){
            print(names(vars_l))
             print(paste("logpv",logpv,jj1))
             jj1 = jj1+1
          }
        }
        }
        }
      }
       
        #lapply(datas1, function(d)d$saveParquet())
        #lapply(vars_l, function(v) v$var)
          ##just take the top
       # if(length(vars_l[[1]]$var)>0) print(vars_l)
        vars_l[[1]]#$var
        })
         res2[unlist(lapply(res2,function(xx) length(xx[[1]]$var)))>0]
    })
    self$post_process(variables,flags_out)
  })
    if(!is.null(sigDB) ){
      sigDB$saveVars(vars_combined, flags,phens)
    }
 
   vars_combined
  },
  extractPredictions=function(all_models,phens, flags, CV = FALSE, liab=T, data_nme  = names(self$datas)){
    #datas = self$datas
    names(data_nme) = data_nme
    nmes_p = names(phens); names(nmes_p) = nmes_p
    #all_models_ = all_models[[1]]; nme_p = nmes_p[[1]]; d_nme = data_nme[[1]];phens1 = phens[[nme_p]]
    res3 = lapply(all_models, function(all_models_){
      res2 = lapply(nmes_p, function(nme_p){
        phens1 = phens[[nme_p]]
        res2=lapply(data_nme, function(d_nme){
          d = self$datas[[d_nme]]
          if(is.null(d)){
            print(d_nme)
            stop("!!")
          }
          d$extractPredictions(all_models_, phens1, nme_p,flags, CV=CV, liab=liab)
        })
      })
   
    res2[lapply(res2,length)>0]
    })
    res3
 },
  evaluateAllModels=function(all_models, phens,flags, db=NULL){ ## different folds with same variables
##                          inds = as.numeric(names(all_models))){
    #self = all_models
    sigDB = self$getSigDB(db)
    
    if(!is.null(sigDB) ){
      eval1 = sigDB$loadEval(flags,phens, user="")
      if(!is.null(eval1) && nrow(eval1)>0){
        return(eval1)
      }
    }
    if(length(all_models)==0) return(NULL)
    nme_d = .readFlag(flags,"test",names(self$datas))
    names(nme_d) = nme_d
    print(nme_d)
    func_strs = fromJSON(.readFlag(flags,"transform_y",'{"y":"function(y) y"}'))
    inverse_func_strs = fromJSON(.readFlag(flags,"transform_y_inverse",'{"y":"function(y) y"}'))
    #k=1
    for(k in 1:length(func_strs)){ #CHECK INVERSE CORRECT
      func0 = eval(str2lang(func_strs[[k]]))  ## should be inverse
     func1 = eval(str2lang(inverse_func_strs[[k]]))  ## should be inverse
     xx = -10:10  
     if(max(abs(apply(cbind(xx,func1(func0(xx))),1,diff)),na.rm=T)>1e-5) stop("please provide inverse function")
    }
    
    mod_nmes = names(all_models); names(mod_nmes) = mod_nmes
    #mod_nme = mod_nmes[[1]]; nme1 = nme_d[[1]] 
    eval1=.merge1_new( lapply(mod_nmes, function(mod_nme){
      inverse_func_str = inverse_func_strs[[mod_nme]]
      all_models_y = all_models[[mod_nme]]
     .merge1_new(lapply(nme_d, function(nme1){
     d = self$datas[[nme1]]
      resd = try(d$evaluateAllModels(all_models_y,phens,inverse_func_str = inverse_func_str,flags))
      if(inherits(resd,"try-error")) {
        print(paste("problem", nme1))
        return(NULL)
      }
      resd
    }),addName="data")
    }),addName="transform_y")
    if(is.null(eval1)) return(NULL)
  #  eval1 = subset(eval1, model!="avg")
    eval2 = eval1%>% pivot_wider(names_from="submeasure")
  #  isfull=eval2$model %in% full_model_nmes
  #  eval2%>%tibble::add_column(isfull=isfull)
    if(!is.null(sigDB) ){
      sigDB$saveEval(eval2, flags,phens, user="")
      eval1 = sigDB$loadEval(flags,phens, user="")
      return(eval1)
    }
    .calcEval1(eval2)
  },
  pvalues=function(vars,phens,transform_y,flags){
    datas = self$datas
    transform_func = eval(str2lang(transform_y))
    
    k = .readFlag(flags, 'rep',length(datas[[1]]$train))
    lapply(vars,function(vars1) {
      lapply(datas, function(d){
        b_i1 = vars1[[length(vars1)]]
        prev_var =  vars1[-length(vars1)]
        
        sig_res = d$calcBetaProj1(phens,k, b_i1, prev_var, transform_func, convert=T )
        unlist(sig_res$pvs)
        #d$train[[k]]$getPvs(prev)
      }) 
    })
    #.getPvsAll(phens,datas,vars1,k))
  }
))
