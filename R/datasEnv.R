
.getPvsAll=function(phens,datas, vars,k){
  pvs_all= lapply(datas, function(d){
    b_i1 = vars[[length(vars)]]
    prev_var =  vars[-length(vars)]
    sig_res = d$calcBetaProj1(phens,k, b_i1, prev_var , convert=T)
    unlist(sig_res$pvs)
    #d$train[[k]]$getPvs(prev)
  })
  pvs_all
}
.getFamily<-function(y_mat){
  types = attr(y_mat, "types")
  if(!is.null(types)){
    return(lapply(types, function(typ){
      if(typ=="double") return("gaussian")
      if(typ=="boolean") return("binomial")
      if(typ=="integer") return("ordinal")
      if(typ=="character") return("multinomial")
    }))
  }
    if(typeof(y_mat)=="list"){
    return(lapply(y_mat,function(y){
      if(!is.numeric(y)){
        y = as.factor(y);
      }
      if(is.factor(y)){
        return(if(length(levels(y))==2 )"binomial" else "multinomial")
      }
      if(length(unique(y[!is.na(y)]))==2) return("binomial")
      vals = unique(y)
      if(sum(abs(vals-round(vals)))<1e-9) return("ordinal")
      return("gaussian")
    }))
  }else{
   
  return(apply(y_mat,2,function(y){
    if(!is.numeric(y)){
      y = as.factor(y);
    }
    if(is.factor(y)){
       return(if(length(levels(y))==2 )"binomial" else "multinomial")
    }
    if(length(unique(y[!is.na(y)]))==2) return("binomial")
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
  initialize=function(
              datasets,
              ys = lapply(datasets, function(d) d$y),
              flags = list(),
              mats = lapply(datasets, function(d) lapply(d$data, function(d1).getSparseMatrices(d1))),
              families=lapply(ys, function(d) .getFamily(d)),
              
                      memDir=NULL){
    datas = lapply(names(mats), function(nme){
      dataObj$new(mats[[nme]], incl_full=T,seed = getOption("seed",42), memDir=if(is.null(memDir)) NULL else paste(memDir, nme,sep="/"))
    })
    names(datas) = names(mats)
    self$type="slow1"
    types_all = getOption("types_all",names(datas[[1]]$data))
    names(types_all) = types_all
    var_threshs=  lapply(types_all, function(v) .readFlag(flags, "var_thresh",0.00))
    genes_incls=.readFlag(flags,"genes_incls",NULL) #,getOption("genes_incls",NULL)
    batch=.readFlag(flags, "batch",0)
    
    nrep = .readFlag(flags,"nrep",if(batch>0) 0 else 1)
    invisible(lapply(datas, function(data) data$init1(var_threshs, nrep=nrep,  batch = batch,genes_incls=genes_incls)))
    varn = getOption("varn",c())
   
    invisible(lapply(1:length(datas), function(ik) {
      y1 = ys[[ik]] #dists[[ik]]$updateYdb(cats[['cats']])
      family = families[[ik]]
      datas[[ik]]$updateY(y1, family=family, CHECK=T)
  #    datas[[ik]]$updateYdb(dists[[ik]]$mydb, cats[['cats']])
  #    missing_vals = self$updateY(y1, family=family, CHECK=T)
      
      datas[[ik]]$initTrain(varn=varn)
    #data[[ik]]$initY()
     return(NULL)
    }))
    self$datas = datas
  },
  
  update=function( flags = list()){
    types_all = getOption("types_all",names(datas[[1]]$data))
    names(types_all) = types_all
    var_threshs=  lapply(types_all, function(v) .readFlag(flags, "var_thresh",0.00))
    genes_incls=.readFlag(flags,"genes_incls",NULL) #,getOption("genes_incls",NULL)
    batch=.readFlag(flags, "batch",0)
    
    nrep = .readFlag(flags,"nrep",if(batch>0) 0 else 1)
    invisible(lapply(datas, function(data) data$init1(var_threshs, nrep=nrep,  batch = batch,genes_incls=genes_incls)))
    varn = getOption("varn",c())
    for(i in 1:length(self$data)){
      datas[[ik]]$initTrain(varn=varn)
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
  pheno=function(maxpheno=1e9){
   self$datas[[1]]$pheno(maxpheno=maxpheno);
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
      comb=.combineAngles(angles,topn=topn)
      as.list(comb)       
    })
    angles_all
  },
  makeAllModels=function(variables, phens, flags, all_models = list()){
    if(length(variables)==0) return(list())
    variables =variables[order(unlist(lapply(variables, function(v) length(v$var))),decreasing=T)]
  
    verbose=.readFlag(flags,"verbose",F)
    for(v_nme in names(variables)){
      if(verbose)print(v_nme)
      v = variables[[v_nme]]
      inds =v$inds
      vars1 = v$var
      nme_ = paste(names(vars1),collapse=";")
      models1 = all_models[[nme_]]
      if(is.null(models1)){
        models1 = self$makeModels( vars1, inds,phens,flags)
        for(k in 1:length(models1)){
          all_models[[names(models1)[[k]]]] = models1[[k]]
        }
      }else{
        new_inds = inds[which(!(inds %in% names(models1)))]
        if(length(new_inds)>0){
          models11 = self$makeModels( vars1, new_inds,phens,flags)
     
         for(k in 1:length(models11)){
           all_models[[names(models11)[[k]]]]= c(all_models[[names(models11)[[k]]]], models11[[k]])
          }
        }
      }
    }
    all_models
  },
  makeModels=function(vars1, inds, phens,flags){
    datas=self$datas
    train_nme = .readFlag(flags,'train', names(datas))
    if(length(which(train_nme %in% names(self$datas)))==0)train_nme = names(self$datas)[[1]]
    verbose=.readFlag(flags,"verbose",FALSE)
    #k = .readFlag(flags, 'rep',length(datas[[1]]$train))
    models = lapply(inds, function(k){
     # fold = (k<length(datas[[1]]$train))
     # lapply(var_eg, function(vars1){
        lapply(datas[names(datas) %in% train_nme], function(d){
         # fold_inds = if(fold) k else 1:length(d$train)
          mods = d$makeModels(phens,vars1, k,verbose)
          mods
        })
      #})
    })
    vars = names(models[[1]][[1]])
    names(vars) = vars
    models1 = lapply(vars, function(v){
      m3 = lapply(models, function(m){
        m2 = lapply(m, function(m1) m1[[v]])
        m2[unlist(lapply(m2, length))>0]
      })
      m3[unlist(lapply(m3, length))>0]
    })
    models1
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
  select=function(phens,flags,verbose=F){#, all_reps=F ## used to be train
    datas1 = self$datas
    topn = .readFlag(flags,'topn', 20)
 #   return_type = .readFlag(flags,'return','model') ##model,eval,plot
    train_nme = .readFlag(flags,'train', names(datas1))
    if(!(train_nme %in% names(datas1))) train_nme = names(datas1)[[1]]
    test_nme = .readFlag(flags,'test', names(datas1))
    maxsize=.readFlag(flags,'max',50)
    num_pvals = min(topn, 10)
    #type='slow1'  #.readFlag(flags, 'type', 'slow1') 
    incl = .readFlag(flags,'data_types',names(datas1[[1]]$data))
    names(incl )= incl
    logpvthresh = log(.readFlag(flags,"pthresh",1e-5))
    vars_l = list(mStateObj$new(c(),c(), NULL))
    logpv=-100
    
 #   k = .readFlag(flags, 'rep',length(datas1[[1]]$train))
    nreps1 =unlist(lapply(datas1, function(x) ncol(x$looc$incl)))
    if(length(table(nreps1))>1) stop(" must have same number of reps in each dataset, use nrep instead of bach")
    nreps = 1:nreps1[[1]]
   names(nreps) = nreps
   names(nreps)[length(nreps)]="full"
   #if(!all_reps) 
    beam = .readFlag(flags,"beam",2)
   
    
    variables =lapply(nreps, function(k){
      print(paste("cv",k,"of",length(nreps)))
      jj1=0
      invisible(lapply(train_nme, function(data_nme) datas1[[data_nme]]$update(k))); ### update training object
    while(logpv<logpvthresh && length(vars_l[[1]]$var)<maxsize){
      angles_all = lapply(vars_l, function(vars){
        angles=lapply(train_nme, function(data_nme) datas1[[data_nme]]$getAngles1(phens,vars$var,incl=incl,k=k, type=self$type))
        comb=.summariseAngles(.combineAngles(angles,topn=topn),topn)
        nxt_vars = lapply(1:num_pvals, function(ik){
          #print(ik)
          pv =  .getPvsAll(phens,datas1[names(datas1) %in% train_nme], c(vars$var,comb[ik]),k)
          mStateObj$new(comb[ik],  .sumChisq(pv) , prev_i=vars)
        })
        names(nxt_vars) = names(comb)[1:length(nxt_vars)]  
        nxt_vars[order(unlist(lapply(nxt_vars, function(nv)nv$cumpv)))]          
      })
      angles_all1 = unlist(angles_all, rec=F)
      ord = order(unlist(lapply(angles_all1, function(nv)nv$cumpv)))
      logpv = angles_all1[[ord[1]]]$logpv
      #print(logpv)
      if(logpv<=logpvthresh){
        vars_l = angles_all1[ord][1:beam]
      }
      if(verbose){
        nmes=names(vars_l[[1]]$var)
       print(nmes[length(nmes)])
         print(paste("logpv",logpv,jj1))
         jj1 = jj1+1
      }
    }
    #lapply(datas1, function(d)d$saveParquet())
    #lapply(vars_l, function(v) v$var)
      ##just take the top
    vars_l[[1]]$var
    })
    
    lens = unlist(lapply(variables, function(v) length(names(v))))
    if(max(lens)==0) return(list())
    sizes = 1:max(lens)
    names(sizes)=sizes
    vars1 = lapply(sizes, function(j){
      indsj = which(lens>=j)
      variables_j = variables[indsj]
      var_names =unlist(lapply(variables_j, function(v) paste(names(v)[1:j],collapse=";"))) ## just take the top one for now
      vars = names(sort(table(unlist(var_names)),decr=T))
      names(vars)=vars
      vars1 = lapply(vars, function(v) {
        inds=as.list(indsj[which(var_names==v)])
        list(inds = inds ,var = variables[[inds[[1]]]][1:j])
      })
    })
    names(vars1) = NULL
    vars2 = unlist(vars1, rec=F)
    flags_out = list(train=train_nme,  max=maxsize, 
                     nreps = nreps, beam=beam,
                     topn = num_pvals, pthresh =exp(logpvthresh), data_types = incl)
    attr(vars2,"flags_out")=toJSON(flags_out, simplifyVector=T, flatten=T)
    vars2
  },
  extractPredictions=function(all_models,phens, flags, CV = FALSE){
    #datas = self$datas
    res2=lapply(self$datas, function(d){
      d$extractPredictions(all_models, phens, falgs, CV=CV)
    })
   
    res2[lapply(res2,length)>0]
 },
  evaluateAllModels=function(all_models, phens,flags){ ## different folds with same variables
##                          inds = as.numeric(names(all_models))){
    if(length(all_models)==0) return(NULL)
    full_models = lapply(all_models, function(all_model1) all_model1[['full']])
    full_models = full_models[unlist(lapply(full_models, length))>0]
    full_model_nmes = names(full_models)
    nme_d = .readFlag(flags,"test",names(self$datas))
    names(nme_d) = nme_d
    print(nme_d)
    eval1=.merge1_new(lapply(nme_d, function(nme1){
     d = self$datas[[nme1]]
      resd = try(d$evaluateAllModels(all_models,phens,flags))
      if(inherits(resd,"try-error")) {
        print(paste("problem", nme1))
        return(NULL)
      }
      resd
    }),addName="data")
 
    if(is.null(eval1)) return(NULL)
    eval2 = eval1%>% pivot_wider(names_from="submeasure")
    isfull=eval2$model %in% full_model_nmes
    eval2%>%tibble::add_column(isfull=isfull)
  },
  pvalues=function(vars,phens,flags){
    datas = self$datas
    k = .readFlag(flags, 'rep',length(datas[[1]]$train))
    lapply(vars,function(vars1) {
      lapply(datas, function(d){
        b_i1 = vars1[[length(vars1)]]
        prev_var =  vars1[-length(vars1)]
        
        sig_res = d$calcBetaProj1(phens,k, b_i1, prev_var, convert=T )
        unlist(sig_res$pvs)
        #d$train[[k]]$getPvs(prev)
      }) 
    })
    #.getPvsAll(phens,datas,vars1,k))
  }
))
