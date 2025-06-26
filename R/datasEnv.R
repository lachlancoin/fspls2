
.getPvsAll<-function(subphens,datas, var_new,k, transform_y){
  transform_func = eval(str2lang(transform_y))
  pvs_all= lapply(datas, function(d){
    b_i1 = var_new[[length(var_new)]]
    prev_var =  var_new[-length(var_new)]
    sig_res = d$calcBetaProj1(subphens,k, b_i1, prev_var , transform_func,convert=T)
    unlist(sig_res$pvs)
    #d$train[[k]]$getPvs(prev)
  })
  pvs_all
}
.getFamily<-function(y_mat, family1=NULL){
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
    var_threshs=  lapply(types_all, function(v) .readFlag(flags, "var_quantile",0.00))
    genes_incls=.readFlag(flags,"genes_incls",NULL) #,getOption("genes_incls",NULL)
    batch=.readFlag(flags, "batch",0)
    all_v_all = .readFlag(flags,"all_v_all",T)
    nrep = .readFlag(flags,"nrep",if(batch>0) 0 else 1)
    invisible(lapply(datas, function(data) data$init1(var_threshs, nrep=nrep,  batch = batch,genes_incls=genes_incls)))
    varn = getOption("varn",c())
   
    invisible(lapply(1:length(datas), function(ik) {
      y1 = ys[[ik]] #dists[[ik]]$updateYdb(cats[['cats']])
      family = families[[ik]]
      datas[[ik]]$updateY(y1, family=family, CHECK=T, all_v_all=all_v_all)
  #    datas[[ik]]$updateYdb(dists[[ik]]$mydb, cats[['cats']])
  #    missing_vals = self$updateY(y1, family=family, CHECK=T)
      
      datas[[ik]]$initTrain(varn=varn)
      datas[[ik]]$initTrain1();
    #data[[ik]]$initY()
     return(NULL)
    }))
    self$datas = datas
  },
  
  update=function( flags = list()){
    types_all = getOption("types_all",names(self$datas[[1]]$data))
    names(types_all) = types_all
    var_threshs=  lapply(types_all, function(v) .readFlag(flags, "var_quantile",0.00))
    genes_incls=.readFlag(flags,"genes_incls",NULL) #,getOption("genes_incls",NULL)
    batch=.readFlag(flags, "batch",0)
    
    nrep = .readFlag(flags,"nrep",if(batch>0) 0 else 1)
    invisible(lapply(self$datas, function(data) data$init1(var_threshs, nrep=nrep,  batch = batch,genes_incls=genes_incls)))
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
  makeAllModels=function(vars_all, phens, flags){
    verbose=.readFlag(flags,"verbose",F)
    nmes_vars_all = names(vars_all); names(nmes_vars_all) = nmes_vars_all
    func_strs = fromJSON(.readFlag(flags,"transform_y",'{"y":"function(y) y"}'))
    
    lapply(nmes_vars_all, function(nme_v_all){
      vars = vars_all[[nme_v_all]]
      func_str = func_strs[[nme_v_all]]
      all_models = list()
      variables = vars$variables
      var_inds = vars$inds
      if(length(variables)==0) return(list())
      ord = order(unlist(lapply(variables, length)),decreasing=T)
      variables = variables[ord]
      var_inds = var_inds[ord]
    
     for(v_nme in names(variables)){
       if(verbose)print(v_nme)
       vars2 = variables[[v_nme]]
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
  },
  makeModels=function(vars2, inds, phens,func_str, flags){
    datas=self$datas
    train_nme = .readFlag(flags,'train', names(datas)[1])
    if(length(which(train_nme %in% names(self$datas)))==0)train_nme = names(self$datas)[[1]]
    verbose=.readFlag(flags,"verbose",FALSE)
    #k = .readFlag(flags, 'rep',length(datas[[1]]$train))
    nmes_inds = names(inds);names(nmes_inds)=nmes_inds
    
    models = lapply(nmes_inds, function(nmes_inds1){
      inds1 = inds[[nmes_inds1]]
      phens1 = phens[[nmes_inds1]]
        mods1 = lapply(inds1, function(k){
            lapply(datas[names(datas) %in% train_nme], function(d){
              mods = d$makeModels(phens1, vars2,k,func_str, verbose)
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
        var2 = var1[[phenn]]$var
        cumpv = var1[[phenn]]$cumpv
        varn = paste(names(var2), collapse=";")
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
   self$convert1(variables, phens)
  },
  select=function(phens,flags,verbose=F ){
    nreps1 =ncol(self$datas[[1]]$looc$incl)
    datas1 = self$datas
    topn = .readFlag(flags,'topn', 20)
    train_nme = .readFlag(flags,'train', names(datas1))
    train_nme = train_nme[train_nme %in% names(datas1)]
    if(length(train_nme)==0) train_nme = names(datas1)[[1]]
    maxsize=.readFlag(flags,'max',50)
    num_pvals = min(topn, 10)
    incl = .readFlag(flags,'data_types',names(datas1[[1]]$data))
    names(incl )= incl
    logpvthresh = log(.readFlag(flags,"pthresh",1e-5))
    logpv=-100
    nreps = 1:nreps1
    names(nreps) = nreps
    beam = .readFlag(flags,"beam",1)
    func_str = fromJSON(.readFlag(flags,"transform_y",'{"y":"function(y) y"}'))
    flags_out = list(train=train_nme,  max=maxsize, 
                     nreps = nreps, beam=beam,func_str = func_str,
                     topn = num_pvals, pthresh =exp(logpvthresh), data_types = incl)
  
    vars_combined=lapply(func_str,function(funcst){
      print(funcst)
     variables=lapply(nreps, function(k){
      print(paste("cv",k,"of",length(nreps)))
      jj1=0
      invisible(lapply(train_nme, function(data_nme) datas1[[data_nme]]$update(k, funcst))); ### update training object
      phens_index = 1:length(phens)
      names(phens_index) = names(phens)
      res2=  lapply(phens_index, function(p_index){
      subphens = phens[[p_index]]
        if(FALSE) cat(p_index); cat("\t");
        vars_l = list(mStateObj$new(c(),c(), NULL))  ## initialise vars_l
        while(logpv<logpvthresh && length(vars_l[[1]]$var)<maxsize){
          angles_all = lapply(vars_l, function(vars){
            angles=lapply(train_nme, function(data_nme) datas1[[data_nme]]$getAngles1(subphens,vars$var,incl=incl,k=k, type=self$type))
            names(angles) = train_nme
            cols_incl = lapply(train_nme, function(data_nme)datas1[[data_nme]]$cols_incl)
             comb1=.combineAngles(angles,cols_incl,topn=topn)
            comb=.summariseAngles(comb1,topn)
            num_pvals1 = min(num_pvals, length(comb))
            nxt_vars = lapply(1:num_pvals1, function(ik){
           #   print(ik)
              var_new = c(vars$var,comb[ik])
              pv =  .getPvsAll(subphens,datas1[names(datas1) %in% train_nme], var_new,k, funcst)
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
       # if(length(vars_l[[1]]$var)>0) print(vars_l)
        vars_l[[1]]#$var
        })
         res2[unlist(lapply(res2,function(xx) length(xx$var)))>0]
    })
    self$post_process(variables,flags_out)
  })
  
 
   vars_combined
  },
  extractPredictions=function(all_models,phens, flags, CV = FALSE, liab=T, data_nme  = names(self$datas)){
    #datas = self$datas
    names(data_nme) = data_nme
    res3 = lapply(all_models, function(all_models_){
      nmes_p = names(phens); names(nmes_p) = nmes_p
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
  evaluateAllModels=function(all_models, phens,flags){ ## different folds with same variables
##                          inds = as.numeric(names(all_models))){
    #self = all_models
    if(length(all_models)==0) return(NULL)
    nme_d = .readFlag(flags,"test",names(self$datas))
    names(nme_d) = nme_d
    print(nme_d)
    func_strs = fromJSON(.readFlag(flags,"transform_y",'{"y":"function(y) y"}'))
    inverse_func_strs = fromJSON(.readFlag(flags,"transform_y_inverse",'{"y":"function(y) y"}'))
    for(k in 1:length(func_strs)){ #CHECK INVERSE CORRECT
      func0 = eval(str2lang(func_strs[[k]]))  ## should be inverse
     func1 = eval(str2lang(inverse_func_strs[[k]]))  ## should be inverse
     xx = -10:10  
     if(max(abs(apply(cbind(xx,func1(func0(xx))),1,diff)),na.rm=T)>1e-5) stop("please provide inverse function")
    }
    
    mod_nmes = names(all_models); names(mod_nmes) = mod_nmes
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
    eval2 = eval1%>% pivot_wider(names_from="submeasure")
  #  isfull=eval2$model %in% full_model_nmes
  #  eval2%>%tibble::add_column(isfull=isfull)
    eval2
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
