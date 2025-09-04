.extractFullVars<-function(vars_all){
  subinds = which(unlist(lapply(vars_all$inds, function(x) length(grep('full',names(x)))))>0)
  list(variables = vars_all$variables[subinds], inds = vars_all$inds[subinds], cumpv = vars_all$cumpv[subinds],
       transf = vars_all$transf[subinds],
       flags = vars_all$flags,
       phens=vars_all$phens, transform_y = vars_all$transform_y)
}
.extractFullModels<-function(all_models, full_model_only=F){
  subinds = which(unlist(lapply(all_models$models, function(x) length(grep('full', names(x)))))>0)
  models1 = all_models$models[subinds]
  res1 =list(flags = all_models$flags, phens=all_models$phens, transform_y = all_models$transform_y)
  if(!full_model_only) {
    res1$models = models1
  }else{
  res1$models = lapply(models1, function(mod1){
    mod1[grep('full', names(mod1))]
  })
  }
  res1
}
.getPvsAll<-function(subphens,datas1, vars_l1, b_i_name,k,func_str1, project=T, useoffset=T){
  
  #var_new = c(prev_i$var,b_i)
  nmesd = names(datas1); names(nmesd)=nmesd
  #nmed = nmesd[[1]]
  pvs_all= lapply(nmesd, function(nmed){
    d = datas1[[nmed]]
    prev_i = vars_l1[[nmed]]
    family = strsplit(names(subphens)[[1]],"\\.")[[1]][1]
    if(family=="multinomial") useoffset=F
    prev_i1 = d$makeNextModel(prev_i,b_i_name,subphens,k, func_str1,family, ypred=NULL, project=project, useglm=T, logpthresh =0, useoffset=useoffset)
  #  if(inherits(prev_i1,"try-error")) {
    #  print(paste("problem", nmed))
   #   return(NULL)
  #  }
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
            mat3[which(is.infinite(mat3))]=NA
            na_ind = which(is.na(mat3))
            na_m = mats2$matrixNA
            if(length(na_ind)>0){
              na_m[na_ind] = 1
              mat31 = apply(mat3,2,function(v){ ## sets mean to col mean
                mv = mean(v, na.rm=T)
                v[is.na(v)]=mv
                v
              })
              mat3 = Matrix(mat31, sparse=TRUE)
             # na_ind1 = which(is.na(mat3), arr.ind=T)
              #na_cols = unique(na_ind1[,2])
            }
            list(matrix = mat3, matrixNA = na_m)
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
    pheno_balance=.readFlag(flags,"pheno_balance",NULL)
  #  invisible(lapply(datas, function(data) data$init1(pheno_balance = pheno_balance, nrep=nrep,  batch = batch)))
    varn = getOption("varn",c())
   
    invisible(lapply(1:length(datas), function(ik) {
      y1 = ys[[ik]] #dists[[ik]]$updateYdb(cats[['cats']])
      family = families[[ik]]
     
      ##need to work on all_v_all
      datas[[ik]]$updateY(y1, family=family, CHECK=T, all_v_all=all_v_all, one_v_rest = one_v_rest)
  #    datas[[ik]]$updateYdb(dists[[ik]]$mydb, cats[['cats']])
  #    missing_vals = self$updateY(y1, family=family, CHECK=T)
      
  #    datas[[ik]]$initTrain(varn=varn)
    #data[[ik]]$initY()
     return(NULL)
    }))
    self$datas = datas
  },
 plotData=function(vars_all1, phens1 = vars_all1$phens, all_types=F, transform_x = NULL, violin=F, assoc=F){
  df4= .merge1_new( lapply(self$datas, function(d) d$plotData(vars_all1, phens1 = phens1, all_types=all_types, transform_x = transform_x, violin=violin, assoc=assoc)), 
                addName="dataset")
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
   getSigDB=function(nme1="",reload= F, clear=F, user=""){
     if(is.null(nme1)) return(NULL)
   if(reload || is.null(self$sigs[[nme1]])){
     self$sigs[[nme1]]=   sigEnv$new(self$sigsdir,nme1, clear=clear, user=user)
     self$sigs[[nme1]]$updateData(data_flags = self$flags, 
                                  data_names =names(self$datas), 
                                  data_types = names(self$datas[[1]]$data),
                                  dims = self$dims(),
                                  phenos =self$datas[[1]]$pheno())
     #self$sigs[[nme1]]$data_id
   }
   self$sigs[[nme1]]
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
  pheno=function(maxpheno=1e9,sep=F, sep_group = F, exclude=NULL, code=NULL){
   res = self$datas[[1]]$pheno(maxpheno=maxpheno, sep=sep, sep_group = sep_group, code = code);
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
  makeAllModels=function(vars_all, phens=vars_all$phens, flags=vars_all$flags, transform_y=vars_all$transform_y,verbose=F, max = 1e6,
                         user="",
                         db=vars_all$db){
    if(typeof(transform_y)=="character"){
      if(length(transform_y)==1) transform_y = fromJSON(transform_y)   ## this is to fix problem with JSON
      transform_y = list(x=transform_y)
    }
    sigDB = self$getSigDB(db, user=user)
    if(!is.null(sigDB) ){
      all_models =try( sigDB$loadModels(flags,phens, transform_y=transform_y))
      if(inherits(all_models,"try-error")) {
        print(paste("problem reading from DB .. recalculating"))
      }else if(!is.null(all_models) && length(all_models$models)>0){
        return(all_models)
      }
    }
    transf=vars_all$transf
    #self$update(phens, flags, transform_y)
    self$updateLOOC(phens,flags)
   # nmes_vars_all = names(vars_all); names(nmes_vars_all) = nmes_vars_all
    #func_strs = fromJSON(.readFlag(flags,"transform_y",'{"y":"function(y) y"}'))
    #nme_v_all = nmes_vars_all[[1]]
    
    logpthresh= log(.readFlag(flags,"pthresh",1e-3))
    project=.readFlag(flags,"project",TRUE)
   # combined_models=lapply(nmes_vars_all, function(nme_v_all){
    #  if(verbose) print(nme_v_all)
      vars = vars_all#[[nme_v_all]]
      func_str = lapply(transform_y, function(xx) xx[[1]]) #func_strs[[nme_v_all]]
      all_models = list()
      variables = vars$variables
      var_inds = vars$inds
   
    #  print(var_inds);stop("!!")
      rem_inds = 1:self$datas[[1]]$nreps()
       names(rem_inds) = as.character(rem_inds)
       names(rem_inds)[which(rem_inds==self$datas[[1]]$nreps())]="full"
       #inds1 = 1:2
       all_models = self$makeModels(list(),names(func_str)[[1]],rem_inds , phens, func_str, flags)
       
      if(length(variables)==0) return(list(models=all_models, flags = flags, phens = phens, transform_y = transform_y, db=db))
      ord = order(unlist(lapply(variables, length)),decreasing=T)
      variables = variables[ord]
      var_inds = var_inds[ord]
     #v_nme = names(vars_all$variables)[1]; max=10; verbose=T; k=1;variables =vars_all$variables; 
      
     for(v_nme in names(variables)){
       print(v_nme)
       if(verbose)print(v_nme)
       vars2 = variables[[v_nme]]
       var_transf=strsplit(names(transf[[v_nme]])[[1]],"_")[[1]]
       if(verbose)print(var_transf)
       vars2 = vars2[1:min(length(vars2), max)]
       inds =var_inds[[v_nme]]
       nme_ = paste(names(vars2),collapse=";")
       models1 = all_models[[nme_]]
       if(is.null(models1)){
        models1 = self$makeModels( vars2, var_transf,inds,phens,func_str,flags)
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
            models2 = self$makeModels( vars2, var_transf,inds#[p_i]
                                       ,phens#[which(names(phens) %in% p_nme)]
                                       ,func_str,flags)
            for(k in 1:length(models2)){
              all_models[[names(models2)[[k]]]] = models2[[k]]#[[p_nme]]#[[p_nme]]
                     
            }
          }else{
            subinds= inds[which(is.na(match(names(inds), names(all_models[[nme_]]))))]
            if(length(subinds)>0){ ## missing inds
              models3 = self$makeModels( vars2, var_transf,subinds,phens#[which(names(phens) %in% p_nme)]
                                         ,func_str, flags)
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
    all_models_=list(models=all_models, flags = flags, phens = phens, transform_y = transform_y, db=db)
    #})
    if(!is.null(sigDB) ){
      sigDB$saveModels (all_models_)
    }
    #combined_models
   all_models_
  },
  makeModels=function(vars2, var_transf,inds, phens,func_str1, flags){
    datas=self$datas
#    print(inds)
    logpthresh= log(.readFlag(flags,"pthresh",1e-3))
    project=.readFlag(flags,"project",TRUE)
    useoffset=.readFlag(flags,"useoffset",TRUE)
    train_nme = .readFlag(flags,'train', names(datas)[1])
    if(length(which(train_nme %in% names(self$datas)))==0)train_nme = names(self$datas)[[1]]
    verbose=.readFlag(flags,"verbose",FALSE)
    useglm=.readFlag(flags,"useglmnet",TRUE)
  #  print(paste("useglm", useglm))
    #k = .readFlag(flags, 'rep',length(datas[[1]]$train))
   # nmes_inds = names(inds);names(nmes_inds)=nmes_inds
    #nmes_inds1 = nmes_inds[[1]];d =self$datas[[1]]
    #models = lapply(nmes_inds, function(nmes_inds1){
      inds1 = inds#[[nmes_inds1]]
      phens1 = phens#[[nmes_inds1]]
      #k=inds1[[1]]
        mods1 = lapply(inds1, function(k){
         # print(k)
            lapply(datas[names(datas) %in% train_nme], function(d){
              mods = d$makeModels(phens1, vars2,var_transf,k,logpthresh = logpthresh,project=project, func_str1 = func_str1,useglm=useglmnet, useoffset=useoffset)
          
            })
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
  post_process=function(variables){
  #  variables = variables[unlist(lapply(variables, function(x) length(x[[1]]$var)))>0]
    full_index = length(variables)
    #lens = unlist(lapply(variables, function(x) length(x$var)))
    #if(max(lens)==0) return(list())
    vars_all = list()
    vars_all1 = list()
    vars_all2 = list()
    vars_all3 = list() #funcstr
    names(variables) = 1:length(variables)
    func_inds = lapply(variables, function(vv) attr(vv,"func_ind"))
    for(repn in names(variables)){
      full = repn==full_index
      var1 = variables[[repn]]
      func_ind=attr(var1,"func_ind")
      func_str1 = paste(names(func_ind), collapse="_")
    #  for(phenn in names(var1)){
        var2 = var1[[1]]$var_names
        if(length(var2)>0){
            names(var2) = var1[[1]]$varnames
            cumpv = attr(var1,"cumpv")#lapply(var1, function(vv) attr(vv,"cumpv"))
            varn = paste(names(var2),collapse=";") #paste(names(var2), collapse=";")
            if(is.null(vars_all[[varn]])){
              vars_all[[varn]] = list()
              vars_all1[[varn]] = var2
              vars_all2[[varn]] = list()
              vars_all3[[varn]] =list(repn) 
              names(vars_all3[[varn]]) = func_str1
            }else{
              vars_all3[[varn]][[func_str1]]=c( vars_all3[[varn]][[func_str1]],repn)
            }
            repn1 = as.list(as.numeric(repn))
            repn2 = as.list(cumpv)
            names(repn1) = if(full) "full" else repn
            names(repn2) = if(full) "full" else repn
            if(is.null(vars_all[[varn]])){
              vars_all[[varn]] =repn1
              vars_all2[[varn]] =repn2
             # vars_all3[[varn]]
            }else{
              vars_all[[varn]] = c(vars_all[[varn]] , repn1)
              vars_all2[[varn]] = c(vars_all2[[varn]] , repn2)
            }
           
            
                       
      }
    }
    vars_combined = list(variables = vars_all1, inds = vars_all,cumpv=vars_all2 ,transf= vars_all3) 
  
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
updateTrain=function( phens, flags, transform_y, verbose=F){
  train_nme = .readFlag(flags,"train", names(self$datas))
  datas1 = self$datas
  invisible(lapply(train_nme, function(data_nme) {
    if(verbose)print(data_nme)
    datas1[[data_nme]]$updateTrain( phens,flags,transform_y, verbose=verbose)
  })); ### update training object - updates all
},
updateLOOC=function( phens, flags,varn=c(),force=F, verbose=F){
  train_nme = names(self$datas) # .readFlag(flags,"train", names(self$datas))
  datas1 = self$datas
  invisible(lapply(train_nme, function(data_nme) datas1[[data_nme]]$updateLOOC( phens,flags,varn=varn,force=force, verbose=verbose))); ### update training object - updates all
  
},

  select=function(phens1,flags,transform_y=list(x=c("function(y) y","function(y) y")), verbose=F, db=NULL,user="" ){#c(y="function(y) y","function(y) y")
    if(typeof(transform_y)=="character"){
      transform_y = list(x=transform_y)
    }
    
    sigDB = self$getSigDB(db,user=user)
    if(!is.null(sigDB) ){
      vars_all = sigDB$loadVars(flags, phens1, transform_y = transform_y)
    
      if(!is.null(vars_all)) return(vars_all)
      vars_all$db=db
    }
    self$updateLOOC(phens1, flags, verbose=verbose)
    self$updateTrain(phens1, flags, transform_y, verbose=verbose)
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
    minsize=.readFlag(flags,'min',0)
    num_pvals = min(topn, 10)
    incls = fromJSON(.readFlag(flags,'data_types',"{}")) 
    if(length(incls) == 0 )incls = list("all"=names(datas1[[1]]$data))
    incls_all = unique(unlist(incls))
    logpvthresh = log(.readFlag(flags,"pthresh",1e-5))
    logpv=-100
    nreps = 1:nreps1
    names(nreps) = nreps
    beam = .readFlag(flags,"beam",1)
    funcst1 =lapply(transform_y, function(xx) xx[1])
    stop_y = .readFlag(flags, 'stop_y',"rand")
    stop_random=F
   # funcst = transform_y[[1]]
   # phens_index = 1:length(phens)
   # names(phens_index) = names(phens)
    var_thresh = lapply(train_nme, function(data_nme){
      lapply(datas1[[data_nme]]$vars, function(v) quantile(v, quantiles))
    })
     #k=1; g_incl  = genes_incls[[1]]; qq =1; incl = incls[[1]]; data_nme = train_nme[[1]]
     variables=lapply(nreps, function(k){
      if(verbose) print(paste("cv",k,"of",length(nreps)))
      jj1=0
      func_ind=c()
      #invisible(lapply(train_nme, function(data_nme) datas1[[data_nme]]$update(k, funcst, phens,incls_all))); ### update training object
     # res2=  lapply(phens_index, function(p_index){
        if(verbose) print(paste(k,length(nreps)))
       
#        invisible(lapply(train_nme, function(data_nme) datas1[[data_nme]]$update(k, funcst, subphens,incls_all))); ### update training object
        if(FALSE) cat(p_index); cat("\t");
        vars_l =list(lapply(train_nme,function(xx) stateObj$new(phens1, NULL,NULL,NULL,NULL,k, var=c(), varnames=c(), W_all = NULL)))
        #vars_l1 = vars_l[[1]]
        for(incl in incls){
          if(verbose) print(incl)
        for(g_incl in genes_incls){
         for(qq in 1:length(quantiles)){
           
         while( (length(vars_l[[1]][[1]]$var) < minsize || logpv<logpvthresh) && length(vars_l[[1]][[1]]$var)<maxsize && ! stop_random){
          angles_all = lapply(vars_l, function(vars_l1){
            nxt_vars1 =  tryCatch({
              varnames = vars_l1[[1]]$var_names; type = self$type
            angles=lapply(train_nme, function(data_nme) {
              #print(data_nme); 
              datas1[[data_nme]]$getAngles1(phens1,varnames,incl=incl,k=k, type=type)
              })
            cols_incl = lapply(train_nme, function(data_nme)datas1[[data_nme]]$cols_incl(var_thresh[[data_nme]],incl, g_incl,qq)) ### fix 
             comb_=.combineAngles(angles,cols_incl,incl,topn=topn, onlyAll = onlyAll)
             nme_comb = names(comb_); names(nme_comb) = nme_comb
            res_inner=lapply(nme_comb, function(nme_c1){
               #print(nme_c1)
                comb = comb_[[nme_c1]]
                if(nrow(comb)==0) return(NULL)
                num_pvals1 = min(num_pvals, nrow(comb))
                inds1p = 1:num_pvals1; names(inds1p) = comb$names[1:length(inds1p)]
                nxt_vars = lapply(inds1p, function(ik){
                   b_i_name = c(comb$data_type[[ik]], comb$names[[ik]])
                   nv = .getPvsAll(phens1,datas1[names(datas1) %in% train_nme], vars_l1, b_i_name,k, funcst1[[nme_c1]], project = project, useoffset=useoffset)
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
            })
            res_inner
            },error=function(w){
              print(w)
              print("error")
              return(NULL)
            })
            nxt_vars1
            
          })
          nme_func = names(funcst1); names(nme_func)=nme_func
          angles_all2=lapply(nme_func, function(xx){
            angles_all_ = lapply(angles_all, function(angles_all1){
              angles_all1[[xx]]
            })
            angles_all_ = angles_all_[unlist(lapply(angles_all_,length))>0]
            angles_all1 = unlist(angles_all_, rec=F)
            ord = order(unlist(lapply(angles_all1, function(nv)attr(nv,"cumpv"))))
            logpv1 = attr(angles_all1[[ord[1]]],"cumpv")
            #print(logpv)
           # if(logpv<=logpvthresh || length(vars_l[[1]][[1]]$var) < minsize){
              vars_l_nxt = angles_all1[ord][1:beam]
            #}
            
            list(vars_l = vars_l_nxt, logpv = logpv1)
          })
          logpvs=unlist(lapply(angles_all2, function(xx) xx$logpv))
          min_ind=which.min(logpvs)
         
          logpv =min(logpvs)
          if(!is.null(stop_y)){
            stop_random = length(grep(stop_y,names(angles_all2)[[min_ind]]))>0
          }
          if(stop_random){
            print(exp(logpvs))
            print(paste("stopping due to random", exp(logpv)))
          }
          #logpv<=logpvthresh || length(vars_l[[1]][[1]]$var) < minsize 
          if((!stop_random && logpv<=logpvthresh) || length(vars_l[[1]][[1]]$var)<minsize  ){
            func_ind = c(func_ind,min_ind)
            
            vars_l =lapply(angles_all2[[min_ind]]$vars_l , function(aa){
              aa$transf = names(angles_all2)[[min_ind]]
              aa
            })
            
          }
          if(verbose){
            print(names(vars_l))
            print(paste("logpv",nme_func[min_ind],logpv,jj1))
            jj1 = jj1+1
          }
          #angles_all = angles_all[unlist(lapply(angles_all,length))>0]
         
        }
        }
        }
      }
       
        #lapply(datas1, function(d)d$saveParquet())
        #lapply(vars_l, function(v) v$var)
          ##just take the top
       # if(length(vars_l[[1]]$var)>0) print(vars_l)
        attr(vars_l[[1]],"func_ind")=func_ind#$var
        vars_l[[1]]
    #     res2[unlist(lapply(res2,function(xx) length(xx[[1]]$var)))>0]
    })
    vars_all=self$post_process(variables)
    vars_all$flags = flags; vars_all$phens = phens1; vars_all$transform_y = transform_y ; vars_all$db=db
    if(length(vars_all$variables)==0) return(vars_all)
    if(!is.null(sigDB) ){
#      saveVars = function(vars_all,flags,phens,transform_y,user="", replace=T){
      sigDB$saveVars(vars_all,replace=T)
    }
    vars_all
  },
  extractPredictions=function(all_models,phens=all_models$phens,transform_y = all_models$transform_y, flags=all_models$flags, CV = FALSE, liab=T, data_nme  = names(self$datas)){
    #datas = self$datas
    inverse_func_str =lapply(transform_y, function(t_y)  t_y[[2]])
    
    self$updateLOOC(phens, flags)
    names(data_nme) = data_nme
    #nmes_p = names(phens); names(nmes_p) = nmes_p
    all_models_y = all_models$models
    #all_models_ = all_models_y[[1]];d_nme = data_nme[[1]];# nme_p = nmes_p[[1]]; phens1 = phens[[nme_p]];
    
    res3 = lapply(all_models_y, function(all_models_){
      #res2 = lapply(nmes_p, function(nme_p){
       # phens1 = phens[[nme_p]]
        res2=lapply(data_nme, function(d_nme){
          d = self$datas[[d_nme]]
          if(is.null(d)){
            print(d_nme)
            stop("!!")
          }
          d$extractPredictions(all_models_, phens,inverse_func_str, flags, CV=CV, liab=liab)
        })
      #})
   
    res2[lapply(res2,length)>0]
    })
    predictions0=res3[unlist(lapply(res3, function(x) length(x[[1]][[1]])))>0]
    predictions0 # 
 },
  evaluateAllModels=function(all_models, phens=all_models$phens,flags=all_models$flags,verbose=F,
                             transform_y=all_models$transform_y, db=all_models$db, user=""){ ## different folds with same variables
##                          inds = as.numeric(names(all_models))){
    #self = all_models
    if(typeof(transform_y)=="character"){
      transform_y = list(x=transform_y)
    }
    sigDB = self$getSigDB(db,user=user)
    
    if(!is.null(sigDB) ){
      eval1 = sigDB$loadEval(flags,phens, transform_y)
      if(!is.null(eval1) && nrow(eval1)>0){
        return(eval1)
      }
    }
    self$updateLOOC(phens, flags)
    if(length(all_models$models)==0) return(NULL)
    nme_d2 = .readFlag(flags,"test",names(self$datas))
    names(nme_d2) = nme_d2
    print(nme_d2)
   # func_strs = fromJSON(.readFlag(flags,"transform_y",'{"y":"function(y) y"}'))
  #  inverse_func_strs = fromJSON(.readFlag(flags,"transform_y_inverse",'{"y":"function(y) y"}'))
    #k=1
      func0 = lapply(transform_y, function(t_y) eval(str2lang(t_y[[1]])))  ## should be inverse
     func1 = lapply(transform_y, function(t_y) eval(str2lang(t_y[[2]])))  ## should be inverse
     xx = -10:10  
     if(max(abs(apply(cbind(xx,func1[[1]](func0[[1]](xx))),1,diff)),na.rm=T)>1e-5) stop("please provide inverse function")
    
#    mod_nmes = names(all_models); names(mod_nmes) = mod_nmes
    #mod_nme = mod_nmes[[1]]; nme1 = nme_d[[1]] 
    #eval1=.merge1_new( lapply(mod_nmes, function(mod_nme){
      inverse_func_str =lapply(transform_y, function(t_y)  t_y[[2]])
      all_models_y = all_models$models#[[mod_nme]]
    eval1 =  .merge1_new(lapply(nme_d2, function(nme1){
      print(nme1)
     d = self$datas[[nme1]]
      resd = d$evaluateAllModels(all_models_y,phens,inverse_func_str = inverse_func_str,flags, verbose=verbose)
      #if(inherits(resd,"try-error")) {
      #  print(resd)
      #  print(paste("problem", nme1))
      #  stop("!!")
      #  return(NULL)
      #}
      resd
    }),addName="data")
    #}),addName="transform_y")
    if(is.null(eval1)) return(NULL)
  #  eval1 = subset(eval1, model!="avg")
    eval2 = eval1%>% pivot_wider(names_from="submeasure") %>% tibble::add_column(transform_y=strsplit(transform_y[[1]]," ")[[1]][2])
  #  isfull=eval2$model %in% full_model_nmes
  #  eval2%>%tibble::add_column(isfull=isfull)
    if(!is.null(sigDB) ){
      sigDB$saveEval(eval2, flags,phens, transform_y = transform_y)
      eval1 = sigDB$loadEval(flags,phens,transform_y = transform_y)
      return(eval1)
    }
    #print("HH")
    .calcEval1(eval2, rename=F)
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
