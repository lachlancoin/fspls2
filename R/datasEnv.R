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

.nonZero<-function(am2){
  am2[unlist(lapply(am2, length))>0]
}

##if(length(cohort)>1){
#       
#        lapply(cohort, function(c){
#          mi1 =  match(rownames(c$matrix), rn)
#          mi0 =  match(rn,rownames(c$matrix))
#          
#         
#          if(typeof(c$matrix)!="S4") stop(" should be sparse matrix")
#          
#          
#          })
#    }


datasEnv<-R6Class("datasEnv", public = list(
  datasH = "list",
  type="character",
  sigsdir="character",
  sigs="environment",
  flags ="list",
  transform_y="character",
  initialize=function(
              datasH,
              flags = list(),
              dbDir="./",
                      memDir=NULL){
    self$flags = flags
    self$datasH = datasH
    transform_y =.readFlag(flags, "transform_y",toJSON(list(x=list(unvfunc="function(y,param) y",func="function(y,param) y", param=1))))
    ### MAKE SIGNATURE DIRECTORY
    self$sigsdir=paste(dbDir,"fspls_signatures",sep="/")
    self$datasH = datasH 
    self$sigs=   sigEnv$new(self$sigsdir,nme1)
    self$transform_y= self$sigs$updateData(data_flags = flags, 
                                 data_names =names(datasH), 
                                 data_types = names(datasH[[1]]$data$data),
                                 dims = self$dims(),
                                 transform_y = transform_y)
  },
 
  clear_db=function(drop=F,exclude="vars", recursive=F){
    if(drop){
      self$sigs$drop_all(exclude=exclude)
      if(recursive){
        lapply(self$datasH, function(dh) dh$clear_db(drop=drop, exclude=exclude))
      }
    }else{
      warning("need to set drop=T if you are sure, this will delete all saved signatures")
    }
    
  },  
  
  dims=function(){
    lapply(self$datasH, function(data) data$dims())
  },
 cats = function(maxpheno = 1e9){
    lapply(self$datasH, function(d) d$cats(maxpheno))
 },
  pheno=function(maxpheno=1e9,sep=F, sep_group = F, code=NULL, memb=NULL){
   res = self$datasH[[1]]$pheno(maxpheno=maxpheno, sep=sep, sep_group = sep_group, code = code,memb=memb);
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
      angles=lapply(train_nme, function(data_nme) datas[[data_nme]]$getAngles1(phens,vars1,incl=incl,k=k, type=self$type, 
                                                                              ))
      cols_incl = lapply(train_nme, function(data_nme)datas1[[data_nme]]$cols_incl)
      comb=.combineAngles(angles,cols_incl,topn=topn)
      as.list(comb)       
    })
    angles_all
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
updateTrain=function( phens, flags,  verbose=F){
  train_nme = .readFlag(flags,"train", names(self$datas))
  mi = match(train_nme,names(self$datas))
  if(length(which(is.na(mi)))>0) stop("train flag is wrong")
  datas1 = self$datas
  invisible(lapply(train_nme, function(data_nme) {
    if(verbose)print(data_nme)
    datas1[[data_nme]]$updateTrain( phens,flags,verbose=verbose)
  })); ### update training object - updates all
},
updateLOOC=function( phens, flags,varn=c(),force=F, verbose=F){
  train_nme = names(self$datas) # .readFlag(flags,"train", names(self$datas))
  datas1 = self$datas
  invisible(lapply(train_nme, function(data_nme) datas1[[data_nme]]$updateLOOC( phens,flags,varn=varn,force=force, verbose=verbose))); ### update training object - updates all
  
},

  select=function(phens,flags, verbose=F, db=NULL,user="",useDB=T ){#c(y="function(y) y","function(y) y")
   
    sigDB = if(useDB) self$sigs else NULL
    
    if(!is.null(sigDB) ){
      vars_all = sigDB$loadVars(flags, phens)
      if(!is.null(vars_all)) return(vars_all)
      vars_all$db=db
    }
    for(dh in self$datasH){
      dh$updateLOOC(phens, flags, verbose=verbose)
      dh$updateTrain(phens, flags,verbose=verbose)
    }
    datas1 = self$datasH
    nreps1 =self$datasH[[1]]$nreps()
    project=.readFlag(flags,"project",T)
    useoffset=.readFlag(flags,"useoffset",T)
    topn = .readFlag(flags,'topn', 20)
    useglm = .readFlag(flags,'useglmnet',T)
    onlyAll = .readFlag(flags,'only_all',F)
    train_nme = .readFlag(flags,'train', names(datas1))
    train_nme = train_nme[train_nme %in% names(datas1)]
    quantiles = sort(fromJSON(.readFlag(flags, "quantiles","[0]")),decreasing=T)
    genes_incls=fromJSON(.readFlag(flags,"genes_incls",'{"all":["all"]}')) #,getOption("genes_incls",NULL)
    if(!is.list(genes_incls)) stop("genes incls should be list")
    if(length(train_nme)==0) train_nme = names(datas1)[[1]]
    names(train_nme) = train_nme
    maxsize=.readFlag(flags,'max',500)
    minsize=.readFlag(flags,'min',0)
    num_pvals = min(topn, 10)
    incls = fromJSON(.readFlag(flags,'data_types',"{}")) 
    if(length(incls) == 0 )incls = list("all"=names(self$datasH[[1]]$data$data))
    incls_all = unique(unlist(incls))
    logpvthresh = log(.readFlag(flags,"pthresh",1e-5))
    logpv=-100
    nreps = 1:nreps1
    names(nreps) = nreps
    beam = .readFlag(flags,"beam",1)
    stop_y = .readFlag(flags, 'stop_y',"rand")
    stop_random=F
   # funcst = transform_y[[1]]
   # phens_index = 1:length(phens)
   # names(phens_index) = names(phens)
    var_thresh = lapply(train_nme, function(data_nme){
      lapply(self$datasH[[data_nme]]$data$vars, function(v) quantile(v, quantiles))
    })
    Wall0 =lapply(phens, function(f) matrix(nrow=0,ncol=0)) 
     #k=1;  qq =1; incl = incls[[1]]; data_nme = train_nme[[1]];g_incl  = genes_incls[[1]];#nme_c1 = names(transform_y)[[1]]
     variables=lapply(nreps, function(k){
      if(verbose) print(paste("cv",k,"of",length(nreps)))
      jj1=0
      func_ind=c()
      #invisible(lapply(train_nme, function(data_nme) datas1[[data_nme]]$update(k, funcst, phens,incls_all))); ### update training object
     # res2=  lapply(phens_index, function(p_index){
        if(verbose) print(paste(k,length(nreps)))
       
#        invisible(lapply(train_nme, function(data_nme) datas1[[data_nme]]$update(k, funcst, subphens,incls_all))); ### update training object
        if(FALSE) cat(p_index); cat("\t");
      
        vars_l =list(lapply(train_nme,function(xx) stateObj$new(phens, NULL,NULL,NULL,NULL,k, var=c(), varnames=c(), W_all =Wall0)))
      names(vars_l) = "empty"
      #vars_l1 = vars_l[[1]]
        for(incl in incls){
         if(verbose) print(incl)
         for(g_incl in genes_incls){
         for(qq in 1:length(quantiles)){
           
         while( (length(vars_l[[1]][[1]]$var) < minsize || logpv<logpvthresh) && length(vars_l[[1]][[1]]$var)<maxsize && ! stop_random){
          angles_all = lapply(vars_l, function(vars_l1){
            nxt_vars1 =  tryCatch({
              varnames = vars_l1[[1]]$var_names; type = self$datasH[[1]]$type
              comb_all1=lapply(train_nme, function(data_nme) {
              #print(data_nme); 
                angleH=list(angles=
                  self$datasH[[data_nme]]$data$getAngles1(phens,varnames,incl=incl,k=k, type=type),
                    cols_incl = self$datasH[[data_nme]]$data$cols_incl(var_thresh[[data_nme]],incl, g_incl,qq, excl=varnames)) ### fix 
                  .combineAngles1(angleH, incl, topn=topn, onlyAll=onlyAll, excl=varnames)
              })
            if(length(comb_all1)>1) stop("need to write a merge function for comb across data")
             comb_ =  comb_all1[[1]]
             #comb_=.combineAngles(anglesH,incl,topn=topn, onlyAll = onlyAll, excl=varnames)
             nme_comb = names(comb_); names(nme_comb) = nme_comb
             Wall = vars_l1[[1]]$W_all
            res_inner=lapply(nme_comb, function(nme_c1){
               nmesp1 = names(comb_[[nme_c1]]); names(nmesp1) = nmesp1
               lapply(nmesp1, function(nme_p1){
                comb = comb_[[nme_c1]][[nme_p1]]
                if(nrow(comb)==0) return(NULL)
                num_pvals1 = min(num_pvals, nrow(comb))
                inds1p = 1:num_pvals1; names(inds1p) = comb$names[1:length(inds1p)]
                nxt_vars = lapply(inds1p, function(ik){
                   b_i_name = c(comb$data_type[[ik]], comb$names[[ik]], nme_c1,nme_p1)
                   nv = try(lapply(train_nme, function(tn){
                     dh = self$datasH[[tn]]
                     prev_i = vars_l1[[tn]]
                     dh$getPvsAll(phens,prev_i, b_i_name,k,  vars_l1[[tn]]$W_all,project = project, 
                                  useglm=useglm,
                                  inv_transform=.readFlag(flags,"x_transform",T),
                                  useoffset=useoffset)
                   }))
                   if(inherits(nv,"try-error")) {
                     print(paste(nme_c1, "error"))
                     return(NULL)
                   }
                  attr(nv,"cumpv")= .sumChisq(unlist(lapply(nv, function(nv1){
                     unlist(nv1$pvs)
                   })))
                  attr(nv,"cumpv_all")= .sumChisq(unlist(lapply(nv, function(nv1){
                    unlist(nv1$pvs_all)
                  })))
                  
                  nv
                  #mStateObj$new(comb[ik],  .sumChisq(pv) , prev_i=prev_i)
                })
                nxt_vars = nxt_vars[unlist(lapply(nxt_vars, length))>0]
                if(length(nxt_vars)==0) return(NULL)
                pvs_list = unlist(lapply(nxt_vars, function(nv) attr(nv,"cumpv")))
               # pvs_list_all = unlist(lapply(nxt_vars, function(nv) attr(nv,"cumpv_all")))
                #print(pvs_list)
                #subinds1 = pvs_list<=logpvthresh
                #if(length(which(subinds1)==0)) return(NULL)
                #nxt_vars = nxt_vars[subinds1]
                #pvs_list = pvs_list[subinds1]
                nxt_vars[order(pvs_list)] 
            })
            })
            res_inner = res_inner[unlist(lapply(res_inner, length))>0]
            if(length(res_inner)==0) return(NULL)
            res_inner
            },error=function(w){
              print(w)
              #print("error")
              return(NULL)
            })
            nxt_vars1 =nxt_vars1[unlist(lapply(nxt_vars1, length))>0]
            if(length(nxt_vars1)==0) return(NULL)
            nxt_vars1
            
          })
          angles_all = angles_all[unlist(lapply(angles_all, length))>0]
          if(length(angles_all)==0){
            stop_random=T;
            next;
          }
          
          ang1 = unlist(unlist(unlist(angles_all, rec=F),rec=F),rec=F)
          logpvs = unlist(lapply(ang1, function(a1)attr(a1,"cumpv")))
          logpvs_all = unlist(lapply(ang1, function(a1)attr(a1,"cumpv_all")))
          ord = order(logpvs)
          names(ord) = names(logpvs)
          ord_all = order(logpvs_all)
          ang1 = ang1[ord_all]
          logpvs = logpvs[ord_all]
          logpvs_all = logpvs_all[ord_all]
          if(!.readFlag(flags,"x_transform",T)){
            if(length(func_ind)>0 ){
              finds = unlist(lapply(ang1, function(a1)a1[[1]]$var[[1]][3]))
              ang1 = ang1[finds==func_ind[[1]]]
            }else{
              func_ind = c(func_ind,ang1[[1]][[1]]$var[[1]][3])
            }
          }
          if(!is.null(stop_y)){
            gp1=grep(stop_y, names(logpvs))
            gp=grep(stop_y, names(logpvs), inv=T)
         
           # print("HERE")
            #print(unlist(list(rand= min(logpvs[gp1]),nonrand=min(logpvs[gp]))))
            stop_random= min(logpvs[gp1])<=min(logpvs[gp])
            #print(head(sort(ord[stop_ind])))
          }
         
     
          #print("HERE ALL")
          if(stop_random){
            if(verbose) print(paste("stopping due to random", exp(logpv), names(logpvs)[which.min(logpvs)]))
          }
          logpv =min(logpvs)
          
         
          #logpv<=logpvthresh || length(vars_l[[1]][[1]]$var) < minsize 
          if((!stop_random && logpv<=logpvthresh) || length(vars_l[[1]][[1]]$var)<minsize  ){
            dupls=(unlist(lapply(ang1, function(a1) paste(unlist(lapply(a1[[1]]$var_names, function(vv1)paste(vv1[1:2],collapse="::"))), collapse=";;"))))
            ang1 = ang1[!duplicated(dupls)]
            vars_l = ang1[1:min(length(ang1),beam)]
          }
          if(verbose){
            print(head(sort(logpvs_all[gp])))
            print(names(vars_l))
            print(paste("logpv",logpv,min(logpvs_all), jj1))
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
    vars_all$flags = flags; vars_all$phens = phens;  vars_all$db=db #vars_all$transform_y = transform_y ;
    if(length(vars_all$variables)==0) return(vars_all)
    if(!is.null(sigDB) ){
#      saveVars = function(vars_all,flags,phens,transform_y,user="", replace=T){
      sigDB$saveVars(vars_all,replace=T)
    }
    vars_all
  },
  extractPredictions=function(all_models,phens=all_models$phens, flags=all_models$flags, CV = FALSE, liab=T, data_nme  = names(self$datas)){
    #datas = self$datas
  #  inverse_func_str =lapply(transform_y, function(t_y)  t_y[[2]])
    
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
          d$extractPredictions(all_models_, phens, flags, CV=CV, liab=liab)
        })
      #})
   
    res2[lapply(res2,length)>0]
    })
    predictions0=res3[unlist(lapply(res3, function(x) length(x[[1]][[1]])))>0]
    predictions0 # 
 },
updateTransforms = function(transform_y){
  self$flags[['transform_y']] = transform_y
  for(k in 1:length(self$datasH)){
    self$datasH[[k]]$updateTransforms(transform_y)
  }
},
  pvalues=function(vars,phens,transform_y,flags){
    stop("not working")
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
