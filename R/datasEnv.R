.extractFullVars<-function(vars_all){
  subinds = which(unlist(lapply(vars_all$inds, function(x) length(grep('full',names(x)))))>0)
  list(variables = vars_all$variables[subinds], inds = vars_all$inds[subinds], cumpv = vars_all$cumpv[subinds],
       transf = vars_all$transf[subinds],
       flags = vars_all$flags,
       phens=vars_all$phens, transform_y = vars_all$transform_y)
}
.mergeComb<-function(comb_all1, flags){
   topn = .readFlag(flags,'topn', 20);num_pvals = min(topn, 20)
  if(length(comb_all1)>1) stop("need to write a merge function for comb across data.  This needs to restrict to num_pvals")
   lapply(comb_all1[[1]], function(c1){
     lapply(c1, function(p1){
       head(p1, num_pvals)
     })
   })
}
.mergeResInner<-function(res_inner1){
  nme_comb = names(res_inner1[[1]]);names(nme_comb) = nme_comb
  res_inner=lapply(nme_comb, function(nme_c1){
    #nme_c1 = nme_comb[[1]]; nmesp1 = names(comb_[[nme_c1]]); names(nmesp1) = nmesp1
    nmesp1 = names(res_inner1[[1]][[nme_c1]]);names(nmesp1) = nmesp1
    lapply(nmesp1, function(nme_p1){
    #  comb = comb_[[nme_c1]][[nme_p1]]
     # if(nrow(comb)==0) return(NULL)
#      nme_comb = names(comb_); names(nme_comb) = nme_comb
      varnames = unique(unlist(lapply(res_inner1, function(ri) names(ri[[nme_c1]][[nme_p1]]))))
     
      #varnames = comb$names; names(varnames) = varnames
      if(length(varnames)==0) return(NULL)
      names(varnames) = varnames
      nxt_vars = lapply(varnames, function(vn){
        nv = lapply(res_inner1, function(ri){
          ri[[nme_c1]][[nme_p1]][[vn]]
        })
        nv1 = list(var_names =nv[[1]]$var_names,
            cum_pv= .sumChisq(unlist(lapply(nv, function(nv1){
      unlist(nv1$pvs)
    }))),
       cumpv_all= .sumChisq(unlist(lapply(nv, function(nv1){
      unlist(nv1$pvs_all)
    }))))
    nv1
      })
  nxt_vars = nxt_vars[unlist(lapply(nxt_vars, length))>0]
  if(length(nxt_vars)==0) return(NULL)
  pvs_list = unlist(lapply(nxt_vars, function(nv) nv$cum_pv))
  nxt_vars[order(pvs_list)] 
    })
  })
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
    nme1="combined"
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
 
  clear_db=function(drop=F,exclude="vars", recursive=T){
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
 post_process=function(variables,beam=1){
   full_index = length(variables)
   vars_all = list()
   vars_all1 = list()
   vars_all2 = list()
   #vars_all3 = list() #funcstr
   names(variables) = 1:length(variables)
  # func_inds = lapply(variables, function(vv) attr(vv,"func_ind"))
   for(repn in names(variables)){
     full = repn==full_index
     
     var1 = variables[[repn]][[beam]]   ### only taking the top1
     var2 = var1$var_names
     if(length(var2)>0){
       names(var2) = names(var1$var_names)
       cumpv = var1$cum_pv#lapply(var1, function(vv) attr(vv,"cumpv"))
       varn = paste(names(var2),collapse=";") #paste(names(var2), collapse=";")
       if(is.null(vars_all[[varn]])){
         vars_all[[varn]] = list()
         vars_all1[[varn]] = var2
         vars_all2[[varn]] = list()
        # vars_all3[[varn]] =list(repn) 
        # names(vars_all3[[varn]]) = func_str1
       }else{
         #vars_all3[[varn]]=c( vars_all3[[varn]],repn)
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
   vars_combined = list(variables = vars_all1, inds = vars_all,cumpv=vars_all2)# ,transf= vars_all3) 
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
 
 nextVars=function(vars_l_todo, expt_id, k1,logpvthresh,stop_y="rand", verbose=F){
   vars_l = vars_l_todo$vars_l
   todo1 = vars_l_todo$todo1
   angles_all = lapply(vars_l, function(vars_l1){
     varnames = vars_l1$var_names; 
     res_inner2 = self$sigs$loadPvals(expt_id, varnames,k1) ## reconstruct ri
     nxt_vars1=.mergeResInner(res_inner2)
     nxt_vars1 =nxt_vars1[unlist(lapply(nxt_vars1, length))>0]
     if(length(nxt_vars1)==0) return(NULL)
     nxt_vars1
   })
   
   angles_all = angles_all[unlist(lapply(angles_all, length))>0]
   if(length(angles_all)==0){
  
     vars_l_todo = list(stop=length(todo1)==1, vars_l = vars_l, todo1 = todo1[-1])
     return(vars_l_todo)
   }
      ang1 = unlist(unlist(unlist(angles_all, rec=F),rec=F),rec=F)
      logpvs = unlist(lapply(ang1, function(a1)a1[["cum_pv"]]))
      logpvs_all = unlist(lapply(ang1, function(a1)a1[["cumpv_all"]]))
      ord = order(logpvs)
      names(ord) = names(logpvs)
      ord_all = order(logpvs_all)
      ang1 = ang1[ord_all]
      logpvs = logpvs[ord_all]
      logpvs_all = logpvs_all[ord_all]
      if(!is.null(stop_y)){
        gp1=grep(stop_y, names(logpvs))
        gp=grep(stop_y, names(logpvs), inv=T)
        
        # print("HERE")
        #print(unlist(list(rand= min(logpvs[gp1]),nonrand=min(logpvs[gp]))))
        stop_random= min(logpvs[gp1])<=min(logpvs[gp])
        #print(head(sort(ord[stop_ind])))
      }
      logpv =min(logpvs)
      
  if(stop_random){
    if(verbose) print(paste("stopping due to random", exp(logpv), names(logpvs)[which.min(logpvs)]))
  }
  ##ADD MORE RESTRICTIONS .. eg maxsize
      #     while( (length(vars_l[[1]]$var_names) < minsize || logpv<logpvthresh) && length(vars_l[[1]]$var_names)<maxsize && ! vars_l_todo$stop_random){
  if((!stop_random && logpv<=logpvthresh  ) ){
    if(verbose){
      print(head(sort(logpvs_all[gp])))
      print(names(vars_l))
    }
    dupls=(unlist(lapply(ang1, function(a1) paste(unlist(lapply(a1[[1]]$var_names, function(vv1)paste(vv1[1:2],collapse="::"))), collapse=";;"))))
    ang1 = ang1[!duplicated(dupls)]
    vars_l_todo = list(stop=F, vars_l = ang1, todo1 = vars_l_todo$todo1)
    return(vars_l_todo)
  }
      vars_l_todo = list(stop=length(todo1)==1, vars_l = vars_l, todo1 = vars_l_todo$todo1[-1], jj = vars_l_todo$jj+1)
      return(vars_l_todo)
},
savePvals=function(expt_id,k1, data_nme, vars_l, comb_){
  for(varn1 in names(vars_l)){ 
    self$sigs$savePvals(expt_id, data_nme, comb_[[varn1]], vars_l[[varn1]]$var_names,k1, useCurrVarnames = F)
  }
},
saveAngles=function(expt_id,k1, data_nme, vars_l, comb_){
  for(varn1 in names(vars_l)){ 
    self$saveAngles(expt_id, data_nme, comb_[[varn1]], vars_l[[varn1]]$var_names,k1)
  }
},

getTodo=function(flags, phens, logpv = -100){
  incls = fromJSON(.readFlag(flags,'data_types',toJSON(names(self$datasH[[1]]$data$data))))
  genes_incls=fromJSON(.readFlag(flags,"genes_incls",'{"all":["all"]}')) #,getOption("genes_incls",NULL)
  quantiles = sort(fromJSON(.readFlag(flags, "quantiles","[0]")),decreasing=T)
  names(incls) = incls;
  todo1 = unlist(unlist(lapply(incls, function(incl){
    lapply(genes_incls, function(g_incl){
      lapply(quantiles, function(qq){
        list(incl = incl, g_incl = g_incl, qq = qq)    
      })
    })
  }), rec=F), rec=F)
  Wall0 =lapply(phens, function(f) matrix(nrow=0,ncol=0))
  # var_thresh = var_thresh[match(names(var_thresh), train_nme)]
  vars_l_todo =
    list(
      todo1 = todo1,
      jj=0,
      logpv=logpv,
      vars_l = list(empty=stateObj$new(phens, NULL,NULL,NULL,NULL, var=c(), varnames=c(), Wall =Wall0)),
      stop=F
    )
  
  vars_l_todo
},
# self$select_k(phens,flags, k1,   var_thresh, quantiles)
 select_k=function(phens,flags, k1,   expt_id,
                     verbose=F){
  
   train_nme = .readFlag(flags,'train', names(self$datasH))
   names(train_nme)=train_nme
   beam = .readFlag(flags,"beam",1)
   maxsize=.readFlag(flags,'max',500)
   minsize=.readFlag(flags,'min',0)
   stop_y = .readFlag(flags, 'stop_y',"rand")
   logpvthresh = log(.readFlag(flags,"pthresh",1e-5))
  # var_thresh = var_thresh[match(names(var_thresh), train_nme)]
   saveAngles=F
 #  loadPv=.readFlag(flags, "loadPV",F)
  
  vars_l_todo = self$getTodo(flags, phens)
   #vars_l1 = vars_l[[1]];   k1=1;  qq =1; incl = incls[[1]]; data_nme = train_nme[[1]];g_incl  = genes_incls[[1]];#nme_c1 = names(transform_y)[[1]]
   while(length(vars_l_todo$todo1)>0){
     todo = vars_l_todo$todo1[[1]]
     g_incl = todo$g_incl; incl = todo$incl; qq_t = todo$qq
     vars_l = vars_l_todo$vars_l
     logpv = vars_l_todo$logpv
              invisible(lapply(train_nme, function(data_nme){
                 comb_ = self$datasH[[data_nme]]$multiAnglesAndPv(phens,k1,flags,expt_id, vars_l_todo , saveAngles=saveAngles, verbose=verbose)
                  if(saveAngles){
                     self$saveAngles(expt_id, k1, data_nme, vars_l, comb_ )  ## how do we get varnames
                   }else{
                     self$savePvals(expt_id,k1, data_nme, vars_l, comb_)
                   }
              }))
             if(saveAngles){  ## this is if separating angles from pval calc
               invisible(lapply(vars_l, function(vars_l1){
                 varnames = vars_l1$var_names
                 comb_= .mergeComb( self$sigs$loadAngles(expt_id,varnames), flags) ## reconstruct comb_angs1
                 res_inner1 = lapply(train_nme, function(data_nme){
                   prev_i = vars_l1[[data_nme]]
                   ri = self$datasH[[data_nme]]$res_inner(comb_,prev_i,flags,k1,expt_id)
                   self$sigs$savePvals(expt_id, data_nme, ri, varnames,useCurrVarnames=F)
                   ri
                 })
               }))
             }
             vars_l_todo = self$nextVars(vars_l_todo, expt_id, k1,logpvthresh,stop_y=stop_y, verbose=verbose)
           }
   vars_l_todo$vars_l
 },
update=function(phens, flags, verbose=F){
  nreps = lapply(self$datasH, function(dh){
    dh$update(phens, flags)
  })
  nreps[[1]]
},
  select=function(phens,flags, verbose=F, useDB=T ){#c(y="function(y) y","function(y) y")
    nreps = self$update(phens, flags, verbose=verbose);
    sigDB =self$sigs
    if(!is.null(self$sigs) && useDB ){
      vars_all = self$sigs$loadVars(flags, phens)
      if(!is.null(vars_all)) return(vars_all)
    }
    expt_id=self$sigs$getExpt(flags=flags, phens = phens, add_new=T)
     variables=lapply(nreps, function(k1){
       if(verbose) print(paste("cv",k1,"of",length(nreps)))
      self$select_k(phens,flags, k1,   expt_id, verbose=verbose)
    })
    vars_all=self$post_process(variables,beam=1)
    vars_all$flags = flags; vars_all$phens = phens; # vars_all$transform_y = transform_y ;
    if(length(vars_all$variables)==0) return(vars_all)
    if(useDB){
      self$sigs$saveVars(vars_all,replace=T)
    }
    vars_all
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
