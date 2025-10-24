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
.mergeResInner<-function(res_inner1, comb_){
  nme_comb = names(comb_); names(nme_comb) = nme_comb
  res_inner=lapply(nme_comb, function(nme_c1){
    nmesp1 = names(comb_[[nme_c1]]); names(nmesp1) = nmesp1
    lapply(nmesp1, function(nme_p1){
      comb = comb_[[nme_c1]][[nme_p1]]
      if(nrow(comb)==0) return(NULL)
      nme_comb = names(comb_); names(nme_comb) = nme_comb
      varnames = comb$names; names(varnames) = varnames
      nxt_vars = lapply(varnames, function(vn){
        nv = lapply(res_inner1, function(ri){
          ri[[nme_c1]][[nme_p1]][[vn]]
        })
    attr(nv,"cumpv")= .sumChisq(unlist(lapply(nv, function(nv1){
      unlist(nv1$pvs)
    })))
    attr(nv,"cumpv_all")= .sumChisq(unlist(lapply(nv, function(nv1){
      unlist(nv1$pvs_all)
    })))
    nv
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
   full_index = length(variables)
   vars_all = list()
   vars_all1 = list()
   vars_all2 = list()
   #vars_all3 = list() #funcstr
   names(variables) = 1:length(variables)
  # func_inds = lapply(variables, function(vv) attr(vv,"func_ind"))
   for(repn in names(variables)){
     full = repn==full_index
     var1 = variables[[repn]]
     var2 = var1[[1]]$var_names
     if(length(var2)>0){
       names(var2) = var1[[1]]$varnames
       cumpv = attr(var1,"cumpv")#lapply(var1, function(vv) attr(vv,"cumpv"))
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

  select=function(phens,flags, verbose=F, user="",useDB=T ){#c(y="function(y) y","function(y) y")
    sigDB =self$sigs
    if(!is.null(sigDB) && useDB ){
      vars_all = sigDB$loadVars(flags, phens)
      if(!is.null(vars_all)) return(vars_all)
    }
    expt_id=sigDB$getExpt(flags=flags, phens = phens, add_new=T)
    for(dh in self$datasH){
      dh$updateLOOC(phens, flags, verbose=verbose)
      dh$updateTrain(phens, flags,verbose=verbose)
    }
    datas1 = self$datasH
    train_nme = .readFlag(flags,'train', names(self$datasH))
    train_nme = train_nme[train_nme %in% names(self$datasH)]
    quantiles = sort(fromJSON(.readFlag(flags, "quantiles","[0]")),decreasing=T)
    genes_incls=fromJSON(.readFlag(flags,"genes_incls",'{"all":["all"]}')) #,getOption("genes_incls",NULL)
    if(!is.list(genes_incls)) stop("genes incls should be list")
    if(length(train_nme)==0) train_nme = names(datas1)[[1]]
    names(train_nme) = train_nme
    maxsize=.readFlag(flags,'max',500)
    minsize=.readFlag(flags,'min',0)
  
    incls = fromJSON(.readFlag(flags,'data_types',"{}")) 
    if(length(incls) == 0 )incls = list("all"=names(self$datasH[[1]]$data$data))
    incls_all = unique(unlist(incls))
    logpvthresh = log(.readFlag(flags,"pthresh",1e-5))
    logpv=-100
    nreps1 =self$datasH[[1]]$nreps()
    nreps = 1:nreps1
    names(nreps) = nreps
    beam = .readFlag(flags,"beam",1)
    stop_y = .readFlag(flags, 'stop_y',"rand")
    stop_random=F
    var_thresh = lapply(train_nme, function(data_nme){
      lapply(self$datasH[[data_nme]]$data$vars, function(v) quantile(v, quantiles))
    })
    Wall0 =lapply(phens, function(f) matrix(nrow=0,ncol=0)) 
    # k1=1;  qq =1; incl = incls[[1]]; data_nme = train_nme[[1]];g_incl  = genes_incls[[1]];#nme_c1 = names(transform_y)[[1]]
     variables=lapply(nreps, function(k1){
      if(verbose) print(paste("cv",k1,"of",length(nreps)))
      jj1=0
      #invisible(lapply(train_nme, function(data_nme) datas1[[data_nme]]$update(k, funcst, phens,incls_all))); ### update training object
     # res2=  lapply(phens_index, function(p_index){
        if(verbose) print(paste(k1,length(nreps)))
       
#        invisible(lapply(train_nme, function(data_nme) datas1[[data_nme]]$update(k, funcst, subphens,incls_all))); ### update training object
      
        vars_l =list(lapply(train_nme,function(xx) stateObj$new(phens, NULL,NULL,NULL,NULL,k1, var=c(), varnames=c(), Wall =Wall0)))
      names(vars_l) = "empty"
      #vars_l1 = vars_l[[1]]
        for(incl in incls){
         if(verbose) print(incl)
         for(g_incl in genes_incls){
         for(qq in 1:length(quantiles)){
         while( (length(vars_l[[1]][[1]]$var) < minsize || logpv<logpvthresh) && length(vars_l[[1]][[1]]$var)<maxsize && ! stop_random){
           updateV = T  ## whether to update vars each iteration  , should only be F for debugging
           {
          angles_all = lapply(vars_l, function(vars_l1){
            nxt_vars1 =  tryCatch({
              varnames = vars_l1[[1]]$var_names; type = self$datasH[[1]]$type
              for(data_nme in train_nme){
                comb_angs1=self$datasH[[data_nme]]$combinedAngles(phens, varnames, incl, k1, type,var_thresh[[data_nme]], g_incl, qq, flags)
                sigDB$saveAngles(expt_id, data_nme, comb_angs1,varnames ) 
              }
              comb_all1=  sigDB$loadAngles(expt_id,varnames) ## reconstruct comb_angs1
               comb_ =  .mergeComb(comb_all1, flags) 
                res_inner1 = lapply(train_nme, function(data_nme){
                 prev_i = vars_l1[[data_nme]]
                 ri = self$datasH[[data_nme]]$res_inner(comb_,prev_i,flags,k1)
                 sigDB$savePvals(expt_id, data_nme, ri, varnames)
                 ri
               })
                if(.readFlag(flags, "loadPV",F)){
                  res_inner2 = sigDB$loadPvals(expt_id, varnames) ## reconstruct ri
                  res_inner=.mergeResInner(res_inner2, comb_)
                }else{
             res_inner=.mergeResInner(res_inner1, comb_)
                }
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
            if(updateV) vars_l = ang1[1:min(length(ang1),beam)]
          }
          if(verbose){
            print(head(sort(logpvs_all[gp])))
            print(names(vars_l))
            print(paste("logpv",logpv,min(logpvs_all), jj1))
            if(updateV) jj1 = jj1+1
          }
         }
         
        }
        }
        }
      }
       
        #lapply(datas1, function(d)d$saveParquet())
        #lapply(vars_l, function(v) v$var)
          ##just take the top
       # if(length(vars_l[[1]]$var)>0) print(vars_l)
      #attr(vars_l[[1]],"func_ind")=func_ind#$var
        vars_l[[1]]   ## return the best one .. perhaps we should allow to return more? 
    #     res2[unlist(lapply(res2,function(xx) length(xx[[1]]$var)))>0]
    })
    vars_all=self$post_process(variables)
    vars_all$flags = flags; vars_all$phens = phens;  vars_all$transform_y = transform_y ;
    if(length(vars_all$variables)==0) return(vars_all)
    if(useDB){
      sigDB$saveVars(vars_all,replace=T)
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
