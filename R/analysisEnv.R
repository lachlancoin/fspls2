


#' main function for variable selection
#' @param datasH a list of dataH objects
#' @param flags list of options
#' @param transform_y transformation object 
#' @param phens list of phenotypes
#' @param dbDir database dir, can be NULL
#' @export
fspls.select<-function(datasH, flags,
                transform_y, 
                phens=datasH[[1]]$pheno()$all,
                dbDir = NULL
               ){#c(y="function(y) y","function(y) y")
  analysis =analysisEnv$new(flags=flags, dbDir=dbDir) ;
  analysis$updateExpt(phens1, flags)
  phens1 = phens
  mc.cores = .readFlag(flags, "mc.cores",1)
  flags$transform_y = toJSON(transform_y);
  verbose=.readFlag(flags,'verbose',F);
  if(flags$topn<flags$beam) stop("beam should be less than topn")
  if(is.null(flags[['data_types']]) || flags[['data_types']]=="{}")flags[['data_types']]=toJSON(datasH[[1]]$data_types())
  
  nreps_all = lapply(datasH, function(dh){
    dh$update(phens1, flags, transform_y);
  })
  nreps = nreps_all[[1]]
  
  ## load if already calculated
    vars_all = analysis$loadVars(phens);
    if(!is.null(vars_all)) return(vars_all)
  
  vars_l_todo = analysis$getTodo(flags, phens1)

 
  
  variables1 <-lapply(nreps, function(k1) {
    vars_l_todo =   analysis$select_k(datasH, phens1, flags, k1, vars_l_todo)
    vars_l_todo$vars_l 
  })
  
  vars_all=post_process(variables1,flags,phens1)
  analysis$saveVars(vars_all);
  vars_all
  
}


.extractFullVars<-function(vars_all0){
  lapply(vars_all0, function(vars_all){
  subinds = which(unlist(lapply(vars_all$inds, function(x) length(grep('full',names(x)))))>0)
  list(variables = vars_all$variables[subinds], inds = vars_all$inds[subinds], cumpv = vars_all$cumpv[subinds],
       transf = vars_all$transf[subinds],
       flags = vars_all$flags,
       phens=vars_all$phens, transform_y = vars_all$transform_y)
  })
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
  #nme_c1 = nme_comb[[1]]; nme_p1 = names(res_inner1[[1]][[nme_c1]])[1]
  
  res_inner=lapply(nme_comb, function(nme_c1){
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
                 angle = sum(unlist(lapply(nv, function(nv1) unlist(nv1$angle)))),
              cum_angle =sum(unlist(lapply(nv, function(nv1) unlist(nv1$angles)))),
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


#' Analysis  Class
#'
#' @description
#' A class that encapsulates analysis of one or more datasets. Only expert users wlil need access to this
#'
#' @details
#' Available methods:
#' \itemize{
#'   \item \code{new()} - Create a new instance
#'   \item \code{name()} - Get dataset name
#'   \item \code{clone()} - Clone the object
#' }
#'  @export
analysisEnv<-R6::R6Class("analysisEnv", 
            inherit = analysisBase,
  private = list(
    type="character",
    sigsdir="character",
    sigs="environment",
    flags ="list",
    updateData=function(datasH, data_names =names(datasH), 
                        data_types = names(datasH[[1]]$data$data),
                        dims   = lapply(datasH, function(data) data$dims())
    ){
      private$transform_y= private$sigs$updateData(data_flags = private$flags, 
                                                   data_names =data_names, 
                                                   data_types = data_types,
                                                   dims = dims,
                                                   transform_y = private$transform_y)
    },
    nextVars=function(flags, phens, vars_l_todo,  k1,
                                          logpvthresh,beam,  comb2_news = NULL,stop_y="rand", verbose=F){
      vars_l = vars_l_todo$vars_l
      todo1 = vars_l_todo$todo1
     # expt_id=super$getExpt(flags, phens, add_new=T)

      useAngles = !is.null(flags$angles_only) && flags$angles_only
      nme_l = names(vars_l); names(nme_l) = nme_l
      angles_all = lapply(nme_l, function(nme_l1){
        vars_l1 = vars_l[[nme_l1]]
        varnames = vars_l1$var_names; 
        #if(!is.null(varnames) && is.null(names(varnames))) names(varnames) = lapply(varnames, paste, collapse=".")
        if(length(comb2_news)>0){
         #  if(length(varnames)==0) varnames = "empty"
           res_inner1 = lapply(comb2_news, function(c2n) c2n[[nme_l1]])
        }else{
            res_inner1 = private$sigs$loadPvals(private$expt_id(), varnames,k1) ## reconstruct ri
        }
        nxt_vars1=.mergeResInner(res_inner1)
        nxt_vars1 =nxt_vars1[unlist(lapply(nxt_vars1, length))>0]
        if(length(nxt_vars1)==0) return(NULL)
        nxt_vars1
      })
      
      angles_all = angles_all[unlist(lapply(angles_all, length))>0]
      if(length(angles_all)==0){
        vars_l_todo = list(stop=length(todo1)==1, vars_l = vars_l, todo1 = todo1[-1])
        return(vars_l_todo)
      }
      ang1 = unlist(unlist(unlist(angles_all, recursive=FALSE),recursive=FALSE),recursive=FALSE)
      angles_ = unlist(lapply(ang1, function(a1)a1[["angle"]]))
      
      vn = unlist(lapply(ang1, function(a1)paste(names(a1[["var_names"]]), collapse=";")), recursive=FALSE)
      names(ang1) = vn 
      
      logpvs =if(useAngles) angles_ else   unlist(lapply(ang1, function(a1)a1[["cum_pv"]]))
      logpvs_all = if(useAngles) unlist(lapply(ang1, function(a1)a1[["cum_angle"]])) else  unlist(lapply(ang1, function(a1)a1[["cumpv_all"]]))
      
      
      ord = order(logpvs)
      names(ord) = names(logpvs)
      ord_all = order(logpvs_all)
      ang1 = ang1[ord_all]
      logpvs = logpvs[ord_all]
      logpvs_all = logpvs_all[ord_all]
      angles_ = angles_[ord_all]
      if(!is.null(stop_y)){
        gp1=grep(stop_y, names(logpvs))
        gp=grep(stop_y, names(logpvs), inv=T)
        
        # print("HERE")
        stop_random =  min(logpvs[gp1]) < min(logpvs[gp]) 
        stop_random1= min(angles_[gp1]) < min(angles_[gp])
        stop_random2= min(logpvs_all[gp1]) < min(logpvs_all[gp])
        
        #print(unlist(list(rand= min(logpvs[gp1]),nonrand=min(logpvs[gp]))))
        # stop_random= min(gp1)<=min(gp)
        #print(head(sort(ord[stop_ind])))
        if(verbose && FALSE){
              print(paste("COMPARING TO RANDOM!!!!! useAngles=", useAngles))
              print(unlist(list(rand=min(logpvs[gp1]),nonrand= min(logpvs[gp]))))
              print("cumulative ")
              print(unlist(list(rand=min(logpvs_all[gp1]),nonrand= min(logpvs_all[gp]))))
        }
        if(!useAngles && verbose) print(unlist(list(rand=min(angles_[gp1]),nonrand= min(angles_[gp]))))
        
      }
      logpv =min(logpvs)
      
      
      if(stop_random || stop_random1 || stop_random2){
        if(verbose) print(paste("stopping due to random", exp(logpv), names(logpvs)[which.min(logpvs)]))
      }
      
      ##ADD MORE RESTRICTIONS .. eg maxsize
      #     while( (length(vars_l[[1]]$var_names) < minsize || logpv<logpvthresh) && length(vars_l[[1]]$var_names)<maxsize && ! vars_l_todo$stop_random){
      if((!stop_random &&  ! stop_random1 && !stop_random2 && logpv<=logpvthresh  ) ){
        if(verbose){
          print(head(sort(logpvs_all[gp])))
          print(names(vars_l))
        }
        dupls=(unlist(lapply(ang1, function(a1) paste(unlist(lapply(a1$var_names, function(vv1)paste(vv1[1:2],collapse="::"))), collapse=";;"))))
        ang1 = ang1[!duplicated(dupls)]
        ang1 = ang1[grep("rand", names(ang1), inv=T)]
        ang1 = ang1[1:min(length(ang1),beam)]
        vars_l_todo = list(stop=F, vars_l = ang1, todo1 = vars_l_todo$todo1)
        if(length(grep("rand", names(vars_l_todo$vars_l)))>0) stop("!!");
        
        return(vars_l_todo)
      }
     # if(length(todo1)==1){
    #    print("could consider saving the vars at this point to the DB.  Maybe also need to record dataset included")
    #  }
     # print("shortening _todo")
      vars_l_todo = list( stop=length(todo1)==1, vars_l = vars_l, todo1 = vars_l_todo$todo1[-1], jj = vars_l_todo$jj+1)
      return(vars_l_todo)
    },
    savePvals=function(flags,phens,k1, data_nme, vars_l, comb2){
      for(varn1 in names(vars_l)){ 
        ri=comb2[[varn1]]
        varnames=vars_l[[varn1]]$var_names
        super$savePvals(flags, phens, ri, varnames, k1, data_nme=data_nme, useCurrVarnames=F)

      }
    },
    saveAngles=function(flags,phens,k1, data_nme, vars_l, comb_){
      for(varn1 in names(vars_l)){ 
        private$saveAngles(flags,phens, data_nme, comb_[[varn1]], vars_l[[varn1]]$var_names,k1)
      }
    }
    
   
    
  ),                      
  public = list(
  
    #' @description Create a new instance
    #' @param flags a list object specifying options
    #' @param dbDir location for database to be written. If NULL then no DB 
  initialize=function(
              flags = list(),
              dbDir=tempdir()
            
                    ){
    super$initialize("combined",flags=flags, dbDir=dbDir);
       
  },
  
  #' selection for a single fold
  #' @param datasH a list of dataH objects
  #' @param phens1 list of phenotypes
  #' @param flags list of options
  #' @param k1  the fold
  #' @param vars_l_todo an object representing what is left to do
  #' @returns vars_l_todo object
  select_k=function(datasH,phens1,flags, k1,
                    vars_l_todo 
                  ){
    verbose = .readFlag(flags,"verbose",F);
    print(paste("fold",k1))
    expt_id = private$expt_id();
    # get_plots=.readFlag(flags, "get_plots",F) 
    stop_y = .readFlag(flags, 'stop_y',"rand")
    logpvthresh = log(.readFlag(flags,"pthresh",0.1))
    beam= .readFlag(flags,"beam",1)
    saveAngles=F
    #plot_results = list()
    # vars_l = analysis$nextVars(expt_id, flags)
    nvar=0;
    nmesH = lapply(datasH, function(dh)dh$name());
    names(datasH) = nmesH;
    names(nmesH) = nmesH
    comb2 = lapply(datasH, function(x) return(NULL))
    #  comb2_old = comb2
    useDB = !is.null(private$sigs)
    # nmesH = names(datasH); names(nmesH) = nmesH;
    while(length(vars_l_todo$todo1)>0 ){
      comb2_news =try(lapply(nmesH, function(nmeh){
        dh = datasH[[nmeh]]
        comb20=comb2[[nmeh]]
        comb_new  = dh$multiAnglesAndPv(comb20 , phens1, k1,flags,expt_id, vars_l_todo)
        comb21 =  lapply(comb_new, function(x) x$pvs)
        private$savePvals(flags,phens1,k1, dh$name(), vars_l_todo$vars_l,comb21)# no need to save here, just keep
        comb21
      }))
      if(inherits(comb2_news,"try-error")) break;
      
      vars_l_todo_new= private$nextVars(flags,phens1, vars_l_todo,  k1,logpvthresh,beam, 
                                        comb2_news=if(useDB) NULL else comb2_news, 
                                        stop_y = stop_y, verbose=verbose)
      if(length(vars_l_todo$todo1)==length(vars_l_todo_new$todo1)){
        comb2 = comb2_news
      }
      vars_l_todo = vars_l_todo_new
      
      nvar = length(vars_l_todo$vars_l[[1]]$var)
      if(verbose) print(names(vars_l_todo$vars_l))
      if(length(vars_l_todo$vars_l[[1]]$var_names)>=flags$max) break;
    }
    vars_l_todo
  },
  
  #' Get todo 
  #'
  #' @param flags list of options
  #' @param phens phenotypes
  #' @return object outlining what is left to do for next iteration
  getTodo=function(flags, phens){
    logpv = -100
    incls = fromJSON(.readFlag(flags,'data_types','{}'))
    genes_incls=fromJSON(.readFlag(flags,"genes_incls",'{"all":["all"]}')) #,getOption("genes_incls",NULL)
    quantiles = sort(fromJSON(.readFlag(flags, "quantiles","[0]")),decreasing=T)
    # names(incls) = incls;
    todo1 = unlist(unlist(lapply(incls, function(incl){
      lapply(genes_incls, function(g_incl){
        lapply(quantiles, function(qq){
          list(incl = incl, g_incl = g_incl, qq = qq)    
        })
      })
    }), recursive=FALSE), recursive=FALSE)
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
  
  
  #' Clear the databases.  Only necessary if using a database and re-running analyses
  #'
  #' @param drop completely drop tables and start from scratch
  #' @param exclude tables to exclue from clearing
  #' @param datasH  list of dataH objects to clear
   clear_db=function(drop=F,exclude="vars", datasH = NULL){
    if(drop){
      if(!is.null(private$sigs))private$sigs$drop_all(exclude=exclude)
      if(!is.null(datasH)){
        lapply(datasH, function(dh) dh$clear_db(drop=drop, exclude=exclude))
      }
    }else{
      warning("need to set drop=T if you are sure, this will delete all saved signatures")
    }
    
  },  
  
 
  #' Get the next vars in iteration.  Internal function
  #'
  #' @param flags list of options
  #' @param phens phenotypes
  #' @param vars_l_todo vars_l_todo
  #' @param comb2 results
  #' @param data_nme name of dataset
  #' @param k1 which repetition
  #' @return object outlining what is left to do

 savePvalsAndNextVars=function(flags,phens, vars_l_todo,comb2,data_nme, k1){
   beam=.readFlag(flags,'beam',1);
   stop_y=.readFlag(flags,"stop_y","rand")
   verbose=.readFlag(flags,"verbose",F);
   logpvthresh = log(.readFlag(flags,'pthresh',0.05));
   angles_only = .readFlag(flags,'angles_only',FALSE);
   if(angles_only) logpvthresh =0;
      private$savePvals(flags,phens,k1, data_nme, vars_l_todo$vars_l,comb2)
    private$nextVars(flags,phens, vars_l_todo,  k1,logpvthresh,beam, stop_y = stop_y, verbose=verbose)
 },
 
 
 #' main function for variable selection
 #' @param datasH a list of dataH objects
 #' @param phens list of phenotypes
 #' @param flags list of options
 #' @param transform_y transformation object 
 #' @param useDB boolean to indicate if results should be saved to database
 select=function(datasH, flags,
                 transform_y, 
                 phens=datasH[[1]]$pheno()$all,
                 useDB=F){#c(y="function(y) y","function(y) y")
   phens1 = phens
   super$updateExpt( phens1, flags)
   mc.cores = .readFlag(flags, "mc.cores",1)
   ##if(mc.cores>1 && )
   flags$transform_y = toJSON(transform_y);
   verbose=.readFlag(flags,'verbose',F);
   if(flags$topn<flags$beam) stop("beam should be less than topn")
   if(is.null(flags[['data_types']]) || flags[['data_types']]=="{}")flags[['data_types']]=toJSON(datasH[[1]]$data_types())
   
   nreps_all = lapply(datasH, function(dh){
     dh$update(phens1, flags, transform_y);
   })
   nreps = nreps_all[[1]]
   
   
     vars_all = super$loadVars(phens)
     if(!is.null(vars_all)) return(vars_all)
  
   vars_l_todo = self$getTodo(flags, phens1)

   
   variables1 <- lapply(nreps, function(k1) {  ## can use mclapply here
      vars_l_todo =   self$select_k(datasH, phens1, flags, k1,  vars_l_todo)
      vars_l_todo$vars_l 
   })
   

   vars_all=post_process(variables1,flags,phens1)
   self$saveVars(vars_all)
   vars_all
   
 }
 


  


))
