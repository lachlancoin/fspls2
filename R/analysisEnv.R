



plot_traj<-function(comb_plot, y="value"  ,facet="data~maxsig", keep_best=10, txtsize=5, step=2){ #y="cumulative";
  if(!is.null(comb_plot$nrep)){
    if(length(unique(comb_plot$nrep))>1) stop(" need to subset on nrep first")
  }
  comb_plot$maxsig = gsub("\n",";", comb_plot$maxsig)
  maxl = min(comb_plot[[y]])
  
  maxsigl1 = unlist(lapply(comb_plot$maxsig, function(x) paste(sort(strsplit(x,";")[[1]]), collapse=";")))
  maxsigl1_u = sort(table(maxsigl1), decreasing=TRUE)
  conv = lapply(names(maxsigl1_u), function(x){
    names(sort(table(comb_plot$maxsig[which(maxsigl1==x)]),decreasing=TRUE))[1]
  })
  names(conv) = names(maxsigl1_u)
  maxsig = unlist(lapply(maxsigl1, function(x) conv[x]))
  comb_plot$maxsig = maxsig
  maxsigl = sort(table(maxsig),decreasing=TRUE)
  
  if(length(maxsigl1_u)>keep_best){
    
    ms1   = factor(maxsig, names(maxsigl)[maxsigl>=maxsigl1_u[keep_best]])
    comb_plot1 =comb_plot[!is.na(ms1),,drop=FALSE]
  }else{
    comb_plot1 = comb_plot
  }
  ms1 = unique(comb_plot$maxsig); names(ms1) = ms1
  levs = sort(unlist(lapply(ms1, function(ms){
    inds1 = comb_plot$maxsig==ms
    min(comb_plot$cumulative[inds1])
  })),decreasing=FALSE)
  labels = unlist(lapply(names(levs), function(lev){
    maxs = strsplit(lev,";")[[1]];
    paste(unlist(lapply(seq.int(1, length(maxs), by=step), function(st){
      
      mi = min(length(maxs), st+step)
      paste(maxs[st:mi], collapse=";")
    })), collapse="\n")
  }))
  
  comb_plot1$maxsig = factor(comb_plot1$maxsig, levels = names(levs), labels=labels)
  comb_plot1[[y]] = -1*comb_plot1[[y]]
  ggp=ggplot(comb_plot1, aes_string("nvar" ,y, color="maxsig"))+geom_point()+geom_line(alpha=.1)+facet_grid(facet)
  ggp=ggp+guides(color = "none")+theme(
    strip.text = element_text(size = txtsize, face = "bold", color = "black")
  )+geom_hline(yintercept = -1*maxl)+scale_y_log10()
  ggp
}


## expands to replace proportion with NA
expandData<-function(dataset, mult = 100 
){
  na_inds = which(dataset$certainty<1)
    non_na_inds = 1:nrow(dataset$y)
  non_na_inds = non_na_inds[-na_inds]
      
  y_alt = 1- dataset$y[na_inds,,drop=FALSE]
  if(max(abs(apply(cbind(sort(unique(y_alt[,1])),c(0,1)),1,diff)))>0) stop("need to have 0 1 values for y")
  if(max(dataset$certainty)>1.0) stop("!!")
  if(min(dataset$certainty<0.0)) stop("!!")
  rownames(y_alt) = paste0(rownames(dataset$y)[na_inds], ".alt")
  
  y_new = rbind(dataset$y, y_alt)
  weights = c(dataset$certainty, 1-dataset$certainty[na_inds])* mult
  names(weights) = c(rownames(dataset$y), paste0(rownames(dataset$y)[na_inds],".alt"))
  d_new =  lapply(dataset$dataset, function(d){
    d1 = d[na_inds,,drop=FALSE];
    rownames(d1) = paste0(rownames(d)[na_inds],".alt")
    d_out = rbind(d,d1)
    d_out
  })
  list(y = y_new, dataset=d_new, na_inds=na_inds,non_na_inds = non_na_inds, 
       orig_inds = 1:nrow(dataset$y), alt_inds = nrow(dataset$y)+1:length(na_inds),
              weights =weights)
  
}



#' main function for variable selection with unknown values
#' @param dataset a list with y and data value
#' @param flags list of options
#' @param transform_x transformation object 
#' @param phens list of phenotypes
#' @param dbDir y_orig the original y (for testing only)
fspls.iterative<-function(dataset,flags, transform_x ){
  options(flags);
  if(is.null(dataset$y) || length(dataset$dataset)==0 || is.null(dataset$certainty)) stop(" dataset not well defined")
  if(!is.factor(dataset$y[[1]])) stop("y should be a factor")
  flags$nfold=1; flags$verbose=F;
  mult = .readFlag(flags, "mult",100)
  dh = dataH$new(dataset$dataset, y = dataset$y, 
                 certainty = dataset$certainty,
                 nme="iterative", flags=flags, transform_x = transform_x, dbDir=NULL)
 vars_all = NULL; all_models = NULL; preds1 = NULL;
 updates = list();
 select_each_iteration = .readFlag(flags, "select_each_iteration",TRUE);
 datasH = list(dh);
 k=0
 phens = dh$pheno()$all
 analysis =analysisEnv$new(flags=flags, dbDir=NULL) ;
  repeat  {
    if(k==0 || select_each_iteration) {
            
          vars_all = analysis$select(datasH,  flags, transform_x, phens = phens)
      
              #vars_all = fspls.select(datasH,  flags, transform_x, phens = dh$pheno()$all)
    }
    k = k+1
   all_models =dh$makeAllModels(vars_all,useDB=FALSE)
   updated =  dh$updateWeights(all_models)
     updates[[k]] = updated
   print(paste("sumdiff",k,abs(updated$sumdiff)))
  if(abs(updated$sumdiff)<=1e-5 && k>2) {
        message("breaking since probs ordered not improving")
        break;
        
      }
  }
  y_new =     dh$getYNew();

   list(dh = dh, vars_all = vars_all, all_models = all_models,updates = updates, 
        y_new = y_new$y, certainty_new = y_new$certainty,error_rate = y_new$error_rate
        
         )
#  plot(roc1)
  
}
#' extracts the full model variables from a cross validation run
#' @param vars_all and object returned by fspls.select
#' @export
extractFullVars<-function(vars_all){
  lapply(vars_all, function(vars_all1){
    subinds = which(unlist(lapply(vars_all1$inds, function(x) length(grep('full',names(x)))))>0)
    list(variables = vars_all1$variables[subinds], inds = vars_all1$inds[subinds], cumpv = vars_all1$cumpv[subinds],
         transf = vars_all1$transf[subinds],
         flags = vars_all1$flags,
         phens=vars_all1$phens, transform_x = vars_all1$transform_x)
  })
}

#' main function for variable selection
#' @param datasH a list of dataH objects
#' @param flags list of options
#' @param transform_x transformation object 
#' @param phens list of phenotypes
#' @param dbDir database dir, can be NULL
#' @export
fspls.select<-function(datasH, flags,
                transform_x, 
                phens=datasH[[1]]$pheno()$all,
                dbDir = NULL
               ){#c(y="function(y) y","function(y) y")
  options(flags);
  analysis =analysisEnv$new(flags=flags, dbDir=dbDir) ;
  analysis$updateExpt(phens, flags)
  #phens1 = phens
  mc.cores = .readFlag(flags, "mc.cores",1)
  flags$transform_x = toJSON(transform_x);
  verbose=.readFlag(flags,'verbose',FALSE);
  if(flags$topn<flags$beam) stop("beam should be less than topn")
  if(is.null(flags[['data_types']]) || flags[['data_types']]=="{}")flags[['data_types']]=toJSON(datasH[[1]]$data_types())
  
  nreps_all = lapply(datasH, function(dh){
    dh$update(phens, flags, transform_x);
  })
  nreps = nreps_all[[1]]
  
  ## load if already calculated
    vars_all = analysis$loadVars(phens);
    if(!is.null(vars_all)) return(vars_all)
  
  vars_l_todo = analysis$getTodo(flags, phens)
  variables1 <-lapply(nreps, function(k1) {
    vars_l_todo =   analysis$select_k(datasH, phens, flags, k1, vars_l_todo)
    vars_l_todo$vars_l 
  })
  
  vars_all=post_process(variables1,flags,phens)
  analysis$saveVars(vars_all);
  vars_all
  
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
.extractFullModels<-function(all_models, full_model_only=FALSE){
  subinds = which(unlist(lapply(all_models$models, function(x) length(grep('full', names(x)))))>0)
  models1 = all_models$models[subinds]
  res1 =list(flags = all_models$flags, phens=all_models$phens, transform_x = all_models$transform_x)
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
      private$transform_x= private$sigs$updateData(data_flags = private$flags, 
                                                   data_names =data_names, 
                                                   data_types = data_types,
                                                   dims = dims,
                                                   transform_x = private$transform_x)
    },
    nextVars=function(flags, phens, vars_l_todo,  k1,
                                          logpvthresh,beam,  comb2_news = NULL,stop_y="rand", verbose=FALSE){
      vars_l = vars_l_todo$vars_l
      todo1 = vars_l_todo$todo1
     # expt_id=super$getExpt(flags, phens, add_new=TRUE)

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
        gp=grep(stop_y, names(logpvs), inv=TRUE)
        
      
        stop_random =  min(logpvs[gp1], na.rm=TRUE) < min(logpvs[gp], na.rm=TRUE) 
       # stop_random1= min(angles_[gp1]) < min(angles_[gp])
        stop_random2= min(logpvs_all[gp1], na.rm=TRUE) < min(logpvs_all[gp], na.rm=TRUE)
        

        # stop_random= min(gp1)<=min(gp)

        if(verbose){
          print(paste(stop_random, stop_random2))
              print(paste("COMPARING TO RANDOM!!!!! useAngles=", useAngles))
              print(unlist(list(rand=min(logpvs[gp1]),nonrand= min(logpvs[gp]))))
              print("cumulative ")
              print(unlist(list(rand=min(logpvs_all[gp1]),nonrand= min(logpvs_all[gp]))))
        }
        if(!useAngles && verbose) {
          print(unlist(list(rand=min(logpvs[gp1]),nonrand= min(logpvs[gp]))))
          print(unlist(list(rand=min(angles_[gp1]),nonrand= min(angles_[gp]))))
        }
        
      }
      logpv =min(logpvs)
      
      
      if(stop_random || stop_random2){
        if(verbose) print(paste("stopping due to random", exp(logpv), names(logpvs)[which.min(logpvs)]))
      }
      
      ##ADD MORE RESTRICTIONS .. eg maxsize
      #     while( (length(vars_l[[1]]$var_names) < minsize || logpv<logpvthresh) && length(vars_l[[1]]$var_names)<maxsize && ! vars_l_todo$stop_random){
      if((!stop_random && !stop_random2 && logpv<=logpvthresh  ) ){
        if(verbose){
          print(head(sort(logpvs_all[gp]),beam))
          print(names(vars_l))
        }
        dupls=(unlist(lapply(ang1, function(a1) paste(unlist(lapply(a1$var_names, function(vv1)paste(vv1[1:2],collapse="::"))), collapse=";;"))))
        ang1 = ang1[!duplicated(dupls)]
        last_non_rand = grep("rand", names(ang1))[1]-1
        if(is.na(last_non_rand)) last_non_rand = beam;      
        ang1 = ang1[1:min(length(ang1),beam, last_non_rand)]
        vars_l_todo = list(stop=FALSE, vars_l = ang1, todo1 = vars_l_todo$todo1)
        if(length(grep("rand", names(vars_l_todo$vars_l)))>0) stop("!!");
        
        return(vars_l_todo)
      }

      vars_l_todo = list( stop=length(todo1)==1, vars_l = vars_l, todo1 = vars_l_todo$todo1[-1], jj = vars_l_todo$jj+1)
      return(vars_l_todo)
    },
    savePvals=function(flags,phens,k1, data_nme, vars_l, comb2){
      for(varn1 in names(vars_l)){ 
        ri=comb2[[varn1]]
        varnames=vars_l[[varn1]]$var_names
        super$savePvals(flags, phens, ri, varnames, k1, data_nme=data_nme, useCurrVarnames=FALSE)

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
  #' gets the trajectory
  #' @param datasH a list of dataH objects after variable selection
  plot_trajectory=function(datasH){
    
    traj = lapply(datasH, function(dh){
      dh$get_trajectory();
    })
    names(traj)=    unlist(lapply(datasH, function(dh) dh$name()))
    comb_plot=.merge1_new(traj, addName="data")
    #plot_traj(comb_plot)
    plot_traj(comb_plot, y="cumulative",keep_best = 5,txtsize=8,step=10,facet="data~maxsig")
  
    },
  
  #' selection for a single fold
  #' @param datasH a list of dataH objects
  #' @param phens list of phenotypes
  #' @param flags list of options
  #' @param k1  the fold
  #' @param vars_l_todo an object representing what is left to do
  #' @returns vars_l_todo object
  select_k=function(datasH,phens,flags, k1,
                    vars_l_todo 
                  ){
    verbose = .readFlag(flags,"verbose",FALSE);
    if(verbose)print(paste("fold",k1))
    expt_id = private$expt_id();
    # get_plots=.readFlag(flags, "get_plots",FALSE) 
    stop_y = .readFlag(flags, 'stop_y',"rand")
    logpvthresh = log(.readFlag(flags,"pthresh",0.1))
    beam= .readFlag(flags,"beam",1)
    saveAngles=FALSE
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
        comb_new  = dh$multiAnglesAndPv(comb20 , phens, k1,flags,expt_id, vars_l_todo)
        comb21 =  lapply(comb_new, function(x) x$pvs)
        private$savePvals(flags,phens,k1, dh$name(), vars_l_todo$vars_l,comb21)# no need to save here, just keep
        comb21
      }))
      if(inherits(comb2_news,"try-error")) break;
      
      vars_l_todo_new= private$nextVars(flags,phens, vars_l_todo,  k1,logpvthresh,beam, 
                                        comb2_news=if(useDB) NULL else comb2_news, 
                                        stop_y = stop_y, verbose=verbose)
      if(length(vars_l_todo$todo1)==length(vars_l_todo_new$todo1)){
        comb2 = comb2_news
      }
      vars_l_todo = vars_l_todo_new
      
      nvar = length(vars_l_todo$vars_l[[1]]$var)
      if(verbose) print(names(vars_l_todo$vars_l))
      
      if(length(vars_l_todo$vars_l[[1]]$var_names)>=flags$max  ) break;
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
    quantiles = sort(fromJSON(.readFlag(flags, "quantiles","[0]")),decreasing=TRUE)
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
        stop=FALSE)
      
    
    vars_l_todo
  },
  
  
  #' Clear the databases.  Only necessary if using a database and re-running analyses
  #'
  #' @param drop completely drop tables and start from scratch
  #' @param exclude tables to exclue from clearing
  #' @param datasH  list of dataH objects to clear
   clear_db=function(drop=FALSE,exclude="vars", datasH = NULL){
    if(drop){
      if(!is.null(private$sigs))private$sigs$drop_all(exclude=exclude)
      if(!is.null(datasH)){
        lapply(datasH, function(dh) dh$clear_db(drop=drop, exclude=exclude))
      }
    }else{
      warning("need to set drop=TRUE if you are sure, this will delete all saved signatures")
    }
    
  },  
  
 
  #' Get the next vars in iteration.  Internal function
  #'
  #' @param flags list of options
  #' @param phens phenotypes
  #' @param vars_l_todo vars_l_todo
  #' @param comb2_new results
  #' @param data_nme name of dataset
  #' @param k1 which repetition
  #' @return object outlining what is left to do

 savePvalsAndNextVars=function(flags,phens, vars_l_todo,comb2_new,data_nme, k1){
   beam=.readFlag(flags,'beam',1);
   stop_y=.readFlag(flags,"stop_y","rand")
   verbose=.readFlag(flags,"verbose",FALSE);
   logpvthresh = log(.readFlag(flags,'pthresh',0.05));
   angles_only = .readFlag(flags,'angles_only',FALSE);
   useDB = !is.null(private$sigs)
   if(angles_only) logpvthresh =0;
    if(useDB)  private$savePvals(flags,phens,k1, data_nme, vars_l_todo$vars_l,comb2_new)
   vars_l_todo_new= private$nextVars(flags,phens, vars_l_todo,  k1,logpvthresh,beam, 
                                     comb2_news=if(useDB) NULL else comb2_new, 
                                     stop_y = stop_y, verbose=verbose)
   vars_l_todo_new
   # private$nextVars(flags,phens, vars_l_todo,  k1,logpvthresh,beam, stop_y = stop_y, verbose=verbose)
 },
 
 
 #' main function for variable selection
 #' @param datasH a list of dataH objects
 #' @param phens list of phenotypes
 #' @param flags list of options
 #' @param transform_x transformation object 
 #' @param useDB boolean to indicate if results should be saved to database
 select=function(datasH, flags,
                 transform_x, 
                 phens=datasH[[1]]$pheno()$all,
                 useDB=FALSE){#c(y="function(y) y","function(y) y")
   phens1 = phens
   super$updateExpt( phens1, flags)
   #mc.cores = .readFlag(flags, "mc.cores",1)
   ##if(mc.cores>1 && )
   flags$transform_x = toJSON(transform_x);
   verbose=.readFlag(flags,'verbose',FALSE);
   if(flags$topn<flags$beam) stop("beam should be less than topn")
   if(is.null(flags[['data_types']]) || flags[['data_types']]=="{}")flags[['data_types']]=toJSON(datasH[[1]]$data_types())
   
   nreps_all = lapply(datasH, function(dh){
     dh$update(phens1, flags, transform_x);
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
