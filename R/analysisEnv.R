random_identifier = "rand" #\\.rand\\."  ## this identifies random variables;


#' Extracts all the variables selected from the variables object. Useful when using beam>1 
#' @param variables a variables object
#' @export
extract_all_variables<-function(variables){
df2 = .merge1_new(lapply(variables, function(v){
    .merge1_new(v)
}), addName="nrep")
return(tibble::as_tibble(df2))  
}



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
fspls.iterative<-function(dataset,flags, transform_x){
  flags$keep_all = FALSE;
  options(flags);
  if(is.null(dataset$y) || length(dataset$dataset)==0 || is.null(dataset$certainty)) stop(" dataset not well defined")
  if(!is.factor(dataset$y[[1]])) stop("y should be a factor")
  flags$nfold=1; flags$verbose=FALSE;
  mult = .readFlag(flags, "mult",100)
  dh = dataH$new(dataset$dataset, y = dataset$y, 
                 certainty = dataset$certainty,
                 nme="iterative", flags=flags)
 vars_all = NULL; all_models = NULL; preds1 = NULL;
 updates = list();
 select_each_iteration = .readFlag(flags, "select_each_iteration",FALSE);
 datasH = list(dh);
 k=0
 phens = dh$pheno()$all

 analysis =analysisEnv$new(flags=flags, dbDir=NULL) ;
 quantiles = list(); 
 sumdiffs = c();
  repeat  {
    {
    y_new =     dh$getYNew();
    q1 = c(y_new$error_rate, quantile(y_new$certainty[y_new$na_inds]));
    names(q1)[[1]] = "error"
    quantiles[[k+1]] =q1
    names(quantiles) = 1:(k+1)
   
  # print(data.frame(quantiles))
    }
    
    if(k==0 || select_each_iteration) {
            
          vars_all = analysis$select(datasH, phens, transform_x)
      
              #vars_all = fspls.select(datasH,  flags, transform_x, phens = dh$pheno()$all)
    }
    k = k+1
   all_models =dh$makeAllModels(vars_all)
   updated =  dh$updateWeights(all_models)
     updates[[k]] = updated
   message(paste("sumdiff",k,abs(updated$sumdiff)))
    sumdiffs[[k]] = abs(updated$sumdiff)
  if(abs(updated$sumdiff)<=1e-5 && k>2) {
        message("breaking since probs ordered not improving")
        break;
        
  }
    
    if(length(sumdiffs)>2 && sumdiffs[[k]] >= sumdiffs[[k-1]]) {
      message("breaking since sumdiff not changing");
      break;
      
    }
  }
 
 {
   y_new =     dh$getYNew();
   q1 = c(y_new$error_rate, quantile(y_new$certainty[y_new$na_inds]));
   names(q1)[[1]] = "error"
   quantiles[[k+1]] =q1
   names(quantiles) = 1:(k+1)
   
 }
  
   list(dh = dh, vars_all = vars_all, all_models = all_models,updates = updates, quantiles = quantiles,
        y_new = y_new$y, certainty_new = y_new$certainty,error_rate = y_new$error_rate
        
         )
#  plot(roc1)
  
}


#' main function for variable selection
#' @param datasH a list of dataH objects
#' @param flags list of options
#' @param transform_x transformation object 
#' @param phens list of phenotypes
#' @param data_types a data type object
#' @param dbDir database dir, can be NULL
#' @export
fspls.select<-function(datasH, flags,
                transform_x, 
                phens=datasH[[1]]$pheno()$all,
                data_types=datasH[[1]]$data_types(),
                                dbDir = NULL
               ){#c(y="function(y) y","function(y) y")
  flags$keep_all = FALSE;
  options(flags);
  analysis =analysisEnv$new(flags=flags, dbDir=dbDir) ;
  vars_all = analysis$select( datasH,phens,transform_x, data_types)
 
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
.mergeResInner<-function(res_inner1){ ## still used? 
  nme_comb = names(res_inner1[[1]]);names(nme_comb) = nme_comb
#  nme_c1 = nme_comb[[1]]; nme_p1 = names(res_inner1[[1]][[nme_c1]])[1] ; ri = res_inner1[[1]]
  
  nxt_vars = .merge1_new(lapply(nme_comb, function(nme_c1){
    nmesp1 = names(res_inner1[[1]][[nme_c1]]);names(nmesp1) = nmesp1
    .merge1_new(lapply(nmesp1, function(nme_p1){
      r_comb = .merge1_new(lapply(res_inner1, function(ri){
        ri[[nme_c1]][[nme_p1]] #|> unite("name1",data_type, names, remove=F,sep=".")
      }), addName="dataset")
      r_comb1 = r_comb |> pivot_wider( id_cols=c(names, data_type), names_from="dataset", values_from=value) #
      r_comb1[,1:2] |> tibble::add_column(value= apply(r_comb1[,-(1:2)],1,sum))
      
      }), addName="param")
      # cumpv_all= .sumChisq(unlist(lapply(nv, function(nv1)unlist(nv1$pvs_all))))
  
    }), addName="transf")
  nxt_vars 
 
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
    means_y=function(datasH){ ## averages over the datasH
      means_y_all =lapply(datasH, function(dh) dh$means_y());
      counts_y_all = lapply(datasH, function(dh) dh$counts_y());
      #means_y = list();
    #  means_y[[colk]][[f_k]][[g_k]]
      nreps = datasH[[1]]$nreps();
      nmes = names(means_y_all); names(nmes) = nmes
  # l1 = 1
      res=lapply(1:2, function(jk){  ##iterates over means_y0 and means_y1
        #print(jk)
          means_y_comb = lapply(nreps, function(k){
          #  print(k)
            means_y = means_y_all[[1]][[k]][[jk]]
            
           nmes1=names(means_y); names(nmes1)=nmes1;
          # nme1  = nmes1[[1]];  nmes2 = names(means_y[[nme1]]); nme2 = nmes2[[1]];  nmes3 = names(means_y[[nme1]][[nme2]]);  nme3 = nmes3[[1]]
      
          lapply(nmes1, function(nme1){
             nmes2 = names(means_y[[nme1]]); names(nmes2) =nmes2
              lapply(nmes2, function(nme2){
               
                  nmes3 = names(means_y[[nme1]][[nme2]]); names(nmes3) = nmes3;
                  lapply(nmes3, function(nme3){
                  df_mean = data.frame(lapply(nmes, function(nme){
                    means_y_all[[nme]][[k]][[jk]][[nme1]][[nme2]][[nme3]]
                    
                  })  )
                  df_counts = data.frame(lapply(nmes, function(nme){
                    counts_y_all[[nme]][[k]][[nme1]][[nme2]][[nme3]]
                  })  )
                  weights = (apply(df_counts,1,function(v) v/sum(v)))
                  apply(df_mean * weights,1,sum)
                 
                })
              })
            })
         })
      })
     
      res
    },
    topAngles = function(vars_l_todo, k1, comb2_news, verbose=FALSE,stop_y=random_identifier,minvars =0){
      vars_l = vars_l_todo$vars_l
      angles_all = .merge1_new(comb2_news, addName="dataset") 
      angles_all =   pivot_wider(angles_all, id_cols=c(partialname, model, varname, signature,names, data, transform, param,prev), names_from="dataset", values_from=value)
      angles_all = angles_all[,1:9] |> tibble::add_column(value= apply(angles_all[,-(1:9)],1,sum))
       o = order(angles_all$value, decreasing=FALSE)
      ang1=angles_all[o,] 
      #unite("partialname", data_type, names, transf, sep=".", remove=FALSE) |> unite("varname", data_type,names, transf,param,sep="." ,remove=FALSE)  |> unite("model",prev,varname,sep=";", remove=FALSE) 
     # ang1$model =  gsub("^empty;", "",ang1$model)

    
      ang1 = subset(ang1,!duplicated(ang1$partialname))
      
     
    },
   
    savePvals=function(k1, data_nme, vars_l, comb2){
      flags = private$flags;
      phens = private$phens;
      for(varn1 in names(vars_l)){ 
        ri=comb2[[varn1]]
        varnames=vars_l[[varn1]]$var_names
        super$savePvals(ri, varnames, k1, data_nme=data_nme, useCurrVarnames=FALSE)

      }
    },
    saveAngles=function(k1, data_nme, vars_l, comb_){
      for(varn1 in names(vars_l)){ 
        private$saveAngles( data_nme, comb_[[varn1]], vars_l[[varn1]]$var_names,k1)
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
  #' @param k1  the fold
  #' @param vars_l_todo an object representing what is left to do
  #' @returns vars_l_todo object
  select_k=function(datasH, k1,
                   means_y_k = NULL
                                 ){
    force=FALSE
    vars_l_todo = self$getTodo(private$flags, private$phens);
    
    invisible(lapply(datasH, function(dh){
      dh$updateTrain(k1,means_y_k = means_y_k, force=force);
    }))
    phens = private$phens;
    flags = private$flags;
    verbose = .readFlag(flags,"verbose",FALSE);
    if(verbose)print(paste("fold",k1))
    expt_id = private$expt_id();
    maxvars = .readFlag(flags,"max",1e6)
    minvars = .readFlag(flags,"min",0)
    stop_y = random_identifier
    logpvthresh = log(.readFlag(flags,"pthresh",0.5))
    beam= .readFlag(flags,"beam",1)
    topn = .readFlag(flags,'topn',20) ## number to calculat pvals for
    useAngles = .readFlag(flags,"angles_only", TRUE)
    
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
    vars_all = list();
    while(length(vars_l_todo$todo1)>0 ){
      comb2_news =try(lapply(nmesH, function(nmeh){
        dh = datasH[[nmeh]]
        comb20=comb2[[nmeh]]
        dh$multiAnglesAndPv(comb20 , k1,expt_id, vars_l_todo)
       
      }))
      if(inherits(comb2_news,"try-error")){
        print(comb2_news) #error
        warning("error so not continuing to extend model")
        break;
      } 
      if(!useAngles){
        ang1 = private$topAngles( vars_l_todo, k1,comb2_news, verbose=verbose );
        rand_ind=which(ang1$transform=="rand")
        non_rand_ind = which(ang1$transform!="rand")
        if(length(non_rand_ind)==0){
          if(nvar <minvars) warning("no non-random indices, breaking. You can try to increase topn or useAngles=TRUE")
          break;
        }
        if(nvar < minvars && min(rand_ind)==1){
          rand_ind = rand_ind[rand_ind>min(non_rand_ind)]
        }
          minv = min(nrow(ang1), topn, rand_ind)
          ang1 = head(ang1,minv)
        
          comb2_filtered = lapply(comb2_news, function(xx){
            xx[match(ang1$model,xx$model),]
#            subset(xx, model %in% ang1$model)
          })
          comb2_news =lapply(nmesH, function(nmeh){
            dh = datasH[[nmeh]]
            comb_filtered=comb2_filtered[[nmeh]];
            comb_new=dh$calcPvs(comb_filtered, k1, expt_id)
           comb_new
          })
       #   ang1 = private$topAngles( vars_l_todo, k1,comb2_new=if(useDB) NULL else comb2_news, verbose=verbose );
       #   cum_pvs = unlist(lapply(ang1, function(x) x$cum_pv))
          #ang2 = ang2[order(cum_pvs)]
        #  if(verbose)print(head(cum_pvs))
          
      }
    # grep(random_identifier, names(ang2))
      
      vars_l_todo_new= self$nextVars( vars_l_todo,  k1,logpvthresh,beam, 
                                        comb2_news,
                                      #  ang1=ang1,
                                        
                                        stop_y = if(nvar<minvars) NULL else stop_y, 
                                      verbose=verbose)
      if(verbose) print(vars_l_todo_new$vars_l)
      if(verbose) print(paste("stop",vars_l_todo_new$stop, "todo",length(vars_l_todo_new$todo1)))
      if(vars_l_todo_new$moveNext){
       
        comb2 = comb2_news
        
      }
     
      vars_l_todo = vars_l_todo_new
      
      nvar = vars_l_todo$vars_l$nvar[[1]]
      
      if(is.na(nvar)) break;
      
      if(verbose) print(names(vars_l_todo$vars_l))
      if(nvar>0) vars_all[[nvar]] = vars_l_todo$vars_l
      if(nvar>=maxvars ) break;
    }
    if(length(vars_all)>0) names(vars_all) = paste("nvar",1:length(vars_all),sep="_")
     return(vars_all)  ## this keeps the vars from each stage of iteration
    #vars_l_todo$vars_l
  },
  
  #' Get todo 
  #'
  #' @param flags list of options 
  #' @param phens phenotypes
  #' @return object outlining what is left to do for next iteration
  getTodo=function(flags, phens){
  
    if(is.null(phens)) stop("!!")
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
    
    # var_thresh = var_thresh[match(names(var_thresh), train_nme)]
    #empty=initialStateObj(phens)
    
    empty = data.frame(list(partialname="empty", model="empty", sumAngle=0, varname="empty", signature="", names="empty", data="", transform="", param="", prev="", value=0, nvar=0, beam=1))
    
    vars_l_todo =
      list(
        todo1 = todo1,
        jj=0,
        logpv=logpv,
        vars_l = empty,
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
  #' @param vars_l_todo vars_l_todo
  #' @param comb2_new results
  #' @param data_nme name of dataset
  #' @param k1 which repetition
  #' @return object outlining what is left to do

 savePvalsAndNextVars=function(vars_l_todo,comb2_new,data_nme, k1){
   flags = private$flags;
 #  phens = private$phens;
   beam=.readFlag(flags,'beam',1);
   stop_y=random_identifier #.readFlag(flags,"stop_y",random_identifier)
   verbose=.readFlag(flags,"verbose",FALSE);
   logpvthresh = log(.readFlag(flags,'pthresh',0.05));
   angles_only = .readFlag(flags,'angles_only',TRUE);
   verbose=.readFlag(flags,'verbose',TRUE);
   useDB = !is.null(private$sigs)
   if(angles_only) logpvthresh =0;
    if(useDB)  private$savePvals(k1, data_nme, vars_l_todo$vars_l,comb2_new)
   #private$topAngles(vars_l_todo, k1, comb2_new,v)
   ang1 = private$topAngles(vars_l_todo, k1, comb2_new, verbose=verbose) ###
   return(ang1);
  # vars_l_todo_new= private$nextVars(vars_l_todo,  k1,logpvthresh,beam, 
  #                                   comb2_news=if(useDB) NULL else comb2_new, 
  #                                   stop_y = stop_y, verbose=verbose)
  # vars_l_todo_new
   # private$nextVars(flags,phens, vars_l_todo,  k1,logpvthresh,beam, stop_y = stop_y, verbose=verbose)
 },
 
 #' Get the next vars in iteration.  Internal function
 #'
 #' @param vars_l_todo vars_l_todo
 #' @param comb2_news results
 #' @param data_nme name of dataset
 #' @param k1 which repetition
 #' @param stop_y should be random_identifier
 #' @return object outlining what is left to do
 nextVars=function(vars_l_todo,  k1,
                   logpvthresh,beam,  comb2_news,stop_y=random_identifier, verbose=FALSE
                                  ){
#   value_string = if(useAngles) "value" else "pvals"
   
   angles_all = .merge1_new(comb2_news, addName="dataset") 
   angles_all$value = angles_all$value+angles_all$sumAngle
   #subset(angles_all, model=="log_pre.ABCD2.pow.1;log_pre.AL391280_2.pow.1")
  # angles_all |> dplyr::summarise(n = dplyr::n(), .by = c(partialname, model, varname, signature, names, data_type, transf, param, prev, nvar, dataset)) |>
  #   dplyr::filter(n > 1L)
   
    angles_all =  pivot_wider( angles_all ,id_cols=c(partialname, model, varname,names, data,  signature,transform, param,prev, nvar), names_from="dataset", values_from=value)
   angles_all = angles_all[,1:10] |> tibble::add_column(value= apply(angles_all[,-(1:10)],1,sum))
   o = order(angles_all$value, decreasing=FALSE)
   ang1=angles_all[o,] 
   
  # ang1 = private$topAngles(vars_l_todo, k1, comb2_news, verbose=verbose)
  # angles_ = unlist(lapply(ang1, function(a1)a1[["angle"]]))
   flags=private$flags; 
   vars_l = vars_l_todo$vars_l
   todo1 = vars_l_todo$todo1
   #logpvs =if(useAngles) angles_ else  -1* unlist(lapply(ang1, function(a1)a1[["cum_ll"]]))
   #angles_all=unlist(lapply(ang1, function(a1)a1[["cum_angle"]]))
   #logpvs_all = if(useAngles) angles_all else  -1*unlist(lapply(ang1, function(a1)a1[["cum_ll"]]))
  # ord = order(logpvs)
  # names(ord) = names(logpvs)
   #ord_all = order(logpvs_all)
   #ang1 = ang1[ord_all]
   #logpvs = logpvs[ord_all]
   #logpvs_all = logpvs_all[ord_all]
   #angles_ = angles_[ord_all]
   #angles_all = angles_all[ord_all]
   stop_random1 = FALSE; stop_random2 = FALSE;
   if(!is.null(stop_y)){
     gp1=grep(stop_y,ang1$transform)
     gp=grep(stop_y, ang1$transform, invert=TRUE)
     
     
       stop_random1 = min(gp1)<min(gp); stop_random2 = stop_random1
    
     
     
     
     # stop_random= min(gp1)<=min(gp)
     
   
     if( verbose) {
       #print(unlist(list(rand=min(logpvs[gp1]),nonrand= min(logpvs[gp]))))  #verbose
       print((list(rand=ang1[min(gp1),],nonrand= ang1[min(gp),]))) #verbose
     }
     
   }
   #logpv =min(logpvs)
  # if(stop_random || stop_random2 || stop_random1){
  #   if(verbose) print(paste("stopping due to random", exp(logpv), names(logpvs)[which.min(logpvs)]))
  # }
 #  if(logpv>logpvthresh){
#     if(verbose) print(paste("stopping due to pvalue", exp(logpv), exp(logpvthresh)))
     
#   }
   
   ##ADD MORE RESTRICTIONS .. eg maxsize
   #     while( (length(vars_l[[1]]$var_names) < minsize || logpv<logpvthresh) && length(vars_l[[1]]$var_names)<maxsize && ! vars_l_todo$stop_random){
   if(( !stop_random2 && !stop_random1  ) ){
    # if(verbose){
    #   print(head(sort(logpvs_all[gp]),beam))
     #  print(names(vars_l))
    # }
     ##remove signatures which are same in different order
     #non_rand = grep(random_identifier, names(ang1),inv=TRUE)
     #dupls = rep(FALSE, length(ang1))
     
     exclusion_in_beam = strsplit(as.character(getOption("exclusion_in_beam","0")),":")[[1]]  ## can be a number of max
     if(exclusion_in_beam[1]=="max" && beam>1){  ## excludes any overlapping signature
       #### REDO
       
    # duplicated=which(duplicated(unlist(vn_all[inds_to_keep])))   
    # if(length(duplicated)>0) stop("not supposed to have duplicates") 
      
     }else if(exclusion_in_beam[1] !="max" && !is.na(as.numeric(exclusion_in_beam[1])) && as.numeric(exclusion_in_beam[1])>0 && beam>1){  ##KEEPS DIFFERENT N-1 signatures
       
     }
    # last_non_rand = grep(random_identifier, names(ang1))[1]-1
    # if(is.na(last_non_rand)) last_non_rand = beam;     
  #   topn = .readFlag(flags,'topn', 20);
   #  num_pvals = min(topn, 20)
     ang1 = subset(ang1,transform!="rand")
     tokeep=1:min(nrow(ang1), beam)
     ang1 = ang1[tokeep,,drop=F] |> tibble::add_column(beam = 1:length(tokeep))#last_non_rand,
     
     
     todo1 = vars_l_todo$todo1;
     todo1[[1]]$incl$nvar = todo1[[1]]$incl$nvar+1;
     if(todo1[[1]]$incl$nvar>= todo1[[1]]$incl$max){
       todo1 = todo1[-1]
     }
     #ang1$sumAngle = ang1$sumAngle+ang1$value
     vars_l_todo = list(stop=FALSE, vars_l = ang1, todo1 = todo1, moveNext=TRUE)
     #if(length(grep(random_identifier, names(vars_l_todo$vars_l)))>0) stop("!!");
     
     return(vars_l_todo)
   }
   todo1=vars_l_todo$todo1[-1]
   vars_l_todo_new = list( stop=length(todo1)==0, vars_l = vars_l, todo1 = todo1, jj = vars_l_todo$jj+1, moveNext=FALSE)
   return(vars_l_todo_new)
 },
 
 #' main function for variable selection
 #' @param datasH a list of dataH objects
 #' @param phens list of phenotypes
 #' @param transform_x transformation object
 #' @param data_types data_types object
 select=function(datasH,phens,transform_x, data_types=datasH[[1]]$data_types(), avg_mean_y = FALSE){#c(y="function(y) y","function(y) y")
   flags = private$flags;
   flags = super$updateExpt( phens, flags, transform_x, data_types)
   verbose=.readFlag(flags,'verbose',FALSE);
   #if(flags$topn<flags$beam) stop("beam should be less than topn")
   
   nreps_all = lapply(datasH, function(dh){
     dh$update(phens=phens, flags=flags, transform_x=transform_x, data_types=data_types);
     dh$nreps();
   })
   nreps = nreps_all[[1]]
   
   
     vars_all = super$loadVars()
     if(!is.null(vars_all)) return(vars_all)
  
    means_y = if(avg_mean_y) private$means_y(datasH) else NULL;
   
   variables1 <- lapply(nreps, function(k1) {  ## can use mclapply here
     
     means_y_k =if(avg_mean_y) means_y[[k1]] else NULL
      self$select_k(datasH,k1,   means_y_k =means_y_k)
     
   })
   
#   vars_all=super$integrate(variables1)

   private$saveVars(variables1)
   variables1
   
 }
 


  


))
