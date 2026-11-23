random_identifier = "rand" #\\.rand\\."  ## this identifies random variables;

#' main function for variable selection
#' @param min minimum number of variables
#' @param max maximum number of variables, number of selected features 
#' @param nfold number of folds.  10 fold implies a 10 non-overlapping 10% splits
#' @param batchsize batch size can be specified alternatively to nfold, but only one can be nonzero
#' @param angles_only whether to use purely angles to find best variables (FALSE faster but less accurate)
#' @param topn if using model fit to rank variables, how many top angles to take  (larger takes more time)
#' @param beam  how many best combinations to take through to next round of iteration. minimum value is 1, but higher explores more combinations
#' @param all_v_all  multi-class outcomes, this builds all vs all classifier, instead of a one vs all #here!!!
#' @param project whether to use the projection approach in model building (TRUE recommended)
#' @param useglmnet whether to use glmnet for model fitting (after variables selected) doesnt work well for randomisation.  Note that useglmnet is TRUE by default for multinomial
#' @param verbose print verbose
#' @param show_warnings  print extra warnings
#' @param binomial.link link parameter for binomial, one is default
#' @param gaussian.link  link parameter for gaussian, zero is default
#' @param fspls.types  a list indicate the output evaluations to consider for different families
#' @param firth=TRUE,
#' @export
getFlags<-function(min=1,  #minimum number of variables
         max=100, #maximum number of variables, number of selected features 
         nfold=1,#number of folds.  10 fold implies a 10 non-overlapping 10% splits
         batchsize=0, #batch size can be specified alternatively to nfold, but only one can be nonzero
         angles_only=TRUE, #whether to use purely angles to find best variables (RECOMMENDED FALSE)
         topn=100,  #if using model fit to rank variables, how many top angles to take  (larger takes more time)
         beam=1, #how many best combinations to take through to next round of iteration. minimum value is 1, but higher explores more combinations
         all_v_all=TRUE,  #for multi-class outcomes, this builds all vs all classifier, instead of a one vs all #here!!!
         project=TRUE,   #whether to use the projection approach in model building (TRUE recommended)
         useglmnet=FALSE,  #whether to use glmnet for model fitting (after variables selected) doesnt work well for randomisation.  Note that useglmnet is TRUE by default for multinomial
         verbose=TRUE,
         lambda = NULL, #for using glmnet
         show_warnings=FALSE,
         binomial.link= c(0.1, 0.3,0.5,1.0),
         gaussian.link = c(0),
         fspls.types= list(gaussian = c("correlation","rms"),binomial=c("AUC","AUPRC","area"),
                              multinomial=c("AUC"),
                              ordinal = c("AUC_all")),
         
         firth=TRUE){
      
   list(min=min,  #minimum number of variables
                 max=max, #maximum number of variables, number of selected features 
                 nfold=nfold,#number of folds.  10 fold implies a 10 non-overlapping 10% splits
                 batchsize=batchsize, #batch size can be specified alternatively to nfold, but only one can be nonzero
                 angles_only=angles_only, #whether to use purely angles to find best variables (RECOMMENDED FALSE)
                 topn=topn,  #if using model fit to rank variables, how many top angles to take  (larger takes more time)
                 beam=beam, #how many best combinations to take through to next round of iteration. minimum value is 1, but higher explores more combinations
                 all_v_all=all_v_all,  #for multi-class outcomes, this builds all vs all classifier, instead of a one vs all #here!!!
                 project=project,   #whether to use the projection approach in model building (TRUE recommended)
                 useglmnet=useglmnet,  #whether to use glmnet for model fitting (after variables selected) doesnt work well for randomisation.  Note that useglmnet is TRUE by default for multinomial
                 verbose=verbose,
                 lambda = lambda, #for using glmnet
                 show_warnings=show_warnings,
                 binomial.link= binomial.link,
                 gaussian.link = gaussian.link,
                 fspls.types = fspls.types,
                 firth=firth
    )
}
#' @description
#' this reads pre-saved results and makes a new analysis object as well as updating datasH objects
#' it re-establishes the variables selected in analysis and the models in dataH objects
#' @param out_dir directory with results saved via save_all_results
#' @param datasH list of dataH objects to update
#' @export
read_all_results<-function(out_dir, datasH = NULL){
  variables0=jsonlite::read_json(paste(out_dir,"variables.json",sep="/"), simplifyVector=TRUE)
  modelfile = paste(out_dir,"all_models.json",sep="/")
  all_models0 = jsonlite::unserializeJSON(paste(readLines(modelfile), collapse="\n"))

  folds = jsonlite::read_json(paste(out_dir,"folds.json",sep="/"), simplifyVector= TRUE);
  params = jsonlite::read_json(paste(out_dir,"params.json",sep="/"), simplifyVector = TRUE)
  analysis = analysisEnv$new(flags=params$flags, dbDir=NULL)
  analysis$variables = variables0
  if(!is.null(datasH)){
    
    for(nme in names(datasH)) {
      dh = datasH[[nme]];
      dh$variables = variables0
      dh$update(phens = params$phens, flags = params$flags, transform_x = params$transform_x, data_types = params$data_types, 
                                folds = folds, 
                                force=TRUE)
                    dh$all_models = all_models0[[nme]]
    }
  }
  analysis
  }


plotInternal<-function(eval5,sep_by,  title1, point, dotsize, shape_nme, color_nme, 
                       linetype_nme, text_nme,  labelsize, showranges,
                       grid0_nme, grid1_nme, scales, txtsize,  line, transformy, showtext, legend){
  
  if(length(sep_by)==0){
    # ph3 = paste(sort(unique(apply(eval5[,names(eval5) %in% title1, drop=FALSE],1,paste,collapse=","))), collapse=" ")
    ggp<-ggplot(eval5);
    if(point) {
      if(dotsize %in% names(eval5)){
        ggp<-ggp+geom_point(aes_string(x="numvars", y="mid",  shape=shape_nme,size=dotsize, color=color_nme))
      }else{
        ggp<-ggp+geom_point(aes_string(x="numvars", y="mid",  shape=shape_nme, color=color_nme),size=dotsize)
        
      }
    }
    if(line)  ggp<-ggp+geom_line(aes_string(x="numvars", y="mid", linetype=linetype_nme,  color=color_nme)) #+ggtitle(paste(title,ph, ph3))
    if(showtext) ggp<-ggp+geom_text_repel(aes_string(x="numvars", y="mid", label=text_nme,color=color_nme),size=labelsize,max.overlaps = 1000) #+ggtitle(paste(title,ph, ph3))
    
    if(showranges){ ## geom_ribbon vs geom_errorbar
      
      ggp<-ggp+ geom_ribbon(aes_string(x = "numvars", ymin="low", ymax="high",linetype=linetype_nme,color=color_nme, fill = color_nme ), alpha = 0.1)
    }
    if(nchar(grid0_nme)>0){
      if(nchar(grid1_nme)>0){
        ggp<-ggp+facet_grid(paste(grid0_nme, grid1_nme,sep="~"),scales=scales)
      }else{
        ggp<-ggp+facet_wrap(grid0_nme, scales=scales)
      }
    }
    legend_position=if(legend) "bottom" else "none";
    ggp<- ggp+ theme(legend.position = legend_position,legend.title = element_text(size = txtsize));#+theme(,    legend.text = element_text(size = 3))
    ggp<-ggp+ scale_y_continuous(transform=transformy);#scale_y_log10() 
    return(ggp+ggtitle(title1))
    
  }
  #print(sep_by)
  levs1 = unique(eval5[[sep_by[[1]]]]); names(levs1) = levs1;
  final_res=lapply(levs1,function(ph){
    subinds1 = eval5[[sep_by[[1]]]]==ph
    
    if(length(sep_by)>0){
      res =plotInternal(eval5[subinds1,,drop=F], sep_by[-1], 
                        paste(title1,paste(sep_by[[1]],ph, sep="=")), 
                        point, dotsize, shape_nme, color_nme,
                        linetype_nme, text_nme, labelsize, showranges,
                        grid0_nme, grid1_nme, scales, txtsize, line, transformy, showtext, legend);
      return(res);
    }
    
    
    
    
  })
  final_res
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
      angles_all = angles_all[,1:9] |> tibble::add_column(value= apply(angles_all[,-(1:9)],1,sum,na.rm=T))
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
              flags = getFlags(),
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
  #' @description extractPredictions from datasH objects
  #' @param datasH
  extractPredictions=function(datasH, prob=TRUE){
    predictions = lapply(datasH, function(dh) dh$extractPredictions(prob=prob)) 
    self$predictions = predictions#  names(datasH), function(nmeh) datasH[[nmeh]]$extractPredictions(prob=TRUE))
    invisible(predictions);
  },
  plotPredictions = function(datasH, show="area"){
    pred_plots = lapply(datasH, function(dh) dh$plotPredictions(show="area"))
    self$pred_plots = pred_plots
    invisible(pred_plots);
  },
  
  #' @description saves all results to a directory
  #' @param datasH a list of dataH objects 
  #' @param out_dir  output directory
  #' @param savePlots whether to save plots
  #' @returns assigned name of dataset
  save_all_results=function(
    datasH,
    out_dir =paste("FSPLS", gsub("\\s+","_",date()), sep="_" ),
    savePlots=TRUE, width = 20 , height = 20
  ){
    working_dir = getwd();
    dir.create(out_dir,recursive = TRUE);
    setwd(out_dir);
     if(!is.null(self$evals_pretty)) write_xlsx(self$evals_pretty, "evaluation.xlsx");
    betas_table = unlist(self$betas_table, recursive=FALSE)
   
    toprint = c(list(variables=self$vars_table), betas_table)
    toprint = toprint[unlist(lapply(toprint, length))>0]
    if(length(toprint)>0){
      write_xlsx(toprint, "models.xlsx");
    }
     toDUMP = list(
      flags = private$flags, phens=private$phens, transform_x = private$transform_x, data_types = private$data_types)
    
    jsonlite::write_json(toDUMP,"params.json", simplifyVector=TRUE, pretty = TRUE);
    
    if(!is.null(self$variables)) jsonlite::write_json(self$variables,"variables.json", simplifyVector=TRUE, pretty = TRUE);
    # if(!is.null(self$all_models)) jsonlite::write_json(self$all_models,"all_models.json", simplifyVector=TRUE, pretty = TRUE);
    if(length(datasH)>0){
      all_models = lapply(datasH, function(dh)dh$all_models)
      folds = lapply(datasH, function(dh)dh$folds())
      jsonlite::write_json(folds, "folds.json", simplifyVector=TRUE, pretty = TRUE);
  
     writeLines(jsonlite::serializeJSON(all_models), "all_models.json", useBytes = TRUE)
    }
    if(savePlots){
      .print("data_plot.pdf", self$data_plots, width = width);
      .print("eval_plot.pdf", self$eval_plot, width = width);
      .print("pred_plots.pdf", self$pred_plots, width = width)
    }
    setwd(working_dir)
    out_dir
  },
  
  
  evaluate=function(datasH){
    evals = lapply(datasH, function(dh) dh$evaluate())
    eval_combined = .merge1_new(evals, addName="dataset")
    self$evals = eval_combined;
    invisible(eval_combined)
  },
  
  tabulateBetas=function(datasH, max_beam=10){
  betas_all0 = lapply(datasH, function (dh) dh$tabulateBetas(max_beam = max_beam))
  nvars = 1:length(betas_all0[[1]]); names(nvars) = nvars;
  betas_all = lapply(nvars, function(nv){
   beams = names(betas_all0[[1]][[nv]]); names(beams) = 1:length(beams)
    lapply(beams, function(beam){
      df1 = betas_all0[[1]][[nv]][[beam]]
      colind = 5:ncol(df1)
      df = data.frame(lapply(betas_all0, function(ba) ba[[nv]][[beam]][,5]))
       cbind(df1[,1:4],df)
   
   })
  })
  self$betas_table = betas_all
  betas_all
  },
  
  
  #' @description extracts and plots the selected variables
  #' @param datasH is a list of dataH objects
  #' @param violin violin plots
  #' @returns ggplots 
  plotData=function(datasH, variables = self$variables, facet="gene~pheno_data",flip  = TRUE,violin=TRUE){
    
    plot_data_all =lapply(datasH, function(dh) dh$plotData(variables ))
                          
                          families = names(plot_data_all[[1]]); names(families) = families
                          
                          plot_all = lapply(families, function(family){
                            df5 = .merge1_new(lapply(plot_data_all, function(pd)pd[[family]])) |> unite('pheno_data', pheno, dataset,sep="\n", remove=FALSE)
                            
                            fam = strsplit(family,"\\.")[[1]][[1]];
                            #df5=subset(plot_data, family==fam)
                            if(flip) facet ="pheno_data~gene";
                            names(df5) =sub("name","gene",names(df5))
                            
                            df5$gene = gsub("_"," ",df5$gene)
                            df5$gene = gsub(" ","\n",df5$gene)
                            
                            df5$gene = factor(df5$gene, levels = df5$gene[!duplicated(df5$gene)])
                            
                            df6=subset(df5, family=="gaussian" | family=="ordinal")
                            df7=(subset(df5, family!="gaussian"))
                            pv6 = getPvals(df6)
                            pv7 = if(fam=="binomial") getPvals(df7) else pv6;
                            # color= if(length(unique(df4$dataset))>1) "dataset" else #"pheno"
                            # facet= if(!is.null(df4[['transform']]) ) "transform~pheno" else "pheno"
                            
                            lens = c(length(table(df5$gene)), length(table(df5$pheno_data)))
                            wrap=TRUE;
                            title=""
                            if(lens[1]>1 & lens[2]>1){
                              
                              wrap = FALSE
                            }else if(lens[2]>1){
                              title= df5$gene[1]
                              facet="pheno_data"
                            }else{
                              title = df5$dataset[[1]]
                              facet="gene"
                            }
                            ggp = NULL; ggp1 = NULL;ggp2=NULL
                            if(nrow(df7)>0){
                              df7$y = factor(df7$y, ordered=TRUE)
                              
                              ggp<-ggplot(df7, aes(x=y, y=value, color=pheno));#+facet_wrap(facet, scales="free");#+ggtitle(unlist(phens1))
                              if(violin){
                                ggp<-ggp+ggplot2::geom_violin()+geom_point()
                              }else{
                                ggp<-ggp+ggplot2::geom_boxplot()
                              }
                              if(!is.null(pv7)){
                                pv7$y = factor(pv7$y, levels = levels(df7$y))
                                
                                ggp<-ggp+ggplot2::geom_text(data=pv7,
                                                            mapping = aes(x=y, y=value,  label=pv), color="black")
                              }
                            }
                            
                            if(nrow(df6)>0){
                              separator="______________"
                              prbs=c(0.33,0.5,0.66)
                              df6_1 =df6 |> unite("comb",gene,family,pheno_data,sep=separator)
                              comb1=unique(df6_1$comb); names(comb1)=comb1
                              quants=.merge1_new(lapply(comb1,function(c1){
                                df6_2 = subset(df6_1, comb==c1)
                                df6_2$y = as.numeric(as.character(df6_2$y))
                                df6_2 = df6_2[order(df6_2$y),]
                                yv = unique(df6_2$y)
                                names(yv)=yv
                                .merge1_new(lapply(yv, function(yv1){
                                  df6_3 = subset(df6_2, y<=yv1)
                                  q1=quantile(df6_3$value,na.rm=TRUE,probs = prbs)
                                  df_=data.frame(t(data.frame(q1)))
                                  df_
                                }),addName="y")
                              }), addName="comb")
                              
                              quants1 = quants |> separate("comb", sep=separator, into=c("gene","family","pheno_data"))
                              names(quants1) = gsub("\\.","",names(quants1))
                              quants1$y = as.numeric(quants1$y)
                              #gsub("X50", "median expression", names(quants1))
                              quants1$gene=factor(quants1$gene, levels = levels(df6$gene))
                              ggp1<-ggplot(quants1, aes(x=y, y=X50, ymin=X33, ymax = X66,color=pheno_data, fill=gene))
                              ggp1<-ggp1+geom_line()+geom_ribbon(alpha=0.1)
                              ggp1<-ggp1+ggplot2::geom_text(data=pv6,
                                                            mapping = aes(x=y, y=value, ymin=value, ymax=value, label=pv), color="black")
                              # ggp1<-ggp1+ggplot2::geom_text(mapping = aes(x=y, y=value, label=pv), color="black",data = pv6)
                              
                              ggp1<-ggp1+labs(y="median gene expression for y<=yval", x="yval")
                              
                              ggp2<-ggplot(df6, aes(x=y, y = value, color=pheno, fill=pheno, shape = dataset))+geom_point()
                              ggp2<-ggp2+ggplot2::geom_text(data=pv6,
                                                            mapping = aes(x=y, y=value, label=pv), color="black")
                              
                            }
                            if(wrap){
                              ggp<-ggp+facet_wrap(facet,scales = "free")+ggtitle(title)
                              ggp1<-ggp1+facet_wrap(facet,scales = "free")+ggtitle(title)
                              ggp2<-ggp2+facet_wrap(facet,scales = "free")+ggtitle(title)
                              
                            }else{
                              ggp<-ggp+facet_grid(facet,scales = "free")+ggtitle(title)
                              ggp1<-ggp1+facet_grid(facet,scales = "free")+ggtitle(title);#+scale_y_continuous(transform="log")
                              ggp2<-ggp2+facet_grid(facet,scales = "free")+ggtitle(title);
                            }
                            
                            resu =  list("binomial"=ggp,"gaussian"=ggp1, "gaussian1"=ggp2)
                            if(flip) resu = lapply(resu, function(r1)r1+ggplot2::coord_flip())
                            resu
                          })
                          
          self$data_plots=plot_all
          invisible(plot_all)
  },
  #' selection for a single fold
  #' @param datasH a list of dataH objects
  #' @param k1  the fold
  #' @param vars_l_todo an object representing what is left to do
  #' @returns vars_l_todo object
  select_k=function(datasH, k1,
                   means_y_k = NULL,
                   minvars = .readFlag(private$flags,"min",0),
                   maxvars = .readFlag(private$flags,"max",1e6)
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
      if(verbose)print(paste("nvar",nvar))
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
      
        
          comb2_filtered = .lapply_nme(comb2_news, function(xx_nme){
            xx = comb2_news[[xx_nme]]
            mi2=match(ang1$model,xx$model);
            
            if(length(which(is.na(mi2)))>0) {
              model_nme=ang1$model[which(is.na(mi2))]
              print(model_nme)
              stop(paste("problem", xx_nme))
            }
            xx[mi2,]
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
   angles_all = angles_all[,1:10] |> tibble::add_column(value= apply(angles_all[,-(1:10)],1,sum,na.rm=T))
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
    #full first
    variables1 =lapply(nreps[length(nreps)], function(k1) { 
      self$select_k(datasH,k1,   means_y_k =if(avg_mean_y) means_y[[k1]] else NULL)
      
    })
    maxvars = length(variables1[[1]]); minvars = maxvars;
    print(maxvars)
    if(length(nreps)>1){
          variables_cvs <- lapply(nreps[-length(nreps)], function(k1) { 
            print(paste("doing ",k1,"................"))
            ## can use mclapply here
            self$select_k(datasH,k1,   means_y_k =if(avg_mean_y) means_y[[k1]] else NULL, 
                          minvars = minvars, maxvars = maxvars)
           
         })
        variables1 = c(variables_cvs, variables1);
    }
   
#   vars_all=super$integrate(variables1)

   private$saveVars(variables1)
   invisible(variables1)
   
 },


#'  Plot the results from evaluateAll
#' @param evals a tibble or data frame from evaluate method 
#' @param shape_color a list of column names to encode as both shape and color
#' @param shape defaults to shape_color but can be specified separately
#' @param text a list of column names to include as text labels
#' @param dotsize a numerical column to encode dotsize, or a number for constant sizes
#' @param color defaults to shape_color but can be specified separately 
#' @param linetype column controlling the type of line
#' @param showranges whether to show 95\% CI as ribbon
#' @param txtsize size of text
#' @param transformy transformation on y axis, from ggplot2
#' @param legend whether to show a legend
#' @param sep_by list of column names to separate into diff plots
#' @param scales used in faceting can be free or fixed
#' @param point whether to show dots
#' @param line whether to include line
#' @param labelsize the size of the text labels
#' @param grid_y how to facet on y 
#' @param grid_x how to facet on x
#' @param title the title
#' @param title1 columns to include in the title
#' @returns nested ggplot2 plots
plotEval=function(eval3=self$evals,
                   shape_color=c("pheno","subpheno"),
                   shape=shape_color,
                   text =c(), ## variable
                   dotsize="nsamps",
                   color=shape_color,
                   linetype=shape,
                   showranges=TRUE, ## linetype="fullmodel"
                   txtsize=1,transformy = "identity",legend=FALSE,sep_by="",scales="free",point=TRUE,line=TRUE,
                   labelsize=2,
                   grid_x = c("cohort","measure"),grid_y = "cv_full",title="", title1=""
){
  grid1 = grid_y; grid0 = grid_x;
  grid0= grid0[grid0 %in% names(eval3)];  grid1= grid1[grid1 %in% names(eval3)]
  beam_levs = c(sort(unique(as.numeric(as.character(eval3$beam)))),"combined"); 
  
  eval3$beam = factor(eval3$beam, levels =beam_levs  )
  
  #eval3$sign = as.character(eval3$sign)
  #eval3$sign = factor(eval3$sign, levels = c(-1,0,1), labels=c("-","","+"))
  #  eval3 = eval3[,names(eval3) %in% c(shape_color,sep_by, linetype, text, dotsize, color, title1, grid0, grid1, "numvars","mid","low","high"),drop=FALSE]
  
  
  l1 = apply(eval3,1,paste, collapse="::");
  #print(which(duplicated(l1)))
  eval3 = eval3[!duplicated(l1),,drop=FALSE]
  showtext = length(text)>0
  linetype_nme = paste(linetype, collapse="_");
  grid0_nme = paste(grid0, collapse="_")
  grid1_nme = if(length(grid1)==0) "" else  paste(grid1,collapse="_");
  shape_nme = paste(shape,collapse="_")
  color_nme = paste(color,collapse="_")
  text_nme = paste(text,collapse="_")
  
  # sep_by_nme = "sep_by"
  subphens = table(eval3$subpheno)
  subphens = subphens[order(as.numeric(unlist(lapply(names(subphens), function(str)strsplit(str,"\\|")[[1]][1]))))]
  eval3$subpheno = factor(eval3$subpheno, levels = names(subphens))
  if(!is.null(eval3$measure))eval3$measure = factor(eval3$measure)
  eval3 = .modify(eval3, color, color_nme)
  eval3 = .modify(eval3, shape, shape_nme)
  eval3 = .modify(eval3, text, text_nme)
  
  eval3 = .modify(eval3, linetype, linetype_nme)
  #  eval3 = .modify(eval3, sep_by, sep_by_nme)
  eval3 = .modify(eval3,grid0, grid0_nme,collapse="\n")
  if(length(grid1)>0){
    eval3 = .modify(eval3,grid1, grid1_nme, collapse="\n")
  }
  eval2 = eval3
  
  phenos = unique(eval2$sep_by)
  names(phenos)=phenos
  
  
  
  
  eval2$numvars = as.numeric(eval2$numvars)
  eval2 = subset(eval2, !is.na(mid))
  # eval2$isfull = (eval2$isfull+1)/2.0
  
  
  ggps =  plotInternal(eval2, sep_by, title1, point, dotsize, shape_nme, color_nme,
                       linetype_nme, text_nme, labelsize, showranges,
                       grid0_nme, grid1_nme, scales, txtsize, line, transformy, showtext, legend);
  self$eval_plot = ggps;
  invisible(ggps)
}





  


))
