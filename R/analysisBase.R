

.print<-function(file, x, width = 30,height  = 20){
  if(is.null(x)) return;
  pdf(file ,width=60, height=20);
  print(x);
  dev.off()
}




check_flags<-function(flags){
  angles_only = .readFlag(flags, 'angles_only', TRUE)
  topn = .readFlag(flags, "topn", 20)
  beam = .readFlag(flags, "beam",1)
  if(topn<beam&& !angles_only) stop("beam should be less than topn")
 
  if(!is.null(.readFlag(flags,"nrep",NULL))){
    flags[['nfold']] = flags[['nrep']]
    warning("replaced nrep with nfold")
  }
  if(!is.null(.readFlag(flags,"batch",NULL))){
    flags[['batchsize']] = flags[['batch']]
    
    warning("replace batch with batchsize")
  }
  
  invisible(flags)
}



#' Analysis Base Class
#'
#' @description
#' A class that encapsulates a dataset holder and/or analysis of datasets
#'
#' @details
#' Available methods:
#' \itemize{
#'   \item \code{new()} - Create a new instance
#'   \item \code{name()} - Get dataset name
#'   \item \code{clone()} - Clone the object
#' }
#'
#' @export
analysisBase<-R6::R6Class("analysisBase", 
                   private = list(
                     sigsdir="character",
                     sigs="environment",
                     nme="character",
                     dbDir="character",
                     exptid="character",
                     flags ="list",
                     data_types = "list",
                     transform_x="list",
                     phens="list",
                     useDB="boolean",
                    
                     
                     savePvals=function(ri , varnames, k1, data_nme=private$nme, useCurrVarnames=TRUE){
                       
                       if(!is.null(private$sigs)) {
                         #                         private$sigs$savePvals(flags, phens, data_nme, ri, varnames,k1, useCurrVarnames =FALSE)
                         #super$savePvals(flags,phens1, ri, varnames,k1,useCurrVarnames=TRUE)
                         
                         private$sigs$savePvals(private$flags,private$phens, data_nme, ri, varnames,k1,useCurrVarnames=useCurrVarnames)
                       }
                     },
                     saveVars = function(variables){
                       if(!is.null(private$sigs)){
                         private$sigs$saveVars(variables,replace=TRUE)   #saving local
                       }
                       self$variables = variables;
                     },
                     savePredictions=function(predictions){
                       self$predictions = predictions
                     },
                     
                     loadVars = function(){
                       phens = private$phens;
                       if(is.null(private$sigs)) return (NULL)
                       private$sigs$loadVars(private$flags, phens)
                     },
                    loadModels=function(){
                       sigDB = private$sigs 
                       if(!is.null(sigDB) ){
                         all_models =try( sigDB$loadModels(private$flags,private$phens))
                         if(inherits(all_models,"try-error")) {
                           warning(paste("problem reading from DB .. recalculating"))
                         }else if(!is.null(all_models) && length(all_models$models)>0){
                           return(all_models)
                         }
                       }
                     },
                     saveModels = function(all_models_){
                       sigDB = private$sigs 
                       if(!is.null(sigDB)){
                         sigDB$saveModels(all_models_)
                       }
                       self$all_models = all_models_;
                       
                     },
                     loadEval=function(){
                       sigDB = private$sigs 
                       
                       if(!is.null(sigDB) ){
                         
                         eval1 = sigDB$loadEval(private$flags,private$phens)
                         if(!is.null(eval1) && nrow(eval1)>0){
                           return(eval1)
                         }
                       }
                       
                     },
                     saveEval = function(eval2){
                       sigDB = private$sigs 
                       if(!is.null(sigDB) ){
                        sigDB$saveEval(eval2, private$flags,private$phens)
                         
                       }
                      self$evals=eval2;
                     },
                     expt_id = function(){
                       return(private$exptid);
                     },
                    
                     
                     updateExpt=function(phens, flags, transform_x, data_types){
                       private$flags = flags;
                       private$transform_x = transform_x;
                       private$flags$transform_x = toJSON(transform_x);
                       private$flags$data_types = toJSON(data_types);
                       private$phens = phens;
                       private$data_types = data_types;
                       #if(is.null(private$flags[['data_types']]) || private$flags[['data_types']]=="{}")private$flags[['data_types']]=toJSON(data_types);
                    
                       if(!is.null(private$sigs)){
                         private$exptid =  private$sigs$getExpt(flags, phens, add_new=TRUE)
                       }
                       private$flags
                       check_flags(private$flags)
                       options(private$flags);
                     },
                    saveDataPlots=function(data_plots){
                      self$data_plots = data_plots
                    },
                     savePlotResults = function(plot_results_eval){
                       self$plot_results_eval = plot_results_eval
                     }
                   
                     
                    
                   ),
                   public = list(
                                       
                     variables=NULL, ##selected variables
                     all_models=NULL,
                     predictions=NULL,
                     evals=NULL,
                     data_plots =NULL,
                     pred_plots=NULL,
                     evals_pretty =NULL,   
                     eval_plot=NULL,
                     betas_table = NULL,
                     vars_table = NULL,
                     #' @description Create a new instance of base class
                     #' @param nme name
                     #' @param dims dims
                     #' @param flags a list object specifying options
                     #' @param dbDir location for database to be written, default to tempdir
                     initialize=function(
                              nme,
                              dims = NULL,
                              flags = list(),
                              dbDir=tempdir()
                            
                     ){
                       
                       private$phens = NULL;
                       private$transform_x = NULL
                       
                       useDB = !is.null(dbDir);
                       private$useDB =useDB;
                       if(is.null(nme)) stop("nme should not be null")
                      private$dbDir = dbDir;
                      private$exptid=0;
                      private$nme=sub("/",".",nme)
                      private$flags = flags
                      private$sigsdir=paste(dbDir,paste0("fspls_signatures__",private$nme),sep="/")
                   #   private$dims = dimsl
                      # private$sigsdir=paste(dbDir,"fspls_signatures",sep="/")
                      # private$sigs=    sigEnv$new(private$sigsdir,nme1,flags, NULL, clear=FALSE)
                       private$sigs=   if(useDB) sigEnv$new(private$sigsdir,private$nme,private$flags, dims, clear=FALSE) else NULL;
                       
                     },
                   #' @description get dataset name;
                   #' @returns assigned name of dataset
                   name = function(){
                     return (private$nme)
                   },
                   
                  
                
                   #' @description integrates variables from different CV runs
                   #' @param max_model_length length maximum model length
                   #' @param variables a variables object
                   #' @param order a new ordering
                   integrate=function(variables, max_model_length=1000, max_beam = 1000, order=NULL){
                    # if(names(variables[['full']][[1]]=="nvar_1")) {
                       ##if variables is all variables, we take the last of each
                      
                     #}
                     flags = private$flags;
                  #   keepAll = .readFlag(flags,"keep_all", TRUE);
                   #  if(keepAll){
                         variables2 = lapply(variables, function(x) {
                           x1=x[[min(length(x), max_model_length)]]
                          
                           lapply(x1, function(x2){
                             if(length(order)==length(x2$var_names)) x2$var_names = x2$var_names[order]
                             x2
                           })
                           
                         
                         })
                      variables = variables2
                    # }
                     phens = private$phens;
                    useDB=private$useDB;
                     useAngles = .readFlag(flags, 'angles_only', TRUE)
                    # useAngles = !is.null(flags$angles_only) && flags$angles_only
                     full_index = length(variables)
                     beam_m = min(length(variables[[full_index]]), max_beam)
                     beams = 1:beam_m
                     names(beams)=beams
                     names(variables) = 1:length(variables)
                     
                     vars_combined=lapply(beams, function(beam){
                       vars_all = list()
                       vars_all1 = list()
                       vars_all2 = list()
                       for(repn in names(variables)){
                         full = repn==full_index
                         if(beam>length(variables[[repn]])) {
                           message("skipping")
                           next;
                         }
                         var1 = variables[[repn]][[beam]]   ### only taking the top1
                         var2 = var1$var_names
                         if(length(var2)>0){
                           names(var2) = names(var1$var_names)
                           cumpv = if(useAngles) var1$cum_angle else  var1$cum_pv#lapply(var1, function(vv) attr(vv,"cumpv"))
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
                       result = list(variables = vars_all1, inds = vars_all,cumpv=vars_all2, beam=beam,flags = flags,phens = phens )# ,transf= vars_all3) 
                       attr(result,"variables")=variables;
                       result
                     })
                     
                     
                     vars_combined
                   },

                  #' @description extracts the variables into a table
                  tabulate_variables=function(variables =self$variables){
                    df2 = .merge1_new(lapply(variables, function(v){
                      .merge1_new(v)
                    }), addName="nrep")
                    res1 = (tibble::as_tibble(df2))  
                    self$vars_table = res1;
                    invisible(res1)
                  },
                  
                  #' @description  tabulates the evals
                  #' @param eval1 a tibble or data frame from evaluateAll method from dataH
                  #' @param cv_index can be full, avg or merged, or a numerical index
                  #' @param measure the measure to be included, eg AUC
                  #' @param values_from which values to include, can be low, mid, high or mid_perc
                  #' @returns formatted table
                   tabulateEvals=function(eval1 = self$evals, minvars =1, cvs = c("avg","merged","full"), beams = 1:10
                  ){
                     if(is.null(eval1)) {
                       warning("need to use evluate function first")
                       return(NULL)
                     }
                    id_cols = c("subpheno","pheno","cv_index","beam","measure","family","data")
                    id_cols = id_cols[id_cols %in% names(eval1)]
                    eval1 = subset(eval1, cv_index %in% cvs)
                    eval1 = subset(eval1, numvars>=minvars)
                    eval1 = subset(eval1, beam %in% beams)
                    eval3 = eval1[,match(c("mid","low","high"),names(eval1))];
                    range=apply(eval3, 1,function(xx) {
                      
                      if(is.na(xx[1])) return("NA")
                      
                      #if(xx[3]>100*xx[1] || xx[2}<100*xx[1]) return("NA(-Inf, Inf)")
                      if(is.na(xx[2])) return(gsub("\\s+","",sprintf("%4.2g",xx[1])));
                      gsub("\\s+","",sprintf("%4.2g (%4.2g-%4.2g)", xx[1], xx[2], xx[3]))
                      #    paste(xx1[1],paste0("(",paste(xx1[-1], collapse="-"),")"))
                    })
                    eval1 = tibble::add_column(eval1, range=range)
                    res = pivot_wider(eval1,id_cols=id_cols,values_from=range, names_from="numvars")
                    res2 = pivot_wider(eval1,id_cols=id_cols,values_from=variable, names_from="numvars")
                    model = apply(res2[,-(1:length(id_cols))],1, function(xx) paste(xx[!is.na(xx)],collapse=";"))
                    res3 = res |> tibble::add_column(model,.after = "family")
                    #cvs = unique(eval1$cv_index); names(cvs) = cvs; cvs = cvs[cvs%in% c("avg","merged","full")]
                    families = unique(res$family); names(families) = families;
                    res_out=lapply(families, function(fam){
                      res4= subset(res3, family==fam)
                      res5 = res4 |> dplyr::arrange(cv_index, pheno, subpheno,measure)
                      res5
                    })
                    self$evals_pretty = res_out
                    invisible(res_out)
                  }
                  
                  
                 
                 
                     
                   )
                  
                 
);
                     