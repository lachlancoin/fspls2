

check_flags<-function(flags){
  if(flags$topn<flags$beam) stop("beam should be less than topn")
  
  if(!is.null(.readFlag(flags,"nrep",NULL))){
    flags[['nfold']] = flags[['nrep']]
    warning("replaced nrep with nfold")
  }
  if(!is.null(.readFlag(flags,"batch",NULL))){
    flags[['batchsize']] = flags[['batch']]
    
    warning("replace batch with batchsize")
  }
  if(flags$topn<flags$beam) stop("beam should be less than topn")
  
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
                     }
                     
                   
                     
                    
                   ),
                   public = list(
                                          
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
                   #' @param variables a variables object
                   integrate=function(variables){
                     flags = private$flags;
                     phens = private$phens;
                    useDB=private$useDB;
                     useAngles = !is.null(flags$angles_only) && flags$angles_only
                     full_index = length(variables) 
                     beams = 1:length(variables[[full_index]])
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
                   }
                 
                 
                     
                   )
                  
                 
);
                     