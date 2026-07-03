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
                     flags ="list",
                     nme="character",
                     dbDir="character",
                     
                     savePvals=function(flags, phens1, ri , varnames, k1, data_nme=private$nme, useCurrVarnames=T){
                       if(!is.null(private$sigs)) {
#                         private$sigs$savePvals(flags, phens, data_nme, ri, varnames,k1, useCurrVarnames = F)
                         #super$savePvals(flags,phens1, ri, varnames,k1,useCurrVarnames=T)
                         
                        private$sigs$savePvals(flags,phens1, data_nme, ri, varnames,k1,useCurrVarnames=useCurrVarnames)
                       }
                     },
                      getExpt=function(flags, phens, add_new=T){
                            if(is.null(private$sigs)) return(0);
                        private$sigs$getExpt(flags, phens, add_new=T)
                      },
                     
                     post_process=function(variables, flags, phens, useDB=FALSE){
                       useAngles = !is.null(flags$angles_only) && flags$angles_only
                       full_index = length(variables) 
                       beams = 1:length(variables[[full_index]])
                       names(beams)=beams
                       names(variables) = 1:length(variables)
                       
                       vars_combined=lapply(beams, function(beam){
                         print(beam)
                         vars_all = list()
                         vars_all1 = list()
                         vars_all2 = list()
                         #vars_all3 = list() #funcstr
                         # func_inds = lapply(variables, function(vv) attr(vv,"func_ind"))
                         
                         for(repn in names(variables)){
                           full = repn==full_index
                           if(beam>length(variables[[repn]])) {
                             print("skipping")
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
                         list(variables = vars_all1, inds = vars_all,cumpv=vars_all2, beam=beam,flags = flags,phens = phens )# ,transf= vars_all3) 
                       })
                       if(useDB  && !is.null(private$sigs)){
                         private$sigs$saveVars(vars_combined,replace=T)   #saving local
                       }
                       
                       vars_combined
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
                       
                       mc.cores = .readFlag(flags, "mc.cores",2)
                       if(mc.cores>1){
                         dbDir= NULL;
                         
                       }
                       useDB = !is.null(dbDir);
                       if(is.null(nme)) stop("nme should not be null")
                      private$dbDir = dbDir;
                      private$nme=sub("/",".",nme)
                      private$flags = flags
                      private$sigsdir=paste(dbDir,paste0("fspls_signatures__",private$nme),sep="/")
                   #   private$dims = dimsl
                      # private$sigsdir=paste(dbDir,"fspls_signatures",sep="/")
                      # private$sigs=    sigEnv$new(private$sigsdir,nme1,flags, NULL, clear=F)
                       private$sigs=   if(useDB) sigEnv$new(private$sigsdir,private$nme,private$flags, dims, clear=F) else NULL;
                       
                     },
                   #' @description get dataset name;
                   #' @returns assigned name of dataset
                   
                   name = function(){
                     return (private$nme)
                   }
                     
                   )
                  
                 
);
                     