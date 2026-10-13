##collection of models

modelObj<-R6::R6Class("modelObj",
                       public = list(
                       
                         models="list",
                         nvar="list",
                         beams="list", ## probably don't need beam
                         k='integer',
                         initialize = function(k ){
                           self$k = k;
                           self$models  =list()
                           self$nvar =list();
                           self$beams = list();
                         },
                         set=function(prev_i1,beam){
                           varlen = length(prev_i1$var_names);
                           
                           nme = if(varlen==0) "empty" else  paste(unlist(lapply(prev_i1$var_names, paste, collapse=".")),collapse=";")
                          
                          # print(paste("setting1", nme, varlen))
                           self$models[[nme]]=prev_i1;
                           self$nvar[[nme]] = varlen
                           self$beams[[nme]]= beam
                         },
                         get = function(prev_i2, b_i_name){
                           pn = if(is.null(prev_i2)) c() else names(prev_i2$var_names);
                           nme1 = paste(c(pn, paste(b_i_name, collapse=".")), collapse=";")
                          if(nme1=="")nme1="empty"
                            self$models[[nme1]]
                         },
                         separate=function(){
                           nvar = unlist(self$nvar)
                           nvars = as.numeric(names(table(nvar)));
                           names(nvars) = nvars
                           lapply(nvars, function(nv){
                           #  print(nv)
                                beams = self$beams[nvar==nv]
                                mods = self$models[nvar==nv]
                                lapply(mods, function(x) x$simplify())
                           })
                          # 
                           
                         }
                         
                         
                       )
)


    

#models = modelObj$new(1)
