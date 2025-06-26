


.calcWall_2<-function(data,  var){
  prev_i = list(W_all =matrix(nrow=0, ncol=0) , var = c() )
  for(k in 1:length(var)){
    W_all2 = .calcWall_1(data, var[[k]], prev_i)
    prev_i = list(W_all = W_all2, var = var[1:k])
  }
  prev_i$W_all
}
#colk is family col


.sumMatrices<-function(matrices){
  if(length(matrices)==0) return(NULL)
  m=matrices[[1]]
  varn = dimnames(m)[[2]]
  if(length(matrices)>1){
    for(k in 2:length(matrices)){
      mi1 = match(varn,dimnames(matrices[[k]])[[2]])
      m[!is.na(mi1)] =m[!is.na(mi1)]+ matrices[[k]][mi1[!is.na(mi1)]]
      m[is.na(mi1)] = 9999
      #m = m+matrices[[k]]
    }
  }
  m
}


mStateObj<-R6Class("mStateObj",
                   public = list(
                     var="list",
                     logpvs="vector",
                     logpv="double",
                     cumpv="double",
                     nme="character",
                     initialize = function(var1, logpv1, prev_i = NULL   ){
                       if(!is.null(prev_i)){
                         self$var = c(prev_i$var, var1)
                         self$logpvs = c(prev_i$logpvs, logpv1)
                       }else{
                         self$var = var1
                         self$logpvs = logpv1
                       }
                       self$logpv = logpv1
                       self$cumpv = .sumChisq(self$logpvs)
                       self$nme = paste(sort(names(self$var)))
                     }
                  
                   )
)


stateObj1<-R6Class("stateObj1", ## cut down stateObj
                   public = list(
                     var="list",
                     name="character",
                     name_prev="character",
#                     varnames="list",
                     initialize = function(
                                           prev_i = NULL,
                                           b_i = NULL,
                                           var = list()
                     ){
                       #  names(data_nmes) = data_nmes
                       if(is.null(prev_i)){
                         self$var = var
                            self$name=paste(unlist(lapply(self$var, paste, collapse=".")),collapse="-")

                         self$name_prev=""
                       }else{
                         self$name_prev = prev_i$name
                         var = prev_i$var
                         var_st = unlist(lapply(var, paste, collapse="."))
                         if(paste(b_i,collapse=".") %in% var_st) {
                           print(paste("selected same variable twice"))#,names(todoInds[[b_i[1]]])[b_i[2]]))
                           warning(paste(" duplicated variable")) #,names(todoInds[[b_i[1]]])[b_i[2]]))
                           return(NULL)
                         }
                         vars = c(var,list(b_i))
                         self$name=paste(unlist(lapply(vars, paste, collapse=".")),collapse="-")
                         vars0=vars
#                         varnames=lapply(vars0, function(vv) c(names(data$data)[vv[1]],dimnames(data$data[[vv[1]]])[[2]][vv[2]]))
                         self$var = vars
#                         self$varnames=varnames
                       }
                     },
                      to_json=function(){
                       names(self$var) =1:length(self$var) 
                      json= toJSON(list(var = self$var,name = self$name, name_prev=self$name_prev))
                      json
                      }
                   )
                   )

stateObj<-R6Class("stateObj",##represents a state of the model
                  public = list(
                    var="list", var_names="list",
                    varnames="list",
                    name="character",
                    name_prev="character",
                    betas_proj="vector", #list of vector
                    W_all="matrix",  ##list of matrix,
                    constants_proj="list",
                    betas="list",
                    nonNA="logical",
                    mean_x="numeric",
                    initialize = function(phensi,data,
                                          betas_new, 
                                          constants_proj,
                                         k, #mean_y,
                                          prev_i = NULL,
                                          b_i = NULL,b_i_name=NULL,
                                          var = list(), varnames = list(),
                                          W_all = matrix(nrow=0, ncol=0),
                                         mean_x = NULL,
                                         proj=F
                                         ){
                      self$constants_proj=constants_proj
                      self$nonNA = data$looc$incl[,k];
                      train = data$train
                      
                      #if(length(phensi)>1) stop("!!")
                      family = unlist(lapply(names(phensi), function(x) getOption("fspls.family",strsplit(x,"\\.")[[1]][1])))
                      names(family)=names(phensi)
                      if(is.null(prev_i)){
                        self$var = var
                        self$var_names = var
                        self$varnames=var
                        self$mean_x=c()
                        self$name=""
                   #     self$cum_pvs_proj = list()
                        self$betas_proj= if(length(var)==0) c() else lapply(data$y, function(yy)(matrix(0,nrow=1, ncol= length(var)))) #lapply(datas, function(x) return(c()))
                        self$betas = if(length(var)==0) NULL else lapply(data$y, function(yy)(matrix(0,ncol=1, nrow= length(var)))) 
                        self$name=if(length(var)==0) "" else paste(unlist(lapply(var, paste, collapse=".")),collapse="-")
                        self$W_all = W_all
                      }else{
                        self$name_prev=prev_i$name
                        var = prev_i$var
                        self$mean_x = c(prev_i$mean_x, mean_x)
                        var_st = unlist(lapply(var, paste, collapse="."))
                       
                        if(paste(b_i,collapse=".") %in% var_st) {
                          print(paste("selected same variable twice"))#,names(todoInds[[b_i[1]]])[b_i[2]]))
                          warning(paste(" duplicated variable")) #,names(todoInds[[b_i[1]]])[b_i[2]]))
                          return(NULL)
                        }
                        vars = c(var,list(b_i))
                        
                        self$name=paste(unlist(lapply(vars, paste, collapse=".")),collapse="-")
                        vars0=vars
                        varnames=paste(b_i_name, collapse=".")
                        W_all_new = NULL
                        if(proj){
                          W_all = prev_i$W_all
                          W_all_new =  data$calcWall(b_i, prev_i$var, prev_i$W_all)
                        }
                        self$var = vars
                        self$var_names = c(prev_i$var_names, list(b_i_name))
                        self$varnames=c(prev_i$varnames, varnames)
                        beta_nme = names(betas_new)
                        names(self$var) = names(self$var_names)
                        nme_betas_new = names(betas_new); names(nme_betas_new) = nme_betas_new
                        self$betas_proj=lapply(nme_betas_new, function(b_n){
                          b_n1 = betas_new[[b_n]]
                          if(family[[b_n]]=="multinomial"){
                            b_n1 =b_n1[[1]]
                          }else{
                            b_n1 = as.matrix(data.frame(b_n1), nrow = 1)
                          }
                          if(proj){
                          rbind(prev_i$betas_proj[[b_n]],b_n1)
                          }else{
                            b_n1
                          }
                        })
                        if(proj){
                        self$betas=lapply(nme_betas_new, function(b_n) {
                          bp = self$betas_proj[[b_n]]
                          W_all_new %*% bp
                        })
                        self$W_all = W_all_new
                        }else{
                          self$betas = self$betas_proj
                          self$W_all = NULL
                        }
                        #if(family=="multinomial"){
                        #  self$betas = list(self$betas)
                        #  names(self$betas) = beta_nme
                        #}else{
                        #  self$betas = lapply(data.frame(self$betas), as.matrix, nrow = nrow(self$betas), ncol=1)
                        #}
                       
                      }
                    },
                  simplify=function(){
                    names(self$var_names)=self$varnames
                    list(betas=self$betas, constants_proj = self$constants_proj,
                         mean_x = self$mean_x,
                         var_names = self$var_names) #, pvs = self$pvs_proj)
                  },
                  updateConst=function(phensi,ypred, data,  k, transform_func, useglm=F, verbose=F){
                    #if(length(phensi)>1) stop("!!")
                    family = unlist(lapply(names(phensi), function(x) getOption("fspls.family",strsplit(x,"\\.")[[1]][1])))
                    names(family)=names(phensi)
                    non_na_x = data$getNonNA(self$var) 
                    self$constants_proj = vector("list", length(phensi))
                    names(self$constants_proj) = names(phensi)
                    na_k=non_na_x & self$nonNA
                  for(kk1 in 1:length(phensi)){
                    #kk1 = 1;  
                    kk = names(phensi)[[kk1]]
                    phensi1 = phensi[[kk1]]
                    y=  if(family[[kk]]=="multinomial")   attr(data$y[[kk]],"factor")[na_k] else data$y[[kk]][na_k,,drop=F]
                    y = transform_func(y)
                    yp1 =ypred$ypreds[[kk1]][na_k,,drop=F]
                    if(family[[kk]]=="multinomial"){
                        levs = levels(y)# c(0,1:length(self$betas[[k]]))
                        m1 = try(multinom(y~as.matrix(yp1),trace=F))
                        sm  = summary(m1, digits=3)
                        const_term = sm$coefficients[,1]
                        self$constants_proj[[kk1]] = const_term
                    }else if(family[[kk]]=="ordinal"){
                        for(kk_1 in 1:length(phensi1)){
                          kk_ = phensi1[kk_1]
                          y1c = y[,kk_]
                          levs1 = min(y1c, na.rm=T):max(y1c,na.rm=T)
                          consts = rep(0, length(levs1)-1)
                          names(consts) = levs1[-length(levs1)]
                          df = data.frame(list(y=factor(y1c,levels=sort(unique(y1c))),x= yp1[,1]))
                          m1=try(polr(y~x,  data=df,Hess=T, method="logistic"))  
                          if(abs(m1$coefficients[[1]]-1)>0.5){
                            print(m1$coefficients)
                            print(self$var)
                            stop(" something gone wrong")
                          }
                          if(inherits(m1,"try-error")){
                            stop("polr error in updating constant in stateObj")
                            ##prob not a great way to do this but avoids errors
                           # gl = glm(data$y[non_na_x & nonNA,kk]~ yp1[,1], family="gaussian")
                            #self$constants_proj[[kk]][[kk_1]] = consts_prev[[kk]] #rep(NA, length(levs1)-1) #gl$coefficients[1]
                          }else{
                            self$constants_proj[[kk]][[kk_1]]= m1$zeta
                          }
                        }
                     }else {
                       spike_slab_iter = getOption("spike_slab_iter",0)
                       ones = rep(1, length(yp1[,1]))
                       for(kk_1 in 1:length(phensi1)){
                         kk_ = phensi1[kk_1]
                         y1c = y[,kk_]
                        if(spike_slab_iter>1){
                         if(family[[kk]]=="binomial"){
                           ab=logit.spike(y1c~yp1[,kk_1], niter= spike_slab_iter, ping =-1)
                         }else{
                           if(getOption("logprint",F)) print("gaussian spike slab")
                           ab=lm.spike(y1c~yp1[,kk_1], niter= spike_slab_iter, ping =-1)#,weights=w)
                         }
                         sm = summary(ab)
                         x_ind = match(c("x","(Intercept)"),dimnames(sm$coefficients)[[1]])
                         coeff = sm$coefficients[x_ind[1],c(1,2,3,5)]
                         self$constants_proj[[kk1]][kk_1]= sm$coefficients[x_ind[2],1]
                         #const_term 
                       }else{
                         if(useglm){
                           self$constants_proj[[kk1]][kk_1] = tryCatch({
                             nonNAy = !is.na(y1c)
                           ridge=glmnet(cbind(ones[nonNAy],yp1[nonNAy,kk_1]),y1c[nonNAy],family=family[[kk]], alpha = 0)
                           rbeta <- coef(ridge,s=min(ridge$lambda))
                           if(abs(rbeta[3,1]-1)>0.1) stop("problem with weights")
                        if(verbose)   print(rbeta)
                          rbeta[1,1]
                           }, error=function(w) {
                             gl = glm(y1c[nonNAy]~ yp1[nonNAy,kk_1], family=family[[kk]])
                             if(abs(gl$coefficients[2]-1)>0.1) stop("problem with weights")
                             if(verbose)   print(gl$coefficients)
                             gl$coefficients[1]
                           
                           })
                         }else{
                        self$constants_proj[[kk1]][kk_1]  <- tryCatch({
                          gl = glm(y1c~ yp1[,kk_1], family=family[[kk]])
                         # print(gl$coefficients)
                          gl$coefficients[1]
                        }, warning=function(w) {
                          print("using glmnet to regularise ")
                          #ones = rep(1, length(yp1[,1]))
                          nonNAy = !is.na(y1c)
                          ridge=glmnet(cbind(ones[nonNAy],yp1[nonNAy,kk_1]),y1c[nonNAy],family=family[[kk]], alpha = 0)
                          rbeta <- coef(ridge,s=min(ridge$lambda))
                          if(verbose)   print(rbeta)
                          rbeta[1,1]
                          
                        })
                         }
                       }
                       }
                     }
                  
                   
          names(  self$constants_proj) = names(phensi)
                  }
                  }   
                  )
)
