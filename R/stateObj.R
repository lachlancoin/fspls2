


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
                    initialize = function(phensi,data,
                                          betas_new, 
                                          constants_proj,
                                         k, #mean_y,
                                          prev_i = NULL,
                                          b_i = NULL,b_i_name=NULL,
                                          var = list(), varnames = list(),
                                          W_all = matrix(nrow=0, ncol=0)
                                         ){
                      self$constants_proj=constants_proj
                      self$nonNA = data$looc$incl[,k];
                      train = data$train
                      
                      if(length(phensi)>1) stop("!!")
                      family = unlist(lapply(names(phensi), function(x) strsplit(x,"\\.")[[1]][1]))
                      names(family)=names(phensi)
                      if(is.null(prev_i)){
                        self$var = var
                        self$var_names = var
                        self$varnames=var
                        self$name=""
                   #     self$cum_pvs_proj = list()
                        self$betas_proj= if(length(var)==0) c() else lapply(data$y, function(yy)(matrix(0,nrow=1, ncol= length(var)))) #lapply(datas, function(x) return(c()))
                        self$betas = if(length(var)==0) NULL else lapply(data$y, function(yy)(matrix(0,ncol=1, nrow= length(var)))) 
                        self$name=if(length(var)==0) "" else paste(unlist(lapply(var, paste, collapse=".")),collapse="-")
                        self$W_all = W_all
                      }else{
                        self$name_prev=prev_i$name
                        var = prev_i$var
                        var_st = unlist(lapply(var, paste, collapse="."))
                        W_all = prev_i$W_all
                        if(paste(b_i,collapse=".") %in% var_st) {
                          print(paste("selected same variable twice"))#,names(todoInds[[b_i[1]]])[b_i[2]]))
                          warning(paste(" duplicated variable")) #,names(todoInds[[b_i[1]]])[b_i[2]]))
                          return(NULL)
                        }
                        vars = c(var,list(b_i))
                        self$name=paste(unlist(lapply(vars, paste, collapse=".")),collapse="-")
                        vars0=vars
                        varnames=paste(b_i_name, collapse=".")
                        W_all_new =  data$calcWall(b_i, prev_i$var, prev_i$W_all)
                        self$var = vars
                        self$var_names = c(prev_i$var_names, list(b_i_name))
                        self$varnames=c(prev_i$varnames, varnames)
                        beta_nme = names(betas_new)
                        if(family=="multinomial"){
                          betas_new =betas_new[[1]]
                        }else{
                          betas_new = as.matrix(data.frame(betas_new), nrow = 1)
                        }
                        self$betas_proj=rbind(prev_i$betas_proj,betas_new)
                        self$betas=W_all_new %*% self$betas_proj 
                        if(family=="multinomial"){
                          self$betas = list(self$betas)
                          names(self$betas) = beta_nme
                        }else{
                          self$betas = lapply(data.frame(self$betas), as.matrix, nrow = nrow(self$betas), ncol=1)
                        }
                        self$W_all = W_all_new
                      }
                    },
                  simplify=function(){
                    names(self$var_names)=self$varnames
                    list(betas=self$betas, constants_proj = self$constants_proj, var_names = self$var_names) #, pvs = self$pvs_proj)
                  },
                  updateConst=function(phensi,ypred, data,  k){
                    if(length(phensi)>1) stop("!!")
                    family = unlist(lapply(names(phensi), function(x) strsplit(x,"\\.")[[1]][1]))
                    names(family)=names(phensi)
                    non_na_x = data$getNonNA(self$var) 
                    self$constants_proj = vector("list", length(phensi))
                    names(self$constants_proj) = names(phensi)
                    na_k=non_na_x & self$nonNA
                    kk1 = 1;  kk = names(phensi)[[kk1]]
                    phensi1 = phensi[[kk1]]
                    y=  if(family[[kk]]=="multinomial")   attr(data$y[[kk]],"factor")[na_k] else data$y[[kk]][na_k,,drop=F]
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
                          df = data.frame(list(y=factor(y1c,levels=sort(unique(y1c))),x= yp1))
                          m1=try(polr(y~x,  data=df,Hess=T, method="logistic"))  
                          if(inherits(m1,"try-error")){
                            warning("polr error in updating constant")
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
                         if(TRUE){
                           self$constants_proj[[kk1]][kk_1] = tryCatch({
                           ridge=glmnet(cbind(ones,yp1[,kk_1]),y1c,family=family[[kk]], alpha = 0)
                           rbeta <- coef(ridge,s=min(ridge$lambda))
                          rbeta[1,1]
                           }, error=function(w) {
                             gl = glm(y1c~ yp1[,kk_1], family=family[[kk]])
                             gl$coefficients[1]
                           })
                         }else{
                        self$constants_proj[[kk1]][kk_1]  <- tryCatch({
                          gl = glm(y1c~ yp1[,kk_1], family=family[[kk]])
                          gl$coefficients[1]
                        }, warning=function(w) {
                          print("using glmnet to regularise ")
                          ones = rep(1, length(yp1[,1]))
                          ridge=glmnet(cbind(ones,yp1[,kk_1]),y1c,family=family[[kk]], alpha = 0)
                          rbeta <- coef(ridge,s=min(ridge$lambda))
                          rbeta[1,1]
                        })
                         }
                       }
                       }
                     }
                  
                   
          names(  self$constants_proj) = names(phensi)
                  }
                    
                  )
)
