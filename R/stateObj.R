


.calcWall_2<-function(data,  var){
  prev_i = list(W_all =matrix(nrow=0, ncol=0) , var = c() )
  for(k in 1:length(var)){
    W_all2 = .calcWall_1(data, var[[k]], prev_i)
    prev_i = list(W_all = W_all2, var = var[1:k])
  }
  prev_i$W_all
}
#colk is family col


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
                   tbls="list",
                   pvs="list",
                   sumPv="numeric",
                    initialize = function(phensi,data,
                                          betas_proj, 
                                          constants_proj,
                                          tbls,
                                         k, #mean_y,
                                          prev_i = NULL,
                                          b_i = NULL,b_i_name=NULL,
                                          var = list(), varnames = list(),
                                          W_all = matrix(nrow=0, ncol=0),
                                         mean_x = NULL, pvs = c(),
                                         useoffset=T
                                         ){
                      self$pvs = pvs
                      self$sumPv = .sumChisq(unlist(pvs))
                      self$constants_proj=constants_proj
                      self$nonNA = data$looc$incl[,k];
                      train = data$train
                      self$tbls = tbls
                      self$W_all = W_all
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
                        self$var = vars
                        self$var_names = c(prev_i$var_names, list(b_i_name))
                        self$varnames=c(prev_i$varnames, varnames)
                        self$betas_proj=betas_proj
                        if(FALSE){
                        beta_nme = names(betas_new)
                        names(self$var) = names(self$var_names)
                        nme_betas_new = names(betas_new); names(nme_betas_new) = nme_betas_new
                       
                        self$betas_proj=lapply(nme_betas_new, function(b_n){
                          b_n1 = betas_new[[b_n]]
                          if(family[[b_n]]=="multinomial"){
                            b_n1 = b_n1[[1]]
                            if(!is.matrix(b_n1)){
                              b_n1 =t(as.matrix(b_n1))
                            }
                            #print(b_n1)
                          }else{
                            b_n1 = as.matrix(data.frame(b_n1)) #, nrow = 1)
                          }
                          if(!useoffset) return(b_n1)
                          if(nrow(b_n1)==1){ 
                            return( rbind(prev_i$betas_proj[[b_n]], b_n1))
                          }
                          rescale =   b_n1[1,]
                          b_n1 = b_n1[-1,,drop=F]
                          bet_n = prev_i$betas_proj[[b_n]]
                          for(jj in 1:nrow(bet_n)){
                            bet_n[jj,] = bet_n[jj,]*rescale
                          }
                          rbind(bet_n, b_n1)
                         
                        })
                        }
                        #self$betas=self$betas_proj  ## now these are the same
                        
#                          self$betas=lapply(self$betas_proj, function(bp){
 #                           b2 = self$W_all %*% bp
  #                        } )
                      }
                    },
                  setOffset=function(){
                    nme_betas_new = names(self$betas_proj)
                    #b_n= nme_betas_new[[1]]
                    for(b_n in nme_betas_new){
                      offset =self$mean_x%*%  self$betas[[b_n]]
                      vv = unlist(self$constants_proj[[b_n]])-offset[1,]
                      self$constants_proj[[b_n]] = vv
                    }
                  },
                   # self$betas=lapply(nme_betas_new, function(b_n) {
                   #    bp = self$betas_proj[[b_n]]
                   #    W_all_new %*% bp
                   #  })
##simplifies and translates into original space
                  simplify=function(transf){
                    nmebp= names(self$betas_proj); names(nmebp)=nmebp
                    self$betas=lapply(nmebp, function(bp){
                      self$W_all[[bp]] %*% self$betas_proj[[bp]]
                    })
                    self$setOffset() 
                    names(self$var_names)=self$varnames
                    list(betas=self$betas, constants_proj = self$constants_proj,
                         mean_x = self$mean_x,tbls = self$tbls,
                         transf=transf,
                     #    W_all = self$W_all,
                         var_names = self$var_names) #, pvs = self$pvs_proj)
                  },
                  updateConst=function(phensi,ypred, data,  transform_func, useglm=F, verbose=F,update=F
                                       ){
                    #if(length(phensi)>1) stop("!!")
                    family = unlist(lapply(names(phensi), function(x) getOption("fspls.family",strsplit(x,"\\.")[[1]][1])))
                    names(family)=names(phensi)
                    non_na_x = data$getNonNA(self$var) 
                    constants_proj = vector("list", length(phensi))
                    names(constants_proj) = names(phensi)
                    na_k=non_na_x & self$nonNA
                  for(kk1 in 1:length(phensi)){
                    if(verbose)print(names(phensi)[[kk1]])
                    #kk1 = 1;  
                    kk = names(phensi)[[kk1]]
                    phensi1 = phensi[[kk1]]
                    y=  if(family[[kk]]=="multinomial")   attr(data$y[[kk]],"factor")[na_k] else data$y[[kk]][na_k,,drop=F]
                    y = transform_func(y)
                    yp1 =ypred$ypreds[[kk1]][na_k,,drop=F]
                    if(family[[kk]]=="multinomial"){
                        levs = levels(y)# c(0,1:length(self$betas[[k]]))
                        m1 = try(multinom(y~as.matrix(yp1[,-1]),trace=F))
                        sm  = summary(m1, digits=3)
                        if(verbose){
                          print(sm$coefficients)
                        }
                        const_term = sm$coefficients[,1]
                        constants_proj[[kk1]] = const_term
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
                            warning(" something gone wrong")
                          }
                          if(inherits(m1,"try-error")){
                            stop("polr error in updating constant in stateObj")
                            ##prob not a great way to do this but avoids errors
                           # gl = glm(data$y[non_na_x & nonNA,kk]~ yp1[,1], family="gaussian")
                            #self$constants_proj[[kk]][[kk_1]] = consts_prev[[kk]] #rep(NA, length(levs1)-1) #gl$coefficients[1]
                          }else{
                            if(verbose){
                              print(m1$coefficients)
                              print(m1$zeta) 
                            }
                            constants_proj[[kk]][[kk_1]]= m1$zeta
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
                         constants_proj[[kk1]][kk_1]= sm$coefficients[x_ind[2],1]
                         #const_term 
                       }else{
                         if(useglm){
                           constants_proj[[kk1]][kk_1] = tryCatch({
                             nonNAy = !is.na(y1c)
                           ridge=glmnet(cbind(ones[nonNAy],yp1[nonNAy,kk_1]),y1c[nonNAy],family=family[[kk]], alpha = 0)
                           rbeta <- coef(ridge,s=min(ridge$lambda))
                         if(verbose)  print(rbeta)
                       #    if(abs(rbeta[3,1]-1)>diff_thresh) stop("problem with weights")
                        if(verbose)   print(rbeta)
                          sum(rbeta[1:2,1])
                           }, error=function(w) {
                             gl = glm(y1c[nonNAy]~ yp1[nonNAy,kk_1], family=family[[kk]])
                            # if(abs(gl$coefficients[2]-1)>diff_thresh) stop("problem with weights")
                             if(verbose)   print(gl$coefficients)
                             gl$coefficients[1]
                           
                           })
                         }else{
                        constants_proj[[kk1]][kk_1]  <- tryCatch({
                          gl = glm(y1c~ yp1[,kk_1], family=family[[kk]])
                         if(verbose) print(gl$coefficients)
                          gl$coefficients[1]
                        }, warning=function(w) {
                          print("using glmnet to regularise ")
                          #ones = rep(1, length(yp1[,1]))
                          nonNAy = !is.na(y1c)
                          ridge=glmnet(cbind(ones[nonNAy],yp1[nonNAy,kk_1]),y1c[nonNAy],family=family[[kk]], alpha = 0)
                          rbeta <- coef(ridge,s=min(ridge$lambda))
                          if(verbose)   print(rbeta)
                          sum(rbeta[1:2,1])
                          
                        })
                         }
                       }
                       }
                     }
                  
                   
          names(  constants_proj) = names(phensi)
                  }
                    if(update)self$constants_proj = constants_proj
                    constants_proj
                  }   
                  )
)
