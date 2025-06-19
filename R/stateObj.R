


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
                   # var1="list",
                    pvs_proj="list",
                    cum_pvs_proj="list",
                    betasP="list",
#                    betas_offset="list", #list of vector
                    betas_proj="vector", #list of vector
                    W_all="matrix",  ##list of matrix,
                    const_off_proj="list",
 #                   const_off_offset="list",
                    constants_proj="list",
                  constants_proj_all="list",
                    betas="list",
                    const="list",
                    mean_x = "list",
                    tbls="list",
                    mean_y="list",
                    initialize = function(phensi,data,train, k,mean_y,
                                          prev_i = NULL,
                                          b_i = NULL,b_i_name=NULL,
                                          var = list(), varnames = list(),
                                          W_all = matrix(nrow=0, ncol=0)
                                         ){
                      #if(is.null(train)) train = data$train[[k]]
                      #  names(data_nmes) = data_nmes
                  #   mean_y = train$means_y[[k]]
                      family = unlist(lapply(names(phensi), function(x) strsplit(x,"\\.")[[1]][1]))
                      self$constants_proj=lapply(1:length(phensi), function(ik){
                        nme = names(phensi)[[ik]]
                        phensi1 = phensi[[ik]]
                        if(family[[ik]]=="multinomial") {
                          fact = attr(data$y[[ik]],"factor")
                          return(rep(0,length(levels(fact))-1 ))
                        }else if(family[[ik]]=="ordinal"){
                          levs = sort(unique(data$y[[ik]][,1]))
                          return(rep(0,length(levs )-1))
                        }else{
                          return(mean_y[[ik]][phensi1])
                        }
                      })
                      names(self$constants_proj) = names(phensi)
                      self$constants_proj_all = list()
                      if(is.null(prev_i)){
                        self$var = var
                        self$var_names = var
                        self$varnames=var
                     #  self$constants_proj = as.list(self$train[[1]]$mean_y)
                        self$name=""
                        self$cum_pvs_proj = list()
                        
                       # self$var1 = list()
                        self$betas_proj= if(length(var)==0) c() else lapply(data$y, function(yy)(matrix(0,nrow=1, ncol= length(var)))) #lapply(datas, function(x) return(c()))
                        self$betas = if(length(var)==0) NULL else lapply(data$y, function(yy)(matrix(0,ncol=1, nrow= length(var)))) 
                        #self$constants_proj = if(length(var)==0) NULL else lapply(data$y, function(yy) (rep(0, length(var))))
                        #self$betas_proj = self$betas_offset
                        self$name=if(length(var)==0) "" else paste(unlist(lapply(var, paste, collapse=".")),collapse="-")
                        self$W_all = W_all
                      }else{
                       # train = data$train
                        self$name_prev=prev_i$name
                        #family = train$family #unlist(lapply(train, function(xx) xx$family))
                        var = prev_i$var
                        var_st = unlist(lapply(var, paste, collapse="."))
                        #betas_offset=prev_i$betasoffset
                     
                        W_all = prev_i$W_all
                        if(paste(b_i,collapse=".") %in% var_st) {
                          print(paste("selected same variable twice"))#,names(todoInds[[b_i[1]]])[b_i[2]]))
                          #print(names(todoInds)[var])
                          warning(paste(" duplicated variable")) #,names(todoInds[[b_i[1]]])[b_i[2]]))
                          return(NULL)
                        }
                        vars = c(var,list(b_i))
                        self$name=paste(unlist(lapply(vars, paste, collapse=".")),collapse="-")
                       # print(paste("HERE", self$name))
                        vars0=vars
                        W_all_new = NULL
                        betas_new_proj=NULL
                        betas_new_offset = NULL
                        #ang=angle[b_i]
                        vv = b_i
                        varnames=paste(names(data$data)[vv[1]],dimnames(data$data[[vv[1]]])[[2]][vv[2]],sep=".")
                        betas_newP = NULL
                       # constants_proj=as.list(self$train[[1]]$mean_y)
                        pvs_proj = NULL
                          W_all_new =  data$calcWall(b_i, prev_i$var, prev_i$W_all)
                          
                          b_new_proj = data$calcBetaProj1(phensi,k,b_i,prev_i$var, convert=F)
                          betas_new_proj = lapply(1:length(b_new_proj$betas), function(kk){
                            dff=data.frame(cbind(prev_i$betas_proj[[kk]], b_new_proj$betas[[kk]]) )
                            names(dff) = c(names(prev_i$betas_proj)[[kk]], varnames)
                            dff
                          })
                          names(betas_new_proj) = names(b_new_proj$betas)
                          tbls = b_new_proj$tbls
                          constants_proj = b_new_proj$constants
                          pvs_proj = b_new_proj$pvs
                          
                          #betas_newP = lapply(1:length(train), function(i) {
                            bet_i = betas_new_proj
                            #    if(DRS>0) b_i = .DRS(b_i, DRS)
                            b_n1=lapply(1:length(bet_i), function(i1){
                              b_i1 = bet_i[[i1]]
                                bn=as.matrix((apply(b_i1,1,function(ikk) (W_all_new%*%ikk)[,1] )))
                                if(nrow(bn)!=ncol(b_i1))bn = t(bn)
                              bn
                            })
                            names(b_n1) = names(betas_new_proj)
                            betas_newP =  b_n1
#                          names(betas_newP) = names(train)
                       
                        mean_x = unlist(lapply(vars, function(vv) data$mean_x[[vv[1]]][vv[2]]))
                          #lapply(train, function(t) unlist(lapply(vars, function(vv) t$mean_x[[vv[1]]][vv[2]])))
                        #names(mean_x) = names(train)
                        #mean_y = train$mean_y #lapply(train, function(t) t$mean_y)
                        #names(mean_y) =names(train)
                        ##if family is binomial we dont offset by mean y
                        const_off_proj=NULL
                        #const_off_offset=NULL
#                          const_off_proj = lapply(1:length(train), function(i1){
                            beta_newPi = betas_newP
                            bnp_indices = 1:length(beta_newPi)
                            names(bnp_indices) = names(beta_newPi)
                            const_off_proj=NULL
                          # const_off_proj= lapply(bnp_indices, function(bnp_i){
                          #    m =  mean_x %*%beta_newPi[[bnp_i]]
                          #    if(!(family %in% c("binomial","ordinal")))  mean_y[bnp_i]-1*m[1,1] else -1*m[1,1]
                          #  })
                          #})
#                          names(const_off_proj) = names(train)
                        #      const_off = train[[i]]$mean_x[x$var] %*% betas
                        
                        const =  const_off_proj 
                        self$var = vars
                        self$var_names = c(prev_i$var_names, list(b_i_name))
                        self$tbls = tbls
                        self$varnames=c(prev_i$varnames, varnames)
#                        self$betas_offset=betas_new_offset
                        self$betas_proj=betas_new_proj
                        self$const_off_proj=const_off_proj
                        self$pvs_proj = pvs_proj
                        self$cum_pvs_proj=c(prev_i$cum_pvs_proj,pvs_proj)
                        #self$constants_proj = constants_proj
                        self$constants_proj_all =  prev_i$constants_proj_all
                        self$constants_proj_all[[length( self$constants_proj_all)+1]] = constants_proj
                             self$betas= betas_newP 
                             self$const=const
                             self$W_all = W_all_new
                             self$betasP = betas_newP
                             self$mean_x = mean_x
                             self$mean_y = mean_y
                      }
                    },
                  simplify=function(){
                    names(self$var_names)=self$varnames
                    list(betas=self$betas, constants_proj = self$constants_proj, var_names = self$var_names, pvs = self$pvs_proj)
                  },
                  updateConst=function(phensi,data, train, k,consts_prev=NULL,CHECK=F){   #REVISIT THIS FOR MULTINOM
                    ##need to consider non_na_x
                    #family = data$family
                    nonNA = data$looc$incl[,k]
                    family = unlist(lapply(names(phensi), function(x) strsplit(x,"\\.")[[1]][1]))
                    names(family)=names(phensi)
                    non_na_x = data$getNonNA(self$var) 
                    self$constants_proj = vector("list", length(phensi))
                    names(self$constants_proj) = names(phensi)
                    na_k=non_na_x & nonNA
                    for(kk1 in 1:length(phensi)){
                      phensi1 = phensi[[kk1]]
                      kk = names(phensi)[[kk1]]
#                      nonNA = train$nonNA[[kk]]
                      if(family[[kk]]=="multinomial"){
                        y = attr(data$y[[kk]],"factor")
                        #names(y) = kk
                       levs = levels(y)# c(0,1:length(self$betas[[k]]))
                       yp1 = .calcYpred_multinom(self,  data$data, non_na_x & nonNA,levs, numvar=NULL,kk=kk1,constants=rep(0, ncol(self$betas[[kk]])), liab=F) 
#                       ypred$ypreds[[j]][[kk]][ind_1,] =  .calcYpred_multinom(prev_kj,  self$data, ind_1, levs, numvar = numvar,kk=kk)  ## for multi-prediction
 #                      ypred$ypreds[[j]][[kk]][na_x,] = rep(NA, ncol(ypred$ypreds[[j]][[kk]]))
                       m1 = try(multinom(y[non_na_x & nonNA]~as.matrix(yp1),trace=F))
                       sm  = summary(m1, digits=3)
                       const_term = sm$coefficients[,1]
                       self$constants_proj[[kk1]] = const_term
                      }else if(family[[kk]]=="ordinal"){
                        y = data$y[[kk]]
                        for(kk_1 in 1:length(phensi1)){
                          kk_ = phensi1[kk_1]
                          levs1 = min(y[,kk_], na.rm=T):max(y[,kk_],na.rm=T)
                          consts = rep(0, length(levs1)-1)
                          names(consts) = levs1[-length(levs1)]
                       yp1 = .calcYpred_ord(self,  data$data, non_na_x & nonNA, levs = levs1,numvar=NULL,kk=kk_1,
                                            constants=consts, liab=F) 
                       #yp1 = .calcYpred_1(self,  data$data, non_na_x, kk=1,numvar=NULL,constants=rep(0, ncol(data$y)))
                      
                          y1c = (y[non_na_x & nonNA,kk_])
                          
                          
                        df = data.frame(list(y=factor(y1c,levels=sort(unique(y1c))),x= yp1))
                       # df = data.frame(list(y = y1c,x=yp1))
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
                       #y = data$y[[kk]]
                       spike_slab_iter = getOption("spike_slab_iter",0)
                       for(kk_1 in 1:length(phensi1)){
                         kk_ = phensi1[kk_1]
                         yp1 = .calcYpred_1(self,  data$data, non_na_x &nonNA, kk=kk_1,constants=0) #rep(0, nrow(train$yTr[[kk]]))) 
                          y=data$y[[kk]][non_na_x & nonNA,kk_1]
                       if(spike_slab_iter>1){
                        # print("using spike slab for const")
                         if(family[[kk]]=="binomial"){
                           ab=logit.spike(y~yp1, niter= spike_slab_iter, ping =-1)
                         }else{
                           if(getOption("logprint",F)) print("gaussian spike slab")
                           ab=lm.spike(y~yp1, niter= spike_slab_iter, ping =-1)#,weights=w)
                           
                         }
                         sm = summary(ab)
                         x_ind = match(c("x","(Intercept)"),dimnames(sm$coefficients)[[1]])
                         coeff = sm$coefficients[x_ind[1],c(1,2,3,5)]
                         self$constants_proj[[kk1]][kk_1]= sm$coefficients[x_ind[2],1]
                         #const_term 
                       }else{
                       
                         if(TRUE){
                           self$constants_proj[[kk1]][kk_1] = tryCatch({
                           ones = rep(1, length(yp1[,1]))
                           ridge=glmnet(cbind(ones,yp1[,1]),data$y[[kk]][non_na_x & nonNA,kk_],family=family[[kk]], alpha = 0)
                           rbeta <- coef(ridge,s=min(ridge$lambda))
                          rbeta[1,1]
                           }, error=function(w) {
                             gl = glm(data$y[[kk]][non_na_x & nonNA,kk_]~ yp1[,1], family=family[[kk]])
                             gl$coefficients[1]
                           })
                         }else{
                        self$constants_proj[[kk1]][kk_1]  <- tryCatch({
                          gl = glm(data$y[[kk]][non_na_x & nonNA,kk_]~ yp1[,1], family=family[[kk]])
                          gl$coefficients[1]
                        }, warning=function(w) {
                          print("using glmnet to regularise ")
                          ones = rep(1, length(yp1[,1]))
                          ridge=glmnet(cbind(ones,yp1[,1]),data$y[[kk]][non_na_x & nonNA,kk_],family=family[[kk]], alpha = 0)
                          rbeta <- coef(ridge,s=min(ridge$lambda))
                          rbeta[1,1]
                        })
                         }
                        
                        
                        
                       }
                       }
                     }
                  
                    if(CHECK && FALSE ){
                      yp2 = .calcYpred_1(self,  data$data, non_na_x, numvar=NULL,constants=unlist(self$constants_proj) )
                        gl = glm(data$y[[kk]][non_na_x,kk_]~ yp2[,kk_], family=train$family)
                        if(abs(gl$coefficients[1])>1e-10) stop(" coefficients not right")
                    }
          names(  self$constants_proj) = names(phensi)
                    }
                  }
                    
                  )
)
