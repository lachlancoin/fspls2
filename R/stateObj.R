
initialStateObj<-function(phens){
  Wall0 =lapply(phens, function(f) matrix(nrow=0,ncol=0))
  stateObj$new(phens,NULL,  NULL,NULL,NULL,NULL,NULL, NULL,NULL, var=c(), varnames=c(), Wall =Wall0)
}

#.calcWall_2<-function(data,  var){
#  prev_i = list(Wall =matrix(nrow=0, ncol=0) , var = c() )
#  for(k in 1:length(var)){
#    Wall2 = .calcWall_1(data, var[[k]], prev_i)
#    prev_i = list(Wall = Wall2, var = var[1:k])
#  }
#  prev_i$Wall
#}
#colk is family col



stateObj<-R6::R6Class("stateObj",##represents a state of the model
                  public = list(
                    var="list", var_names="list",
                    varnames="list",
                    name="character",
                    angle="double",
                    angles="list",
                    logLiks ="list",
                    degsf="list",
                    degsf0="list",
                    lambdas="list",
                    logLiks0="list", ## logLiks of nested model
                    sumAngle="numeric",
                    name_prev="character",
                    betas_proj="vector", #list of vector
                    Wall="matrix",  ##list of matrix,
                    constants_proj="list",
                    constants = "list",
                   betas="list",
                    #nonNA="logical",
                    mean_x="numeric",
                   tbls="list",
                  # pvs="list",pvs_sep="list",
                  # pvs_all ="list",
                  # sumPv="numeric",
                  
                    initialize = function(phensi,data,
                                          betas_proj, 
                                          betas,
                                          constants_proj,
                                          tbls,
                                          logLiks,degsf,
                                          lambdas,
                                      #   k, #mean_y,
                                          prev_i = NULL,
                                          b_i = NULL,b_i_name=NULL,
                                          var = list(), varnames = list(),
                                          Wall = matrix(nrow=0, ncol=0),
                                         mean_x = NULL, 
                                      
                                         useoffset=TRUE
                                         ){
                      self$angle=NULL;
                      self$sumAngle=0;
                      self$angles = c();
                      #self$pvs = pvs
                      self$logLiks = logLiks;
                    
                   #   print("H");print(self$degsf);
                      self$constants_proj=constants_proj
                      self$constants = NULL
                    self$lambdas = lambdas
                      train = data$train
                      self$tbls = tbls
                      self$Wall = Wall
                                       family = unlist(lapply(names(phensi), function(x) getOption("fspls.family",strsplit(x,"\\.")[[1]][1])))
                                       self$degsf=degsf;
                                       
                        
                                       
                      names(family)=names(phensi)
                      self$betas_proj=betas_proj
                      self$betas = betas;
               #       print("HHHHH");print(self$betas)
                      if(is.null(prev_i)){
                        self$logLiks0 = NULL
                        self$degsf0=NULL
                        self$var = var
                        self$var_names = var
                        self$varnames=var
                        self$mean_x=c()
                        self$name=""
                   #     self$cum_pvs_proj = list()
                    #    self$betas_proj= betas_proj #if(length(var)==0) c() else lapply(data$y, function(yy)(matrix(0,nrow=1, ncol= length(var)))) #lapply(datas, function(x) return(c()))
                      #  self$betas = if(length(var)==0) NULL else lapply(data$y, function(yy)(matrix(0,ncol=1, nrow= length(var)))) 
                        self$name=if(length(var)==0) "" else paste(unlist(lapply(var, paste, collapse=".")),collapse="-")
                      }else{
                        self$logLiks0=prev_i$logLiks;
                        self$degsf0=prev_i$degsf;
                        self$name_prev=prev_i$name
                        var = prev_i$var
                        self$mean_x = c(prev_i$mean_x, mean_x)
                        var_st = unlist(lapply(var, paste, collapse="."))
                       
                        if(paste(b_i,collapse=".") %in% var_st) {
                          #print(paste("selected same variable twice"))#,names(todoInds[[b_i[1]]])[b_i[2]]))
                          warning(paste(" duplicated variable")) #,names(todoInds[[b_i[1]]])[b_i[2]]))
                          return(NULL)
                        }
                        varnames=paste(b_i_name, collapse=".")
                        
                        vars = c(var,list(b_i))
                        names(vars)[length(vars)] = varnames
                        self$name=paste(unlist(lapply(vars, paste, collapse=".")),collapse="-")
                        vars0=vars
                        self$var = vars
                        self$var_names = c(prev_i$var_names, list(b_i_name))
                        self$varnames=c(prev_i$varnames, varnames)
                        names(self$var_names)[[length(self$var_names)]] = varnames
                        names(self$varnames)[[length(self$varnames)]] = varnames
                        
                        if(length(self$degsf)==0){
                          self$degsf = unlist(.lapply_nme(phensi, function(x){
                            fam = strsplit(x,"\\.")[[1]][1]
                            if(fam %in% c("multinomial","ordinal")) return(0);
                            return(rep(0, length(phensi[[x]])))
                          }))
                          
                        }
                        if(length(self$logLiks)==0){
                          self$logLiks = unlist(.lapply_nme(phensi, function(x){
                            fam = strsplit(x,"\\.")[[1]][1]
                            if(fam %in% c("multinomial","ordinal")) return(-Inf);
                            return(rep(-Inf, length(phensi[[x]])))
                          }))
                          
                        }
                        

                      }
#print(self$betas)
                    },
                  calcPvalue=function(log.p=TRUE){
                    if(is.null(self$logLiks0)) return(NULL);
                    .lapply_nme(self$logLiks, function(nme){
                      .lapply_nme(self$logLiks[[nme]], function(nme1){
                        ll = self$logLiks[[nme]][[nme1]];
                        ll0 = self$logLiks0[[nme]][[nme1]];
                        df = self$degsf[[nme]][[nme1]];
                        df0 = self$degsf0[[nme]][[nme1]]
                        .lrt(ll,ll0,df,df0, log.p=log.p)
                      })
                    })
                  },
    
                  setOffset=function(){
                    nme_betas_new = names(self$betas_proj)
                    self$constants = self$constants_proj
                    #b_n= nme_betas_new[[1]]
                    for(b_n in nme_betas_new){
                      offset =self$mean_x%*%  self$betas[[b_n]]
                      vv = unlist(self$constants_proj[[b_n]])-offset[1,]
                      self$constants[[b_n]] = vv  ## this was constants_proj
                    }
                  },
                   # self$betas=lapply(nme_betas_new, function(b_n) {
                   #    bp = self$betas_proj[[b_n]]
                   #    Wall_new %*% bp
                   #  })
##simplifies and translates into original space
                  simplify=function(){
                     
                   # self$setOffset() 
                    names(self$var_names)=self$varnames
                    list(betas=self$betas, constants_proj = self$constants_proj,
                         constants = self$constants,
                         logLiks = self$logLiks, 
                         logLiks0 = self$logLiks0,
                         degsf = self$degsf,
                         degsf0 = self$degsf0,
                         lambdas = self$lambdas,
                         #pvs = self$pvs, pvs_all=self$pvs_all,pvs_sep = self$pvs_sep,
                      #   cum_pv= .sumChisq(unlist(self$pvs)),
                       #  cumpv_all= .sumChisq(unlist(self$pvs_all)),
                         betas_proj = self$betas_proj,
                         Wall = self$Wall,
                         mean_x = self$mean_x,tbls = self$tbls,
                     #    Wall = self$Wall,
                         var_names = self$var_names) #, pvs = self$pvs_proj)
                  },
                  updateConst=function(phensi,ypred, data,  k,useglm=FALSE, verbose=FALSE,update=FALSE
                                       ){
                    #if(length(phensi)>1) stop("!!")
                    family = unlist(lapply(names(phensi), function(x) getOption("fspls.family",strsplit(x,"\\.")[[1]][1])))
                    names(family)=names(phensi)
                    non_na_x = data$getNonNA(self$var) 
                    constants_proj = vector("list", length(phensi))
                    names(constants_proj) = names(phensi)
                    na_k=non_na_x & data$looc$incl[,k];
                  for(kk1 in 1:length(phensi)){
                    if(verbose)print(names(phensi)[[kk1]])
                    #kk1 = 1;  
                    kk = names(phensi)[[kk1]]
                    phensi1 = phensi[[kk1]]
                    y=  if(family[[kk]]=="multinomial")   attr(data$y[[kk]],"factor")[na_k] else data$y[[kk]][na_k,,drop=FALSE]
                    yp1 =ypred$ypreds[[kk1]][na_k,,drop=FALSE]
                    if(family[[kk]]=="multinomial"){
                        levs = levels(y)# c(0,1:length(self$betas[[k]]))
                        m1 = try(multinom(y~as.matrix(yp1[,-1]),trace=FALSE))
                        sm  = summary(m1, digits=3)
                        if(verbose)   print(sm$coefficients)
                        
                        
                        const_term = sm$coefficients[,1]
                        constants_proj[[kk1]] = const_term
                    }else if(family[[kk]]=="ordinal"){
                        for(kk_1 in 1:length(phensi1)){
                          kk_ = phensi1[kk_1]
                          y1c = y[,kk_]
                          levs1 = min(y1c, na.rm=TRUE):max(y1c,na.rm=TRUE)
                          consts = rep(0, length(levs1)-1)
                          names(consts) = levs1[-length(levs1)]
                          df = data.frame(list(y=factor(y1c,levels=sort(unique(y1c))),x= yp1[,1]))
                          m1=try(polr(y~x,  data=df,Hess=TRUE, method="logistic"))  
                 
                          if(abs(m1$coefficients[[1]]-1)>0.5){
                           # print(m1$coefficients)
                           # print(self$var)
                            warning(" something gone wrong")
                          }
                          if(inherits(m1,"try-error")){
                            stop("polr error in updating constant in stateObj")
                            ##prob not a great way to do this but avoids errors
                           # gl = glm(data$y[non_na_x & nonNA,kk]~ yp1[,1], family="gaussian")
                            #self$constants_proj[[kk]][[kk_1]] = consts_prev[[kk]] #rep(NA, length(levs1)-1) #gl$coefficients[1]
                          }else{
                            if(verbose) print(m1$coefficients)
                             
                            if(verbose)  print(m1$zeta) 
                            
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
                          # if(getOption("logprint",FALSE)) print("gaussian spike slab")
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
                           ridge=glmnet(cbind(ones[nonNAy],yp1[nonNAy,kk_1]),y1c[nonNAy],family=family[[kk]], alpha = 0,
                                        standardize=FALSE,
                                        lambda = getOption("lambda"))
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
                         warning("using glmnet to regularise ")
                          #ones = rep(1, length(yp1[,1]))
                          nonNAy = !is.na(y1c)
                          ridge=glmnet(cbind(ones[nonNAy],yp1[nonNAy,kk_1]),y1c[nonNAy],family=family[[kk]], alpha = 0,lambda = getOption("lambda"))
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
