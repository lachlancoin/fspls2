##this is cut down version for storing 
trainObj1<-R6Class("trainObj1",
                   public = list(
                     rmsv_ ="matrix",
                     prev="list",
                     prev_old="list",
                     looc_incl_k_ij = "vector",
                     # ypred="ypredObj",
                     
                     initialize = function( train){
                       self$rmsv_ = train$rmsv_
                       self$prev = train$prev
                       self$looc_incl_k_ij = train$looc_incl_k_ij
                       self$prev_old = train$prev_old
                       self$ypred = train$ypred$clone()
                     }
                   )
)

trainObj<-R6Class("trainObj",
                  public = list(
                    initialize = function( data,k,var=list(), varnames = list(), W_all = matrix(nrow=0, ncol=0)){
                      types_ =     getOption("fspls.types", fromJSON('{"gaussian": "rank_correlation","binomial" : "AUC"}'))
                      family = data$family
                      looc_incl_k_ij=data$looc$incl[,k]
                      if(is.null(looc_incl_k_ij)) stop("looc is null")
                      #data = datas[[ij]]
                      yTr = data$y
                      self$family=family
                      nonNA = lapply(1:ncol(yTr), function(ik){
                        n_na = !is.na(yTr[,ik]) #apply(yTr[,ik], 1, function(v) length(which(is.na(v)))==0)
                       #na_y = length(which(n_na))>0
                       n_na & looc_incl_k_ij
                      })
                      nmes1 = names(data$data); names(nmes1) = nmes1
                      if(length(nmes1)==0) stop("no nmes")
                      
                      #}
                      n = nrow(data$y)
                      #self$ypred=ypredObj$new(data,NULL, family)
                    #  phensi=1:length(self$ypred$ypreds)
                     # self$rmsv_=  .calcRMSV_1(phensi,self$ypred$ypreds, data,nonNA, length(self$ypred$ypreds), types_=types_,family=family)
                      
                                   
                      ##should move cols_incl to dataObj
                      if(FALSE){
                        yvar=apply(yTr,2,function(v1) all(duplicated(v1[!is.na(v1)])))
                        ycols = which(!yvar)
                        if(length(ycols)==0){
                          warning(" all phenotype cols have zero variance")  
                          return (NULL)
                        }
                      }else{
                        ycols = 1:ncol(yTr)
                      }
#                      y1 = yTr[nonNA,ycols,drop=F]
                     
                      #self$y2 = yTr[,ycols,drop=F]
                      y1 = vector("list", ncol(yTr))
                      names(y1) = names(yTr)
                      maxy1 =vector("list",ncol(yTr))
                      mean_y=vector("list",ncol(yTr))
                      for(colk in 1:length(family)){
                       
                       if(family[[colk]]=="multinomial"){
                        levs1 = levels(yTr[,colk])
                        names(levs1) = paste(names(y1)[colk],levs1)
                        y1[[colk]] = data.frame(lapply(levs1, function(lv){
                          y3 = rep(0,nrow(yTr))
                          y3[is.na(yTr[,colk])]=NA
                          y3[yTr[,colk]==lv]=1
                          y3
                        }))
                       
                       }else{
                        y1[[colk]] = yTr[,colk,drop=F]
                       }
                        maxy1[[colk]]=apply(y1[[colk]][nonNA[[colk]],,drop=F],2,max)
                        maxy1[[colk]][maxy1[[colk]]==0]=1  
                        d_w = data$weights[nonNA[[colk]],colk]
                        mean_y[[colk]] = rep(0, ncol(y1[[colk]]))
                        for(jk in 1:ncol(y1[[colk]])){
                          mean_y[[colk]][[jk]] = (y1[[colk]][nonNA[[colk]],jk]%*% d_w)/sum(d_w) 
                        }
                      }
                    
#                      d_weights = lapply(nonNA, function(n_na) data$weights[n_na,ycols,drop=F])
                      #if(family=="binomial" && max(maxy1)>1) stop(" not binomial y data")
                      
                      #3if(family!="multinomial"){
                     
                      #}
                      
                      y=y1
#                      ymod = vector("list", ncol(yTr))   
                      for(colk in 1:length(y1)){
                        for(j in 1:length(mean_y[[colk]])){
                          y[[colk]][,j] = data$weights[,colk]*(y1[[colk]][,j]  - mean_y[[colk]][j]) #yTr[,j]  - mean_y[j]
                        }
                        y[[colk]][!nonNA[[colk]],] = rep(0, ncol(y[[colk]]))
                        if(family[[colk]] %in% c("binomial","multinomial")){
                        ##so 0 1 values maxy should be 1
                        ymean = mean_y[[colk]] #mean(y, na.rm=T)  ##NEED TO ADD BACK Y MEAN HERE
                        ymean1 = log(ymean) - log(maxy1[[colk]][[1]]-ymean)
                        #        ymean1 = log(ymean) - log(1-ymean)
 #                       ymod[[colk]] =  t(as.matrix(data$weights[,colk]*((y[[colk]]+ymean)-1/(1+exp(-ymean1)))))
                      }else{
  #                      ymod[[colk]]=t(as.matrix(data$weights[,colk]*(y[[colk]])))
                      }
                      #if(family!="multinomial"){
   #                   mean_check= apply(ymod[[colk]][,nonNA[[colk]],drop=F],1,mean)
  #                    if(max(abs(mean(mean_check)))>1e-5){
   #                     stop(paste("assumption of mean 0 ymod", mean_check))
    #                  }
                      }
                      
                      self$yTr = lapply(y, t)  #should this be ymod??
                      ymean = lapply(self$yTr,function(yTr1) apply(yTr1,1,mean)) ## should be mean 0
                      if(max(abs(unlist(ymean)))>1e-5)stop("problem")
                     # self$y1 = y1
                    #  self$ymod = ymod
                    #  self$weights = d_weights
                  #    self$norm = norm
                      self$nonNA = nonNA
                     # self$min_na = min(which(nonNA))
                      #self$max_na = max(which(nonNA))
                    #  self$na_y =na_y
                      #self$non_na_inds = lapply(nonNA, function(n_na) rep(T,length(which(nonNA)))
                      self$looc_incl_k_ij = looc_incl_k_ij
                   #   self$mean_x = mean_x
                      self$mean_y = mean_y
                      self_ynme = dimnames(data$y)[[2]]
                      #                      list(y=y, y1=y1, ymod=ymod,weights = d_weights,  norm = norm,nonNA = nonNA, 
                      #                           min_na = min(which(nonNA)),max_na = max(which(nonNA)),na_y = na_y,
                      #                           non_na_inds =rep(T,length(which(nonNA))),
                      #                           looc_incl_k =looc_incl_k, mean_x = mean_x, mean_y = mean_y, ynme=dimnames(data$y)[[2]])
                      
                      
                      
                      self$prev = stateObj$new(phensi, data, self,k,var=var, varnames=varnames, W_all = W_all)
                    
                    #  self$prev= prevs1
                    #  self$prev_old=NULL
                      
#                      yTr = lapply(self$y, t)   #precalculate products
                    #  print("precalc products")
                      self$products= lapply(1:length(data$data), function(ik){
                        x = data$data[[ik]]
                        lapply(self$yTr, function(yTr1){
                         # print(dim(yTr1))
                        #  print(dim(x))
                          resu1=if(isbigmatrix(x)) dgemm(A=(yTr1),B=x) else yTr1 %*% x
                          dimnames(resu1) = list(dimnames(yTr1)[[1]],dimnames(x)[[2]])
                          resu1
                        })
                      })
                      names(self$products) = names(data$data)
                      
                    },
getPvs=function(prev){
  # lapply(self$prevs, .getWeights11_1, pvs=T)
  .getWeights11_1(self$prev,pvs=T)
},
reorder=function(o){
  #o1=match(names(o), names(self$prev))
  self$prev = self$prev[o]
},
keep=function(tokeep1){
  if(length(tokeep1)==0){
    #print("winding back")
    self$prev = self$prev_old;
    self$prev_old=NULL
  }else{
    self$prev=self$prev[tokeep1]
  }
},
## we keep this in case we need to wind back 


getMaxBetaProj=function(){
  #     lapply(self$prevs, function(prevk){
  mabv = lapply(self$prev, function(pk) max(abs(unlist(pk$betas_proj))))
  #   })
  mabv
},                    
                    yTr="list",
          #ymod="list", 
                    #y1="list",
                   # y2="matrix",
                    #weights = "matrix",
               
                    nonNA="list",
                    family="character",
#                    min_na = "numeric",max_na="numeric",
#na_y = "logical",
#                    non_na_inds = "vector",
looc_incl_k_ij = "vector",
                 #   ypred="ypredObj",
                  #  rmsv_ ="matrix",
                    prev="list",
products="list",
                  #  prev_old="list",
          
#                    mean_x = "vector",
mean_y="vector"
                  )
)



