








trainObj<-R6Class("trainObj",
                  public = list(
                    yTr="list",
                    k="numeric",
                    means_y = "list",
                    y1="list",
                    #nonNA="list",
                    family="character",
                    looc_incl = "vector",
                    prev="list",
                    products="list",
                    
                    initialize = function(y, looc, family=names(y)){  #y is a list of sparse matrices
                    #  types_ =     getOption("fspls.types", fromJSON('{"gaussian": "rank_correlation","binomial" : "AUC"}'))
                      self$family=family
                      self$prev=list()
                      self$means_y = list()
                      self$looc_incl=looc$incl
                      self$y1=y
#                      self$yNA = yNA
                      self$yTr = lapply(self$y1,function(y11){
                        t(y11)
                      })
                      self$k=NA
                    },
                    transform=function(weights, k){
                      y1 = self$y1
                      looc_incl_k_ij = self$looc_incl[,k]
                      means_y = lapply(y1, function(y1)rep(0, ncol(y1)))
                      for(colk in 1:length(y1)){
                        ncols=ncol(y1[[colk]])
                        meansy = rep(0, ncols)
                        for(j in 1:ncols){
                             v = y1[[colk]][,j]
                             nonNA1 = !is.na(v) & looc_incl_k_ij
                             d_w = weights[nonNA1]
                             meansy[[j]] = (v[nonNA1]%*% d_w)/sum(d_w) 
                            self$yTr[[colk]][j,] = weights*(v  - meansy[j]) #y[,j]  - mean_y[j]
                            if(length(which(!nonNA1))>0){
                              self$yTr[[colk]][j,!nonNA1] =0  
                            }
                        }
                    #    print(meansy)
                        means_y[[colk]] = meansy
                      }
                      self$means_y[[k]] = means_y
                    },
                   # nonNA1=function(k,y){  ##more efficient way to do this with sparse Matrices
                  #    looc_incl_k_ij=self$looc_incl[,k]
                  #    nonNA = lapply(y, function(yTr){
                  #        t(apply(yTr,2, function(v){
                  #            !is.na(v) & looc_incl_k_ij
                  #          }))
                  #    })
                  #    nonNA
                  #  },
                    update=function(data,k,
                                    var=list(), varnames = list(), W_all = matrix(nrow=0, ncol=0)){
                        ycols = 1:length(self$y1)
                        self$k = k
                      #phensi = ycols
                        #.transform<-function(y1,  nonNA, weights, mean_y){
                      self$transform(data$weights, k)
                      ymean = lapply(self$yTr,function(yTr1) apply(yTr1,1,mean)) ## should be mean 0
                      if(max(abs(unlist(ymean)))>1e-5)stop("problem")
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
                      phens = data$pheno()
                      phensi = data$phensi(phens)
                      self$prev[[k]] = stateObj$new(phensi, data, self,k, self$mean_y,var=var, varnames=varnames, W_all = W_all)
                    },
                    
                    calcMean=function( looc_incl_k_ij, weights){
                      lapply(self$y1, function(yTr){
                        apply(yTr,2, function(v){
                          nonNA = !is.na(v) & looc_incl_k_ij
                          d_w = weights[nonNA]
                          #mean_y[[colk]] = rep(0, ncol(y1[[colk]]))
                          #for(jk in 1:ncol(y1[[colk]])){
                          (v[nonNA]%*% d_w)/sum(d_w) 
                          #}
                        })
                      })
                      
                    },          
                    
getPvs=function(prev){
  # lapply(self$prevs, .getWeights11_1, pvs=T)
  .getWeights11_1(prev,pvs=T)
},
reorder=function(o,k){
  #o1=match(names(o), names(self$prev))
  self$prev = self$prev[[k]][o]
},
#keep=function(tokeep1){
#  if(length(tokeep1)==0){
    #print("winding back")
#    self$prev = self$prev_old;
##    self$prev_old=NULL
#  }else{
#    self$prev=self$prev[tokeep1]
#  }
#},
## we keep this in case we need to wind back 


getMaxBetaProj=function(k){
  #     lapply(self$prevs, function(prevk){
  mabv = lapply(self$prev[[k]], function(pk) max(abs(unlist(pk$betas_proj))))
  #   })
  mabv
}                    
 
                  )
)



