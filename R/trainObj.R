








trainObj<-R6::R6Class("trainObj",
                  public = list(
                    yTr="list",
                    k="numeric",
                    means_y0 = "list", ## before transform
                    means_y1 = "list", ## before transform
                    
                    counts_y="list",
                    y1="list",
                    #nonNA="list",
                    family="character",
                    looc_incl_k_ij = "vector",
                    products="list",
                    transforms="list",
                    phens1 = "list",
                    incl="list",
                    subphens ="list",
                  #  funcs="closure",
                  #  func_str = "character",
                  #  ifuncs="list",
                    
                    initialize = function(y, weights, looc, incl, transforms,family=names(y)
                                          ){  #y is a list of sparse matrices
                      if(length(names(transforms))==0) stop("no names on transform")
                     # print(y)
                      
                    #  types_ =     getOption("fspls.types", fromJSON('{"gaussian": "rank_correlation","binomial" : "AUC"}'))
                      self$family=family
                      self$transforms = transforms
                      self$means_y0 = NULL;  self$means_y1 = NULL;
                      self$counts_y = NULL;
                     # self$looc_incl_k_ij=  self$looc_incl[,k]
                      self$looc_incl_k_ij=looc;
                      self$y1=y
                      self$yTr =
                        lapply(self$y1,function(y11){
                          lapply(transforms, function(f_k){
                            params = f_k$params;names(params)=params
                            lapply(params, function(p1){
                               Matrix::t(y11)
                            })
                          })
                          })
                      self$k=NA
                      self$incl = incl
                      self$calcMeans(weights)
                    },
                    diffTransforms=function(transforms){
                      toJSON(names(transforms))!= toJSON(names(self$transforms))
                    },
                    product=function(ik,ii, phensi1){ #  self$train$products[[ik]][[ii]][phensi1,,drop=FALSE]  #[,self$cols_incl[[ik]],drop=FALSE]
                     produ=lapply(self$products[[ik]][[ii]], function(p0){
                       lapply(p0, function(p1){
                          p1[match( names(phensi1),dimnames(p1)[[1]]),,drop=FALSE]
                     })
                     })
                     produ
                    },
                  calcMeans=function(weights){  ## has same patter as transform but split out so that we can average across datasets
                    #print("calc means")
                    funcst=self$transforms
                    y1 = self$y1
                    looc_incl_k_ij = self$looc_incl_k_ij;
                  
                    means_y0 = lapply(y1, function(y1_) lapply(funcst, function(fst) lapply(fst$params,function(f3) rep(0, ncol(y1_)))));
                    means_y1 = lapply(y1, function(y1_) lapply(funcst, function(fst) lapply(fst$params,function(f3) rep(0, ncol(y1_)))));
                    
                    counts_y = lapply(y1, function(y1_) lapply(funcst, function(fst) lapply(fst$params,function(f3) rep(0, ncol(y1_)))))
                 #   inds_to_do =1:length(self$y1);## which(names(self$y1) %in% names(phens1))
                  #  if(length(inds_to_do)==0) stop("inds_to_do is empty")
                    for(colk in 1:length(self$y1)){
                      #colk = inds_to_do[colk1]
                      #print(colk)
                      ncols=ncol(y1[[colk]])
                      inds_to_do_1 = 1:ncol(y1[[colk]])
#                      inds_to_do_1 = which(dimnames(y1[[colk]])[[2]] %in% phens1[[colk1]])
                      for(f_k in 1:length(funcst)){
                        funcs =  funcst[[f_k]]$invfunc #could also be invfunc?
                        params = funcst[[f_k]][[3]]; names(params)=params
                        for(g_k in 1:length(params)){
                          pow1 = params[[g_k]]
                          meansy = rep(0, ncols)
                          meansy1 = rep(0, ncols)
                          county = rep(0, ncols)
                          for(j in inds_to_do_1){
                            nonNA1 =  looc_incl_k_ij
                            if(getOption("mean_before_transf",TRUE)){
                              v = y1[[colk]][nonNA1,j] 
                              nonNA2 = !is.na(v)
                              d_w = weights[nonNA1][nonNA2]
                              meansy[[j]] = (v[nonNA2]%*% d_w)/sum(d_w) 
                            }
                              v = funcs(y1[[colk]][nonNA1,j]-meansy[[j]], pow1) 
                              nonNA2 = !is.na(v)
                              d_w = weights[nonNA1][nonNA2]
                               meansy1[[j]] = (v[nonNA2]%*% d_w)/sum(d_w) 
                            
                            county[[j]] = length(which(nonNA2));
                          }
                          means_y0[[colk]][[f_k]][[g_k]] = meansy
                          means_y1[[colk]][[f_k]][[g_k]] = meansy1
                          counts_y[[colk]][[f_k]][[g_k]] = county
                        }
                      }
                    }
                    self$means_y0 = means_y0
                    self$means_y1 = means_y1
                    self$counts_y = counts_y;
                  },
                    transform=function(weights,   means_y0 = self$means_y0, means_y1 = self$means_y1){
                      funcst=self$transforms
                      y1 = self$y1
                      if(is.null(means_y0)) {
                        stop("need to calc means first")
                       
                      }
                     
                      
                     # print("transform")
                      looc_incl_k_ij = self$looc_incl_k_ij
                    #  inds_to_do = which(names(self$y1) %in% names(phens1))
                  #    if(length(inds_to_do)==0) stop("inds_to_do is empty")
                     # for(colk1 in 1:length(inds_to_do)){
                      for(colk in 1:length(self$y1)){
                      # colk = inds_to_do[colk1]
                       ncols=ncol(y1[[colk]])
                       inds_to_do_1 = 1:ncol(y1[[colk]])#which(dimnames(y1[[colk]])[[2]] %in% phens1[[colk1]])
                       for(f_k in 1:length(funcst)){
                          funcs =  funcst[[f_k]]$invfunc 
                          params = funcst[[f_k]][[3]]; names(params)=params
                          for(g_k in 1:length(params)){
                            pow1 = params[[g_k]]
                              meansy0 = means_y0[[colk]][[f_k]][[g_k]]
                              meansy1 = means_y1[[colk]][[f_k]][[g_k]]
                              
                              for(j in inds_to_do_1){
                                nonNA1 =  looc_incl_k_ij
                               
                                    v = funcs(y1[[colk]][nonNA1,j]-meansy0[[j]], pow1)  ## should we subtract mean before transformation  
                                    nonNA2 = !is.na(v)
                                    v2 =  weights[nonNA1]*(v  - meansy1[[j]])
                                
                                   #  d_w = weights[nonNA1][nonNA2]
                                    # meansy[[j]] = (v[nonNA2]%*% d_w)/sum(d_w) 
                                    
                                     v2 = v2/sd(v2,na.rm=TRUE)
                                 ### new line to avoid giving advantage to transformations 
                                 self$yTr[[colk]][[f_k]][[g_k]][j,nonNA1] =v2 #y[,j]  - mean_y[j]
                                    if(length(which(!nonNA1))>0){
                                      self$yTr[[colk]][[f_k]][[g_k]][j,!nonNA1] =0   ## will not contribute to dot product
                                   
                                    }
                                 if(length(which(!nonNA2))>0){
                                   self$yTr[[colk]][[f_k]][[g_k]][j,which(nonNA1)[!nonNA2]] =0
                                 }
                                 vars1 = apply(self$yTr[[colk]][[f_k]][[g_k]][j,,drop=FALSE],1,var, na.rm=TRUE)
                                 if(min(vars1)==0) {
                                   if(min(apply(self$y1[[colk]][,j,drop=FALSE],2,var,na.rm=TRUE))>0){
                                     if(getOption("verbose",FALSE))print(self$y1[[colk]][,j])
                                     warning(paste(" transformations gave raise to zero variance, choose diff transformations",toJSON(funcs),pow1,
                                                   colnames(y1[[colk]])[[inds_to_do_1[j]]],
                                                   sep="\n"))
                                   }

                                 }
                                 }
                            
                          }
                       }
                      }
                      #self$means_y[[k]] = means_y
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
                  
                    update=function(data,subphens, means_y = list(y0=self$means_y0, y1=self$means_y1), force=FALSE){
                      
                      if(!force && toJSON(subphens)==toJSON(self$subphens)){
                        if(length(which(unlist(lapply(self$products, is.null))))==0){
                        ##print("not updating"); 
                          return(NULL)
                        }
                      }
                      self$subphens = subphens
                      incl = self$incl
                      phen_fams = names(subphens) #unlist(lapply(phens, function(ph) names(ph)))
                      phen_fam = unique(phen_fams)
                      names(phen_fam) = phen_fam
                      phens1 = lapply(phen_fam, function(pf){
                        ab = unique(unlist(subphens[which(phen_fams==pf)])  )
                        names(ab)=ab
                        ab
                      })
                      self$phens1 = phens1
                        ycols = which(names(self$y1) %in% names(phens1))
                        funcst=self$transforms
#                        self$k = k
                     
                      self$transform(data$weights,  means_y0 = means_y[[1]], means_y1 = means_y[[2]])
                      inds1 = which(names(self$yTr) %in% names(phens1))
                   #   for(i in 1:length(self$funcs)){
                      nmes_phens1 = names(phens1); names(nmes_phens1) = nmes_phens1;
                      nmes_funcst = names(funcst);names(nmes_funcst)=nmes_funcst
                    #  ymean = lapply(nmes_phens1, function(nme_p1){
                    #    lapply(nmes_funcst, function(nme_f1){
                    #      nme_t1=names(self$transforms[[nme_f1]]$params); names(nme_t1) = nme_t1
                    #      lapply(nme_t1,function(p1){
                           # print(paste(nme_p1, nme_f1, p1))
                    #        yTr1 = self$yTr[[nme_p1]][[nme_f1]][[p1]]
                    #        subinds = dimnames(yTr1)[[1]] %in% unlist(phens1) #[[nmes_phens1]]
                    #        if(length(which(subinds))==0) subinds = dimnames(yTr1)[[1]] %in% phens1[nmes_phens1]
                    #        apply(yTr1[subinds,,drop=FALSE],1,mean,na.rm=TRUE)
                    #      })
                    #    })
                    #  })
                       incl1 = names(data$data); names(incl1) = incl1
                      #ik = incl1[[1]]; nme_p1 = nmes_phens1[[1]]; nme_f1 = nmes_funcst[[1]]; p1 = names(self$transforms[[nme_f1]]$params)[[1]]
                      self$products= lapply(incl1, function(ik){
                        if(!(ik %in% incl)) return(NULL)
                        x = data$data[[ik]]
                        lapply(nmes_phens1, function(nme_p1){
                          lapply(nmes_funcst, function(nme_f1){
                            nme_t1=names(self$transforms[[nme_f1]]$params); names(nme_t1) = nme_t1
                            lapply(nme_t1,function(p1){
                          yTr1 = self$yTr[[nme_p1]][[nme_f1]][[p1]]
                          subinds = dimnames(yTr1)[[1]] %in% unlist(phens1)#[[nmes_phens1]]
                          if(length(which(subinds))==0) subinds = dimnames(yTr1)[[1]] %in% phens1[[nmes_phens1]]
                          yTr1 = yTr1[subinds,,drop=FALSE]
                          resu1=if(isbigmatrix(x) && typeof(yTr1)!="S4") dgemm(A=yTr1,B=x) else yTr1 %*% x
                            dimnames(resu1) = list(dimnames(yTr1)[[1]],dimnames(x)[[2]])
                            resu1
                            })
                          })
                        })
                      })
                    },
                    
                           
                    
getPvs=function(prev){
  # lapply(self$prevs, .getWeights11_1, pvs=TRUE)
  .getWeights11_1(prev,pvs=TRUE)
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



