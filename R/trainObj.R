








trainObj<-R6::R6Class("trainObj",
                  public = list(
                    yTr="list",
                    k="numeric",
                    means_y = "list",
                    y1="list",
                    #nonNA="list",
                    family="character",
                    looc_incl = "vector",
                    products="list",
                    transforms="list",
                    phens1 = "list",
                    incl="list",
                    subphens ="list",
                  #  funcs="closure",
                  #  func_str = "character",
                  #  ifuncs="list",
                    
                    initialize = function(y, looc, incl, transforms,family=names(y)
                                          ){  #y is a list of sparse matrices
                      if(length(names(transforms))==0) stop("no names on transform")
                     # print(y)
                      
                    #  types_ =     getOption("fspls.types", fromJSON('{"gaussian": "rank_correlation","binomial" : "AUC"}'))
                      self$family=family
                      self$transforms = transforms
                      self$means_y = list()
                      self$looc_incl=looc$incl
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
                    transform=function(weights, k,  phens1){
                       funcst=self$transforms
                     # func_str1 = paste("function(x)",func_str)
                      y1 = self$y1
                      looc_incl_k_ij = self$looc_incl[,k]
                      means_y = lapply(y1, function(y1_) lapply(funcst, function(fst) lapply(fst$params,function(f3) rep(0, ncol(y1_)))))
                      inds_to_do = which(names(self$y1) %in% names(phens1))
                      #colk1=1; f_k = 1; g_k = 1; colk = inds_to_do[colk1]; j= which(dimnames(y1[[colk]])[[2]] %in% phens1[[colk1]])[1]
                      
                      for(colk1 in 1:length(inds_to_do)){
                       colk = inds_to_do[colk1]
                       ncols=ncol(y1[[colk]])
                       inds_to_do_1 = which(dimnames(y1[[colk]])[[2]] %in% phens1[[colk1]])
                       for(f_k in 1:length(funcst)){
                          funcs =  funcst[[f_k]]$invfunc #could also be invfunc?
                          params = funcst[[f_k]][[3]]; names(params)=params
                          for(g_k in 1:length(params)){
                            pow1 = params[[g_k]]
                              meansy = rep(0, ncols)
                              for(j in inds_to_do_1){
                              #  print("APPLYING FUNC TO Y HERE")
                              #  print(funcs);
                                nonNA1 =  looc_incl_k_ij
                                v = funcs(y1[[colk]][nonNA1,j], pow1) 
                                nonNA2 = !is.na(v)
                               # print(head(v))
                               #print(head(nonNA2));
                                # if(length(which(is.na(v)))> 0.3 *length(v)) warning(paste(" more than 30% NA in transformation", toJSON(funcs)))
                                 d_w = weights[nonNA1][nonNA2]
                                 meansy[[j]] = (v[nonNA2]%*% d_w)/sum(d_w) 
                                 v2 =  weights[nonNA1]*(v  - meansy[j])
                                 v2 = v2/sd(v2,na.rm=TRUE)  ### new line to avoid giving advantage to transformations
                                 self$yTr[[colk]][[f_k]][[g_k]][j,nonNA1] =v2 #y[,j]  - mean_y[j]
                                    if(length(which(!nonNA1))>0){
                                      self$yTr[[colk]][[f_k]][[g_k]][j,!nonNA1] =0   ## will not contribute to dot product
                                   
                                    }
                                 if(length(which(!nonNA2))>0){
                                   self$yTr[[colk]][[f_k]][[g_k]][j,which(nonNA1)[!nonNA2]] =0
                                 }
                                 vars1 = apply(self$yTr[[colk]][[f_k]][[g_k]][j,,drop=FALSE],1,var, na.rm=TRUE)
                                 if(min(vars1)==0) stop(paste(" transformations gave raise to zero variance, choose diff transformations",toJSON(funcs)))
                                 }
                            
                        #    print(meansy)
                            means_y[[colk]][[f_k]][[g_k]] = meansy
                          }
                       }
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
                    update=function(data,k,subphens, force=FALSE){
                      
                      if(!force && toJSON(subphens)==toJSON(self$subphens) && k == self$k){
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
                        self$k = k
                      self$transform(data$weights, k, phens1)
                      inds1 = which(names(self$yTr) %in% names(phens1))
                   #   for(i in 1:length(self$funcs)){
                      nmes_phens1 = names(phens1); names(nmes_phens1) = nmes_phens1;
                      nmes_funcst = names(funcst);names(nmes_funcst)=nmes_funcst
                      ymean = lapply(nmes_phens1, function(nme_p1){
                        lapply(nmes_funcst, function(nme_f1){
                          nme_t1=names(self$transforms[[nme_f1]]$params); names(nme_t1) = nme_t1
                          lapply(nme_t1,function(p1){
                           # print(paste(nme_p1, nme_f1, p1))
                            yTr1 = self$yTr[[nme_p1]][[nme_f1]][[p1]]
                            subinds = dimnames(yTr1)[[1]] %in% unlist(phens1) #[[nmes_phens1]]
                            if(length(which(subinds))==0) subinds = dimnames(yTr1)[[1]] %in% phens1[nmes_phens1]
                            apply(yTr1[subinds,,drop=FALSE],1,mean,na.rm=TRUE)
                          })
                        })
                      })
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
                    
                    calcMean=function( looc_incl_k_ij, weights){
                      if(TRUE) stop(" not used")
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



