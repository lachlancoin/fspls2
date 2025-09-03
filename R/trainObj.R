








trainObj<-R6Class("trainObj",
                  public = list(
                    yTr="list",
                    k="numeric",
                    means_y = "list",
                    y1="list",
                    #nonNA="list",
                    family="character",
                    looc_incl = "vector",
                    products="list",
                    func_str="character",
                    phens1 = "list",
                    incl="list",
                    subphens ="list",
                  #  funcs="closure",
                  #  func_str = "character",
                  #  ifuncs="list",
                    
                    initialize = function(y, looc, incl, funcst,family=names(y)
                                          ){  #y is a list of sparse matrices
                      if(length(names(funcst))==0) stop("no names on transform")
                      
                    #  types_ =     getOption("fspls.types", fromJSON('{"gaussian": "rank_correlation","binomial" : "AUC"}'))
                      self$family=family
                      self$func_str = lapply(funcst, function(xx) xx[[1]])
                      self$means_y = list()
                      self$looc_incl=looc$incl
                      self$y1=y
                      self$yTr =
                        lapply(self$y1,function(y11){
                          lapply(funcst, function(f_k){
                            y12 = t(y11)
                          })
                          })
                      self$k=NA
                      self$incl = incl
                      
                    },
                    product=function(ik,ii, phensi1){ #  self$train$products[[ik]][[ii]][phensi1,,drop=F]  #[,self$cols_incl[[ik]],drop=F]
                     lapply(self$products[[ik]][[ii]], function(p1){
                      p1[match( names(phensi1),dimnames(p1)[[1]]),,drop=F]
                     })
                    },
                    transform=function(weights, k,  phens1){
                       funcst=self$func_str
                     # func_str1 = paste("function(x)",func_str)
                      y1 = self$y1
                      looc_incl_k_ij = self$looc_incl[,k]
                      means_y = lapply(y1, function(y1_) lapply(funcst, function(fst) rep(0, ncol(y1_))))
                      inds_to_do = which(names(self$y1) %in% names(phens1))
                      for(colk1 in 1:length(inds_to_do)){
                       colk = inds_to_do[colk1]
                       ncols=ncol(y1[[colk]])
                       inds_to_do_1 = which(dimnames(y1[[colk]])[[2]] %in% phens1[[colk1]])
                       for(f_k in 1:length(funcst)){
                          funcs =  eval(str2lang(funcst[[f_k]]))
                          meansy = rep(0, ncols)
                          for(j in inds_to_do_1){
                             v = funcs(y1[[colk]][,j])
                             nonNA1 = !is.na(v) & looc_incl_k_ij
                             d_w = weights[nonNA1]
                             meansy[[j]] = (v[nonNA1]%*% d_w)/sum(d_w) 
                             self$yTr[[colk]][[f_k]][j,] = weights*(v  - meansy[j]) #y[,j]  - mean_y[j]
                                if(length(which(!nonNA1))>0){
                                  self$yTr[[colk]][[f_k]][j,!nonNA1] =0  
                                }
                             }
                        
                    #    print(meansy)
                        means_y[[colk]][[f_k]] = meansy
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
                    update=function(data,k,subphens){
                      if(toJSON(subphens)==toJSON(self$subphens) && k == self$k){
                        print("not updating"); return(NULL)
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
                     
                        funcst=self$func_str
                      #  if(!is.na(k)){
                      #   if(k[1]==self$k[1])  return(NULL)
                      #  }
                        self$k = k
                          # no need to update
                      #phensi = ycols
                        #.transform<-function(y1,  nonNA, weights, mean_y){
                      self$transform(data$weights, k, phens1)
                      inds1 = which(names(self$yTr) %in% names(phens1))
                   #   for(i in 1:length(self$funcs)){
                      nmes_phens1 = names(phens1); names(nmes_phens1) = nmes_phens1;
                      nmes_funcst = names(funcst);names(nmes_funcst)=nmes_funcst
                      ymean = lapply(nmes_phens1, function(nme_p1){
                        lapply(nmes_funcst, function(nme_f1){
                            yTr1 = self$yTr[[nme_p1]][[nme_f1]]
                            subinds = dimnames(yTr1)[[1]] %in% unlist(phens1) #[[nmes_phens1]]
                            if(length(which(subinds))==0) subinds = dimnames(yTr1)[[1]] %in% phens1[nmes_phens1]
                            apply(yTr1[subinds,,drop=F],1,mean,na.rm=T)
                        })
                      })
                       incl1 = names(data$data); names(incl1) = incl1
                      self$products= lapply(incl1, function(ik){
                        if(!(ik %in% incl)) return(NULL)
                        x = data$data[[ik]]
                        lapply(nmes_phens1, function(nme_p1){
                          lapply(nmes_funcst, function(nme_f1){
                          yTr1 = self$yTr[[nme_p1]][[nme_f1]]
                          subinds = dimnames(yTr1)[[1]] %in% unlist(phens1)#[[nmes_phens1]]
                          if(length(which(subinds))==0) subinds = dimnames(yTr1)[[1]] %in% phens1[[nmes_phens1]]
                          yTr1 = yTr1[subinds,,drop=F]
                         # lapply(yTr_, function(yTr1){ ##NEED TO ADD POWS
                            # print(dim(yTr1))
                            #  print(dim(x))
                            resu1=if(isbigmatrix(x) && typeof(yTr1)!="S4") dgemm(A=yTr1,B=x) else yTr1 %*% x
                            dimnames(resu1) = list(dimnames(yTr1)[[1]],dimnames(x)[[2]])
                            resu1
                          })
                        })
                      })
                     # names(self$products) = names(data$data)
                     # phens = data$pheno()
                    #  phensi = data$phensi(phens)
                   #   self$prev[[k]] = stateObj$new(phensi, data, self,k, self$mean_y,var=var, varnames=varnames, W_all = W_all)
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



