
.expandR<-function(U, nonNA){
  U1 = matrix(0,nrow =length(nonNA), ncol = ncol(U))
  U1[nonNA,]=U
  dimnames(U)
  U1
}

UDVPObj<-R6::R6Class("UVDPObj", public = list(
  P="matrix",
  inds = "vector",
  centralise="logical",
  alias = "vector",
  U="matrix",
  Ut="matrix",
  Dinv="matrix",
  Vinv="matrix",
  VDU = "matrix",
  I="matrix",
  P2="matrix",
  var="vector",
  initialize=function(data,var1,
                      P= data$extractData(var1, adjust=TRUE),
                      check=FALSE,centralise=FALSE
                      ){  ## d = train[[i]] ## already centralised
    #d = data$train
    nonNA = rep(TRUE, data$nrow) # use all  data$train$nonNA
   # var_df = data.frame(var1)
    self$var = var1
    #var=var1
  #  x_all = lapply( 1:length(data$data), function(k){
  #    var1_inds = which(var_df[1,]==k)
  ##    var1_ = as.numeric(var_df[2,var1_inds])
  #    x = if(length(var1_)==0) data$data[[k]][, 1,drop=FALSE][,c()] else as.matrix(data$data[[k]][, var1_,drop=FALSE])
  #    if(length(var1_)>0){
  #      for(j in 1:ncol(x)){
  #        x[,j] = x[,j]- data$mean_x[[k]][var1_[j]]  ## changed this so that projection is based on full matrix
  #      }
  #    }
  #    x
  #  })
  #  x_all = x_all[unlist(lapply(x_all, function(ab) !is.null(dim(ab))))]
  #  x11 = as.matrix(data.frame(x_all))
    
    inds = which(nonNA)
    
    if(length(inds)==0) stop(" should not be zero")
    #  P = data[,colinds,drop=FALSE]
    #P = x11
    lengthvars =ncol(P) #ncol(P)
    if(centralise) P= apply(P, 2, .centralise, inds)  ##NEED TO REVISIT CENTRALISATION ON FLY
    
    if(lengthvars==0){
      self$P = P
      self$inds = inds
      self$centralise = centralise
      self$alias = NULL
      self$U = NULL
      self$Ut = NULL
      self$Dinv = NULL
      self$Vinv = NULL
      self$VDU = NULL
      self$P2 = NULL
#      return(list(P=P, inds = inds, alias = alias, centralise = centralise))
    }else{
      P1 = P[inds,,drop=FALSE] 
    
      ##double check column means are zero
      if(check){
        me = apply(P1,2,mean, na.rm=TRUE);
        if(max(abs(me))>1e-5){
         # print(me)
          stop('err')
        }
      }
    #vars = apply(P1,2,var)
    #  if(!is.null(inds) & length(which(inds))>0)  else P
      alias = which(apply(P1,1,.cntNA)==0)
    
    #  svd = dgesdd(A=P1[alias,,drop=FALSE])
      svd = svd( P1[alias,,drop=FALSE])
     U = svd$u
     V = t(svd$v)
    
    if(lengthvars==1) {
      Dinv = as.matrix(1/svd$d)
    }else {
      Dinv = diag(1/svd$d)
    }
    Vinv = solve(V)
    #Dinv = solve(D)
    U_exp =  .expandR(U, nonNA)
    inv1= which(svd$d>0)  #to avoid infinite values
    VDU = Vinv %*% Dinv[,inv1] %*% t(U_exp[,inv1,drop=FALSE])
    
    
    P2 =NULL #diag(nrow(P)) - P%*%VDU  somewhat expensive to calcualte, so we avoid
    self$U =U_exp
    self$Ut = t(self$U)
    self$Dinv = Dinv
    self$Vinv = Vinv
    self$VDU = VDU
    self$P = P #as.big.matrix(P)
    self$P2 = P2
    self$inds = inds
    self$alias = alias
    self$centralise = centralise
    
    }
  },
getW=function(){
  UDV=self
  alias = UDV$alias
  U = UDV$U
  Dinv = UDV$Dinv
  Vinv = UDV$Vinv
  W = Vinv %*% Dinv %*% t(U)
  W
},
getWall=function(x, Wall1){
  UDV=self
  W = Wall1
  if(is.null(Wall1)){
    Wall=matrix(nrow=0, ncol=0)
  }else if(nrow(W)==0){
    Wall=matrix(1)
  }else{
    
    W_h_best_i= UDV$VDU %*%  x    ##UDV$VDU[,nonNA,drop=FALSE] %*%  x  #d$x[,b_i]
    Wall3 = cbind(W,-W_h_best_i[,1])
    Wall=rbind(Wall3,c(rep(0,ncol(Wall3)-1),1))
  }
  Wall
},
P2_old=function(){
  P = self$P
  VDU = self$VDU
  I = diag(nrow(P))
 # P2=diag(nrow(P)) - P%*%VDU
  P2 = dgemm(ALPHA = -1.0, A=P, B = VDU, C = I)
  P2
}
)
)
#.calcUDVP<-function(datas,train, var1){
#  data_inds = 1:length(datas)
#  names(data_inds) = names(datas)
#  UDVP_h=lapply(data_inds, function(i) UDVPObj$new(train[[i]], datas[[i]], var1))
 
#}
