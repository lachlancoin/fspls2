#.mergeRMSV<-function(datas){
#  rmsv_=.merge1(lapply(datas, function(d)d$rmsv_),num_cols="value",addName="data")
#  rmsv1 = .findMinRMSV(rmsv_, params$mult, fspls.sum=getOption("fspls.sum",T))
#  min(rmsv1)
#}

calcVar<-function(yp,y,family="gaussian"){
  aa=try(ci_var(yp-y))
  vy = var(y)
  if(inherits(aa,"try-error")){
    res= (c(NA,var(yp-y)/vy,NA))
  }else{
  res = c(aa$interval[1], aa$estimate, aa$interval[2])/vy
  }
  names(res) =c("low","mid","high")
  res
}

.corCI<-function(x,y, method="pearson", probs = c(0.025,0.975)){
  type=if(method=="spearman") "bootstrap" else "normal"
  if(length(x)<5){
      res = c(-1, cor(x,y,method=method),1)
  }else{
    cor1 = try(ci_cor(x,y,method=method, probs=probs, type=type))
  
      if(inherits(cor1,"try-error")){
        res = c(-1, cor(x,y,method=method),1)
      }else{
        range=cor1$interval
        res = c(range[1], cor1$estimate, range[2])
      }
  }
  names(res) =c("low","mid","high")
  res
}

.prob<-function(vec)vec/sum(vec)
.samp<-function(pr, levs){
  ind = which(pr==max(pr,na.rm=T))
  return(levs[ind[1]])
  #sample(levs,1,pr=pr)
  ##levs[which(pr==max(pr))[1]]
}
#head(.logistic(-1*m1$coeff*x + m1$zeta[1]))
## this currently returns probablity less than or equal to
liability1<-function(yp,const){
  ab = as.matrix(data.frame(lapply(const,function(v){
    .logistic(v-yp)  ## shouldnt this be plus??
  })) )
  dimnames(ab)[[2]] = names(const)
  if(TRUE) return(ab)
  ac = t(apply(ab,1,function(x)diff(c(0, x, 1))))
  ac
}
liability<-function(xM){
  levs = attr(xM,"levs")
  alias = apply(xM,1,.cntNA)==0
  M = xM[alias,,drop=F]
  len =dim(M)[1]
  M1 = cbind(rep(1,len),exp(M))
  dimnames(M1)[[2]] = levs
  su = apply(M1,1,sum)
  pr = t(apply(M1,1,.prob))
  pr
  # res = rep(NA,length(alias))
  #  res[alias] = apply(pr,1,.samp,levs)
  #  attr(xM,"levs") = levs
  #  res
}



.calcYpred_multinom<-function(prev_kj, data, ind_1,levs,
                              numvar=NULL, DRS_thresh = 1e-5,kk=1,
                              liab=T,
                              vars1 = prev_kj$var,
                              betas =  prev_kj$betas,
                              # levs = names(prev_kj$tbls[[kk]]),
                              constants= prev_kj$constants_proj[[kk]]){
  #  if(getOption("fspls.DRS",F)){
  #    betas1 = apply(betas1,c(1,2), function(b) if(b>DRS_thresh) 1 else if(b<-1*DRS_thresh) -1 else 0)
  #  }
  
  if(length(vars1)==0 ){ 
    yp = .rep(rep(1/length(levs), length(levs)), length(which(ind_1)))
    return (yp)
  }
  yp = .rep(constants, length(which(ind_1)))
  
  if(!is.null(numvar) && numvar==0){
    attr(yp, "levs") =levs
    return(liability(yp))  
  }
  
  # betas1 = t(as.matrix(data.frame(betas)))
  #print(names(prev_))
  # print(prev_kj$constants_proj)
  if(!is.null(numvar))vars1 = vars1[1:numvar]
  #  const = prev_kj$const[[i]][[ycol]]
  df1 = as.matrix(data.frame(vars1))
  
  
  if(nrow(df1)>0){
    for(ki in 1:length(data)){
      inds_11 = which(df1[1,]==ki)
      if(length(inds_11)>0){
        ##CHECK IF CORRECT
        yp = yp + data[[ki]][ind_1,df1[2,inds_11], drop=F]%*% betas[[kk]][inds_11,,drop=F]
      }
    }
    #  yp = yp+ unlist(prev_kj$constants_proj)
  }
  # print(yp)
  attr(yp, "levs") =levs
  if(!liab) return(yp)
  liability(yp)
}



.calcYpred_ord<-function(prev_kj, data, ind_1,levs,
                         numvar=NULL, kk=1,
                         liab=T,
                         constants= prev_kj$constants_proj[[kk]]){
  #  if(getOption("fspls.DRS",F)){
  #    betas1 = apply(betas1,c(1,2), function(b) if(b>DRS_thresh) 1 else if(b<-1*DRS_thresh) -1 else 0)
  #  }
  vars1 = prev_kj$var
  #  yp = .rep(zeros, length(which(ind_1)))
  yp = .rep(0, length(which(ind_1)))
  
  if(length(vars1)==0 || !is.null(numvar) && numvar==0){
    attr(yp, "levs") =levs
    return(if(liab) liability1(yp, constants) else yp)
  }
  
  betas =  prev_kj$betas
  betas1 = betas[[kk]]
  #print(names(prev_))
  # print(prev_kj$constants_proj)
  if(!is.null(numvar))vars1 = vars1[1:numvar]
  #  const = prev_kj$const[[i]][[ycol]]
  df1 = as.matrix(data.frame(vars1))
  
  
  if(nrow(df1)>0){
    for(ki in 1:length(data)){
      inds_11 = which(df1[1,]==ki)
      if(length(inds_11)>0){
        ##CHECK IF CORRECT
        yp = yp + data[[ki]][ind_1,df1[2,inds_11], drop=F]%*% betas1[inds_11,,drop=F]
      }
    }
    #  yp = yp+ unlist(prev_kj$constants_proj)
  }
  # print(yp)
  attr(yp, "levs") =levs
  if(typeof(yp)=="S4") yp = as.matrix(yp)
  if(!liab) return(yp)
  liability1(yp, constants)
}




.calcYpred_binomial<-function(prev_kj, data, ind_1,
                              numvar=NULL, kk=1,
                              liab=T,
                              constants= prev_kj$constants_proj[[kk]]){
  #  if(getOption("fspls.DRS",F)){
  #    betas1 = apply(betas1,c(1,2), function(b) if(b>DRS_thresh) 1 else if(b<-1*DRS_thresh) -1 else 0)
  #  }
  vars1 = prev_kj$var
  #  yp = .rep(zeros, length(which(ind_1)))
  yp = .rep(0, length(which(ind_1)))
  
  if(length(vars1)==0 || !is.null(numvar) && numvar==0){
    #print(dim(yp)); print(constants)
    return(if(liab) .logistic(yp+constants[[1]]) else yp)
  }
  
  betas =  prev_kj$betas
  betas1 = betas[[kk]]
  #print(names(prev_))
  # print(prev_kj$constants_proj)
  if(!is.null(numvar))vars1 = vars1[1:numvar]
  #  const = prev_kj$const[[i]][[ycol]]
  df1 = as.matrix(data.frame(vars1))
  
  
  if(nrow(df1)>0){
    for(ki in 1:length(data)){
      inds_11 = which(df1[1,]==ki)
      if(length(inds_11)>0){
        ##CHECK IF CORRECT
        yp = yp + data[[ki]][ind_1,df1[2,inds_11], drop=F]%*% betas1[inds_11,,drop=F]
      }
    }
    yp=as.matrix(yp)
    #  yp = yp+ unlist(prev_kj$constants_proj)
  }
  # print(yp)
  if(!liab) return(yp)
  .logistic(yp+constants[[1]])
}

.calcYpred_1<-function(prev_kj, data, ind_1,
                       kk=1,
                       constants= prev_kj$constants_proj[[kk]]
                       #constants= unlist(prev_kj$constants_proj)
){
 
  yp = .rep(constants, length(which(ind_1)))
  betas =  prev_kj$betas
  vars1 = prev_kj$var
  if(length(vars1)==0) return(yp)
  betas1 = betas[[kk]]
  
  vars1 = prev_kj$var
  #print(names(prev_))
  # print(prev_kj$constants_proj)
  #  const = prev_kj$const[[i]][[ycol]]
  df1 = as.matrix(data.frame(vars1))
  
  
  
  
  if(nrow(df1)>0){
    for(kj in 1:length(data)){
      inds_11 = which(df1[1,]==kj)
      if(length(inds_11)>0){
        ##CHECK IF CORRECT
       # print(dim(data[[kk]]))
      #  print(length(which(ind_1)))
      #  print(df1[2,inds_11])
      #  d1 = (data[[kj]][,,drop=F])
        prod=data[[kj]][ind_1,df1[2,inds_11], drop=F]%*% betas1[inds_11,,drop=F]
      #  print(dim(prod))
      #  print(dim(yp))
        yp = yp + prod
      }
    }
    #  yp = yp+ unlist(prev_kj$constants_proj)
  }
  as.matrix(yp) 
}


rmse_interval <- function(rmse, deg_free, p_lower = 0.025, p_upper = 0.975){
  unlist(list(low = sqrt(deg_free / qchisq(p_upper, df = deg_free)) * rmse,
         mid = rmse,
         high = sqrt(deg_free / qchisq(p_lower, df = deg_free)) * rmse))
}

calcRMS<-function( predy,yTs, family , CI = F, rmsea=T, rel=F){
  if(length(yTs)==0) return(c(NA,NA,NA))
  rms = NA
  alias = which(!(is.na(predy) | is.na(yTs)))
  y = predy
  if(length(y)!=length(yTs)) stop(paste("error",length(y),length(yTs)));
  
  rms= sqrt(sum((yTs[alias] - y[alias])^2)/length(yTs[alias]))
  if(rel){
   scale =  sqrt(sum((yTs[alias] - mean(yTs[alias]))^2)/length(yTs[alias]))
  }else{
    scale = 1
  }
  rmse_interval(rms, length(predy))/scale
}
.areaBetween<-function(yp,y1,family="binomial"){
  #if(family=="ordinal"){
  #  y1[y1>=1]=1
  #}
  ap1=try(getAreaPlot(yp, y1))
  if(inherits(ap1,"try-error")){
    return(rep(NA,3))
  }
  #print(attr(ap1,"area"))
  return(c(NA,attr(ap1,"area")[[1]],NA))
}

  
.areaBetweenOld<-function(yp, y,family="binomial"){
  if(family=="binomial"){
    #    print("logistic")
    m = glm(y~1,offset=yp,family=family)
    if(abs(m$coefficients[1])>1e10){
      m = glm(y~yp,family=family)
    }
    yp2 = .logistic(yp+m$coefficients[1])
  }else{
    yp2 = yp
  }
  ty= table(y)  ## assumes ty is ordered, which it should be
  nmey = as.numeric(names(ty))
  range = c(nmey[1], nmey[length(nmey)])
  names(nmey) = nmey
  cdfs=lapply(nmey,function(yv){
    ecdf(yp2[y==yv])
  })
  kns = lapply(cdfs,  stats::knots)
  kn = unique(sort(unlist(kns)))
  matr0 = cbind(c(range[1],kn),c(kn,range[2]))
  midp = apply(matr0,1,mean)
  diffp = apply(matr0,1,function(v)v[2]-v[1])
  matr = as.matrix(data.frame(lapply(cdfs, function(cdf) cdf(midp))))
  dimnames(matr)[[2]] = names(cdfs)
  dimnames(matr)[[1]] = c(range[1], kn)
  #print(matr)
  diff = apply(matr,1,function(v) .diffs_all(v, ty))  ## should this be weighted by num samples? 
  return(c(NA,sum( diffp * diff),NA))
}
.diffs_all<-function(v,ty){
  m2=cbind(ty[-length(ty)], ty[-1])
  ty2 =apply(m2,1,function(v)v[2]+v[1])
  ty3 = ty2/sum(ty2)
  m1 = cbind(v[-length(v)], v[-1], ty3)
  sum(apply(m1,1,function(v2) (v2[1] - v2[2])*v2[3]))
}
.calcAUCW<-function(ypred,y, w,
                    conf.level=getOption("conf.level",0.95)
                    ){
  if(TRUE){
    if(length(which(y==1))==0 || length(which(y==0))==0 ) return(c(NA,0.5,NA))
    ##NOTE THE WEIGHTED VERSION SEEMS NOT TO WORK WELL
  #  print(y)
  #  print(ypred)
    roc1=try(roc(y,ypred, quiet=T))
    if(inherits(roc1,"try-error")){
      return(rep(NA,3))
    }
    #print(ci(roc1)[1:3])
   cir = ci(roc1, conf.level=conf.level)[1:3]
   if(is.na(cir[2])) cir[2] = roc1$auc
   return(cir)
  }
  
  tw = table(w)
  nw = as.numeric(names(tw))
  tw = tw * nw
  tw = tw/sum(tw)
  #   print("H")
  #  print(w)
  aucs=data.frame(lapply(nw, function(w1) .calcAUCW1(y,ypred, w==w1,q=q)))
  apply(aucs,1,function(auc) sum(auc*tw))
}
.logistic<-function(y ) 1.0 / ( 1.0 + exp(-y))
.slug<-function(x){
 gsub("/",".", gsub(" ",".",gsub("^X","",x)))
}


.misclass<-function(y_,y1_, weights_, auc=T){
  y = as.numeric(y_)
  y1 = as.numeric(y1_)
  if(length(weights_)!=length(y)) stop("problem")
  nonNA =       !is.na(y1) & !is.na(y)
  default = if(auc) NA else list()
  if(length(which(nonNA))==0) return (default)
  y = y[nonNA]
  y1 = y1[nonNA]
  weights = weights_[nonNA]
  o= order(y)
  leny = length(y)
  m =  cbind(y1,o,1:length(y), weights)
  y2 = m[o,]
  # print(head(y))
  indsk = 1:(leny-1)
  c1 = c(0,cumsum(y2[indsk,1]*y2[indsk,4]))  #fourth column is weights
  c2 = rev(cumsum(rev(y2[,4]*(1-y2[,1]))))
  # cbind(c1,c2,y2)  
  #  cbind(c1,c2)
  error =  c1+c2
  # print(error)
  minerr = min(error)
  ind = which(error == minerr)
  thresh =  y[ m[ind,3]] + 0.0000001
  # minerr_check= sum(y1[which(y<=thresh)]) + sum(1-y1[which(y>thresh)])
  corr = length(y) - minerr
  if(auc) return (corr/length(y))
  c(minerr,corr, minerr/length(y), corr/length(y), thresh)
}


##flip is doing NAs nont nonNAs


.scoreInternal=function(yp,y1, w1, type_i, fam,thresh){
  minlength = 1;
  if(type_i %in% c("correlation","rank_correlation","area","area_full","AUC","AUC_full","AUC_all")) minlength=2
  if(length(y1)<minlength){
   # print('h')
    return(unlist(list(low=NA, mid=NA, high=NA)))
  }
  if(type_i=="misclass"){
    if(fam=="gaussian") y1=as.numeric(factor(y1, sort(unique(y1))))-1
    yp=if(fam=="binomial") plogis(yp[,1]) else yp[,1]
    rms = -1*.misclass(yp,y1, w1,auc=T)
    names(rms) =c("low","mid","high")
  }else if(type_i=="area"){
    rms = 1*.areaBetween(yp[,1], y1, family=fam)[2]
    names(rms)=names(y)[[ycol]]
  }else if(type_i=="area_full"){
    rms = 1*.areaBetween(yp[,1], y1, family=fam)
    # names(rms)=names(y)[[ycol]]
    names(rms)=c("low","mid","high")
  }else if(type_i=="AUC"){
    rms = 1*(.calcAUCW(yp[,1],y1, w1))
    names(rms)=c("low","mid","high")
  }else if(type_i=="AUC_full"){
    rms = 1*(.calcAUCW(yp[,1],y1, w1))
    names(rms)=c("low","mid","high")
  }else if(type_i=="AUC_all"){
    if(fam=="ordinal"){
      levs = min(y1, na.rm=T):(max(y1,na.rm=T)-1) 
      lev_inds = 1:length(levs)
      names(lev_inds)=unlist(lapply(levs, function(l)paste(l,l+1,sep="|")))
      rms_l = .merge1_new(lapply(lev_inds, function(kj){
        gt =y1>levs[kj]
        if(length(which(gt))==0) return(NA)
        y2 = ifelse(gt,1,0)
        data.frame(value= 1*(.calcAUCW(yp[,kj],y2, w1)[2]), submeasure=c("low","mid","high"))
      }),addName="subpheno")
      #names(rms_l )=paste(names(y)[[ycol]],levs,sep=".")
      rms = rms_l #+0.5
      #            rms = c(rms_l, sum(rms_l, na.rm=T))
      #            names(rms) = c(levs,"sum")
    }else{
      #.calcAUCW(ypred, y, w)[2] for weighted AUC
      levs = unique(as.character(y1))
      mi4 = match(.slug(levs), .slug(dimnames(yp)[[2]]))
      levs = levs[!is.na(mi4)]
      mi4 = match(.slug(levs), .slug(dimnames(yp)[[2]]))
      lev_inds = 1:length(levs)
      names(lev_inds)=levs
      rms_l = .merge1_new(lapply(lev_inds, function(kj){
        y2 = ifelse(y1==levs[kj],1,0)
        data.frame(list(value=1*(.calcAUCW(yp[,mi4[kj]],y2, w1))),submeasure=c("low","mid","high"))
      }),addName="subpheno")
      rms = rms_l
    }
  }else if(type_i=="correlation" || type_i=="rank_correlation"){
    method = if(type_i=="correlation") "pearson" else "spearman"
    if(var(yp[,1])==0)warning("no variance in y")
    rms =  if(nrow(yp)==0 || var(yp[,1])==0) c(NA,NA,NA) else  .corCI(yp[,1], y1,method=method)#weights=w1[nonNA]
    names(rms)=c("low","mid","high")
  }else if(type_i=="rms"){
    rms=calcRMS(yp[,1],y1)
  }else if(type_i=="var"){
    rms=calcVar(yp[,1],y1)
  }else if(type_i=="youden"){
    rms=1*.youden(yp,y1)
  } else if(type_i=="youden_full"){
    rms=1*.youden(yp,y1)
  }else if(type_i=="youden_sens"){
    rms=1*.youden(yp,y1,typ="sens")
  }else if(type_i=="youden_spec"){
    rms=1*.youden(yp,y1,typ="spec")
  }else if(type_i=="youden_sens_full"){
    rms=1*.youden(yp,y1,typ="sens")
  }else if(type_i=="youden_spec_full"){
    rms=1*.youden(yp,y1,typ="spec")
  }else if(type_i=="sens_spec"){
    rms = .calcSensSpec(yp,y1,w1)
  }else if(type_i=="sens"){
    rms = .calcSens(yp,y1,w1, thresh=thresh)
  }else if(type_i=="npv"){
    rms = .calcNPV(yp,y1,w1, thresh=thresh)
  }else if(type_i=="f1"){
    rms = .calcF1(yp,y1,w1, thresh=thresh)
  }else if(type_i=="ppv"){
    rms = .calcPPV(yp,y1,w1, thresh=thresh)
  }else if(type_i=="spec"){
    rms = .calcSpec(yp,y1,w1, thresh=thresh)
  }else{
    stop(paste("type not recognised:",type_i))
  }
  rms
}
.initYpred1<-function(y, phensi){
  #inds1 = 1:maxn
  #names(inds1)=inds1
  nmey =names(phensi)
  names(nmey) = nmey
  family = unlist(lapply(nmey, function(x) strsplit(x,"\\.")[[1]][1]))
  
  #lapply(inds1, function(x)  {
    res = lapply(nmey, function(i){
      yi = y[[i]]
      if(family[[i]]=="multinomial"){
        yi = attr(yi,"factor")
         levs1 = levels(yi)
         names(levs1) = levs1
         rr1 = unlist(lapply(levs1, function(l1){
           vv = ifelse(yi ==l1, 1, 0)
           mean(vv,na.rm=T)
#           rep(mean(vv, na.rm=T), length(vv))
         }))
         rr = do.call(rbind, replicate(length(yi),rr1, simplify=FALSE))
         dimnames(rr) = list(dimnames(y)[[1]],levs1)
       }else if(family[[i]]=="ordinal"){
           levs1 = min(yi,na.rm=T):max(yi,na.rm=T)
           names(levs1) = levs1
          rr0= unlist(lapply(levs1[-length(levs1)], function(l1){
             vv = ifelse(yi ==l1, 1, 0)
             mean(vv,na.rm=T)
#             rep(mean(vv, na.rm=T), length(vv))
           }))
          rr1 = cumsum(rr0)
          rr = do.call(rbind, replicate(length(yi),rr1, simplify=FALSE))
          dimnames(rr) = list(dimnames(y)[[1]],levs1[-length(levs1)])
    }else{
      subinds = phensi[[i]]
      dimn = list(rownames(y[[i]]), colnames(y[[i]])[subinds])
      rr = Matrix(0,nrow = length(dimn[[1]]), ncol = length(dimn[[2]]),dimnames = dimn ,sparse=T)
        #as.matrix(data.frame(rep(mean(yi,na.rm=T),length(yi))))
    }
    rr
    })
    res
  #})
}

ypredObj<-R6Class("ypredObj", public = list(
  wname="char",
  ypreds="list",
  weights="list",
  family="character",
  #rms_prev="numeric",
  #rmsv="vector",
  rmsv_="matrix",
  types_="list",
  phensi="numeric",
  nrow="numeric",
  
  initialize=function(data, phensi,  family = data$family,
                     
                      types_=getOption("fspls.types", 
                                     default_types)){
   # rms_prev = 99999
  wname=NULL
  rmsv=NULL
  self$phensi = phensi
  self$family=family
  self$types_ = types_
    ypreds=.initYpred1(data$y, phensi)
 # self$params=params
  self$wname=wname
  self$ypreds=ypreds
  self$weights=data$weights
  self$nrow = data$nrow
  #self$rmsv_ = rmsv_
 # self$rms_prev = rms_prev
#  self$rmsv=rmsv
#  list(params=params,wname=wname,ypreds = ypreds, weights = weights,rms_prev=rms_prev, rmsv=rmsv)
  },
#updateYpredsInds=function(data,prev,flip nonNA  ){ #= self$train$looc_incl[,k2]
#  ypred = self
#  phensi = self$phensi
#  within=(k2==self$nreps())
#  ypred$updateYP(data, prev, nonNA, !within)  
#},

updateYP=function(data,prev,  nonNA, flip=T, ignore.na=F){
  ypred = self
  prev_kj = prev 
  phensi=self$phensi
  prev_kj$var = lapply(prev_kj$var_names, data$convert)  ## this would not be threadsafe
  vars_to_incl = which(unlist(lapply(prev_kj$var, length))==2)
  prev_kj$var = prev_kj$var[vars_to_incl]
  prev_kj$varnames = prev_kj$varbanes[vars_to_incl]
  prev_kj$betas = lapply(prev_kj$betas, function(xx)xx[vars_to_incl,,drop=F])
  na_x = if(ignore.na) rep(F, self$nrow) else data$getNA(prev_kj$var)
  for(kk1 in names(phensi)){ #} 1:length( ypred$ypreds)){
    kk = phensi[[kk1]]
    if(is.null(nonNA)){
      ind_1 = if(flip) rep(T, self$nrow) else rep(F,self$nrow)
    }else{
      ind_1 = if(flip) !nonNA else nonNA
    }
    levs1 = NULL
    family = strsplit(kk1,"\\.")[[1]][1]
    # if(family=="multinomial") levs1=dimnames(self$y[[kk1]])[[2]]
    #  if(family=="ordinal")levs1 = min(self$y[[kk1]][,kk], na.rm=T):max(self$y[[kk1]][,kk],na.rm=T) 
    self$calcYpred(prev_kj,data$data,ind_1,kk1, kk,na_x, family=family)
  }
  #   names(ypred$ypreds) = names(self$y)
  
  #}
},

#calcYpred(prev_kj,self$data,ind_1,levs,numvar,kk1, kk)
  calcYpred=function(prev_kj, data, ind_1,  kk1, kk,na_x,   family = self$family[[kk]]){  ## kk1 in model space 
  
    #      ypred$ypreds[[kk]]$calcYpred(prev_kj,self$data,ind_1,levs,numvar,kk1, self$family[[kk]])
    if(family=="multinomial"){
      levs = dimnames( self$ypreds[[kk1]])[[2]]
      ab=.calcYpred_multinom(prev_kj,  data, ind_1, levs, kk=kk1)  ## for multi-prediction
      mi22 = match(dimnames( self$ypreds[[kk1]])[[2]], levs)
      self$ypreds[[kk1]][ind_1,is.na(mi22)]=0
      self$ypreds[[kk1]][ind_1,!is.na(mi22)] =  ab[,mi22[!is.na(mi22)]] 
    }else if(family=="ordinal"){
      constants = prev_kj$constants_proj[[kk1]]
      for(kk_1 in 1:length(kk)){
        levs1 = 0:length(constants[[kk_1]])
         ab= .calcYpred_ord(prev_kj,  data, ind_1, levs = levs1,
                                                 kk=kk_1, constants = constants[[kk_1]])  ## for multi-prediction
           self$ypreds[[kk1]][ind_1,] =  ab
      }
    }else if(family=="binomial"){
      constants = prev_kj$constants_proj[[kk1]]
      for(kk_1 in 1:length(kk)){
        ab =   .calcYpred_binomial(prev_kj,  data, ind_1, 
                                                    kk=kk_1, constants = constants[kk_1])  ## for multi-prediction
        self$ypreds[[kk1]][ind_1,kk_1] =ab
      }
    }else{
     # self$ypreds[[kk]]= array(dim=c(length(ind_1),length(kk)))
      constants = prev_kj$constants_proj[[kk1]]
      for(kk_1 in 1:length(kk)){
        ab = .calcYpred_1(prev_kj,  data, ind_1,kk=kk_1, constants=constants[kk_1]) 
        self$ypreds[[kk1]][ind_1,kk_1] = ab[,1]
      }
    }
    if(length(which(na_x))>0){
      self$ypreds[[kk1]][na_x,] = rep(NA, ncol(self$ypreds[[kk1]]))
    }
    
  },
calcRMSV=function(y, nonNA,       flip=F){
  phensi = self$phensi
  ypreds = self$ypreds
  types_ =self$types_
  nme_phens = names(phensi)
  names(nme_phens) = nme_phens
  ind_1 = if(flip) !nonNA else nonNA
  w1= if(is.null(ind_1)) self$weights else self$weights[ind_1]
  nsamps = length(which(ind_1))
  .merge1_new(lapply(nme_phens, function(nme_p1){
    family= strsplit(nme_p1,"\\.")[[1]][1]
    fam = family
    mtype = match(family,names(types_))
    if(length(which(is.na(mtype)))>0) stop("could not find match in types_")
    types_i = types_[[mtype]]
    names(types_i) = types_i
    ycols = phensi[[nme_p1]]
    ypreds1 = ypreds[[nme_p1]]
    names(ycols) = names(ypreds1)
    #  print(names(ycols))
    phens = dimnames(ypreds1)[[2]]
    ycol_inds = 1:length(ycols)

    y2 = if(is.null(ind_1)) y[[nme_p1]] else  y[[nme_p1]][ind_1,,drop=F]
    yp2 = if(is.null(ind_1))ypreds1 else  ypreds1[ind_1,,drop=F] 
    names(ycol_inds) = dimnames(y2)[[2]][ycols]
    rms_1=lapply(ycol_inds, function(ycol_ind){
      ycol = ycols[ycol_ind]
      y1 =  y2[,ycol]
      yp =if(fam=="ordinal") yp2 else yp2[,ycol_ind,drop=F] 
      nonNA = !is.na(y1)
      nonNA = nonNA & !is.na(yp[,1]) 
      .merge1_new(
        lapply(types_i, function(type_i1){
          type_i1s = strsplit(type_i1,"\\.")[[1]]
          type_i= type_i1s[[1]]
          thresh = if(length(type_i1s)>1) type_i1s[2] else NA
          rms =.scoreInternal(yp[nonNA,,drop=F], y1[nonNA],w1[nonNA], type_i, fam,thresh)
          if(is.null(rms)) stop(paste(type_i,"rms NULL"))
          if(is.null(names(rms))) names(rms) = 1:length(rms)
          if(typeof(rms)=="list"){
            return(rms)
          }else{
            df12=data.frame(cbind(names(rms), rms)) #, type_i) #as.matrix(rms)
            names(df12)=c("submeasure","value")#,"measure")
            df12[,2]=as.numeric(df12[,2])
          }
          df12
        }), num_cols="value", addName="measure")%>% tibble::add_column(nsamps = nsamps, cv=flip)
    })
    rms_2 = .fixBeforeMerge(rms_1)
    df3=.merge1_new(rms_2,num_cols="value", addName="pheno")
    if(!("subpheno" %in% names(df3))){
      df3 = df3 %>% tibble::add_column(subpheno="")
    }
    df3    
  }),addName="family")
  #.merge1_new(rms_3, num_cols = "value", addName="beam")
  # dimnames(res_df)[[1]] = names(rms_3)
  # as.matrix(res_df)
}

)
)





###funcs
.youden<-function(yp, y1,w, typ=NULL){
  #  print(yp)
  threshv = sort(unique(yp[,1]))
  names(threshv) = threshv
  df1 = unlist(lapply(threshv, function(t){
     .calcSens(yp, y1,w,t)[2]+.calcSpec(yp, y1,w,t)[2]-1
  }))
  t = threshv[which.max(df1)][1]
  if(is.null(typ)){
   return( .calcSens(yp, y1,w,t)+.calcSpec(yp, y1,w,t)-1)
  }
  if(typ=="sens") return( .calcSens(yp, y1,w,t))
  .calcSpec(yp, y1,w,t)
}

.calcSensSpec<-function(yp, y1,w){
#  print(yp)
 threshv = sort(unique(yp[,1]))
 names(threshv) = threshv
 df1 = t(data.frame(lapply(threshv, function(t){
   c(.calcSens(yp, y1,w,t), .calcSpec(yp, y1,w,t))
 })))
 dimnames(df1) = list(threshv, c("sens.low", "sens.mid","sens.high", "spec.low","spec.mid","spec.high"))
 i=which.max(apply(df1, 1, function(v) v[2] +v[5]))
res1 = df1[i,]
names(res1)  = paste(names(res1),signif(threshv[i], digits =2),sep=':')
  res1
}
.calcMissing<-function(yp, lower, upper){
  M = length(which(yp>=lower & yp<=upper))
  T = length(yp)
  a=binom.confint(M,T, method="probit",conf.level=getOption("conf.level",0.95))
  res = c(a$lower, a$mean, a$upper) 
  names(res)=c("low","mid","high")
  res
}
.calcSens<-function(yp, y1,w, thresh=0.5){
  P =length( which( y1==1))
  N =length( which( y1==0))
  TP = length(which(y1[yp>=thresh]==1))
  
  #  TN = spec*N
  #  o = order(yp, decreasing=F)
  #  cum=cumsum(1-y1[o])
  #  pos=which(cum>=TN)[1]
  
  #  TN1 = length(which(y1[o][1:pos]==0))
  
  #  TP=length(which(y1[o][(pos+1):length(y1)]==1))
  
  if(P==0) return (c(NA,NA,NA))
  a=binom.confint(TP,P, method="prop.test",conf.level=getOption("conf.level",0.95))
  res=c(a$lower, a$mean, a$upper)
  names(res)=c("low","mid","high")
  res
}
.calcPrev<-function(yp,y1,w,thresh =0.5){
  P =length( which( y1==1))
  N =length( which( y1==0))
  if(P+N==0) return (c(NA,NA,NA))
  a=binom.confint(P,P+N, method="prop.test", conf.level=getOption("conf.level",0.95))
  
  res=c(a$lower, a$mean, a$upper)
  names(res)=c("low","mid","high")
  res
}

.calcF1<-function(yp, y1,w, thresh=0.5){
  TP = length(which(y1[yp>=thresh]==1))
  FP = length(which(y1[yp>=thresh]==0))
  FN = length(which(y1[yp<thresh]==1))
  TN = length(which(y1[yp<thresh]==0))
  if(2*TP+FP+FN==0) return(rep(NA,3))
  a=binom.confint(2*TP,2*TP+FP+FN, method="probit",conf.level=getOption("conf.level",0.95))
  res=c(a$lower, a$mean, a$upper)
  names(res)=c("low","mid","high")
  res
}


.calcNPV<-function(yp, y1,w, thresh=0.5){
  FN = length(which(y1[yp<thresh]==1))
  TN = length(which(y1[yp<thresh]==0))
  if(TN+FN==0) return(c(NA,NA,NA))
  #print(paste(TN, FN))
  a=binom.confint(TN,TN+FN, method="probit",conf.level=getOption("conf.level",0.95))
  res=c(a$lower, a$mean, a$upper)
  names(res)=c("low","mid","high")
  res
}

.calcPPV<-function(yp, y1,w, thresh=0.5){
  TP = length(which(y1[yp>=thresh]==1))
  FP = length(which(y1[yp>=thresh]==0))
  if(TP+FP==0) return(c(NA,NA,NA))
  a=binom.confint(TP,TP+FP, method="probit",conf.level=getOption("conf.level",0.95))
  res=c(a$lower, a$mean, a$upper)
  names(res)=c("low","mid","high")
  res
}
.calcSpec<-function(yp, y1,w, thresh =0.5){
  P =length( which( y1==1))
  N =length( which( y1==0))
  TN = length(which(y1[yp<thresh]==0))
  #TP = sens*P
  
  #o = order(yp, decreasing=T)
  #cum=cumsum(y1[o])
  #pos=which(cum>=TP)[1]
  
  #TN = length(which(y1[o][(pos+1):length(y1)]==0))
  #TP1 = length(which(y1[o][1:pos]==1))
  if(N==0) return (c(NA,NA,NA))
  a=binom.confint(TN,N, method="prop.test",conf.level=getOption("conf.level",0.95))
  
  res=c(a$lower, a$mean, a$upper)
  names(res)=c("low","mid","high")
  res
}


.calcPPV<-function(yp, y1,w, thresh =0.5){ #TP/TP+FP
  TP = length(which(y1[yp>=thresh]==1))
  FP = length(which(y1[yp>=thresh]==0))
  if(TP+FP==0) return( rep(NA,3))
  a=binom.confint(TP,TP+FP, method="prop.test",conf.level=getOption("conf.level",0.95))
  res=c(a$lower, a$mean, a$upper)
  names(res)=c("low","mid","high")
  res
}

.calcNPV<-function(yp, y1,w, thresh =0.5){ #TN/TN+FN
  TN = length(which(y1[yp<thresh]==0))
  FN = length(which(y1[yp<thresh]==1))
  if(TN+FN==0) return (rep(NA,3))
  a=binom.confint(TN,TN+FN, method="prop.test",conf.level=getOption("conf.level",0.95))
 res= c(a$lower, a$mean, a$upper)
 names(res)=c("low","mid","high")
 res
}

