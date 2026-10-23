#.mergeRMSV<-function(datas){
#  rmsv_=.merge1(lapply(datas, function(d)d$rmsv_),num_cols="value",addName="data")
#  rmsv1 = .findMinRMSV(rmsv_, params$mult, fspls.sum=getOption("fspls.sum",T))
#  min(rmsv1)
#}
default_rms = unlist(list(low=NA, mid=NA, high=NA))

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
  ind = which(pr==max(pr,na.rm=TRUE))
  return(levs[ind[1]])
  #sample(levs,1,pr=pr)
  ##levs[which(pr==max(pr))[1]]
}
#head(.logistic(-1*m1$coeff*x + m1$zeta[1]))
## this currently returns probablity less than or equal to
liability1<-function(yp,const, liab=TRUE){
  ab = as.matrix(data.frame(lapply(const,function(v){
    if(liab) .logistic(yp-v) else yp-v  ## according to description of polr in R help
   # if(liab) .logistic(v-yp) else v-yp  ## according to description of polr in R help
  })) )
  dimnames(ab)[[2]] = names(const)
  if(TRUE) return(ab)
  ac = t(apply(ab,1,function(x)diff(c(0, x, 1))))
  ac
}
liability_old<-function(xM){
  levs = attr(xM,"levs")
  alias = apply(xM,1,.cntNA)==0
  M = xM[alias,,drop=FALSE]
  len =dim(M)[1]
  M1 = cbind(rep(1,len),exp(M))
  dimnames(M1)[[2]] = levs
  su = apply(M1,1,sum)
  pr = t(apply(M1,1,.prob))
  pr
}
liability<-function(xM){ ## use with glmnet output
  levs =dimnames(xM)[[2]]
  alias = apply(xM,1,.cntNA)==0
  M = xM[alias,,drop=FALSE]
  len =dim(M)[1]
  M = M/max(M,na.rm=TRUE)
  M1 = exp(M)
  dimnames(M1)[[2]] = levs
  su = apply(M1,1,sum)
  pr = t(apply(M1,1,.prob))
  pr
}


.calcYpred_multinom<-function(prev_kj, data, ind_1,levs,
                              kk=1,
                              liab=TRUE,
                              vars1 = prev_kj$var,
                              betas =  prev_kj$betas,
                              transforms = NULL,
                              # levs = names(prev_kj$tbls[[kk]]),
                              constants= prev_kj$constants_proj[[kk]]){
  #  if(getOption("fspls.DRS",F)){
  #    betas1 = apply(betas1,c(1,2), function(b) if(b>DRS_thresh) 1 else if(b<-1*DRS_thresh) -1 else 0)
  #  }
  yp = .rep(constants, length(which(ind_1)))
  
  if(length(vars1)==0 ){ 
    return( liability(yp))
#    yp = .rep(rep(1/length(levs), length(levs)), length(which(ind_1)))
#    return (yp)
  }
  beta1 = betas[[kk]]
  incl1 = unlist(lapply(vars1, length))>1
  df1 = as.matrix(data.frame(vars1[incl1]))
  beta1 = beta1[incl1,,drop=FALSE]
  
  
  if(nrow(df1)>0){
    for(ki in 1:length(data)){
      inds_11 = which(df1[1,]==ki)
      if(length(inds_11)>0){
        ##CHECK IF CORRECT
#        d2 =data[[ki]][ind_1,df1[2,inds_11], drop=FALSE]
        d2 = .extract1(data[[ki]],ind_1, df1, inds_11 ,transforms)
        #print(beta1)
        extra=d2%*%  beta1[inds_11,,drop=FALSE]
        yp = yp + extra
      }
    }
    #  yp = yp+ unlist(prev_kj$constants_proj)
  }
  # print(yp)
  levs = names(prev_kj$tbls[[kk]])
  attr(yp, "levs") =levs
  #print(yp)
  if(!liab) return(yp)
  liability(yp)
}



.calcYpred_ord<-function(prev_kj, data, ind_1,levs,
                         kk=1,
                         liab=TRUE,
                         transforms = NULL,
                         betas1=prev_kj$betas[[kk]],
                         const =prev_kj$constants_proj[[kk]][[1]]){
  #  if(getOption("fspls.DRS",F)){
  #    betas1 = apply(betas1,c(1,2), function(b) if(b>DRS_thresh) 1 else if(b<-1*DRS_thresh) -1 else 0)
  #  }
  vars1 = prev_kj$var
  #  yp = .rep(zeros, length(which(ind_1)))
  yp = .rep(0, length(which(ind_1)))
  
  if(length(vars1)==0){
    attr(yp, "levs") =levs
    return(liability1(yp, const, liab=liab))
#    return(if(liab) liability1(yp, constants) else yp)
  }
  
#  betas =  prev_kj$betas
#  betas1 = betas[[kk]]
  #print(names(prev_))
  # print(prev_kj$constants_proj)
  #if(!is.null(numvar))vars1 = vars1[1:numvar]
  #  const = prev_kj$const[[i]][[ycol]]
  incl1 = unlist(lapply(vars1, length))>1
  df1 = as.matrix(data.frame(vars1[incl1]))
  betas1 = betas1[incl1,,drop=FALSE]
  
  if(nrow(df1)>0){
    for(ki in 1:length(data)){
      inds_11 = which(df1[1,]==ki)
      if(length(inds_11)>0){
        ##CHECK IF CORRECT
        d2 = .extract1(data[[ki]],ind_1, df1, inds_11 ,transforms)
        yp = yp + d2%*% betas1[inds_11,,drop=FALSE] # d2 = .extract1(data1[[kj]],ind_1, df1, inds_11 ,transforms)
      }
    }
    #  yp = yp+ unlist(prev_kj$constants_proj)
  }
  # print(yp)
  attr(yp, "levs") =levs
  if(typeof(yp)=="S4") yp = as.matrix(yp)
  #if(!liab) return(yp)
  liability1(yp, const, liab=liab)
}




.calcYpred_binomial<-function(prev_kj, data, ind_1,
                               kk=1,kk1="",
                              liab=TRUE,
                              transforms = NULL,
                              betas1 = prev_kj$betas[[kk1]],
                              constants= prev_kj$constants_proj[[kk1]]){
  #  if(getOption("fspls.DRS",F)){
  #    betas1 = apply(betas1,c(1,2), function(b) if(b>DRS_thresh) 1 else if(b<-1*DRS_thresh) -1 else 0)
  #  }
  vars1 = prev_kj$var
  #  yp = .rep(zeros, length(which(ind_1)))
#  yp = .rep(constants, length(which(ind_1)))
  yp = .rep(constants, length(which(ind_1)))
  
  if(length(vars1)==0 ){
    #print(dim(yp)); print(constants)
    return(if(liab) .logistic(yp+constants[[1]]) else yp)
  }
  
  #betas =  prev_kj$betas
  #betas1 = betas[[kk1]]
  
 # mean_x = prev_kj$mean_x
  #print(mean_x)
  #print(betas1)
  #y1_off = mean_x %*% betas1
  #print(y1_off)
  
#  yp = .rep(constants+y1_off[1,1], length(which(ind_1)))
  
  
  
  incl1 = unlist(lapply(vars1, length))>1
  df1 = as.matrix(data.frame(vars1[incl1]))
  betas1 = betas1[incl1,,drop=FALSE]
  
  #print(names(prev_))
  # print(prev_kj$constants_proj)
  #  const = prev_kj$const[[i]][[ycol]]
  df1 = as.matrix(data.frame(vars1))
  
  
  if(nrow(df1)>0){
    for(ki in 1:length(data)){
      inds_11 = which(df1[1,]==ki)
      if(length(inds_11)>0){
        ##CHECK IF CORRECT
        #d2 =data[[ki]][ind_1,df1[2,inds_11], drop=FALSE]
        d2 = .extract1(data[[ki]],ind_1, df1, inds_11 ,transforms)
        #d2 =t( t(data[[ki]][ind_1,df1[2,inds_11], drop=FALSE]) - prev_kj$mean_x[inds_11])
        yp = yp + d2%*% betas1[inds_11,,drop=FALSE]
      }
    }
    yp=as.matrix(yp)
    #  yp = yp+ unlist(prev_kj$constants_proj)
  }
  # print(yp)
  if(!liab) return(yp)
  .logistic(yp)
}

.extract1<-function(data2,ind_1, df1, inds_11, transforms){
  d2 = data2[ind_1,df1[2,inds_11], drop=FALSE]
  if(!is.null(transforms)){
    tr = transforms[df1[3,inds_11]]
    for(jkk in 1:ncol(d2)){
      d2[,jkk] = tr[[jkk]][[1]](d2[,jkk])
    }
  }
  d2
}



rmse_interval <- function(rmse, deg_free, p_lower = 0.025, p_upper = 0.975){
  unlist(list(low = sqrt(deg_free / qchisq(p_upper, df = deg_free)) * rmse,
         mid = rmse,
         high = sqrt(deg_free / qchisq(p_lower, df = deg_free)) * rmse))
}

calcRMS<-function( predy,yTs, family , CI =FALSE, rmsea=TRUE, rel=FALSE){
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
  ap1=try(summariseAreaPlot(getAreaPlot(yp, y1)))
  
  #ggplot(ap1, aes(x=knots, y=value, color=subpheno, size=counts))+geom_point()
  if(inherits(ap1,"try-error")){
   resu = (rep(NA,3))
  }
  #print(attr(ap1,"area"))
  resu=(c(NA,attr(ap1,"area")[[1]],NA))
  names(resu)=c("low","mid","high")
  resu
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

.logistic<-function(y ) 1.0 / ( 1.0 + exp(-y))
.slug<-function(x){
 gsub("/",".", gsub(" ",".",gsub("^X","",x)))
}


.misclass<-function(y_,y1_, weights_, auc=TRUE){
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

.calcArea<-function(yp_kj,y2, type_i = "area"){
  ap1 = getAreaPlot(yp_kj, y2)
  ap2=summariseAreaPlot(ap1)
  ap2[[type_i]]
}
##flip is doing NAs nont nonNAs
#   .scoreInternal(yp[nonNA,,drop=FALSE], y1[nonNA],w1[nonNA], type_i, fam,thresh)

.scoreInternal<-function(yp,y1, w1, type_i, fam,thresh,liab=TRUE){
  
  minlength = 1;
  if(type_i=="AUC" && fam=="ordinal") type_i = "AUC_all"
  if(type_i=="AUPRC" && fam=="ordinal") type_i = "AUPRC_all"
  if(type_i=="area" && fam=="ordinal") type_i = "area_all"
  
 # if(type_i=="area" && fam=="ordinal") type_i = "AUC_all"
  if(type_i %in% c("correlation","rank_correlation","area","area_full","AUC","AUC_full","AUC_all","var")) minlength=2
  if(length(y1)<minlength){
   # print('h')
    return(unlist(list(low=NA, mid=NA, high=NA)))
  }
 
  if(type_i %in% c("area","max_diff","max_diff_x")){
    rms=tryCatch({
      ap1 = getAreaPlot(yp, y1)
      ap2=summariseAreaPlot(ap1)
      ap2[[type_i]]
    },error=function(w){
      return(rep(NA,3))
    })
     names(rms)=c("low","mid","high")
  }else if(type_i=="misclass"){
    if(fam=="gaussian") y1=as.numeric(factor(y1, sort(unique(y1))))-1
    yp=if(fam=="binomial") plogis(yp[,1]) else yp[,1]
    rms = -1*.misclass(yp,y1, w1,auc=TRUE)
    names(rms) =c("low","mid","high")
 } else if(type_i=="AUPRC"){
   tbl = table(y1);
   tbl = sort(table(y1), decreasing=TRUE);
   if(length(tbl)<=1 || tbl[2]<=1) {
    
       rms = default_rms;
     }else{
       rms = 1*(.calcAUPRC(yp[,1],y1, w1))
       names(rms)=c("low","mid","high")
     }
   }else if(type_i=="AUC"){
    rms = 1*(.calcAUCW(yp[,1],y1, w1))
    names(rms)=c("low","mid","high")
  }else if(type_i=="AUC_full"){
    rms = 1*(.calcAUCW(yp[,1],y1, w1))
    names(rms)=c("low","mid","high")
  }else if(type_i %in% c("AUC_all","AUPRC_all","area_all")){
    if(fam=="ordinal"){
      levs = if(liab) unlist(lapply(colnames(yp), function(xi)as.numeric(strsplit(xi,"\\|")[[1]][1]))) else sort(unique(y1))[-1]
      lev_inds = 1:length(levs)
      names(lev_inds)=levs
      rms_l = .merge1_new(lapply(lev_inds, function(kj){
        gt =y1<=levs[kj]  ##interpret as probability of being less or equal to  count
        if(length(which(gt))==0) return(NULL)
        y2 = ifelse(gt,1,0)
        ypkj = if(liab)yp[,kj] else yp[,1]
       
        
        if(type_i=="area_all") return(data.frame(value= 1*(.calcArea(ypkj,y2)), submeasure=c("low","mid","high")))
        
        if(type_i=="AUC_all") return(data.frame(value= 1*(.calcAUCW(ypkj,y2, w1)), submeasure=c("low","mid","high")))
        tbl = sort(table(y2), decreasing=TRUE);
        if(length(tbl)<=1 || tbl[2]<=1) return(data.frame(list(value=default_rms), submeasure=names(default_rms)))
        
        data.frame(value= 1*(.calcAUPRC(ypkj,y2, w1)), submeasure=c("low","mid","high"))
      }),addName="subpheno")
      #names(rms_l )=paste(names(y)[[ycol]],levs,sep=".")
      rms = rms_l #+0.5
      #            rms = c(rms_l, sum(rms_l, na.rm=TRUE))
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
  }else if(type_i=="mad"){
    rms=c(NA,median(abs(yp[,1]-y1)),NA)
    names(rms)=c("low","mid","high")
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
.initYpred1<-function(y, phensi, family, liab){
  #inds1 = 1:maxn
  #names(inds1)=inds1
  nmey =names(phensi)
  names(nmey) = nmey
  indsy = 1:length(nmey); names(indsy)=nmey
  #family = unlist(lapply(nmey, function(x) strsplit(x,"\\.")[[1]][1]))
  
  #lapply(inds1, function(x)  {
    res = lapply(indsy, function(i){
      i1 = nmey[[i]]
      yi = y[[i1]]
      if(family[[i]]=="multinomial"){
        yi = attr(yi,"factor")
         levs1 = levels(yi)
         names(levs1) = levs1
         rr1 = unlist(lapply(levs1, function(l1){
           vv = ifelse(yi ==l1, 1, 0)
           mean(vv,na.rm=TRUE)
#           rep(mean(vv, na.rm=TRUE), length(vv))
         }))
         rr = do.call(rbind, replicate(length(yi),rr1, simplify=FALSE))
         dimnames(rr) = list(dimnames(y[[1]])[[1]],levs1)
       }else if(family[[i]]=="ordinal" && liab){
         yi= yi[,1]
           levs1 = levels(factor(yi, ordered=T))
           levs2 = unlist(lapply(2:length(levs1), function(i3) paste(levs1[i3-1], levs1[i3], sep="|")))
         
          rr0= unlist(lapply(levs1[-length(levs1)], function(l1){
             vv = ifelse(yi ==l1, 1, 0)
             mean(vv,na.rm=TRUE)
#             rep(mean(vv, na.rm=TRUE), length(vv))
           }))
          rr1 = 1-cumsum(rr0)
          rr = do.call(rbind, replicate(length(yi),rr1, simplify=FALSE))
          dimnames(rr) = list(dimnames(y[[1]])[[1]],levs2)
    }else{
      subinds = phensi[[i]]
      dimn = list(rownames(y[[i]]), colnames(yi)[subinds])
      rr = Matrix(0,nrow = length(dimn[[1]]), ncol = length(dimn[[2]]),dimnames = dimn ,sparse=TRUE)
     ## colnames(rr) = colnames(yi)[subinds]
        #as.matrix(data.frame(rep(mean(yi,na.rm=TRUE),length(yi))))
    }
    rr
    })
    res
  #})
}

ypredObj<-R6::R6Class("ypredObj", public = list(
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
  liab="boolean",
  
  initialize=function(data, phensi,  liab=TRUE,
                      family = unlist(lapply(names(phensi), function(str)getOption("fspls.family",strsplit(str,"\\.")[[1]][1]))),
                     
                      types_=getOption("fspls.types", 
                                     default_types)){
   # rms_prev = 99999
  wname=NULL
  rmsv=NULL
  self$liab = liab
  self$phensi = phensi
  self$family=family
  self$types_ = types_
  self$ypreds=.initYpred1(data$y, phensi, family, liab)
 # self$params=params
  self$wname=wname
  #self$ypreds=ypreds
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
## use inv_transform on y  
updateEnsemble=function(ypreds1, coeff1, nonNA, colnames,liab=FALSE,flip=TRUE){
  nmes = names(coeff1); names(nmes)=nmes
  
  if(is.null(nonNA)){
    ind_1 = if(flip) rep(TRUE, self$nrow) else rep(FALSE,self$nrow)
  }else{
    ind_1 = if(flip) !nonNA else nonNA
  }
  for(nme in nmes){
    lambda = attr(coeff1[[nme]],"lambda")
  family = strsplit(nme,"\\.")[[1]][1]
  #if(is.null(lambda) && family!="gaussian" && family!="multinomial") stop("lambda is NULL")
  
  if(family=="ordinal"){
    
    interc = attr(coeff1[[nme]],"zeta")
    coeff2 = as.matrix(coeff1[[nme]])
  }else{
    interc = unlist(coeff1[[nme]][1,,drop=TRUE])
    coeff2 = as.matrix(coeff1[[nme]][-1,,drop=F]); 
  }
#  print(nme)
    if(nrow(coeff2)==0){
      self$ypreds[[nme]][ind_1,] = matrix(interc, nrow = length(ind_1), ncol = length(interc), byrow = TRUE)
    }else if(family=="multinomial"){
      indy = 1:ncol(coeff2); names(indy) = names(coeff2)
      names(ypreds1) = colnames
       li1 = unlist( lapply(ypreds1, function(yp){
         df0=data.frame(yp$ypreds[[nme]][ind_1,])
        }), recursive=FALSE)
       mat=as.matrix(data.frame(li1))
      # match(rownames(coeff2), colnames(mat))
        # coeff1=coeff[[nme]][[i]]
       m1 =   matrix(unlist(interc), nrow = nrow(mat), ncol = ncol(coeff2), byrow = TRUE)
        
       self$ypreds[[nme]][ind_1,] =   mat %*% coeff2 +m1
   
    }else{
      if(is.null(lambda)) lambda = rep(1,ncol(coeff2))
         # mi1 = match( colnames,rownames(coeff2));
        #  mi2 = match( rownames(coeff2), colnames);
        #  if(length(which(is.na(mi2)))>0 || length(which(is.na(mi1)))>0) warning("coeffs dont match");
          
         # coeff2 = as.matrix(coeff2[mi1,,drop=F])
      #    if(length(which(is.na(mi1)))>0) stop("problem")
          indy = 1:ncol(coeff2); names(indy) = names(coeff2)
          li1 = lapply(indy, function(i){
            mat=as.matrix(data.frame(lapply(ypreds1, function(yp){
              yp$ypreds[[nme]][ind_1,i]
            })))
            if(family=="gaussian"){
              fam = gaussian(link = sinh_asinh_link(lambda[[i]]))
              return(fam$linkinv(mat %*% coeff2[,i]+interc[i]))
            }
            if(family=="ordinal"){
              if(!liab) return( mat %*% coeff2[,i])
              aop=1-ao_probs(mat,interc, coeff2[,i],  lambda[[i]], type="cum")
              return(aop); 
            }
            if(family=="binomial"){
              ab =  mat %*% coeff2[,i]+interc[i]
              if(!liab) return(ab)
              return(pao_(ab, lambda[[i]]))
            }
            
            stop("unrecognised family")
           # coeff1=coeff[[nme]][[i]]
           
          })
          result=if(length(li1)==1) li1[[1]] else as.matrix(data.frame(li1))
          #colnames(ab)
          self$ypreds[[nme]][ind_1,] = result
    }
    
  }
  preds = self$predictions(nonNA, flip=flip)
  preds
},
updateYP=function(d,full_model,  nonNA,flip=TRUE, ignore.na=FALSE){
  prev_i1=full_model
  liab = self$liab
  ypred = self
  prev_kj = prev_i1 
  phensi=self$phensi
  prev_kj$var = lapply(prev_kj$var_names, d$convert)  ## this would not be threadsafe
  na_x = if(ignore.na) rep(FALSE, self$nrow) else d$getNA(prev_kj$var)
  #kk1 = names(phensi)[[1]]
#  inv_func = NULL
  len = length(prev_kj$var)
  if(is.null(nonNA)){
    ind_1 = if(flip) rep(TRUE, self$nrow) else rep(FALSE,self$nrow)
  }else{
    ind_1 = if(flip) !nonNA else nonNA
  }
   for(kk1 in names(phensi)){ #} 1:length( ypred$ypreds)){
    kk = phensi[[kk1]]
    levs1 = NULL
    family =getOption("fspls.family",strsplit(kk1,"\\.")[[1]][1])
    self$calcYpred(prev_kj,d,ind_1,kk1, kk,na_x, family=family)
  }
  #   names(ypred$ypreds) = names(self$y)
  
  #}
},

#calcYpred(prev_kj,self$data,ind_1,levs,numvar,kk1, kk)
# if inv_func is NULL we should forward transform the x otherwise we transform y
  calcYpred=function(prev_kj, d, ind_1,  kk1, kk,na_x, 
                     family = self$family[[kk]]){  ## kk1 in model space 
    transforms = d$getTransforms(prev_kj$var )
    liab=self$liab
    x_ = d$extractData(prev_kj$var, adjust=FALSE)
    betas=prev_kj$betas_proj[[kk1]]
   
    if(length(prev_kj$var)==0){
          x_trans =matrix(0, nrow = length(which(ind_1)) ,ncol = nrow( betas))
    }else{
   
      x_trans=.eval1_(x_[ind_1,,drop=FALSE], prev_kj$Wall[[kk1]], transforms, family)
      
          #ab = x_trans %*% betas ## mean_x = prev_kj$mean_x)
          #ab = as.matrix(ab)
    }
    if(nrow(betas)!=ncol(x_trans)) betas = t(betas)
    #print(cbind(dim(x_trans),dim(betas) ))
#    print(dim(x_))
    lambdas = prev_kj$lambdas[[kk1]]
    constants = prev_kj$constants_proj[[kk1]]
    if(family=="multinomial"){
      ab =   x_trans %*% betas
      levs = names(prev_kj$tbls[[kk1]][kk])
      if(is.null(dim(ab)) ){ 
        ab = matrix(1, ncol=length(levs), nrow = length(which(ind_1)))
        dimnames(ab)[[2]] = levs
      }
      if(length(lambdas)>1){
        for(jkk in 1:ncol(ab)){
          ab[,jkk] = if(liab) pao_(ab[,jkk]+ constants[[jkk]], lambdas[[jkk]]) else ab[,jkk]+ constants[[jkk]]
        }
      #  ab = t(apply(ab,1,.prob))
        ab = ab/apply(ab,1,sum)
      }else{
     
        attr(ab, "levs") =levs
        for(jkk in 1:ncol(ab)){
          ab[,jkk] = ab[,jkk]+ constants[[jkk]]
        }
        if(max(abs(ab))==0) ab =ab+1
        #print(yp)
        if(liab) ab=liability(ab)
      }
    
      
      
      
      levs = dimnames( self$ypreds[[kk1]])[[2]]
      mi22 = match(dimnames( self$ypreds[[kk1]])[[2]], dimnames(ab)[[2]])
      self$ypreds[[kk1]][ind_1,is.na(mi22)]=0
      self$ypreds[[kk1]][ind_1,!is.na(mi22)] =  as.matrix(ab[,mi22[!is.na(mi22)]]) 
    }else if(family=="ordinal"){
     
      ncols = ncol(self$ypreds[[kk1]])
      constk = constants[[1]]
      
      if(length(constk)==0 || (length(constk!=ncols) && max(abs(constk))<1e-9)){
        constk = rep(0, ncols); names(constk) = colnames(self$ypreds[[kk1]])      
      }
      
     
      
      
      
      if(liab){
        #mi1 = match(names(constk), colnames(self$ypreds[[kk1]]))
        mi2 = match(colnames(self$ypreds[[kk1]]),names(constk))
        constk1 = constk[mi2] ##can give rise to NA columns
        
        nas = which(is.na(mi2))
        if(length(nas)>0){
          st_end_constk = data.frame(lapply(names(constk), function(xz) as.numeric(strsplit(xz,"\\|")[[1]])))
          
          constk1[nas] = unlist(lapply(nas, function(xz) {
          
            spl = as.numeric(strsplit(  colnames(self$ypreds[[kk1]])[[xz]],"\\|")[[1]])
            if(spl[[1]]>max(st_end_constk[1,])){
              return (max(constk))
            }
            if(spl[[2]] <min(st_end_constk[2,])){
              return(min(constk))
            }
            st1 = which(st_end_constk[1,]==spl[[1]]); end1 = which(st_end_constk[2,]==spl[[2]])
         
            if(length(st1)==0 || length(end1)==0) return(NA)
            mean(constk[st1:end1])
          }))
          
        }
#       ab1= liability1(ab, constk, liab=liab) 
       aop=1-ao_probs(x_trans,constk1, betas,  lambdas[[1]], type="cum")
   
        self$ypreds[[kk1]][ind_1,] =  aop;
        #colnames(self$ypreds[[kk1]]) = names(constk)
      }else{
        ab = x_trans %*% betas
        levs1 = names(constk)
        attr(ab, "levs") =levs1
        self$ypreds[[kk1]][ind_1,] = ab
        #colnames(self$ypreds[[kk1]]) = colnames(betas);
      }
    }else if(family=="binomial"){
      ab = x_trans %*% betas
      for(jkk in 1:ncol(ab)){
        ab[,jkk] = if(liab) pao_(ab[,jkk]+ constants[[jkk]], lambdas[[jkk]]) else ab[,jkk]+ constants[[jkk]]
      }
      
        self$ypreds[[kk1]][ind_1,] =ab
        colnames(self$ypreds[[kk1]]) = colnames(ab)

    }else{
      kk_1 = 1

      constants = unlist(prev_kj$constants_proj[[kk1]])
      ab = x_trans %*% betas
      for(jkk in 1:ncol(ab)){
        fam = gaussian(link = sinh_asinh_link(lambdas[[jkk]]))
        ab[,jkk] = fam$linkinv(ab[,jkk]+ constants[jkk])
      }
  
        self$ypreds[[kk1]][ind_1,] = ab
        colnames(self$ypreds[[kk1]]) = colnames(ab)

    }
    
    if(!is.null(na_x) && length(which(na_x))>0){
      self$ypreds[[kk1]][na_x,] = rep(NA, ncol(self$ypreds[[kk1]]))
    }
    
  },

 predictions=function(nonNA, flip=FALSE){
   ypreds = self$ypreds
   ind_1 = if(flip) nonNA else !nonNA  ## because we setting to NA
   i1 = which(ind_1);
   l1 = length(i1)
  # if(length(which(ind_1))==0) return ypreds
   lapply(ypreds, function(yp){
     yp1 = yp;
     if(l1>0){
       for(kk in 1:ncol(yp1)){
         yp1[i1,kk]=NA
       }
     }
     yp1
   })
 },
calcRMSV=function(y, nonNA,    flip=FALSE, betas=NULL,
                  nvar = NA, cv_index = NA, beam =NA
                  ){
  phensi = self$phensi
  ypreds = self$ypreds
  types_ =self$types_
  nme_phens = names(phensi)
  names(nme_phens) = nme_phens
  ind_1 = if(flip) !nonNA else nonNA
  w1= if(is.null(ind_1)) self$weights else self$weights[ind_1]
  nsamps = length(which(ind_1))
  #nme_p1 = nme_phens[[1]]
  res_phens=lapply(nme_phens, function(nme_p1){
    #print(nme_p1)
    beta = if(is.null(betas))NULL else betas[[nme_p1]]
    family=  getOption("fspls.family",strsplit(nme_p1,"\\.")[[1]][1])
    fam = family
    mtype = match(family,names(types_))
    if(length(which(is.na(mtype)))>0) stop("could not find match in types_")
    types_i = types_[[mtype]]
    names(types_i) = types_i
    ycols = phensi[[nme_p1]]
    ypreds1 = ypreds[[nme_p1]]
  #  names(ycols) = names(ypreds1)  ##SAFE TO REMOVE
    #  print(names(ycols))
    phens = dimnames(ypreds1)[[2]]
    ycol_inds = 1:length(ycols)

    y2 = if(is.null(ind_1)) y[[nme_p1]] else  y[[nme_p1]][ind_1,,drop=FALSE]
    yp2 = if(is.null(ind_1))ypreds1 else  ypreds1[ind_1,,drop=FALSE] 
    names(ycol_inds) = dimnames(y2)[[2]][ycols]
    rms_1=lapply(ycol_inds, function(ycol_ind){
      beta1 = if(is.null(beta)) NA else beta[nrow(beta),ycol_ind]
      #print(ycol_ind)
      ycol = ycols[ycol_ind]
      y1 =  y2[,ycol]
      yp =if(fam=="ordinal" ) yp2 else yp2[,ycol_ind,drop=FALSE] 
    #  yp = inverse_func(yp)
      
      nonNAy = nonNA & !is.na(yp[,1])  &  !is.na(y1)
      names(types_i) = types_i
     aa2= .merge1_new(
        lapply(types_i, function(type_i1){
          #print(type_i1)
          type_i1s = strsplit(type_i1,"\\.")[[1]]
          type_i= type_i1s[[1]]
          thresh = if(length(type_i1s)>1) type_i1s[2] else NA
          if(length(which(nonNAy))<=1){
            return(NULL)
#            rms = c(NA,NA,NA);   names(rms) = c("low", "mid","high")
            
          }else{
          rms =tryCatch({
            .scoreInternal(yp[nonNAy,,drop=FALSE], y1[nonNAy],w1[nonNAy], type_i, fam,thresh, liab=self$liab)
          },error=function(w){
           print("warning");
            print(w);
            
            return(NULL)
          })
          }
         #  print("HERE");
        #   print(rms)
          if(is.null(rms)) {
            warning(paste("rms NULL typei1, nme_p1, ycol_ind, cv_index, nvar, beam", type_i1, nme_p1, ycol_ind, cv_index, nvar, beam))
            return(NULL) ;
          }
          if(is.null(names(rms))) names(rms) = 1:length(rms)
          if(typeof(rms)=="list"){
            return(rms)
          }else{
            df12=data.frame(cbind(names(rms), rms)) #, type_i) #as.matrix(rms)
            names(df12)=c("submeasure","value")#,"measure")
            df12[,2]=as.numeric(df12[,2])
          }
          df12
        }), num_cols="value", addName="measure")
     if(is.null(aa2)) return(NULL)
      aa2|> tibble::add_column(nsamps = nsamps, cv=flip, beta = beta1, sign = sign(beta1))
    })
    
    rms_2 = .fixBeforeMerge(rms_1)
    df3=.merge1_new(rms_2,num_cols="value", addName="pheno")
    if(!("subpheno" %in% names(df3)) && !is.null(df3)){
      df3 = df3 |> tibble::add_column(subpheno="")
    }
    df3    
  })
  #aa=readRDS("aa.rds")
 # nme_aa = names(aa[[1]])
#  aa_new = lapply(aa, function(aa1){
#     aa1[,match(nme_aa,names(aa1))]
#  })
  res_aa = .merge1_new(res_phens,addName="family")
  #lapply(aa, function(x)names(x))
  if(is.null(res_aa)) return(res_aa)
  res_aa |> pivot_wider(names_from="submeasure", values_from="value")
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
  TR = length(yp)
  a=binom.confint(M,TR, methods="probit",conf.level=getOption("conf.level",0.95))
  res = c(a$lower, a$mean, a$upper) 
  names(res)=c("low","mid","high")
  res
}
.calcSens<-function(yp, y1,w, thresh=0.5){
  P =length( which( y1==1))
  N =length( which( y1==0))
  TP = length(which(y1[yp>=thresh]==1))
  
  #  TN = spec*N
  #  o = order(yp, decreasing=FALSE)
  #  cum=cumsum(1-y1[o])
  #  pos=which(cum>=TN)[1]
  
  #  TN1 = length(which(y1[o][1:pos]==0))
  
  #  TP=length(which(y1[o][(pos+1):length(y1)]==1))
  
  if(P==0) return (c(NA,NA,NA))
  a=binom.confint(TP,P, methods="prop.test",conf.level=getOption("conf.level",0.95))
  res=c(a$lower, a$mean, a$upper)
  names(res)=c("low","mid","high")
  res
}
.calcPrev<-function(yp,y1,w,thresh =0.5){
  P =length( which( y1==1))
  N =length( which( y1==0))
  if(P+N==0) return (c(NA,NA,NA))
  a=binom.confint(P,P+N, methods="prop.test", conf.level=getOption("conf.level",0.95))
  
  res=c(a$lower, a$mean, a$upper)
  names(res)=c("low","mid","high")
  res
}
.better<-function(rmsv2, rmsv){
  subinds = which(rmsv2$submeasure=="mid" & !is.na(rmsv2$value) & !is.na(rmsv$value))
  rmsv2 = rmsv2[subinds,,drop=FALSE]
  rmsv =  rmsv[subinds,,drop=FALSE]
  if(nrow(rmsv2)!=nrow(rmsv)) stop("need same dims")
  mult = rep(1, nrow(rmsv2)) 
  mult[tolower(rmsv$measure) %in% c("correlation","auc","sens","spec","f1","youden")] =  -1  ## things where bigger is better, mult by -1   
  
  rmsv2$value= rmsv2$value* mult
  rmsv$value = rmsv$value*mult
  better = rmsv2$value<rmsv$value
  names(better) = rmsv2$measure
  better
}

.calcF1<-function(yp, y1,w, thresh=0.5){
  TP = length(which(y1[yp>=thresh]==1))
  FP = length(which(y1[yp>=thresh]==0))
  FN = length(which(y1[yp<thresh]==1))
  TN = length(which(y1[yp<thresh]==0))
  if(2*TP+FP+FN==0) return(rep(NA,3))
  a=binom.confint(2*TP,2*TP+FP+FN, methods="probit",conf.level=getOption("conf.level",0.95))
  res=c(a$lower, a$mean, a$upper)
  names(res)=c("low","mid","high")
  res
}


.calcNPV<-function(yp, y1,w, thresh=0.5){
  FN = length(which(y1[yp<thresh]==1))
  TN = length(which(y1[yp<thresh]==0))
  if(TN+FN==0) return(c(NA,NA,NA))
  #print(paste(TN, FN))
  a=binom.confint(TN,TN+FN, methods="probit",conf.level=getOption("conf.level",0.95))
  res=c(a$lower, a$mean, a$upper)
  names(res)=c("low","mid","high")
  res
}

.calcPPV<-function(yp, y1,w, thresh=0.5){
  TP = length(which(y1[yp>=thresh]==1))
  FP = length(which(y1[yp>=thresh]==0))
  if(TP+FP==0) return(c(NA,NA,NA))
  a=binom.confint(TP,TP+FP, methods="probit",conf.level=getOption("conf.level",0.95))
  res=c(a$lower, a$mean, a$upper)
  names(res)=c("low","mid","high")
  res
}
.calcSpec<-function(yp, y1,w, thresh =0.5){
  P =length( which( y1==1))
  N =length( which( y1==0))
  TN = length(which(y1[yp<thresh]==0))
  #TP = sens*P
  
  #o = order(yp, decreasing=TRUE)
  #cum=cumsum(y1[o])
  #pos=which(cum>=TP)[1]
  
  #TN = length(which(y1[o][(pos+1):length(y1)]==0))
  #TP1 = length(which(y1[o][1:pos]==1))
  if(N==0) return (c(NA,NA,NA))
  a=binom.confint(TN,N, methods="prop.test",conf.level=getOption("conf.level",0.95))
  
  res=c(a$lower, a$mean, a$upper)
  names(res)=c("low","mid","high")
  res
}


.calcPPV<-function(yp, y1,w, thresh =0.5){ #TP/TP+FP
  TP = length(which(y1[yp>=thresh]==1))
  FP = length(which(y1[yp>=thresh]==0))
  if(TP+FP==0) return( rep(NA,3))
  a=binom.confint(TP,TP+FP, methods="prop.test",conf.level=getOption("conf.level",0.95))
  res=c(a$lower, a$mean, a$upper)
  names(res)=c("low","mid","high")
  res
}

.calcNPV<-function(yp, y1,w, thresh =0.5){ #TN/TN+FN
  TN = length(which(y1[yp<thresh]==0))
  FN = length(which(y1[yp<thresh]==1))
  if(TN+FN==0) return (rep(NA,3))
  a=binom.confint(TN,TN+FN, methods="prop.test",conf.level=getOption("conf.level",0.95))
 res= c(a$lower, a$mean, a$upper)
 names(res)=c("low","mid","high")
 res
}

