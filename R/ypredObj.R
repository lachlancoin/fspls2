#.mergeRMSV<-function(datas){
#  rmsv_=.merge1(lapply(datas, function(d)d$rmsv_),num_cols="value",addName="data")
#  rmsv1 = .findMinRMSV(rmsv_, params$mult, fspls.sum=getOption("fspls.sum",T))
#  min(rmsv1)
#}


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
                       numvar=NULL, DRS_thresh = 1e-5,kk=1,
                       constants= prev_kj$constants_proj[[kk]]
                       #constants= unlist(prev_kj$constants_proj)
){
  yp = .rep(constants, length(which(ind_1)))
  betas =  prev_kj$betas
  vars1 = prev_kj$var
  if(length(vars1)==0) return(yp)
  if(!is.null(numvar) && numvar==0) return(yp) 
  betas1 = betas[[kk]]
  
  vars1 = prev_kj$var
  #print(names(prev_))
  # print(prev_kj$constants_proj)
  if(!is.null(numvar))vars1 = vars1[1:numvar]
  #  const = prev_kj$const[[i]][[ycol]]
  df1 = as.matrix(data.frame(vars1))
  
  
  
  
  if(nrow(df1)>0){
    for(kk in 1:length(data)){
      inds_11 = which(df1[1,]==kk)
      if(length(inds_11)>0){
        ##CHECK IF CORRECT
        yp = yp + data[[kk]][ind_1,df1[2,inds_11], drop=F]%*% betas1[inds_11,,drop=F]
      }
    }
    #  yp = yp+ unlist(prev_kj$constants_proj)
  }
  as.matrix(yp) 
}



calcRMS<-function( predy,yTs, family , CI = F, rmsea=T){
  if(length(yTs)==0) return(c(NA,NA,NA))
  rms = NA
  alias = which(!(is.na(predy) | is.na(yTs)))
  y = predy
  if(length(y)!=length(yTs)) stop(paste("error",length(y),length(yTs)));
  
  ##
  #  sort(unlist(lapply(1:10,function(i) {
  #     rands=alias[sample.int(length(alias), replace=T)]
  #    sqrt(sum((yTs[rands] - y[rands])^2)/length(yTs[rands]))
  #    })))
  rms= sqrt(sum((yTs[alias] - y[alias])^2)/length(yTs[alias]))
  if(CI){
    if(rmsea){
      cis = try(ci.rmsea(rms, df=1, N=length(yTs[alias]), conf.level = 0.95))
    }else{
      cis = c(NA, rms,NA)
    }
    return(cis) 
    
  } 
  rms*  (length(yTs)/length(alias))
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
.calcRMSV_1<-function(phensi,ypreds, data, nonNAs, 
                    types_ = getOption("fspls.types", 
                                       fromJSON('{"gaussian": "rank_correlation","binomial":"AUC","multinomial":"AUC_all","ordinal" : "AUC_all"}')),
                                                           
                    flip=F,
                    family=data$family){
 # inds1 = 1:maxl
  mtype = match(family,names(types_))
  if(length(which(is.na(mtype)))>0) stop("could not find match in types_")
 # names(mtype)=family
#  print(mtype)
#  names(inds1)=names(ypreds)[inds1]
 # rms_3=lapply(inds1, function(jj){
      ycols = phensi
      names(ycols) = names(ypreds)[phensi]
      if(is.null(names(ypreds))) names(ycols) =dimnames(ypreds)[[2]]
    #  print(names(ycols))
      phens = names(ypreds) 
      rms_1=lapply(ycols, function(ycol){
        ind_1 = if(flip) !nonNAs[[ycol]] else nonNAs[[ycol]]
        nsamps = length(which(ind_1))
        types_i = types_[[mtype[[ycol]]]]
        names(types_i) = types_i
        .merge1_new(lapply(types_i, function(type_i1){
        # print(type_i1)
          type_i1s = strsplit(type_i1,"\\.")[[1]]
          type_i= type_i1s[[1]]
          thresh = if(length(type_i1s)>1) type_i1s[2] else NA
        y1 =  if(is.null(ind_1)) data$y[,ycol] else  data$y[ind_1,ycol]
        fam = family[[ycol]]
        yp = if(is.null(ind_1))ypreds[[ycol]] else  ypreds[[ycol]][ind_1,,drop=F] 
        w1= if(is.null(ind_1)) data$weights[,ycol] else data$weights[ind_1,ycol]
        nonNA = !is.na(y1)
        nonNA = nonNA & !is.na(yp[,1])
      #  print(type_i)
        if(type_i=="misclass"){
          if(fam=="gaussian") y1=as.numeric(factor(y1, sort(unique(y1))))-1
          yp=if(fam=="binomial") plogis(yp[,1]) else yp[,1]
          rms = -1*.misclass(yp,y1, w1,auc=T)
        }else if(type_i=="area"){
          rms = -1 *.areaBetween(yp[nonNA,1], y1[nonNA], family=fam)[2]
          names(rms)=names(data$y)[[ycol]]
        }else if(type_i=="area_full"){
          rms = -1 *.areaBetween(yp[nonNA,1], y1[nonNA], family=fam)
         # names(rms)=names(data$y)[[ycol]]
          names(rms)=paste(names(data$y)[[ycol]],c("low","mid","high"),sep=".")
        }else if(type_i=="AUC"){
          #.calcAUCW(ypred, y, w)[2] for weighted AUC
          rms = -1*(.calcAUCW(yp[,1],y1, data$weights[nonNA,ycol]))
          names(rms)=c("low","mid","high")
        #  print(rms)
          #rms[i] = -1*calcAUC(yp,data$y[,ycol],T)
        }else if(type_i=="AUC_full"){
          #.calcAUCW(ypred, y, w)[2] for weighted AUC
          rms = -1*(.calcAUCW(yp[,1],y1, data$weights[,ycol]))
          names(rms)=paste(names(data$y)[[ycol]],c("low","mid","high"),sep=".")
          #  print(rms)
          #rms[i] = -1*calcAUC(yp,data$y[,ycol],T)
        }else if(type_i=="AUC_all"){
          if(fam=="ordinal"){
          levs = min(y1, na.rm=T):(max(y1,na.rm=T)-1) 
          lev_inds = 1:length(levs)
          names(lev_inds)=unlist(lapply(levs, function(l)paste(l,l+1,sep="|")))
            rms_l = .merge1_new(lapply(lev_inds, function(kj){
              gt =y1>levs[kj]
              if(length(which(gt))==0) return(NA)
              y2 = ifelse(gt,1,0)
             data.frame(value= -1*(.calcAUCW(yp[,kj],y2, data$weights[,ycol])[2]), submeasure=c("low","mid","high"))
            }),addName="subpheno")
            #names(rms_l )=paste(names(data$y)[[ycol]],levs,sep=".")
            rms = rms_l #+0.5
#            rms = c(rms_l, sum(rms_l, na.rm=T))
#            names(rms) = c(levs,"sum")
          }else{
          #.calcAUCW(ypred, y, w)[2] for weighted AUC
          levs = unique(as.character(y1[nonNA]))
          mi4 = match(.slug(levs), .slug(dimnames(yp)[[2]]))
          levs = levs[!is.na(mi4)]
          mi4 = match(.slug(levs), .slug(dimnames(yp)[[2]]))
          lev_inds = 1:length(levs)
          names(lev_inds)=levs
          
          rms_l = .merge1_new(lapply(lev_inds, function(kj){
            y2 = ifelse(y1==levs[kj],1,0)
            data.frame(list(value=-1*(.calcAUCW(yp[nonNA,mi4[kj]],y2[nonNA], data$weights[,ycol]))),submeasure=c("low","mid","high"))
          }),addName="subpheno")
          rms = rms_l#+0.5
#          rms = c(rms_l, sum(rms_l, na.rm=T))
#          names(rms) = c(levs,"sum")
#          rms = sum(rms_l, na.rm=T)
          }
          #rms[i] = -1*calcAUC(yp,data$y[,ycol],T)
        }else if(type_i=="correlation"){
          if(length(type_i1s)==1){
              rms =  if(length(yp[nonNA,1])==0 || var(yp[nonNA,1])==0) c(NA,0,NA) else -1* weightedCorr(yp[nonNA,1], y1[nonNA], weights=w1[nonNA],method="Pearson")
              #names(rms)=dimnames(data$y)[[2]][[ycol]]
              }else{
            if(var(yp[nonNA,1])==0) {
            rms = c(NA,0,NA) 
          }else{
            if(length(yp[nonNA,1])<3 || var(yp[nonNA,1])==0) {
             rms = c(NA,0,NA) 
            }else{
           ab=   ci_cor(yp[nonNA,1], y1[nonNA], method="pearson", probs = c(0.05, 0.95))
         rms = -1*c(ab$interval[1], ab$estimate, ab$interval[2])
            }
          }
            names(rms)=paste(dimnames(data$y)[[2]][[ycol]],c("low","mid","high"), sep=".")
          }
           #                                use="pairwise.complete.obs")
          
        }else if(type_i=="rank_correlation"){
          rms = if(var(yp[nonNA,1])==0) rep(0, ncol(data$y)) else  -1* weightedCorr(yp[nonNA,1], y1[nonNA],weights=w1[nonNA], method="Spearman") 
         # print(var(yp[nonNA]))
          #, use="pairwise.complete.obs")
          if(is.null(rms)){
            print(yp[nonNA,1])
            print(y1[nonNA])
          }
         # print(rms)
        #  print(names(data$y))
          names(rms)=names(data$y)[[ycol]]
        }else if(type_i=="rms"){
          rms=calcRMS(yp[,1],y1)
          names(rms)=names(data$y)[[ycol]]
          }else if(type_i=="youden"){
            rms=-1*.youden(yp,y1)[2]
            names(rms)=names(data$y)[[ycol]]
           } else if(type_i=="youden_full"){
            rms=-1*.youden(yp,y1)
            names(rms)=paste(names(data$y)[[ycol]],c("low","mid","high"),sep=".")
#            names(rms)=names(data$y)[[ycol]]
          }else if(type_i=="youden_sens"){
            rms=-1*.youden(yp,y1,typ="sens")[2]
            names(rms)=names(data$y)[[ycol]]
          }else if(type_i=="youden_spec"){
            rms=-1*.youden(yp,y1,typ="spec")[2]
            names(rms)=names(data$y)[[ycol]]
          }else if(type_i=="youden_sens_full"){
            rms=-1*.youden(yp,y1,typ="sens")
            names(rms)=paste(names(data$y)[[ycol]],c("low","mid","high"),sep=".")
          }else if(type_i=="youden_spec_full"){
            rms=-1*.youden(yp,y1,typ="spec")
            names(rms)=paste(names(data$y)[[ycol]],c("low","mid","high"),sep=".")
          }else if(type_i=="sens_spec"){
            #.calcAUCW(ypred, y, w)[2] for weighted AUC
            rms = .calcSensSpec(yp,y1,w1)
            #print(rms)
            #rms[i] = -1*calcAUC(yp,datas[[i]]$y[,ycol],T)
         }else if(type_i=="sens"){
            #.calcAUCW(ypred, y, w)[2] for weighted AUC
            rms = .calcSens(yp,y1,w1, thresh=thresh)
            #print(rms)
            #rms[i] = -1*calcAUC(yp,datas[[i]]$y[,ycol],T)
          }else if(type_i=="npv"){
            #.calcAUCW(ypred, y, w)[2] for weighted AUC
            rms = .calcNPV(yp,y1,w1, thresh=thresh)
            #rms[i] = -1*calcAUC(yp,datas[[i]]$y[,ycol],T)
          }else if(type_i=="f1"){
            #.calcAUCW(ypred, y, w)[2] for weighted AUC
            rms = .calcF1(yp,y1,w1, thresh=thresh)
            #rms[i] = -1*calcAUC(yp,datas[[i]]$y[,ycol],T)
          }else if(type_i=="ppv"){
            #.calcAUCW(ypred, y, w)[2] for weighted AUC
            rms = .calcPPV(yp,y1,w1, thresh=thresh)
            #rms[i] = -1*calcAUC(yp,datas[[i]]$y[,ycol],T)
          }else if(type_i=="spec"){
            #.calcAUCW(ypred, y, w)[2] for weighted AUC
            
            rms = .calcSpec(yp,y1,w1, thresh=thresh)
            #rms[i] = -1*calcAUC(yp,datas[[i]]$y[,ycol],T)
        }else{
          stop(paste("type not recognised:",type_i))
        }
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
      df3=.merge1_new(rms_1,num_cols="value", addName="pheno")
    df3    
  ##})
  #.merge1_new(rms_3, num_cols = "value", addName="beam")
  # dimnames(res_df)[[1]] = names(rms_3)
  # as.matrix(res_df)
}

.initYpred1<-function(data, family){
  #inds1 = 1:maxn
  #names(inds1)=inds1
  ycols = 1:ncol(data$y)
  names(ycols) = names(data$y)
  #lapply(inds1, function(x)  {
    res = lapply(ycols, function(i){
      yi = data$y[,i]
      if(family[[i]]=="multinomial"){
         if(!is.factor(yi)) stop(" yi needs to be a factor for multinomial")
         levs1 = levels(yi)
         names(levs1) = levs1
         rr1 = unlist(lapply(levs1, function(l1){
           vv = ifelse(yi ==l1, 1, 0)
           mean(vv,na.rm=T)
#           rep(mean(vv, na.rm=T), length(vv))
         }))
         rr = do.call(rbind, replicate(length(yi),rr1, simplify=FALSE))
         dimnames(rr) = list(dimnames(data$y)[[1]],levs1)
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
          dimnames(rr) = list(dimnames(data$y)[[1]],levs1[-length(levs1)])
    }else{
      rr = as.matrix(data.frame(rep(mean(yi,na.rm=T),length(yi))))
      dimnames(rr) = list(dimnames(data$y)[[1]],"X0")
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
  
  initialize=function(data,weights,   family = data$family){
   # rms_prev = 99999
  wname=NULL
  rmsv=NULL
  self$family=family
 # mult1 = params$mult
 # type =  unlist(lapply(family, function(f1) params$types_[[f1]]))
  if(is.null(weights)){
   
    ypreds=.initYpred1(data, family)
  }else{
    stop("not sure about this")
    beam2=10000  ##IF USING WEIGHTS WE DONT WANT THIS TO EFFECT
    weights=lapply(weights, function(w){
      w[,c(T,!is.na(match( dimnames(w)[[2]][-1],names( todoInds)))),drop=F]
    })
    
    weights=weights[unlist(lapply(weights, function(w) ncol(w)))>1]
    names(weights)=lapply(weights, function(w) paste(dimnames(w)[[2]][-1], collapse=","))
    if(length(weights)==0) stop("no weights left")
    ##expand the names to include all substrings
    wname=unique(unlist(lapply(names(weights), function(x){
      strspl=strsplit(x,",")[[1]]
      unlist(lapply(1:length(strspl), function(xi){
        paste(strspl[1:xi],collapse=",")
      }))
    }),rec=F))
    ##find max branching
    mn=max(.len(.branch(wname)))
    ypreds=.initYpred1(data, mn*2, family)
    
    #    print(wname)
    params$max=max(unlist(lapply(weights,function(x) ncol(x)-1 )))
    #print(paste("max", params$max))
  }
 # self$params=params
  self$wname=wname
  self$ypreds=ypreds
  self$weights=weights
  #self$rmsv_ = rmsv_
 # self$rms_prev = rms_prev
#  self$rmsv=rmsv
#  list(params=params,wname=wname,ypreds = ypreds, weights = weights,rms_prev=rms_prev, rmsv=rmsv)
  },
#calcYpred(prev_kj,self$data,ind_1,levs,numvar,kk1, kk)
  calcYpred=function(prev_kj, data, ind_1,  kk1, kk,na_x, levs1){  ## kk1 in model space 
    family = self$family[[kk]]
    #      ypred$ypreds[[kk]]$calcYpred(prev_kj,self$data,ind_1,levs,numvar,kk1, self$family[[kk]])
    if(family=="multinomial"){
      levs = dimnames( self$ypreds[[kk]])[[2]]
      ab=.calcYpred_multinom(prev_kj,  data, ind_1, levs, kk=kk1)  ## for multi-prediction
      mi22 = match(dimnames( self$ypreds[[kk]])[[2]], levs)
      self$ypreds[[kk]][ind_1,is.na(mi22)]=0
      self$ypreds[[kk]][ind_1,!is.na(mi22)] =  ab[,mi22[!is.na(mi22)]] 
    }else if(family=="ordinal"){
      constants = prev_kj$constants_proj[[kk1]]
      if(length(constants)<length(levs1)-1){
        constants = c(constants, rep(constants[length(constants)],length(levs1)-1-length(constants)))
      }
      constants = constants[1:length(levs1)-1]  ##if too long.  We assume that the lower counts are common but might not have higer counts
   aa= .calcYpred_ord(prev_kj,  data, ind_1, levs = levs1,
                                                 kk=kk1, constants = constants)  ## for multi-prediction
    self$ypreds[[kk]][ind_1,] =  aa
    }else if(family=="binomial"){
      constants = prev_kj$constants_proj[[kk1]]
      self$ypreds[[kk]][ind_1,] =  .calcYpred_binomial(prev_kj,  data, ind_1, 
                                                    kk=kk1, constants = constants)  ## for multi-prediction
    }else{
      ab = .calcYpred_1(prev_kj,  data, ind_1,kk=kk1) 
      self$ypreds[[kk]][ind_1,] =  unlist(ab[,1]) ##why need unlist here??
    }
    self$ypreds[[kk]][na_x,] = rep(NA, ncol(self$ypreds[[kk]]))
    
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

