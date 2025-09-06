default_types=fromJSON('{"gaussian": "correlation","binomial" : "AUC","multinomial" : "AUC","ordinal" :"AUC"}')

.getZetaBinary<-function(x,y, w){
  levsy =sort(unique(y[!is.na(y)]))
  yn = .mkBinary(y)
  m1=glm(yn~x,  family="binomial", weights=w)
  ll1 = logLik(m1)
  ll2 =  logLik(update(m1, ~1))
  pv = pchisq((2*(ll1 - ll2)),attr(ll1,"df")[1]-attr(ll2,"df")[1],lower.tail=FALSE,log.p=F)
  beta_new1 = m1$coefficients[-1]
  inds = 1:(length(levsy)-1); 
  names(inds) = unlist(lapply(inds, function(kk) paste(levsy[inds[kk]], levsy[inds[kk]+1],sep="|")))
  levsy1=levsy[inds]; names(levsy1)= names(inds)
  if(is.vector(x)) {
    prod= x *beta_new1 
    }else{
     prod =  x%*% beta_new1
     prod = prod[,1]
    }
  const_term=unlist(lapply(levsy1, function(l1){
    y2 = .mkBinary(y, thresh = l1)
  
    m1=glm(y2~1, offset=prod, family="binomial", weights=w)
    -m1$coefficients[[1]]
  }))
  list(beta_new1 = beta_new1, const_term = const_term, pv = pv)
}

.mkBinary<-function(y, thresh=mean(y,na.rm=T)){
  yn = rep(NA, length(y))
  yn[y<=thresh]=0
  yn[y>thresh]=1
  yn
}
.calcPvalue<-function(x,y, beta_new1, yp1,w, family){   ## this seems to not work anymore for multinomial
  if(length(which(!is.na(y)))==0) return (0)
  if(is.matrix(x) || typeof(x)=="S4"){
    yp_new = x %*% beta_new1
  }else{
    yp_new = as.matrix(x * beta_new1, ncol=1)
    yp1 = matrix(yp1, ncol=1, nrow= nrow(yp_new))
  }
  if(family=="multinomial"){
    m1 = multinom_ridge(yp1,y,w)
    m2 = multinom_ridge(yp_new,y,w)
    ll2 = -0.5 * m2$dev  
    ll1 = -0.5 *m1$dev        
    pv1 = .lrt(ll2,ll1,length(levels(y)), 1, log.p=T)
  }else if(family=="ordinal"){
    pv1 =tryCatch({
      df1 = data.frame(cbind(y,as.matrix(yp1 )));
      df1$y=factor(df1$y, levels = sort(unique(df1$y)))  
      func = paste0("y~",paste(colnames(df1)[-1], collapse="+"))
      m1=polr(func,  data=df1,weights=w,Hess=T, method="logistic")
      df1[,2] =  yp_new[,1]
      m2=polr(func,  data=df1,weights=w,Hess=T, method="logistic")
      ll1 = logLik(m1)
      ll2 =  logLik(m2)
      .lrt(ll2,ll1,2,1,log.p=T)  
    },error=function(ew){
      m1=glm(y~yp1[,1],weights=w,family="gaussian")
      m2=glm(y~yp_new[,1], weights=w,family="gaussian")
      ll2 = logLik(m2)
      ll1 =  logLik(m1)
      .lrt(ll2,ll1,2,1, log.p=T)
    })
  }else{
    m1=glm(y~yp1[,1],weights=w,family=family)
    m2=glm(y~yp_new, weights = w, family=family)
    ll2 = logLik(m2)
    ll1 =  logLik(m1)
    pv1 = .lrt(ll2,ll1,2,1, log.p=T)
  }
  pv1
}

multinom_ridge<-function(x,y,w,lambda=NULL){
  if(is.null(dim(x))){
    one_inds = 1:2
    x = cbind(1,x)
  }else{
    one_inds = 1
  }
#  m1=(multinom(y~x,weights=w, trace=F))
#  sm1  = summary(m1)
  
  ridge=glmnet(x,y, family="multinomial", alpha=0, weights=w, lambda=lambda)
  mins = min(ridge$lambda)
#predy=predict(ridge,x,s=min(ridge$lambda))[,,1]
  rbeta <- coef(ridge,s=min(ridge$lambda))
  rdf = as.matrix(data.frame(lapply(rbeta, as.matrix)))
  dimnames(rdf)[[2]]= names(rbeta)
  devs = deviance(ridge)
  list(const = apply(rdf[one_inds,,drop=F],2,sum), beta = rdf[-one_inds,],dev= devs[length(devs)])
}

## gets matrices without NAs
.oneVRest<-function(y){
  res = data.frame(lapply(1:ncol(y), function(i){
    y1 = y[,i]
    if(!is.factor(y1)) y1 = factor(y1, levels = sort(unique(y1)))
    levs = levels(y1)
    names(levs) = levs
    inds = 1:length(levs)
    names(inds) = levs[inds]
    lapply(levs, function(lev){
      y2= rep(0, length(y1))
      y2[y1==lev]=1
      y2[is.na(y1)]=NA
      y2
    })
  }))
  rownames(res) = rownames(y)
 Matrix(as.matrix(res),sparse=T)
}
.expandAllvAll<-function(y){
  res = data.frame(lapply(1:ncol(y), function(i){
    y1 = y[,i]
    if(!is.factor(y1)) y1 = factor(y1, levels = sort(unique(y1)))
    levs = levels(y1)
    inds = 1:length(levs)
    names(inds) = levs[inds]
    unlist(lapply(inds[-1], function(ind){
      inds2 = 1:(ind-1)
      names(inds2) = levs[inds2]
      lapply(inds2, function(ind2){
        y1 = rep(NA,  length(y))
        y1[y==levs[[ind]]]=0
        y1[y==levs[[ind2]]]=1
        y1
      })
    }),rec=F)
  }))
  rownames(res) = rownames(y)
  m1 = Matrix(as.matrix(res),sparse=T)
}

.expandFactors<-function(y, max_cats=50){
  res = lapply(1:ncol(y), function(i){
    print(i)
    y1 = y[,i]
    if(!is.factor(y1)) y1 = factor(y1, levels = sort(unique(y1)))
    levs1 = levels(y1)
    if(length(levs1) > max_cats) {
      warning(paste("more than max_cats", names(y)[[i]],length(levs1), max_cats))
      return(NULL)
    }
    mat = Matrix(0, nrow = length(y1), ncol = length(levs1), dimnames = list(rownames(y), levs1),sparse=T)
    if(length(which(is.na(y1)))>0){
      mat[is.na(y1),]=rep(NA, length(levs1))
    }
    for(j in 1:length(levs1)){
      indsy1 = which(y1==levs1[[j]])
      if(length(indsy1)>0){
        mat[indsy1,j]=1
      }
    }
    attr(mat,"factor")=y1
    mat
  })

  names(res) = dimnames(y)[[2]]
  res = res[unlist(lapply(res, length))>0]
  res
}

.sumMatrices<-function(matrices, onlyAll=F){
  if(length(matrices)==0) return(NULL)
  if(onlyAll){
    nmes_m = names(matrices[[1]])
    if(length(matrices)>1){
      for(k in 2:length(matrices)){
        nmes_m = nmes_m[which(nmes_m %in% names(matrices[[k]]))]
      }
    }
  }else{
    nmes_m =unlist(lapply(matrices, function(m) names(m)))
    nmes_m = nmes_m[!duplicated(nmes_m)]
  }
  m=rep(0, length(nmes_m)); names(m) = nmes_m
  for(k in 1:length(matrices)){
      mi1 = match(nmes_m,names(matrices[[k]]))
    #  print(length(which(is.na(mi1))))
    #   test =(cbind(names(matrices[[k]][mi1[!is.na(mi1)]]), names((m[!is.na(mi1)]))))
    #   print(min(apply(test,1, function(x) x[1] ==x[2])))
      m[!is.na(mi1)] =m[!is.na(mi1)]+ matrices[[k]][mi1[!is.na(mi1)]]
     if(onlyAll){
       m[is.na(mi1)] = 9999
     }
      #m = m+matrices[[k]]
  }
  m
}
.combineAngles1<-function(angles1,cols_incl1, incl, topn=100){ 
  nme_trans = names(angles1[[1]][[1]]); names(nme_trans) = nme_trans
  nmes_angs1 = names(angles1); names(nmes_angs1)=nmes_angs1
  lapply(nme_trans, function(nmet){
    df= .merge1_new(lapply(nmes_angs1, function(inc1){
      angles2 = angles1[[inc1]]
      col_incl = cols_incl1[[inc1]]
      .merge1_new(lapply(angles2, function(angles3){
        mat= angles3[[nmet]]
        .merge1_new(apply(mat,1, function(x){
            xp = sort(x, partial = topn)[topn]
            inds = which(x <= xp & col_incl)
            score = x[inds]
            data.frame(cbind(score,inds))
        }), addName="pheno")
      }), addName="family")
    }), addName="data_type")
    df[order(df$score),]
  })
}
.combineAngles<-function(angles,cols_incl, incl,onlyAll=F, topn=100){
  names(incl)=incl
  nme_trans = names(angles[[1]][[1]][[1]]); names(nme_trans) = nme_trans
  comb_all1=lapply(nme_trans, function(nme_t1){
  #nme_t1 = nme_trans[[1]]
      comb_all = lapply(incl,function(inc1){
        matrices = lapply(1:length(angles), function(i){
          ang1 = angles[[i]][[inc1]]
          if(is.null(ang1)) return(NULL)
          col_incl = cols_incl[[i]][[inc1]]
          #ang1=angle1[[inc1]]
          ang2=ang1[[1]][[nme_t1]]
          cs = colSums(ang2)
          if(length(ang1)>1){
            for(jk in 1:length(ang1)){
              cs = cs+colSums(ang1[[jk]][[nme_t1]])
            }
          }
          cs[col_incl]
        })
        names(matrices) = names(angles)
        matrices[unlist(lapply(matrices, length))>0]
        if(length(matrices)==0) return(NULL)
        .sumMatrices(matrices, onlyAll=onlyAll)
      })
 
  top_angles=whichpart1(comb_all, n=topn, return_scores=T)
  t1 = .merge1_new(lapply(top_angles, function(ta){
    data.frame(list(names = names(ta), value=ta))
  }),addName="data_type")
  t1 = t1[order(t1$value),]
  subset(t1, value<999)
  })
  comb_all1
}


.check<-function(var1,var2){
  if(length(var1)!=length(var2)) stop("not match var")
  if(length(var1)>0){
    if(length(which(unlist(var1)!=unlist(var2)))>0) stop("problem with check !!") ## should also check contenst
  }
  return(T)
}

fitModel<-function(yTr,x1,offset=NULL, family=getOption("family","binomial"),
                   weights = rep(1, length(yTr)),
                   useBF = getOption("useBF",FALSE)){ #
  ll=0
  #n = dim(x)[2]
  
  if(is.null(offset)){
    if(family!="multinomial"){
      if(is.null(x1)) lm<-glm(yTr~1,family=family, weights=weights)
      else lm<-glm(yTr~x1,family=family, weights=weights)
    }else{
      #offset = cbind(rep(0,dim(offset)[1]),offset)  not sure about cbind  here
      if(is.null(x1)) invisible(capture.output({lm<-multinom(yTr~1)}))
      else invisible(capture.output({lm<-multinom(yTr~x1)}))
    }
    if(useBF){
      coe = summary(lm)$coeff
      ll =ll+abf(coe[2,])
    }else{
      ll=ll+logLik(lm)
    }
  }else{
  
      if(family!="multinomial"){
        if(is.null(x1)) lm<-glm(yTr~offset(offset),family=family, weights=weights)
        else lm<-glm(yTr~x1+offset(offset),family=family, weights=weights)
      }else{
        #offset = cbind(rep(0,dim(offset)[1]),offset)  not sure about cbind  here
        if(is.null(x1)) invisible(capture.output({lm<-multinom(yTr~offset(offset))}))
        else invisible(capture.output({lm<-multinom(yTr~x1+offset(offset))}))
      }
      if(useBF){
        coe = summary(lm)$coeff
        ll =ll+abf(coe[2,])
      }else{
        ll=ll+logLik(lm)
      }
  }
  ll
}

.convertOrdinalRanges<-function(f2){
  v2=unlist(lapply(names(f2), function(f3) {
    rnge=as.numeric(strsplit(f3,"\\|")[[1]])
    val = f2[[f3]]
    if(length(rnge)==1){
      v2 = val;
      names(v2) = rnge[1]
      return(v2)
    }
    v3 = rep(val, rnge[2]-rnge[1]+1)
    names(v3) = rnge[1]:rnge[2]
    v3
  }))
  minp = min(as.numeric(names(v2)))
  if(minp>1){
    ##extend teh first bback to 1
    v2_1 = rep(v2[1], minp-1)
    names(v2_1) = 1:(minp-1)
    v2 = c(v2_1,v2)
  }
  v2
}

.getWeights11_1<-function(prev,const=F, pvs=F){
  betas=lapply(prev, function(prev_i){
    phen_nmes = names(prev_i$betas)
    if(pvs){
      pv_proj = prev_i$pvs_proj
      #      dimn =  list( pheno=datanmes, subpheno=phen_nmes)
      #      dims = unlist(lapply(dimn, length))
      v1 = rep(NA, length(phen_nmes)) #array(dimnames =dimn, dim = dims)
      names(v1) = phen_nmes
      for(i in 1:length(v1)){
        v1[i]    = prev_i$pvs_proj[[i]]
      }
      
      
    }else if(const){
      const_proj = prev_i$constants_proj
      const_proj1 = 
        lapply(const_proj, function(f2){
          if(is.null(names(f2))){
            return(f2);
          }
          v2 = .convertOrdinalRanges(f2)
          
        })
      
      nconst=unlist(lapply(const_proj1, length),rec=T)
      dimn =  list(subpheno=phen_nmes , variables = 1:max(nconst))
      dims = c(length(phen_nmes), max(nconst))
      names(dims) = names(dimn)
      v1 = array(dimnames =dimn, dim = dims)
      for(i in 1:dims[[1]]){
        for(j in 1:dims[[2]]){
          v_offset=const_proj1[[i]][[j]] 
          v1[i,j]=const_proj1[[i]][[j]] 
          #vec1= prev_i$const_proj1[[j]] -(v_offset - v_offset[1]) ####CHECK THIS DOES NOT SEEM RIGHT
          ## old - not entirely sure about this constants, but seems ok 
          #prev_i$constants_proj[[i]][[j]] - prev_i$const_off_proj[[i]][[j]]  
          #v1[j,1:length(vec1)]    =vec1    
        }
      }
    }else{
      #      dimn =  list( pheno=datanmes, subpheno=phen_nmes,variables=prev_i$varnames )
      dimn =  list( subpheno=phen_nmes,variables=unlist(lapply(prev_i$varnames, paste, collapse=".")) )
      dims = unlist(lapply(dimn, length))
      v1 = array(dimnames =dimn, dim = dims)
      for(i in 1:dims[[1]]){
        v1[i,]    = prev_i$betas[[i]]      
      }
      
      
    }
    v1
    #    v2=data.frame(t(v1))
    #    names(v2) = dimnames(v1)[[1]]
    #    v2
  })
  betas
}

.checkDatas<-function(datas){
  for(k in 1:length(datas[[1]]$data)){
    dn1 = dimnames(datas[[1]]$data[[k]])[[2]]
    
    
    for(j in 2:length(datas)){
      dn2 = dimnames(datas[[j]]$data[[k]])[[2]]
      mi1 = match(dn2, dn1)
      if(mi1[[1]]!=1 || mi1[[length(mi1)]]!=length(mi1)) stop("problem")
    }
  }
}



.convertToTrainingData1<-function(depmapData,  drug,
                                  nme="all",
                                  small_class_thresh = 5,
                                   incl_data=NULL
                                        , phenotypes = NULL){
  .conv1=function(x) strsplit(x,"::")[[1]][1]
  
  if(is.null(depmapData$nme2Treatment) || !is.null(phenotypes)){
    drug_code = drug
    names(drug_code) = drug
  }else{
    drug = drug[which(drug %in% names(depmapData$nme2Treatment))]
    mi1 = match(drug,names(depmapData$nme2Treatment))
    drug_code =depmapData$nme2Treatment[mi1[!is.na(mi1)]]
  }
  
  dat = lapply(depmapData$input_files, function(inp){
    zd = inp$zdesc[[nme]]
      lapply(zd,function(zd1){
        m1 = attach.big.matrix(zd1[['vals']])
        attr(m1,"norm")=attr(zd1,"norm")
        attr(m1,"mean_x")=attr(zd1,"mean_x")
        isna = attach.big.matrix(zd1[['isNA']])
        list(vals=m1, isNA=isna)
      })
  })
  pheno=NULL
  if(is.null(phenotypes)){
    zd = depmapData$pheno_files$pheno$zdesc[[nme]]
   pheno = 
     lapply(zd, function(zd1)attach.big.matrix(zd1[['vals']]))
    pheno = lapply(pheno, function(p) {
      p1=p$x; 
      dimnames(p1)[[2]] =unlist(lapply(dimnames(p1)[[2]], .conv1));
      p1
      })
  }
  beam=getOption("fspls.bean",c(1,1))
  #maxn= beam[1] *beam[2]
  
  trans_names = names(depmapData$transformations)
  names(subset_nmes) = subset_nmes
  if(!is.null(incl_data)) subset_nmes = subset_nmes[incl_data]
  types1 = names(dat)[names(dat)!="pheno"]
  names(types1) =types1
  names(trans_names) = trans_names
    data2 = lapply(types1, function(t1){
      lapply(trans_names, function(nme_t){
          d1 = dat[[t1]][[nme_t]][['vals']]
          if(is.null(d1)) {
            d1 = dat[[t1]][[nme_t]][['vals']]
          }
      d1
    })
    })
    data2Na = lapply(types1, function(t1){
      lapply(trans_names, function(nme_t){
        d1 = dat[[t1]][[nme_t]][['isNA']]
        if(is.null(d1)) {
          d1 = dat[[t1]][[nme_t]][['isNA']]
        }
        d1
      })
    })
    
    data = unlist(data2, recursive=F)
    dataNA = unlist(data2Na, recursive=F)
    names(data) = gsub("\\.",":", names(data))
    names(dataNA) = gsub("\\.",":", names(dataNA))
    y_matr = if(is.null(phenotypes)) pheno[[nme]] else phenotypes
    if(is.null(rownames(y_matr)[[1]]))stop(" row names cannot be null ")
    mi2 =match(rownames(data[[1]]), rownames(y_matr))
    
#    y = as.matrix(y_matr[,grep(drug_code, dimnames(y_matr)[[2]]),drop=F])
    mi=match(drug_code, dimnames(y_matr)[[2]])
    if(length(which(is.na(mi)))>0){
    #  stop(paste("did not find phenotype", paste(drug_code, collapse=",")))
      y = data.frame(array(0,dimnames =list(rownames(data[[1]]), drug_code) , dim = c(nrow(data[[1]]),length(drug_code))))
    
    } else{   
    #print(mi2)
      y = y_matr[mi2,mi,drop=F]
      dimnames(y)[[2]] = names(drug_code)
    }
    fam = getOption("family",NULL)
    for(k in 1:ncol(y)){
       if(!is.null(fam) && fam[k]=="multinomial"){
          y [,k] = factor(y[,k])
        }
      if(is.factor(y[,k])){
        if(length(levels(y[,k]))>2){
        ##levels from most to least abundant
        y[y[,k] %in% names(which(table(y[,k])<small_class_thresh)),]=NA
        
        y[,k] = factor(y[,k], levels = names(sort(table(as.character(y[!is.na(y[,k]),k])),decr=T)))
        }else{
        #  y[,k] = as.numeric(y[,k])-1
        }
      }
    }
    #print(head(y))
    if(length(drug)==1 && drug_code!=drug) dimnames(y)[[2]]=names(drug_code)
    weights = apply(y, c(1,2), function(c) 1)
    ones = weights
    y = data.frame(y)
    family=if(!is.null(getOption("family",NULL))) getOption("family") else .inferFamily(y)
    #print(paste("family",family))
    res=dataObj$new(data,dataNA,y,family, weights)
  res
}





dataObj<-R6Class("dataObj", public = list(
  dist="environment",
  data="list",
  dataNA="list",
  types="vector",
  nrow="integer",
  vars="list",
#  var_thresh="double",
  min_minor="double",
  prev="list",
  y="matrix",
 # ones="matrix",
#  ones_x="matrix",
  #ymod="matrix",
  weights="matrix",
  train="list", #trainObj",
  subset="logical",  ### used to restrict to subsets in looc Obj
  looc="loocObj",
  mean_x="list",
  norm="vector",
  norms_list="list",
  norm_done="list",   ##dataframe of prev calculated
  mem_dirs="list",  ##outfile to write new norms calculated in this session
 # cols_incl="list",
  ypreds_all="ypredObj",
  family="character",
 
  UDVP="UDVPObj",


## converts between text and indices
convert=function(b_i1){
  #print(b_i1)
  data_ind = which(names(self$data)==b_i1[[1]])
  if(length(data_ind)==0){
    print(names(self$data))
    stop(paste("could not find",b_i1[[1]]))
  }
  var_ind = which(dimnames(self$data[[data_ind]])[[2]]==b_i1[2])
  if(length(var_ind)==0) return(NULL)
  c(data_ind, var_ind)
},
cats = function(maxpheno = 1e9){
  fam=self$family
  phens =  lapply(1:length(fam), function(i){
    y1 = self$y[[i]]
    if(fam[[i]]=="multinomial"){
      y1 = attr(y1, "factor")
      l1 = list(levels(y1))
    }else{
      cn = colnames(y1)[1:min(maxpheno, ncol(y1))]
      l1 = list(apply(y1,2,table))
    }
    names(l1) = names(self$y)[[i]]
    l1
  })
  names(phens) = names(self$y)
},
pheno = function(maxpheno=1e9,sep=F, sep_group=F,code = NULL){ 
 phens =  lapply(self$y, function(y1) colnames(y1)[1:min(maxpheno, ncol(y1))])
 if(!is.null(code)){
   phens =  lapply(phens, function(p){
    mi2=match(code,p)
    p2 = p[mi2]
    names(p2) = names(code)
    p2
   })
 }
 nmes = names(phens)
 names(nmes)=nmes
 if(sep_group){
   return(lapply(nmes, function(nme){
    phens[nme]
   }))
 }
mult = grep("multinomial",names(phens))
  if(sep){
   
    nmes1 = if(length(mult)>0) nmes[-mult] else nmes
   l1= lapply(nmes1, function(nme){
     ph = phens[[nme]]
     inds = 1:length(ph)
     names(inds)= inds
    lapply(inds, function(i){
       l3 = list(ph[i])
       names(l3) = nme
       l3
     })
   })
   l12= lapply(nmes[mult], function(nme){
     ph = phens[[nme]]
       l3 = list(ph)
       names(l3) = nme
       l3
   })
   return( c(l12,(unlist(l1,rec=F))))
  }else{
    return(list(all=phens))
  }
},
####does regression just on orthogonal component
calcBetaProj1=function(subphens,k,b_i,prev_var, transform_func,convert=T, betas = list(), strict=F,project=F,  useglm=getOption("glmnet",T), useoffset=F){
  if(convert){
    b_i = self$convert(b_i)
    if(length(b_i)<2) {
      pvs = unlist(lapply(phensi,function(xx) 999)) ## returning large positive value since its not in dataset
       return(list(pvs=pvs))
    }
    prev_var = lapply(prev_var, self$convert)
    subphens = self$phensi(subphens)
  }
  nmesi = names(subphens)
  names(nmesi) = nmesi
  #nme=nmesi[[1]]
  res = lapply(nmesi, function(nme){
  #  print(nme)
    phensi1 =  match(nme, names(self$y))
    family=getOption("fspls.family",self$family[phensi1])
    phensi_=subphens[[nme]]
    betas1 = betas[[nme]]
    if(family=="multinomial") useoffset=F
    if(!is.null(betas1) && length(prev_var)>0 ){
      res1=self$calcBetaProjAll(nme,phensi_,family,  k, b_i, prev_var, transform_func,betas1,project=project, strict=strict,useglm=useglm, useoffset=useoffset)
    }else{
      res1 = self$calcBetaProj(nme,phensi_,family,  k, b_i, prev_var, transform_func,strict=strict, useglm=useglm)
    }
    res1
  })
  nme = names(res[[1]])
  
  names(nme)=nme
  res1 = lapply(nme, function(n){
    lapply(res, function(r) r[[n]])
  })
  res1
},

calcBetaProj=function(nme,phensi_,family, k,b_i,prev_var, transform_func, strict=F, useglm=getOption("glmnet",T)){
  #b_i = b_i1
  data = self
  train = data$train
  d = train
  # print(family)fcalz
  nonNA = self$looc$incl[,k]
  vars1 = c(prev_var,list(b_i))
  data$updateUDVP(prev_var)
  
  #train = data$train
  ys =self$y[[nme]]##  if(family %in% c("binomial","ordinal")) (d$y1) else if(family =="multinomial") d$y2 else d$y
  if(family=="multinomial"){
    ys = data.frame(list(attr(ys,"factor")))
    names(ys) = nme
  }else if(length(phensi_) < ncol(ys)){
    ys = ys[,phensi_,drop=F] 
  }
  ncoly = ncol(ys)
  betas = apply(ys,2,function(v1) list())# apply(train$y,2, function(v1) list())
  tbls =  apply(ys,2, function(v1) list())
  
  #  lapply(train, function(v) apply(v$y,2, function(v1) list()))
  constants = apply(ys,2, function(v1) list())
  #lapply(train, function(v) apply(v$y,2, function(v1) list()))
  pvs =  apply(ys,2, function(v1) list())
  
  #lapply(train, function(v) apply(v$y,2, function(v1) list()))
  
  #for(i in 1:length(train)){
  
  j = length(vars1)
  #UDV = UDVP_h[[i]]
  #     x = d$x[,vars1[j]]
  if(length(vars1[[j]])==0){
    non_na_x = rep(T, self$nrow)
    useglm=F
  }else{
    non_na_x = if(is.null(data$dataNA[[vars1[[j]][1]]])) rep(T,nrow(ys) ) else !(data$dataNA[[vars1[[j]][1]]][, vars1[[j]][2]] )
  }
  if(length(vars1[[j]])==0){
    x_=rep(1, self$nrow)
    }else{
  x1 =  data$data[[vars1[[j]][1]]][, vars1[[j]][2]] 
  mean_adj = data$mean_x[[vars1[[j]][1]]][vars1[[j]][2]]
  
  x_ = x1- mean_adj
    }
  UDV = data$UDVP ## should check it corresponds to prev_i
  if(!is.null(UDV$VDU)){
    if(!is.null(UDV$P2)) {
      x_ = UDV$P2 %*% x_
    }else{
      #print("opp")
      ## revisit this
      # print("calculating betas on projection")
      #x_ = UDV$P2() %*%x_ ## AFTER PROJECTION
      x_ = x_ - UDV$P %*% (UDV$VDU %*% x_)
      if(FALSE){
        ##this confirms orthogonality
       unlist(lapply( 1:length(UDV$var), function(jk){
         x2_ = self$data[[UDV$var[[jk]][1]]][,UDV$var[[jk]][2]]- self$mean_x[[UDV$var[[jk]][1]]][UDV$var[[jk]][2]]
         x2_%*%x_
      }))
      }
    }
  }
  # x = x[non_na_x]
  varx = var(x_[non_na_x],na.rm=T)
  if(is.na(varx) )stop("variance  NA")
  #if(varx<1e-10)stop("variance  too small")
  #k =NULL ##delete this later
  for(kk in 1:ncoly){
    nonNAk = nonNA & non_na_x
    x = x_[nonNAk ]
    y = transform_func(ys[,kk])[nonNAk]
    w = data$weights[nonNAk]
    beta_new1=0;
    
    if(family=="multinomial"){
      ty=as.list(table(y))
      tbls[[kk]] = ty[ty>0]
      
      rdf = try(multinom_ridge(x,y,w))
      if(inherits(rdf,"try-error")){  ## setting null values
        levsy = levels(y)
        beta_new1 = rep(0, length(levsy)); names(beta_new1) = levsy
        const_term = beta_new1
      }else{
      #logpv1 = pchisq(rdf$dev/2,df=1,low=F,log=T)
        beta_new1=rdf$beta
        const_term = rdf$const
      }
    }else if(family=="ordinal" ){
      ty=as.list(table(y))
      tbls[[kk]] = ty[ty>0]
      use_bin = length(unique(y,na.rm=T))<=2
      if(!use_bin){
        df = data.frame(cbind(y,x ))
        df$y = factor(y, levels = sort(unique(y,na.rm=T)))  
        # print("using polr") 
        m1=try(polr(y~x,  data=df,weights=w,Hess=T, method="logistic"),silent=T)
        #predict(m1, df, type = "p")
        
        if(inherits(m1,"try-error")){
          if(strict) warning("!!!")
          resz = .getZetaBinary(x,y, w)
          beta_new1=resz$beta_new1
          const_term =resz$const_term
          pv1 = resz$pv1
          #const_term = -1*sm$coefficients[1,1]
        }else{
          ll1 = logLik(m1)
          ll2 =  logLik(update(m1, ~1))
          pv = pchisq((2*(ll1 - ll2)),attr(ll1,"df")[1]-attr(ll2,"df")[1],lower.tail=FALSE,log.p=F)
          coeff = c(m1$coefficients,1,1,pv[1])
          const_term = m1$zeta
          beta_new1=coeff[1]
          pv1 = coeff[4] # pt(abs(coeff[3]),df=1,lower.tail=F)  ## not sure this correct
        }
        
      }else{
        stop("here")
        m1=glm(y~x, family="binomial", weights=w)
        sm  = summary(m1)
        coeff = sm$coeff[2,]
        pv1 = coeff[4]
        beta_new1=coeff[1]
        const_term = sm$coefficients[1,1]
      }
    }else{
      spike_slab_iter = getOption("spike_slab_iter",0)
      nonNAy = !is.na(y)
      if(length(which(nonNAy))==0){
        const_term =0; beta_new1 = 0
      }else if(spike_slab_iter>1){
        # print("using spike slab")
        #ab2 = lm.spike(x~y, niter=spike_slab_iter, ping -1)
        ab2=lm.spike(x~y, niter= spike_slab_iter, ping =-1)#,weights=w)
        sm = summary(ab2)
        pv1=max(1e-20,1-sm$coefficients[1,5])
        if(family=="binomial"){
          ab=logit.spike(y~x, niter= spike_slab_iter, ping =-1)
        }else{
          #if(getOption("logprint",F)) print("gaussian spike slab")
          ab=lm.spike(y~x, niter= spike_slab_iter, ping =-1)#,weights=w)
          
        }
        sm = summary(ab)
        x_ind = match(c("x","(Intercept)"),dimnames(sm$coefficients)[[1]])
        coeff = sm$coefficients[x_ind[1],c(1,2,3,5)]
        const_term = sm$coefficients[x_ind[2],1]
        beta_new1=coeff[1]
        pv1 = 1-coeff[4]
      }else{
        if(useglm){
          sm2<- tryCatch({
          ones = rep(1, length(x))
          x_mod = cbind(ones,x)
       
          ridge=glmnet(x_mod[nonNAy,,drop=F],y[nonNAy],weights = w[nonNAy], family=family, alpha = 0)
          
          rbeta <- coef(ridge,s=min(ridge$lambda))
          const_term = sum(rbeta[1:2,1])
          beta_new1 = rbeta[3,1]
          if(FALSE){
            ypred = predict(ridge,x_mod, s=min(ridge$lambda))
            ll1 = .logLik(y,ypred, family=family)
            ll2 =.logLik(y,family=family) 
            pv1 = .lrt(ll1,ll2,2,1, log.p=F)
          }else{
            nulldev=ridge$nulldev
            dev= ridge$dev.ratio[which.min(ridge$lambda)] *nulldev
            #       pv1 = pchisq(-1*(dev-nulldev),df=1,low=F)
          }
          list(const_term = const_term, beta_new1 = beta_new1)
          }, error=function(errw) {
            m1=glm(y~x, family=family, weights=w) ## including weights lead to non-convergence
            sm  = summary(m1)
            #print(var(x))
            if(nrow(sm$coeff)<2){
              coeff = rep(0,4)
              const_term=0
              pv1=1
              beta_new1 = 0
            }else{
              coeff = sm$coeff[2,]
              
              const_term = sm$coefficients[1,1]
              beta_new1=coeff[1]
              pv1 = coeff[4]
            }
            list(const_term = const_term,const_term, beta_new1 = beta_new1, pv1=pv1)    
          })
        }else{
        sm2 <- tryCatch({
          m1=glm(y~x, family=family, weights=w) ## including weights lead to non-convergence
          sm  = summary(m1)
          #print(var(x))
          if(nrow(sm$coeff)<2){
            coeff = rep(0,4)
            const_term=sm$coefficients[1,1]
            pv1=1
            beta_new1 = 0
          }else{
            coeff = sm$coeff[2,]
            
            const_term = sm$coefficients[1,1]
            beta_new1=coeff[1]
            pv1 = coeff[4]
          }
          
          if(FALSE){
            ypred=beta_new1*x + const_term
            ll1 = .logLik(y,ypred, family=family)
            ll2 =.logLik(y, family=family) 
            pv1 = .lrt(ll1,ll2,2,1, log.p=F)
          }
          # print(paste(pv1_1, pv1, family))
          # print(pv1)
          #pchisq((2*(ll1 - ll2)),attr(ll1,"df")[1]-attr(ll2,"df")[1],lower.tail=FALSE,log.p=F)
          list(const_term = const_term, beta_new1 = beta_new1, pv1=pv1)
        }, warning=function(errw) {
          print("using glmnet 1")
          ones = rep(1, length(x))
          x_mod = cbind(ones,x)
          nonNAy = !is.na(y)
          ridge=glmnet(x_mod[nonNAy,,drop=F],y[nonNAy],family=family,weights=w[nonNAy], alpha = 0)
          #ridge=glmnet(x_mod,y,family=family, alpha = 0)
          rbeta <- coef(ridge,s=min(ridge$lambda))
          const_term = sum(rbeta[1:2,1])
          beta_new1 = rbeta[3,1]
          if(FALSE){
            ypred = predict(ridge,x_mod, s=min(ridge$lambda))
            ll1 = .logLik(y,ypred, family=family)
            ll2 =.logLik(y,family=family) 
            pv1 = .lrt(ll1,ll2,2,1, log.p=F)
          }else{
            nulldev=ridge$nulldev
            dev= ridge$dev.ratio[which.min(ridge$lambda)] *nulldev
            #       pv1 = pchisq(-1*(dev-nulldev),df=1,low=F)
          }
          list(const_term = const_term,const_term, beta_new1 = beta_new1)
        })
        }
        const_term = sm2$const_term
        beta_new1 = sm2$beta_new1
        #  pv1 = sm2$pv1
        
      }
      
    }
    
    betas[[kk]] =beta_new1
    
    
    ##
    #yp = beta_new1%*%x1
    #a = lm(y~yp[1,])
    #summary(a)
    if(family=="multinomial" ){  ##this consistent for pvs?? easiest way to get pvalue for multinomial and single variable x
      m1=lm(x~y)
      ll1 = logLik(m1)
      ll2 =  logLik(update(m1, ~1))
      pv1 = .lrt(ll1,ll2,length(levels(y)),1, log.p=T)
    }else{
       pv1 = .calcPvalue(x,y, beta_new1, 1,w, family)
       
    }
    
    pvs[[kk]] = pv1
    constants[[kk]] = const_term  #-mean_adj*beta_new1
  }
  list(betas=betas, constants = constants,pvs = pvs,tbls=tbls)
},
##project variable controls whether the projection of the variable is fitted, or the variable itself
calcBetaProjAll=function(nme,phensi_,family, k,b_i,prev_var, transform_func,betas1, project=F, strict=F, useoffset = F,useglm=getOption("glmnet",T)){
  data = self
  if(!useoffset) project=F
  nonNA = self$looc$incl[,k]
  vars1 = c(prev_var,list(b_i)) 
  data$updateUDVP(prev_var)
  ys =self$y[[nme]]##  if(family %in% c("binomial","ordinal")) (d$y1) else if(family =="multinomial") d$y2 else d$y
  if(family=="multinomial"){
    ys = data.frame(list(attr(ys,"factor")))
    names(ys) = nme
  }else if(length(phensi_) < ncol(ys)){
    ys = ys[,phensi_,drop=F] 
  }
  ncoly = ncol(ys)
  betas = apply(ys,2,function(v1) list())# apply(train$y,2, function(v1) list())
  tbls =  apply(ys,2, function(v1) list())
  #Walls =  apply(ys,2, function(v1) list())
  
  constants = apply(ys,2, function(v1) list())
  pvs =  apply(ys,2, function(v1) list())
  j = length(vars1)
  non_na_x = if(is.null(data$dataNA[[vars1[[j]][1]]])) rep(T,nrow(ys) ) else !(data$dataNA[[vars1[[j]][1]]][, vars1[[j]][2]] )
  x_ =self$extractData(vars1, adjust=T)
  
  UDV = data$UDVP ## should check it corresponds to prev_i
  x1_ = x_[,ncol(x_), drop=F]
  if(useoffset && project){
        x1_ = x1_ - UDV$P %*% (UDV$VDU %*% x1_)
        if(FALSE){ ##THIS DEMONSTRATE ORTHOGONALITY
              unlist(lapply( 1:length(UDV$var), function(jk){
                x2_ = self$data[[UDV$var[[jk]][1]]][,UDV$var[[jk]][2]]- self$mean_x[[UDV$var[[jk]][1]]][UDV$var[[jk]][2]]
                x2_%*%x1_
              }))
        }
  }
  
  #mean_x = apply(x_,2,mean, na.rm=T)
  if(length(which(duplicated(colnames(x_)))>0)){
    print(vars1);
    stop("problem")
  }
  non_na_x = apply(x_,1,function(v) length(which(is.na(v))))==0
  for(kk in 1:ncoly){
    nonNAk = nonNA & non_na_x
    y = transform_func(ys[,kk])[nonNAk]
    w = data$weights[nonNAk]
    beta_new1=0;
    const_term=0
    Wall1 = NULL
    yp1 =  if(family=="multinomial")  x_[,-ncol(x_),drop=F ] %*%  betas1  else x_[,-ncol(x_),drop=F ] %*%  betas1 [,kk,drop=F]
    
    if(useoffset){
      if(project){ 
        UDVP1=UDVPObj$new(self, NULL,yp1) #,check=F, centralise=F)
        Wall1 = UDVP1$getWall(x_[,ncol(x_)], diag(1))
     # Walls[[kk]] = Wall1
      }else{
        Wall1 = diag(ncol(yp1))
      }
      #x = cbind(yp1, x[,ncol(x)])
      x = cbind(yp1[nonNAk,], x1_[nonNAk,1])
      dimnames(x)[[2]] = c(paste0("A",1:(ncol(x)-1)),"x")
    }else{
      x = x_[nonNAk,,drop=F]
    }
    yp1 = yp1[nonNAk,,drop=F]  
    if(family=="multinomial"){
       ty=as.list(table(y))
      tbls[[kk]] = ty[ty>0]
      
      rdf=multinom_ridge(x,y,w)
      const_term = rdf$const
      beta_new1 =  rdf$beta
    }else if(family=="ordinal" ){
      use_bin = length(unique(y))<=2
       ty=as.list(table(y))
      tbls[[kk]] = ty[ty>0]
      if(!use_bin){
        colnames(x)=paste("x", 1:ncol(x), sep="")
        df = data.frame(cbind(y,as.matrix(x) ))
        df$y = factor(y, levels = sort(unique(y)))  
        # print("using polr")
        func = paste0("y~",paste(colnames(x), collapse="+"))
        m1=try(polr(func,  data=df,weights=w,Hess=T, method="logistic"),silent=T)
        #predict(m1, df, type = "p")
        
        if(inherits(m1,"try-error")){
          if(strict) warning("polr!!!")
          resz = .getZetaBinary(x,y, w)  ##need as.matrix(x) ??
          beta_new1=resz$beta_new1
          const_term =resz$const_term
   #       pv1 = resz$pv1
          #  print(df)
          # print(m1)
#          yn = .mkBinary(y)
         
 #         m1=glm(yn~as.matrix(x), family="binomial", weights=w)
          #warning("polr error .. using gaussian!!")
          #m1=glm(y~x, family="gaussian", weights=w)
  #        sm  = summary(m1)
  #        beta_new1=sm$coefficients[-1,1]
  #        const_term = -1*sm$coefficients[1,1]
        }else{
          coeff = m1$coefficients
          const_term = m1$zeta
          beta_new1=coeff
        }
        
      }else{
        stop("here")
        m1=glm(y~as.matrix(x), family="binomial", weights=w)
        sm  = summary(m1)
        beta_new1=sm$coefficients[-1,1]
        const_term = -1*sm$coefficients[1,1]
      }
    }else{
      nonNAy = !is.na(y)
      spike_slab_iter = getOption("spike_slab_iter",0)
      if(length(which(nonNAy))==0){
        const_term=0; beta_new1 = rep(0,ncol(x))
      }else if(spike_slab_iter>1){
        # print("using spike slab")
        #ab2 = lm.spike(x~y, niter=spike_slab_iter, ping -1)
        ab2=lm.spike(as.matrix(x)~y, niter= spike_slab_iter, ping =-1)#,weights=w)
        sm = summary(ab2)
        if(family=="binomial"){
          ab=logit.spike(y~as.matrix(x), niter= spike_slab_iter, ping =-1)
        }else{
          #if(getOption("logprint",F)) print("gaussian spike slab")
          ab=lm.spike(y~x, niter= spike_slab_iter, ping =-1)#,weights=w)
          
        }
        stop("redo")
        sm = summary(ab)
        x_ind = match(c("x","(Intercept)"),dimnames(sm$coefficients)[[1]])
        coeff = sm$coefficients[x_ind[1],c(1,2,3,5)]
        const_term = sm$coefficients[x_ind[2],1]
        beta_new1=coeff[1]
      }else{
       
        if(useglm){
          sm2<- tryCatch({
            #if(ncol(x)==1){
           # ones = rep(1, nrow(x))
          #  x_mod = cbind(ones,x)
            #}else{
            #  x_mod = x
            #}
           
         
                ridge=glmnet(x[nonNAy,,drop=F],y[nonNAy],family=family,weights=w[nonNAy], alpha = 0)
                rbeta <- coef(ridge,s=min(ridge$lambda))
                const_term = rbeta[1,1]
                 beta_new1 =  rbeta[-1,1]
         
            list(const_term = const_term, beta_new1 = beta_new1)
          }, error=function(errw) {
            
                m1=glm(y~as.matrix(x), family=family, weights=w) ## including weights lead to non-convergence
            sm  = summary(m1)
            #print(var(x))
            if(nrow(sm$coeff)<1+ncol(x)){
              coeff = rep(0,4)
              const_term=0
              beta_new1 = rep(0, ncol(x))
            }else{
#              coeff = sm$coeff[2,]
              
              const_term = sm$coefficients[1,1]
              beta_new1=sm$coeff[-1,2]
              if(length(beta_new1)!=ncol(x)){
                print(beta_new1)
                print(colnames(x))
              }
              names(beta_new1)=colnames(x)
           #   pv1 = coeff[4]
            }
            list(const_term = const_term,const_term, beta_new1 = beta_new1)    
          }
)
          beta_new1 = sm2$beta_new1
          const_term = sm2$const_term
        }else{
          sm2 <- tryCatch({
          
             m1=glm(y~as.matrix(x), family=family, weights=w) ## including weights lead to non-convergence
            sm  = summary(m1)
            #print(var(x))
            if(nrow(sm$coeff)<1+ncol(x)){
              coeff = rep(0,4)
              const_term=0
              beta_new1 = rep(0, ncol(x))
            }else{
              coeff = sm$coeff[2,]
              
              const_term = sm$coefficients[1,1]
              beta_new1=sm$coeff[-1,2]
              # print(beta_new1)
              #  print(colnames(x))
              if(length(beta_new1)!=ncol(x)){
                print(beta_new1)
                print(colnames(x))
              }
              names(beta_new1)=colnames(x)
            }
            
            list(const_term = const_term, beta_new1 = beta_new1)
          }, warning=function(errw) {
            print("using glmnet 1")
            #ones = rep(1, nrow(x))
            #x_mod = cbind(ones,x)
            nonNAy = !is.na(y)
           
              ridge=glmnet(x[nonNAy,,drop=F],y[nonNAy],family=family, alpha = 0, weights=w[nonNAy])
              #ridge=glmnet(x_mod,y,family=family, alpha = 0)
              rbeta <- coef(ridge,s=min(ridge$lambda))
              #const_term = rbeta[1,1]
              const_term = rbeta[1,1]
              beta_new1 = rbeta[-(1),1]
            list(const_term = const_term,beta_new1 = beta_new1)
          })
          beta_new1 = sm2$beta_new1
          const_term = sm2$const_term
        }
      }
      
    }
      pv1 = .calcPvalue(x,y, beta_new1, yp1,w, family)
    
      
        
      
      pvs[[kk]] = pv1
      if(project){
       betas[[kk]] = if(family=="multinomial") Wall1 %*% beta_new1 else   (Wall1 %*% beta_new1)[,1]
      }else{
        betas[[kk]] = beta_new1
      }
      constants[[kk]] = const_term  #-mean_adj*beta_new1
  }
  list(betas=betas, constants = constants,tbls = tbls, pvs = pvs)
},
##var and W_all1 are from one smaller model
calcWall=function(b_i, var, W_all1){
  data = self
 # var = prev_i$var
  data$updateUDVP(var)
  UDV = data$UDVP ## should check it corresponds to prev_i
  x = data$data[[b_i[1]]][, b_i[2]] - data$mean_x[[b_i[1]]][b_i[2]]
  UDV$getWall(x, W_all1)
},
getConstantsProj=function(phensi){  ##default constants
  data = self
  family =unlist(lapply(names(phensi),function(x) getOption("fspls.family",  strsplit(x,"\\.")[[1]][1])))
  constants_proj = lapply(1:length(phensi), function(ik){
    nme = names(phensi)[[ik]]
    phensi1 = phensi[[ik]]
    #  if(length(phensi1)>1) stop("this should just have length 1")
    if(family[[ik]]=="multinomial") {
      fact = attr(data$y[[ik]],"factor")
      return(rep(0,length(levels(fact))-1 ))
    }else if(family[[ik]]=="ordinal"){
      levs = sort(unique(data$y[[ik]][,1]))
      return(rep(0,length(levs )-1))
    }else{
      return(rep(0,length(phensi1)))
    }
  })
  names(constants_proj) = names(phensi)
  constants_proj
},
extractData=function(var, adjust=T,convert=F){
  if(convert)var = lapply(var,self$convert)
  Dall = Matrix( 0,nrow = self$nrow, ncol = length(var), dimnames = list(dimnames(self$data[[1]])[[1]], names(var)), sparse=T)
  nme = rep("", length(var))
  if(length(var)==0)return(Dall)
  for(jk in 1:length(var)){
    v1 = var[[jk]]
    Dall[,jk] = self$data[[v1[1]]][,v1[2]];
    nme[[jk]] = paste(names(self$data)[[v1[1]]], dimnames(self$data[[v1[1]]])[[2]][v1[2]])
    if(adjust)Dall[,jk] = Dall[,jk]- self$mean_x[[v1[1]]][v1[2]]
  }
  dimnames(Dall)[[2]] = nme
  Dall
},

#vars2 = vars_all[[1]]$variables[[1]];phens1 = phens[[1]]; k=1; useoffset=T
makeModels=function(phens1, vars2,var_transf, k, 
                    project=T,logpthresh = -5,useglm=T,useoffset=T,
                  func_str1="function(y) y"
                  ){
  phen2 = phens1
  #if(length(phens1)>1) stop("assuming just one phenotype here")
  ypred=self$ypred(phen2)
  phensi = self$phensi(phen2)
  subphens = phensi
 # prev_is = lapply(fold_inds, function(k) stateObj$new(phensi,self,self$train[[k]], k))
 prev_i= NULL; #self$prev[[k]]  #lapply(fold_inds, function(k) self$train$prev[[k]])
  #if(length(prev_i$var)>0) stop("problem")
  nmes = c()
  #mean_y = self$train$means_y[[k]]
 
  data = self;
  len = length(vars2)
  models = vector("list", len)
  useglm=getOption("glmnet",T)
  fams1 = lapply(names(phens1), function(st1)strsplit(st1,"\\.")[[1]][1])
  family = getOption("fspls.family",fams1) #ypred$family[[1]]
  if(family[[1]]=="multinomial") useoffset=F
  nme1 = ""
  b_i_name=c()
#  func_str1=func_str[1]## could think about updating this for new iteration, but for moment stick
  #transform_func = eval(str2lang(func_str[[1]]))
  
  prev_i = self$makeNextModel(prev_i,b_i_name,subphens,k, func_str1[[1]],
                              family, ypred=ypred, project=project, useglm=F, logpthresh =logpthresh,useoffset=useoffset)
  models[[1]] = prev_i$simplify(names(func_str1)[[1]])
  nmes[[1]] = "empty"
  jk=1
  while(jk<=len){
    b_i_name = vars2[[jk]]
#    transform_func=eval(str2lang(func_str1[var_transf[[jk]]]))
    prev_i1 = self$makeNextModel(prev_i,b_i_name,subphens,k,func_str1[[var_transf[[jk]]]],
                                 family, ypred=ypred, project=project, useglm=useglm, logpthresh =logpthresh,useoffset=useoffset)
    if(is.null(prev_i1)) break;
    nme2 = paste(vars2[[jk]], collapse=".")
    nme1 = if(jk==1)  nme2 else paste(nme1, nme2,sep=";")
    nmes[[jk+1]] = nme1
    models[[jk+1]] =prev_i1$simplify(var_transf[[jk]])
    prev_i = prev_i1
    jk = jk+1
  }
  models = models[1:length(nmes)]
  names(models) = nmes
  models
},
updateWeights=function(subphens=self$pheno()[[1]][1]){
  y1 = self$y[[names(subphens)[1]]]
  fam = strsplit(names(subphens)[[1]],"\\.")[[1]][1]
  if(fam=='multinomial'){
    fact=attr(y1,"factor")
  }else if(fam!="gaussian"){
    colind = which(dimnames(y1)[[2]]==subphens[[1]][1])
    if(length(colind)==0){
      print('h')
      colind=which(dimnames(y1)[[2]]==names(subphens[[1]])[1]) 
    }
    fact = apply(y1[,colind ,drop=F  ],1,paste,collapse="_")
  }
  if(fam=="gaussian"){
   warning('not renormalising gaussian')
    return(NULL)
  }else{
    naInds = fact=="NA" | is.na(fact)
    tbl = table(fact[!naInds])
    w1   =(1.0/tbl)
    #w1 = w1/sum(w1)
    self$weights = rep(0, length(self$weights))
    for(kk in 1:length(w1)){
      self$weights[which(fact==names(tbl)[[kk]])] = w1[[kk]]
    }
    self$weights[naInds] = min(w1)
  }
},
 makeNextModel=function(prev_i, b_i_name, subphens, k,funcst, family, ypred=NULL, project=T, useglm=T,    logpthresh = -5,
                        useoffset=T,inv_funcst=NULL,
                        CHECK=getOption("fspls.check",F),
                        verbose=getOption("fspls.verbose1",F)) {
   data =self
   transform_func = eval(str2lang(funcst))
    b_i = if(length(b_i_name)==0) c() else self$convert(b_i_name)
    #if(is.null(b_i)) return(NULL)
    mean_x =if(length(b_i_name)==0) c() else  self$mean_x [[b_i[[1]]]][b_i[2]]
    prev_var = if(is.null(prev_i)) list() else prev_i$var
    self$updateUDVP(prev_var)
    #W_all = data$calcWall(b_i, prev_i$var, prev_i$W_all) ## WALL not important, we can get rid of it later
    #prev_var = if(jk==1) prev_i$var  else lapply(vars2[1:(jk-1)], self$convert)
    betas = prev_i$betas
  
    b_new_proj = self$calcBetaProj1(subphens,k,b_i,prev_var, transform_func, betas = betas, project=project,convert=F,  
                                    strict=T, useglm=useglm, useoffset=useoffset) 
   # b_new_proj1 = self$calcBetaProj1(subphens,k,b_i,prev_var, transform_func, betas = betas, project=!project,convert=F,    strict=T, useglm=useglm, useoffset=useoffset) 
    
    betas_new = b_new_proj$betas
      constants_proj =if(family[[1]]=="multinomial") b_new_proj$constants[[1]] else b_new_proj$constants
      pvs =if(family[[1]]=="multinomial") b_new_proj$pvs[[1]] else b_new_proj$pvs
   #   print(pvs)
    #  if(.sumChisq(pvs)>logpthresh ){
        
     # return(NULL)
    #  }
    tbls = if(family[[1]]=="multinomial") b_new_proj$tbls[[1]]  else b_new_proj$tbls
    prev_i1=stateObj$new(subphens,data, betas_new,constants_proj, tbls, k,prev_i , b_i,b_i_name=b_i_name, mean_x = mean_x, W_all = NULL,pvs =pvs, useoffset=useoffset)
  #  prev_i1$setOffset() 
    if(is.null(prev_i)) return(prev_i1)
    prev_i1$setOffset()
    if(FALSE){ ##JUST TO CHECK WHAT GLM VARIABLES  WOULD BE IF WE JUST FIT ALL
      extractd0 = self$extractData(c(prev_i$var, list(b_i)), adjust=F)
      means_x = apply(extractd0,2,mean)
      extractd = self$extractData(c(prev_i$var, list(b_i)), adjust=T)
      subnme = names(subphens)[[1]]
      family = strsplit(subnme,"\\.")[[1]]
      df2 = data.frame(lapply(subphens[[1]], function(colind){
       
      #  family = getOption("fspls.family",self$family)
       # phensi = self$phensi(subphens)
        if(family=="multinomial"){
          y_m = attr(self$y[[1]],"factor")
          m1 = try(multinom(y~as.matrix(extractd0),weights=w, trace=F))
          if(inherits(m1,"try-error")){
            stop("multinom error")
          }
          sm  = summary(m1)
          betas_n1 = t(sm$coefficients[,-1])
          xM = t(t(extractd0 %*% betas_n1)+ sm$coefficients[,1])
          attr(xM, "levs") = levels(y_m)
          x_p = liability(xM)
          roc(self$y[[1]][,1],x_p[,1])
          roc(self$y[[1]][,2],x_p[,2])
          roc(self$y[[1]][,3],x_p[,3])
          list(betas_n1, prev_i1$betas[[1]])
        }else{
          y11 = self$y[[subnme]][,subphens[[subnme]],drop=F]
          self$updateWeights(subphens)
          nonNAy = !is.na(y11[,1])
        ridge=glmnet(cbind(1,extractd0[nonNAy,]),y11[nonNAy,,drop=F] ,family=family, alpha = 0, weights =self$weights[nonNAy])
        rbeta <- coef(ridge,s=min(ridge$lambda))
        aa = predict(ridge,cbind(1,extractd0[nonNAy,]),s=min(ridge1$lambda), family=family)
        ridge1 = glmnet(cbind(1, extractd0[nonNAy,] %*%  rbeta[-(1:2),,drop=F]), y11[nonNAy,,drop=F], family = family, alpha=0)
        rbeta1 <- coef(ridge1,s=min(ridge1$lambda))
        print(rbeta1)
        return(cbind(rbeta[-(1:2),],prev_i1$betas[[1]]))
        }
#        (rbeta[-(1:2),] - prev_i1$betas[[1]][,colind,drop=F])
      }))
      print(df2)
#      dimnames(df2)[[2]] = 1:ncol(self$y[[1]])
      }
  
    if(verbose  && !is.null(ypred) && !is.null(inv_funcst)){
      inv_func = eval(str2lang(inv_funcst))
#      inv_transform_func = lapply(inverse_funcstr, function(str1) eval(str2lang(str1)))  ## should be inverse
      ypred$updateYP(data, prev_i1, nonNA, inv_func=inv_func,flip=FALSE, liab=F)
      new_const = prev_i1$updateConst(subphens,ypred, data,k, transform_func, useglm=getOption("glmnet",T),verbose=verbose, update=F)
      
      ypred$updateYP(data, prev_i1, nonNA,  inv_func=inv_func,flip=FALSE, liab=T)
      rmsv=(ypred$calcRMSV(self$y, nonNA,     flip=FALSE))
      print(quantile(rmsv$value))
    #  median(ypred$ypreds$binomial)
    }
  return(prev_i1)
   
},
nreps=function(){
  ncol(self$looc$incl)
},


getNonNAInds=function(inds){
    df=data.frame(lapply(inds, function(k){
      self$looc$incl[,k]#  $nonNA1(k, self$y)
      #nonNAs[[kk]] 
    }))
    apply(df,1,min)!=0
},

plotData=function(vars_all1, phens1 = vars_all1$phens, all_types=F, transform_x = NULL, violin=F, assoc=F){
  phensi = self$phensi(phens1)
  nmei = names(phensi); names(nmei) = nmei
  df = data.frame(lapply(nmei, function(nmei1){
    as.matrix(self$y[[nmei1]][,phensi[[nmei1]],drop=F])
  }))
  
  variables = vars_all1$variables
  names(variables)=NULL
  vars = unlist(variables, recursive=F)
  incls = names(self$data); names(incls)=incls
  
  if(all_types){
    genes = unique(unlist(lapply(vars, function(x) x[[2]])))
    names(genes) = genes
    vars1 = lapply(incls, function(incl){
      lapply(genes, function(g){
       c(incl,g)
      })
    })
    vars = unlist(vars1,rec=F)
  }
  into=c("data","gene")
  
  if(!is.null(transform_x) && all_types){
    nmev = names(vars)
    transform_x1 = fromJSON(transform_x)
    nme_t = names(transform_x1)
    to_repl=paste0("_",names(transform_x1))
    for(kk in nme_t){
      nmev = gsub(paste0("_",kk), paste0(".",kk), nmev)
    }
    names(vars) = nmev
    into=c("data","transform","gene")
  #  facet="transform~data"
  }
 
  
  df2 = data.frame(lapply(vars, function(vark){
   inds = self$convert(vark)
   x = self$data[[inds[1]]][,inds[2]]
   na_x = self$dataNA[[inds[1]]][,inds[2]]
   x[which(na_x)] = NA
   x
  }))
  names(df2)=unlist(lapply(vars, function(x)paste(x,collapse="__")))
 
  
  if(assoc){
    nmes_df = names(df);names(nmes_df) = nmes_df;
    nmes_df2 = names(df2);names(nmes_df2) = nmes_df2;
    
    pvs=lapply(nmes_df, function(nme_df_){
      lapply(nmes_df2, function(nme_df2_){
          m2=lm(df2[[nme_df2_]] ~df[[nme_df_]])
          m1=glm(df2[[nme_df2_]] ~1)
        ll2 = logLik(m2)
        ll1 =  logLik(m1)
        pv1 = .lrt(ll2,ll1,2,1, log.p=F)
        pv1
      })
    })
    pv_res = sort(unlist(pvs))
    print(pv_res)
    return(pv_res)
  }
  
  nme_df = names(df); names(nme_df) = nme_df
  df4 = .merge1_new(lapply(nme_df, function(nmedf1){
    df_k = df[[nmedf1]]
    df3=df2 %>% tibble::add_column(y=df_k) %>% pivot_longer(names(df2)) %>% separate("name",sep="__", into=into)
  }),addName="pheno")
  df4$y = factor(df4$y)
  df4 = subset(df4, !is.na(y))
  
  df4
},

extractPredictions=function(all_models_,phens1,inverse_func_str, flags, CV = FALSE,liab=T,
                            ypred = self$ypred(phens1)){
  d = self
  transform_func = lapply(inverse_func_str, function(str1) eval(str2lang(str1)))  ## should be inverse
  
  nmesp = names(phens1)
  names(nmesp)= nmesp
  #all_models1_ = all_models_[[1]]
  all_models1 = all_models_
 # res_all = lapply(all_models_, function(all_models1_){
  #  all_models1 = all_models1_[[nme_p]]
 # evals = .merge1_new(lapply(all_models, function(all_models1){
    full_ind = names(all_models1)=="full"
    full_model = all_models1[["full"]]
    nmesm = names(all_models1)[!full_ind];
    inds=as.numeric(nmesm)
    nmes_models = names(all_models1[[1]]);
    names(nmes_models) = nmes_models
  #  ypred = ypredObj$new(d,NULL, d$family)
    evals = lapply(nmes_models, function(nmes1){
    #  dim(ypred$ypreds[[1]])
      if(!is.null(full_model) && ! CV){
        nonNA =self$looc$incl[,self$nreps()]
        prev_i1 = full_model[[nmes1]]
        inv_func=transform_func[[prev_i1$transf]]
        ypred$updateYP(d,prev_i1, nonNA, inv_func=inv_func, flip=FALSE, liab=liab )

        return(lapply(nmesp, function(nmesp1){
          phens1[[nmesp1]]
          yy2 = d$y[[nmesp1]]
          mi  = match(phens1[[nmesp1]],dimnames(yy2)[[2]])
          list(y=d$y[[nmesp1]][,mi,drop=F], ypred= ypred$ypreds[[nmesp1]])
        } )    )   
#        d$updateYpredsInds(phens,full_model[[nmes1]], d$nreps(), ypred )
      }
      if(length(nmesm)>0 && CV){
        inds = as.numeric(nmesm)
        for(j in 1:length(nmesm)){
          nonNA =self$looc$incl[,inds[[j]]]
          prev_i1 = all_models1[[j]][[nmes1]]
          inv_func=transform_func[[prev_i1$transf]]
          ypred$updateYP(self, prev_i1, nonNA,inv_func=inv_func, flip=TRUE, liab=liab)
#          d$updateYpredsInds(phens,all_models1[[j]][[nmes1]], inds[[j]], ypred)
        }
        nonNA=self$getNonNAInds(inds)
        return(lapply(nmesp, function(nmesp1){
          yy2 = d$y[[nmesp1]]
          mi  = match(phens1[[nmesp1]],dimnames(yy2)[[2]])
        list(y=d$y[[nmesp1]][!nonNA,mi,drop=F], ypred=ypred$ypreds[[nmesp1]][!nonNA,,drop=F]) 
        }))
       
      }
      return(NULL)
      #     rbind(res1,res2)
      #    }),addName="variable")
    })
    #,addName="trainedOn")
    #}),addName="fullmodel")
  
  
 # })
#res_all[lapply(res_all,length)>0]
evals
},
ypred=function(phens1){
  family = unlist(lapply(names(phens1), function(str)getOption("fspls.family",strsplit(str,"\\.")[[1]][1])))
  ypr = ypredObj$new(self,self$phensi(phens1),family=family)
  ypr
},# inverse_func_strs = fromJSON(.readFlag(flags,"transform_y_inverse",'{"y":"function(y) y"}'))
#all_models_y = all_models$y; inverse_func_str = fromJSON(flags1$transform_y_inverse)[[1]]; self = datasAll$datas[[1]]
evaluateAllModels=function(all_models_y,phens,inverse_func_str, flags,
                           ypred = self$ypred(phens), #lapply(phens, function(phens1) self$ypred(phens1)),
                           verbose=F
                         ){
  d = self
  self$updateLOOC(phens,flags)
  liab = .readFlag(flags,"liab",T)  ## whether to evaluate with liability , default is true
  transform_func = lapply(inverse_func_str, function(str1) eval(str2lang(str1)))  ## should be inverse
  
#  ypred = self$ypred(phens)
  group_names= names(all_models_y); names(group_names)=group_names
  numvars = unlist(lapply(group_names, function(x) if(x=="empty") 0 else length(strsplit(x,";")[[1]])))
  numvars1 = sort(unique(numvars))
  names(numvars1) = numvars1
  #pheno_nmes = names(phens); names(pheno_nmes)=pheno_nmes
  if(length(all_models_y)==0) return(NULL)
  nmes_models = names(all_models_y[[1]][[1]]) #[[1]]);
  names(nmes_models) = nmes_models
  #numvar = numvars[[1]]; nmes1 = nmes_models[[1]];group_name = group_names[[1]]
  evals_all = .merge1_new(lapply(numvars1, function(numvar){
    if(verbose)print(paste("numvar",numvar))
    #.merge1_new(lapply(pheno_nmes, function(pheno_nme){
      #ypred = ypreds[[pheno_nme]]; 
      if(is.null(ypred)) stop("ypred is null")
      .merge1_new(lapply(nmes_models, function(nmes1){
        group_names2 = group_names[numvars==numvar]
        evals = .merge1_new(lapply(group_names2, function(group_name){
          if(verbose) print(paste(numvar,nmes1,group_name))
          all_models1 = all_models_y[[group_name]]#[[pheno_nme]]   
              full_ind = names(all_models1)=="full"
              full_model = all_models1[["full"]][[nmes1]]
              all_models2 = lapply(all_models1[!full_ind], function(am) am[[nmes1]])
              if(length(all_models2)>0){
                all_models2 = all_models2[!unlist(lapply(all_models2, is.null))]
              }
              nmesm = names(all_models2)
              inds=as.numeric(nmesm)
              res1 = NULL; res2 = NULL
              if(!is.null(full_model)){
                #ypredObj$updateYP(self, phens, )#= self$looc$incl[,k2]
                nonNA =self$looc$incl[,self$nreps()]
                ypred$updateYP(d, full_model, nonNA, inv_func=transform_func[[full_model$transf]], flip=FALSE, liab=liab)
                res1 = ypred$calcRMSV(self$y, nonNA,      flip=FALSE)%>% tibble::add_column(isfull=T, transf=full_model$transf)
                #res1 = self$getRMSVInds(phens, d$nreps(), ypred)  
              }
              if(length(nmesm)>0){
                transf=c()
                for(j in 1:length(nmesm)){
                  nonNA =self$looc$incl[,inds[[j]]]
                  prev_i1 = all_models2[[j]]
                  transf = c(transf,prev_i1$transf)
                  ypred$updateYP(d, prev_i1, nonNA,inv_func=transform_func[[prev_i1$transf]], flip=TRUE)
                  #          self$updateYpredsInds(phens,all_models1[[j]][[nmes1]], inds[[j]], ypred)
                }
                nonNA=self$getNonNAInds(inds)
                res2 = ypred$calcRMSV(self$y,nonNA, flip=TRUE)%>% tibble::add_column(isfull=F, transf=paste(unique(transf), collapse="_"))
              }
              rbind(res1,res2)
        }),addName="model")
        if(is.null(evals) ) return(NULL)
        ## this calculates average scores across cv
        evals1 = unite(evals,"phenomodel","pheno","model")
        evals$isfull[evals1$phenomodel %in% evals1$phenomodel[evals1$isfull]] = T
        evals3 = subset(evals, cv==T)%>% pivot_wider(names_from="pheno", names_prefix="pivoted_")
        if(nrow(evals3)==0) return(evals)
        evals4 = unite(evals3,"comb","submeasure","measure","subpheno","family", remove=F)
        combs = unique(evals4$comb)
        mi = match(c("comb","submeasure","measure","subpheno","family","cv","isfull","model","transf"),names(evals4))
        avg_inds =( 1:ncol(evals4))[-mi]
       # print(evals4)
        avgs = .merge1_new(lapply(combs, function(comb1){
          s1 = subset(evals4, comb==comb1)
          s2 = s1[1,,drop=F]
          avg_v = apply(s1[,avg_inds,drop=F],2,mean,na.rm=T)
          s2[avg_inds] = as.list(avg_v)
          s2$model = "avg"
          s2
        }))
        avgs1 = avgs[,-1]%>%pivot_longer(cols=starts_with("pivoted_"),names_prefix="pivoted_",names_to="pheno")
        mi2 = match(names(evals),names(avgs1))
      
        rbind(evals, avgs1[,mi2])
      }),addName="trainedOn")
   # }),addName="pheno_group")
  }),addName="numvars")
  if(!is.null(evals_all$numvars)){
  evals_all$numvars = as.numeric(evals_all$numvars)
  }
  evals_all
},



  updateModel=function(k,best_all_i, model_prev,to_keep,CHECK=T){
    if(TRUE) stop("not updating")
    sprev = self$train$prev[[k]]
    best_all=lapply(1:length(best_all_i), function(j){
      prev_i= sprev[[to_keep[j]]]
      best_i = best_all_i[[j]]
      if(CHECK){
        .check(prev_i$var, model_prev[[to_keep[j]]]$var)
      }
      #data_ind = which(names(self$data)==b_i[[1]])
      #b_i1 = c(data_ind, which(dimnames(self$data[[data_ind]])[[2]]==b_i[2]))
      nxt_v= lapply(best_i, function(b_i) stateObj$new(self, self$train[[k]], k,prev_i,self$convert(b_i1)))
         lapply(nxt_v, function(nv) nv$updateConst(self, self$train[[k]],k))
      #                       nxt_v = nxt_v[!unlist(lapply(nxt_v, is.null))]
      names(nxt_v) = unlist(lapply(nxt_v,  function(nv)paste(unlist(lapply(nv$var, function(vv){
        paste(names(self$data)[vv[1]],dimnames(self$data[[vv[1]]])[[2]][vv[2]],sep=".")
      })), collapse=",")))
      nxt_v
    })
    self$train[[k]]$prev_old = self$train[[k]]$prev
    self$train[[k]]$prev=unlist(best_all, rec=F)
  },
  getMaxBetaProj=function(){
    #     lapply(self$prevs, function(prevk){
    mabv = lapply(self$prev, function(pk) max(abs(unlist(pk$betas_proj))))
    #   })
    mabv
  },
  updateUDVP=function(var){
    UDV=self$UDVP
    if(!is.null(UDV) && !is.null(UDV$var)){
      if(length(UDV$var)==0 &&length(var)==0) return(NULL)
      if(length(UDV$var)==length(var) &&  length(which(unlist(UDV$var)!=unlist(var)))==0) return(NULL) ## dont need to update
    }
    Dall = self$extractData(var, adjust=T)
    self$UDVP=UDVPObj$new(self, var,Dall)
  },
projOut1=function(ik){
  UDV=self$UDVP
#  d = self$train[[k]]
 # nonNA = d$nonNA
  x = self$data[[ik]] ## would be better not to do this to keep as bigmatrix
#  mean_x=self$mean_x[[ik]]
  #Dall = t(t(x[d$nonNA,])-mean_x)
  if(length(UDV$var)>0){
    alias = UDV$alias
#    if(length(alias)<length(which(nonNA))) stop(" problem")
    Ut = UDV$Ut
    Dinv = UDV$Dinv
    Vinv = UDV$Vinv
    #I = UDV$I
    #P2 = UDV$P2
   # Ut=t(U)
    #W = (Vinv %*% Dinv) %*% (t(U) %*% Dall[alias,,drop=F])
    #UD = Ut %*% Dall[d$nonNA,]  ##t(U) %*% Dall[alias,,drop=F]
    #mean_x1 = .rep(mean_x,nrow(x) )
    ## assumes rows of Ut sum to zero
    #Ut[,!d$nonNA]=rep(0, ncol(UD))
    UD = if(isbigmatrix(x)) dgemm(A=Ut,B= x) else Ut %*%x
    VinvDinv = Vinv %*% Dinv
    #UD2 = Ut %*% mean_x1
    #for(kk in 1:nrow(Ut)) UD1[kk,] = UD1[kk,]-Ut[kk,]*mean_x
    W =  if(isbigmatrix(x)) dgemm(A=VinvDinv, B=UD) else VinvDinv%*% UD
 #PW = dgemm(A=UDV$P, B=W)
    ##P is dim nx1 and W is 1xm
    #means of PW are zero - should check
    #so norm is PW^2
 #   normPW2 = -1 *biganalytics::apply(PW,2, function(g) sqrt(sum(g^2)))
      
   # x3 = x-UDV$P %*% W
    
   # x3 = dgemm(ALPHA = -1.0, A=UDV$P, B = W, BETA=1.0,C = x) #x - PW
#       normx3 = -1 *biganalytics::apply(gens,2, function(g) sqrt(sum((g-mean(g, na.rm=T))^2),na.rm=T))

   # mean_x = dgemm(A=self$ones_x, B=x3)[1,]
   # mean_x = self$mean_x[[ik]]
   # x4 =t(t(x3[])-mean_x)
  #  x2 = P2 %*% (t(t(x[])-mean_x))
   return(W)  #x-UDV$P %*% W
    
  }else{
    return (NULL)
  }
},
getProjectedData=function(varnames=NULL,
                          vars=lapply(varnames, self$convert)
                          ){
  self$updateUDVP(vars)
  inds = 1:length(self$data)
  names(inds) = names(self$data)
  lapply(inds, function(ik){
      self$projOut(ik)
  })
},
getVariance=function(){
  inds = 1:length(self$data)
  names(inds) = names(self$data)
  lapply(inds, function(ik){
    if(typeof(self$data[[ik]])=="S4"){
    return(sparse_variance(self$data[[ik]]))
    }else{
      return(biganalytics::apply(self$data[[ik]], 2,var, na.rm=T))
    }
  })
},
getProjectedData1=function(prev_var, b_i){
  self$updateUDVP(prev_var)
  UDV=self$UDVP
  D_all = self$extractData(c(prev_var,list(b_i)), adjust=T)
  a= UDV$P %*% UDV$VDU # %*% D_all
  d2=D_all[,ncol(D_all)] - a%*% D_all[,ncol(D_all)]
  d2
},
projOut=function(ik){ 
  UDV=self$UDVP
  x = self$data[[ik]] #[d$nonNA,,drop=F]
  mean_x=self$mean_x[[ik]]
  Dall = t(t(x)-mean_x)
  if(length(UDV$var)>0){
    alias = UDV$alias
    U = UDV$U
    Dinv = UDV$Dinv
    Vinv = UDV$Vinv
    W = Vinv %*% Dinv %*% t(U) %*% Dall[alias,,drop=F]
    x2 = Dall - UDV$P %*% W
    return(x2)
  }else{
    return(Dall)
  }
},
saveParquet=function(){
  mem_dirs = self$mem_dirs
  if(length(mem_dirs)==0) return(NULL)
  nme_d = names(self$norms_list)
  lapply(nme_d, function(nme1){
    print(nme1)
    xx = mem_dirs[[nme1]]
    if(length(self$norms_list[[nme1]])>0){
      files = grep(".pq",dir(xx, full=T),v=T)
      out_F= paste(xx, paste0(length(files)+1,".gz.pq"),sep="/")
      dat2=as_tibble(self$norms_list[[nme1]])
   # print(names(dat2))
     write_parquet(dat2, out_F, compression = "gzip", compression_level = 5)
    }
  })
  return(NULL)
},
getNorm=function(W,var,ik,type){
  if( length(var)==0){
    return(self$norm[[ik]])
  }
 
  #var_st= paste(unlist(lapply(var, paste,collapse=",")),collapse=":")
  var_list1= lapply(var, function(v1){
    paste(names(self$data)[v1[1]], dimnames(self$data[[v1[1]]])[[2]][ v1[2]], sep=",")
  })

  var_st = paste(sort(unlist(var_list1)),collapse=":")
  ##this is where we store norms to avoid calculating multiple times
  ## another option could be in a sqldb
#  print(paste(ik,"HH", var_st,length(self$norms_list)))
  #print(length(self$norms_list))
  #print(length(self$norms_list)[[ik]])
 
  if(ik>length(self$norms_list)) stop(paste("problem", ik, length(self$norms_list)))
  
  norm = NULL
  if(!is.null(self$norm_done)){
    norm = self$norm_done[[names(self$data)[ik]]][[var_st]]
    if(!is.null(norm)){
   #   print(paste("reusing prev",var_st))
      return(norm)
    }
  }
  sn1 = self$norms_list[[ik]]
   if(length(sn1)>0 && getOption("store.norm",TRUE)){
    norm = sn1[[var_st]]  #if precalculated
    if(!is.null(norm)) {
    #  print(paste("reusing",var_st))
      return(norm)
    }
   }

  if(type =="slow"){
    x = W
    warning(" this happens if we use projOut")
    #x2 = t(t(x[])-mean_x)
    ##mean_x already subtracted
    norm=-1*apply(x[]^2,2,sum,na.rm=T)^.5
    #          product = abs(yTr%*%x[,self$cols_incl[[ik]]])
  }else if(TRUE){ ## working with W instead of x
    mean_x=self$mean_x[[ik]] 
    P = self$UDVP$P  
    x_ = self$data[[ik]]#note x = x_ - P%*%W  Also note P%*W columns have mean zero
    ncolw = ncol(W)
    nrowx = nrow(x_)
    #if(TRUE){ ## calculates as a difference
      stepsize=getOption("fspls.stepsize",10000)
      steps=seq(1,ncolw,stepsize)
      #print(paste("norm ", length(steps)))
      norm2=lapply(steps, function(ii){
        ii1 = min(ii+stepsize-1, ncolw)
        apply(x_[,ii:ii1]- P %*% W[,ii:ii1],2,function(g) sum((g-mean(g))^2))
      }) ##, mc.cores = getOption("mc.cores.calc_norm",1))
      norm = -1*sqrt(unlist(norm2))
     # print("done calc")
    #}
    #else{  ##this very slow using hadamard product
#      norm12 = self$norm[[i]]^2
#      norm2 = sapply(1:ncolw, function(ii){
#        g = x_[,ii]
 #       pw = P%*% W[,ii]
#        norm12[[ii]]+sum(pw^2) - 2*sum(g*pw)
#      })
 #     norm = -1*sqrt(norm2)
  #    print("done calc")
  }else{  ## not currently used
    x3 =  self$data[[ik]]-  self$UDVP$P %*%W
      norm = -1*biganalytics::apply(x3,2,function(g) sqrt(sum((g-mean(g))^2)))
  }
 # print(paste("updating", length(self$norms_list[[ik]])))
  if(getOption("store.norm",TRUE)){
    self$norms_list[[ik]][[var_st]]=norm
  }
  return(norm)
},
getAngleInnerOld=function(phensi,ik,k,var, type="slow",var_thresh=1e-5){
  
  assoc=(type %in% c("assoc","assoc1"))
  W = if(type %in% c("slow","assoc"))self$projOut(ik) else self$projOut1(ik)
  ##NOTE PROJOUT1 ALSO SUBTRACTS MEAN, BUT FOR PROJOUT WE HAVE TO ADJUST FOR MEAN
  ##IF WE USE PROJOUT1 then x is actually W
  #mean_x = if(type %in% c("slow","assoc")) NULL else self$mean_x[[ik]]#  x_s$mean_x
  yTr =self$train[[k]]$yTr  #this is zero in the NA positions of d,ie yTr[,d$nonNA]  should be all zero
  ## should we store this as transpose??
 
  #y = self$y
 #d = self$train[[k]]
 # nonNA = lapply(d$nonNA,t)
#  yTr[,!nonNA]=0
     if(assoc){
       stop("rethink this")
       if(!is.null(self$dataNA)) stop(" should have NA in matrix for this")
       #const_term = 0 ##NEED TO THINK ABOUT THIS
       offset = NULL #rep(const_term,nrow(y))
       angle = matrix(1e9, nrow=nrow(yTr),ncol=ncol(x))
       #x1
       #ll0 = fitModel(y[,1],NULL,offset)#,indices)
       for(kk in 1:ncol(y)){
         ##this assumes NAs are NAs not in diff matrix
         lls = apply(x[nonNA,self$cols_incl[[ik]]],2,function(v){
           nonNAv = !is.na(v)
           fitModel(y[nonNA,kk][nonNAv],v[nonNAv],offset)#,indices)
         })
         
         lls = -1*lls
        
         angle[kk,self$cols_incl[[ik]]] = lls
       }
     }else{
        norm = self$getNorm(W,var, ik, type)
        norm_sel = abs(norm[unlist(lapply(var, function(v) if(v[1]==ik) v[2] else NULL))])
        if(length(which(norm_sel>var_thresh))>0) warning(" not projecting out properly .. possibly due to correlated vars !")
        to_rem = which(norm>-var_thresh)
        P = self$UDVP$P
        angles1 = lapply(names(phensi), function(ii){
         #phi = phensi[[nme_i]]
          nme_i = ii
          phensi1 = phensi[[ii]]
       #  lapply(phi, function(ii){
          
          products =self$train[[k]]$product(ik,ii,phensi1);
          nmes_products = names(products); names(nmes_products) = nmes_products
          angle=lapply(nmes_products, function(nmes_prod){
            product = products[[nmes_prod]]
            yTr1 = yTr[[nme_i]][[nmes_prod]]
            if(length(var)==0){
             
              #  self$train$products[[ik]][[ii]][phensi1,,drop=F]  #[,self$cols_incl[[ik]],drop=F]
            }else{
              #PW = P %*% W
              #diff = dgemm(A=yTr1,B=PW)
              PY = yTr1[phensi1,,drop=F] %*% P
              #diff1 = PY %*% W
              product=product-  PY %*% W  #[,self$cols_incl[[ik]],drop=F]
             # dimnames(product) = dimnames(self$train$products[[ik]][[ii]])
            }
            angle_1= t(abs(product[]))/(norm)
            if(nrow(angle_1) !=nrow(yTr1)) angle_1 = t(angle_1)
            if(length(to_rem)>0){
                angle_1[,to_rem]=999  #after we project out the projected out columns have zero norm
            }
            angle_1
          })
        #if(FALSE){
          ##JUST TO SHOW THAT SUBTRACTING MEAN DOESNT MAKE DIFF BECAUSE yTR has mean zero
         # x1 = t(t(x[]) - self$mean_x[[ik]])
        #  product1 = abs(yTr1%*%x1[,self$cols_incl[[ik]]])
        #  product2 = abs(yTr1%*%x1[,self$cols_incl[[ik]]])
        #}
        #the product should not include information from the NA see comment above
        
        ##note yTr has zeros where nonNA is false
        #print(table(yTr[1,!nonNA]))
        ##NOW CALCULATE NORM,WHICH WE NOT ADJUSTING FOR NA
      #  print(dim(product))
      #  print(length(self$cols_incl[[ik]]))
        
        # angle = matrix(999,nrow=nrow(yTr1),ncol = ncol(W))  ## we could try not to calculate this each time
       #   if(length(which(!self$cols_incl))>0){
      #      angle[,!self$cols_incl[[ik]]] = 999 #product2#  [self$cols_incl[[ik]]]
        #  }
       
         ##because norm is actually negative number
        

         angle
        })
     }
  names(angles1) = names(phensi)
     return(angles1)
},
getAngleInner=function(phensi,ik,k,var){
  stop("problems with this")
  d = self$train[[k]]
  UDV=self$UDVP
  var = UDV$var
  dnorm =self$norm[[ik]]  ##USING THE SAME NORM EVEN AFTER PROJECTION, NEED TO FIX
  x=self$data[[ik]]
  ymod = self$ymod[[ik]] #matrix(0, nrow=ncol(y), ncol = nrow(x))
  ymod[,!d$nonNA]=rep(0, nrow(ymod)) ## set the NA values to zero so not included in matrix mult
  #UDV=UDVP_h[[i]]
  
  if(length(var)>0){
    ymod[,d$nonNA] = d$ymod - (d$ymod %*% UDV$P) %*% UDV$VDU
  }else{
    ymod[,d$nonNA] = d$ymod
  }
  #angle[[ik]]
  angle1=if(isbigmatrix(x)) dgemm( A=ymod, B=x) else ymod %*% x
  # if(run_separate){  ##this is true
  ##this does the angle normalisation within each sample.  Whereas using sqnorm does it overall.Within each sample is more stable
  for(ki in 1:nrow(angle1)){
    angle1[ki,] = abs(angle1[ki,])/dnorm
    angle1[ki,which(!self$cols_incl[[ik]])]=999
    angle1[ki,which(abs(dnorm)<1e-5)]=999
  }
  angle1
  #}
},
phensi=function(subphens){
  phensi=lapply(names(subphens), function(nme){
    mi2=match(subphens[[nme]], colnames(self$y[[nme]]))
    names(mi2) = subphens[[nme]]
    mi2
  })
  names(phensi) = names(subphens)
  phensi
},
##incl is which of the layers to include
getAngles1=function(subphens,varnames,incl=names(self$data), k=1,type="slow1"){
  var=lapply(varnames, self$convert)
  var = var[unlist(lapply(var, length))>0]
  phensi = self$phensi(subphens)
  names(phensi)=names(subphens)
  self$getAngles(phensi,var ,incl=incl, k=k, type=type)
},
  getAngles=function(phensi, var,incl=names(self$data), k=1,type="slow1"){  ## type can be fast, slow, assoc, slow1 .. fast gives wrong results
    run_separate=T
    self$updateUDVP(var)
    ntrans=attr(self$data,"ntrans") 
    if(is.null(ntrans)) ntrans=1
    #d = self$train
    nmesd = names(self$norm)
    angles_d= vector('list', length(self$data)) ## need angle object
    names(angles_d) = names(self$data)
    for(ik in 1:length(self$norm)){
    #  print(ik)
      if(names(self$norm)[[ik]] %in% incl){
        if(type=="fast"){
          stop("not working")
          angles_d[[ik]] = self$getAngleInner(phensi,ik,k,var)
        }else{
          angles_d[[ik]] =  tryCatch({
            self$getAngleInnerOld(phensi,ik, k,var,type)
          },error=function(errw){
            print(w)
            NULL
          })
            
        }
      }
    }    
    # angle[which(dnorm==0)]=999
    #  if(FALSE){
    ##SO LONG AS YMOD IS MEAN ZERO, COLUNS OF X DONT HAVE TO BE
    #    if(abs(mean(ymod[1,]))>1e-10) stop("!!")
    #    ones = as.matrix(rep(1, nrow(x)))
    #    means=t(as.matrix(d$mean_x[todoInds]))
    #    x2=ones %*% means
    #    aa1 = (ymod %*% ones) %*% means
    # if (max(abs(diff))>1e-5) stop("!!")     
    #  }
    return(angles_d)
    
  },

storeWeights = function(flags, remote=""){
  weights = self$getWeights()
  jsonw = lapply(weights, function(w) lapply(w,toJSONM))
  df1=.merge1(lapply(jsonw, function(json){
      data.frame(types = toJSON(self$types), pheno = toJSON(dimnames(self$y)[[2]]), 
                 family = toJSON(self$family),
                 offset= json$offset, weights = json$weights, date = date(), flags = toJSON(flags), ip=remote)
  }), addName="name")
  df1 = df1 %>% tibble::add_column(size= unlist(lapply(df1$name, function(x) length(strsplit(as.character(x),",")[[1]]))))
  try(dbWriteTable(dist$mydb, "weights", df1,overwrite=F,append=T)
)
},




translate=function(prev1){
  prev3=lapply(prev1, function(pr){
   pr$clone()
  })
  for(kk in 1:length(prev3)){
    prev3[[kk]]$var=lapply(prev1[[kk]]$varnames, function(vn){
      vns=strsplit(vn,"\\.")[[1]]
      
      i1 = which(names(self$norm)==paste(vns[-length(vns)],collapse="."))
      if(length(i1)==0){
        i1 = which(names(self$norm)==vns[1])
        i2 = which(names(self$norm[[i1]])==paste(vns[-1], collapse="."))
      }else{
        i2 = which(names(self$norm[[i1]])==vns[[length(vns)]]) #paste(vns[-1], collapse="."))
      }
      if(length(i2)==0) stop ("could not translate")
      vv=c(i1,i2)
        ## print(vv)
      vv
    })
   
  }
  prev3
},



updateYpredsAll=function( phensi, inverse_func_str,
                          ypred = self$ypreds_all,
                           prevsk = lapply(1:self$nreps(), function(k) self$train[[k]]$prev),
                           nme = names(prevsk[[length(prevsk)]])
                           ){ ## within=T  means predict on samples trained on
  transform_func = lapply(inverse_func_str, function(str1) eval(str2lang(str1)))  ## should be inverse
  
  for(k in 1:(self$nreps()-1)){
     ind_1=     self$train[[k]]$nonNA
     #d = self$train[[k]]
     inv_func = transform_funct[[prevsk[[k]]$transf]]
     self$updateYP(phensi,prevsk[[k]], ypred, ind_1, inv_func=inv_func, TRUE)  
  }
},

reorder=function(o,k){
  self$train[[k]]$reorder(o)
},
getNA=function(var){
  if(length(var)==0) return(c())
li1 = lapply(var, function(vars1){
   if(length(vars1)<2) return(NULL)
    if(is.null(self$dataNA[[vars1[1]]])) rep(F, self$nrow)  else  (self$dataNA[[vars1[1]]][, vars1[2]] )
  })
df1 = data.frame(li1[unlist(lapply(li1, length))>0])
  apply(df1,1,max)==1
},
 getNonNA=function(var){
    if(length(var)==0) return(rep(T, self$nrow))
    df1 = data.frame(lapply(var, function(vars1){
      if(is.null(self$dataNA[[vars1[1]]])) rep(T, self$nrow)  else  !(self$dataNA[[vars1[1]]][, vars1[2]] )
    }))
    apply(df1,1,min)==1
  },
  extractVar=function(varn){
    lapply(varn, function(vn){
      deets=strsplit(vn,"\\.")[[1]]
      ind1=which(names(self$data)==deets[[1]])
      ind2 = which(dimnames(self$data[[ind1]])[[2]]==deets[2])
      c(ind1,ind2)
    })
  },

## gets ready for training - updates train, prev looc
updateTrain=function(phens,flags, transform_y, verbose=F){ ## this updates the reps and train  ## called after updateLOOC
  nrep = ncol(self$looc$incl)
  if(verbose) print("update train")
  incls = fromJSON(.readFlag(flags,'data_types',"{}")) 
  if(length(incls) == 0 )incls = list("all"=names(self$data))
  incls_all = unique(unlist(incls))
  if(length(self$train)!=nrep  ||  is.null(self$train[[1]]) || toJSON(self$train[[1]]$func_str)!=toJSON(transform_y)){
    self$train = lapply(1:nrep, function(k)trainObj$new( self$y,self$looc , incls_all, transform_y,family=self$family)) #lapply(1:ncol,function(k)
  }
  if(!is.null(self$subset)){
    ## apply subset via the looc object to avoid subsetting big matrix
    self$looc$incl[!self$subset,] = rep(FALSE, ncol(self$looc$incl))
  }
  within=T
  for(k in 1:length(self$train)){
    if(verbose) print(paste("update",k))
    self$train[[k]]$update(self,k,phens)
  }
  reweight=.readFlag(flags,"reweight",FALSE)
  if(reweight){
    self$updateWeights(phens)
  }
},
updateLOOC=function(phens,flags,varn=c(), force=F, verbose=F){
  seed=.readFlag(flags,"seed",42)
  #incl = incls_all
  batch=.readFlag(flags, "batch",0)
  nrep = .readFlag(flags,"nrep",if(batch>0) 0 else 1)
  nrows=nrow(self$data[[1]])
  if(!force && !is.null(self$looc) && nrep ==self$looc$nrep && batch == self$looc$batch && seed == self$looc$seed && nrows ==self$looc$nrows){
    if(verbose)print("no need to update, although we should probably check if the phens changed")
    return(NULL)
  }
  pheno_balance = if(.readFlag(flags,"pheno_balance",FALSE)) unlist(phens) else NULL
  self$looc=loocObj$new(self, nrep = nrep, batch = batch,
                        incl_full = T,
                        seed = seed,
                        nrow = nrows,
                        pheno_balance = pheno_balance)
  if(nrep==1){
    #should really do this inside looc obj
    self$looc$incl = self$looc$incl[,2,drop=F]
  }
  phensi = self$phensi(phens)
  self$prev = list()
  nrep = ncol(self$looc$incl)
  var = self$extractVar(varn)
  W_all = if(length(var)==0) matrix(nrow=0, ncol=0) else  .calcWall_2(self, var)# lapply(datas, function(x) return(matrix(nrow=0, ncol=0)))
  for(k in 1:nrep){
    #  print(k)
    #    phensi,data, betas_new, constants_proj,  k, #mean_y,
    self$prev[[k]] = stateObj$new(phensi, self,NULL,NULL,NULL,k, var=var, varnames=varn, W_all = W_all)
  }
 
  self$train = NULL #lapply(1:nrep, function(k) return(NULL)) #lapply(1:ncol,function(k)
  
},
randomise=function(n= nrow(self$y[[1]]),
                   indices=sample.int(n,n)){
 self$y = lapply(self$y, function(y1){
   y2 = y1[indices,,drop=F]
   fact = attr(y1,"factor")
   if(!is.null(fact)){
     attr(y2,"factor")=fact[indices]
   }
   y2
 })
  if(!is.null(self$train)){
    self$initTrain();  
  }
},
  
cols_incl =function(var_threshs, incl = names(self$norm),g_incl = NULL,qq=1){ #incl, g_incl,qq
  norm = self$norm
 names(incl) = incl
  lapply(incl, function(norm_nme){
     varthreshs = var_threshs[[norm_nme]]
      norm1 = norm[[norm_nme]]
      varthresh = var_threshs[[norm_nme]][qq]
      variance = self$vars[[norm_nme]]
     # varthresh = var_thresh[[norm_nme]]
      var_res = variance>=varthresh
      if(!is.null(g_incl) && length(g_incl)>0 && g_incl!="all"){
        var_res = var_res & (names(norm1) %in% g_incl)
      }
      var_res
  })
},
  updateY=function(y1,    family=NULL,CHECK=T, all_v_all=F, one_v_rest=F){ ## updates y
    if(is.null(rownames(y1))){
      if(nrow(y1)==nrow(self$data[[1]])){
        rownames(y1) = rownames(self$data[[1]])
      }
    }
    mi1 = match(dimnames(self$data[[1]])[[1]],dimnames(y1)[[1]])
    matching_vals = dimnames(self$data[[1]])[[1]][!is.na(mi1)]
    missing_vals = dimnames(self$data[[1]])[[1]][is.na(mi1)]
    if(length(matching_vals)==0) stop("could not match pheno with matrix")
    family=if(!is.null(family)) family else .inferFamily(y1)
    subinds = unlist(lapply(family, length))>0
  #  family = family]
    fams = unique(family[subinds])
    names(fams) = fams
    y = lapply(fams, function(fam){
      #print(fam)
      inds = which(family==fam)
      names(inds) = dimnames(y1)[[2]][inds]
      if(fam=="multinomial"){
        if(all_v_all){
          return(.expandAllvAll(y1[mi1,inds,drop=F]))
        }else if(one_v_rest){
          return(.oneVRest(y1[mi1,inds,drop=F]))
        }
        
        return(.expandFactors(y1[mi1,inds,drop=F]))
      }else if(fam=="ordinal"){
        l1 = lapply(inds, function(ind){
        #  print(ind)
          tryCatch({
          Matrix(as.matrix(y1[mi1,ind,drop=F]))
          },error=function(ew){
            return(NULL)
          })
        })
        l1 [unlist(lapply(l1, length))>0]
        
      }else{
        if(typeof(y1)=="S4"){
          return (list(y1[mi1,inds,drop=F]))
         # return(list(getSparseSubMatrix(getSparseSubMatrix(y1, inds, by="col"), mi1, by="row")))
        }
        if(fam!="gaussian"){
       for(jj in inds){
         
          
         levs = levels(y1[,jj])
         if(is.null(levs)){
           levs = sort(unique(y1[,jj]))
           y1[,jj] = factor(y1[,jj],levels=levs)
         }
           y1[,jj] = as.numeric(y1[,jj])-1
       }
       }
        
        return(list(Matrix(as.matrix(y1[mi1,inds,drop=F]), sparse=T)))
      }
    })
    if(all_v_all) names(y) = gsub("multinomial","binomial.multiway",names(y))
    if(one_v_rest) names(y) = gsub("multinomial","binomial",names(y))
   self$family = unlist(lapply(names(y), function(nme) rep(gsub(".multiway","",nme),
                                                           if(is.list(y[[nme]])) length(y[[nme]]) else ncol(y[[nme]]))),rec=F)
   self$y = unlist(y, rec=F)
   self$weights = rep(1, length(mi1))
    return(NULL)
   # list(missing=missing_vals, matching=matching_vals)
   #self$weights = weights
  },
  initialize_dist=function(dist,types){
    self$subset=NULL
    mem_dirp = getOption("fspls.mem_dir", NULL)
    if(!is.null(mem_dirp)){
      nme_data = names(data)
      names(nme_data) = nme_data
      mem_dirs = lapply(nme_data,function(xx){
        d1 = paste(mem_dirp,nme,xx,sep="/") 
        dir.create(d1, rec=T)
        d1
      })
      self$norm_done = lapply(mem_dirs,function(xx){
        open_dataset(xx)%>% collect()
      })
      self$mem_dirs = mem_dirs
    }else{
      self$mem_dirs = NULL
      self$norm_done = NULL
    }
    
    self$types = types
    self$dist = dist
    types_m = subset(dist$types, type %in%types)
    matrices_all = lapply(types_m$rds, function(rds){
      fi = paste0("^",rev(strsplit(rds,"./")[[1]])[1],".")
      fi1= grep(fi,dir(dist$dir),v=T)
      names(fi1) = unlist(lapply(fi1, function(fi2){
        rev(strsplit(fi2,"\\.")[[1]])[1]
      }))
      matrices =lapply(fi1,function(f1){
        readRDS(paste(dist$dir,f1,sep="/"))
      })
      matrices
    })
    names(matrices_all) = types_m$type
    matricesNA = unlist(lapply(matrices_all,function(m1){
      nonNA = which(names(m1)!="NA")
      na_r = vector("list",length(nonNA))
      names(na_r) = names(m1)[nonNA]
      for(j in 1:length(na_r))na_r[[j]] = m1[['NA']]
      na_r
    }), recursive=F)
    matrices = unlist(lapply(matrices_all,function(m1){
      m1[names(m1)!="NA"]
    }), recursive=F)
    self$data = matrices
    self$dataNA = matricesNA  
   self$nrow = nrow(self$data[[1]])
    self$train=NULL
    nrowd =  nrow(self$data[[1]])
    ym = array(1,dim = c(1, nrowd))
    self$mean_x = lapply(self$data, function(d1){
      mx = attr(d1,"mean_x")
      if(is.null(mx)){
        mx =dgemm( A=ym, B=d1)[1,]/sum(ym)
      }
      mx
    })
    nmes1 = names(self$data)
    self$norm = lapply(self$data, function(d1) {
      normm = attr(d1,"norm")
      if(is.null(normm)){
        normm=   -1 *biganalytics::apply(d1,2, function(g) sqrt(sum((g-mean(g))^2)))
        
      }
      #      print(normm)
      normm
      
    })
    self$norms_list = lapply(self$norm, function(n) list())
    
  },
  initialize=function( cohort, 
                      incl_full=T,seed = 42, memDir = NULL) { ## mem_dirp is for saving scores
   
    if(!is.null(memDir)){
      mem_dirp = memDir
      mem_dirp = paste(memDir,nrow(cohort$rna$matrix),sep="/")  #getOption("fspls.mem_dir", NULL)
      dir.create(mem_dirp, recursive = T);
     nme_data = names(data)
     names(nme_data) = nme_data
     mem_dirs = lapply(nme_data,function(xx){
       d1 = paste(mem_dirp,xx,sep="/") 
       dir.create(d1, rec=T)
       d1
     })
     self$norm_done = lapply(mem_dirs,function(xx){
       open_dataset(xx)%>% collect()
     })
     self$mem_dirs = mem_dirs
    }else{
      self$mem_dirs = NULL
      self$norm_done = NULL
    }
    self$prev=list()
#    self$family=if(!is.null(family)) family else .inferFamily(y)
  self$looc=NULL
    self$subset = NULL
    #nrowy = nrow(data[[1]])
    #ncoly = 1
#    self$ones =     matrix(1, nrow=nrowy, ncol = ncoly)
#    self$ones_x =     matrix(1, ncol=nrowy, nrow= 1)
    self$data = lapply(cohort, function(c) {
      cm = c$matrix
      if(is.null(rownames(cm))) rownames(cm) = 1:nrow(cm)
      cm
    })
    self$dataNA = lapply(cohort, function(c) {
      cm = c$matrixNA
      if(is.null(rownames(cm))) rownames(cm) = 1:nrow(cm)
      cm
      })

    self$nrow = nrow(self$data[[1]])
    self$train=NULL
   
     
    #self$ypreds_all= ypredObj$new(self,NULL, maxn,family)
    nrowd =  nrow(self$data[[1]])
    ym = array(1,dim = c(1, nrowd))
    self$mean_x = lapply(self$data, function(d1){
      mx = attr(d1,"mean_x")
      if(is.null(mx)){
        mx =if(typeof(d1)=="S4")  (ym %*% d1)[1,]/sum(ym) else    dgemm( A=ym, B=d1)[1,]/sum(ym)
      }
      mx
      })
    
  
    nmes1 = names(self$data)
    self$norm = lapply(self$data, function(d1) {
      normm = attr(d1,"norm")
      if(is.null(normm)){
        if(typeof(d1)=="S4"){ ## this not efficient
          means = rowMeans(d1)
          
         # normm=   -1 *apply(d1,2, function(g) sqrt(sum((g-mean(g))^2)))
          normm = -1*sparse_norm(d1)
        }else{
        normm=   -1 *biganalytics::apply(d1,2, function(g) sqrt(sum((g-mean(g))^2)))
        }
      }
#      print(normm)
      normm
      
      })
    self$norms_list = lapply(self$norm, function(n) list())
    nmes1 = names(self$data); names(nmes1) = nmes1
    norm=self$norm
    n = nrow(self$data[[1]])
    
    vars = lapply(names(norm), function(norm_nme){
      norm1 = norm[[norm_nme]]
      (norm1*norm1)/(n-1)
    })
    names(vars) = names(norm)
    
    
    
    self$vars = vars
    self$UDVP=NULL;
  }
)
)


