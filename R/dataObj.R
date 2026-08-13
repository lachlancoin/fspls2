
.better1<-function(rmsv, rmsv2){
  measure = rmsv$measure[[1]]
  subpheno = rmsv$subpheno[[1]]
  cv = rmsv$cv[[1]]
  rmsv_ = subset(rmsv, submeasure=="mid"& cv==cv & measure==measure & subpheno == subpheno) |> pivot_wider(names_from="pheno", values_from="value")
  rmsv2_ = subset(rmsv2, submeasure=="mid"& cv==cv & measure==measure & subpheno == subpheno) |> pivot_wider(names_from="pheno", values_from="value")
a1 = rbind(rmsv_[1,-(1:6)], rmsv2_[1,-(1:6)])

a2 = cbind(c("bef","aft"),a1)
names(a2)[1] = "nme"
a2
}

default_types=jsonlite::fromJSON('{"gaussian": "correlation","binomial" : "AUC","multinomial" : "AUC","ordinal" :"AUC"}')
.convertToTransform<-function(transform_x, fromJSON=FALSE){
  str1 = if(fromJSON) jsonlite::fromJSON(transform_x) else transform_x 
lapply(str1, function(t_y){
  t_y1 =  lapply(t_y[1:2], function(t_y1) eval(str2lang(t_y1)))
  t_y1$params = t_y[[3]]
  .checkInverse(t_y1)
  t_y1 
})
}





## assumes Wall1 is upper diagonal
## uses data with mean subtracted
## does not add any constant term
.eval1_<-function(x_,  Wall2, transf, params, family, CHECK=FALSE){
  if(ncol(Wall2) > 1){
    if( max(abs(Wall2[lower.tri(Wall2)]))>0) stop("!!")
  }
  #if(is.null(dim(beta_new2))) beta_new2 =as.matrix(beta_new2, nrow = ncol(x_), ncol = 1)
 # print(dim(x_)); print(dim(Wall2));
  if(ncol(x_)==0) return(x_);
  t2 = data.frame(vapply(1:ncol(x_), function(jj){
    t1 = transf[[jj]]$func(x_[,1:jj] %*% Wall2[1:jj,jj,drop=FALSE], transf[[jj]]$param) 
   as.vector(t1)
  }, rep(0, nrow(x_))))
  names(t2) = names(transf)
  as.matrix(t2)
  ##t2 %*% beta_new2
}
.getZetaBinary<-function(x,y, w){
  levsy =sort(unique(y[!is.na(y)]))
  yn = .mkBinary(y)
  m1=glm(yn~as.matrix(x),  family="binomial", weights=w)
  ll1 = logLik(m1)
  ll2 =  logLik(update(m1, ~1))
  pv = pchisq((2*(ll1 - ll2)),attr(ll1,"df")[1]-attr(ll2,"df")[1],lower.tail=FALSE,log.p=FALSE)
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

.mkBinary<-function(y, thresh=mean(y,na.rm=TRUE)){
  yn = rep(NA, length(y))
  yn[y<=thresh]=0
  yn[y>thresh]=1
  yn
}


#pv1 = .calcPvalue(x_trans[nonNA ,,drop=FALSE]%*% betas[[kk]],yp1[,kk,drop=FALSE],y, w, family)

.calcPvalue<-function(yp_new,yp1k,y,  w, family){   ## this seems to not work anymore for multinomial
  if(length(which(!is.na(y)))==0) return (0)
  #if(family=="multinomial"){
  #  yp_new = .eval1_(x_, Wall2, transf, family)
  #}else{
  #  yp_new = .eval1_(x_, Wall2, transf, family) %*% beta_new2
  #}
  if(family=="multinomial"){
    m1 = multinom_ridge(yp1k,y,w)
    m2 = multinom_ridge(yp_new,y,w)
    ll2 = -0.5 * m2$dev  
    ll1 = -0.5 *m1$dev       
    df2  = (ncol(yp_new)-1) * (length(levels(y))-1)
    df1  = (ncol(yp1k)-1) * (length(levels(y))-1)
    
    pv1 = .lrt(ll2,ll1,df2, df1, log.p=TRUE)
  }else if(family=="ordinal"){
    pv1 =tryCatch({
      df1 = data.frame(cbind(y,as.matrix(yp1k )));
      df1$y=factor(df1$y, levels = sort(unique(df1$y)))  
      func = paste0("y~",paste(colnames(df1)[-1], collapse="+"))
      m1=suppressWarnings(polr(func,  data=df1,weights=w,Hess=TRUE, method="logistic"))
      df1[,2] =  yp_new[,1]
      m2=suppressWarnings(polr(func,  data=df1,weights=w,Hess=TRUE, method="logistic"))
      ll1 = logLik(m1)
      ll2 =  logLik(m2)
      .lrt(ll2,ll1,2,1,log.p=TRUE)  
    },error=function(ew){
      m1=glm(y~yp1k[,1],weights=w,family="gaussian")
      m2=glm(y~yp_new[,1], weights=w,family="gaussian")
      ll2 = logLik(m2)
      ll1 =  logLik(m1)
      .lrt(ll2,ll1,2,1, log.p=TRUE)
    })
  }else{
    if(ncol(yp_new)>1 || ncol(yp1k)>1) stop("not expecting")
   
    m1=suppressWarnings(glm(y~yp1k[,1],family=family,weights=w)) ##, weights should be integer
    m2=suppressWarnings(glm(y~yp_new[,1], family=family, weights = w) )
    ll2 = logLik(m2)
    ll1 =  logLik(m1)
    if(ll1==0)warning("problem, zero likelihood")
    pv1 = .lrt(ll2,ll1,2,1, log.p=TRUE)
  }
  as.vector(pv1)
}

multinom_ridge<-function(x,y,w,lambda=getOption("lambda")){
  if(is.null(dim(x)) || ncol(x)==1){
    one_inds = 1:2
    x = cbind(1,x)
  }else{
    one_inds = 1
  }
#  m1=(multinom(y~x,weights=w, trace=FALSE))
#  sm1  = summary(m1)
  
  
  ridge<- tryCatch(
    glmnet(x,y, family="multinomial", alpha=0, weights=w, lambda=lambda, standardize=FALSE), #, lambda.min.ratio = 1e-10),
    error = function(e) {
      if (grepl("zero variance", conditionMessage(e))) {
        NULL
      } else {
        stop(e)
      }
    }
  )
    
    if(is.null(ridge)) return(ridge);
  
  mins = min(ridge$lambda)
#predy=predict(ridge,x,s=min(ridge$lambda))[,,1]
  rbeta <- coef(ridge,s=min(ridge$lambda))
  rdf = as.matrix(data.frame(lapply(rbeta, as.matrix)))
  dimnames(rdf)[[2]]= names(rbeta)
  devs = deviance(ridge)
  list(const = apply(rdf[one_inds,,drop=FALSE],2,sum), beta = rdf[-one_inds,],dev= devs[length(devs)])
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
 Matrix(as.matrix(res),sparse=TRUE)
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
    }),recursive=FALSE)
  }))
  rownames(res) = rownames(y)
  m1 = Matrix(as.matrix(res),sparse=TRUE)
}

.expandFactors<-function(y, max_cats=50){
  res = lapply(1:ncol(y), function(i){
  #  print(i)
    y1 = y[,i]
    if(!is.factor(y1)) y1 = factor(y1, levels = sort(unique(y1)))
    levs1 = levels(y1)
    if(length(levs1) > max_cats) {
      warning(paste("more than max_cats", names(y)[[i]],length(levs1), max_cats))
      return(NULL)
    }
    mat = Matrix(0, nrow = length(y1), ncol = length(levs1), dimnames = list(rownames(y), levs1),sparse=TRUE)
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

.sumMatrices<-function(matrices, onlyAll=FALSE){
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
  

.combineAngles<-function(anglesH, incl,onlyAll=FALSE, topn=100, excl = list()){
  names(incl)=incl
  angles = lapply(anglesH, function(aa) aa$angles)
  cols_incl = lapply(anglesH, function(aa) aa$cols_incl)
  
  nme_trans = names(angles[[1]][[1]][[1]]); names(nme_trans) = nme_trans
  comb_all1=lapply(nme_trans, function(nme_t1){
  #nme_t1 = nme_trans[[1]]
     nme_pow = names(angles[[1]][[1]][[1]][[nme_t1]]); names(nme_pow)=nme_pow
     lapply(nme_pow, function(nme_p1){
      comb_all = lapply(incl,function(inc1){
        matrices = lapply(1:length(angles), function(i){
          ang1 = angles[[i]][[inc1]]
          if(is.null(ang1)) return(NULL)
          col_incl = cols_incl[[i]][[inc1]]
          #ang1=angle1[[inc1]]
          ang2=ang1[[1]][[nme_t1]][[nme_p1]]
          cs = Matrix::colSums(ang2)
          if(length(ang1)>1){
            for(jk in 1:length(ang1)){
              cs = cs+Matrix::colSums(ang1[[jk]][[nme_t1]][[nme_p1]])
            }
          }
          excl1 = excl[unlist(lapply(excl, function(ex) ex[3]==nme_t1 && ex[1] == inc1 && ex[4] ==nme_p1))]
          if(length(excl1)>0){
          col_incl[which(names(col_incl) %in% unlist(lapply(excl1, function(ex) ex[2])))]=FALSE
          }
          cs[col_incl]
        })
        names(matrices) = names(angles)
        matrices[unlist(lapply(matrices, length))>0]
        if(length(matrices)==0) return(NULL)
        .sumMatrices(matrices, onlyAll=onlyAll)
      })
 
  top_angles=whichpart1(comb_all, n=topn, return_scores=TRUE)
  t1 = .merge1_new(lapply(top_angles, function(ta){
    data.frame(list(names = names(ta), value=ta))
  }),addName="data_type")
  t1 = t1[order(t1$value),]
  subset(t1, value<999)
  })
  })
  comb_all1
}


.check<-function(var1,var2){
  if(length(var1)!=length(var2)) stop("not match var")
  if(length(var1)>0){
    if(length(which(unlist(var1)!=unlist(var2)))>0) stop("problem with check !!") ## should also check contenst
  }
  return(TRUE)
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

.getWeights11_1<-function(prev,const=FALSE, pvs=FALSE){
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
      
      nconst=unlist(lapply(const_proj1, length),recursive=TRUE)
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







dataObj<-R6::R6Class("dataObj", public = list(
  dist="environment",
  data="list",
  db_name="character",
  dataNA="list",
  types="vector",
  nrow="integer",
  vars="list",
  direction="character",
#  var_thresh="double",
transforms = "list",  ## this is the functions
#transform_x="list", ## this is the string
#default_transform="closure",
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
  if(length(b_i1)==0) return(NULL)
  data_ind = which(names(self$data)==b_i1[[1]])
  if(length(data_ind)==0){
   # print(names(self$data))
    stop(paste("could not find",b_i1[[1]]))
  }
  var_ind = which(dimnames(self$data[[data_ind]])[[2]]==b_i1[2])
  if(length(var_ind)==0) return(NULL)
  trans_ind=match(b_i1[[3]],names(self$transforms))
  pow_ind = match(b_i1[[4]], names(self$transforms[[trans_ind]]$params))
  c(data_ind, var_ind, trans_ind, pow_ind)
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
dims=function(){
  data = self
  res = list(
    nonNA=lapply(data$data, function(d){
      list(dim = dim(d),bigmatrix=is.big.matrix(d))
    }),
    na=lapply(data$dataNA, function(d){
      list(dim = dim(d),bigmatrix=is.big.matrix(d))
    }))
  res
},
pheno = function(maxpheno=1e9,sep=FALSE, sep_group=FALSE,code = NULL,memb=NULL){ 
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
 if(!is.null(memb)){
   levs = sort(unique(memb));names(levs)=levs
   resu1 = lapply(levs, function(lev){
     m1 = names(memb[memb==lev])
     li = lapply(self$y, function(y1){
       colnames(y1)[which(colnames(y1) %in% m1)]
     })
     li[unlist(lapply(li, length))>0]
   })
   return(resu1)
 }
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
   return( c(l12,(unlist(l1,recursive=FALSE))))
  }else{
    return(list(all=phens))
  }
},
####does regression just on orthogonal component
calcBetaProj1=function(subphens,k,b_i,b_i_name, prev_var, Wall,convert=TRUE, betas = list(), strict=FALSE,project=FALSE, 
                                             useglm=getOption("useglmnet",TRUE), useoffset=FALSE){
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
    #print(nme)
    phensi1 =  match(nme, names(self$y))
    family=getOption("fspls.family",self$family[phensi1])
    phensi_=subphens[[nme]]
    betas1 = betas[[nme]]
    Wall1 = Wall[[nme]]
    if(family=="multinomial") useoffset=FALSE
    if(length(prev_var)>0 ){
      if(is.null(betas1) ) stop("no betas");
      res1=self$calcBetaProjAll(nme,phensi_,family,  k, b_i, b_i_name,prev_var, Wall1,betas1,project=project, strict=strict,useglm=useglm, useoffset=useoffset)
    }else{
      res1 = self$calcBetaProj(nme,phensi_,family,  k, b_i, b_i_name, prev_var, Wall1,strict=strict, useglm=useglm)
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
getTransforms=function(vars1){
  default_transform = eval(str2lang("function(x,pow) x"))
  nme1 ="func"; ## always using func since we use invfunc in trainObj  if(inv_transform)"invfunc" else "func"
  deft = list(func=default_transform, param=1)
  transf = list(deft)
  if(length(vars1)>0 && !is.null(vars1[[1]])){
  #  print("h")
   # print(vars1)
    transf = lapply(vars1, function(v1) {
      tf = self$transforms[[v1[[3]]]];
      param = tf$params[[v1[[4]]]]
    
      list(func=tf$func, param =param)
    })
  }
 
  transf
},
calcBetaProj=function(nme,phensi_,family, k,b_i,b_i_name, prev_var,Wall, strict=FALSE, 
                    
                      useglm=getOption("useglmnet",TRUE)){
  #b_i = b_i1
  if(length(prev_var)>0) stop("problem")
  data = self
  train = data$train
  d = train

  # print(family)fcalz
  nonNA = self$looc$incl[,k]
  vars1 = c(prev_var,list(b_i))
  names(vars1)[length(vars1)]  = paste(b_i_name,collapse=".")
  
  data$updateUDVP(prev_var)
  transf = self$getTransforms(vars1)
  ys =self$y[[nme]]##  if(family %in% c("binomial","ordinal")) (d$y1) else if(family =="multinomial") d$y2 else d$y
  if(family=="multinomial"){
    ys = data.frame(list(attr(ys,"factor")))
    names(ys) = nme
  }else if(length(phensi_) < ncol(ys)){
    ys = ys[,phensi_,drop=FALSE] 
  }
  ncoly = ncol(ys)
  betas = apply(ys,2,function(v1) list())# apply(train$y,2, function(v1) list())
  tbls =  apply(ys,2, function(v1) list())
    constants = apply(ys,2, function(v1) list())
  
  pvs =  apply(ys,2, function(v1) list())
    j = length(vars1)
    if(length(vars1[[j]])==0){
   # stop(" this doesnt happen")
    non_na_x = rep(TRUE, self$nrow)
    useglm=FALSE
  }else{
    non_na_x = if(is.null(data$dataNA[[vars1[[j]][1]]])) rep(TRUE,nrow(ys) ) else !(data$dataNA[[vars1[[j]][1]]][, vars1[[j]][2]] )
  }
  if(length(which(non_na_x))==0){
  #  print(vars1);
    stop("no non_zero");
  }
  if(length(vars1[[j]])==0){
    #stop(" this doesnt happen")
    x_=matrix(1, nrow=self$nrow, ncol=1)
    meanx = 1
    }else{
      x_ =self$extractData(vars1, adjust=FALSE, convert=FALSE)
      meanx = apply(x_,2,mean)
      
      
    }
  UDV = data$UDVP ## should check it corresponds to prev_i
  Wall1 =self$UDVP$getWall(x_[,ncol(x_)]-meanx[ncol(x_)], Wall)
  
  varx = var(x_[non_na_x],na.rm=TRUE)
  if(is.na(varx) ){
    stop("variance  NA")
  }
  #if(varx<1e-10)stop("variance  too small")
  #k =NULL ##delete this later
 # transform_func_y = self$default_transform
 # if(length(b_i)>0 && !inv_transform) {
#    transform_func_y= self$transforms[[b_i[[3]]]][[1]]
#  }
 # print(transform_func_y)
#  .print_verbose("transform_func_y",transform_func_y,1)
  x_trans  =  .eval1_(x_,  Wall1, transf, family)
  for(kk in 1:ncoly){
    nonNAk = nonNA & non_na_x
    nonNAk1 = non_na_x[nonNA]; ## second subset
    x = x_trans[nonNAk,];##transf[[1]]$func(x_[nonNA], transf[[1]]$param)[nonNAk1];
    
   # x = transf[[1]]$func(x_, transf[[1]]$param)[nonNAk ]
    y = ys[,kk][nonNAk] ##transform_func_y(ys[,kk])[nonNAk]
    w = data$weights[nonNAk]
    beta_new1=0;
    #x = x_t[nonNAk]
    if(family=="multinomial"){
      ty=as.list(table(y))
      tbls[[kk]] = ty[ty>0]
      
     # rdf = try(multinom_ridge(x,y,w))
      
      rdf=multinom_ridge(x, y, w)
      if (is.null(rdf)) {
        levsy = levels(y)
        beta_new1 = rep(0, length(levsy)); names(beta_new1) = levsy
        const_term = beta_new1
      }else{
      #logpv1 = pchisq(rdf$dev/2,df=1,low=FALSE,log=TRUE)
        beta_new1=rdf$beta
        const_term = rdf$const
      }
    }else if(family=="ordinal" ){
      ty=as.list(table(y))
      tbls[[kk]] = ty[ty>0]
      use_bin = length(unique(y,na.rm=TRUE))<=2
      if(!use_bin){
        df = data.frame(cbind(y,x ))
        df$y = factor(y, levels = sort(unique(y,na.rm=TRUE)))  
        # print("using polr") 
        m1=try(suppressWarnings(polr(y~x,  data=df,weights=w,Hess=TRUE, method="logistic")),silent=TRUE)
        #predict(m1, df, type = "p")
        
        if(inherits(m1,"try-error")){
          if(strict) warning("error with fitting polr, refitting as many binary")
          resz = .getZetaBinary(x,y, w)
          beta_new1=resz$beta_new1
          const_term =resz$const_term
          pv1 = resz$pv1
          #const_term = -1*sm$coefficients[1,1]
        }else{
          ll1 = logLik(m1)
          ll2 =  logLik(update(m1, ~1))
          pv = pchisq((2*(ll1 - ll2)),attr(ll1,"df")[1]-attr(ll2,"df")[1],lower.tail=FALSE,log.p=FALSE)
          coeff = c(m1$coefficients,1,1,pv[1])
          const_term = m1$zeta
          beta_new1=coeff[1]
          pv1 = coeff[4] # pt(abs(coeff[3]),df=1,lower.tail=FALSE)  ## not sure this correct
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
     # spike_slab_iter = getOption("spike_slab_iter",0)
      nonNAy = !is.na(y)
      if(length(which(nonNAy))==0){
        const_term =0; beta_new1 = 0
      }else{
        if(useglm){
          sm2<- tryCatch({
        
          ones = rep(1, length(x))
          x_mod = cbind(ones,x)
       
          ridge=glmnet(x_mod[nonNAy,,drop=FALSE],y[nonNAy],weights = w[nonNAy], family=family, alpha = 0,
                       standardize=FALSE,
                       lambda = getOption("lambda"), lambda.min.ratio = 1e-10)
          
          rbeta <- coef(ridge,s=min(ridge$lambda))
          const_term = sum(rbeta[1:2,1])
          beta_new1 = rbeta[3,1]
          if(FALSE){
            ypred = predict(ridge,x_mod, s=min(ridge$lambda))
            ll1 = .logLik(y,ypred, family=family)
            ll2 =.logLik(y,family=family) 
            pv1 = .lrt(ll1,ll2,2,1, log.p=FALSE)
          }else{
            nulldev=ridge$nulldev
            dev= ridge$dev.ratio[which.min(ridge$lambda)] *nulldev
            #       pv1 = pchisq(-1*(dev-nulldev),df=1,low=FALSE)
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
          withCallingHandlers({
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
            pv1 = .lrt(ll1,ll2,2,1, log.p=FALSE)
          }
          # print(paste(pv1_1, pv1, family))
          # print(pv1)
          #pchisq((2*(ll1 - ll2)),attr(ll1,"df")[1]-attr(ll2,"df")[1],lower.tail=FALSE,log.p=FALSE)
          list(const_term = const_term, beta_new1 = beta_new1, pv1=pv1)
          }, warning = function(w) {
            if(grepl("fitted probabilities numerically 0 or 1", conditionMessage(w))){
              invokeRestart("muffleWarning")   # ignore this one, keep going
            }
            # anything else falls through and propagates up to tryCatch's handler
          })
        }, warning=function(errw) {
         # print(errw)
       #   message("using glmnet 810")
          ones = rep(1, length(x))
          x_mod = cbind(ones,x)
          nonNAy = !is.na(y)
          ridge=glmnet(x_mod[nonNAy,,drop=FALSE],y[nonNAy],family=family,weights=w[nonNAy], alpha = 0, 
                       standardize=FALSE,
                       lambda = getOption("lambda"), lambda.min.ratio = 1e-10)
          #ridge=glmnet(x_mod,y,family=family, alpha = 0)
          rbeta <- coef(ridge,s=min(ridge$lambda))
          const_term = sum(rbeta[1:2,1])
          beta_new1 = rbeta[3,1]
          if(FALSE){
            ypred = predict(ridge,x_mod, s=min(ridge$lambda))
            ll1 = .logLik(y,ypred, family=family)
            ll2 =.logLik(y,family=family) 
            pv1 = .lrt(ll1,ll2,2,1, log.p=FALSE)
          }else{
            nulldev=ridge$nulldev
            dev= ridge$dev.ratio[which.min(ridge$lambda)] *nulldev
            #       pv1 = pchisq(-1*(dev-nulldev),df=1,low=FALSE)
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
      pv1 = .lrt(ll1,ll2,length(levels(y))-1,1, log.p=TRUE)
    }else{
      yp_new=x_trans[nonNA ,,drop=FALSE]%*% beta_new1
      pv1 = .calcPvalue(yp_new,matrix(rep(1, length(which(nonNA)))),ys[nonNA,kk], data$weights[nonNA], family)
      
#       pv1 = .calcPvalue(matrix(x_[nonNA, drop=FALSE]), ys[nonNA,kk], beta_new1, matrix(rep(1, length(which(nonNA)))),
 #                        data$weights[nonNA], family,Wall1, transf)
    }
    pvs[[kk]] = pv1
    constants[[kk]] = const_term  #-mean_adj*beta_new1
  }
  list(betas=betas, constants = constants,pvs = pvs,tbls=tbls,Wall1=Wall1)
},

testBetaProj=function(vars_all, transform_x){
  self$updateTransform(transform_x);
  b_i_name = vars_all[[1]]$variables[[1]][[1]]
  b_i = self$convert(b_i_name)
  nme = names(self$y)[[1]]
  phensi_ = colnames(self$y[[1]])
  family=nme;
  k = 1
  prev_var = list(); 
  vars1 = c(prev_var,list(b_i))
  names(vars1)[length(vars1)]  = paste(b_i_name,collapse=".")
  self$updateUDVP(prev_var)
  transf = self$getTransforms(vars1)
  Wall =diag(0); 
  useglm = TRUE; strict=FALSE;
  ys =self$y[[nme]];
  y = ys[,1]
  non_na_x = rep(TRUE, self$nrow)
  x_ =self$extractData(vars1, adjust=FALSE, convert=FALSE)
  meanx = apply(x_,2,mean)
  params = seq(1, 200)/100; names(params) = params
  ab =  lapply(params, function(p){
    x1 = transf[[1]]$func(x_,p);
    s = summary(  glm(y ~x1[,1], family=nme))
    s$coeff[2,4]
  })
  ac=unlist(lapply(ab, log))
  plot(params, ac)
},
#getWall=function(prev_var,x_){
#  Wall1 = diag(1)
 # for(kkj in 1:length(prev_var)){
#    UDVP1=UDVPObj$new(self, prev_var[kkj]) #,check=FALSE, centralise=FALSE)
 #   Wall1 = UDVP1$getWall(x_[,kkj+1], Wall1)
#  }
#},

##project variable controls whether the projection of the variable is fitted, or the variable itself
## Wall1 is projection from previous
calcBetaProjAll=function(nme,phensi_,family, k,b_i,b_i_name, prev_var, Wall1,betas1, project=FALSE, strict=FALSE, 
                         CHECK=getOption("fspls.check",TRUE),
                                              useoffset =FALSE,useglm=getOption("useglmnet",TRUE)){
  data = self 
  if(!useoffset) project=FALSE
  nonNA = self$looc$incl[,k]
  vars1 = c(prev_var,list(b_i)) 
  names(vars1)[length(vars1)]  = paste(b_i_name,collapse=".")
  data$updateUDVP(prev_var)
 
  ys =self$y[[nme]]##  if(family %in% c("binomial","ordinal")) (d$y1) else if(family =="multinomial") d$y2 else d$y
  if(family=="multinomial"){
    ys = data.frame(list(attr(ys,"factor")))
    names(ys) = nme
  }else if(length(phensi_) < ncol(ys)){
    ys = ys[,phensi_,drop=FALSE] 
  }
  ncoly = ncol(ys)
  betas = apply(ys,2,function(v1) list())# apply(train$y,2, function(v1) list())
  tbls =  apply(ys,2, function(v1) list())
  constants = apply(ys,2, function(v1) list())
  pvs =  apply(ys,2, function(v1) list())
  j = length(vars1)
  non_na_x = if(is.null(data$dataNA[[vars1[[j]][1]]])) rep(TRUE,nrow(ys) ) else !(data$dataNA[[vars1[[j]][1]]][, vars1[[j]][2]] )
  x_ =self$extractData(vars1, adjust=FALSE)
  meanx = t(replicate(nrow(x_), apply(x_,2,mean,na.rm=TRUE)))
  
  transf = self$getTransforms(vars1) #lapply(vars1, function(v1) self$transforms[[v1[[3]]]][[2]])
  UDV = data$UDVP ## should check it corresponds to prev_i
  Wall2 = UDV$getWall(x_[,ncol(x_)]-meanx[,ncol(x_)], Wall1) ## updated projection  may not need to subtract meanx  .. leaving it for now
  if(!project) Wall2 = diag(ncol(Wall2))
  if(CHECK){
    x1_ = x_[,ncol(x_), drop=FALSE]  
 
  if(project){ ## PROBABLY DONT NEED THIS ANYMORE SINCE USING WALL
      #####  x1_ = (x_ %*% Wall2-x_)[,ncol(x_), drop=FALSE] # should be zero for all except last column 
    
         x1_ = x1_ - UDV$P %*% (UDV$VDU %*% x1_)   
       
        if(CHECK){ ##THIS xDEMONSTRATE ORTHOGONALITY
             d3 = self$extractData(UDV$var, adjust=TRUE)
          #   var_x1 = var(x1_[,1])
              x5=unlist(lapply( 1:length(UDV$var), function(jk){
                
                x2_ =d3[,jk] # self$data[[UDV$var[[jk]][1]]][,UDV$var[[jk]][2]]- self$mean__x(UDV$var[[jk]]) # [[UDV$var[[jk]][1]]][UDV$var[[jk]][2]]
                x3=x2_%*%x1_[,1]
                x3[1,1]
              }))
              names(x5) = dimnames(x_)[[2]][-ncol(x_)]
              inds111=which(abs(x5)>1e-5 )#*sd(x1_[,1]))
            if( length(inds111)){
              warning(paste("no longer orthogonal but might be due to transformation giving large SD", max(abs(x5)), dimnames(x_)[[2]][ncol(x_)], " vs " ,paste(names(inds111),collase=",")), sd(x1_[,1], sd(d3[,inds111[1]])))
            }
          x3=x_ %*%Wall2 #- meanx %*% Wall2
          #plot(x3[,ncol(x3)], x1_[,1])
          chck=sum(abs(x3[,ncol(x3)]-mean(x3[,ncol(x3)])-(x1_[,1] - mean(x1_[,1]))))
          if(chck>1e-5){
            warning(paste("problem with Wall",chck))
          }
         # cbind(x_[,1:ncol(x_)-1], x1)
        }
  }
  
  if(!project){
    Wall1 = diag(1, ncol(Wall1))
  }
  }
  if(length(which(duplicated(colnames(x_)))>0)){
    
    stop("problem .. duplicated colnames")
  }
  non_na_x = apply(x_,1,function(v) length(which(is.na(v))))==0
  #transform_func_y = self$default_transform
  #if(!inv_transform){
  #    transform_func_y =self$transforms[[b_i[[3]]]][[1]]
  #}
 # 
   x_trans  =  .eval1_(x_,  Wall2, transf, family)
   x_trans_prev  =  .eval1_(x_[,-ncol(x_),drop=F],  Wall1, transf[-ncol(x_)], family)
   
#x_trans=  t(t(x_trans) -apply(x_trans,2,mean)  )
#   sum(abs(transf1$func(x1_, transf1$param[1])-x_trans[,ncol(x_trans)] - ( mean(x1_[,1]) - mean(x_trans[,ncol(x_trans)]))))
  
   if(family=="multinomial"){
  yp1 = x_trans_prev[nonNA,,drop=FALSE]
  }else{
    yp1 =   x_trans_prev[nonNA,,drop=FALSE] %*% betas1
  }
  
  for(kk in 1:ncoly){
    nonNAk = nonNA & non_na_x
    nonNAk1 = non_na_x[nonNA]; ## second subset
    y=ys[,kk][nonNAk]
    w = data$weights[nonNAk]
    beta_new1=0;
    const_term=0
   # Wall1 = NULL
   
    #  if(family=="multinomial")  (x_1 %*%  Wall1) %*% betas1  else (x_[,-ncol(x_),drop=FALSE ] %*%  Wall1) %*% betas1 [,kk,drop=FALSE]
    
   if(useoffset){
      x = cbind(yp1[,kk], x_trans[nonNA,ncol(x_trans)])
      dimnames(x)[[2]] = c(paste0("A",1:(ncol(x)-1)),"x")
    }else{
     # warning("no transformation here without offset!")
      x1_ = x_[,ncol(x_), drop=FALSE]  
      transf1 = transf[[ncol(x_)]]
      x =cbind(yp1,transf1$func( x1_[nonNAk,,drop=FALSE], transf1$param))
     
     # for(kk3 in 1:ncol(x)){
    #    transf1$func(x[,], transf1$param)
    #  }
    }
   # yp1 = yp1[nonNAk1,kk,drop=FALSE]  
    if(family=="multinomial"){
       ty=as.list(table(y))
      tbls[[kk]] = ty[ty>0]
      #mean_x = apply(x,2, mean,na.rm=T)
      #x = t(t(x) - mean_x)
      #print(apply(x,2,mean)) 
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
      
        func = paste0("y~",paste(colnames(x), collapse="+"))
        m1=try(suppressWarnings(polr(func,  data=df,weights=w,Hess=TRUE, method="logistic")),silent=TRUE)
        #predict(m1, df, type = "p")
        
        if(inherits(m1,"try-error")){
        #  if(strict) warning("polr!!!")
          resz = .getZetaBinary(x,y, w)  ##need as.matrix(x) ??
          beta_new1=resz$beta_new1
          const_term =resz$const_term
  
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
      #spike_slab_iter = getOption("spike_slab_iter",0)
      if(length(which(nonNAy))==0){
        const_term=0; beta_new1 = rep(0,ncol(x))
      }else {
       
        if(useglm){
          #if(CHECK && FALSE){
          #  print(summary(glm(y~x[,1])))
          #  ridge=glmnet(cbind(1,x[nonNAy,1,drop=FALSE]),y[nonNAy],family=family,weights=w[nonNAy], alpha = 0)
          #  rbeta <- coef(ridge,s=min(ridge$lambda))
          #  print(rbeta)
          #}
          
          sm2<- tryCatch({
            #if(ncol(x)==1){
           # ones = rep(1, nrow(x))
          #  x_mod = cbind(ones,x)
            #}else{
            #  x_mod = x
            #}
                ridge=glmnet(x[nonNAy,,drop=FALSE],y[nonNAy],family=family,weights=w[nonNAy], alpha = 0,
                             standardize=FALSE,
                             lambda = getOption("lambda"), lambda.min.ratio = 1e-10)
                rbeta <- coef(ridge,s=min(ridge$lambda))
                const_term = rbeta[1,1]
                 beta_new1 =  rbeta[-1,1]
         
            list(const_term = const_term, beta_new1 = beta_new1)
          }, error=function(errw) {
            warning("reverting to glmnet")
                m1=glm(y~as.matrix(x), family=family, weights=w, standardize=FALSE) ## including weights lead to non-convergence
            sm  = summary(m1)
        
            if(nrow(sm$coeff)<1+ncol(x)){
              coeff = rep(0,4)
              const_term=0
              beta_new1 = rep(0, ncol(x))
            }else{
#              coeff = sm$coeff[2,]
              
              const_term = sm$coefficients[1,1]
              beta_new1=sm$coeff[-1,2]
              #if(length(beta_new1)!=ncol(x)){
              #  print(beta_new1)
              #  print(colnames(x))
              #}
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
            withCallingHandlers({
             m1=glm(y~as.matrix(x), family=family, weights=w) ## including weights lead to non-convergence
            sm  = summary(m1)
            #print(var(x))
            if(nrow(sm$coeff)<1+ncol(x)){
              coeff = rep(0,4)
              const_term=0
              beta_new1 = rep(0, ncol(x))
            }else{
              #coeff = sm$coeff[2,]
              
              const_term = sm$coefficients[1,1]
              beta_new1=sm$coefficients[-1,1]
           
              names(beta_new1)=colnames(x)
            }
            
            list(const_term = const_term, beta_new1 = beta_new1)
            }, warning = function(w) {
              if(grepl("fitted probabilities numerically 0 or 1", conditionMessage(w))){
                invokeRestart("muffleWarning")   # ignore this one, keep going
              }
              # anything else falls through and propagates up to tryCatch's handler
            })
          }, warning=function(errw) {
         #   print(errw)
          #  warning("using glmnet 1132")
            #ones = rep(1, nrow(x))
            #x_mod = cbind(ones,x)
            nonNAy = !is.na(y)
           
              ridge=glmnet(x[nonNAy,,drop=FALSE],y[nonNAy],family=family, alpha = 0, weights=w[nonNAy],lambda = getOption("lambda"), 
                           standardize=FALSE,
                           lambda.min.ratio = 1e-10)
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
     
      if(useoffset){ #project){
       ##betas[[kk]] = if(family=="multinomial") Wall1 %*% beta_new1 else   (Wall1 %*% beta_new1)[,1]
        betas[[kk]] = rbind(betas1[,kk,drop=FALSE] * beta_new1[1], beta_new1[2])
      }else{
        betas[[kk]] = beta_new1
      }
     # print("H")
    #  print(dim(x_));
    #  print(nonNA);
    #  print(dim(y))
    #  print(dim(yp1));
    if(family=="multinomial"){
      yp_new=x_trans[nonNA ,,drop=FALSE]
      pv1 = .calcPvalue(yp_new,yp1,y, w, family)
      
    }else{
      yp_new=x_trans[nonNA ,,drop=FALSE]%*% betas[[kk]]
    
       pv1 = .calcPvalue(yp_new,yp1[,kk,drop=FALSE],y, w, family)
    }
    pvs[[kk]] = pv1
      constants[[kk]] = const_term  #-mean_adj*beta_new1
  }
  list(betas=betas, constants = constants,tbls = tbls, pvs = pvs, Wall = Wall2)
},
##var and Wall1 are from one smaller model
#calcWall=function(b_i, var, Wall1, inv_transform=getOption("x_transform",TRUE)){
#  data = self
#  data$updateUDVP(var, inv_transform=inv_transform)
#  UDV = data$UDVP ## should check it corresponds to prev_i
#  x = data$data[[b_i[1]]][, b_i[2]] - data$mean_x[[b_i[1]]][b_i[2]]
#  UDV$getWall(x, Wall1)
#},
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
mean__x = function(v1){
  if(length(v1)==0) return(c())
#  self$mean_x[[v1[[1]]]][v1[[2]]]
  mean(self$data[[v1[1]]][,v1[2]])
},
extractData=function(var, adjust=TRUE,convert=FALSE){
  if(convert)var = lapply(var,self$convert)
  Dall = Matrix( 0,nrow = self$nrow, ncol = length(var), dimnames = list(dimnames(self$data[[1]])[[1]], names(var)), sparse=TRUE)
  nme = rep("", length(var))
  if(length(var)==0)return(Dall)
  for(jk in 1:length(var)){
    v1 = var[[jk]]
    Dall[,jk] = self$data[[v1[1]]][,v1[2]];
    t_ind = v1[[3]]
    nme[[jk]] = paste(names(self$data)[[v1[1]]], dimnames(self$data[[v1[1]]])[[2]][v1[2]])
  
  
    #print(c( mean(Dall[,jk]),self$mean_x[[v1[1]]][v1[2]]))
    if(adjust)Dall[,jk] = Dall[,jk]- mean(Dall[,jk]) #self$mean_x[[v1[1]]][v1[2]]
    #if(inv_transform) {
    #  Dall[,jk] = self$transforms[[t_ind]][[2]](Dall[,jk]) ## applies inverse transform
    #  nme[[jk]] = paste(nme[[jk]], names(transforms)[[t_ind]])
    #}
  }
  dimnames(Dall)[[2]] = nme
  Dall
},

#vars2 = vars_all[[1]]$variables[[1]];phens1 = phens[[1]]; k=1; useoffset=TRUE
makeModels=function(phens1, vars2,k, 
                    project=TRUE,logpthresh = -5,useglm=TRUE,useoffset=TRUE,
                    flags = list(),checkRMSV=TRUE, verbose=FALSE
                  
                  ){
    nonNA = self$looc$incl[,k]
  if(is.null(self$looc)){
    self$updateLOOC( phens1,flags,varn=c(),force=FALSE, verbose=FALSE)
  }
  phen2 = phens1
  ypred=self$ypred(phen2)
  phensi = self$phensi(phen2)
  subphens = phensi
  nmes = c()
  data = self;
  len = length(vars2)
  models = vector("list", len)
  useglm=getOption("useglmnet",TRUE)
  fams1 = lapply(names(phens1), function(st1)strsplit(st1,"\\.")[[1]][1])
  family = getOption("fspls.family",fams1) #ypred$family[[1]]
  if(family[[1]]=="multinomial") {
    useoffset=FALSE
    project=FALSE
  }
  nme1 = ""
  b_i_name=c()
  names(family)=family
  #transform_x1=transform_x[[var_transf[[1]]]]
 # nme_c1 = var_transf[[1]]
  Wall =lapply(subphens, function(f) matrix(nrow=0, ncol=0))
  prev_i = self$makeNextModel(NULL,b_i_name,subphens,k, 
                              family, ypred=ypred, project=project, useglm=FALSE, 
                      Wall=Wall,
                              logpthresh =logpthresh,useoffset=useoffset)
  models[[1]] = prev_i$simplify(useoffset)
  nmes[[1]] = "empty"
  jk=1

 #Wall = prev_i$Wall
 #rmsv=NULL
 rmsv=if(checkRMSV) self$checkRMSV(subphens,prev_i, ypred, nonNA,verbose=verbose) else NULL
  #print(rmsv1)
  while(jk<=len){
    b_i_name = vars2[[jk]]
#    transform_func=eval(str2lang(func_str1[var_transf[[jk]]]))
   # transform_x1=transform_x[[var_transf[[jk]]]]
    prev_i1 = self$makeNextModel(prev_i,b_i_name,subphens,k,Wall=Wall,
                                 family, ypred=ypred, project=project, useglm=useglm, logpthresh =logpthresh,useoffset=useoffset)
    if(is.null(prev_i1)) break;
    nme2 = paste(vars2[[jk]], collapse=".")
    nme1 = if(jk==1)  nme2 else paste(nme1, nme2,sep=";")
    nmes[[jk+1]] = nme1
    models[[jk+1]] =prev_i1$simplify(useoffset)
    if(checkRMSV){
      rmsv2=self$checkRMSV(subphens,prev_i1, ypred, nonNA,verbose=verbose)
     # print(subset(rmsv2, submeasure=="mid"))
      if(!is.null(rmsv)){
      better = .better1(rmsv, rmsv2)
       improvement = apply(better[,-1,drop=FALSE],2, function(x)x[2]-x[1])
       #print(better)
       worse = improvement[which(improvement < -1e-3)]
       if(length(worse)>0){
        
         warning(paste("not improving!", toJSON(as.list(sort(worse)))))
       }
       
       
   
      }
      rmsv = rmsv2
    }
   
    prev_i = prev_i1
    Wall = prev_i1$Wall
  
    jk = jk+1
   # print(paste(jk,nrow(Wall[[1]])))
  }
  ## next we need to adjust projection and set offset
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
      
      colind=which(dimnames(y1)[[2]]==names(subphens[[1]])[1]) 
    }
    fact = apply(y1[,colind ,drop=FALSE ],1,paste,collapse="_")
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
checkRMSV=function(subphens, prev_i1, ypred, nonNA,verbose=FALSE, useglm=TRUE){
   data=self
  ypred$updateYP(data, prev_i1, nonNA, flip=FALSE, liab=FALSE)
  rmsv=(ypred$calcRMSV(self$y, nonNA,     flip=FALSE))
   subset( rmsv, submeasure=="mid")
},
 makeNextModel=function(prev_i2, b_i_name, phens, k, family, ypred=NULL, project=TRUE, useglm=TRUE,    logpthresh = -5,
                        useoffset=TRUE,
                        Wall = prev_i2$Wall,
                        CHECK=getOption("fspls.check",FALSE),
                        verbose=getOption("fspls.verbose1",FALSE)) {
    data =self
    b_i = if(length(b_i_name)==0) c() else self$convert(b_i_name)
    prev_var = if(is.null(prev_i2)) list() else lapply(prev_i2$var_names, self$convert)
    self$updateUDVP(prev_var)
    #Wall = data$calcWall(b_i, prev_i$var, prev_i$Wall) ## WALL not important, we can get rid of it later
    #prev_var = if(jk==1) prev_i$var  else lapply(vars2[1:(jk-1)], self$convert)
    betas = prev_i2$betas_proj;
  
    b_new_proj <- withCallingHandlers(
      self$calcBetaProj1(phens, k, b_i, b_i_name, prev_var, Wall,
                         betas = betas, project = project, convert = FALSE,
                         strict = TRUE, useglm = useglm, useoffset = useoffset),
      warning = function(w) {
        if (grepl("fitted probabilities numerically 0 or 1 occurred", conditionMessage(w))) {
          invokeRestart("muffleWarning")
        }
      }
    )  
    if(length(b_i)>0  && sum(unlist(b_new_proj$pvs))==0 && getOption("show_warnings",FALSE)) {
   
      warning(paste("log pvalues should not reach 0 at ", paste(b_i_name, collapse=",")))
    }
    Wall2 = b_new_proj$Wall
   # b_new_proj1 = self$calcBetaProj1(subphens,k,b_i,prev_var,  betas = betas, project=!project,convert=FALSE,    strict=TRUE, useglm=useglm, useoffset=useoffset) 
    
    betas_new = b_new_proj$betas
      constants_proj =if(family[[1]]=="multinomial") b_new_proj$constants[[1]] else b_new_proj$constants
      pvs =if(family[[1]]=="multinomial") b_new_proj$pvs[[1]] else b_new_proj$pvs
   #   print(pvs)
    #  if(.sumChisq(pvs)>logpthresh ){
        
     # return(NULL)
    #  }
    tbls = if(family[[1]]=="multinomial") b_new_proj$tbls[[1]]  else b_new_proj$tbls
    
    if(family[[1]]=="multinomial"){
      betas_proj=lapply(betas_new, function(b_n){
        if(is.matrix(b_n[[1]])) return(b_n[[1]])
        dw = t(as.matrix(data.frame(b_n)))
      })
    }else{
      betas_proj=lapply(betas_new, function(b_n){
        as.matrix(data.frame(b_n))
      })
    }
    mean_x = self$mean__x(b_i) #[[b_i[[1]]]][b_i[2]]
    
    prev_i1=stateObj$new(phens,data, betas_proj,constants_proj, tbls, prev_i2 , b_i,b_i_name=b_i_name, mean_x = mean_x, Wall = Wall2,
                         pvs =pvs, useoffset=useoffset)
  #  prev_i1$setOffset() 
    #if(is.null(prev_i2)) return(prev_i1)
    #prev_i1$setOffset()
    if(FALSE){ ##JUST TO CHECK WHAT GLM VARIABLES  WOULD BE IF WE JUST FIT ALL
      extractd0 = self$extractData(c(prev_i2$var, list(b_i)), adjust=FALSE)
      means_x = apply(extractd0,2,mean)
      extractd = self$extractData(c(prev_i2$var, list(b_i)), adjust=TRUE)
      subnme = names(phens)[[1]]
      family = strsplit(subnme,"\\.")[[1]][1]
      df2 = data.frame(lapply(phens[[1]], function(colind){
        if(family=="multinomial"){
          y_m = attr(self$y[[1]],"factor")
          m1 = try(multinom(y~as.matrix(extractd0),weights=w, trace=FALSE))
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
          y11 = self$y[[subnme]][,phens[[subnme]],drop=FALSE]
          self$updateWeights(phens)
          nonNAy = !is.na(y11[,1])
          ab=glm(y11[,1] ~as.matrix(extractd0))
        ridge=glmnet(cbind(1,extractd0[nonNAy,]),y11[nonNAy,,drop=FALSE] ,family=family, alpha = 0,
                     standardize=FALSE,
                     weights =self$weights[nonNAy],lambda = getOption("lambda"), lambda.min.ratio = 1e-10)
        rbeta <- coef(ridge,s=min(ridge$lambda))
        aa = predict(ridge,cbind(1,extractd0[nonNAy,]),s=min(ridge$lambda), family=family)
        ridge1 = glmnet(cbind(1, extractd0[nonNAy,] %*%  rbeta[-(1:2),,drop=FALSE]), y11[nonNAy,,drop=FALSE], family = family, alpha=0,lambda = getOption("lambda"), 
                        standardize=FALSE,
                        lambda.min.ratio = 1e-10)
        rbeta1 <- coef(ridge1,s=min(ridge1$lambda))
      
        return(cbind(rbeta[-(1:2),],prev_i1$betas[[1]],ab$coefficients))
        }
#        (rbeta[-(1:2),] - prev_i1$betas[[1]][,colind,drop=FALSE])
      }))
     
#      dimnames(df2)[[2]] = 1:ncol(self$y[[1]])
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

plotData=function(vars_all1, phens1 = vars_all1$phens, all_types=FALSE, transform_x = NULL, violin=FALSE, assoc=FALSE){
  phensi = self$phensi(phens1)
  nmei = names(phensi); names(nmei) = nmei
  df = data.frame(lapply(nmei, function(nmei1){
    as.matrix(self$y[[nmei1]][,phensi[[nmei1]],drop=FALSE])
  }))
  
  variables = vars_all1$variables
  names(variables)=NULL
  vars = unlist(variables, recursive=FALSE)
  incls = names(self$data); names(incls)=incls
  
  if(all_types){
    genes = unique(unlist(lapply(vars, function(x) x[[2]])))
    names(genes) = genes
    vars1 = lapply(incls, function(incl){
      lapply(genes, function(g){
       c(incl,g)
      })
    })
    vars = unlist(vars1,recursive=FALSE)
  }
  into=c("data","gene")
  
  if(!is.null(transform_x) && all_types){
    nmev = names(vars)
    transform_x1 = jsonlite::fromJSON(transform_x)
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
        pv1 = .lrt(ll2,ll1,2,1, log.p=FALSE)
        pv1
      })
    })
    pv_res = sort(unlist(pvs))
    
    return(pv_res)
  }
  
  nme_df = names(df); names(nme_df) = nme_df
  df4 = .merge1_new(lapply(nme_df, function(nmedf1){
    df_k = df[[nmedf1]]
    df3=df2 |> tibble::add_column(y=df_k) |> pivot_longer(names(df2)) |> separate("name",sep="__", into=into)
  }),addName="pheno")
  df4$y = factor(df4$y)
  df4 = subset(df4, !is.na(y))
  
  df4
},

ypred=function(phens1){
  family = unlist(lapply(names(phens1), function(str)getOption("fspls.family",strsplit(str,"\\.")[[1]][1])))
  ypr = ypredObj$new(self,self$phensi(phens1),family=family)
  ypr
},# inverse_func_strs = jsonlite::fromJSON(.readFlag(flags,"transform_x_inverse",'{"y":"function(y) y"}'))
#all_models_y = all_models$y; inverse_func_str = jsonlite::fromJSON(flags1$transform_x_inverse)[[1]]; self = datasAll$datas[[1]]

extractPredictions=function(all_models_y,phens, flags, 
                            transform_x = jsonlite::fromJSON(flags$transform_x),
                            ypred = self$ypred(phens), liab=TRUE
                                   ){
  
  self$updateTransforms(transform_x   )      
  d = self
  self$updateLOOC(phens,flags)
    ## whether to evaluate with liability , default is true
  #  ypred = self$ypred(phens)
  minv = .readFlag(flags,"min",0)
  maxv = .readFlag(flags,"max",1e6)
  verbose=.readFlag(flags, "verbose",FALSE)
  group_names= names(all_models_y); names(group_names)=group_names
  numvars = unlist(lapply(group_names, function(x) if(x=="empty") 0 else length(strsplit(x,";")[[1]])))
  numvars1 = sort(unique(numvars))
  names(numvars1) = numvars1
  numvars1 = numvars1[numvars1>=minv & numvars1<=maxv]
  #pheno_nmes = names(phens); names(pheno_nmes)=pheno_nmes
  if(length(all_models_y)==0) return(NULL)
  #
  #
  #nmes_models = names(all_models_y[[1]]);names(nmes_models) = nmes_models;  numvar = numvars1[[3]]; nmes1 = nmes_models[[1]];  
  #group_names2 = group_names[numvars==numvar]; group_name = group_names2[[1]]
  evals_all =lapply(numvars1, function(numvar){
    if(verbose)cat(paste("numvar",numvar))
    #.merge1_new(lapply(pheno_nmes, function(pheno_nme){
    #ypred = ypreds[[pheno_nme]]; 
    if(is.null(ypred)) stop("ypred is null")
    # .merge1_new(lapply(nmes_models, function(nmes1){
    group_names2 = group_names[numvars==numvar]
    # evals = .merge1_new(lapply(group_names2, function(group_name){
    #    if(verbose) cat(paste(numvar,nmes1,group_name))
    all_models1 = all_models_y[names(all_models_y) %in% group_names2]#[[pheno_nme]]   
    all_models2 = lapply(all_models1, function(am) lapply(am, function(am1) am1))
    all_models2_full = all_models2[unlist(lapply(all_models2, function(am)"full" %in% names(am) ))]
    if(length(all_models2)>0){
      all_models2 = all_models2[!unlist(lapply(all_models2, is.null))]
    }
    names(all_models2) = NULL
    names(all_models2_full) = NULL
    all_models3 = unlist(all_models2,recursive=FALSE)
    all_models3_full = unlist(all_models2_full,recursive=FALSE)
    full_model = all_models3[["full"]]
    full_model_nme=paste(names(full_model$var_names), collapse=";")
    nmesm = grep("full",names(all_models3),inv=TRUE,value=TRUE);
    nmesm_full = grep("full",names(all_models3_full),inv=TRUE,value=TRUE);
    inds=as.numeric(nmesm); 
    inds_full = as.numeric(nmesm_full)
    res1 = list(); #res2 = list(); res3 = list()
    if(!is.null(full_model)){
      #ypredObj$updateYP(self, phens, )#= self$looc$incl[,k2]
      nonNA =self$looc$incl[,self$nreps()]
      ypred$updateYP(d, full_model, nonNA, flip=FALSE, liab=liab)
      #res1 = ypred$calcRMSV(self$y, nonNA,      flip=FALSE)
      # print(res1);
      res1 [["full"]] =  ypred$predictions(nonNA, flip=FALSE)
      for(lk in 1:length(res1[['full']])){
        attr(res1[['full']][[lk]],"betas") = full_model$betas[[lk]][numvar,]
      }
      #res1 = self$getRMSVInds(phens, d$nreps(), ypred)  
    }
    if(length(nmesm)>0){
      #transf=c()
      for(j in 1:length(nmesm)){
        nonNA =self$looc$incl[,inds[[j]]]
        prev_i1 = all_models3[[j]]
        #   transf = c(transf,prev_i1$transf)
        ypred$updateYP(d, prev_i1, nonNA, flip=TRUE)
        #          self$updateYpredsInds(phens,all_models1[[j]][[nmes1]], inds[[j]], ypred)
      }
      nonNA=self$getNonNAInds(inds)
      #res2 = ypred$calcRMSV(self$y,nonNA, flip=TRUE)|> tibble::add_column(isfull=FALSE,model="cv")
      res1[["cv"]] =  ypred$predictions(nonNA, flip=TRUE)
    }
    if(length(nmesm_full)>0){
      #transf=c()
      for(j in 1:length(nmesm_full)){
        nonNA =self$looc$incl[,inds_full[[j]]]
        prev_i1 = all_models3_full[[j]]
        #   transf = c(transf,prev_i1$transf)
        ypred$updateYP(d, prev_i1, nonNA, flip=TRUE)
        #          self$updateYpredsInds(phens,all_models1[[j]][[nmes1]], inds[[j]], ypred)
      }
      nonNA=self$getNonNAInds(inds_full)
     # res3 = ypred$calcRMSV(self$y,nonNA, flip=TRUE)|> tibble::add_column(isfull=TRUE,model=full_model_nme)
      res1[['cv_full']] = ypred$predictions(nonNA, flip=TRUE)
    }
    res1
  #  rbind(res1,res2,res3)
    #        }),addName="model")
    #     }),addName="trainedOn")
    # }),addName="pheno_group")
  })
  #if(!is.null(evals_all$numvars)){
  #  evals_all$numvars = as.numeric(evals_all$numvars)
  #}
  evals_all

},
evaluateAllModels=function(all_models_y,phens,flags,
                           ypred = self$ypred(phens), #lapply(phens, function(phens1) self$ypred(phens1)),
                           verbose=FALSE
                         ){

  d = self
  self$updateLOOC(phens,flags)
  liab = .readFlag(flags,"liab",TRUE)  ## whether to evaluate with liability , default is true
#  ypred = self$ypred(phens)
  group_names= names(all_models_y); names(group_names)=group_names
  numvars = unlist(lapply(group_names, function(x) if(x=="empty") 0 else length(strsplit(x,";")[[1]])))
  numvars1 = sort(unique(numvars))
  names(numvars1) = numvars1
  #pheno_nmes = names(phens); names(pheno_nmes)=pheno_nmes
  if(length(all_models_y)==0) return(NULL)
  
  #nmes_models = names(all_models_y[[1]]);names(nmes_models) = nmes_models;  numvar = numvars1[[1]]; nmes1 = nmes_models[[1]];  group_names2 = group_names[numvars==numvar]; group_name = group_names2[[1]]
  evals_all = .merge1_new(lapply(numvars1, function(numvar){
    if(verbose)cat(paste("numvar",numvar))
          if(is.null(ypred)) stop("ypred is null")
            group_names2 = group_names[numvars==numvar]
              all_models1 = all_models_y[names(all_models_y) %in% group_names2]#[[pheno_nme]]   
              all_models2 = lapply(all_models1, function(am) lapply(am, function(am1) am1))
              all_models2_full = all_models2[unlist(lapply(all_models2, function(am)"full" %in% names(am) ))]
              if(length(all_models2)>0){
                all_models2 = all_models2[!unlist(lapply(all_models2, is.null))]
              }
              names(all_models2) = NULL
              names(all_models2_full) = NULL
              all_models3 = unlist(all_models2,recursive=FALSE)
              all_models3_full = unlist(all_models2_full,recursive=FALSE)
              full_model = all_models3[["full"]]
               full_model_nme=paste(names(full_model$var_names), collapse=";")
              nmesm = grep("full",names(all_models3),inv=TRUE,value=TRUE);
              nmesm_full = grep("full",names(all_models3_full),inv=TRUE,value=TRUE);
              inds=as.numeric(nmesm); 
              inds_full = as.numeric(nmesm_full)
              res1 = NULL; res2 = NULL; res3 = NULL
              if(!is.null(full_model)){
                #ypredObj$updateYP(self, phens, )#= self$looc$incl[,k2]
                nonNA =self$looc$incl[,self$nreps()]
                ypred$updateYP(d, full_model, nonNA,  flip=FALSE, liab=liab)
                res1 = ypred$calcRMSV(self$y, nonNA,      flip=FALSE)
               
                beta = unlist(lapply(1:nrow(res1), function(ii){
               
                  b1 = full_model$betas[[res1$family[ii]]]
                  ik2 = which(colnames(b1)==res1$pheno[[ii]])
                  if(length(ik2)==0)ik2 = 1
                  b1[nrow(b1),ik2]
                }))
                
                res1 = res1 |> tibble::add_column(isfull=TRUE, model=full_model_nme, beta, sign = sign(beta))
                #res1 = self$getRMSVInds(phens, d$nreps(), ypred)  
              }
              if(length(nmesm)>0){
                #transf=c()
                for(j in 1:length(nmesm)){
                  nonNA =self$looc$incl[,inds[[j]]]
                  prev_i1 = all_models3[[j]]
               #   transf = c(transf,prev_i1$transf)
                  ypred$updateYP(d, prev_i1, nonNA, flip=TRUE)
                  #          self$updateYpredsInds(phens,all_models1[[j]][[nmes1]], inds[[j]], ypred)
                }
                nonNA=self$getNonNAInds(inds)
              
                res2 = ypred$calcRMSV(self$y,nonNA, flip=TRUE)|> tibble::add_column(isfull=FALSE,model="cv", beta=NA, sign = NA)
              }
              if(length(nmesm_full)>0){
                #transf=c()
                for(j in 1:length(nmesm_full)){
                  nonNA =self$looc$incl[,inds_full[[j]]]
                  prev_i1 = all_models3_full[[j]]
                  #   transf = c(transf,prev_i1$transf)
                  ypred$updateYP(d, prev_i1, nonNA, flip=TRUE)
                  #          self$updateYpredsInds(phens,all_models1[[j]][[nmes1]], inds[[j]], ypred)
                }
                nonNA=self$getNonNAInds(inds_full)
                res3 = ypred$calcRMSV(self$y,nonNA, flip=TRUE)|> tibble::add_column(isfull=TRUE,model=full_model_nme, beta=NA, sign = NA)
              }
              rbind(res1,res2,res3)
#        }),addName="model")
#     }),addName="trainedOn")
   # }),addName="pheno_group")
  }),addName="numvars")
  if(!is.null(evals_all$numvars)){
  evals_all$numvars = as.numeric(evals_all$numvars)
  }
  evals_all
},



  updateModel=function(k,best_all_i, model_prev,to_keep,CHECK=TRUE){
    if(TRUE) stop("not updating..")
    sprev = self$train$prev[[k]]
    best_all=lapply(1:length(best_all_i), function(j){
      prev_i= sprev[[to_keep[j]]]
      best_i = best_all_i[[j]]
      if(CHECK){
        .check(prev_i$var, model_prev[[to_keep[j]]]$var)
      }
      #data_ind = which(names(self$data)==b_i[[1]])
      #b_i1 = c(data_ind, which(dimnames(self$data[[data_ind]])[[2]]==b_i[2]))
      nxt_v= lapply(best_i, function(b_i) stateObj$new(self, self$train[[k]], prev_i,self$convert(b_i1)))
         lapply(nxt_v, function(nv) nv$updateConst(self, self$train[[k]],k))
      #                       nxt_v = nxt_v[!unlist(lapply(nxt_v, is.null))]
      names(nxt_v) = unlist(lapply(nxt_v,  function(nv)paste(unlist(lapply(nv$var, function(vv){
        paste(names(self$data)[vv[1]],dimnames(self$data[[vv[1]]])[[2]][vv[2]],sep=".")
      })), collapse=",")))
      nxt_v
    })
    self$train[[k]]$prev_old = self$train[[k]]$prev
    self$train[[k]]$prev=unlist(best_all, recursive=FALSE)
  },
  getMaxBetaProj=function(){
    #     lapply(self$prevs, function(prevk){
    mabv = lapply(self$prev, function(pk) max(abs(unlist(pk$betas_proj))))
    #   })
    mabv
  },
  updateUDVP=function(var){
    inv_transform=FALSE
    UDV=self$UDVP
    if(!is.null(UDV) && !is.null(UDV$var)){
      if(length(UDV$var)==0 &&length(var)==0) return(NULL)
      if(length(UDV$var)==length(var) &&  length(which(unlist(UDV$var)!=unlist(var)))==0) return(NULL) ## dont need to update
    }
    Dall = self$extractData(var, adjust=TRUE)
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
    #W = (Vinv %*% Dinv) %*% (t(U) %*% Dall[alias,,drop=FALSE])
    #UD = Ut %*% Dall[d$nonNA,]  ##t(U) %*% Dall[alias,,drop=FALSE]
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
#       normx3 = -1 *biganalytics::apply(gens,2, function(g) sqrt(sum((g-mean(g, na.rm=TRUE))^2),na.rm=TRUE))

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
      return(apply(self$data[[ik]], 2,var, na.rm=TRUE))
      
      #return(biganalytics::apply(self$data[[ik]], 2,var, na.rm=TRUE))
    }
  })
},
getProjectedData1=function(prev_var, b_i){
  self$updateUDVP(prev_var)
  UDV=self$UDVP
  D_all = self$extractData(c(prev_var,list(b_i)), adjust=TRUE)
  a= UDV$P %*% UDV$VDU # %*% D_all
  d2=D_all[,ncol(D_all)] - a%*% D_all[,ncol(D_all)]
  d2
},
projOut=function(ik){ 
  #warn("!! mean_x might not be right")
  UDV=self$UDVP
  x = self$data[[ik]] #[d$nonNA,,drop=FALSE]
  #mean_x= self$mean_x[[ik]]
  mean_x = apply(x,2,mean)  ## replacement
  Dall = t(t(x)-mean_x)
  if(length(UDV$var)>0){
    alias = UDV$alias
    U = UDV$U
    Dinv = UDV$Dinv
    Vinv = UDV$Vinv
    W = Vinv %*% Dinv %*% t(U) %*% Dall[alias,,drop=FALSE]
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
    
    xx = mem_dirs[[nme1]]
    if(length(self$norms_list[[nme1]])>0){
      files = grep(".pq",dir(xx, full=TRUE),value=TRUE)
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
    norm=-1*apply(x[]^2,2,sum,na.rm=TRUE)^.5
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
    norm = -1*apply(x3,2,function(g) sqrt(sum((g-mean(g))^2)))
    
#      norm = -1*biganalytics::apply(x3,2,function(g) sqrt(sum((g-mean(g))^2)))
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
      stop("not supported")
      
     }else{
        norm = self$getNorm(W,var, ik, type)
        norm_sel = abs(norm[unlist(lapply(var, function(v) if(v[1]==ik) v[2] else NULL))])
        if(length(which(norm_sel>var_thresh))>0 ) warning(" not projecting out properly .. possibly due to correlated vars !")
        to_rem = which(norm>-var_thresh)
        P = self$UDVP$P
        angles1 = lapply(names(phensi), function(ii){
                 nme_i = ii
          phensi1 = phensi[[ii]]
                 products =self$train[[k]]$product(ik,ii,phensi1);
          nmes_products = names(products); names(nmes_products) = nmes_products
          #nmes_prod = "rand"; nmes_prod1 = names(products[[nmes_prod]])[1]
          angle=lapply(nmes_products, function(nmes_prod){
            product0 = products[[nmes_prod]]
            nmes_products1 = names(product0); names(nmes_products1) = nmes_products1
            lapply(nmes_products1, function(nmes_prod1){
              product = products[[nmes_prod]][[nmes_prod1]]
            yTr1 = yTr[[nme_i]][[nmes_prod]][[nmes_prod1]]
            if(length(var)==0){
             
              #  self$train$products[[ik]][[ii]][phensi1,,drop=FALSE]  #[,self$cols_incl[[ik]],drop=FALSE]
            }else{
              #PW = P %*% W
              #diff = dgemm(A=yTr1,B=PW)
              PY = yTr1[phensi1,,drop=FALSE] %*% P
              #diff1 = PY %*% W
              product=product-  PY %*% W  #[,self$cols_incl[[ik]],drop=FALSE]
             # dimnames(product) = dimnames(self$train$products[[ik]][[ii]])
            }
            direction=self$direction
           # print(paste("dir",direction));
            if(direction=="+"){
              angle_1= t((product[]))/(norm)
              angle_1[angle_1>0]=999
            }else if(direction=="-"){
              angle_1= -1* (t((product[]))/(norm))
              angle_1[angle_1>0]=999
              
            }else{
              angle_1= t(abs(product[]))/(norm)
            }
            if(ncol(angle_1) != length(norm) || colnames(angle_1)[1]!=names(norm)[1]) angle_1 = t(angle_1)
            #if(nrow(angle_1) !=nrow(yTr1)) angle_1 = t(angle_1)
            if(length(to_rem)>0){
                angle_1[,to_rem]=999  #after we project out the projected out columns have zero norm
            }
            angle_1
          })
          })
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
    if(is.null(self$y[[nme]])) stop("subphens is wrong")
    mi2=match(subphens[[nme]], colnames(self$y[[nme]]))
    names(mi2) = subphens[[nme]]
    na_ind= which(is.na(mi2))
    if(length(na_ind)>0) stop("should not have NA")
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
    run_separate=TRUE
    self$updateUDVP(var)
    ntrans=attr(self$data,"ntrans") 
    if(is.null(ntrans)) ntrans=1
    #d = self$train
    nmesd = names(self$norm)
    angles_d= vector('list', length(self$data)) ## need angle object
    names(angles_d) = names(self$data)
    for(ik in 1:length(self$norm)){
    #  (ik)
      if(names(self$norm)[[ik]] %in% incl){
        if(type=="fast"){
          stop("not working")
          angles_d[[ik]] = self$getAngleInner(phensi,ik,k,var)
        }else{
          angles_d[[ik]] =  tryCatch({
            self$getAngleInnerOld(phensi,ik, k,var,type)
          },error=function(errw){
            message(errw)
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
  df1 = df1 |> tibble::add_column(size= unlist(lapply(df1$name, function(x) length(strsplit(as.character(x),",")[[1]]))))
  try(dbWriteTable(dist$mydb, "weights", df1,overwrite=FALSE,append=TRUE)
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





reorder=function(o,k){
  self$train[[k]]$reorder(o)
},
getNA=function(var){
  if(length(var)==0) return(c())
li1 = lapply(var, function(vars1){
   if(length(vars1)<2) return(NULL)
    if(is.null(self$dataNA[[vars1[1]]])) rep(FALSE, self$nrow)  else  (self$dataNA[[vars1[1]]][, vars1[2]] )
  })
df1 = data.frame(li1[unlist(lapply(li1, length))>0])
  apply(df1,1,max)==1
},
 getNonNA=function(var){
    if(length(var)==0) return(rep(TRUE, self$nrow))
    df1 = data.frame(lapply(var, function(vars1){
      if(is.null(self$dataNA[[vars1[1]]])) rep(TRUE, self$nrow)  else  !(self$dataNA[[vars1[1]]][, vars1[2]] )
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
updateTransform=function(transform_x){
  transforms =  .convertToTransform(transform_x)
  update_trans = toJSON(self$transforms)!= toJSON(transforms)
  self$transforms = transforms
  update_trans
},
## gets ready for training - updates train, prev looc
updateTrain=function(phens,flags, transform_x, verbose=FALSE, force=FALSE){ ## this updates the reps and train  ## called after updateLOOC
 update_trans = self$updateTransform(transform_x);
  nrep = ncol(self$looc$incl)
  if(verbose) cat("update train")
  incls = jsonlite::fromJSON(.readFlag(flags,'data_types',"{}")) 
  if(length(incls) == 0 )incls = list("all"=names(self$data))
  incls_all = unique(unlist(incls))
  if(force || length(self$train)!=nrep  ||  is.null(self$train[[1]]) || update_trans){ # || toJSON(self$train[[1]]$func_str)!=toJSON(transform_x)){
    self$train = lapply(1:nrep, function(k)trainObj$new( self$y,self$looc , incls_all, self$transforms,family=self$family)) #lapply(1:ncol,function(k)
  }
  if(!is.null(self$subset)){
    ## apply subset via the looc object to avoid subsetting big matrix
    self$looc$incl[!self$subset,] = rep(FALSE, ncol(self$looc$incl))
  }
  within=TRUE
  for(k in 1:length(self$train)){
    if(verbose) cat(paste("update",k))
    self$train[[k]]$update(self,k,phens, force=force)
  }
  reweight=.readFlag(flags,"reweight",FALSE)
  if(reweight){
    self$updateWeights(phens)
  }
},
updateLOOC=function(phens,flags,varn=c(), force=FALSE, verbose=FALSE){
  seed=.readFlag(flags,"seed",42)
  #incl = incls_all
  batch=.readFlag(flags, "batchsize",0)
  nrep = .readFlag(flags,"nfold",if(batch>0) 0 else 1)
  if(batch>0 && nrep>0) warning("only one of batch or nrep greater than zero")
  nrows=nrow(self$data[[1]])
  if(!force && !is.null(self$looc) && nrep ==self$looc$nrep && batch == self$looc$batch && seed == self$looc$seed && nrows ==self$looc$nrows){
    if(verbose)cat("no need to update, although we should probably check if the phens changed")
    return(NULL)
  }
  pheno_balance = if(.readFlag(flags,"pheno_balance",FALSE)) unlist(phens) else NULL
  self$looc=loocObj$new(self, nrep = nrep, batch = batch,
                        incl_full =TRUE,
                        seed = seed,
                        nrow = nrows,
                        pheno_balance = pheno_balance)
  if(nrep==1){
    #should really do this inside looc obj
    self$looc$incl = self$looc$incl[,2,drop=FALSE]
  }
  phensi = self$phensi(phens)
  self$prev = list()
  nrep = ncol(self$looc$incl)
  var = self$extractVar(varn)
  if(length(var)>0) stop("need to reimplement calcWall_2 for this")# .calcWall_2(self, var)# lapply(datas, function(x) return(matrix(nrow=0, ncol=0)))
  Wall =  matrix(nrow=0, ncol=0) 
  for(k in 1:nrep){
    #  print(k)
    #    phensi,data, betas_new, constants_proj,  k, #mean_y,
    self$prev[[k]] = stateObj$new(phensi, self,NULL,NULL,NULL, var=var, varnames=varn, Wall = Wall)
  }
 
  self$train = NULL #lapply(1:nrep, function(k) return(NULL)) #lapply(1:ncol,function(k)
  
},
randomise=function(n= nrow(self$y[[1]]),
                   indices=sample.int(n,n)){
 self$y = lapply(self$y, function(y1){
   y2 = y1[indices,,drop=FALSE]
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
  
cols_incl =function(var_threshs, incl = names(self$norm),g_incl = NULL, excl = list()){ #incl, g_incl,qq
  norm = self$norm
 names(incl) = incl
 useall = length(g_incl)==1 & g_incl[[1]]=="all"
  lapply(incl, function(norm_nme){
     varthreshs = var_threshs[[norm_nme]]
      norm1 = norm[[norm_nme]]
      varthresh = var_threshs[[norm_nme]]
      variance = self$vars[[norm_nme]]
     # varthresh = var_thresh[[norm_nme]]
      var_res = variance>=varthresh
      if(!is.null(g_incl) && length(g_incl)>0 && !useall){
        var_res = var_res & (names(norm1) %in% g_incl)
      }
      excl1 = excl[unlist(lapply(excl, function(e1)e1[[1]]==norm_nme))]
      var_res[(names(var_res) %in% unlist( lapply(excl1, function(e1) e1[[2]])))]=FALSE
      
      #var_res = var_res[!(names(var_res) %in% unlist( lapply(excl1, function(e1) e1[[2]])))]
      var_res
  })
},
  updateY=function(y1,weights, preprocessed=FALSE,family=NULL,CHECK=TRUE, all_v_all=FALSE, one_v_rest=FALSE){ ## updates y
    if(preprocessed){
      self$y = y1;
      self$family = family;
      self$weights = rep(1, nrow(y1[[1]]))
    }
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
          return(.expandAllvAll(y1[mi1,inds,drop=FALSE]))
        }else if(one_v_rest){
          return(.oneVRest(y1[mi1,inds,drop=FALSE]))
        }
        
        return(.expandFactors(y1[mi1,inds,drop=FALSE]))
      }else if(fam=="ordinal"){
        l1 = lapply(inds, function(ind){
        #  print(ind)
          tryCatch({
          Matrix(as.matrix(y1[mi1,ind,drop=FALSE]))
          },error=function(ew){
            return(NULL)
          })
        })
        l1 [unlist(lapply(l1, length))>0]
        
      }else{
        if(typeof(y1)=="S4"){
          return (list(y1[mi1,inds,drop=FALSE]))
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
        
        return(list(Matrix(as.matrix(y1[mi1,inds,drop=FALSE]), sparse=TRUE)))
      }
    })
    if(all_v_all) names(y) = gsub("multinomial","binomial.multiway",names(y))
    if(one_v_rest) names(y) = gsub("multinomial","binomial",names(y))
   self$family = unlist(lapply(names(y), function(nme) rep(gsub(".multiway","",nme),
                                                           if(is.list(y[[nme]])) length(y[[nme]]) else ncol(y[[nme]]))),recursive=FALSE)
   self$y = unlist(y, recursive=FALSE)
   self$weights = weights
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
        dir.create(d1, recursive=TRUE, showWarnings = FALSE)
        d1
      })
      self$norm_done = lapply(mem_dirs,function(xx){
        open_dataset(xx)|> collect()
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
      fi1= grep(fi,dir(dist$dir),value=TRUE)
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
    }), recursive=FALSE)
    matrices = unlist(lapply(matrices_all,function(m1){
      m1[names(m1)!="NA"]
    }), recursive=FALSE)
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
        #normm=   -1 *biganalytics::apply(d1,2, function(g) sqrt(sum((g-mean(g))^2)))
        normm=   -1 *apply(d1,2, function(g) sqrt(sum((g-mean(g))^2)))
        
      }
      #      print(normm)
      normm
      
    })
    self$norms_list = lapply(self$norm, function(n) list())
    
  },
updateTransforms=function(transform_x){
  self$transforms =.convertToTransform(transform_x)
},
split=function(proportions){
 
  db_name = self$db_name
  #    p = proportions[[1]]
  prop1 = cumsum(proportions)
  nrow = self$nrow
  prop2 = c(0,round(prop1*nrow))
  inds_new = lapply(1:(length(prop2)-1), function(i){
    start = prop2[i]+1
    end = prop2[i+1]
    inds = prop2[i]:prop2[i+1]
  }); 
  names(inds_new) =  paste(db_name,prop2[-1],sep=".")
  mats = lapply(inds_new, function(inds){
    nme_d = names(self$data); names(nme_d) = nme_d
      lapply(nme_d, function(nme){
          list(
          matrix = self$data[[nme]][inds,,drop=FALSE],
          matrixNA =  self$dataNA[[nme]][inds,,drop=FALSE]
          )
        })
  })
    ys =  lapply(inds_new, function(inds){
      lapply(self$y, function(y1){
      y1[inds,,drop=FALSE]
    })
    })
  list(mats = mats, ys = ys);  
},
  initialize=function( cohort,  db_name,dbDir,flags,
                      incl_full=TRUE,seed = 42, memDir = NULL) { ## mem_dirp is for saving scores
  
    self$direction = .readFlag(flags,'direction','none');
    self$db_name=db_name;
    if(!is.null(memDir)){
      mem_dirp = memDir
      mem_dirp = paste(memDir,nrow(cohort$rna$matrix),sep="/")  #getOption("fspls.mem_dir", NULL)
      dir.create(mem_dirp, recursive =TRUE, showWarnings=FALSE);
     nme_data = names(data)
     names(nme_data) = nme_data
     mem_dirs = lapply(nme_data,function(xx){
       d1 = paste(mem_dirp,xx,sep="/") 
       dir.create(d1, recursive=TRUE)
       d1
     })
     
     self$norm_done = lapply(mem_dirs,function(xx){
       open_dataset(xx)|> collect()
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
    #if(nrow(self$y[[1]])!=self$nrow)stop("problem")
    if(nrow(self$dataNA[[1]])!=self$nrow){
        lens = unlist(lapply(self$dataNA, function(d) length(which(d>0)))  )
        if(max(lens)>0) stop("!! NA matrix has different dimensions")
        ##otherwise we fix
        warning("NA matrix has diff dimensions")
        self$dataNA = lapply(self$data, function(d1) {
          dimnames = dimnames(d1); dims= unlist(lapply(dimnames, length))
        sparseMatrix(c(), c(),  dims = dims,dimnames = dimnames)
        })
    }
    self$train=NULL
   
     
    #self$ypreds_all= ypredObj$new(self,NULL, maxn,family)
    nrowd =  nrow(self$data[[1]])
    ym = array(1,dim = c(1, nrowd))
    self$mean_x = lapply(self$data, function(d1){
      mx = attr(d1,"mean_x")
      if(is.null(mx)){
        mx =if(typeof(d1)=="S4")  (ym[,1:nrow(d1),drop=FALSE] %*% d1)[1,]/sum(ym) else    dgemm( A=ym[,1:nrow(d1),drop=FALSE], B=d1)[1,]/sum(ym)
        #mx =if(typeof(d1)=="S4")  (ym %*% d1)[1,]/sum(ym) else    dgemm( A=ym, B=d1)[1,]/sum(ym)
      }
      mx
      })
    
  
    nmes1 = names(self$data)
    #print("HHHH");
    self$norm = lapply(self$data, function(d1) {
      
      normm = attr(d1,"norm")
      if(is.null(normm)){
        if(typeof(d1)=="S4"){ ## this not efficient
         # print(dim(d1))
        #  print(typeof(d1))
          means = Matrix::rowMeans(d1)
          
         # normm=   -1 *apply(d1,2, function(g) sqrt(sum((g-mean(g))^2)))
          normm = -1*sparse_norm(d1)
        }else{
          #normm=   -1 *biganalytics::apply(d1,2, function(g) sqrt(sum((g-mean(g))^2)))
          
        normm=   -1 *apply(d1,2, function(g) sqrt(sum((g-mean(g))^2)))
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
    
   
    
   
    #self$default_transform = eval(str2lang("function(x,pow) x"))
  }
)
)


