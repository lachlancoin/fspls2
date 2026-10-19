nullLikelihood <- function(y, w = NULL, family = c("binomial", "gaussian", "ordinal","multinomial")) {
  family <- match.arg(family)
  if (is.null(w)) w <- rep(1, length(y))
  ok <- !is.na(y) & !is.na(w)
  y1 <- y[ok]; w1 <- w[ok]
  n  <- sum(w1)
  if(family=="multinomial"){
    y1 <- droplevels(factor(y1, ordered = TRUE))
    levsy = levels(y1)
    tab <- tapply(w1, y1, sum)
    tab <- tab[tab > 0]                        # drop zero-weight classes
    m   <- sum(tab)
    ll  <- sum(tab * log(tab / m))
    attr(ll, "df") <- length(tab) - 1
  }else if (family == "ordinal") {
    y1 <- droplevels(factor(y1, ordered = TRUE))
    nj <- tapply(w1, y1, sum)
    nj <- nj[nj > 0]
    ll = (sum(nj * log(nj / sum(nj))))
    attr(ll,"df") = length(levels(y1))-1
  }else if (family == "binomial") {
    if (!all(y1 %in% c(0, 1))) stop("binomial needs a 0/1 response")
    mu <- sum(w1 * y1) / n
    if (mu <= 0 || mu >= 1){
      ll = 0;attr(ll,"df")=0
              # all one class
    }else{
      ll = (n * (mu * log(mu) + (1 - mu) * log1p(-mu)))
      attr(ll,"df") = 1;
    }
  }else{
  # gaussian
  mu  <- sum(w1 * y1) / n
  rss <- sum(w1 * (y1 - mu)^2)
  ll = -n / 2 * (log(2 * pi * rss / n) + 1)         # frequency weights
  attr(ll, "df")=2
  # for precision weights, add: + 0.5 * sum(log(w1))
  }
  ll
}

fitOrdinal<-function(y,x, lambs,family, firth=TRUE, steps = getOption("clm.steps",1)){
  colnames(x)=paste("x", 1:ncol(x), sep="")
  df = data.frame(cbind(y,as.matrix(x) ))
  df$y = factor(y, levels = sort(unique(y)))  
  func = paste0("y~",paste(colnames(x), collapse="+"))
  tryCatch({
  fits = clm_ao_path(func, df, lambdas = lambs, firth = FALSE, steps = steps)
  
  cfs= clm_ao_path(func, df,lambs, firth=FALSE, steps=steps)
  max_ind = which.max(unlist(lapply(cfs, function(x) x$logLik)))
  if(firth){
    m1r = clm_ao(func, df, lambda = cfs[[max_ind]]$lambda, start = cfs[[max_ind]]$par, steps=1, firth=TRUE)
  }else{
    m1r = cfs[[max_ind]]
  }
  const_term = if(family=="binomial") -1*m1r$alpha else m1r$alpha;
  list(beta_new1 = m1r$beta, const_term =const_term, lambda = m1r$lambda, logLik = m1r$logLik,
               extra_df = 1)
 
   
 
 
  }, error=function(errw) {
    return(NULL);
    tryCatch({
    resz = .getZetaBinary(x,y, w)  ##need as.matrix(x) ??
    return(list(beta_new1=resz$beta_new1,const_term =resz$const_term, lambda=1, logLik =resz$logLik ,
               extra_df =1))
    }, error=function(errw1) {
      
    
      return(NULL)
    })
  })
}

fitGLMnet<-function(y,x,family,w, fams){
  nonNAy = !is.na(y)
#  names(lambs) = lambs;
  n = length(which(nonNAy));
  ncolx = ncol(x)
  if(ncolx==1){
    ones = rep(1,nrow(x))
    x = cbind(ones, x)
  }
  tryCatch({
    rbetas = lapply(fams, function(fam){
   
    
    ridge <-  glmnet::cv.glmnet(x[nonNAy,,drop=FALSE],y[nonNAy],family=fam,
      #df1, y,family = family,
                                   alpha = 0,
                                   standardize = FALSE,
      weights=w[nonNAy], #weights = self$weights[nonNAy],
                                relax=TRUE,
                                   lambda.min.ratio = 1e-10,
                              
                                   nfolds = 10) 
      #  ridge=glmnet(x[nonNAy,,drop=FALSE],y[nonNAy],family=linkfun,weights=w[nonNAy], alpha = 0,
      #           standardize=FALSE,
      #           lambda = getOption("lambda"), lambda.min.ratio = 1e-10)
    
    rbeta <- coef(ridge,s=min(ridge$lambda),gamma = "gamma.min")
          mu <- predict(ridge,s="lambda.min", newx = x[nonNAy,,drop=FALSE], type = "response")     # n x nlambda matrix
         
           if(family=="binomial"){
            logLik <- colSums(y * log(mu) + (1 - y) * log1p(-mu))
         }else{
          # gaussian
            rss      <- colSums((y - mu)^2)
            logLik <- -n / 2 * (log(2 * pi * rss / n) + 1)
         }
      attr(rbeta, "logLik")=logLik;
      rbeta
    })
     lls = unlist(lapply(rbetas, function(xx)attr(xx,"logLik")))
     rbeta = rbetas[[which.max(lls)]];
    const_term = rbeta[1,1]
    beta_new1 =  if(ncolx==1) rbeta[-(1:2),1] else rbeta[-1,1]
    
   
  
    
    list(const_term = const_term, beta_new1 = beta_new1, lambda = as.numeric(names(lambs)[[which.max(lls)]]),logLik = attr(rbeta,"logLik"), extra_df = 1)
  }, error=function(errw) {
    warning("reverting to normal glm")
    fitGLM(y,x,family,w, fams)
  })
  
}

fitGLM<-function(y,x, family,w,
                 fams,
                 firth=FALSE, type="AS_mean"){ ## type"MPL_Jeffreys
  #if(family!="binomial") stop("needs to be binomial")
 # names(lambs) = lambs;
  if(is.null(colnames(x))){
    colnames(x) = paste("x",1:ncol(x),sep="")
  }
  # if(is.null(dim(x))) return( fitGLM(y,x,family,w, lambda, firthIfFail=FALSE))
  tryCatch({
    func = paste0("y~",paste(colnames(x), collapse="+"))
    df = as.data.frame(x); 
    df$y = if(family=="binomial") factor(y, levels = sort(unique(y), decreasing=T))  else y; 
    rbetas = .lapply_nme(fams, function(lambda){
    
    if(family=="binomial"){
      
        if(firth){
          m_mean <- glm(func, data = df, family = fams[[lambda]], method = "brglm_fit", type = type)
        }else{
          m_mean <- glm(func, data = df, family = fams[[lambda]])
          
        }
        
        coeffs = coef(m_mean)
        beta_new1  = -1*coeffs[-1]
     
        names(beta_new1)= colnames(x)
        const_term  = -1 *coeffs[1];
    }else{
     # fam = gaussian(link = sinh_asinh_link(lambda))
      if(firth && FALSE){ ## no point
        m_mean <- glm(func, data = df, family = fams[[lambda]], method = "brglm_fit", type = type)
      }else{
        m_mean <- glm(func, data = df, family = fams[[lambda]])
        
      }
      coeffs = coef(m_mean)
      beta_new1  = coeffs[-1]
      names(beta_new1)= colnames(x)
      const_term  = coeffs[1];
    }
   
      list(const_term = const_term, beta_new1 = beta_new1, lambda =as.numeric(lambda), logLik= logLik(m_mean), extra_df = 1)
    })
    rbeta = rbetas[[which.max(unlist(lapply(rbetas, function(xx) xx$logLik)))]]
  
    
    return(rbeta)
  } ,error=function(e) {
    message(conditionMessage(e)) 
    if(firth==FALSE && family!="gaussian" && getOption("revertToFirth", TRUE)){
      return(fitGLM(y,x, family,w,fams, firth=TRUE, type="AS_mean")) ## type"MPL_Jeffreys
        
    }else{
      return(NULL);
    }
  })
}


.averageResults<-function(res4){
  if(length(res4)==0) return(NULL)
  if(length(res4)==1)return(res4[[1]])
  joins = lapply(res4, function(a){
    a1=a|> unite("join", subpheno,pheno, family,measure,cv, beam,sep="_")
    a1$join
  })
  join1 = joins[[1]];
  if(length(which(duplicated(join1)))>0) stop("duplicated elements not supposed to occur")
  for(join in joins[-1]){
    join1 = join1[join1 %in% join]
  }
  
  res_final = res4[[1]][match(join1,joins[[1]]),,drop=F]
  if(nrow(res_final)==0) return(res_final)
  error_high = (res_final$high-res_final$mid)^2
  error_low = (res_final$mid-res_final$mid)^2
  
  for(j  in 2:length(res4)){
    indsi2 = match(join1,joins[[j]])
    error_high = error_high +( res4[[j]]$high[indsi2] -  res4[[j]]$mid[indsi2])^2
    error_low = error_low +( res4[[j]]$low[indsi2] -  res4[[j]]$mid[indsi2])^2
    
    res_final$mid =  res_final$mid + res4[[j]]$mid[indsi2]
    
    
    
  }
  error_low = sqrt(error_low) / length(res4);
  error_high = sqrt(error_high) / length(res4);
  
  res_final$mid = res_final$mid/length(res4)
  res_final$low = res_final$mid - error_low;
  res_final$high = res_final$mid +error_high
  
  #res_final$cv_index = rep("avg", nrow(res_final))
  res_final$cv_index = rep("avg", nrow(res_final))
  res_final
  
}

.better1<-function(rmsv, rmsv2){
  measure = rmsv$measure[[1]]
  subpheno = rmsv$subpheno[[1]]
  cv = rmsv$cv[[1]]
  rmsv_ = subset(rmsv, cv==cv & measure==measure & subpheno == subpheno) |> pivot_wider(names_from="pheno", values_from="mid")
  ncol1 = ncol(rmsv)-2
  rmsv2_ = subset(rmsv2, cv==cv & measure==measure & subpheno == subpheno) |> pivot_wider(names_from="pheno", values_from="mid")
a1 = rbind(rmsv_[1,-(1:ncol1)], rmsv2_[1,-(1:ncol1)])

a2 = cbind(c("bef","aft"),a1)
names(a2)[1] = "nme"
a2
}

default_types=jsonlite::fromJSON('{"gaussian": "correlation","binomial" : "AUC","multinomial" : "AUC","ordinal" :"AUC"}')
.convertToTransform<-function(transform_x, fromJSON=is.character(transform_x)){
  str1 = if(fromJSON) jsonlite::fromJSON(transform_x) else transform_x 
lapply(str1, function(t_y){
  t_y1 =  lapply(t_y[1:2], function(t_y1) if(is.character(t_y1)) eval(str2lang(t_y1)) else t_y1)
  t_y1$params = t_y[[3]]
  .checkInverse(t_y1)
  t_y1 
})
}





## assumes Wall1 is upper diagonal
## uses data with mean subtracted
## does not add any constant term

.eval1_<-function(x_,  Wall2, transf,  family, CHECK=FALSE){
  if(ncol(Wall2) > 1){
    if( max(abs(Wall2[lower.tri(Wall2)]))>0) stop("!!")
  }
  #if(is.null(dim(beta_new2))) beta_new2 =as.matrix(beta_new2, nrow = ncol(x_), ncol = 1)
 # print(dim(x_)); print(dim(Wall2));
  if(ncol(x_)==0) return(x_);
  if(getOption("projectAfterTransform",FALSE)){ ## PRETTY SURE THIS IS WRONG
    
    x2 = data.frame(vapply(1:ncol(x_), function(jj){
      t1 = transf[[jj]]$func(x_[,jj], transf[[jj]]$param)
     # t1 = transf[[jj]]$func(x_[,1:jj] %*% Wall2[1:jj,jj,drop=FALSE], transf[[jj]]$param) 
      as.vector(t1)
    }, rep(0, nrow(x_))))
    t2 = data.frame(vapply(1:ncol(x_), function(jj){

       t1 = as.matrix(x2[,1:jj,drop=FALSE]) %*% Wall2[1:jj,jj,drop=FALSE]
      as.vector(t1)
    }, rep(0, nrow(x_))))
   
  }else{ ##projection first and then transformation 
    t2 = data.frame(vapply(1:ncol(x_), function(jj){
      t1 = transf[[jj]]$func(x_[,1:jj] %*% Wall2[1:jj,jj,drop=FALSE], transf[[jj]]$param) 
     as.vector(t1)
    }, rep(0, nrow(x_))))
   
  }
  if(nrow(x_)==1) t2 = t(t2)
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
  list(beta_new1 = beta_new1, const_term = const_term, pv = pv, logLik=ll1)
}

.mkBinary<-function(y, thresh=mean(y,na.rm=TRUE)){
  yn = rep(NA, length(y))
  yn[y<=thresh]=0
  yn[y>thresh]=1
  yn
}


multinom_nnet <- function(x, y, w = NULL, decay = 1e-3, standardise = FALSE,
                          maxit = 500) {
  print("NNET")
  x <- as.matrix(x)
  orig_names <- colnames(x)
  if (is.null(orig_names)) orig_names <- paste0("x", seq_len(ncol(x)))
  if (is.null(w)) w <- rep(1, length(y))
  
  #ok <- !is.na(y) & !is.na(w) & complete.cases(x)
  #x <- x[ok, , drop = FALSE]; y <- y[ok]; w <- w[ok]
  
  tbl <- sort(table(y), decreasing = TRUE)
  tbl <- tbl[tbl > 0]
  lev <- names(tbl)                       # most frequent class first = reference
  K   <- length(lev)
  
  # Standardise so the ridge penalty treats predictors equally
  if (standardise) {
    ctr <- colMeans(x)
    scl <- apply(x, 2, sd); scl[scl == 0] <- 1
    xs  <- sweep(sweep(x, 2, ctr), 2, scl, "/")
  } else {
    ctr <- rep(0, ncol(x)); scl <- rep(1, ncol(x)); xs <- x
  }
  
  xn <- paste0("x", seq_len(ncol(x)))
  df <- as.data.frame(xs); names(df) <- xn
  df$y <- factor(as.character(y), levels = lev)
  df$w <- w
  fml <- reformulate(xn, response = "y")
  
  m <- nnet::multinom(fml, data = df, weights = w, decay = decay, maxit = maxit,
                      MaxNWts = 1e7, trace = FALSE)
  if (m$convergence != 0) warning("multinom did not converge within maxit")
  
  cf <- coef(m)
  if (is.null(dim(cf))) {                 # K = 2: coef() returns a named vector
    cf <- matrix(cf, nrow = 1, dimnames = list(lev[2], names(cf)))
  }
  cf <- cf[, c("(Intercept)", xn), drop = FALSE]
  
  # Reference-class form -> sum-to-zero (symmetric) form, as glmnet reports
  full <- rbind(0, cf)
  rownames(full) <- c(lev[1], rownames(cf))
  full <- full[lev, , drop = FALSE]
  sym  <- sweep(full, 2, colMeans(full))  # K x (p+1)
  
  # Back-transform to the original predictor scale
  beta_orig  <- sweep(sym[, -1, drop = FALSE], 2, scl, "/")
  const_orig <- sym[, 1] - drop(beta_orig %*% ctr)
  coefmat <- cbind(const_orig, beta_orig)
  #colnames(coefmat) <- paste0("X", 0:ncol(x))
  rownames(coefmat) <- lev
  
  # Unpenalised multinomial log-likelihood at the fitted coefficients
  eta <- cbind(1, x) %*% t(coefmat)
  P   <- exp(eta - apply(eta, 1, max)); P <- P / rowSums(P)
  idx <- cbind(seq_len(nrow(P)), match(as.character(df$y), lev))
  ll  <- sum(w * log(pmax(P[idx], 1e-300)))
  mi1 = match(rownames(coefmat), levels(y))
  list(coef = coefmat, constant = coefmat[mi1, 1],
       beta = t(coefmat[mi1, -1, drop = FALSE]),
       extra_df = K - 1, logLik = ll, levels = lev,
       decay = decay, convergence = m$convergence)
}
#pv1 = .calcPvalue(x_trans[nonNA ,,drop=FALSE]%*% betas[[kk]],yp1[,kk,drop=FALSE],y, w, family)
multinom_firth <- function(x, y, w = NULL, type = "AS_mean") {
  print("FIRTH")
  
  x <- as.matrix(x)
  orig_names <- colnames(x)
  if (is.null(orig_names)) orig_names <- paste0("x", seq_len(ncol(x)))
  if (is.null(w)) w <- rep(1, length(y))
  
  ok <- !is.na(y) & !is.na(w) & complete.cases(x)
  x <- x[ok, , drop = FALSE]; y <- y[ok]; w <- w[ok]
  lev = levels(y)
  tbl <-table(y)
                           # most frequent class first = reference
  K   <- length(levels(y))
  
  xn <- paste0("x", seq_len(ncol(x)))            # safe internal names
  df <- as.data.frame(x); names(df) <- xn
  df$y <- y
  df$w <- w
  fml  <- reformulate(xn, response = "y")
  m    <- brglm2::brmultinom(fml, data = df, weights = w, type = type, ref = which.max(tbl))
  
  cf <- coef(m)                                  # assumed (K-1) x (p+1); check with str(cf)
  cf <- cf[, c("(Intercept)", xn), drop = FALSE] # fix the column order by name
  full <- rbind(0, cf)
  rownames(full) <- c(lev[which.max(tbl)], rownames(cf))
  full <- full[lev, , drop = FALSE]
  sym  <- t(sweep(full, 2, colMeans(full)))      # (p+1) x K, sum-to-zero across classes
  rownames(sym) <- c("(Intercept)", orig_names)
  
  eta <- cbind(1, x) %*% sym                     # n x K, columns in `lev` order
  P   <- exp(eta - apply(eta, 1, max)); P <- P / rowSums(P)
  ll  <- sum(w * log(P[cbind(seq_len(nrow(P)), match(as.character(df$y), lev))]))
  
  list(constant = sym[1, ], beta = sym[-1, , drop = FALSE],
       extra_df = K - 1, logLik = ll, levels = lev)
}


multinom_ridge<-function(x,y,w,lambda=getOption("lambda")){
  degf = (ncol(x)-1) * (length(levels(y))-1)
  
  
    # y must be a factor; the first level is the reference class by default
   # d$y <- factor(d$y)
    
    
  
  if(is.null(dim(x)) || ncol(x)==1){
    one_inds = 1:2
    x = cbind(1,x)
  }else{
    one_inds = 1 
  }
#  m1=(multinom(y~x,weights=w, trace=FALSE))
#  sm1  = summary(m1)
  
  
  ridge<- tryCatch({
    glmnet::cv.glmnet(x,y, family="multinomial", alpha=0, weights=w,lambda.min.ratio=1e-10, standardize=FALSE) ## should try relax?
    }, #, lambda.min.ratio = 1e-10),
    error = function(e) {
      if (grepl("zero variance", conditionMessage(e))) {
        NULL
      } else {
        stop(e)
      }
    }
  )
 if(is.null(ridge)) return(NULL)
     
   
    
    
  
#predy=predict(ridge,x,s=min(ridge$lambda))[,,1]
  rbeta <- coef(ridge, s = "lambda.min",gamma = "gamma.min")
  rdf = as.matrix(data.frame(lapply(rbeta, as.matrix)))
  dimnames(rdf)[[2]]= names(rbeta)
  mu <- predict(ridge, s = "lambda.min", newx = x, type = "response")[, , 1]   # n x K matrix
  
  #mu =predict(ridge,s="lambda.min", gamma = "gamma.min",newx = x, type = "response")[,,1]     # n x nlambda matrix
  yi <- as.integer(factor(y, levels = colnames(mu)))     # column index of each observed class
  logLik <- sum(log(mu[cbind(seq_len(nrow(mu)), yi)]))
  attr(logLik ,"df") = degf
  list(const = apply(rdf[one_inds,,drop=FALSE],2,sum), beta = rdf[-one_inds,],
   extra_df=length(levels(y))-1,
     logLik = logLik )
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

  attr(m1,"factor")= data.frame( lapply(data.frame(y), function(y1){
    if(!is.factor(y1)) y1 = factor(y1, levels = sort(unique(y1)))
    y1
  }));
  m1;
}

.expandFactors<-function(y, max_cats=50){
  res = lapply(1:ncol(y), function(i){
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
  folds="list",
  db_name="character",
  dataNA="list",
  types="vector",
  nrow="integer",
  vars="list",
#  direction="character",
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
  fams_all = "list",
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
calcBetaProj1=function(subphens,k,b_i,b_i_name, prev_var, Wall,prev_i2,
                       convert=TRUE, strict=FALSE,project=TRUE,useoffset=TRUE){
  if(convert){
    b_i = self$convert(b_i)
    if(length(b_i)<2) {
      pvs = unlist(lapply(phensi,function(xx) 999)) ## returning large positive value since its not in dataset
       return(list(pvs=pvs, pvs_sep = pvs))
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
    family=self$family[[phensi1]]
    
    fams = self$fams_all[[phensi1]];
    phensi_=subphens[[nme]]
    Wall1 = Wall[[nme]]
    useglm=getOption("useglmnet",FALSE)
    
    if(family=="multinomial" && !useglm){
      family="binomial";
    }
    if(length(prev_var)==0) {
      project=FALSE ;#family!="multinomial"
      useoffset=FALSE;
    }
    
#    useoffset = length(prev_var)>0 && family!="multinomial"
    res1=self$calcBetaProjAll(nme,phensi_,family,  k, b_i, b_i_name,prev_var, Wall1,
                              prev_i2,fams,useglm=useglm,
                              project=project, strict=strict,useoffset=useoffset)
   
   # data.frame(res1$betas)
    if(family=="multinomial"){
     # betas_new = res1$betas
          
        res1$constants = res1$constants[[1]] 
        res1$tbls = res1$tbls[[1]] 
        res1$betas_proj = res1$betas_proj[[1]]
        if(!is.matrix(res1$betas_proj)){
          res1$betas_proj = matrix(res1$betas_proj, nrow =1 , ncol = length(res1$betas_proj))
          colnames(res1$betas_proj) = names(res1$constants)
#        }
        }
        res1$betas = res1$Wall %*%res1$betas_proj  ## assuming no projection
    }else{
      res1$betas_proj = as.matrix( data.frame(res1$betas_proj))
      res1$betas = res1$Wall %*% res1$betas_proj
    }
 
    res1
  })
  nme = names(res[[1]])
  
  names(nme)=nme
  res2 = lapply(nme, function(n){
    lapply(res, function(r) r[[n]])
  })
  res2
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
calcBetaProjAll=function(nme,phensi_,family, k,b_i,b_i_name, prev_var, Wall1,prev_i2,fams,
                         project=TRUE, strict=FALSE, 
                         useglm=getOption("useglmnet",FALSE),
                         
                         CHECK=getOption("fspls.check",TRUE),
                        useoffset =TRUE){
  betas1 = prev_i2$betas_proj[[nme]];
  constants1 = prev_i2$constants_proj[[nme]];
  prevLL =prev_i2$logLiks[[nme]];
  prev_df = prev_i2$degsf[[nme]]
  data = self 
  firth = getOption("firth",FALSE);
  if(length(prev_var)==0){
    useoffset = FALSE;
    project=FALSE;
  }
  
  nonNA = self$looc$incl[,k]
  vars1 = c(prev_var,list(b_i)) 
  names(vars1)[length(vars1)]  = paste(b_i_name,collapse=".")
  data$updateUDVP(prev_var)
 
  ys =self$y[[nme]]##  if(family %in% c("binomial","ordinal")) (d$y1) else if(family =="multinomial") d$y2 else d$y
  if(family=="multinomial" && useglm){
    useoffset=FALSE;
    ys = data.frame(list(attr(ys,"factor")))
    names(ys) = nme
  }else if(length(phensi_) < ncol(ys)){
    ys = ys[,phensi_,drop=FALSE] 
  }
  ncoly = ncol(ys)
  betas_proj = apply(ys,2,function(v1) list())# apply(train$y,2, function(v1) list())
  betas = apply(ys,2,function(v1) list())# apply(train$y,2, function(v1) list())
  
  lambdas = rep(NA, length(betas))
  logLiks_new = apply(ys,2, function(v1) list())
  degsf = apply(ys,2, function(v1) list())
  
  tbls =  apply(ys,2, function(v1) list())
  constants = apply(ys,2, function(v1) list())
  pvs =  apply(ys,2, function(v1) list())
  pvs_sep =  apply(ys,2, function(v1) list())
  j = length(vars1)
 
  zero_model=length(vars1[[j]])==0
  
  if(zero_model){
    # stop(" this doesnt happen")
    non_na_x = rep(TRUE, self$nrow)
    useglm=FALSE
  }else{
   # non_na_x = if(is.null(data$dataNA[[vars1[[j]][1]]])) rep(TRUE,nrow(ys) ) else !(data$dataNA[[vars1[[j]][1]]][, vars1[[j]][2]] )
    non_na_x = if(is.null(data$dataNA[[vars1[[j]][1]]])) rep(TRUE,nrow(ys) ) else !(data$dataNA[[vars1[[j]][1]]][, vars1[[j]][2]] )
  }
  if(length(which(non_na_x))==0){
    #  print(vars1);
    stop("no non_zero");
  }
  if(zero_model){
    #stop(" this doesnt happen")
    zero_model = TRUE
    x_=matrix(1, nrow=self$nrow, ncol=1)
    colnames(x_)="ones"
    meanx = 1
  }else{
    
  #  x_ =self$extractData(vars1, adjust=FALSE, convert=FALSE)
  #  meanx = apply(x_,2,mean)
    
    x_ =self$extractData(vars1, adjust=FALSE, convert=FALSE)
    meanx = t(replicate(nrow(x_), apply(x_,2,mean,na.rm=TRUE)))
    colnames(x_) = names(vars1)
  }
    
  
  transf = self$getTransforms(vars1) #lapply(vars1, function(v1) self$transforms[[v1[[3]]]][[2]])
  UDV = data$UDVP ## should check it corresponds to prev_i
  #Wall2 = UDV$getWall(x_[,ncol(x_)]-meanx[,ncol(x_)], Wall1) ## updated projection  may not need to subtract meanx  .. leaving it for now
  Wall2 = UDV$getWall(x_[,ncol(x_)], Wall1) ## updated projection  may not need to subtract meanx 
  
    if(!project) Wall2 = diag(ncol(Wall2))
  if(CHECK && !is.null(UDV$VDU) && !project){
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
   x_trans_prev  =  x_trans[,-ncol(x_trans),drop=F];
   #x_trans_prev = .eval1_(x_[,-ncol(x_),drop=FALSE],  Wall1, transf[-ncol(x_)], family)
   
#x_trans=  t(t(x_trans) -apply(x_trans,2,mean)  )
#   sum(abs(transf1$func(x1_, transf1$param[1])-x_trans[,ncol(x_trans)] - ( mean(x1_[,1]) - mean(x_trans[,ncol(x_trans)]))))
  
  if(!useoffset ){
    yp1 = x_trans_prev[nonNA,,drop=FALSE]
  }else{
    yp1 =   x_trans_prev[nonNA,,drop=FALSE] %*% betas1
  }
 #  lambs=if(family=="gaussian") getOption("gaussian.link", 0.0) else  getOption("binomial.link", 1.0);
#   names(lambs)=lambs
   
  for(kk in 1:ncoly){
  #  print(paste(kk,ncoly))
    nonNAk = nonNA & non_na_x
    nonNAk1 = non_na_x[nonNA]; ## second subset
    y=ys[,kk][nonNAk]
    nonNAy = !is.na(y)
    w = data$weights[nonNAk]
    beta_new1=0;
    const_term=0
    
    if(useoffset ){
      ## in this case yp1 is in dimensions of y since multiplied by beta
      x = cbind(yp1[,kk], x_trans[nonNA,ncol(x_trans)])
      
        dimnames(x)[[2]] = c(paste0("A",1:(ncol(x)-1)),"x")
      
    }else{
      x = cbind(yp1, x_trans[nonNA,ncol(x_trans)])
      
    }
    varx = apply(x,2,var,na.rm=T)
#    print(family); print(head(y));
    if(family=="gaussian" && var(y,na.rm=TRUE)==0){
      sm2 = NULL
    }else if(length(table(y))<=1){
      sm2 = NULL
    }else if(family=="multinomial" && !useoffset){
       ty=as.list(table(y))
      tbls[[kk]] = ty[ty>0]
      #mean_x = apply(x,2, mean,na.rm=TRUE)
      #x = t(t(x) - mean_x)
      #print(apply(x,2,mean)) 
      
      rdf= multinom_ridge(x,y,w)
      if(is.null(rdf)){
        sm2 = NULL
      }else{
      beta1 = rdf$beta
      if(length(which(is.na(match(colnames(beta1), names(rdf$const)))))>0) beta1 = t(beta1)
     # print("HHHH")
    #  print(beta1);
    #  print(Wall2 %*% beta1)
        sm2 = list(const_term = rdf$const, beta_new1 =  beta1, lambda=1, logLik = rdf$logLik, 
                   extra_df = length(tbls[[kk]])-1)
      }
    }else if(ncol(x)==1 && var(x)==0 ){
      sm2 = NULL
    
    }else if(length(which(nonNAy))==0){
      warning("no nonNA y")
      sm2 = NULL;
    }else if(family=="ordinal" ){
      ty=as.list(table(y))
      tbls[[kk]] = ty[ty>0]
      lambs= getOption("binomial.link", 1.0);
        sm2 = fitOrdinal(y,x,  lambs, family, firth = firth);
    }else if(useglm){
      
      sm2 =  fitGLMnet(y[nonNAy],x[nonNAy,,drop=F],family,w[nonNAy], fams)
    
    }else{
      sm2 = fitGLM(y,x,family,w, fams, firth=firth)
     
    }
    if(zero_model ){
      lambdas[[kk]]=1.0;
      logLiks_new[[kk]] = nullLikelihood(y, w, family);
      degsf[[kk]] =attr( logLiks_new[[kk]],"df")
      if(is.null(degsf[[kk]])) stop("cannot be NULL");
      beta_new1 = rep(0, ncol(x))
      const_term = 0;  
       if(family=="multinomial"){
         const_term = rep(0, length(unique(y[nonNAy])))
         names(const_term) = unique(y[nonNAy]);
         beta_new1 = matrix(0, nrow=1, ncol = length(const_term));
         colnames(beta_new1) = names(const_term);
       }else if(family=="ordinal"){
         const_term = rep(0, length(unique(y[nonNAy]))-1)
         
       }
    }else if(is.null(sm2)){
      
        lambdas[[kk]]=1.0;
        logLiks_new[[kk]] = prevLL[[kk]];
        degsf[[kk]] = prev_df[[kk]];
        if(family=="multinomial"){
         beta_new1 = cbind(betas1, rep(0, nrow(betas1)))  ## add zero row
        }else{
          beta_new1 = rep(0, ncol(x))
        }
      
          const_term = constants1[[kk]];
        
        if(useoffset && family!="multinomial"){
          beta_new1[1] = 1;
        }
      }else{
             lambdas[[kk]] = sm2$lambda
             logLiks_new[[kk]] =as.vector(sm2$logLik)
             degsf[[kk]] = prev_df[[kk]]+sm2$extra_df
             beta_new1= sm2$beta_new1;
             const_term = sm2$const_term
      }
      if(useoffset ){ #project){

       ##betas[[kk]] = if(family=="multinomial") Wall1 %*% beta_new1 else   (Wall1 %*% beta_new1)[,1]
          betas_proj[[kk]] = rbind(betas1[,kk,drop=FALSE] * beta_new1[1], beta_new1[2])
        
      }else{
        betas_proj[[kk]] = beta_new1
      }
    
    
         constants[[kk]] = const_term  #-mean_adj*beta_new1
  
  }
   if(length(degsf)==0) stop("cannot be zero!");
  list( constants = constants, 
       betas_proj = betas_proj,
       tbls=tbls, Wall = Wall2, lambdas = lambdas, degsf = degsf, logLiks = logLiks_new)
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


fitEnsemble=function(preds, phens, CHECK=FALSE, verbose=FALSE,  steps = getOption("clm.steps",1)){
  nvars = 1:length(preds);names(nvars) = names(preds)
  nmes_phens = names(phens); names(nmes_phens) =nmes_phens;
firth = getOption("firth",FALSE);
coeffs =  
  lapply(nvars, function(i){
    nrep_inds = 1:length(preds[[i]]);
    names(nrep_inds)=names(preds[[i]]);
    a2=lapply(nrep_inds, function(k){
      k_name = names(nrep_inds)[k];
      if(length(preds[[i]][[k_name]])<=1) return(NULL)
    
    lambdas= lapply(preds[[i]][[k_name]], function(x) attr(x,"lambdas"))
  
    a0=  lapply(nmes_phens, function(nme_p){
        #print(nme_p)
        family = strsplit(nme_p,"\\.")[[1]][1]
        y = self$y[[nme_p]]
       if(verbose) print(paste('fit ens',i,k_name, nme_p))
        if(family=="multinomial"){
          preds1 =  (lapply(preds[[i]][[k_name]], function(predi){
            predi[[nme_p]]
          }))
          y = attr(y,"factor")
          include = !(unlist(lapply(preds1,is.null)))
          include_ind = which(include)
          nonNAy = self$looc$incl[,k] & !is.na(y)
          weights = self$weights[nonNAy]
          if(length(include_ind)==0){
            return(NULL)
           # res0 = matrix0,nrow = 1, ncol)
            #res0 = 1; names(res0)="(Intercept)"
            #return(res0)
          }
          df1 = as.matrix(data.frame(unlist(lapply(preds1[include_ind], function(x) data.frame(x[nonNAy,, drop=FALSE])), recursive =FALSE)))
         ###SHOULD CHECK THE NAMES
         
          
          rbeta = tryCatch({
            ridge_cv <-  glmnet::cv.glmnet(df1, y,
                                        family = family,
                                        alpha = 0,
                                        standardize = FALSE,
                                        weights = self$weights[nonNAy],
                                        lambda.min.ratio = 1e-10,
                                        lower.limits = 0 ,
                                        relax=TRUE,
                                        nfolds = 10)               # or however many folds you want
          
          # extract coefficients at the CV-optimal lambda
          rbets <- coef(ridge_cv, s = "lambda.min",gamma = "gamma.min")  
          rbeta = data.frame(lapply(rbets,function(rb) rb[,1,drop=TRUE]))
          return(rbeta)
          },error=function(ew){
            rbeta1 = matrix(0, nrow = 1+ncol(df1), ncol =length(levels(y)))
           
            dimnames(rbeta1)=list(c("(Intercept)",colnames(df1)), levels(y))
           rbeta1[1+(1:ncol(rbeta1)),1:ncol(rbeta1)]=diag(1, ncol(rbeta1))
            return (rbeta1)
          })
          return(rbeta)
          
        }
       
        if(family=="ordinal"){
          preds1 =  (lapply(preds[[i]][[k_name]], function(predi){
            predi[[nme_p]]
          }))

          include = !(unlist(lapply(preds1,is.null)))
          include_ind = which(include)
          nonNAy = self$looc$incl[,k] & !is.na(y[,1])
          weights = self$weights[nonNAy]
          if(length(include_ind)==0){
            return(NULL)
          }
          lambdas1 =unlist(lapply(lambdas[include_ind], function(x) x[[nme_p]]))
          lambda1 = mean(lambdas1)
          rbeta = tryCatch({
          df0 = as.matrix(data.frame(unlist(lapply(preds1[include_ind], function(x) data.frame(as.matrix(x[nonNAy,, drop=FALSE]))), recursive =FALSE)))
          ###SHOULD CHECK THE NAMES
          nonNAx = !is.na(df0[,1])
          
          df1 = data.frame(as.matrix(cbind(y[nonNAy,1,drop=F][nonNAx,,drop=F],df0[nonNAx,,drop=F])));
          colnames(df1)[1]="y"
          
          df1$y=factor(df1$y,ordered=TRUE) # levels = sort(unique(df1$y)))
          func = paste0("y~",paste(colnames(df1)[-1], collapse="+"))
          #m1 <- ordinal::clm(func, data = df1,  link = "logit")
                      
          m1 = clm_ao(func, df1,lower=0,lambda = lambda1, firth = firth, steps  =steps)
          
          
          rbeta1 =  m1$beta
          attr(rbeta1,"zeta") = m1$alpha
          attr(rbeta1,"lambda")=m1$lambda
          rbeta1
          
          },error=function(ew){
            print(ew)
            return(NULL);
#            rbeta1 = matrix(0, nrow = 1+ncol(df1), ncol =length(levels(y)))
#            
#            dimnames(rbeta1)=list(c("(Intercept)",colnames(df1)), levels(y))
#            rbeta1[1+(1:ncol(rbeta1)),1:ncol(rbeta1)]=diag(1, ncol(rbeta1))
            return (rbeta1)
          })
          return(rbeta)
          
        }
        ncoly = phens[[nme_p]]; names(ncoly) = ncoly
        ncoly_ind = 1:length(ncoly); names(ncoly_ind) = names(ncoly);
        li2 = lapply(ncoly_ind, function(y_ind){
          yname = ncoly[[y_ind]]
         # print(paste(i,k, nme_p, yname))
          ik = which(colnames(y)==yname)
          preds1 =  (lapply(preds[[i]][[k_name]], function(predi){
                  predi[[nme_p]][,ik]
          }))
          include = unlist(lapply(preds1,length))>0
          include_ind = which(include)
          nonNAy = self$looc$incl[,k] & !is.na(y[,ik])
          weights = self$weights[nonNAy]
          if(length(include_ind)==0){
            res0 = 1; names(res0)="(Intercept)"
            return(res0)
          }
          if(i==1 ){
            res0 = rep(0, 1+length(include_ind)); res0[2]=1; ## default all weight on first
            names(res0)[1] = "(Intercept)";names(res0)[-1] =names(preds1)[include_ind]
            return(res0)
          }
           df1 = data.frame(preds1[include_ind])[nonNAy,,drop=F];
           names(df1) = names(preds1)[include_ind]
           lambda1 =mean(unlist(lapply(lambdas[include_ind], function(x) x[[nme_p]][[y_ind]])))
           if(family=="binomial") {
             
             linkfun=  binomial(link = ao_link(lambda1, reflected = TRUE))
           }else{
             linkfun=  gaussian(link=sinh_asinh_link(lambda))
             
           }
           
          if(length(include_ind)==1){
            
            summ = summary(glm(y[nonNAy,ik]~df1[[1]], family=linkfun, weights = weights))
            rbeta = summ$coefficients[,1]
            if(length(rbeta)==1){ ##coeff only
              rbeta = c(0,1)
            }
            names(rbeta) = c("(Intercept)",names(df1)[[1]])
            return(rbeta);
          
          }else{
           
            df2 = as.matrix(df1);
            nonNAx = apply(df2,1,function(v) length(which(is.na(v)))==0)

         
           rbeta<- tryCatch({
              ridge_cv <- glmnet::cv.glmnet(df2[nonNAx,,drop=F], y[nonNAy, ik][nonNAx],
                                    family =linkfun ,
                                    alpha = 0,
                                    standardize = FALSE,
                                    weights = self$weights[nonNAy][nonNAx],
                                    lambda.min.ratio = 1e-10,
                                    relax=TRUE,
                                    lower.limits = 0 ,
                                    nfolds = 10)               # or however many folds you want
              
              # extract coefficients at the CV-optimal lambda
              rbets <- coef(ridge_cv, s = "lambda.min",gamma = "gamma.min")  
              rbets[,1]
             
            
            },error=function(errw){
              message(errw) 
         
              rbeta = rep(0, ncol(df2)+1); rbeta[2] = 1
              names(rbeta) = c("(Intercept)",names(df1))
              rbeta
            })
           attr(rbeta,"lambda")=lambda1
          
                 
          }
          if(CHECK && family %in% c("binomial","gaussian") && FALSE){
           
            
            if(family=="binomial"){
              aucs1 = apply(df1, 2, function(aa) roc(y[nonNAy,ik], aa, quiet=TRUE)$auc)
              aa = predict(ridge_cv,as.matrix(df1),s="lambda.min",family=family)
              aucs_new = apply(aa,2, function(aa2) roc(y[nonNAy, ik], aa2, quiet=TRUE)$auc)
            }else if(family=="gaussian"){
              aucs1 = apply(df1, 2, function(aa) cor(y[nonNAy,ik], aa))
              aa = predict(ridge_cv,as.matrix(df1),family=family)
              aucs_new = apply(aa,2, function(aa2) cor(y[nonNAy,ik], aa2))
              
            }    
            
            
            #print(paste(max(aucs1), max(aucs_new)))
            if(max(aucs1)> max(aucs_new)) {
              
              warning(paste("refitting did not improve auc",i,k,nme_p, yname,max_ind,nme_p, yname, max(aucs1), max(aucs_new)))
            }
          }
          rbeta
          })
        li2 = li2[unlist(lapply(li2, length))>0]
        data.frame(li2)
        })
    a0[unlist(lapply(a0, length))>0]
    #  })
  })
  a2[unlist(lapply(a2, length))>0]
  
  })
  coeffs[unlist(lapply(coeffs, length))>0]
  
},


#vars2 = vars_all[[1]]$variables[[1]];phens1 = phens[[1]]; k=1; useoffset=TRUE
makeModels=function(models, phens1, vars2,k, beam,
                    flags = list()
                  ){
  checkRMSV = .readFlag(flags,"checkRMSV",FALSE)
  logpthresh= log(.readFlag(flags,"pthresh",1e-3))
  project=.readFlag(flags,"project",TRUE)
  useoffset=.readFlag(flags,"useoffset",TRUE)
  verbose=.readFlag(flags,"verbose",FALSE)
  useglm1=.readFlag(flags,"useglmnet",FALSE)
    nonNA = self$looc$incl[,k]
   # if(verbose) print("making models")
  if(is.null(self$looc)){
    self$updateLOOC( phens1,flags,varn=c(),force=FALSE, verbose=FALSE)
  }
  phen2 = phens1
  ypred=self$ypred(phen2)
  phensi = self$phensi(phen2)
  subphens = phensi
  #nmes = c()
  data = self;
  len = length(vars2)
  #models = vector("list", len)
  fams1 = lapply(names(phens1), function(st1)strsplit(st1,"\\.")[[1]][1])
  #family = getOption("fspls.family",fams1) #ypred$family[[1]]
  #names(family)=family
  
 # if(family[[1]]=="multinomial") {
#    useoffset=FALSE
#    project=FALSE
#  }
  nme1 = ""
  b_i_name=c()
  #transform_x1=transform_x[[var_transf[[1]]]]
 # nme_c1 = var_transf[[1]]
  Wall =lapply(subphens, function(f) matrix(nrow=0, ncol=0))
  prev_i2 = self$makeNextModel(models, NULL,b_i_name,subphens,k, beam,
                              ypred=ypred, project=project,useoffset=useoffset,
                      Wall=Wall,
                              logpthresh =logpthresh)
  jk=1
 rmsv=if(checkRMSV) self$checkRMSV(subphens,prev_i2, ypred, nonNA,verbose=verbose) else NULL
  while(jk<=len){
    b_i_name = vars2[[jk]]
    if(verbose) cat(paste(jk,"of",len,";"))
    prev_i1 = try(self$makeNextModel(models, prev_i2,b_i_name,subphens,k,beam, Wall=Wall,
                                  ypred=ypred, project=project,  logpthresh =logpthresh))
    if(is.null(prev_i1)) break;
    
     if(inherits(prev_i1,"try-error")){ 
       warning("fitting error");
       break;
     }
    nme2 = paste(vars2[[
      jk]], collapse=".")
    nme1 = if(jk==1)  nme2 else paste(nme1, nme2,sep=";")
    if(checkRMSV){
      rmsv2=self$checkRMSV(subphens,prev_i1, ypred, nonNA,verbose=verbose)
     
      
      if(!is.null(rmsv)){
        ab=pivot_wider(rmsv, names_from="measure", values_from="mid")
        ab1 =as.matrix(ab[,ncol(ab)]); rownames(ab1) = ab$pheno;
        ab2=pivot_wider(rmsv2, names_from="measure", values_from="mid")
        ab21 =as.matrix(ab2[,ncol(ab2)]); rownames(ab21) = ab2$pheno;
        ab3 =cbind(ab1, ab21)
        ab3 = cbind(ab3,apply(ab3,1,function(v)v[2] - v[1]))
        colnames(ab3) = c("AUC1", "AUC2", "diff")
       if(verbose){
         print(t(ab3))
       }
        neginds = which(ab3[,3]<0)
        if(length(neginds)>0){
          warning(paste("AUC decreasing for phenotypes: ", paste(names(neginds))))
        }
    
      }
      rmsv = rmsv2
      if(verbose && !is.null(rmsv)){
        #ncol1 = ncol(rmsv)-2
        rmsv = subset(rmsv, measure==rmsv$measure[[1]])
        ab=pivot_wider(rmsv, names_from="measure", values_from="mid")
        ab1 =as.matrix(ab[,ncol(ab)]); rownames(ab1) = ab$pheno;
        print(t(ab1))
      
        #print(apply(ab, 2, quantile, na.rm=T))
      }
    }
   
    prev_i2 = prev_i1
    Wall = prev_i1$Wall
  
    jk = jk+1
   # print(paste(jk,nrow(Wall[[1]])))
  }
  ## next we need to adjust projection and set offset
 # models = models[1:length(nmes)]
#  names(models) = nmes
#  models
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
  ypred$updateYP(data, prev_i1, nonNA, flip=FALSE)
  rmsv=(ypred$calcRMSV(self$y, nonNA,     flip=FALSE))
  rmsv
},
 makeNextModel=function(models,prev_i2, b_i_name, phens, k,beam, ypred=NULL, project=TRUE,  useoffset=TRUE, 
                        logpthresh = -5,
                        Wall = prev_i2$Wall,
                        CHECK=getOption("fspls.check",FALSE),
                     #   fitting = TRUE,  ## whether we are trying to fit the model or just get a pvalue
                        verbose=getOption("fspls.verbose1",FALSE)) {
    prev_i1 = models$get(prev_i2, b_i_name);
    if(!is.null(prev_i1)) return(prev_i1);
    data =self
    transform_x = self$transforms
    b_i = if(length(b_i_name)==0) c() else self$convert(b_i_name)
    prev_var = if(is.null(prev_i2)) list() else lapply(prev_i2$var_names, self$convert)
    self$updateUDVP(prev_var)
 
    explore=getOption("explore_params",TRUE)
   
    if(explore && length(b_i_name)>1 && length(transform_x[[b_i_name[3]]]$params)>1){
      nme_t1 = b_i_name[[3]]
      params1 = transform_x[[nme_t1]]$params
        b_new_proj1 = lapply(params1, function(p){
          b_i_name1 = b_i_name; b_i_name1[3] = nme_t1;          b_i_name1[4] = p
         

          b2=withCallingHandlers(
            self$calcBetaProj1(phens, k, self$convert(b_i_name1), b_i_name1, prev_var, Wall,
                               prev_i2,
                              
                               project = project, convert = FALSE,useoffset=useoffset,
                             #  fitting = fitting,
                               strict = TRUE),
            warning = function(w) {
              if (grepl("fitted probabilities numerically 0 or 1 occurred", conditionMessage(w))) {
                invokeRestart("muffleWarning")
              }
            }
          )
          b2$b_i_name = b_i_name1;
          b2
        })
        logLik=unlist(lapply(b_new_proj1, function(x) sum(unlist(x$logLiks))))
        #if(verbose)print(pvl);
        best_ind = which.max(logLik)
        b_new_proj=b_new_proj1[[best_ind]]
    
     # pvl=unlist(lapply(b_new_proj_all, function(x) .sumChisq(x$pvs)))
    #  best_ind = which.min(pvl)
    #  b_new_proj= b_new_proj_all[[best_ind]]
   
    }else{
          b_new_proj <- withCallingHandlers(
            self$calcBetaProj1(phens, k, b_i, b_i_name, prev_var, Wall,
                               #fitting = fitting, 
                             prev_i2,
                               project = project, convert = FALSE,useoffset=useoffset,
                               strict = TRUE),
            warning = function(w) {
              if (grepl("fitted probabilities numerically 0 or 1 occurred", conditionMessage(w))) {
                invokeRestart("muffleWarning")
              }
            }
          )
          b_new_proj$b_i_name = b_i_name;
    }
   
    Wall2 = b_new_proj$Wall
   # b_new_proj1 = self$calcBetaProj1(subphens,k,b_i,prev_var,  betas = betas, project=!project,convert=FALSE,    strict=TRUE, useglm=useglm, useoffset=useoffset) 
    
    betas_proj = b_new_proj$betas_proj
    betas = b_new_proj[['betas']]
    
    logLiks = b_new_proj$logLiks;
    lambdas = b_new_proj$lambdas;
    degsf = b_new_proj$degsf;
      constants_proj =b_new_proj$constants
    
    tbls = b_new_proj$tbls
    
   
    mean_x = self$mean__x(b_i) #[[b_i[[1]]]][b_i[2]]
    
    prev_i1=stateObj$new(phens,data, betas_proj,betas,constants_proj, tbls, logLiks,degsf,lambdas,
                         prev_i = prev_i2 , b_i =  self$convert(b_new_proj$b_i_name),
                         b_i_name=b_new_proj$b_i_name, mean_x = mean_x, Wall = Wall2,
                         useoffset=useoffset)
 
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
  
  models$set(prev_i1, beam);
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
  df0 = lapply(nmei, function(nmei1){
    y2=self$y[[nmei1]]
    y3 = attr(y2,"factor")
    if(length(y3)==0) y3 = y2
    as.matrix(y3)
    #y3
  })
  
  df = data.frame(df0)
  names(df) = unlist(lapply(nmei, function(nmei1){
    y2=self$y[[nmei1]]
    y3 = attr(y2,"factor")
    if(length(y3)==0) y3 = y2
    paste(nmei,colnames(y3),sep=".")
  }))
  variables = vars_all1[[1]]$variables
  names(variables)=NULL
  vars = unlist(variables, recursive=FALSE)
  vars = vars[!duplicated(names(vars))]
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
   x[which(na_x==1)] = NA
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
    df3=df2 |> tibble::add_column(y=df_k,.name_repair="minimal") |> pivot_longer(names(df2)) |> separate("name",sep="__", into=into)
  }),addName="pheno")
  if(names(phensi)[[1]]!="gaussian"){
    df4$y = factor(df4$y)
  }
  df4 = subset(df4, !is.na(y))
  
  df4
},

ypred=function(phens1, liab=TRUE){
  family = unlist(lapply(names(phens1), function(str)getOption("fspls.family",strsplit(str,"\\.")[[1]][1])))
  ypr = ypredObj$new(self,self$phensi(phens1),family=family,liab =liab)
  ypr
},# inverse_func_strs = jsonlite::fromJSON(.readFlag(flags,"transform_x_inverse",'{"y":"function(y) y"}'))
#all_models_y = all_models$y; inverse_func_str = jsonlite::fromJSON(flags1$transform_x_inverse)[[1]]; self = datasAll$datas[[1]]

extractPredictions=function(all_models_full,phens, flags, 
                            transform_x = jsonlite::fromJSON(flags$transform_x),
                            liab=FALSE,merge=TRUE,
                            minv = .readFlag(flags,"min",0),
                            maxv = .readFlag(flags,"max",1e6),
                            coeffs = NULL 
                                      ){
  #all_models_y0 = all_models_full
  if(is.null(phens)) stop("phens is null")
  if(!is.null(transform_x))self$updateTransforms(transform_x   )      
  d1 = self
  self$updateLOOC(phens,flags)
  verbose=.readFlag(flags, "verbose",FALSE)
  lens = unlist(lapply(all_models_full, length))
  numvars1 = 0:(max(lens)-1)
  minv = min(max(numvars1), minv)
  numvars1 = as.character(numvars1[numvars1>=minv & numvars1<=maxv])
  
  #pheno_nmes = names(phens); names(pheno_nmes)=pheno_nmes
  if(length(all_models_full)==0) return(NULL)
  numvar_inds=1:length(numvars1); names(numvar_inds) = numvars1
  lapply(numvar_inds, function(numvar_ind){
    numvar = numvars1[[numvar_ind]]
    if(verbose)cat(paste("numvar",numvar))
  #  all_models3 = self$getAllModels(all_models_full, numvar);
    all_models3 = lapply(all_models_full, function(x) x[[numvar]])
    coeff = coeffs[[numvar]]
  
   ab=self$getPredictions(all_models3, phens, numvar,numvar_ind,  liab,
                          ypred = self$ypred(phens, liab),
                          coeff=coeff, merge = merge)
    ab
  })
   
   

},
##beams, nrep
getPredictions=function(all_models3, phens, numvar, numvar_ind,liab,
                        ypred = self$ypred(phens, liab=liab),
                        coeff= NULL,
                        merge=TRUE
                    
){
  flip=merge
nmesm =names(all_models3)
names(nmesm) = nmesm;
#nmem="full"; ij=1
res1 =  lapply(nmesm, function(nmem){
  nmem_ind = if(nmem=="full") ncol(self$looc$incl) else as.numeric(nmem)
   nonNA = self$looc$incl[,nmem_ind] ## only predict on OOF
      len = length(all_models3[[nmem]])
    beam_inds = 1:len; names(beam_inds) =names(all_models3[[nmem]])[beam_inds]
    
  res_mods = lapply(beam_inds, function(ij){
    full_model=  all_models3[[nmem]][[ij]]
    if(is.null(full_model)){
     #warning(paste("no model", "beam:",ij,"fold:" ,nmem, "numvar: ", numvar))
    return(NULL);
    }
    ypred$updateYP(self, full_model, nonNA,   flip=if(nmem=="full") FALSE else flip )
                         
    r1 =  ypred$predictions(nonNA, flip=if(nmem=="full") FALSE else flip)
    attr(r1,"betas") = full_model$betas
    attr(r1,"lambdas") = full_model$lambdas;
    attr(r1,"variable") = if(numvar=="0") rep("",4) else full_model$var_names[[as.numeric(numvar)]]
    
    r1
    #res1 = self$getRMSVInds(phens, d$nreps(), ypred)  
  })
  res_mods[unlist(lapply(res_mods, length))>0]
 
})
res1 = res1[unlist(lapply(res1, length))>0]
nmesm = names(res1)
non_full = nmesm!="full"
nmesm1 = nmesm[non_full]

if(length(nmesm1)>1 && merge){ ##merging only 
  lens_ = unlist(lapply(res1[non_full], function(x) length(x)))
  #ypreds=lapply(beam_inds, function(x) self$ypred(phens))
  nmesp = names(phens); names(nmesp)=nmesp
  nmesm_full = lapply(all_models3, function(x) attr(x,"nmes_full"))   ## this is which folds for each model select same as full model
  beam_inds = 1:max(lens_);  names(beam_inds) = beam_inds

  res_merged= lapply(beam_inds, function(mi){ ## merge for each beam ..  
   # print(mi)
    rm1 = lapply(nmesp, function(nmep){
      family=strsplit(nmep,"\\.")[[1]][1]
      y = self$y[[nmep]]
      
      if(family=="ordinal"){
        
    #    y = attr(y, "factor")
        yi = y[,1]
        levs1 = min(yi,na.rm=TRUE):max(yi,na.rm=TRUE)
        names(levs1) = levs1
        result = matrix(NA,nrow = length(y), ncol = length(levs1)-1) #ncol =length(unique(y[,1]))-1 )
        dimnames(result) = list(rownames(y),levs1[-length(levs1)])
        for(nmem in nmesm1){
        #  print(nmem)
          subinds1 = which(!self$looc$incl[,as.numeric(nmem)])
          if(mi <= length(res1[[nmem]])){
            result[subinds1,] =  res1[[nmem]][[mi]][[nmep]][subinds1,];
          }
        }
        return(data.frame(result));
      }
    
      
      ynames=1:length(phens[[nmep]]); names(ynames)=phens[[nmep]]
    ab= data.frame(lapply(ynames, function(yind){
        df0=(lapply(nmesm1,function(nmem){
#          print(nmem)
          if(mi>length(res1[[nmem]])) return(NULL)
          res1[[nmem]][[mi]][[nmep]][,yind]
        }))
        df=data.frame(df0[lapply(df0, length)>0 ])
        v2=suppressWarnings(apply(df,1,max,na.rm=T))
        v2[is.infinite(v2)]=NA
        v2
      }))
    
       ab
      
    })
    rm1[unlist(lapply(rm1, length))>0]
      })
  res1$merged =  res_merged[unlist(lapply(res_merged, length))>0]
  
}

if(!is.null(coeff)){
  nreps = names(res1); names(nreps) =nreps
  ypreds = lapply(1:max(unlist(lapply(res1, length))),function(x) self$ypred(phens, liab)) 
  res3 = lapply(nreps, function(nrep){
    res2 = res1[[nrep]]
    res2 = res2[unlist(lapply(res2, length))>0]
    coeff1 = coeff[[nrep]]
    if(is.null(coeff1)) return(res2)
   
      nonNA =rep(TRUE, self$nrow)## predict everything
      flip=FALSE; ##use everything
      incl1 = c()
     # lambdas = lapply(res2, function(x) attr(x,"lambda"))
      if(length(res2)>0){
        incl1 = 1:length(res2);
      for(ij in 1:length(res2)){
        ypreds[[ij]]$ypreds = res2[[ij]]
      }
      }
      if(length(incl1)==1){
        preds = ypreds[incl1][[1]]$ypreds
      }else{
        preds = try(ypred$updateEnsemble(ypreds[incl1], coeff1, nonNA,
                                         colnames=names(res2),liab=liab, flip=flip))
      }
      if(inherits(preds,"try-error")){
        warning("error with updating ensemble")
          
      }else{
        res2[['combined']] = preds
      }
      res2
  })
  return(res3)
}else{
  return(res1)
}

},


evaluate=function(predictions,phens,flags,
                          ypred=self$ypred(phens, liab = TRUE),
                          incl_cv_models =TRUE,
                           verbose=FALSE
                         ){
  d = self
  self$updateLOOC(phens,flags)
  nonNA = rep(TRUE, self$nrow)
 # liab = .readFlag(flags,"liab",TRUE)  ## whether to evaluate with liability , default is true
  
  numvars1 = names(predictions); names(numvars1) = numvars1
 evals_all = .merge1_new(lapply(predictions, function(pred_nvar){
   df2 =  .merge1_new(lapply(pred_nvar, function(pred_nvar_cv){
    # print("HH")
     beams = 1:length(pred_nvar_cv); names(beams) = names(pred_nvar_cv)
      .merge1_new(lapply(beams, function(beam){
      #  print(beam)
        pred_nvar_cv_beam = pred_nvar_cv[[beam]]
        pred2 = pred_nvar_cv_beam
        betas=attr(pred2,"betas"); variables =attr(pred2,"variables")
        if(nrow(pred2[[1]])<length(nonNA)) stop(paste("dimension mismatch", nrow(pred2[[1]]), self$nrow))
        ypred$ypreds = pred2
        ##note the NAs are already embedded in the prediction
       res33= ypred$calcRMSV(self$y, nonNA, flip=FALSE, betas=betas) |> tibble::add_column(beam=beam)
       res33
      }), addName ="model" )
    }),addName="cv_index")
   #cv_inds = names(pred_nvar);names(cv_inds) = cv_inds;
   
  #cv_inds = cv_inds[] 
    df3=subset(df2, cv_index!="full" & cv_index!="merged")
    
    cv_inds = names(table(df3$cv_index)); names(cv_inds) = cv_inds
    res4=lapply(cv_inds, function(cvi){
      subset(df3, cv_index==cvi)
    })
   avg= .averageResults(res4)
   
   
     rbind(df2,avg)
  }),addName="numvars")

      evals_all$numvars = as.numeric(evals_all$numvars)
  #}
    evals_all$beam[evals_all$model=="combined"]="combined"
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
projOut1=function(ik){ ## calculates W which is Vinv %*% Dinv %*% Ut %*% x
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
  if( length(var)==0 || !getOption("update_norm_each_iteration", FALSE)){
    return(self$norm[[ik]])
  }
 
  #var_st= paste(unlist(lapply(var, paste,collapse=",")),collapse=":")
  var_list1= lapply(var, function(v1){
    paste(names(self$data)[v1[1]], dimnames(self$data[[v1[1]]])[[2]][ v1[2]], sep=",")
  })

  var_st = paste(sort(unlist(var_list1)),collapse=":")
 
 
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
        res1 = apply(x_[,ii:ii1]- P %*% W[,ii:ii1],2,function(g) sum((g-mean(g))^2))
        names(res1) = colnames(x_)[ii:ii1]
        res1
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
getAngleInnerOld=function(phensi,ik,k,var,type="slow", direction=NULL,var_thresh=1e-5){
#  assoc=(type %in% c("assoc","assoc1"))
  #W = if(type %in% c("slow","assoc"))self$projOut(ik) else 
  W = self$projOut1(ik)  ## depends on x only, not y
  ##NOTE PROJOUT1 ALSO SUBTRACTS MEAN, BUT FOR PROJOUT WE HAVE TO ADJUST FOR MEAN
  ##IF WE USE PROJOUT1 then x is actually W
  #mean_x = if(type %in% c("slow","assoc")) NULL else self$mean_x[[ik]]#  x_s$mean_x
  yTr =self$train[[k]]$yTr  #this is zero in the NA positions of d,ie yTr[,d$nonNA]  should be all zero
  ## should we store this as transpose??
 
  #y = self$y
 #d = self$train[[k]]
 # nonNA = lapply(d$nonNA,t)
#  yTr[,!nonNA]=0
   
        norm = self$getNorm(W,var, ik, type) ## note the norm of previously selected variables is zero, unless maintain_norm=TRUE
        #norm_sel = abs(norm[unlist(lapply(var, function(v) if(v[1]==ik) v[2] else NULL))])
     #   if(length(which(norm_sel>var_thresh))>0 ) warning(" not projecting out properly .. possibly due to correlated vars !")
      
        to_rem1 =  unlist(lapply(var, function(v1)v1[2]))[unlist(lapply(var, function(v1)v1[1]))==ik]
        to_rem = which(norm>-var_thresh)
        P = self$UDVP$P #depends on x only and var
        #ii = names(phensi)[[1]]; nmes_prod = names(products)[[1]]; nmes_prod1 = names(products[[nmes_prod]])[[1]]
        angles1 = lapply(names(phensi), function(ii){
          
                 nme_i = ii
          phensi1 = phensi[[ii]]
          products =self$train[[k]]$product(ik,ii,phensi1);
          nmes_products = names(products); names(nmes_products) = nmes_products
          angle=lapply(nmes_products, function(nmes_prod){
            #print(nmes_prod)
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
              if(length(to_rem1)>0 && max(abs(apply(product[,to_rem1,drop=F],2,sum, na.rm=TRUE)))>var_thresh){
                warning("may not be projecting properly")
              }
             # dimnames(product) = dimnames(self$train$products[[ik]][[ii]])  
            }
           # direction=self$direction
            angle_1= t((product[]))/(norm)
            if(ncol(angle_1) != length(norm) || colnames(angle_1)[1]!=names(norm)[1]) angle_1 = t(angle_1)
            direction1 = direction[[ii]]
            mi1 = match( rownames(angle_1),names(direction1))
            direction1 = direction1[mi1];
           # print(mi1)
           # if(length(direction1) !=nrow(angle_1)) stop("problem")
         #   message(unlist(direction1));
          #  print("direction");
          #  print(unlist(direction1));
          #  print(unlist(rownames(angle_1)))
            for(jk in 1:nrow(angle_1)){
              dir1  = direction1[jk]
           #   print(jk); print(direction1[[jk]])
              if(dir1=="+"){
                angle_1[jk,angle_1[jk,]>0]=999
              }else if( dir1=="-"){
                angle_1[jk,] = -1*angle_1[jk,]
                angle_1[jk,angle_1[jk,]>0]=999
              }else{
                angle_1[jk,] = -1*abs(angle_1[jk,,drop=FALSE])
              }
            }
            if(length(to_rem)>0){
                angle_1[,to_rem]=999  #after we project out the projected out columns have zero norm
            }
            angle_1
          })
          })
         angle
        })
     
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
getDirection=function(inds1){
# 
  pos_inds = unique(inds1$pos_inds);
  neg_inds = unique(inds1$neg_inds);
  dir = inds1$direction
  direction = 
        lapply(self$y, function(y1){
          mi1=match(pos_inds,colnames(y1));
          mi2 = match( neg_inds,colnames(y1))
        
          if(length(pos_inds)+length(neg_inds)>0 && length(which(!is.na(mi1)))+length(which(!is.na(mi2)))==0){
            pos_inds1 = unique(unlist(lapply(pos_inds, function(pi1) grep(paste0(pi1,"$"),colnames(y1),v=TRUE)), recursive=TRUE))
            neg_inds1 = unique(unlist(lapply(neg_inds, function(pi1) grep(paste0(pi1,"$"),colnames(y1),v=TRUE)), recursive=TRUE))
            mi1=match(pos_inds1,colnames(y1));
            mi2 = match( neg_inds1,colnames(y1))
            
          }
          if(length(which(duplicated(c(mi1,mi2)) ) )>0) stop("inconsistent")
          
          res = rep("", ncol(y1))
          names(res)=colnames(y1)
          res[mi1] = "+";
          res[mi2] = "-";
          as.list(res)
        })
  direction    

 # dir1 = rep(list(direction),times)
#  names(dir1) = 1:times
#  dir1
},
##incl is which of the layers to include
getAngles1=function(subphens,varnames,incl, k=1,type="slow1"){
  var=lapply(varnames, self$convert)
  var = var[unlist(lapply(var, length))>0]
  phensi = self$phensi(subphens)
  names(phensi)=names(subphens)
  self$getAngles(phensi,var ,incl, k=k, type=type)
},
  getAngles=function(phensi, var,incl, k=1,type="slow1"){  ## type can be fast, slow, assoc, slow1 .. fast gives wrong results
    run_separate=TRUE
    self$updateUDVP(var)
    ntrans=attr(self$data,"ntrans") 
    if(is.null(ntrans)) ntrans=1
    #d = self$train
    nmesd = names(self$norm)
    angles_d= vector('list', length(self$data)) ## need angle object
    names(angles_d) = names(self$data)
    #incl$direction
    for(ik in 1:length(self$norm)){
    #  (ik)
     # print(ik)
      
      if(names(self$norm)[[ik]] %in% incl$types){
          angles_d[[ik]] =  tryCatch({
            self$getAngleInnerOld(phensi,ik, k,var,type, incl$direction)
          },error=function(errw){
            message(errw)
            NULL
          })
            
        
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
  transforms =  try(.convertToTransform(transform_x))
  if(inherits(transforms,"try-error")){
    print(transforms)
    print(transform_x);
    stop(" could not update transforms") 
  }
  update_trans = toJSON(self$transforms)!= toJSON(transforms)
  self$transforms = transforms
  update_trans
},
means_y=function(k){
  list(self$train[[k]]$means_y0, self$train[[k]]$means_y1)
},

counts_y=function(k){
  self$train[[k]]$counts_y
},
updateTrain =function(k, phens,
                      means_y = self$means_y(k), 
                       force=FALSE){
#  for(k in 1:length(self$train)){
    #if(verbose) cat(paste("update",k))
  
    self$train[[k]]$update(self,phens,  means_y = means_y, force=force)
 # }
},
## gets ready for training - updates train, prev looc
initTrain=function(phens,flags, transform_x, verbose=FALSE, force=FALSE){ ## this updates the reps and train  ## called after updateLOOC
 update_trans = self$updateTransform(transform_x);
  nrep = ncol(self$looc$incl)
  if(verbose) cat("update train")
  incls = jsonlite::fromJSON(.readFlag(flags,'data_types',"{}")) 
  if(length(incls) == 0 )incls = list("all"=names(self$data))
  incls_all = unique(unlist(incls))
  if(force || length(self$train)!=nrep  ||  is.null(self$train[[1]]) || update_trans){ # || toJSON(self$train[[1]]$func_str)!=toJSON(transform_x)){
    self$train = lapply(1:nrep, function(k)trainObj$new( self$y,self$weights, self$looc$incl[,k] , incls_all, self$transforms,family=self$family)) #lapply(1:ncol,function(k)
  }
  if(!is.null(self$subset)){
    ## apply subset via the looc object to avoid subsetting big matrix
    self$looc$incl[!self$subset,] = rep(FALSE, ncol(self$looc$incl))
  }
 # within=TRUE
 
 # reweight=.readFlag(flags,"reweight",FALSE)
#  if(reweight){
#    self$updateWeights(phens)
#  }
},
updateLOOC=function(phens,flags,varn=c(), force=FALSE, verbose=FALSE){
  seed=.readFlag(flags,"seed",42)
  #incl = incls_all
  batch=.readFlag(flags, "batchsize",0)

  nrep = .readFlag(flags,"nfold",if(batch>0) 0 else 1)
  if(batch>0 && nrep>0) warning("only one of batch or nrep greater than zero")
  nrows=nrow(self$data[[1]])
  if(!force && !is.null(self$looc) && nrep ==self$looc$nfold && batch == self$looc$batch && seed == self$looc$seed && nrows ==self$looc$nrows){
    if(verbose)cat("no need to update, although we should probably check if the phens changed")
    return(NULL)
  }
  pheno_balance = if(.readFlag(flags,"pheno_balance",FALSE)) unlist(phens) else NULL
  self$looc=loocObj$new(self, nfold = nrep, batch = batch,
                        folds = self$folds, 
                        incl_full =TRUE,
                        seed = seed,
                        nrow = nrows,
                        pheno_balance = pheno_balance)
  if(nrep==1 && length(self$folds)==0){
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
    self$prev[[k]] = stateObj$new(phensi, self,NULL,NULL,NULL,NULL, NULL,NULL,NULL, var=var, varnames=varn, Wall = Wall)
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
   self$fams_all = .lapply_nme(self$y, function(fam1){
     family = strsplit(fam1,"\\.")[[1]][1];
     fams = NULL;
     if(family=='binomial' || family=="multinomial" || family=="ordinal"){
       lambs= getOption("binomial.link", 1.0);names(lambs)=lambs
       fams = lapply(lambs, function(lambda) if(lambda==1) "binomial" else binomial(link = ao2(lambda)))
     }else if(family=="gaussian"){
       lambs= getOption("gaussian.link", 0.0); names(lambs)=lambs
       fams = lapply(lambs, function(lambda)if(lambda==0) "gaussian" else gaussian(link=sinh_asinh_link(lambda)))
       
     } 
     fams;
   })
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

  initialize=function( cohort,  db_name,dbDir,flags,folds = NULL,
                      incl_full=TRUE,seed = 42, memDir = NULL) { ## mem_dirp is for saving scores
  
    #self$direction = .readFlag(flags,'direction','none');
    self$db_name=db_name;
    self$folds = folds; 
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
  self$fams_all = NULL;
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


