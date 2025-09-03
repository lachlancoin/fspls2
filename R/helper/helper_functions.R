.renameModels<-function(eval,len=3){
    models = lapply(eval$model, function(x)strsplit(x,";")[[1]])
    models1 =unlist(models)
    models1 = models1[!duplicated(models1)]
    m1=lapply(models1,function(x) {
      y=strsplit(x,"\\.")[[1]]
      c(y[[1]], paste(y[-1],collapse="."))
    })
    types=(unlist(lapply(m1, function(m1)m1[[1]])))
    types1 = unique(types)
    names(types1)=lapply(types1,function(str)substr(str,1,min(len, nchar(str))))
    models2 = unlist(lapply(types1, function(typ){
      m2 = lapply(m1[which(types==typ)], function(x) paste(x,collapse="."))
      names(m2) = 1:length(m2)
      m2
    }),rec=T)
    #eval$model1 = eval$model
    eval$model= (unlist(lapply(models, function(m){
      paste(names(models2)[match(m, models2)],collapse=";")
    })))
    attr(eval,"translate")=models2
    eval
}

.avg<-function(eval0){
  nme_cols1 = c("data","subpheno","measure","pheno","trainedOn","transform_y","pheno_group","numvars","cv","family")
  rem_cols = names(eval0)[!(names(eval0) %in% nme_cols1)]
  nme_cols2 = c(nme_cols1, rem_cols)
  mi = match(nme_cols2, names(eval0))
eval1 = eval0[,mi] 
##not working

eval2 = unite(head(eval1),"data:subpheno", remove=F)
length(unique(eval1$`data:family`))
}

.calcEval1<-function(eval0, rename=T, len = 3){ #c("trainedOn","measure","subpheno")
  eval0_1 = unite(eval0, "cohort_measure_pheno_trained", "data","subpheno","measure","pheno","trainedOn","transform_y",remove=F)
  eval0_avg = subset(eval0_1, model=="avg") %>% tibble::add_column("fullmodel"="avg")
  eval =  subset(eval0_1, model!="avg")
  if(rename)eval=.renameModels(eval, len=len)
  cvs = unique(eval$cv)
  cohorts = unique(eval$cohort_measure_pheno_trained)
  eval2 = .merge1_new(lapply(cvs, function(cv1){
    .merge1_new(lapply(cohorts, function(cm){
    eval1 = subset(eval, cv==cv1 & cohort_measure_pheno_trained==cm)
    if(nrow(eval1)==0) return(NULL)
    .merge1_new(lapply(1:nrow(eval1), function(i){
      inds1 = grep(paste0("^",eval1$model[[i]]), eval1$model)
      inds2 = which(eval1$numvars[inds1]==max(eval1$numvars[inds1]))
      do.call(rbind, replicate(length(inds1[inds2]), eval1[i,], simplify=FALSE))  %>% tibble::add_column(fullmodel = eval1$model[inds1[inds2]])
    }))
    }))
  }))
  eval0_avg$cv="CV=avg"
  eval2$cv = paste("CV=",eval2$cv)
  eval21 = rbind(eval2,eval0_avg)
  eval21$isfull = paste("FULL=",eval21$isfull);
  eval3 = unite(eval21,"cv_full", "cv","isfull",remove=F,sep=" ")
  eval3$cv_full[eval3$cv=="CV=avg"]=eval3$cv[eval3$cv=="CV=avg"]
  eval3
}

.modify<-function(eval3, shape_color,
                  shape_color_nme ){
  if(!(shape_color_nme %in% names(eval3))){
      eval3 = eval3 %>% tibble::add_column(shape_color_nme = apply(eval3[,names(eval3) %in% shape_color,drop=F],1,paste, collapse=" "))
        names(eval3) = gsub("shape_color_nme", shape_color_nme, names(eval3))
  }
  eval3
}

.plotEval2<-function(eval1,...){
  if(is.null(eval1[['experiment_id']])){
    return(.plotEval1(eval1,...))
  }
  expt1 = unique(eval1$experiment_id)
  names(expt1) = expt1
  lapply(expt1,function(expt) {
    eval3 = subset(eval1, experiment_id==expt)
    .plotEval1(eval3,...)
  })
}

#      maxn = max(eval2$nsamps)
  #head(eval4)#pivot_wider(eval2, names_from = c("cv", "fullmodel", "numvars"), values_from ="value")
.plotEval1<-function(eval3,
           shape_color=c("pheno","subpheno"),linetype="fullmodel",showranges=T,
           txtsize=1,logy=F,legend=F,sep_by="",scales="free",point=T,line=T,
           grid0 = c("cohort","measure"),grid1 = "cv_full",title="", title1=""
          ){
 
  linetype_nme = paste(linetype, collapse="_");
  grid0_nme = paste(grid0, collapse="_")
  grid1_nme = paste(grid1,collapse="_");
  shape_color_nme = paste(shape_color,collapse="_")
  sep_by_nme = "sep_by"
  eval3 = .modify(eval3, shape_color, shape_color_nme)
  eval3 = .modify(eval3, linetype, linetype_nme)
  eval3 = .modify(eval3, sep_by, sep_by_nme)
  eval3 = .modify(eval3,grid0, grid0_nme)
  eval3 = .modify(eval3,grid1, grid1_nme)
  eval2 = eval3
  
  phenos = unique(eval2$sep_by)
  names(phenos)=phenos
  
  
  subphens = table(eval2$subpheno)
  subphens = subphens[order(as.numeric(unlist(lapply(names(subphens), function(str)strsplit(str,"\\|")[[1]][1]))))]
  eval2$subpheno = factor(eval2$subpheno, levels = names(subphens))
  eval2$numvars = as.numeric(eval2$numvars)
# eval2$isfull = (eval2$isfull+1)/2.0
  ggps=lapply(phenos, function(ph){ 
    eval5 = subset(eval2, sep_by==ph & !is.na(mid))
    ph3 = paste(sort(unique(apply(eval5[,names(eval5) %in% title1, drop=F],1,paste,collapse=","))), collapse=" ")
  ggp<-ggplot(eval5);
  if(point) ggp<-ggp+geom_point(aes_string(x="numvars", y="mid",  shape=shape_color_nme,size="nsamps", color=shape_color_nme))
  if(line)  ggp<-ggp+geom_line(aes_string(x="numvars", y="mid", linetype=linetype_nme,  color=shape_color_nme)) +ggtitle(paste(title,ph, ph3))
  if(showranges){ ## geom_ribbon vs geom_errorbar
    
   ggp<-ggp+ geom_ribbon(aes_string(x = "numvars", ymin="low", ymax="high",linetype=linetype_nme,color=shape_color_nme, fill = shape_color_nme ), alpha = 0.1)
  }
  if(nchar(grid0_nme)>0){
    if(nchar(grid1_nme)>0){
       ggp<-ggp+facet_grid(paste(grid0_nme, grid1_nme,sep="~"),scales=scales)
    }else{
      ggp<-ggp+facet_wrap(grid0_nme, scales=scales)
    }
  }
  legend_position=if(legend) "bottom" else "none";
  ggp<- ggp+ theme(legend.position = legend_position,legend.title = element_text(size = txtsize));#+theme(,    legend.text = element_text(size = 3))
  if(logy)ggp<-ggp+ scale_y_log10() 
  ggp
  })
  ggps
}
.randomize<-function(y){
  y[sample.int(length(y))]
}

.getRandomFuncs<-function(n){ ## although these are same, every invocation will give different results
  inds = 1:n
  names(inds)=paste0("rand_",inds);
  lapply(inds, function(i){
   c( "function(y) .randomize(y)","function(y) .randomize(y)")
  })
}

##this gets transformations for x variable
.getTransformFuncs<-function(pows,offset=1e-10){
  names(pows) = paste("pow", round(pows,2),sep="_")
  transf=lapply(pows, function(pow){ ## 1e-10 avoids problems with zeros
    paste0("function(x) sign(x+",offset,") * abs(x+",offset,")^",pow)#,paste0("function(x) sign(x) * abs(x)^",1/pow))
  })
  toJSON(transf)
}

mergeSparseMatrices<-function(m1,m2, by="row"){
  m1 <- as(m1, "TsparseMatrix")
  m2 <- as(m2, "TsparseMatrix")
  if(by=="row"){
  if("x" %in% slotNames(m1)){
    mat1 <- sparseMatrix(i = 1+c(m1@i, m2@i+nrow(m1)),
                         j = 1+c(m1@j, m2@j),
                         x = c(m1@x, m2@x))
  }else{
    mat1 <- sparseMatrix(i = 1+c(m1@i, m2@i+nrow(m1)),
                         j = 1+c(m1@j, m2@j))
  }
  }else{
    if("x" %in% slotNames(m1)){
      mat1 <- sparseMatrix(j = 1+c(m1@j, m2@j+ncol(m1)),
                           i = 1+c(m1@i, m2@i),
                           x = c(m1@x, m2@x))
    }else{
      mat1 <- sparseMatrix(j = 1+c(m1@j, m2@j+ncol(m1)),
                           i = 1+c(m1@i, m2@i))
    }
  }
  mat1
}

expandSparseMatrix<-function(counts, n,  vec, by="row"){
  mat1=replicate(n,vec,simplify = F)
  if(by=="row"){
     mat =do.call(rbind,mat1 )
  }else{
    mat = do.call(cbind,mat1)
  }
 mergeSparseMatrices(counts, Matrix(mat,sparse=T),by=by)
}

getSparseSubMatrix<-function(counts, inds,by='col'){
  colInds = inds
  rowInds = inds
  if(by=='col' && ncol(counts)==length(inds)){
    #check if we even need submatrix
     if(max(apply(cbind(colInds, 1:ncol(counts)),1,function(v) abs(v[2]-v[1])))==0) return(counts);
  }
  if(by=='row' && nrow(counts)==length(inds)){
    #check if we even need submatrix
    if(max(apply(cbind(rowInds, 1:nrow(counts)),1,function(v) abs(v[2]-v[1])))==0) return(counts);
  }
  if(by=='col'){
        colinds1 = colInds -1
        X2 <- as(counts, "TsparseMatrix")
        subinds = which(X2@j %in% colinds1)
        mi = match( X2@j[subinds],colinds1)
        dimnames = list(dimnames(X2)[[1]], dimnames(X2)[[2]][colInds])
        dims = c(dim(X2)[[1]], length(colInds))
        
        if("x" %in% slotNames(X2)){
          mat1 <- sparseMatrix(i = X2@i[subinds]+1, 
                             j = mi,
                           x = X2@x[subinds],
                           dims = dims, 
                           dimnames = dimnames)
        }else{
          mat1 <- sparseMatrix(i = X2@i[subinds]+1, 
                               j = mi,
                               dims = dims, 
                               dimnames = dimnames)
                           
        }
         
  }else{
    rowinds1 = rowInds -1
    X2 <- as(counts, "TsparseMatrix")
    subinds = which(X2@i %in% rowinds1)
    mi = match( X2@i[subinds],rowinds1)
    dimnames = list(dimnames(X2)[[1]][rowInds], dimnames(X2)[[2]])
    dims = c(length(rowInds), dim(X2)[[2]])
    if("x" %in% slotNames(X2)){
      mat1 <- sparseMatrix(j = X2@j[subinds]+1, 
                           i = mi,
                           x = X2@x[subinds],
                           dims = dims, 
                           dimnames = dimnames)
    }else{
      mat1 <- sparseMatrix(j = X2@j[subinds]+1, 
                           i = mi,
                           dims  =dims, 
                           dimnames = dimnames)
      
    }
  }
    mat1
}
sparse_norm<-function(X){ 
  
 nrow = nrow(X)
  X<-as(X,"CsparseMatrix")
  res<-split(X@x, findInterval(seq_len(nnzero(X)), X@p, left.open=TRUE))
  aa=unlist(vapply(res,
                   function(x) {
                     mean = sum(x)/nrow
                     nzero = nrow -length(x)
                     sqrt(sum((x - mean)^2)+ nzero *(mean^2))
                   },
                   FUN.VALUE = c(1)))
  res1 = rep(0, ncol(X))
  res1[as.numeric(names(aa))] = aa
  names(res1) = colnames(X)
  res1
}
sparse_variance<-function(X){ 
  nrow = nrow(X)
  X<-as(X,"CsparseMatrix")
  res<-split(X@x, findInterval(seq_len(nnzero(X)), X@p, left.open=TRUE))
  aa=unlist(vapply(res,
                   function(x) {
                     mean = sum(x)/nrow
                     nzero = nrow -length(x)
                     (sum((x - mean)^2)+ nzero *(mean^2))/nrow
                   },
                   FUN.VALUE = c(1)))
  res1 = rep(0, ncol(X))
  res1[as.numeric(names(aa))] = aa
  names(res1) = colnames(X)
  res1
}
sparse_levs<-function(X){ 
  nrow = nrow(X)
  X<-as(X,"CsparseMatrix")
  res<-split(X@x, findInterval(seq_len(nnzero(X)), X@p, left.open=TRUE))
  aa=unlist(vapply(res,
                   function(x) {
                    length(unique(x))
                   },
                   FUN.VALUE = c(1)))
  res1 = rep(0, ncol(X))
  res1[as.numeric(names(aa))] = aa
  names(res1) = colnames(X)
  res1
}




.makeMultiClass<-function(y){
  if(ncol(y)>1) stop("only one col")
  y0 = y[,1]
  t = table(y0)
  t = t[order(names(t))]
  if(length(t)==2){
    nme = names(t)
    return(as.matrix(as.numeric(factor(y0, levels = sort(nme)))-1,dims = dims(y)))
  }
  df=data.frame(unlist(lapply(1:(length(t)-1), function(k){
      base = names(t)[k]
      nme = names(t)[(k+1) : length(t)]
      names(nme) = paste(nme,".",base,sep="")
      df=lapply(nme, function(n){
        y1 = rep(NA, length(y))
        y1[y==base]=0
        y1[y==n]=1
        y1
      })
  }),rec=F))
  
  dimnames(df)[[1]] = dimnames(y)[[1]]
  df
}

.sumChisq<-function(pvi){
  if(length(pvi)==0) return(0)
  pvi1=unlist(pvi)
  pchisq(sum(qchisq(pvi1,df=1, lower.tail=F, log.p=T)), df=length(pvi1), lower.tail=F, log.p=T)
}

.readFlag<-function(flags,key, default){
  res = flags[[key]]
  if(is.null(res)) return(default)
  res
}
.mergeC<-function(v){
  res = c()
  if(length(v)>0){
    for(i in 1:length(v)) res = c(res,v[[i]])
  }
  res
}
isbigmatrix<-function(x){
  typeof(x)!="S4"
}

.fixBeforeMerge<-function(t){
  nmes_all = unique(unlist(lapply(t, names)))
  lapply(t, function(t1){
    mi1=match(nmes_all,names(t1))
    na_i = which(is.na(mi1))
    if(length(na_i)>0){
      t1 = cbind(t1,array(NA, dim = c(nrow(t1), ncol=length(na_i)), dimnames =list(NULL, nmes_all[na_i])))
      mi1=match(nmes_all,names(t1))
    }
  
    t1[,mi1,drop=F]
  })
}

.merge1_new<-function(t,num_cols = c(), addName=NULL){
  if(length(t)==0) return(NULL)
  t = t[!unlist(lapply(t, is.null))]
  t = t[unlist(lapply(t,nrow))>0]
  
  if(length(t)==0) return(NULL)
  if(is.null(names(t))) addName=NULL
  nrows = unlist(lapply(t, nrow))
  end = cumsum(nrows)
  start = end -nrows+1
  colnames = names(t[[1]])
  ncol = ncol(t[[1]])
  add = !is.null(addName)
  if(!is.null(addName)) colnames = c(colnames, addName)
  df = data.frame(array(dim = c(sum(nrows), length(colnames)), dimnames = list(NULL, colnames)))
  for(k in 1:length(t)){
    range = start[[k]] : end[[k]]
    df[range,1:ncol] = t[[k]]
    if(add){
      df[range,ncol+1] = names(t)[[k]]
    }
  }
  df
}



.merge1<-function(t,num_cols = c(),uniq_cols=c(), addName=NULL, rowNames=F){
  t=t[unlist(lapply(t, nrow))>0]
  t=t[unlist(lapply(t, length))>0]
  if(length(t)==0) return(NULL)
  t1 = t[[1]]
  if(!is.null(addName)){
    nme = rep(names(t)[[1]], nrow(t1))
    t1 = cbind(t1, nme)
  }
  #if(length(t)==1) return(t1)
  if(length(t)>1){
    for(i in 2:length(t)){
      t_i=t[[i]]
      if(!is.null(addName)){
        nme = rep(names(t)[[i]], nrow(t_i))
        t_i = cbind(t_i, nme)
      }
      t1 = rbind(t1,t_i)
    }
  }
  if(rowNames){
    dimnames(t1)[[1]] = .mergeC(lapply(t, function(t_)dimnames(t_)[[1]]))
  }else{
    dimnames(t1)[[1]]=1:nrow(t1)
  }
  t2 = data.frame(t1)
  names(t2) = dimnames(t1)[[2]]
  if(!is.null(addName)){
    names(t2)[ncol(t2)]=addName
  }
  
  if(length(uniq_cols)>0){
    facts = factor(apply(t2[,which(names(t2) %in% uniq_cols),drop=F],1,paste, collapse="."))
    t2 = t2[!duplicated(facts),,drop=F]
  }
  num_cols = num_cols[which(num_cols %in% names(t2))]
 # lev_cols =names(t2)[!( names(t2)%in% num_cols)]
  for(j in num_cols){
    t2[[j]] = as.numeric(t2[[j]])
  }
#  for(j in lev_cols){
 #   t2[[j]] = factor(t2[[j]])
    
#  }
  
  t2
}



whichpart1<-function(angle, n=10, one_for_each=F, return_scores=F){
  nulli=unlist(lapply(angle, is.null))
  # names(angle) = 1:length(angle)
  wp1=lapply(angle[!nulli], function(x1){
#   x1 = apply(x,2,combine_func)
    #if(!is.null(dim(x1))) {
     # x1 = apply(x1,2,sum)
      wp_ = whichpart(x1,n)
      return(x1[wp_])
#    }else{
#    wp_ = whichpart(x1,n)
#    wp_1 = x1[1,wp_]
#    wp_1
#    }
  })
  if(return_scores) return(wp1)
  wp11 = unlist(wp1)
  wp2 = if(one_for_each)  1:length(wp11) else whichpart(wp11,n=n)
  #names(wp2)
  a1 = lapply(names(wp11[wp2]),function(v1){
    split1 = strsplit(v1,"\\.")[[1]]
    c(split1[1], paste(split1[-1], collapse="."))
  })
  a1
 # lapply(a1, function(a) c(which(names(angle)==a[1]),as.numeric(a[2])))
}

whichpart <- function(x, n=10) {
  nonNA = which(!is.na(x))
  if(n==1) return(which(x==min(x[nonNA]))[1])
  if(length(nonNA)>n) {
    #nx <- length(x)
    #nacnt = nx -length(nonNA)
    # p <- nx-n-nacnt
    xp <- sort(x, partial=n)[n]
    inds = which(x <= xp)
    if(length(inds)<n){
      #print(xp)
      #print(paste(nx,nacnt,p,n))
      ##print(head(xp))
      ##print(head(x))
      #print(inds)
      stop("problem in which part")
    }
    inds = inds[1:n]
  }else{
    inds = nonNA	
  }
  inds[order(x[inds])]
  
}



.keepUniq<-function(vars, removeTrans=T){
  if(removeTrans){
    vals1 = unlist(lapply(names(vars), function(st)paste(sort(
      unlist(lapply(strsplit(st,",")[[1]], function(st1) paste(strsplit(st1,"\\.")[[1]][-1], collapse="."))))
      ,collapse=",")))
  }else{
    vals1 = unlist(lapply(names(vars), function(st)paste(sort( strsplit(st,",")[[1]])
                                                         ,collapse=",")))
  }
  vars[!duplicated(vals1)]
}


.findMinRMSV<-function(rmsv_, mult=rep(1, length(levels(rmsv_$data))), fspls.sum=T){
  
  if(fspls.sum){
    beams = levels(rmsv_$beam)
    names(beams)=beams
    rmsv_allp = lapply(beams, function(p1){
      rmsv_sub=subset(rmsv_, beam==p1)
      rmsv_sub1  =pivot_wider(rmsv_sub,id_cols="phens", values_from="value", names_from="data")
      rmsv=as.matrix(rmsv_sub1[,-1,drop=F])
      dimnames(rmsv) [[1]] = rmsv_sub1$phens
      rmsv[is.na(rmsv)]=0 ## set NA to zero for multiplication
      #rmsv1 =   (rmsv %*% mult)[,1]
      #attr(rmsv1,"phen")=p1
      #sum(rmsv1)
      rmsv[,1]
    })
    ord=order(unlist(rmsv_allp))
    r1=unlist(rmsv_allp[ord])
    return(r1)
    #rmsv_allp[which(min(unlist(rmsv_allp)))]
  }else{
    phens1 = levels(rmsv_$phens)
    names(phens1) = phens1
    rmsv_allp = lapply(phens1, function(p1){
      rmsv_sub=subset(rmsv_, phens==p1)
      rmsv_sub1  =pivot_wider(rmsv_sub,id_cols="beam", values_from="value", names_from="data")
      rmsv=as.matrix(rmsv_sub1[,-1,drop=F])
      dimnames(rmsv) [[1]] = rmsv_sub1$beam
      rmsv1 =   (rmsv %*% mult)[,1]
      attr(rmsv1,"phen")=p1whichpart
      rmsv1
    })
    rmsv_allp[[which.min(unlist(lapply(rmsv_allp, function(v) min(v))))[1]]]
  }
}
.replot_dist<-function(rdsfile,  min=1,max=4, outpdf1=NULL, typed="datas"){
  fi = grep("rds",dir(rdsfile, full=T),v=T)
  if(length(fi)<min) return(NULL)
  names(fi) = unlist(lapply(fi, function(st)gsub(".rds","",rev(strsplit(st,"/")[[1]])[1])))
  fi = fi[order(as.numeric(names(fi)))]
  fi = fi[as.numeric(names(fi))>0]
  if(length(fi)<min) return(NULL)
  fi = fi[1:min(length(fi),max)]
  
  df0 = .merge1(lapply(fi, function(fi1){
    ar = readRDS(fi1)
    .merge1(lapply(ar$datas[[typed]], function(d){
      y1=d$y[,1]
      yp =d$ypreds_all$ypreds[[1]][[1]][,1]
      roc1 =roc(y1,yp,plot=F)
      df2=data.frame(cbind(roc1$sensitivities, roc1$specificities))
      names(df2) = c("sens","spec")
      nme = names(ar$datas[[1]]$train[[1]]$prev)[1]
      print(paste(nme, roc1$auc))
      nme=gsub(",","\n",nme)
      title=paste(nme, collapse="\n")
      title=nme
     
      cbind(df2,title)
      #  if(type=="area") return(plotAreas(yp,y1,title=nme)) else return(ggroc(roc1))
    }), num_cols = c("sens","spec"), addName="lineage")
  }), num_cols = c("sens","spec"), addName="index")
  
  
  df = .merge1(lapply(fi, function(fi1){
    print(fi1)
    ar = readRDS(fi1)
    .merge1(lapply(ar[[typed]], function(d){
      y1=d$y[,1]
      yp =d$ypreds_all$ypreds[[1]][[1]][,1]
      df2=data.frame(cbind(yp,y1))
      names(df2) = c("value","pheno")
      tab = table(df2)
      df3=data.frame(cbind(as.numeric(dimnames(tab)[[1]]), tab))
      
      names(df3) =c("knots","0","1")
     
      df4 = pivot_longer(df3,names(df3[-1]))
      
      names(df4) = c("knots", "pheno", "counts")
      df4=df4[df4$counts>0,,drop=F]
      #ggplot(df4, aes(x=knots, y=value))+geom_point()
      #ggplot(df4, aes(x=knots, y=name, size=value))+geom_point()
      #ggplot(df2, aes(x=factor(pheno), y=value))+geom_jitter(alpha = 0.9, width=0.1, size=.2)
      roc1 =roc(y1,yp,plot=F)
      nme = names(ar$datas[[1]]$train[[1]]$prev)[1]
      print(paste(nme, roc1$auc))
      nme=gsub(",","\n",nme)
      title=paste(nme, collapse="\n")
      title=nme
#      cbind(df2,title)
      cbind(df4,title)
      #  if(type=="area") return(plotAreas(yp,y1,title=nme)) else return(ggroc(roc1))
    }), num_cols = c("knots","value","counts"), addName="lineage")
  }), num_cols = c("knots","value","counts"), addName="index")
 # df$label = apply(cbind(as.character(df$pheno), round(df$value,2)),1,paste,collapse=" ")
  df = df[!is.na(df$pheno),]
  df$pheno=as.numeric(as.character(df$pheno))
  ggp1=ggplot(df, aes(x=knots, y=pheno))#+geom_jitter(alpha = 0.9, width=0.1, size=.2)+ggtitle(rdsfile)#+geom_violin(aes(x=pheno,y=value))
  ggp1<-ggp1+geom_point(aes(size=counts),color="grey")
  ggp1<-ggp1+geom_text_repel(aes(x=knots,y=pheno, label=counts ),size=2)
  if(!is.null(df0)){
  ggp1<-ggp1+geom_line(data=df0, aes(x=1-spec,y=sens), linetype='dashed',colour="grey")
  }
  #ggp1=ggplot(df, aes(x=knots, y=value, color=pheno, linetype=type))+geom_line()+geom_point(aes(size=counts))
#  ggp1= ggp1+geom_text(data=txt_df,nudge_y=0.02,
#                       inherit.aes =T, 
#                       aes(x=knots,y=value,label=label,color=pheno),size=2)
  ggp1<-ggp1+facet_grid(lineage~title)+ggtitle(rdsfile)+scale_x_continuous(limits = c(-0.1,1.1))+scale_y_continuous(limits = c(-0.1,1.1))
  ggp1
if(!is.null(outpdf1))  try(ggsave(outpdf1, plot=ggp1, width =45, height =45, units = "cm",limitsize=F))
  #attr(ggp1,"text_df")=txt_df
  invisible(ggp1)
}

.replot<-function(ar,  typed="datas"){
    #ar = readRDS(rdsfile)
    df=.merge1(lapply(ar[[typed]], function(d){
      y1=d$y[,1]
      yp =d$ypreds_all$ypreds[[1]][[1]][,1]
      roc1 =roc(y1,yp,plot=F)
      nme = names(ar$datas[[1]]$train[[1]]$prev)[1]
      print(paste(nme, roc1$auc))
      nme=gsub(",","\n",nme)
      title=paste(nme, collapse="\n")
      title=nme
      area_plot=getAreaPlot(yp, y1,nme)#, title=title)
      area_plot = addExtraLines(area_plot,yp,y1,nme)
      mat = d$train[[length(d$train)]]$prev[[1]]$betas[[1]]
      if(is.matrix(mat) ){
        beta= try(mat[nrow(mat),1])
        beta = sign(beta)*log10(abs(beta))
        area_plot1 = data.frame(rbind(c(0,beta,0,"beta","summary",nme), c(1,beta,0,"beta","summary",nme)))
        names(area_plot1) = names(area_plot)
        area_plot=rbind(area_plot, area_plot1)
      }
      area_plot
      #  if(type=="area") return(plotAreas(yp,y1,title=nme)) else return(ggroc(roc1))
    }), num_cols = c("knots","value"), addName="lineage")
  df$label = apply(cbind(as.character(df$pheno), round(df$value,2)),1,paste,collapse=" ")
  beta_inds = df$pheno=='beta'
  beta_scale = max(abs(df$value[beta_inds]))
  df$value[beta_inds] = df$value[beta_inds]/(2*beta_scale) + .5
 df
}
.replot_plot<-function(df, title="", outpdf1=NULL){
  txt_df = subset(df ,type=="summary" & knots==0)
  txt_df$knots=0.02
  txt_df$knots[txt_df$pheno=="beta"] = 0.95
  txt_df$label1 = unlist(lapply(as.character(txt_df$label), function(st)strsplit(st," ")[[1]][1]))
  df$label = as.character(df$label)
  ggp1=ggplot(df, aes(x=knots, y=value, color=pheno, linetype=type))+geom_line()+geom_point(aes(size=counts))+ggtitle(title)
  if(!is.infinite(beta_scale)){
  ggp1=ggp1+scale_y_continuous(
    "value", 
    sec.axis = sec_axis(~ (.+.5)*(2*beta_scale), name = "betas")
  )  
  }
  ggp1= ggp1+geom_text(data=txt_df,nudge_y=0.02,
                       inherit.aes =T, 
                       aes(x=knots,y=value,label=label,color=pheno),size=2)
  ggp1<-ggp1+facet_grid(lineage~drug)
  attr(ggp1,"text_df")=txt_df
  if(!is.null(outpdf1))  try(ggsave(outpdf1, plot=ggp1, width =45, height =45, units = "cm",limitsize=F))
  
  ggp2=ggplot(txt_df, aes(x=drug, y=value, color=lineage))+geom_point()+facet_grid(label1~.)
  
  list(ggp1 = ggp1,ggp2=ggp2)
}


addExtraLines<-function(df,yp,y1,title="",input=list()){
  auc = (.calcAUCW(yp,y1))[2]
  youden = (.youden(as.matrix(yp),y1)[2])
  area = attr(df,"area")
  #area =  (.areaBetween(yp,y1))[[2]]
  
  extra = c(auc, youden,area)
  nme_e = c("auc","youden","area")
  extra = data.frame(rbind(cbind(extra,  0), cbind(extra,1)))
  
  extra = cbind(c(nme_e,nme_e), extra)
  extra = cbind(extra, apply(extra[,c(1,2)],1,paste,collapse="="))
  names(extra) = c("pheno","value","knots","label")
  
  extra = cbind(extra,"summary")
  names(extra)[ncol(extra)]="type"
  extra = cbind(extra,0)
  names(extra)[ncol(extra)]="counts"
  
  df=rbind(df, extra[,match(names(df), names(extra))])
  # print(paste(i1,length(kn)))
  cbind(df,title)
}
getModels=function(models, index="full"){
  res=lapply(models, function(x){
    mi=match(index, names(x))
    if(is.na(mi)) return(NULL)
    x[mi]
  })
  res = res[lapply(res, length)>0]
}

getAreaPlot<-function(yp, y1,title = "", input = list()){
  levs = sort( unique(y1[!is.na(y1)]))
  if(length(levs)<=1) return(NULL)
  names(levs)=levs
  knots=sort(unique(unlist(lapply(levs, function(t){
    inds = y1==t
    cdf=ecdf(yp[inds])
    kn=stats::knots(cdf)
  }))))
  if(!is.null(input$range)){
    knots = c(input$range[0], knots, input$range[1])
  }else{
    knots = c(0,knots,1)
  }
  df=.merge1(lapply(levs, function(t){
    inds = y1==t
    tab1 = table(yp[inds])
    cdf=ecdf(yp[inds])
    mi2=match(knots, names(tab1))
    no_dupl = !duplicated(mi2)
    #mi2 = mi2[!duplicated(mi2)]
    tab1_col = tab1[mi2[no_dupl]]
    tab1_col[is.na(tab1_col)]=0
    df1 = data.frame(cbind(knots[no_dupl],cdf(knots[no_dupl]), tab1_col))
    names(df1)= c("knots","value","counts")
    df1
  }), addName="subpheno",num_cols=c("knots","value"))
  pw = pivot_wider(df, names_from="subpheno", id_cols="knots")
  diff=0
  for(k in 2:nrow(pw)){
    delta = (pw[k,1]-pw[k-1,1])
    diff = diff+((pw[k-1,2] + pw[k,2] - (pw[k-1,3] + pw[k,3]))/2) * delta
  }
  area=diff
  
  df = cbind(df, "cumulative")
  
  
  
  names(df)[ncol(df)] = "type"
 
 
  attr(df,"area")=area
  df
}
.combinePv<-function(pvs){
  return(min(pvs))
  chisq = qchisq(pvs, df=1, lower.tail=F)
  pchisq(sum(chisq), df = length(pvs), lower.tail=F)
}

.plotArea1<-function(predictions,family="binomial",rename=T, 
                     subset = 1:length(predictions[[1]][[1]]), #[[1]]),
                     len = 1,grid="model~pheno", max_vars = 100){
  area_p= .getAreaPlot1(predictions,family)
  .plotArea(area_p, rename=rename, len = len,grid=grid,  max_vars = max_vars,title=title)
}

.getAreaPlot1<-function(predictions0, families="binomial"){
  if(is.null(families)) families = names(predictions0[[1]][[1]][[1]])
  names(families)=families
 #j=1; model = predictions0[[2]]; train = model[[1]]; test = train[[1]]; famnme = families[[1]]; fam =test[[famnme]]; family = strsplit(famnme,"\\.")[[1]][1];  phens = dimnames(fam$y)[[2]]; names(phens)=phens
  
  area_p=.merge1_new(lapply(predictions0, function(model){  
    .merge1_new(lapply(model, function(train){
      .merge1_new(lapply(train, function(test){
#        families = names(test); names(families)=families
        .merge1_new(lapply(families, function(famnme){
          fam = test[[famnme]]
          family = strsplit(famnme,"\\.")[[1]][1]
        phens = dimnames(fam$y)[[2]]; names(phens)=phens
       # phens = phens[unlist(lapply(phens, function(p)length(grep(p_incl,p))))>0]
        phensi= 1:length(phens)
        names(phensi) = phens
        .merge1_new(lapply(phensi, function(j){
  #        indsk = 1:ncol(fam$ypred[[phen]]);
  #        names(indsk)=dimnames(fam$ypred[[phen]])[[2]]
         # .merge1_new(lapply(indsk , function(k){
            if(family=="gaussian"){
               ap  = data.frame(list(knots = fam$y[,j],value=fam$ypred[,j]), counts=0, subpheno="", type="points")
               names(ap)=c("knots","value","counts","subpheno", "type")
               ap = ap %>% tibble::add_column(sample=dimnames(fam$y)[[1]])
              
            }else if(family=="ordinal"){
              ap = getAreaPlot(fam$ypred, fam$y[,j])
        }else{
            ap = getAreaPlot(fam$ypred[,j], fam$y[,j])
          
            }
            ap
          #}),addName="subpheno")
        }), addName="pheno")
        }), addName="family")
      }), addName="test")
    }), addName="train")
    }), addName="model")
  #attr(area_p,"family")=family
  lens = unlist(vapply(area_p$model, function(x) if(x=="empty") 0 else length(strsplit(x,";")[[1]]), FUN.VALUE = c(1)))
 
   area_p%>% tibble::add_column(lens = lens)
}
.plotAreaSep=function(area_p, sep="pheno",...){
  phens=unique(area_p[[sep]]);names(phens)=phens
 lapply(phens, function(p){
   p = phens[[1]]
  area_p1=area_p[area_p[[sep]]==p,,drop=F]
  .plotArea(area_p1, ...)
 })  
}
.takeMax1<-function(a1, max_vars=100){
  ab=unite(subset(a1, lens<=max_vars),"comb","sample","CV",remove=F)
  levs = unique(ab$comb); names(levs)=levs
  ab2=.merge1_new(lapply(levs, function(l1){
    sub1 = subset(ab, comb==l1 & !is.na(knots))
    sub1[which.max(sub1$lens),,drop=F]
  }))
  #ab2$pheno = factor(ab2$pheno, labels=levels(ab$pheno))
  ab2
}

.takeMax<-function(area_p1){
  ab=unite(area_p1,"comb","sample","pheno","CV",remove=F)
  ab$pheno = factor(ab$pheno)
  levs = unique(ab$comb); names(levs)=levs
  ab2=.merge1_new(lapply(levs, function(l1){
    sub1 = subset(ab, comb==l1 & !is.na(knots))
    sub1[which.max(sub1$lens),,drop=F]
  }))
  ab2$pheno = factor(ab2$pheno, labels=levels(ab$pheno))
  ab2
}
.getTextLayer<-function(area_p1,p = c(0.1,0.9),r2=F, arrange_by_sample=F, reorder=F){
  ab = if(arrange_by_sample) unite(area_p1,"comb","sample","CV",remove=F) else unite(area_p1,"comb","lens","pheno","CV",remove=F)
  ab$pheno = factor(ab$pheno)
  ab$CV = factor(ab$CV)
  ab$sample=factor(ab$sample)
  levs = unique(ab$comb); names(levs)=levs
  ab2=.merge1_new(lapply(levs, function(l1){
    sub1 = subset(ab, comb==l1 & !is.na(knots))
    model = unique(sub1$model)
    pheno1 = unique(as.character(sub1$pheno))
    sample = unique(as.character(sub1$sample))
    if(length(model)>1)model="avg"
    if(length(sample)>1)sample="avg"
    if(length(pheno1)>1)pheno1="avg"
    corr = if(nrow(sub1)<2) NA else cor(sub1$knots, sub1$value, method="pearson", use="pairwise.complete.obs")
    label = if(r2)   paste0("r2=",round(corr^2,digits=3))   else paste0("cor=",round(corr,digits=3))
    data.frame(list(lens = sub1$lens[1],
                    model=model,
                    CV = as.character(sub1$CV[1]),
                    corr=corr,
                    sign=if(is.na(corr) || sign(corr)<0) "italic" else "bold",
                    sample=sample,
                    pheno = pheno1, knots=quantile(sub1$knots,na.rm=T, probs = p[1]), value=quantile(sub1$value, na.rm=T, probs =p[2]), 
                    label= label))
  }))
    ab2$CV = factor(ab2$CV, levels = levels(ab$CV))
   ab2$pheno = factor(ab2$pheno, levels = levels(ab$pheno))
   ab2$sample = factor(ab2$sample, levels = levels(ab$sample))
  if(arrange_by_sample){
    phens2=levels(ab2$sample);names(phens2)=phens2
    
    levs1 = names(sort(unlist(lapply(phens2, function(p2){
      sb1=subset(ab2, sample==p2)
      max(sb1$cor,na.rm=T)
    })),decr=T))
    if(reorder)ab2$sample = factor(ab2$sample, levels=levs1)
  }else{
    phens2=levels(ab2$pheno);names(phens2)=phens2
    
  levs1 = names(sort(unlist(lapply(phens2, function(p2){
    sb1=subset(ab2, pheno==p2 & CV!="No CV")
    max(sb1$cor)
  })),decr=T))
  if(reorder)ab2$pheno = factor(ab2$pheno, levels=levs1)
  }
  ab2
}
.plotArea<-function(area_p0, CV=F,alpha =0.8, maxphens = 50,p=c(0.1,0.9),
                    arrange_by_sample=F,takeMax=F,
                    rename=T, len = 1,shapes=F,code_len=3,showText=F,
                    max_samps=100,reorder=F,
                    grid=if(!CV) "model~pheno" else "lens ~pheno",  max_vars = 10,scales="free",title="", addText=F,r2=F){
  print("now plotting")
  if(!is.null(area_p0$sample) && !arrange_by_sample){
    area_p0$sample = as.character(as.numeric(factor(area_p0$sample)))
  }
  family = area_p0$family[[1]]
  area_p1=if(rename) .renameModels(area_p0, len=len) else area_p0
 
  area_p1 = subset(area_p1, lens<=max_vars)
  if(arrange_by_sample){
  # if(takeMax) area_p1 = .takeMax(area_p1) 
   if(!showText){
     pheno_code = as.numeric(area_p1$pheno) 
   }else{
    pheno_code= unlist(lapply(area_p1$pheno, substr, 1,code_len)  )
   }
   area_p1 = area_p1%>%tibble::add_column(pheno_code=pheno_code)
  }
  textLayer=NULL
  if(addText){
    textLayer=.getTextLayer(area_p1,r2=r2,p=p, arrange_by_sample=arrange_by_sample, reorder=reorder)
    if(arrange_by_sample){
      levs = levels(textLayer$sample)
      
      if(length(levs)> max_samps){
        
        textLayer=subset(textLayer, sample %in% levs[1:max_samps])
        area_p1 = subset(area_p1, sample %in% levs[1:max_samps])
        levs = levs[1:max_samps]
      }
      
      area_p1$sample = factor(area_p1$sample, levels = levs)
    }else{
      levs = levels(textLayer$pheno)
      
      if(length(levs)>maxphens ){
        
        textLayer=subset(textLayer, pheno %in% levs[1:maxphens])
        area_p1 = subset(area_p1, pheno %in% levs[1:maxphens])
        levs = levs[1:maxphens]
      }
      
      area_p1$pheno = factor(area_p1$pheno, levels = levs)
    }
    if(scales=="fixed"){
      textLayer$knots = mean(textLayer$knots);
      textLayer$value = mean(textLayer$value)
    }
    #print(levels(textLayer$sample))
    #textLayer$sign = factor(textLayer$sign)
  }
  
  
  
  #print(area_p1$subpheno)
  #if(length(unique(area_p1$subpheno))==0) area_p = area_p[,!(names(area_p) %in% "subpheno")]
  if(!is.null(area_p1[['subpheno']])){
    if(family!="multinomial"){
    area_p1$subpheno = factor(area_p1$subpheno, levels = sort(as.numeric(unique(area_p1$subpheno))))
    }else{
      area_p1$subpheno = factor(area_p1$subpheno, levels = sort((unique(area_p1$subpheno))))
      
    }
  }
   
   # if(!is.null(textLayer)) textLayer = subset(textLayer, lens<=max_vars)
    #area_p1$lens = factor(area_p1$lens)
  #  area_p1$CV = factor(area_p1$CV)
    if(arrange_by_sample){
      
      if(!takeMax){
        ggp<-ggplot(area_p1, aes(x=knots, y=value, color=lens))
        
      }else if(shapes ){
        area_p1$lens = factor(area_p1$lens)
      ggp<-ggplot(area_p1, aes(x=knots, y=value, color=pheno, shape=lens))
      }else{
        if(length(unique(area_p1$pheno))>30){
          area_p1$pheno=as.numeric(factor(area_p1$pheno))
          }
        ggp<-ggplot(area_p1, aes(x=knots, y=value, color=pheno))
        
      }
    }else if(is.null(area_p1[['subpheno']]) ||  length(unique(area_p1$subpheno))<=1){
    ggp<-ggplot(area_p1, aes(x=knots, y=value, color=model))
    
  }else{
  ggp<-ggplot(area_p1, aes(x=knots, y=value, color=subpheno))
  }
  if(family!="gaussian") ggp<-ggp+geom_line()
  
  if(family=="gaussian")ggp<-ggp+geom_abline(slope=1, intercept=0, alpha = 0.5, linetype="dashed")
  if(!is.null(area_p1$sample) ){
    if(arrange_by_sample){
    
    #if(!is.null(textLayer))    textLayer =textLayer%>%tibble::add_column(pheno_code= "avg" )
    #area_p1$pheno_code
      if(length(unique(area_p1$pheno))>30){
        ggp=ggp+geom_point( size=1,alpha = alpha)
      }else{
        ggp=ggp+geom_text_repel(aes(label=pheno_code), alpha = alpha) 
      }
    }else{
    ggp=ggp+geom_text(aes(label=sample), alpha = alpha) 
    }
  }else if(max(area_p1$counts,na.rm=T)>5) {
    ggp=ggp+geom_point(aes(size=counts), alpha = alpha) 
  }else {
    ggp=ggp+geom_point( size=2,alpha = alpha)
  }
    
    if(addText){
      ggp<-ggp+geom_text(data=textLayer, aes(x=knots, y=value, label=label, fontface=sign), inherit.aes = FALSE)
      
    }
    if(length(grep("~",grid))>0){
      ggp<-ggp+facet_grid(grid,scales=scales)
    }else{
      ggp<-ggp+facet_wrap(grid, scales=scales)
    }
  
  ggp+ggtitle(title)
}

.cntNA<-function(vec) length(which(is.na(vec)))

.rep<-function(v,n) t(matrix(rep(v,n),ncol=n))


.getRMSPrev=function(rmsv_){
  beam_levs = levels(rmsv_$beam)
  names(beam_levs) = beam_levs
  rmsv1=unlist(lapply(beam_levs, function(beam){
    ab=subset(rmsv_, beam==beam)
    sum(ab$value,na.rm=T)  #could also try min ? 
  }))
  best_beam = names(which.min(rmsv1))[1]
  rmsv2 = subset(rmsv_, beam==best_beam)
  rmsv_prev1 = rmsv2$value
  names(rmsv_prev1) = apply(rmsv2[,names(rmsv2) %in% c("subpheno","pheno","data")],1,paste,collapse=".")
  rmsv_prev1
}


.addColumn<-function(phenot,json1){
  json_nme = names(json1)
  names(json_nme) = json_nme
  res_all=lapply(json_nme,function(nme){
    disease_class = phenot[[nme]]
    res=data.frame(lapply(json1[[nme]],function(tblkk){
#      tblkk = tbl[kk,]
       dc = rep(NA, length(disease_class))
       dc[(disease_class %in% tblkk)] = disease_class[disease_class %in% tblkk]
       tbls1 = names(sort(table(dc), decr=T))
       factor(dc, levels = tbls1)
      #factor(bacviral, levels = tblkk)
    }))
    names(res) = names(json1[[nme]])
    res
  })
  res_all1 = data.frame(unlist(res_all,rec=F))
  cbind(phenot, res_all1)
}
.inferFamily<-function(y){
  unlist(lapply(y, function(y1) {
    ty = table(y1)
    if(is.factor(y1)) {
      return(if(length(ty)>2) "multinomial" else "binomial")
    }else if(length(ty)==2 && names(ty)[[1]]=="0" && names(ty)[[2]] == "1"){
     return("binomial")      
    }else if( max(apply(cbind(y1,round(y1) ),1,function(x) abs(diff(x))),na.rm=T)==0){
      return("ordinal")
    }else{
      return("gaussian")
    }
  }))
}
.getSubset<-function(phenotypes, col){
  levs = levels(phenotypes[[col]])
  names(levs) = levs
  l1=lapply(levs,function(lev){
    r1 = (phenotypes[[col]]==lev)
    r1[is.na(r1)]=FALSE
    r1
  })
  l1[!unlist(lapply(l1, function(xx) length(which(xx))))==0]
}
.subset<-function(dat, subset){
  dat1 = dat$clone(deep=T)
  dat1$subset = subset
  dat1
}


.appendPredictedPheno<-function(phenotypes, data, nmes =c("Bacterial","Viral","NonInfectious") , thresh =0.5){
  ypred_all1 = data$ypreds_all$ypreds[[1]][[1]]
  dimnames(ypred_all1)[[1]] = dimnames(data$y)[[1]]
  names(ypred_all1 ) = nmes

  disease_class_predicted = 
    factor(apply(ypred_all1,1,function(xx) names(ypred_all1)[which.max(xx)]), levels = names(ypred_all1))
  disease_class_predicted [apply(ypred_all1,1, max)<thresh] = NA
  
 print( table(cbind(disease_class_predicted),data$y[,1]))
  phenotypes[['disease_class_predicted']] = disease_class_predicted
#  return(cbind(phenotypes, disease_class_predicted))
  phenotypes   
}


.splitByPheno<-function(datas, phenotypes, nme='disease_class_predicted'){
  subsets = .getSubset(phenotypes, nme)
  datas_=unlist(lapply(subsets, function(subset){
    lapply(datas, function(d) {
      .subset(d, subset)
    })
  }),rec=F)
  c(datas, datas_)  
}
#datas1 is independent validation
.calcRMSVAll<-function(datas,datas1,cv, label, models){
  if(!is.null(datas1)){
    mi1 = match(names(datas), names(datas1))
    if(length(which(!is.na(mi1)))>0){
    datas = datas[!is.na(mi1)]
    datas1 = datas1[mi1[!is.na(mi1)]]
#    print(cbind(names(datas), names(datas1)))
    }
  }
 
#  cnts_all=unlist(lapply(models, function(m) (m$cnt)))
#  print(table(unlist(lapply(models, function(m) (names(m$prev))))))
#  if(length(models)==1){
#    return(NULL)
#  }
  numv = -1 #c(0:min(cnts_all),-1)
  names(numv) = numv
  names(numv)[which(numv<0)]="full"
  if(cv){
  rmsval = .merge1(lapply( numv, function(kk){
    #kk=0 means evaluate the full models
    rms_cv = .merge1(lapply(datas,function(d) d$getRMSVAll(numvar=if(kk<0) NULL else kk)),num_cols="value",addName="dataset")
    rms_cv
  }),num_cols="value",addName="numvars")
 
  }else{
      k = length(datas[[1]]$train)
      prevs = lapply(datas, function(d){
        prevk=models[[k]]$prev[[1]]
        d$makeModel(k,  prevk$vars)
      })
      data1_inds = 1:length(datas1); names(data1_inds) = names(datas1)## because we kept the full model here
      rmsval =.merge1(lapply( numv, function(kk){
        print(kk)
        #kk=-1 means evaluate the full models
        .merge1(

          lapply(data1_inds, function(ij) {
            ij2 = if(length(datas)==1) 1 else ij
            prev1 = prevs[ij2] #datas[[ij2]]$train[[k]]$prev
            prev2 = if(label=="discovery") prev1 else datas1[[ij]]$translate(prev1)
            datas1[[ij]]$initY()
            datas1[[ij]]$importModel(prev2,  datas[[ij2]]$family,numvar=if(kk<0) NULL else kk)
            })
                ,num_cols ="value", addName = "dataset")
      })
      ,num_cols="value",addName="numvars")
  }
  rmsval$value[rmsval$measure=="AUC_all"] =  rmsval$value[rmsval$measure=="AUC_all"]-0.5
  return(cbind(rmsval, label))
}

.rewind<-function(models, datas){
  for(k in 1:length(models))models[[k]]$unwind(datas,k) 
}
.trim<-function(models, datas,to_keep , inds_m=1:length(models)){
  for(k in inds_m){
    models[[k]]$keep(to_keep)
    lapply(datas, function(d) d$train[[k]]$keep(to_keep))
  }
}
.updateModels<-function(models, datas, prev,inds_m, useInternalCVAsStopping){
  for(k in inds_m){
    model = models[[k]]
    if(!useInternalCVAsStopping) model$prev_old=model$prev
    model$prev=prev
    model$cnt = length(prev)
    lapply(datas, function(d) d$updateModel1(k,prev, useOld = useInternalCVAsStopping ))
  }
}

.reorder<-function(models, datas, o, inds_m = 1:length(models)){
  for(k in inds_m){
    model = models[[k]]
    model$reorder(o)
    lapply(datas, function(d) d$reorder(o,k))
  }
}
.summarise<-function(rms_cv, dim1="beam",avg=F){
  #dim1 = 'pheno'
  select1 =   apply(rms_cv[,which(names(rms_cv) %in% dim1),drop=F],1,paste, collapse=".")
  select=factor(select1, lev = select1[!duplicated(select1)])
  
  levs = levels(select)
  .merge1(lapply(levs, function(l1){
    rms_cv1 =rms_cv[select==l1,,drop=F]
  rms_cv2 = rms_cv1[1,,drop=F]
  rms_cv2$value = if(avg) mean(rms_cv1$value,na.rm=T) else sum(rms_cv1$value,na.rm=T)
  rms_cv2$subpheno="comb"
  rms_cv2
  }), num_cols="value")
}
.runAnalysisAll<-function(datas, datas1, params,rdsdir=NULL){
  options(params)
  if(!is.null(rdsdir)) dir.create(rdsdir, rec=T)
  {  #this is setup
      genes_incls=getOption("genes_incls",NULL)
      cv_skips_allowed=getOption("cv_skip",0) ;  signatures = getOption("signatures",NULL)
      incls = getOption("incls",list(names(datas$all$data)))
      ##this is for cross validation evaluation 
      if(length(datas1)>0){
        mi1 = match(names(datas), names(datas1))
        if(length(which(!is.na(mi1)))==0) stop("!!")
        datas = datas[!is.na(mi1)]
        datas1 = datas1[mi1[!is.na(mi1)]]
        # print(cbind(names(datas), names(datas1)))
      }
      num_var =sum(unlist( lapply(datas[[1]]$data, ncol)))
      max_vars = min(num_var,getOption("maxvars",100))
      data1_inds = 1:length(datas1)
      names(data1_inds) = names(datas1)
      types_all = names(datas[[1]]$data)
      names(types_all) = types_all
      var_threshs = lapply(types_all, function(x) getOption("var_thresh",1e-8))
      print("initialising")
      for(ij in 1:length(datas)){print(ij);datas[[ij]]$init1(var_threshs, genes_incls=genes_incls)}
      nreps_all =unlist(lapply(datas, function(d) ncol(d$looc$incl)))
      nreps = table(nreps_all)
      if(length(nreps)>1) {
        ##NEED TO REDO BATCHING
        if(getOption("fspls.batch",0)==0) stop("!!")
        print(paste("new nrep",min(nreps_all)-1 ))
        for(ij in 1:length(datas))datas[[ij]]$init1(var_threshs, nrep = min(nreps_all)-1, batch=0)
        nreps_all =unlist(lapply(datas, function(d) ncol(d$looc$incl)))
        nreps = table(nreps_all)
        if(length(nreps)>1) stop("could not resolve problems with nrep and batch")
      }
      nrep1  = as.numeric(names(nreps))
      nrep=if(nrep1==1)1 else nrep1-1 
      print("building models")
      models=lapply(1:nrep1, function(k) {
        ##following line could be applied on yPreds above
        rmsv_=.merge1(lapply(datas, function(d)d$getRMSV(k)),num_cols="value",addName="data")  #this is call 3 to datas
        rmsv_prev1 =   .getRMSPrev(rmsv_)
        modelObj$new(names(datas), names(datas[[1]]$data),  rmsv_prev1)
        
      })
      names(models) = 1:length(models)
      # jks=1; k=1;
      rms_prev=sum(models[[1]]$rmsv_prev, na.rm=T)
      skipped=0
      useInternalCVAsStopping= getOption("useInternalCVAsStopping",F)
      orderCV = getOption("orderCV",T)
      pv_only=getOption("fspls.pval_only",F)
      beam=getOption("fspls.beam",c(1,1))
      logpthresh=log(getOption("fspls.pthresh1",0.05))
      replaceModel = beam[1] *beam[2] >1  || !useInternalCVAsStopping  ## can only not replace models if no beam
  }
  for(jks in 1:length(incls)){
    print(paste("training", jks))
    for(model in models) model$finished=F
    finished = FALSE
    to_keep = 1  ## which samples to use in next round
    maxvars = getOption("maxvars",100)
    cnt1=1
    maxvars_jks = maxvars[min(jks, length(maxvars))]
    while(!finished && models[[length(models)]]$cnt<sum(maxvars)){
      if(cnt1>maxvars_jks){
        break;
      }else{
        cnt1 = cnt1+1
      }
      incl=incls[[jks]]
      inds_m = if(useInternalCVAsStopping) 1:length(models) else length(models)
      pvslist_all =lapply(inds_m, function(k){
        model = models[[k]]
        nxt_var = if(models[[k]]$cnt<length(signatures) ) signatures[[models[[k]]$cnt+1]] else NULL
        model$simpleForwardTrain(datas,k,  exclude=exclude,nxt_var = nxt_var, incl = incl ,  weights=NULL, to_keep=to_keep)
      })
      pvslist_all1 =pvslist_all[[length(pvslist_all)]]$pv
      pvslist_all2 = pvslist_all[[length(pvslist_all)]]$cumpv
      print(pvslist_all1)
      pvtk1=  pvslist_all1<logpthresh
      if(length(which(pvtk1))==0){
        .trim(models, datas,c(),inds_m  )
        print("breaking on pvalue")
        finished=T
        break;
      }
      
      if(!pv_only &useInternalCVAsStopping & nrep>1){
        rms_cv1 = .merge1(lapply(datas,function(d) .summarise(d$getRMSVAll())),num_cols="value",addName="dataset")
        rms_prev1=min(rms_cv1$value)
        if(rms_prev1>=rms_prev){
          .trim(models, datas,c(),inds_m )
          finished=T
          break;
        }
      }
      vars = names(pvslist_all2)
      prev = models[[length(models)]]$prev
      if(length(models)>1 && replaceModel) {
        .updateModels(models, datas, prev,1:(length(models)-1),useInternalCVAsStopping)
      }
      if(!pv_only){
          if(length(models)==1){
              rms_cv = .merge1(lapply(datas,function(d) .summarise(d$getRMSV(length(models)))),num_cols="value",addName="dataset")
        }else{
              rms_cv = .merge1(lapply(datas,function(d) .summarise(d$getRMSVAll())),num_cols="value",addName="dataset")
          }
          rms_cv = rms_cv[match(vars, rms_cv$beam),,drop=F]
          o = if(orderCV) order(rms_cv$value) else order(pvslist_all2) 
          names(o) = if(orderCV)  rms_cv$beam[o]  else names(pvslist_all2)[o]
      }else{
        o = order(pvslist_all2)
        names(o) = names(pvslist_all2)[o]
      }
         to_keep =o[1:min(length(vars),beam[2])]
         #attr(models,"to_keep")=to_keep
          pvtk=   pvslist_all1[to_keep]<logpthresh 
          if(length(which(pvtk))==0){
            .trim(models, datas,c()  )
            finished=T
            break;
          }
          #if(length(which(!pvtk))>0) pvtk = pvtk[1:which(!pvtk)[1]]  # only keep until first non-sig
          #to_keep = to_keep[pvtk]
#          tokeep1 = rms_cv$beam[to_keep]
          if(!pv_only){
            if(!useInternalCVAsStopping || length(models)==1) rms_prev1=min(rms_cv$value[to_keep],na.rm=T)
            if(rms_prev1>rms_prev){
              skipped = skipped+1
              print(paste("skipped ",skipped))
            }else{
              #reset skipped and reset bar
              print("reset skipped")
              skipped=0
              rms_prev = rms_prev1
            }
          
          if(skipped> cv_skips_allowed){
            if(length(to_keep)>0){ ## need to wind back since rms got worse
              for(k in 1:length(models)){
                models[[k]]$keep(c())
                lapply(datas, function(d) d$train[[k]]$keep(c()))
              }
            }
            finished=T
            print(paste("breaking here because cv detetiorating", models[[1]]$cnt))
            
          #  for(k in which(!finished)){
          #    models[[k]]$unwind(datas,k) 
          #  }
            break;
          }
          }
          if(length(to_keep)==0){
            finished=T
            break
          }
          to_keep =1:min(length(vars),beam[2])  ## really only comes into play if beam>1
          attr(models,"to_keep")=to_keep
          .reorder(models, datas, o)          
#          attr(models,"to_keep")=to_keep
      if(!is.null(rdsdir)){
        outdir1 =  paste(rdsdir,paste(models[[length(models)]]$cnt,"rds",sep="."),sep="/")
        xls_file = paste(rdsdir,paste(models[[length(models)]]$cnt,"xlsx",sep="."),sep="/")
        print(paste("saving data to outdir1", outdir1))
        #sr =  .saveDatas(datas,datas1, params, outdir1,xls_file, types=c("AUC","area","sens_spec"))
        sr =  .saveDatas(datas,datas1, outdir1,xls_file, types=c("AUC","area","youden_sens","youden_spec"))
        if(is.null(sr$rms$validation)){
          if(is.null(sr$rms$crossvalidation)){
            print(head(.summarise(sr$rms$discovery, c("beam","pheno","measure"),avg=T)))
          }else{
            print(head(.summarise(sr$rms$crossvalidation, c("beam","pheno","measure"),avg=T)))
            
          }
        }else{
            print(head(.summarise(sr$rms$validation, c("beam","pheno","measure"),avg=T)))
        }
        
      }
    }
  }
  
  
  invisible(models)
}


.logLik<-function(y,ypred=mean(y), family="binomial"){
  if(family=="binomial"){
    yp1=.logistic(ypred)
    ll1=sum(apply(cbind(y,yp1), 1,function(v) if(v[1]==1) log(v[2]) else log(1-v[2])))
  }else{
    ll1=sum(apply(cbind(y,ypred), 1,function(v) dnorm(v[2], mean = v[1],log=T)))
  }
  
}
.lrt<-function(ll1,ll2,  df1 =  attr(ll1,"df")[1] ,  df2 = attr(ll2,"df")[1],totweight=1,log.p=F){
  if(df2>df1) stop(paste("df2 should be bigger than df1 ",df2,df1))
  pchisq((2*(ll1 - ll2))/(totweight),df1-df2,lower.tail=FALSE,log.p=log.p)
}
.norm1<-function(g){
  sqrt(sum((g - mean(g, na.rm = T))^2, na.rm = T))
}
#.plotAUC<-function(rdsf){
#  
#}


.plotRMS2<-function(tEnvs,type="cv",fam=tEnvs[[1]][[1]]$datas[[1]]$family,
                    measures = NULL){
  #print(tEnv$rms_list_all)
  a7=.merge1(lapply(tEnvs, function(tEnvs1){
    .merge1(lapply(tEnvs1, function(tEnv){
      names(tEnv$rms_list_all)=1:length(tEnv$rms_list_all)
      a1 = .merge1(list(
        cv=.merge1(lapply(tEnv$rms_list_all, function(t) .mod(t$cross)), addName="lens1",num_cols=c("value","lens")),
        val=.merge1(lapply(tEnv$rms_list_all, function(t) .mod(t$val)), addName="lens1",num_cols=c("value","lens")),
        disc=.merge1(lapply(tEnv$rms_list_all, function(t) .mod(t$disc)), addName="lens1",num_cols=c("value","lens"))),
        num_cols = c("value","lens","lens1"),addName="type")
      o = order(a1$dataset)
      #a1[o,]
      a1$value = -1*(a1$value)
      a2 = subset(a1,beam1==1)
      if(fam=="multinomial"){
        return(a2)
      }
      a3 =  a2[grep("\\.mid$",a2$subpheno),,drop=F]
      a4 =  a2[grep("\\.low$",a2$subpheno),,drop=F]
      a5 =  a2[grep("\\.high$",a2$subpheno),,drop=F]
      low = a4$value
      high = a5$value
      a6 = cbind(a3,low, high)
      a6
    }),addName="dataset1",num_cols=c("value","lens","lens1","low","high"))
  }),addName="drug",num_cols=c("value","lens","low","lens1","high"))
  a71 = a7[a7$type==type,,drop=F]
  if(!is.null(measures)){
    a71 = subset(a71, measure %in% measures)
  }
  if(fam %in% c("multinomial","binomial")){
    ggp<-ggplot(a71,aes(x=lens1, y=value, color=subpheno, shape=dataset,linetype=dataset))+geom_line()+geom_point()
    ggp<-ggp+facet_grid(measure ~dataset)
  }else{
    
    ggp<-ggplot(a71,aes(x=lens1, y=value, color=dataset1, linetype=beam1))+geom_line()+geom_point()
    ggp<-ggp+geom_ribbon(aes(ymin=low, ymax=high, fill=dataset),alpha = 0.05)
    ggp<-ggp+facet_grid(measure~drug,scales="free_y")+ggtitle(type) # +ggtitle(paste(a1$beam[a1$len==max(a1$lens)][1],"cross validation",sep="\n"))
  }
  return(ggp)
}



.replot2<-function(tEnvs, type1="cv", type="area", outpdf1=NULL, typed="datas"){
  df = .merge1(lapply(tEnvs, function(tEnvs1){
    .merge1(lapply(tEnvs1, function(tEnv){
      names(tEnv$rms_list_all)=1:length(tEnv$rms_list_all)
      
      .merge1(lapply(tEnv[[typed]], function(d){
        y1=d$y[,1]
        yp =d$ypreds_all$ypreds[[1]][[1]][,1]
        roc1 =roc(y1,yp,plot=F)
        nme = names(ar$datas[[1]]$train[[1]]$prev)[1]
        print(paste(nme, roc1$auc))
        nme=gsub(",","\n",nme)
        title=paste(nme, collapse="\n")
        title=nme
        area_plot=getAreaPlot(yp, y1,nme)#, title=title)
        area_plot = addExtraLines(area_plot,yp,y1,nme)
        mat = d$train[[length(d$train)]]$prev[[1]]$betas[[1]]
        if(is.matrix(mat) ){
          beta= try(mat[nrow(mat),1])
          beta = sign(beta)*log10(abs(beta))
          area_plot1 = data.frame(rbind(c(0,beta,0,"beta","summary",nme), c(1,beta,0,"beta","summary",nme)))
          names(area_plot1) = names(area_plot)
          area_plot=rbind(area_plot, area_plot1)
        }
        area_plot
        #  if(type=="area") return(plotAreas(yp,y1,title=nme)) else return(ggroc(roc1))
      }), num_cols = c("knots","value"), addName="lineage")
    }), num_cols = c("knots","value","counts"), addName="index")
  }),num_cols = c("knots","value","counts"), addName="drug")
  df$label = apply(cbind(as.character(df$pheno), round(df$value,2)),1,paste,collapse=" ")
  beta_inds = df$pheno=='beta'
  beta_scale = max(abs(df$value[beta_inds]))
  df$value[beta_inds] = df$value[beta_inds]/(2*beta_scale) + .5
  txt_df = subset(df ,type=="summary" & knots==0)
  txt_df$knots=0.02
  txt_df$knots[txt_df$pheno=="beta"] = 0.95
  
  ggp1=ggplot(df, aes(x=knots, y=value, color=pheno, linetype=type))+geom_line()+geom_point(aes(size=counts))+ggtitle(rdsfile)
  if(!is.infinite(beta_scale)){
    ggp1=ggp1+scale_y_continuous(
      "value", 
      sec.axis = sec_axis(~ (.+.5)*(2*beta_scale), name = "betas")
    )  
  }
  ggp1= ggp1+geom_text(data=txt_df,nudge_y=0.02,
                       inherit.aes =T, 
                       aes(x=knots,y=value,label=label,color=pheno),size=2)
  ggp1<-ggp1+facet_grid(lineage~drug)
  attr(ggp1,"text_df")=txt_df
  if(!is.null(outpdf1))  try(ggsave(outpdf1, plot=ggp1, width =45, height =45, units = "cm",limitsize=F))
  
  invisible(ggp1)
}



.convert<-function(data,mode="rna", expand=F, max=100, factor=F){
  names(data) =gsub("x_data","data",tolower(names(data)))
  if(is.list(data$y)){
    data$y = as.matrix(data$y)
    if(nrow(data$y)!=nrow(data$data)) stop("!!")
    dimnames(data$y)[[1]] = dimnames(data$data)[[1]]
  }
  if(!is.matrix(data$y)) {
    data$y =  data.frame(matrix(data$y, dimnames = list(dimnames(data$data)[[1]],"y")))
    if(is.character(data$y[1,1])){
      if(expand){
        data$y = .makeMultiClass(data$y)
        
      }else{
        data$y[[1]] = factor(data$y[[1]])
      }
    }
  }
  
  if(is.null(dimnames(data$y))) dimnames(data$y) = list(dimnames(data$data)[[1]], "y")
  if(factor){
    data$y=data.frame(lapply(data.frame(data$y), function(z) factor(paste("X",z,sep=""))))
  }
  if(max < nrow(data$y)){
    inds = sample.int(nrow(data$y),max)
    data$y = data$y[inds,,drop=F]
    data$data = data$data[inds,,drop=F]
  }
  
  dataset = list(data$data);
  names(dataset)=mode
  y = data$y
  list(dataset=dataset,y=y)
}

######from rawlinson paper
## get data from this repo https://github.com/dn-ra/FSPLS-publication-repo via git clone
.readRawlinsonData<-function(filenames,path){
  if(is.null(names(filenames))) names(filenames)=lapply(filenames, function(f)rev(strsplit(f,"/")[[1]])[1])
  
  #files = grep(".Rds",dir(path,f=T,rec=T),v=T)
  datas=lapply(filenames, function(file){
    data_in = readRDS(paste(path,file,sep="/"))
    .convert(data=data_in, mode="rna", expand=F)
  })
  datasets = lapply(datas, function(d) list(dataset=d$dataset, y = d$y))
  datasets
}

.getFinalModel<-function(all_models, pheno_groups=names(all_models[[1]]),target_size=NULL){
  names(pheno_groups) = pheno_groups
  final_models = lapply(pheno_groups, function(pg){
  full_models = lapply(all_models, function(all_models1) all_models1[[pg]][['full']])
  full_models = full_models[unlist(lapply(full_models, length))>0]
  model_size= unlist(lapply(names(full_models), function(x) length(strsplit(x,";")[[1]])))
  if(is.null(target_size) || is.character(target_size)){
    target_size = max(model_size)
  }
  final_model = full_models[[which(model_size==target_size)]]
  final_model
  })
}

#mod1 = final_models[[1]][[1]]; final_model = final_models[[1]]
.extractWeights<-function(final_models){
  weights = lapply(final_models, function(final_model){
model_weights = unlist(lapply(final_model, function(mod1){
  phen1= names(mod1$betas)
  names(phen1) = phen1
  lapply(phen1, function(p1){
    bet=mod1$betas[[p1]]
    varn = data.frame(t(data.frame(lapply(names(mod1$var_names), function(str)strsplit(str,"\\.")[1]))))
    names(varn)=c("type","var")
    df=cbind(varn,bet)
    df
  })
}),rec=F)
cbind(model_weights[[1]][,1:2], data.frame(lapply(model_weights, function(mw)mw[,-(1:2)])))

  })
}
