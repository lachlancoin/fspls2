
rbind_sparse <- function(mats) {
  # sanity check: all matrices must have the same number of columns
  ncols <- vapply(mats, ncol, integer(1))
  stopifnot(length(unique(ncols)) == 1)
  ncol_total <- ncols[1]
  
  # row offset for each matrix = cumulative sum of preceding row counts
  nrows <- vapply(mats, nrow, integer(1))
  offsets <- c(0, cumsum(nrows)[-length(nrows)])
  
  # extract triplets from each matrix and shift row indices
  triplets <- Map(function(m, off) {
    t <- Matrix::mat2triplet(m)
    list(i = t$i + off, j = t$j, x = t$x)
  }, mats, offsets)
  i_all <- unlist(lapply(triplets, `[[`, "i"), use.names = FALSE)
  j_all <- unlist(lapply(triplets, `[[`, "j"), use.names = FALSE)
  x_all <- unlist(lapply(triplets, `[[`, "x"), use.names = FALSE)
  
  sparseMatrix(i = i_all, j = j_all, x = x_all,
               dims = c(sum(nrows), ncol_total))
}


##drop, which are the info cols, not numerical
##assumes first col is rowname
.processBigCSV<-function(filename, rownames = 1, drop=2:6, chunk_size=100000, output = paste0(filename,".mm"), row_mean_thresh=0){
  
        dt <- data.table::fread(filename, drop = drop,
                                skip = 0, nrows = 0, header = TRUE)
        #chunk_size <- 50000
        coln =names(dt)
      nrow_total=NA
      mats <- list()
      mats_type = list()
      i <- 1
      skip <- 1  # header
    max = 3e9
      repeat {
      
        dt <- data.table::fread(filename, drop = drop,
                    skip = skip, nrows = chunk_size, header = FALSE)
      
        if (nrow(dt) == 0) break
        m1 = as(as.matrix(dt[,-1,drop=FALSE]), "sparseMatrix")
        rs = rowSums(m1)/ncol(m1)
    
        subind = rs>=row_mean_thresh
    
        m1 = m1[subind,,drop=FALSE]
        rownames(m1) = dt[[1]][subind]
        mats[[i]] <- m1
        skip <- skip + chunk_size
        i <- i + 1
        if(i>max) break
        if (nrow(dt) < chunk_size) break
      }
     mat1 = rbind_sparse(mats)
    colnames(mat1) = coln[-1]
    Matrix::writeMM(mat1, output)
    
}
#writeMM(m1, "counts.mm")


#' Fit a GLM, suppressing perfect separation warnings
#'
#' @param ... Arguments passed to glm
#' @return A fitted glm object
#' @noRd
glm1 <- function(...) {
  args <- list(...)
  withCallingHandlers(
    do.call(glm, args),
    warning = function(w) {
      if (grepl("fitted probabilities numerically 0 or 1 occurred", conditionMessage(w))) {
        invokeRestart("muffleWarning")
      }
    }
  )
}


##BIG MATRIX NOT SUPPORTED IN THIS VERSION
is.big.matrix<-function(mat){
  if(typeof(mat)=="S4") return(FALSE);
  #print(typeof(mat))
  return (FALSE)
 # return (bigmemory::is.big.matrix(mat));
}

###used strictly for debugging
#.self<-function(dh){
#  assign("self", dh, envir = .GlobalEnv) 
#  assign("private",  dh[[".__enclos_env__"]]$private, envir = .GlobalEnv) 
#  assign("super", dh[[".__enclos_env__"]]$super, envir = .GlobalEnv)
#}

## this from sce single cell format
.convertSCEToSparse<-function(sce){
comb = cbind(sce$donorID, sce$majorCluster)
dimnames(comb)[[2]] = c("donorID","majorCluster")
ncohort = length(levels(sce$donorID))
ncluster = length(levels(sce$majorCluster))

j = (comb[,1]-1)*ncluster + comb[,2]
i = 1:nrow(comb)
dims = c(nrow(comb), ncohort *ncluster)
spM = sparseMatrix(i,j,dims = dims) #dimnames=list(sce$cohortID, sce$majorCluster))
cell_names = levels(sce$majorCluster)
sample_names = levels(sce$donorID)
coln=unlist(lapply(sample_names, function(sn){
  paste(sn, cell_names,sep=".")
}))
colnames(spM) = coln

countsMatrix = sce@assays@data$X %*% spM
countsMatrix
}

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
    }),recursive=TRUE)
    #eval$model1 = eval$model
    eval$model= (unlist(lapply(models, function(m){
      paste(names(models2)[match(m, models2)],collapse=";")
    })))
    attr(eval,"translate")=models2
    eval
}
.print_verbose<-function(txt,str,lev){
  if(lev<=getOption("fspls.verbose",0)) {
    cat(paste(txt,"::"))
    cat(str)
  }
}

#' Plot the results from evaluateAll
#'
#' @param flags a list of options
#' @return corrected flags
#' @export
check_flags<-function(flags){
  if(!is.null(.readFlag(flags,"nrep",NULL))){
    flags[['nfold']] = flags[['nrep']]
    warning("replaced nrep with nfold")
  }
  if(!is.null(.readFlag(flags,"batch",NULL))){
    flags[['batchsize']] = flags[['batch']]
    
   warning("replace batch with batchsize")
  }
  invisible(flags)
}
.avg<-function(eval0){
  nme_cols1 = c("data","subpheno","measure","pheno","trainedOn","pheno_group","numvars","cv","family")
  rem_cols = names(eval0)[!(names(eval0) %in% nme_cols1)]
  nme_cols2 = c(nme_cols1, rem_cols)
  mi = match(nme_cols2, names(eval0))
eval1 = eval0[,mi] 
##not working

eval2 = tidyr::unite(head(eval1),"data:subpheno", remove=FALSE)
length(unique(eval1$`data:family`))
}

.calcEval1<-function(eval0, rename=TRUE, len = 3){ #c("trainedOn","measure","subpheno")
  eval0_1 = unite(eval0, "cohort_measure_pheno_trained", "data","subpheno","measure","pheno","trainedOn",remove=FALSE)
  eval0_avg = subset(eval0_1, model=="avg") |>  tibble::add_column("fullmodel"="avg")
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
      do.call(rbind, replicate(length(inds1[inds2]), eval1[i,], simplify=FALSE))  |>  tibble::add_column(fullmodel = eval1$model[inds1[inds2]])
    }))
    }))
  }))
  eval0_avg$cv="CV=avg"
  eval2$cv = paste("CV=",eval2$cv)
  eval21 = rbind(eval2,eval0_avg)
  eval21$isfull = paste("FULL=",eval21$isfull);
  eval3 = unite(eval21,"cv_full", "cv","isfull",remove=FALSE,sep=" ")
  eval3$cv_full[eval3$cv=="CV=avg"]=eval3$cv[eval3$cv=="CV=avg"]
  eval3
}

.modify<-function(eval3, shape_color,
                  shape_color_nme ){
  if(!(shape_color_nme %in% names(eval3))){
    eval3_sub=eval3[,names(eval3) %in% shape_color,drop=FALSE]
    for(jk in 1:length(eval3_sub)) eval3_sub[[jk]]=factor(eval3_sub[[jk]])
    levs1 = lapply(eval3_sub, function(vv) levels(vv))
    levs_all = levs1[[1]]
    if(length(levs1)>1){
      for(jk in 2:length(levs1)){
      levs_all = unlist(lapply(levs1[[jk]], function(levs11) paste(levs_all, levs11)))
      }
    }
      eval4 = eval3 |>  tibble::add_column(shape_color_nme = apply(eval3_sub,1,paste, collapse=" "))
      
      eval4[['shape_color_nme']]=factor(eval4[['shape_color_nme']], levels = levs_all)
        names(eval4) = gsub("shape_color_nme", shape_color_nme, names(eval4))
        return(eval4)
  }
  eval3
}


plotEval2<-function(eval1,...){
  if(is.null(eval1[['experiment_id']])){
    return(plotEval(eval1,...))
  }
  expt1 = unique(eval1$experiment_id)
  names(expt1) = expt1
  lapply(expt1,function(expt) {
    eval3 = subset(eval1, experiment_id==expt)
    plotEval(eval3,...)
  })
}

#'  Plot the results from evaluateAll
#' @param eval3 a tibble or data frame from evaluateAll method from dataH
#' @param shape_color a list of column names to encode as both shape and color
#' @param shape defaults to shape_color but can be specified separately
#' @param text a list of column names to include as text labels
#' @param dotsize a numerical column to encode dotsize, or a number for constant sizes
#' @param color defaults to shape_color but can be specified separately 
#' @param linetype column controlling the type of line
#' @param showranges whether to show 95\% CI as ribbon
#' @param txtsize size of text
#' @param logy whether to log y axis
#' @param legend whether to show a legend
#' @param sep_by list of column names to separate into diff plots
#' @param scales used in faceting can be free or fixed
#' @param point whether to show dots
#' @param line whether to include line
#' @param labelsize the size of the text labels
#' @param grid0 how to facet on y 
#' @param grid1 how to facet on x
#' @param title the title
#' @param title1 columns to include in the title
#' @returns ggplot2 plot
#' @export
plotEval<-function(eval3,
          shape_color=c("pheno","subpheno"),
          shape=shape_color,
          text ="variable",
          dotsize="nsamps",
           color=shape_color,
          linetype=shape,
          showranges=TRUE, ## linetype="fullmodel"
           txtsize=1,logy=FALSE,legend=FALSE,sep_by="",scales="free",point=TRUE,line=TRUE,
       labelsize=2,
           grid0 = c("cohort","measure"),grid1 = "cv_full",title="", title1=""
          ){
  
  grid0= grid0[grid0 %in% names(eval3)];  grid1= grid1[grid1 %in% names(eval3)]
  
  eval3$beam = factor(eval3$beam, levels =  sort(unique(as.numeric(eval3$beam))))
  
  eval3$sign = as.character(eval3$sign)
  eval3$sign = factor(eval3$sign, levels = c(-1,0,1), labels=c("-","","+"))
      eval3 = eval3[,names(eval3) %in% c(shape_color,sep_by, linetype, text, dotsize, color, title1, grid0, grid1, "numvars","mid","low","high"),drop=FALSE]
  
  
  l1 = apply(eval3,1,paste, collapse="::");
  #print(which(duplicated(l1)))
  eval3 = eval3[!duplicated(l1),,drop=FALSE]
  showtext = length(text)>0
  linetype_nme = paste(linetype, collapse="_");
  grid0_nme = paste(grid0, collapse="_")
  grid1_nme = if(length(grid1)==0) NULL else  paste(grid1,collapse="_");
  shape_nme = paste(shape,collapse="_")
  color_nme = paste(color,collapse="_")
  text_nme = paste(text,collapse="_")
  
  sep_by_nme = "sep_by"
  subphens = table(eval3$subpheno)
  subphens = subphens[order(as.numeric(unlist(lapply(names(subphens), function(str)strsplit(str,"\\|")[[1]][1]))))]
  eval3$subpheno = factor(eval3$subpheno, levels = names(subphens))
  eval3$measure = factor(eval3$measure)
  eval3 = .modify(eval3, color, color_nme)
  eval3 = .modify(eval3, shape, shape_nme)
  eval3 = .modify(eval3, text, text_nme)
  
  eval3 = .modify(eval3, linetype, linetype_nme)
  eval3 = .modify(eval3, sep_by, sep_by_nme)
  eval3 = .modify(eval3,grid0, grid0_nme)
  if(!is.null(grid1_nme)){
  eval3 = .modify(eval3,grid1, grid1_nme)
  }
  eval2 = eval3
  
  phenos = unique(eval2$sep_by)
  names(phenos)=phenos
  
  
  
 
  eval2$numvars = as.numeric(eval2$numvars)
# eval2$isfull = (eval2$isfull+1)/2.0
  ggps=lapply(phenos, function(ph){ 
    eval5 = subset(eval2, sep_by==ph & !is.na(mid))
    ph3 = paste(sort(unique(apply(eval5[,names(eval5) %in% title1, drop=FALSE],1,paste,collapse=","))), collapse=" ")
  ggp<-ggplot(eval5);
  if(point) {
    if(dotsize %in% names(eval5)){
    ggp<-ggp+geom_point(aes_string(x="numvars", y="mid",  shape=shape_nme,size=dotsize, color=color_nme))
    }else{
      ggp<-ggp+geom_point(aes_string(x="numvars", y="mid",  shape=shape_nme, color=color_nme),size=dotsize)
      
    }
  }
  if(line)  ggp<-ggp+geom_line(aes_string(x="numvars", y="mid", linetype=linetype_nme,  color=color_nme)) #+ggtitle(paste(title,ph, ph3))
  if(showtext) ggp<-ggp+geom_text_repel(aes_string(x="numvars", y="mid", label=text_nme,color=color_nme),size=labelsize,max.overlaps = 1000) #+ggtitle(paste(title,ph, ph3))
  
  if(showranges){ ## geom_ribbon vs geom_errorbar
    
   ggp<-ggp+ geom_ribbon(aes_string(x = "numvars", ymin="low", ymax="high",linetype=linetype_nme,color=color_nme, fill = color_nme ), alpha = 0.1)
  }
  if(nchar(grid0_nme)>0){
    if(!is.null(grid1_nme) && nchar(grid1_nme)>0){
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

rotate <- function(x, n) {
  n <- n %% length(x)          # handle n > length or negative n
  x[c((n + 1):length(x), seq_len(n))]
}

permute<-function(x, seed, norm = 1, offset = 0){
  x1=(x+offset)/norm; 
  rotate(x1, seed)
}

invpermute = function(y1, seed, norm=1,offset=0){
 rotate(y1,-seed)*norm-offset 
}

randomize <- function(y, seed, norm = 1, offset=0) {
  y1=(y+offset)/norm; 
  set.seed(seed)
  inds <- sample.int(length(y1))
  y1[inds]
}

invrandomize <- function(y1, seed, norm=1, offset=0) {
  set.seed(seed)
  inds <- sample.int(length(y1))
  y1[order(inds)]*norm-offset   # clean and idiomatic
  # equivalently: y1[match(1:length(y1), inds)]  -- your version, also correct
  # or:           { r <- y1; r[inds] <- y1; r }  -- most explicit
}

# Test
#y <- 1:10
#y_rand <- randomize(y, seed = 42)
#y_back <- invrandomize(y_rand, seed = 42)
#identical(y, y_back)  # TRUE


.calcAUCW1<-function(ypred,y, w,
                    conf.level=getOption("conf.level",0.95)
){
                     tw = table(w)
                       nw = as.numeric(names(tw))
                       tw = tw * nw
                       tw = tw/sum(tw)
                       
                       aucs=data.frame(lapply(nw, function(w1) .calcAUCW(y,ypred, w==w1,conf.level=conf.level)))
                       apply(aucs,1,function(auc) sum(auc*tw))   
    
                     }

.calcAUCW<-function(ypred,y, w,
                    conf.level=getOption("conf.level",0.95)
){
  
    if(length(which(y==1))==0 || length(which(y==0))==0 ) return(c(0,0.5,1))
    ##NOTE THE WEIGHTED VERSION SEEMS NOT TO WORK WELL
    #  print(y)
    #  print(ypred)
    roc1=try(roc(y,ypred, quiet=TRUE))
    if(inherits(roc1,"try-error")){
      return(rep(NA,3))
    }
    #print(ci(roc1)[1:3])
    cir = ci(roc1, conf.level=conf.level)[1:3]
    if(is.na(cir[2])) cir[2] = roc1$auc
    if(is.na(cir[1])) cir[1] = 0
    if(is.na(cir[3])) cir[3] = 1
    return(cir)
 
}



.calcAverageAccuracy<-function(comb){
  comb1=unite(comb,comb,experiment_id, cv_full, measure,sep="__");
  comb1_lev = unique(comb1$comb); names(comb1_lev) = comb1_lev
  medians=.merge1_new(lapply(comb1_lev, function(x){
    aa= subset(comb1, comb==x)
    data.frame(list(mid=median(aa$mid,na.rm=TRUE), comb=x))
  }))
  medians = medians |>  separate(comb,c("experiment_id","cv_full","measure"), sep="__")
  medians=medians[order(medians$cv_full),]
  
  subset(medians, cv_full=="CV=avg")
}

.getPermFuncs<-function(n,norm=1, offset=0, CHECK=FALSE){ ## although these are same, every invocation will give different results
  if(n==0) return(list())
  inds = sample.int(4*n,2*n, replace=FALSE)
 # inds = inds[which(inds %%size !=0)][1:n]
  names(inds) = inds
  transf=list(invfunc=paste0("function(y,seed) invpermute(y,seed, ",norm,", ",offset,")"),
              func=paste0("function(x,seed) permute(x,seed, ",norm,", ",offset,")"),params=as.list(inds))
  transf
}


.getRandomFuncs<-function(n,norm=1, offset=0, CHECK=FALSE){ ## although these are same, every invocation will give different results
  if(n==0) return(list())
  n1 = max(n, 1000)
 inds = sample.int(2*n,n, replace=FALSE)
 names(inds) = inds
 which(duplicated(inds))
 transf=list(invfunc=paste0("function(y,seed) invrandomize(y,seed, ",norm,", ",offset,")"),
             func=paste0("function(x,seed) randomize(x,seed, ",norm,", ",offset,")"),params=as.list(inds))
 if(CHECK){
   ggp= .checkInverse1(transf)
   ggp
 }
 transf
}
.checkInverse1<-function(transf){
  xx=seq(-5,5,by=0.05)
      t_y1= lapply(transf[1:2],function(funcstr1) eval(str2lang(funcstr1)) )
      l1 = lapply(transf[[3]], function(pow){
      ab = data.frame(.checkInverse(t_y1,pow, xx))
  ab
  })
      names(l1)=transf[[3]]
  ab_all=.merge1_new(l1,addName="func")
  ggplot(ab_all, aes(x=xx, y=y_1, color=func))+geom_line();# +scale_y_log10()
}
.checkInverse<-function(t_y1,pow = t_y1[[3]],  xx = -10:10 ){
#  func0 = lapply(transform_y, function(t_y) eval(str2lang(t_y[[1]])))  ## should be inverse
  #func1 = lapply(transform_y, function(t_y) eval(str2lang(t_y[[2]])))  ## should be inverse
  llm=lapply(pow,function(pow1){
         y_1 = t_y1[[1]](xx,pow1)
  y_2 = t_y1[[2]](y_1,pow1)
  m1 = cbind(xx, y_1, y_2)
  #m1 = cbind(t_y[[2]](t_y[[1]](xx)), xx)
  
  diffs = apply(m1[,c(1,3)],1,diff)
  if(max(abs(diffs), na.rm=TRUE)>1e-7){
  #  print(t_y1)
  #  print(m1)
 stop("not inverse")
  }
  invisible(m1)
  })
  invisible(llm)
 # print("ok")
}

#' Get object for transforming the y variables
#'
#' @param pows power to raise to
#' @param offset offset to subtract 
#' @param n_random how many random variables
#' @param perm  whether random is permutation
#' @param norm rescaling factor
#' @param CHECK check whether inverse function works
#' @return transformation object
#' @export
getYTransform<-function(pows = c(1),offset=0,  n_random=1,perm=FALSE, norm=1,CHECK=FALSE){
  if(n_random <1) warning(" recommended to have at least one random permutation");
  if(!( 1 %in% pows)) warning("recommended to have a pow of 1, which is the untransformed y ")
  funcs = list()
  if(length(pows)>0) funcs = c(funcs,list(pow=.getTransformFuncs(pows,  norm = norm, offset=offset, CHECK=CHECK) ))
  if(n_random>0){
    if(perm) {
      funcs = c(funcs,list( rand=.getPermFuncs(n_random,  norm = norm, offset=offset)))
      
    }else{
      funcs = c(funcs,list( rand=.getRandomFuncs(n_random, norm = norm, offset = offset)))
      
    }
  }
 
  funcs
}
getXTransform<-function(pows= c(1),offset=1e-10){
  .getTransformFuncs(pows, offset=offset)
}

##exp is problematic because of neg numbers, particularly after centralisation
##could work with adding back in mean values?? may not generalise to unseen datasets
getExpFunc<-function(pows, rev=FALSE,offset=0.1, CHECK=FALSE){  
  if(length(which(pows<=0))>0) stop("not possible")
  if(length(pows)==0) return (list())
  names(pows)=pows
 
 
                 
  if(rev){
    warning("probably not going to work because x gets centralised before transform")
    warning("this assumes that x+offset is strictly positive")
    
    transf=list(invfunc =  paste0("function(y,pow,norm=",norm,",offset=",offset,") expfunc(y,pow,norm,offset)"),
                func=paste0("function(x,pow,norm=",norm,",offset=",offset,") logfunc(x,pow,norm,offset)"), params = as.list(pows))
    
  }else{
    warning("this assumes that y+offset strictly positive")
    transf=list(
                invfunc=paste0("function(y,pow,norm=",norm,",offset=",offset,") logfunc(y,pow,norm,offset)"),
                func =  paste0("function(x,pow,norm=",norm,",offset=",offset,") expfunc(x,pow,norm,offset)"),
                params = as.list(pows))
   
  }
  if(CHECK){
    ggp= .checkInverse1(transf)
    ggp
  }
  transf
}
#.getRandomFuncs(3)
#c(list(x=c("function(y) y","function(y) y")), random_funcs)
expfunc<-function(y,pow, norm=1,offset=0.1){
  y1 = pow^y
 y1*norm-offset
}
#pow is positive
logfunc<-function(x,pow, norm=1,offset=0.1){
  x1=(x+offset)/norm; 
  log(x1)/log(pow)
}

powfunc<-function(x,pow, norm=1, offset=0.0001){
 # pow=v[1]; norm=v[2]; offset=v[3]
  x1=(x+offset)/norm; 
  sign(x1) * abs(x1)^pow
}
invpowfunc<-function(y,pow, norm=1, offset=0.0001){
  #pow=v[1]; norm=v[2]; offset=v[3];
  y1 = sign(y) * abs(y)^(1/pow); 
  y1*norm-offset
}

##this gets transformations for x variable
##invfunc applied to y; func applied to x
.getTransformFuncs<-function(pows,
                             offset=0.1, norm=1,CHECK=FALSE){
  names(pows)=pows
  transf=list(invfunc=paste0("function(y,seed) invpowfunc(y,seed, ",norm,", ",offset,")"),
              func=paste0("function(x,seed) powfunc(x,seed, ",norm,", ",offset,")"),params=as.list(pows))
  
   
    if(CHECK){
    ggp= .checkInverse1(transf)
    ggp
    }
 transf
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
  mat1=replicate(n,vec,simplify =FALSE)
  if(by=="row"){
     mat =do.call(rbind,mat1 )
  }else{
    mat = do.call(cbind,mat1)
  }
 mergeSparseMatrices(counts, Matrix(mat,sparse=TRUE),by=by)
}


#' Convert matrix into sparse submatrices
#'
#' @param counts a matrix
#' @param inds indices for subsetting
#' @param by can be col or row
#' @return sparse sub matrix
#' @noRd
getSparseSubMatrix<-function(counts, inds,by='col'){
  #newNames = if(by=='col')colnames((counts)) else rownames(counts)
  if(!is.numeric(inds)) {
    tomatch = if(by=="col")colnames(counts) else rownames(counts);
    inds1 = match(inds, tomatch)
    nonna = !is.na(inds1)
    
    inds = inds1[nonna]
    
  }
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
  }),recursive=FALSE))
  
  dimnames(df)[[1]] = dimnames(y)[[1]]
  df
}

.sumChisq<-function(pvi){
  if(length(pvi)==0) return(0)
  pvi1=unlist(pvi)
  pchisq(sum(qchisq(pvi1,df=1, lower.tail=FALSE, log.p=TRUE)), df=length(pvi1), lower.tail=FALSE, log.p=TRUE)
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
  
    t1[,mi1,drop=FALSE]
  })
}

.split<-function(vn1,nme){
  mod_nme=unique(vn1[[nme]]);names(mod_nme)=mod_nme
  lapply(mod_nme, function(mn) vn1[vn1[[nme]]==mn,,drop=FALSE])
}
.splitAll<-function(vn1, nmes, func=identity){
  aa1=.split(vn1, nmes[1])
  if(length(nmes)==1){
    return(lapply(aa1,func))
  }else{
    lapply(aa1, .splitAll,nmes[-1], func)
  } 
}
.merge_all<-function(ri, nmes, func=identity){
  if(length(nmes)==1){
    .merge1_new(lapply(ri,func ), addName=nmes[[1]])
  }else{
    .merge1_new(lapply(ri, .merge_all, nmes[-1], func), addName = nmes[1])
  }
}

.lapply_nme=function(y,FUN){
  
  res = lapply(names(y),FUN)
  names(res) = names(y)
  res
}

.merge_lapply=function(phens,nme, FUN){
   .merge1_new(lapply(phens,FUN), addName=nme)
  
}

.merge_lapply_nme=function(phens,nme, FUN){
  .merge1_new(.lapply_nme(phens,FUN), addName=nme)
  
}

.merge1_new<-function(t,num_cols = c(), addName=NULL, checkNames=TRUE){
  if(checkNames && length(t)>0){
    nme_aa = names(t[[1]])
    t1 = lapply(t, function(aa1){
      aa1[,match(nme_aa,names(aa1)),drop=FALSE]
    })
    t = t1
  }
  t = t[!unlist(lapply(t, is.null))]
  
  if(length(t)==0) return(NULL)
  if(!is.data.frame(t[[1]])) stop("not dataframe")
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
  rownames = 1:sum(nrows)
  for(k in 1:length(t)){
    range = start[[k]] : end[[k]]
    df[range,1:ncol] = t[[k]]
    rn = rownames(t[[k]])
    if(!is.null(rn)){
      rownames[range] = rn
    }
    if(add){
      df[range,ncol+1] = names(t)[[k]]
    }
  }
  if(length(which(duplicated(rownames)))==0){
    #print(rownames)
    rownames(df) = rownames
  }
  df
}



.merge1<-function(t,num_cols = c(),uniq_cols=c(), addName=NULL, rowNames=FALSE){
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
    facts = factor(apply(t2[,which(names(t2) %in% uniq_cols),drop=FALSE],1,paste, collapse="."))
    t2 = t2[!duplicated(facts),,drop=FALSE]
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



whichpart1<-function(angle, n=10, one_for_each=FALSE, return_scores=FALSE){
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



.keepUniq<-function(vars, removeTrans=TRUE){
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


.findMinRMSV<-function(rmsv_, mult=rep(1, length(levels(rmsv_$data))), fspls.sum=TRUE){
  
  if(fspls.sum){
    beams = levels(rmsv_$beam)
    names(beams)=beams
    rmsv_allp = lapply(beams, function(p1){
      rmsv_sub=subset(rmsv_, beam==p1)
      rmsv_sub1  =pivot_wider(rmsv_sub,id_cols="phens", values_from="value", names_from="data")
      rmsv=as.matrix(rmsv_sub1[,-1,drop=FALSE])
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
      rmsv=as.matrix(rmsv_sub1[,-1,drop=FALSE])
      dimnames(rmsv) [[1]] = rmsv_sub1$beam
      rmsv1 =   (rmsv %*% mult)[,1]
      attr(rmsv1,"phen")=p1
      rmsv1
    })
    rmsv_allp[[which.min(unlist(lapply(rmsv_allp, function(v) min(v))))[1]]]
  }
}
.replot_dist<-function(rdsfile,  min=1,max=4, outpdf1=NULL, typed="datas"){
  fi = grep("rds",dir(rdsfile, full.names=TRUE),value=TRUE)
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
      roc1 =roc(y1,yp,plot=FALSE)
      df2=data.frame(cbind(roc1$sensitivities, roc1$specificities))
      names(df2) = c("sens","spec")
      nme = names(ar$datas[[1]]$train[[1]]$prev)[1]
      #print(paste(nme, roc1$auc))
      nme=gsub(",","\n",nme)
      title=paste(nme, collapse="\n")
      title=nme
     
      cbind(df2,title)
      #  if(type=="area") return(plotAreas(yp,y1,title=nme)) else return(ggroc(roc1))
    }), num_cols = c("sens","spec"), addName="lineage")
  }), num_cols = c("sens","spec"), addName="index")
  
  
  df = .merge1(lapply(fi, function(fi1){
    #print(fi1)
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
      df4=df4[df4$counts>0,,drop=FALSE]
      #ggplot(df4, aes(x=knots, y=value))+geom_point()
      #ggplot(df4, aes(x=knots, y=name, size=value))+geom_point()
      #ggplot(df2, aes(x=factor(pheno), y=value))+geom_jitter(alpha = 0.9, width=0.1, size=.2)
      roc1 =roc(y1,yp,plot=FALSE)
      nme = names(ar$datas[[1]]$train[[1]]$prev)[1]
    #  print(paste(nme, roc1$auc))
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
#                       inherit.aes =TRUE, 
#                       aes(x=knots,y=value,label=label,color=pheno),size=2)
  ggp1<-ggp1+facet_grid(lineage~title)+ggtitle(rdsfile)+scale_x_continuous(limits = c(-0.1,1.1))+scale_y_continuous(limits = c(-0.1,1.1))
  ggp1
if(!is.null(outpdf1))  try(ggsave(outpdf1, plot=ggp1, width =45, height =45, units = "cm",limitsize=FALSE))
  #attr(ggp1,"text_df")=txt_df
  invisible(ggp1)
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
                       inherit.aes =TRUE, 
                       aes(x=knots,y=value,label=label,color=pheno),size=2)
  ggp1<-ggp1+facet_grid(lineage~drug)
  attr(ggp1,"text_df")=txt_df
  if(!is.null(outpdf1))  try(ggsave(outpdf1, plot=ggp1, width =45, height =45, units = "cm",limitsize=FALSE))
  
  ggp2=ggplot(txt_df, aes(x=drug, y=value, color=lineage))+geom_point()+facet_grid(label1~.)
  
  list(ggp1 = ggp1,ggp2=ggp2)
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
  df_l=lapply(levs, function(t){
    inds = y1==t
    tab1 = table(yp[inds])
    cdf=ecdf(yp[inds])
    mi2=match(knots, names(tab1))
   # no_dupl = !duplicated(mi2)
    #mi2 = mi2[!duplicated(mi2)]
    tab1_col = tab1[mi2]
    tab1_col[is.na(tab1_col)]=0
   # df1 = data.frame(cbind(knots[no_dupl],cdf(knots[no_dupl]), tab1_col))
    
    df1 = data.frame(cbind(knots,cdf(knots), tab1_col))
    
    names(df1)= c("knots","value","counts")
    df1
  })
  df=.merge1_new(df_l, addName="subpheno",num_cols=c("knots","value"))
  df = cbind(df, "cumulative")
  names(df)[ncol(df)] = "type"
  df
}
summariseAreaPlot<-function(df){
  pw = pivot_wider(df, names_from="subpheno", id_cols="knots")
  diff=0
  diff_v = apply(pw, 1, function(v)abs(v[2]-v[1]))
  for(k in 2:nrow(pw)){
    delta = (pw[k,1]-pw[k-1,1])
    diff = diff+((pw[k-1,2] + pw[k,2] - (pw[k-1,3] + pw[k,3]))/2) * delta
  }
  area=diff
  

 
  max_diff_x = pw$knots[which(diff_v==max(diff_v))]
  
 list(area=c(NA,diff[[1]],NA),
  max_diff=c(NA,max(diff_v),NA),
  max_diff_x= c(min(max_diff_x),median(max_diff_x), max(max_diff_x)))
}
.combinePv<-function(pvs){
  return(min(pvs))
  chisq = qchisq(pvs, df=1, lower.tail=FALSE)
  pchisq(sum(chisq), df = length(pvs), lower.tail=FALSE)
}

.plotArea1<-function(predictions,family="binomial",rename=TRUE, 
                     subset = 1:length(predictions[[1]][[1]]), #[[1]]),
                     len = 1,grid="model~pheno", max_vars = 100){
  area_p= .getAreaPlot1(predictions,family)
  .plotArea(area_p, rename=rename, len = len,grid=grid,  max_vars = max_vars,title=title)
}

.getAreaPlot1<-function(predictions0, families="binomial"){
  if(is.null(families)) families = names(predictions0[[1]][[1]][[1]])
  names(families)=families
 #j=1; model = predictions0[[2]]; train = model[[1]]; test = train[[1]]; famnme = families[[1]]; fam =test[[famnme]]; family = strsplit(famnme,"\\.")[[1]][1];  phens = dimnames(fam$y)[[2]]; names(phens)=phens
  
  area_p=.merge1_new(lapply(predictions0, function(train){  
    .merge1_new(lapply(train, function(model){
      test=model
     # .merge1_new(lapply(model, function(test){
        families1 = names(test); names(families1)=families1
        .merge1_new(lapply(families1, function(famnme){
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
               ap = ap |>  tibble::add_column(sample=dimnames(fam$y)[[1]])
              
            }else if(family=="ordinal"){
              ap = getAreaPlot(fam$ypred, fam$y[,j])
        }else{
            ap = getAreaPlot(fam$ypred[,j], fam$y[,j])
          
            }
            ap
          #}),addName="subpheno")
        }), addName="pheno")
        }), addName="family")
    #  }), addName="test")
    }), addName="model")
    }), addName="train")
  #attr(area_p,"family")=family
  lens = unlist(vapply(area_p$model, function(x) if(x=="empty") 0 else length(strsplit(x,";")[[1]]), FUN.VALUE = c(1)))
 
   area_p|>  tibble::add_column(lens = lens)
}
.plotAreaSep=function(area_p, sep="pheno",...){
  phens=unique(area_p[[sep]]);names(phens)=phens
 lapply(phens, function(p){
   p = phens[[1]]
  area_p1=area_p[area_p[[sep]]==p,,drop=FALSE]
  .plotArea(area_p1, ...)
 })  
}

.takeAvg<-function(eval1){
  #does not deal with subpheno if different
  nme_cols = c("cv_full","measure","numvars","data","cv","isfull","experiment_id","fullmodel","model","transf","subpheno")
  ev2=pivot_wider(eval1,
                  id_cols =nme_cols ,
                  names_from="pheno", values_from="mid")
  meanv=apply(ev2[,-(1:length(nme_cols))],1,mean)
  data.frame(ev2[,1:length(nme_cols)]|>  tibble::add_column(mid=meanv, pheno="avg", nsamps=10))
}
.takeMax1<-function(a1, max_vars=100){
  ab=unite(subset(a1, lens<=max_vars),"comb","sample","CV","pheno",remove=FALSE)
  levs = unique(ab$comb); names(levs)=levs
  ab2=.merge1_new(lapply(levs, function(l1){
    sub1 = subset(ab, comb==l1 & !is.na(knots))
    sub1[which.max(sub1$lens),,drop=FALSE]
  }))
  #ab2$pheno = factor(ab2$pheno, labels=levels(ab$pheno))
  ab2
}


.takeMax<-function(area_p1){
  ab=unite(area_p1,"comb","sample","pheno","CV",remove=FALSE)
  ab$pheno = factor(ab$pheno)
  levs = unique(ab$comb); names(levs)=levs
  ab2=.merge1_new(lapply(levs, function(l1){
    sub1 = subset(ab, comb==l1 & !is.na(knots))
    sub1[which.max(sub1$lens),,drop=FALSE]
  }))
  ab2$pheno = factor(ab2$pheno, labels=levels(ab$pheno))
  ab2
}
.getTextLayer<-function(area_p1,p = c(0.1,0.9),r2=FALSE, arrange_by_sample=FALSE, reorder=FALSE){
  ab = if(arrange_by_sample) unite(area_p1,"comb","sample","CV",remove=FALSE) else unite(area_p1,"comb","lens","pheno","CV",remove=FALSE)
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
                    pheno = pheno1, knots=quantile(sub1$knots,na.rm=TRUE, probs = p[1]), value=quantile(sub1$value, na.rm=TRUE, probs =p[2]), 
                    label= label))
  }))
    ab2$CV = factor(ab2$CV, levels = levels(ab$CV))
   ab2$pheno = factor(ab2$pheno, levels = levels(ab$pheno))
   ab2$sample = factor(ab2$sample, levels = levels(ab$sample))
  if(arrange_by_sample){
    phens2=levels(ab2$sample);names(phens2)=phens2
    
    levs1 = names(sort(unlist(lapply(phens2, function(p2){
      sb1=subset(ab2, sample==p2)
      max(sb1$cor,na.rm=TRUE)
    })),decreasing=TRUE))
    if(reorder)ab2$sample = factor(ab2$sample, levels=levs1)
  }else{
    phens2=levels(ab2$pheno);names(phens2)=phens2
    
  levs1 = names(sort(unlist(lapply(phens2, function(p2){
    sb1=subset(ab2, pheno==p2 & CV!="No CV")
    max(sb1$cor)
  })),decreasing=TRUE))
  if(reorder)ab2$pheno = factor(ab2$pheno, levels=levs1)
  }
  ab2
}

.compareCorrelation<-function(area_p1, maxv=5){
  samps = unique(area_p1$sample); names(samps)=samps
  cvs=unique(area_p1$CV); names(cvs)=cvs
  area_p_max0 = data.frame(.takeMax1(area_p1, max_vars=0))
  
  mvars = sort(unique(area_p1$lens)); names(mvars)=mvars
  mvars = mvars[mvars>0 & mvars<maxv]
  res_all=.merge1_new(lapply(mvars, function(max_vars){
  #  print(max_vars)
  area_p_max = data.frame(.takeMax1(area_p1, max_vars=max_vars))
  data.frame(lapply(cvs, function(cv){
    r2 = unlist(lapply(samps, function(samp){
    # print(paste(cv, samp))
      apm = subset(area_p_max, sample==samp & CV==cv)
      apm0 = subset(area_p_max0, sample==samp & CV==cv)
      if(nrow(apm0)<=1) return(NA)
      resd = try(cor(apm$knots, apm$value,use="pairwise.complete.obs")-cor(apm0$knots,apm0$value,use="pairwise.complete.obs"))
      if(inherits(resd,"try-error")) return(NA)
      resd
    }))
    mean(r2,na.rm=TRUE)
  }))
  }),addName="maxvars")
  res_all
}
.plotArea<-function(area_p0, CV=FALSE,alpha =0.8, maxphens = 50,p=c(0.1,0.9),
                    arrange_by_sample=FALSE,takeMax=FALSE,
                    rename=TRUE, len = 1,shapes=FALSE,code_len=3,showText=FALSE,
                    max_samps=100,reorder=FALSE,
                    grid=if(!CV) "model~pheno" else "lens ~pheno",  max_vars = 10,scales="free",title="", addText=FALSE,r2=FALSE){
  #print("now plotting")
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
   area_p1 = area_p1|> tibble::add_column(pheno_code=pheno_code)
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
        ggp<-ggplot(area_p1, aes(x=knots, y=value, color=pheno))
        
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
    
    #if(!is.null(textLayer))    textLayer =textLayer|> tibble::add_column(pheno_code= "avg" )
    #area_p1$pheno_code
      if(length(unique(area_p1$pheno))>30){
        ggp=ggp+geom_point( size=1,alpha = alpha)
      }else{
        ggp=ggp+geom_text_repel(aes(label=pheno_code), alpha = alpha) 
      }
    }else{
    ggp=ggp+geom_text(aes(label=sample), alpha = alpha) 
    }
  }else if(max(area_p1$counts,na.rm=TRUE)>5) {
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
    sum(ab$value,na.rm=TRUE)  #could also try min ? 
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
       tbls1 = names(sort(table(dc), decreasing=TRUE))
       factor(dc, levels = tbls1)
      #factor(bacviral, levels = tblkk)
    }))
    names(res) = names(json1[[nme]])
    res
  })
  res_all1 = data.frame(unlist(res_all,recursive=FALSE))
  cbind(phenot, res_all1)
}
.inferFamily<-function(y){
  unlist(lapply(y, function(y1) {
    ty = table(y1)
    if(is.factor(y1)) {
      return(if(length(ty)>2) "multinomial" else "binomial")
    }else if(length(ty)==2 && names(ty)[[1]]=="0" && names(ty)[[2]] == "1"){
     return("binomial")      
    }else if( max(apply(cbind(y1,round(y1) ),1,function(x) abs(diff(x))),na.rm=TRUE)==0){
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
  dat1 = dat$clone(deep=TRUE)
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
  
# print( table(cbind(disease_class_predicted),data$y[,1]))
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
  }),recursive=FALSE)
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
       # print(kk)
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
.summarise<-function(rms_cv, dim1="beam",avg=FALSE){
  #dim1 = 'pheno'
  select1 =   apply(rms_cv[,which(names(rms_cv) %in% dim1),drop=FALSE],1,paste, collapse=".")
  select=factor(select1, levels= select1[!duplicated(select1)])
  
  levs = levels(select)
  .merge1(lapply(levs, function(l1){
    rms_cv1 =rms_cv[select==l1,,drop=FALSE]
  rms_cv2 = rms_cv1[1,,drop=FALSE]
  rms_cv2$value = if(avg) mean(rms_cv1$value,na.rm=TRUE) else sum(rms_cv1$value,na.rm=TRUE)
  rms_cv2$subpheno="comb"
  rms_cv2
  }), num_cols="value")
}

.logLik<-function(y,ypred=mean(y), family="binomial"){
  if(family=="binomial"){
    yp1=.logistic(ypred)
    ll1=sum(apply(cbind(y,yp1), 1,function(v) if(v[1]==1) log(v[2]) else log(1-v[2])))
  }else{
    ll1=sum(apply(cbind(y,ypred), 1,function(v) dnorm(v[2], mean = v[1],log=TRUE)))
  }
  
}
.lrt<-function(ll1,ll2,  df1 =  attr(ll1,"df")[1] ,  df2 = attr(ll2,"df")[1],totweight=1,log.p=FALSE){
  if(df2>df1) stop(paste("df2 should be bigger than df1 ",df2,df1))
  pchisq((2*(ll1 - ll2))/(totweight),df1-df2,lower.tail=FALSE,log.p=log.p)
}
.norm1<-function(g){
  sqrt(sum((g - mean(g, na.rm =TRUE))^2, na.rm =TRUE))
}
#.plotAUC<-function(rdsf){
#  
#}







.convert<-function(data,mode="rna", nme="none",expand=FALSE, max=100, factor=FALSE){
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
    data$y = data$y[inds,,drop=FALSE]
    data$data = data$data[inds,,drop=FALSE]
  }
  
  dataset = list(data$data);
  names(dataset)=mode
  y = data$y
  list(dataset=dataset,y=y,nme=nme)
}

######from rawlinson paper
## get data from this repo https://github.com/dn-ra/FSPLS-publication-repo via git clone
.readRawlinsonData<-function(filenames,path){
  if(is.null(names(filenames))) names(filenames)=lapply(filenames, function(f)rev(strsplit(f,"/")[[1]])[1])
  
  #files = grep(".Rds",dir(path,f=TRUE,recursive=TRUERUE),value=TRUE)
  datas=lapply(filenames, function(file){
    data_in = readRDS(paste(path,file,sep="/"))
    .convert(data=data_in, mode="rna", expand=FALSE)
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
}),recursive=FALSE)
cbind(model_weights[[1]][,1:2], data.frame(lapply(model_weights, function(mw)mw[,-(1:2)])))

  })
}
