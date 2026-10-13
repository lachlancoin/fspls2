
#' extracts the full model variables from a cross validation run
#' @param leng length of dataset
#' @param nfold number of reps
#' @param batch batch size (one of batch or nfold is non NA)
#' @param randomize whether to randomize
#' @export
getFolds<-function(leng, nfold=NA, batch=NA, randomize=TRUE, balance = NULL){
  if(!is.null(balance)){
    if(is.na(nfold)) stop("need to specify fold")
    levs = levels(balance);names(levs) = levs
    folds_sep = lapply(levs, function(lev){
       inds1 = which(balance==lev)
       folds = getFolds(length(inds1), nfold = nfold, batch = batch)
       folds1 = lapply(folds, function(fold){
         inds1[fold]
       })
    })
    fols_ = 1:nfold; names(fols_) = fols_;
    return(lapply(fols_,function(fol){
      v = unlist(lapply(folds_sep, function(fs)fs[[fol]])); names(v) = NULL;
      v
    }))
  }
  len_y1 = leng
  if(is.na(nfold) && is.na(batch))stop("either batch or NA needs to be numeric")
  inds = if(randomize) sample.int(len_y1, len_y1) else 1:len_y1
  if(!is.na(nfold) && nfold>0){
      left_over=len_y1 %% nfold
      batch = (len_y1-left_over)/nfold
  }else{
      nfold  = floor(len_y1 / batch)
      left_over=len_y1 %% nfold
    
  }
      len_y2 = len_y1 -left_over
      reps_k = seq(1, len_y2+1, batch)
     # reps_k[length(reps_k)] = len_y1+1
    #  incl = matrix(TRUE, nrow = len_y, ncol =length(reps_k)-1 )
      if(left_over>0){
        for(jj in 1:left_over){
          for(kk in (jj+1):length(reps_k)){
            reps_k[kk]=reps_k[kk]+1
          }
        }
      }
      folds= vector('list', length(reps_k)-1) ## need angle object
      
      for(jj in 2:length(reps_k)){
        folds[[jj-1]] = sort(inds[reps_k[jj-1]:(reps_k[jj]-1)])
      }
    
     if(length(which(duplicated(unlist(folds))))>0) stop("error")
     if(length(unique(unlist(folds))) < len_y1) stop("error")
      names(folds) = 1:length(folds)
#      if(addFull) folds[["full"]] = 1:len_y1
  folds
  

}

loocObj<-R6::R6Class("loocObj", public = list(
  useAll="vector",
  nfold="numeric",
  batch="numeric",
  incl="matrix",
  seed="numeric",
  #reps="vector",
  rowsToDo="vector",
  nrows = "numeric",
  initialize=function(data,   ## data can be null if nrows is not NULL
                      incl_full =TRUE,
                      nrows = nrow(data$data[[1]]),
                      folds = c(),
                       ## randomisation
                      nfold=getOption("nfold",1), 
                      batch=getOption("batchsize",NA),
                      seed = 42,
                      pheno_balance=NULL
                      ){
    
    self$nfold=nfold
    self$batch=batch
    self$seed = seed
    self$nrows = nrows;
    rand = sample(nrows)
    len_y = nrows #length(data$y[,1])
    if(length(folds)==0){
      folds = getFolds(nrows, nfold=nfold, batch = batch, randomize=TRUE);
    }

    useAll = rep(TRUE, len_y)
    len_y1 = len_y #length(non_na_inds)
  
      self$nfold = length(folds);
      self$batch=0;
      incl = matrix(T, nrow = len_y, ncol =length(folds) )
      for(jk in 1:length(folds)){
        incl[folds[[jk]],jk]=FALSE
      }
      if(incl_full) incl = cbind( incl, useAll)
      self$incl = incl
   
  }
))
