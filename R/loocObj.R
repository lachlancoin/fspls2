

loocObj<-R6Class("loocObj", public = list(
  useAll="vector",
  nrep="numeric",
  incl="matrix",
  #reps="vector",
  rowsToDo="vector",
  initialize=function(data,   ## data can be null if nrows is not NULL
                      incl_full = T,
                      nrows = nrow(data$data[[1]]),
                      rand = sample(nrows),  ## randomisation
                      nrep=getOption("fspls.nrep",1), 
                      batch=getOption("fspls.batch",0),
                      pheno_balance=NULL
                      ){
   
    len_y = nrows #length(data$y[,1])
#    nonNA = apply(data$y,c(1,2),function(x) !is.na(x))
#    non_na_inds = which(apply(nonNA, 1,function(x) length(which(x))>0))
    useAll = rep(T, len_y)
    len_y1 = len_y #length(non_na_inds)
    tbl2 = NULL
    
    if(length(pheno_balance)==1){
      phens1 = unlist(lapply(data$y, function(y1) which(dimnames(y1)[[2]]==pheno_balance)))
      y2 = data$y[[names(phens1)[[1]]]][,phens1[[1]]]
      tbl2=table(y2)
    }
    
    if(!is.null(nrep) && nrep!=0){
#      batch  = ceiling((len_y1)/nrep)
     # incl = matrix(T, nrow = len_y1, ncol =nrep)
      left_over=len_y1 %% nrep
      batch = (len_y1-left_over)/nrep
      reps_k = seq(1+left_over, len_y1+1, batch)
      reps_k[length(reps_k)] = len_y1+1
      incl = matrix(T, nrow = len_y, ncol =length(reps_k)-1 )
      for(jj in 2:length(reps_k)){
        incl[reps_k[jj-1]:(reps_k[jj]-1),jj-1] = F
      }
      if(left_over>0){
        for(jk in 1:length(left_over)){
          incl[jk, jk]=FALSE
        }
      }
      if(incl_full) incl = cbind( incl, useAll)
      self$incl = incl[rand,]
      # max=nrep
    }else if(is.null(batch) || batch==0 || is.na(batch) || nrep==1){
      self$useAll = useAll
      self$nrep=0
      self$incl=as.matrix(useAll)
    }else if(length(pheno_balance)==1  && batch==1 && min(tbl2)>5){
      phens1 = unlist(lapply(data$y, function(y1) which(dimnames(y1)[[2]]==pheno_balance)))
      y2 = data$y[[names(phens1)[[1]]]][,phens1[[1]]]
      tbl=table(y2)
      cnt=min(tbl)
      nmet = sort(unique(y2)); names(nmet) = nmet
      v2 = rep(T, length(y2))
      incl2 = unlist(lapply(nmet, function(nmet1){
        s1= as.list(sample(which(y2==nmet1),cnt,replace=T))
      }))
        incl= data.frame(lapply(incl2, function(s2){
            v3 = v2
            v3[s2] = FALSE
            v3
        }))
        if(incl_full) incl = cbind( incl, useAll)
        self$incl = incl
    }else{
      reps_k = seq(1, len_y1+1, batch)
      reps_k[length(reps_k)] = len_y1+1
      incl = matrix(T, nrow = len_y, ncol =length(reps_k)-1 )
      for(jj in 2:length(reps_k)){
        incl[reps_k[jj-1]:(reps_k[jj]-1),jj-1] = F
      }
      if(incl_full) incl = cbind( incl, useAll)
      self$incl = incl[rand,]
    }
  # if(nrep==1 && batch==0){
  #    print("keeping useAll only")
  #    self$incl = self$incl[,2,drop=F]
  #  }
    #print(paste("loocObj", dim(self$incl)))
  }
))
