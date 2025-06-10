## redundant should remove

#default_types=fromJSON('{"gaussian": "rank_correlation","binomial" : "AUC"}')
#default_types=fromJSON('{"gaussian": "rank_correlation","binomial" : "AUC","multinomial" : "AUC_all","ordinal" :"AUC_all"}')


.mod<-function(mat){
  if(is.null(mat)) return(NULL)
  lens=unlist(lapply(as.character(mat$beam), function(str)length(strsplit(str,",")[[1]])))
  beam1 = as.numeric(mat$beam)
  cbind(mat,lens, beam1)
}
.readDir<-function(rdsfile, model_len = getOption("fspls.nrep",1)+1){
  fi = grep(paste(model_len, "rds",sep="."),dir(rdsfile, full=T),v=T)
  #if(length(fi)<min) return(NULL)
  names(fi) = unlist(lapply(fi, function(st)gsub(".rds","",rev(strsplit(st,"/")[[1]])[1])))
  fi = fi[order(as.numeric(names(fi)))]
  fi = fi[as.numeric(names(fi))>0]
}
.makeTenv<-function(rdsfile){
 fi = .readDir(rdsfile)
  #if(length(fi)<min) return(NULL)
  #fi = fi[1:min(length(fi),max)]
  tEnv = list(rms_list_all = lapply(fi, function(fi1){
    ar = readRDS(fi1)
    ar$rms
  }))
}
.makeTenv1<-function(rdsfile){
  fi = .readDir(rdsfile)
  #if(length(fi)<min) return(NULL)
  #fi = fi[1:min(length(fi),max)]
     readRDS(fi)
}


.plotRMS<-function(tEnv){
  #print(tEnv$rms_list_all)
  a1 = .merge1(list(
    cv=.merge1(lapply(tEnv$rms_list_all, function(t) .mod(t$cross)), num_cols=c("value","lens")),
    val=.merge1(lapply(tEnv$rms_list_all, function(t) .mod(t$val)), num_cols=c("value","lens")),
    disc=.merge1(lapply(tEnv$rms_list_all, function(t) .mod(t$disc)), num_cols=c("value","lens"))),
    num_cols = c("value","lens"),addName="type")
  o = order(a1$dataset)
  #a1[o,]
  a1$value = -1*(a1$value)
  a2 = subset(a1,beam1==1)
  if(tEnv$datas[[1]]$family=="multinomial"){
    ggp<-ggplot(a2,aes(x=lens, y=value, color=subpheno, shape=dataset,linetype=dataset))+geom_line()+geom_point()
    ggp<-ggp+facet_grid(measure ~type)
  }else{
    a3 =  a2[grep("\\.mid$",a2$subpheno),,drop=F]
    a4 =  a2[grep("\\.low$",a2$subpheno),,drop=F]
    a5 =  a2[grep("\\.high$",a2$subpheno),,drop=F]
    low = a4$value
    high = a5$value
    a6 = cbind(a3,low, high)
    ggp<-ggplot(a6,aes(x=lens, y=value, color=dataset, linetype=beam1))+geom_line()+geom_point()
    ggp<-ggp+geom_ribbon(aes(ymin=low, ymax=high, fill=dataset),alpha = 0.05)
    ggp<-ggp+facet_grid(measure~type,scales="free_y")+ggtitle(paste(a1$beam[a1$len==max(a1$lens)][1],"cross validation",sep="\n"))
  }
  return(ggp)
}



trainingEnv<-R6Class("trainingEnv", public = list(
  models = "list",
  datas = "list",
  datas1 = "list",
  incls = "list",
  incl="vector",
  print_types="vector",
  default_types="vector",
  rms_list_all="list",
  rms_prev="numeric",
  to_keep="vector",
  jks="numeric",
  rdsdir="character",
  m1="list",
  cnt1="numeric",
  finished="vector",
  rmsv_cv0="matrix",
  num_var="numeric",
  maxvars="vector",
  maxvars_jks="vector",
  skipped="numeric",
  #UDVP="UDVPObj",  
  init1=function(){
    self$cnt1=0
    datas = self$datas
    datas1 = self$datas1
    data1_inds = 1:length(datas1)
    rms_list_all = list()
    names(data1_inds) = names(datas1)
    types_all = names(datas[[1]]$data)
    names(types_all) = types_all
    if(TRUE){
      var_threshs = lapply(types_all, function(x) getOption("var_thresh",1e-8))
      genes_incls=getOption("genes_incls",NULL)
      varn = getOption("varn",c())
      cv_skips_allowed=getOption("cv_skip",0) ;  signatures = getOption("signatures",NULL)
    #print("initialising")
      invisible(lapply(datas, function(d1) d1$init1(var_threshs, genes_incls=genes_incls)))
      invisible(lapply(datas, function(d1) d1$initTrain(varn=varn)))
      invisible(lapply(datas, function(d1) d1$initY()))
      invisible(lapply(datas1, function(d1) d1$initY()))
    }
    #var_threshs1 = lapply(var_threshs, function(xx)1e-9)
    #invisible(lapply(datas1, function(d1) d1$init1(var_threshs1, genes_incls=genes_incls, varn = varn)))
    
    nreps_all =unlist(lapply(datas, function(d) ncol(d$looc$incl)))
    nreps = table(nreps_all)
    if(length(nreps)>1) {
      ##NEED TO REDO BATCHING
      if(getOption("fspls.batch",0)==0) stop("!!")
      #print(paste("new nrep",min(nreps_all)-1 ))
      for(ij in 1:length(datas))datas[[ij]]$init1(var_threshs, nrep = min(nreps_all)-1, batch=0)
      nreps_all =unlist(lapply(datas, function(d) ncol(d$looc$incl)))
      nreps = table(nreps_all)
      if(length(nreps)>1) stop("could not resolve problems with nrep and batch")
    }
    nrep1  = as.numeric(names(nreps))
    nrep=if(nrep1==1)1 else nrep1-1 
    #print("building models")
    var=datas[[1]]$extractVar(varn)
    if(length(getOption("types_", fromJSON('{"gaussian": "rank_correlation","binomial" : "AUC"}'))[[datas[[1]]$family]])!=1) stop(getOption("types_", default_types))
    
    
   
    
    models=lapply(1:nrep1, function(k) {
      ##following line could be applied on yPreds above
      rmsv_=.merge1(lapply(datas, function(d)d$getRMSV(k)),num_cols="value",addName="data")  #this is call 3 to datas
      rmsv_prev1 =   .getRMSPrev(rmsv_)
      
      modelObj$new(names(datas), names(datas[[1]]$data),  rmsv_prev1,var = var)
      
    })
    names(models) = 1:length(models)
    
   
    
    # jks=1; k=1;
    self$models = models   
    self$rms_prev=sum(models[[1]]$rmsv_prev, na.rm=T)
    self$incls = getOption("incls",list(names(datas[[1]]$data)))
    self$num_var =sum(unlist( lapply(datas[[1]]$data, ncol)))
    self$maxvars = getOption("maxvars",self$num_var)
    if(length(models)>1){
      self$rmsv_cv0= .summarise(.calcRMSVAll(datas, NULL,cv=T, label="cv"),c("fspls.beam","pheno"))
    }else{
      self$rmsv_cv0 = NULL
    }
    self$skipped=0
    self$m1=lapply(1:length(models), function(ki){
      list(rms_prev=sum(models[[ki]]$rmsv_prev, na.rm=T), to_keep=1,skipped=0)
    })
    ##NEED TO CHECK THAT THIS OK
    self$finished=rep(F, length(self$m1))
  },
  saveDatas=function(){
    xls_file=NULL
    rdsdir = self$rdsdir
    models = self$models
    k = length(models)
    outdir1 =  paste(rdsdir,paste(1,k,"rds",sep="."),sep="/")
    xls_file = paste(rdsdir,paste(1,k,"xlsx",sep="."),sep="/")
      rms = self$rms_list_all
      data_new = lapply(self$datas, function(d) dataObj1$new(d))
     data1 =   if(is.null(self$datas1))NULL else lapply(self$datas1, function(d) dataObj1$new(d))
     rms = self$rms_list_all[[length(self$rms_list_all)]]
    resu = list(datas = data_new, datas1 = data1, rms =rms, rms_list_all = self$rms_list_all, file = file,
                incls_all =names(self$datas[[1]]$data), finished = self$finished,  nme = names(self$datas[[1]]$train[[1]]$prev))
    file = outdir1
    if(!is.null(file)) try(saveRDS(resu,file ))
    if(!is.null(xls_file)){
      betas = .merge1(lapply(self$datas, function(d) {
        .merge1(lapply(d$train[[length(d$train)]]$prev, function(p){
          .merge1(list("const"=.merge1(lapply(p$constants_proj, function(pb){
            df=data.frame(pb)
            names(df)="estimate"
            df
          }),addName="drug",num_cols="estimate"),
          "beta"=
          .merge1(lapply(p$betas, function(pb){
            df=data.frame(pb)
            names(df)="estimate"
            df
          }), addName="drug", num_cols ="estimate" )
          ),addName="type",num_cols="estimate")
        }), addName="beam",num_cols="estimate")
      }), addName="dataset", num_cols="estimate")
      comb1=list(auc=.merge1(rms[unlist(lapply(rms, length))>0], num_cols = "value"), betas=betas)
      #names(comb1)[1] = names(datas[[1]]$y)[[1]]
      try(WriteXLS(comb1, ExcelFileName =xls_file))
    }
    invisible(resu)
  },
  updateFromRDS=function(fi1){
    models = self$models
    if(!is.null(fi1)){
      ar = readRDS(fi1)
      mi1 = match(names(self$datas), names(ar$datas))
      nonNA = which(!is.na(mi1))
      if(length(nonNA)>0){
        for(ij in nonNA){
          for(k in 1:length(models)){
            prevk = ar$datas[[mi1[ij]]]$train[[k]]$prev
            self$datas[[ij]]$train[[k]]$prev = prevk
            self$models[[k]]$prev=lapply(prevk, function(pki){
              stateObj1$new(prev_i =NULL , b_i =NULL,var =pki$var )
            })
            self$models[[k]]$cnt = length(prevk[[1]]$var)
          }
          incls_all = ar$incls_all
          if(is.null(incls_all)) incls_all= names(datas[[1]]$data)
          vars_to_date = incls_all[sort(unique(unlist(lapply(self$models, function(m) lapply(m$prev, function(p)lapply(p$var, function(yy)yy[1]))))))]
          self$cnt1 = self$models[[1]]$cnt
          return(vars_to_date)
          #self$datas[[ij]]$updateFromRDS()
        }
        
      }
      #datas_prev = ar$datas []
      
    }
    return(NULL)
  },
  initialize=function(datas, datas1,incls){
    if(length(grep("\\.",names(datas)))>0) stop(" should not have . in name")
    if(length(grep("\\.",names(datas1)))>0) stop(" should not have . in name")
   rdsdir = getOption("rdsdir",NULL)
    self$datas = datas
    self$rms_list_all = list()
    self$datas1 =datas1
    types_ = getOption("types_", default_types)
    self$print_types = getOption("print_types",types_[[datas[[1]]$family[[1]]]])
    self$default_types = types_
    if(!is.null(rdsdir)){
      dir.create(rdsdir, rec=T)
      self$rdsdir=rdsdir
    }else{
      self$rdsdir=NULL
    }
    self$incls = incls
   
    ##this is for cross validation evaluation 
    if(length(datas1)>0){
      mi1 = match(names(datas), names(datas1))
      if(length(which(is.na(mi1)))==0){
        self$datas = datas[!is.na(mi1)]
        self$datas1 = datas1[mi1[!is.na(mi1)]]
      }
      # #print(cbind(names(datas), names(datas1)))
    }
    
   
    #    for(ij in 1:length(datas)){#print(ij);datas[[ij]]$init1(var_threshs, genes_incls=genes_incls)}
   
  },
    saveDatas1=function(){
      datas=self$datas; datas1=self$datas1; types=self$print_types
      types1_ = getOption("types_", default_types)
      if(!is.null(types)){
      
        types_ = types1_
        types_[[datas[[1]]$family[[1]]]] = types
        options("types_"=types_)
      }
      rms = list(
        validation = if(is.null(datas1)) NULL else .calcRMSVAll(datas, datas1,cv=F, label="validation",self$models),
        crossvalidation = if(length(datas[[1]]$train)<=1) NULL else .calcRMSVAll(datas, NULL,cv=T, label="cv",self$models),
        discovery = .calcRMSVAll(datas, datas,cv=F, label="discovery", self$models)
      )
     
        options("types_"=types1_)
      
      rms
  },
 
  runAnalysisInner=function(){
    models = self$models
    datas=self$datas
    datas1 = self$datas1
    incl = self$incl
    nrep1 = length(models)
    nrep=if(nrep1==1)1 else nrep1-1 
    cv_skips_allowed=getOption("cv_skip",0) ;  signatures = getOption("signatures",NULL)
    finished=F
    useInternalCVAsStopping= getOption("useInternalCVAsStopping",F)
    orderCV = getOption("orderCV",T)
    pv_only=getOption("fspls.pval_only",T)
    beam=getOption("fspls.beam",c(1,1))
    logpthresh=log(getOption("fspls.pthresh1",0.01))
    replaceModel = beam[1] *beam[2] >1  || !useInternalCVAsStopping  
   # cnt1 =0
    to_keep=1
#    print(paste("maxvars",self$maxvars_jks))
    while(!finished && length(to_keep)>0)
    {
      if(self$cnt1>=self$maxvars_jks){
        break;
      }else{
        self$cnt1 = self$cnt1+1
      }
      
      inds_m = if(useInternalCVAsStopping) 1:length(models) else length(models)
      pvslist_all =lapply(inds_m, function(k){
        model = models[[k]]
        nxt_var = if(models[[k]]$cnt<length(signatures) ) signatures[[models[[k]]$cnt+1]] else NULL
        model$simpleForwardTrain(datas,k,  exclude=exclude,nxt_var = nxt_var, incl = self$incl ,  weights=NULL, to_keep=to_keep)
      })
      pvslist_all1 =lapply(pvslist_all[[length(pvslist_all)]]$pvs_all1, function(x) min(unlist(x)))
      pvslist_all2 = pvslist_all1
    #  pvslist_all2 = pvslist_all[[length(pvslist_all)]]$cumpv  need to find another way
      #print(pvslist_all1)
      pvtk1=  pvslist_all1<logpthresh
      if(length(which(pvtk1))==0){
        .trim(models, datas,c(),inds_m=inds_m  )
       # print("breaking on pvalue")
        finished=T
        return(NULL)
      }
      
      if(!pv_only &useInternalCVAsStopping & nrep>1){
        rms_cv1 = .merge1(lapply(datas,function(d) .summarise(d$getRMSVAll())),num_cols="value",addName="dataset")
        rms_prev1=min(rms_cv1$value)
        if(rms_prev1>=self$rms_prev){
          .trim(models, datas,c(),inds_m )
          finished=T
          return(NULL)
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
        return(NULL)
      }
      #if(length(which(!pvtk))>0) pvtk = pvtk[1:which(!pvtk)[1]]  # only keep until first non-sig
      #to_keep = to_keep[pvtk]
      #          tokeep1 = rms_cv$beam[to_keep]
      if(!pv_only){
        if(!useInternalCVAsStopping || length(models)==1) rms_prev1=min(rms_cv$value[to_keep],na.rm=T)
        if(rms_prev1>self$rms_prev){
          skipped = skipped+1
       #   print(paste("skipped ",skipped))
        }else{
          #reset skipped and reset bar
         # print("reset skipped")
          skipped=0
          self$rms_prev = rms_prev1
        }
        
        if(skipped> cv_skips_allowed){
          if(length(to_keep)>0){ ## need to wind back since rms got worse
            for(k in 1:length(models)){
              models[[k]]$keep(c())
              lapply(datas, function(d) d$train[[k]]$keep(c()))
            }
          }
          finished=T
       #   print(paste("breaking here because cv detetiorating", models[[1]]$cnt))
          
          #  for(k in which(!finished)){
          #    models[[k]]$unwind(datas,k) 
          #  }
          return(NULL)
        }
      }
      if(length(to_keep)==0){
        finished=T
        return(NULL)
      }
      to_keep =1:min(length(vars),beam[2])  ## really only comes into play if beam>1
      attr(models,"to_keep")=to_keep
      self$rms_list_all[[length( self$rms_list_all)+1]] = self$saveDatas1()
      if(!is.null(self$rdsdir)) try(self$saveDatas())
      
      #if(!is.null(self$rdsdir))try(self$saveDatas())
   # print(unlist(models[[1]]$pvs_all))
      .reorder(models, datas, o)    
    }
  #  names(self$rms_list_all) = names(self$models[[1]]$pvs_all)[1:length(self$rms_list_all)]
  },
  updatePheno=function(phenos, drug, phenos1 = NULL,   drug1=NULL){
    invisible(lapply(self$datas, function(d) d$updateY(phenos[,match(drug, names(phenos)),drop=F])))
    if(!is.null(self$datas1)){
      #names(drug1) = drug
      invisible(lapply(self$datas1, function(d) d$updateY(phenos1[,match(drug1, dimnames(phenos1)),drop=F],nme=drug)))
    }
   # self$init1();
  },
  updateJk=function(jks){
    self$jks=jks
    self$incl = self$incls[[jks]]
    self$maxvars_jks = self$maxvars[min(jks, length(self$maxvars))]
    for(model in self$models) model$finished=F
    
  },
  runAnalysisAll=function(outpdf3=NULL, fi1 = NULL){
    print("initialising ...")
    self$init1()  ## initialises everything
    start_pos=1
    if(length(fi1)>0) {
      vars_to_date = self$updateFromRDS(fi1)
      cnts = unlist(lapply(self$incls, function(si) length(which(si %in% vars_to_date))),rec=F)
      start_pos = max(which(cnts>0))
    }
    print(".. done")
#    models = .initialise(datas,datas1)
   # models = self$models
    run_sep = getOption("fspls.run_sep",F) && length(self$models)>1
    cnt1 = self$cnt1
    useInternalCVAsStopping= getOption("useInternalCVAsStopping",F)
    
   
    ## can only not replace models if no beam
    for(jks in start_pos:length(self$incls)){
      self$updateJk(jks)
      if(run_sep){
        inds_run = 1:length(self$models)
        #cnt1=0
        finished_all = F
        while(!finished_all && self$cnt1<self$maxvars_jks)
          {
          #for(ki in 1:length(models)){
         # aa=parLapply(cl, inds_run,function(ki) self$runAnalysisInnerSep(ki) )
          self$m1 =mclapply(inds_run, function(ki){ 
            print(ki)
            self$runAnalysisInnerSep(ki)
          }, mc.cores=getOption("mc.cores.runsep",1))
          self$cnt1 = self$cnt1+1
          self$finished=unlist(lapply(self$m1, function(m2) is.null(m2) || length(m2$to_keep)==0))
          print(self$finished)
          sort(table(unlist(lapply(self$models, function(m)(names(m$prev))))))
          
          
          rmsv_cv= .summarise(.calcRMSVAll(self$datas, NULL,cv=T, label="cv"),c("beam","pheno"))
          if(sum(rmsv_cv$value)>=sum(self$rmsv_cv0$value) &&  useInternalCVAsStopping){
            print("stop cv")
            .trim(self$models, self$datas,c(),inds_run )
            break;
          }else{
            self$rmsv_cv0 = rmsv_cv
          }
          self$rms_list_all[[length( self$rms_list_all)+1]] = self$saveDatas1()
          print(self$rms_list_all[[length(self$rms_list_all)]])
          if(getOption("printModel",FALSE)) try(self$saveDatas())
          if(!is.null(outpdf3) & self$cnt1>1){
            print(paste("plotting",outpdf3, self$cnt1))
            ggp=self$plotRMS()
            try(ggsave(outpdf3, plot=ggp, width =45, height =45, units = "cm",limitsize=F))
          }
          #print(self$rms_list_all[[length(self$rms_list_all)]])
          finished_all = if(length(self$m1)==1) self$finished[1] else length(which(self$finished[-length(self$m1)]))>=0.5*length(self$m1[-length(self$m1)])
        }
        #    .printModel(models, datas, datas1, rdsdir)
      }else{
        self$runAnalysisInner()
      }
    }
    lapply(self$datas, function(d)d$saveParquet())
    invisible(self$models)
  },
  runAnalysisInnerSep=function(ki){
   
    if(self$finished[ki]) return(NULL)
    rms_prev = self$m1[[ki]]$rms_prev
    to_keep=self$m1[[ki]]$to_keep
    skipped = self$m1[[ki]]$skipped
    finished=F
    cv_skips_allowed=getOption("cv_skip",0) ;  signatures = getOption("signatures",NULL)
    useInternalCVAsStopping=  T #getOption("useInternalCVAsStopping",F)
    orderCV = getOption("orderCV",T)
    pv_only=getOption("fspls.pval_only",T)
    beam=getOption("fspls.beam",c(1,1))
    logpthresh=log(getOption("fspls.pthresh1",0.05))
    replaceModel = F #beam[1] *beam[2] >1  || !useInternalCVAsStopping  
    models = self$models
    datas=self$datas
    #cnt1 =models[[ki]]$cnt
    inds_m = ki  #if(useInternalCVAsStopping) 1:length(models) else length(models)
    pvslist_all =lapply(inds_m, function(k){
      model = self$models[[k]]
      nxt_var = if(models[[k]]$cnt<length(signatures) ) signatures[[models[[k]]$cnt+1]] else NULL
      model$simpleForwardTrain(datas,k,  exclude=exclude,nxt_var = nxt_var, incl = self$incl ,  weights=NULL, to_keep=to_keep)
    })
    pvslist_all1 =pvslist_all[[length(pvslist_all)]]$pv
    pvslist_all2 = pvslist_all[[length(pvslist_all)]]$cumpv
    #print(pvslist_all1)
    pvtk1=  pvslist_all1<logpthresh
    if(length(which(pvtk1))==0){
      .trim(models, datas,c(),inds_m  )
      print("breaking on pvalue")
      finished=T
      return(NULL)
    }
    
    if(!pv_only &useInternalCVAsStopping  & FALSE){
      rms_cv1 = .merge1(lapply(datas,function(d) .summarise(d$getRMSVAll())),num_cols="value",addName="dataset")
      rms_prev1=min(rms_cv1$value)
      if(rms_prev1>=rms_prev){
        .trim(models, datas,c(),inds_m )
        finished=T
        return(NULL)
      }
    }
    vars = names(pvslist_all2)
    # prev =self$models[[length(models)]]$prev
    prev =self$models[[ki]]$prev
    #if(length(models)>1 && replaceModel) {
    #  .updateModels(models, datas, prev,1:(length(models)-1),useInternalCVAsStopping)
    #}
    if(!pv_only){
      rms_cv = .merge1(lapply(datas,function(d) .summarise(d$getRMSV(ki,within=T))),num_cols="value",addName="dataset")
      #if(length(models)==1){
      #  rms_cv = .merge1(lapply(datas,function(d) .summarise(d$getRMSV(length(models)))),num_cols="value",addName="dataset")
      #}else{
      #  rms_cv = .merge1(lapply(datas,function(d) .summarise(d$getRMSVAll())),num_cols="value",addName="dataset")
      #}
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
      #.trim(models, datas,c()  )
      .trim(models, datas,c(), inds_m=inds_m  )
      finished=T
      return(NULL)
    }
    if(!pv_only){
      #  if(!useInternalCVAsStopping || length(models)==1) rms_prev1=min(rms_cv$value[to_keep],na.rm=T)
      rms_prev1 = min(rms_cv$value[to_keep], na.rm=T)
      if(rms_prev1>rms_prev){
        skipped = skipped+1
        print(paste("skipped ",skipped))
      }else{
        #reset skipped and reset bar
      #  print("reset skipped")
        skipped=0
        rms_prev = rms_prev1
      }
      
      if(skipped> cv_skips_allowed){
        if(length(to_keep)>0){ ## need to wind back since rms got worse
          # for(k in 1:length(models)){
          for(k in inds_m){
           self$models[[k]]$keep(c())
            lapply(datas, function(d) d$train[[k]]$keep(c()))
          }
        }
        finished=T
        #print(paste("breaking here because cv detetiorating",self$models[[1]]$cnt))
        
        #  for(k in which(!finished)){
        #   self$models[[k]]$unwind(datas,k) 
        #  }
        return(NULL)
      }
    }
    if(length(to_keep)==0){
      finished=T
      return(NULL)
    }
    to_keep =1:min(length(vars),beam[2])  ## really only comes into play if beam>1
    attr(models,"to_keep")=to_keep
    .reorder(models, datas, o, inds_m=inds_m)    
    #return()
    #.printModel(models, datas, datas1, rdsdir)
    
    #          attr(models,"to_keep")=to_keep
    #}
    return(list(to_keep=to_keep, rms_prev=rms_prev, skipped=skipped))
  },
##FOLLOWING FUNCTION GETS LIST OF DRUGS WHICH HAVE lower CI of correlation >0
getList=function(drugs,phenos, phenos1, n=1){
  options("maxvars"=n)
  self$maxvars=n
  datas = self$datas
  datas1 = self$datas1
  rms_all = lapply(1:length(drugs), function(i1){
    drug1 = drugs[i1]
    drug = names(drugs)[i1]
    print(paste(drug,i1))
    names(drug1) = drug
    self$updatePheno(phenos, drug, phenos1 = phenos1, drug1 = drug1)
    v1 = var(datas1[[1]]$y,na.rm=T)
    if(is.na(v1) || v1<1e-9) return(NULL)
    self$runAnalysisAll()   
   
    #sort(table(unlist(lapply(self$models, function(m) names(m$prev)))))
   
    rms11= try(self$saveDatas1())
    if(inherits(rms11,"try-error")){
     return (NULL) 
    }
    prev=datas[[1]]$train[[length(datas[[1]]$train)]]$prev[[1]]
    betas = apply(apply(prev$betas[[1]],c(1,2),round,2),2,paste,collapse=",")
    res0=rbind(
    cbind(rms11$crossvalidation[grep('low',rms11$crossvalidation$subpheno),,drop=F],betas),
    cbind(rms11$validation[grep('mid',rms11$validation$subpheno),,drop=F],betas))
    print(res0)
    res0
  })
  names(rms_all) = drugs
  rms_all = rms_all[!unlist(lapply(rms_all, is.null))]
  rms_all1 = .merge1(rms_all,  num_cols="value")
  rms_all1
},
plotRMS=function(){
  .plotRMS(self)
  },
getFullList=function(drugs, phenos ,phenos1,max=10){
  n=1
  drugl = list()
  drugl_val =list()
  while(length(drugs)>0 && n<=max){
    ab0 = self$getList(drugs,phenos, phenos1, n=n)
    ab = subset(ab0, label=="cv")
    ab_val = subset(ab0, label=="validation")
    o=order(ab$value)
    ab = ab[ o,,drop=F]
    ab_val = ab_val[ o,,drop=F]
    drugl[[n]] = ab
    drugl_val[[n]]=ab_val
    varl=unlist(lapply(ab$beam, function(b)length(strsplit(as.character(b),",")[[1]])))
    drugs2 = ab$pheno [(which(ab$value<0 & varl>=n))]
    head(ab[order(ab$value, decreasing=F),c(2,4,5)])
    if(length(drugs2)==0) break;
    n=n+1
    drugs2 = unique(as.character(drugs2))
    drugs=drugs[which(drugs %in% drugs2 | names(drugs) %in% drugs2)]
  }
  names(drugl) = paste("cv",1:length(drugl),sep=".")
  names(drugl_val) = paste("val",1:length(drugl_val),sep=".")
  c(drugl, drugl_val)
}
  
)
)
