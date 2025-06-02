.refactor<-function(pvs_all1){
  pvn = names(pvs_all1[[1]])
  names(pvn)=pvn
  lapply(pvn, function(pvn1){
    lapply(pvs_all1, function(v) v[[pvn1]])
  })
}


modelObj<-R6Class("modelObj",##represents a state of the model
                  public = list(
                    prev="list",  prev_old="list",
                    #params="list",
                    pvs_all="list",
                    gaussians="vector",
                    weights_all="list", 
                    data_types = "vector",
                    pvs_all_cum="list",
                    rmsv_prev="numeric", rmsv_all="list",rmsv_prev_matr = "matrix",
                    rmsv_="matrix",rmsv_all1="list",
                    cnt="numeric",cnt_skipped="numeric",finished="logical",
                    const_all="list",
                   # sqnorm="list",
                    initialize = function(data_names, data_types,  rmsv_prev,var = list()){
                     # self$sqnorm = lapply(datas, function(d){
                      #    d$train$sqnorm
                      #})
                       self$data_types =  data_types # names(datas[[1]]$data)  #this is pro
                       self$pvs_all_cum = list()
                     
                      self$gaussians=rep("gaussian", length(data_names))
                      #  self$prevs = lapply(datas[[1]]$prevs, function(p) lapply(p, function(p1)list(var=p1$var)))
                      prevs1 = list(stateObj1$new(var=var))
                      names(prevs1) = lapply(prevs1, function(v)v$name)
                      self$prev=prevs1
                      self$prev_old = NULL
                    #  self$params = params
                      self$pvs_all = list()
                      self$weights_all =  list() ;
                      self$const_all = list();
                    
                      self$rmsv_prev = rmsv_prev
                      self$cnt=0;
                      self$rmsv_all = list();
                      self$rmsv_all1 = list();
                      self$finished=F
                      self$cnt_skipped=0
                    },
                    #comb holds weights and rmsv 
                    ## diversity is between 0 and 1 and controls how diverse the sigs are
                    #prep=prepN_; ypreds = prep$ypreds
                   reorder=function(o){
                  #   o1=match(names(o), names(self$prev))
                     self$prev = self$prev[o]
                   },
                    simpleForwardStep=function(datas,k,incl, best_all_i=NULL, 
                                               exclude=list(),
                                               ntrans=1,
                                               to_keep=1,
                                               fspls.sum= T ##evaluate sum of angles over phenos
                    ){
                      model=self
                      {
                        if(is.null(incl)) stop("incl should not be null")
                     # params=model$params
                        beam = getOption("fspls.beam",c(1,1))
                      n=beam[1] ## how much best angle options to keep 
                       mult=rep(1, length(datas))
                      if(length(mult)<length(datas)) stop("!! length mult problem")
                      projection =TRUE
                      betas_projection=TRUE
                      calcProjection=TRUE
                      calcOffset=FALSE
                     # family=datas[[1]]$family
                      best_j = NULL
                     # if(length(family)<length(datas)) stop("!! length family not right")
                      t = Sys.time()
                      data_types=model$data_types
                      names(data_types)=data_types
                      prev = model$prev
                      }
                      if(is.null(best_all_i)){
                          best_all_i=lapply(to_keep, function(j){
                            prev_i = prev[[j]]
                            var = prev_i$var
                            type=getOption("type","slow1")
                            angles = lapply(datas,function(d) d$getAngles1(var,incl,k, type))  #THIS IS CALL 1 to data.  should do all CV calls to data at once
                            names(incl)=incl
                            angle1=lapply(data_types,function(inc1){
                              a1=lapply(angles, function(angle_d) {
                                a2=angle_d[[inc1]]
                                if(is.null(a2)) return(NULL)
                                df2 = data.frame(lapply(a2, function(a3){
                                  apply(a3,2,min) ## minangle over subphenos (e.g. within each multinomial), could also be sum
                                }))
                                apply(df2,1,min)   ## min angle over phenos could also consider sum
                                })
                              .sumMatrices(a1) ## sum over data sets
                            })
                            best_i = whichpart1(angle1,n=n*ntrans,getOption("fspls.one_for_each",F))
                           # bestv=angle1[[best_i[[1]][1]]][best_i[[1]][2]]
                          #  if(is.na(bestv) || bestv==999) stop("something went wrong")
                            best_i
                          })
                      }
                      best_all=lapply(1:length(best_all_i), function(j){
                        prev_i= prev[[to_keep[j]]]
                        best_i = best_all_i[[j]]
                        nxt_v=lapply(best_i,function(b_i) stateObj1$new( prev_i, b_i))
                        names(nxt_v) = unlist(lapply(nxt_v,  function(nv)paste(unlist(lapply(nv$var, function(vv){
                          #      paste(names(todoInds)[[vv[1]]], names(todoInds[[vv[1]]])[vv[2]],sep=".")
                          #varnames=lapply(vars0, function(vv) c(names(datas[[1]]$data)[vv[1]],dimnames(datas[[1]]$data[[vv[1]]])[[2]][vv[2]]))
                          paste(vv[[1]],vv[[2]],sep=".")
                        })), collapse=",")))
                        nxt_v
                      })
                      model_prev = model$prev
                      if(is.null(model_prev)) stop("is null model_prev")
                        if(getOption("update_data",FALSE)) invisible( lapply(datas, function(d) d$updateModel(k,best_all_i, model_prev, to_keep))) #THIS IS CALL 2 to data to update the model 
                        
                      names(best_all) = NULL
                      b_all1 = unlist(best_all, rec=F)
                      #print(b_all1)
                      b_all1
                    },
                    ##FOLLOWING TRAINS THE MODEL ONE STEP
                    simpleForwardTrain=function(datas, k,
                                                exclude=list(),
                                                incl = names(datas[[1]]$data),
                                                nxt_var = NULL,
                                                to_keep = 1,
                                                weights=NULL  #prep$weights
                    ){
                      model=self
                      {
                      fspls.sum = getOption("fspls.sum",T)
                      beta_v= NULL #if(DRS) .makeVMatr(getOption("beta_vp",2),getOption("beta_v",c(0,1,2,3))) else NULL
                      #ycol=prep$ycol
                      time0=Sys.time()
                      prev =model$prev
                      weights_all = model$weights_all
                      const_all = model$const_all
                      pvs_all = model$pvs_all
                      rmsv_all= model$rmsv_all
                      rmsv_all1 = model$rmsv_all1
                      if(is.null(rmsv_all1)) rmsv_all1 = list()
                      cnt=length(rmsv_all)
                      cnt_skipped = model$cnt_skipped
                      if(is.null(cnt_skipped)) cnt_skipped = 0;
                      rmsv_prev =   model$rmsv_prev  
                      wname=names(weights)
                     # params = model$params
                      #if(is.null(params$min)) params$min=0
                        best_all_i = NULL
                     
                        if(length(nxt_var)>0){
                          best_all_i = lapply(nxt_var, function(nxt_var1){  
                          bi1 = strsplit(nxt_var1, "\\.")[[1]]
                          di =which(names(datas[[1]]$data)==bi1[[1]])
                          gne_nme = paste(bi1[-1], collapse=".")
                          gne_nmes = dimnames(datas[[1]]$data[[bi1[[1]]]])[[2]]
                          c1=which(gne_nmes==gne_nme)
                          if(length(c1)==0){
                            c1 = grep(gne_nme,gne_nmes )
                            if(length(c1)>1){
                              c1=c1[which(unlist(lapply(gne_nmes[c1],function(cc) strsplit(cc,"\\s+")[[1]][[1]]))==gne_nme)]
                            }
                          }
                          if(length(di)!=1 || length(c1)!=1) stop(nxt_var)
                          list(c(di,c1))
                        })
                        }
                      }
                        # prev=prevs1[[k]]
                        prev=model$simpleForwardStep(datas,k,incl,
                                                best_all_i = best_all_i,
                                                fspls.sum=fspls.sum,
                                                to_keep=to_keep,
                                                  exclude=exclude)
                       
                        pvs_all1=.refactor(lapply(datas, function(d){
                          lapply(prev, function(prv){
                            b_i1 = prv$var[[length(prv$var)]]
                            prev_var =  prv$var[-length(prv$var)]
                            sig_res = d$calcBetaProj1(k, b_i1, prev_var )
                            sig_res$pvs
                          #d$train[[k]]$getPvs(prev)
                          })
                        }))
                        #  pvs_all1 = lapply(prevs, .getWeights11, names(datas), pvs=T)
                        #pvs_all1 = data.frame(pvs_all1)  
                       
                        if(getOption("spike_slab_iter",0)>0){
                          pvslist = lapply(pvs_all1, function(pvn1){
                            min(unlist(pvn1))
                          })
                          #   print("pvslist")
                          #  print(head(pvslist))
                          #tokeep1 = pvslist>getOption("fspls.pthresh1",1.0) 
                          
                        }else{
                          pvslist = lapply(pvs_all1, function(pvn1){
                            .combinePv(unlist(pvn1))
                          })
                         # tokeep1 = pvslist<getOption("fspls.pthresh1",1.0) 
                        }
                        model$cnt = model$cnt+1
                        model$prev_old=model$prev
                        model$prev=prev
                        nmes1 = names(datas[[1]]$train[[k]]$prev)
                        names(nmes1) = nmes1
                        if(FALSE){
                            pvs_cum=unlist(lapply(nmes1, function(nme1) {
                              sum(unlist(lapply(datas, function(d){
                                p1 = d$train[[k]]$prev[[nme1]]
                                sum(unlist( lapply(p1$cum_pvs_proj, function(p2){
                                  -1*qchisq(p2, df=1, lower.tail=F,log.p=T)
                                })))
                              })))
                            }))
                            pvslist1 = unlist(pvslist)
                            pvs_cum = pvs_cum[match(names(pvslist1), names(pvs_cum))]
                            model$pvs_all_cum [[model$cnt]]=  pvs_cum
                        }
                        model$pvs_all[[model$cnt]]=unlist(pvslist)
                      return( list(pv=pvslist1, pvs_all1 = pvs_all1))
                    },
filter=function(datas, k,
                            exclude=list(),
                            incl = names(datas[[1]]$data),
                            nxt_var = NULL,
                            weights=NULL  #prep$weights
){
  model=self
                      lens=sort(table(unlist(names(prev))), decr=T)
                      if(length(lens)==0){
                        model$finished=T
                        if(getOption("logprint",F)) print("breaking 0")
                        return(model$finished)
                      }
                      
              #        vars=lens #.getRank(prevs,lens)
              #        vars = .keepUniq(vars)
                      beam =getOption("fspls.beam",c(1,1))
                      maxn1= min(length(vars), beam[2]*beam[1])
                      if(FALSE){
                        #diversity=getOption("fspls.diversity",c(1.0,1.0)
                        o1=.applyDiversity(vars,cnt:1,diversity, maxn1, def_val=9999)
                        #print(head(cbind(names(vars), names(vars)[o1]),20))
                        vars = vars[o1]
                      }
                      #,max(unlist( lapply(prevs, length))))
                      if(FALSE &&length(vars)>length(ypreds)){
                        ##this is because hard to precalcute how many we need when use weights
                        # print(vars)
                        print(paste("warning not enough ypreds", length(ypreds), length(vars)))
                        ypreds=.initYpred(datas, length(vars))
                      }
                      #maxl = min(length(vars), length(ypreds))
                      #tokeep=names(vars[1:maxl])
                      tokeep = names(vars)
                      
                      #rmsv_=.merge1(lapply(datas, function(d)d$rmsv_),num_cols="value",addName="data")
                      
                      rmsv_=.merge1(lapply(datas, function(d)d$getRMSV(k)),num_cols="value",addName="data")  #this is call 3 to datas
                      rmsv_prev1 = .getRMSPrev(rmsv_)
                      beam_levs = levels(rmsv_$beam)
                      names(beam_levs) = beam_levs
                      rmsv1=unlist(lapply(beam_levs, function(beam){
                        ab=rmsv_[rmsv_$beam==beam,,drop=F]
                        sum(ab$value,na.rm=T)  #could also try min ? 
                      }))
                      rmsv1 = rmsv1[match(rmsv_$beam, names(rmsv1) )]
                     
                      #rmsv1 = .findMinRMSV(rmsv_,mult, fspls.sum=fspls.sum)
                     # rmsv_prev1 = .findMinRMSV(rmsv_,mult, fspls.sum=fspls.sum)
                      rmsv_prev =    model$rmsv_prev
                      rms_diff = apply(cbind(rmsv_prev1,rmsv_prev), 1, function(v) v[2]-v[1])
                      diff_thresh = getOption("diff_thresh",-0.001)
                      
                      rmsv_all[[cnt]]=rmsv1 
                      pvs_all[[cnt]] = pvs_all1
                      attr(rmsv_all[[cnt]],"phen")=attr(rmsv1,"phen")
                      rmsv_all1[[cnt]] = rmsv_
                    #  rms_prev1 = min(rmsv1) #min( apply(rmsv,1,min))
                      ###WARNING NEXT BIT EXECUTES LOTS OF IFELSE
                      if(inherits(rmsv_prev1,"try-error")){
                        print(paste("problem calculating rms_prev1 ",cnt))
                        print(paste(" need to break"))
                        model$finished=T 
                        lapply(datas, function(d) d$train[[k]]$keep(c()))
                        
                        if(getOption("logprint",F)) print("breaking 1")
                        return(model$finished)
                      }else if(length(which(!is.na(rmsv1)))==0 ){
                        if(getOption("logprint",F)) print("breaking 2")
                        lapply(datas, function(d) d$train[[k]]$keep(c()))
                        model$finished=T
                        return(model$finished)
                      }else if(is.null(nxt_var) && min(rms_diff,na.rm=T)<diff_thresh && cnt>=getOption("min",0) && is.null(weights) && cnt_skipped >= getOption("fspls.max_skipped",0)){
                        model$finished=T
                        lapply(datas, function(d) d$train[[k]]$keep(c()))
                        # print(paste(rms_prev1, rmsv_prev, cnt, params$min))
                        if(getOption("logprint",F)) print("breaking 3")
                        return(model$finished)
                      }else{
                        if(min(rms_diff,na.rm=T)<diff_thresh ){
                          model$cnt_skipped = model$cnt_skipped+1
                        }
                        o = order(rmsv1)
                        o1 = match(names(rmsv1[o]), names(prev))
                       
                        rmsv1 = rmsv1[o]
                        rmsv_ = rmsv_[o,]
                        
                        prev = prev[o1]
                        lapply(datas, function(d) d$reorder(o1,k))
                        beam = getOption("fspls.beam",c(1,1))
                        maxn1= min(length(prev), beam[2])
                        if(cnt>1 && FALSE){
                          #print("checking diversity")
                          ##diversity=getOption("fspls.diversity",c(1.0,1.0)
                          o=.applyDiversity(rmsv1,cnt:1,diversity, maxn1, def_val=9999)
                          
                          maxn1 = min(maxn1, length(o))
                        }
                        
                        # m1=cbind(order(rmsv1[o]),order(lens[o]), order(vars[o]))
                        #  dimnames(m1)=list(names(rmsv1[o]),c("rmsv","len","vars"))
                        #print(head(rmsv1[o]))
                        #print("m1")
                        #   print(head(m1))
                        #  print(head(rmsvo,params$beam[3]))
                        
#                          p1= prev[match(names(rmsv1)[o],names(prev))[1:maxn1]]
                      #    prev=p1[!unlist(lapply(p1, is.null))]
                        #})
                        #prevs = prevs2
                        cnt=cnt+1
                        
                        #print(head(rmsv[o],params$beam[1]))
                       
                        
                        nxt_var_keep = rep(0, length(pvslist))
                        if(!is.null(nxt_var)){
                          #to ensure we keep the nxt_var
                          last_var =  unlist(lapply(names(pvslist), function(st)rev(strsplit(st,",")[[1]])[[1]]))
                          nxt_var_keep[grep(nxt_var,last_var)]=1
                        }
                        tokeep1 = 
                          which(
                          nxt_var_keep | (tokeep1 & names(tokeep1) %in% names(rmsv1)[1:maxn1] ) )#only keep good pvalues and top rmsv
                        #print(paste("nxt_var", nxt_var))
                        if(length(tokeep1)>0){
                          rmsv2 = rmsv1[tokeep1]
                          if(rmsv2[1] >rmsv1[1]) tokeep1 = c()
                        }
                        #print(paste(names(pvslist)))
                        lapply(datas, function(d) d$train[[k]]$keep(tokeep1))  #this is call 4 to datas, saying what to keep
                        if(length(tokeep1)==0){
                        #  print("to keep is length zero")
                        #  print(pvs_all1)
                        #  print(pvslist)
                          # print(pvslist)
                          model$finished=T
                         # if(getOption("logprint",F)) print(pvslist)
                          # print(paste(rms_prev1, rmsv_prev, cnt, params$min))
                          if(getOption("logprint",F)){
                            print("breaking 4")
                            print(pvslist)
                          }
                          return(model$finished)
                        }
                        prev=prev[tokeep1] #lapply(prevs, function(prev) prev[tokeep1])
                        pvs_all1=.refactor(lapply(datas, function(d)d$train[[k]]$getPvs()))
                        
                      #  pvs_all1 = .getWeights11(prev,names(datas), pvs=T) # lapply(prevs, .getWeights11, names(datas), pvs=T)
                        pvs_all1 = lapply(pvs_all1, unlist,rec=F)   
                        
                      if(FALSE){ ## we dont really needs this
                        weights1 = .refactor(lapply(datas, function(d)d$getBetas()))
                        weights1 = lapply(weights1, unlist,rec=F)   
                        
                        const1 = .refactor(lapply(datas, function(d)d$getConstants()))#.getWeights11(prev,names(datas),const=T) # lapply(prevs, .getWeights11,names(datas), const=T)
                        const1 = unlist(const1, rec=F)
                        
                        weights_all[[cnt]] = weights1
                        const_all[[cnt]] = const1
                      }
                       
                        
                        #weights1=unlist(weights_all,rec=F)
                        model$weights_all = weights_all
                        model$const_all= const_all
                        model$pvs_all = pvs_all
                        model$rmsv_prev=rmsv_prev1
                        model$cnt=cnt;
                       # model$cnt_skipped = cnt_skipped
                        model$rmsv_ = rmsv_
                        
                        model$rmsv_all = rmsv_all
                        model$rmsv_all1  =rmsv_all1
                        model$finished=F
                        model$prev_old = model$prev
                        model$prev = prev
                        model$rmsv_prev_matr = rmsv_prev1
                        
                        return(model$finished)
                      }
                      model$finished=T
                      #model=list(rmsv_all = rmsv_all, 
                      #           weights=weights,
                      #           weights_all = weights_all, params=params, prevs=prevs1,aucT=.rmsvToTable(rmsv_all), weights1=unlist(weights_all,rec=F))
                      return(model$finished)
                    },
keep=function(tokeep1){
  if(length(tokeep1)==0){
    #print("winding back")
    self$prev = self$prev_old;
    self$prev_old=NULL
    self$cnt=self$cnt-1
  }else{
    self$prev=self$prev[tokeep1]
  }
},
                    update=function(prevb){
                      self$prev = lapply(prevb, function(x) x$clone())
                      self$cnt = length(prevb[[1]]$var)
                    },
                    unwind=function(datas,k){
                      tokeep1 = c()
                      lapply(datas, function(d) d$train[[k]]$keep(tokeep1))
                      self$prev = self$prev_old;
                      self$prev_old=NULL
                      self$cnt = self$cnt-1
                    } 
                    
                  )
)





##k is which iteration
