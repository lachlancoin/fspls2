#private = self[[".__enclos_env__"]]$private

mergeAll = function(comb2_new, beam){
  comb1 = .merge1_new(lapply(comb2_new, function(c){
    .merge1_new(lapply(c$angles, function(c1){
    
      .merge1_new(c1,addName="pow")
    }),addName="func")
  }), addName="prev") 
   #comb1$prev 
   o =  order(comb1$value)
   comb1[o[1:beam],]
  
}

.check_data<-function(data, y){
  non_na = apply(y,2,function(x) length(x[!is.na(x)]))
  if(min(non_na)==0) warning("one y column has all NA")
  rn_y = rownames(y);
  if(length(rn_y)==0) warning(" define rownames for y");
  rns = lapply(data, function(d){
    rnd = rownames(d);
    if(length(rnd)==0) warning(" define rownames for each data object");
    mi1 = match(rnd, rn_y)  ;mi2 = match(rn_y,rnd)  ;
    if(length(which(is.na(mi1)))>0) warning(" cannot match rownames from data  to y")
    if(length(which(is.na(mi2)))>0) warning(" cannot match rownames from data to y")
    if(mi1[!is.na(mi1)][1]!=1) warning("not a match")
  });
  
}

get_comb_plot<-function(plot_results, nvar){
  if(length(plot_results)==0) return(NULL)
  names(plot_results) = 1:length(plot_results)
  plot_results[[1]]$prev=""
  comb_plot = .merge1_new(plot_results, addName="nvar");
  comb_plot = comb_plot[comb_plot$func!="rand",,drop=F]
  cumulative=comb_plot$value + comb_plot$sumAngle
  
  sigs = comb_plot$signature; 
  maxsig=as.factor(unlist(lapply(sigs, function(sig1){
    inds1=grep(paste0("^",sig1), comb_plot$signature)
    #    inds1 = inds1[]
    inds1 = inds1[which.min(cumulative[inds1])]
    comb_plot$sig[inds1]
  })))
  comb_plot = comb_plot|>tibble::add_column(cumulative, totalvar = nvar, maxsig)
  comb_plot$nvar = as.numeric(comb_plot$nvar)
  comb_plot
}
plot_traj_all<-function(comb_plot, y="value", facet="maxsig~.", keep_best=10, txtsize=5, step=2){
  nreps = unique(comb_plot$nrep); names(nreps)=nreps;
  lapply(nreps,function(nr){
    comb_plot1 = subset(comb_plot, nrep==nr)
    plot_traj(comb_plot1,y,facet, keep_best=keep_best,txtsize=txtsize, step = step)
  })
}
plot_traj<-function(comb_plot, y="value"  ,facet=".~maxsig", keep_best=10, txtsize=5, step=2){ #y="cumulative";
  if(!is.null(comb_plot$nrep)){
    if(length(unique(comb_plot$nrep))>1) stop(" need to subset on nrep first")
  }
  comb_plot$maxsig = gsub("\n",";", comb_plot$maxsig)
  maxl = min(comb_plot[[y]])
  
  maxsigl1 = unlist(lapply(comb_plot$maxsig, function(x) paste(sort(strsplit(x,";")[[1]]), collapse=";")))
  maxsigl1_u = sort(table(maxsigl1), decreasing=T)
  conv = lapply(names(maxsigl1_u), function(x){
    names(sort(table(comb_plot$maxsig[which(maxsigl1==x)]),decreasing=T))[1]
  })
  names(conv) = names(maxsigl1_u)
  maxsig = unlist(lapply(maxsigl1, function(x) conv[x]))
  comb_plot$maxsig = maxsig
  maxsigl = sort(table(maxsig),decreasing=T)
  
  if(length(maxsigl1_u)>keep_best){
    
   ms1   = factor(maxsig, names(maxsigl)[maxsigl>=maxsigl1_u[keep_best]])
    comb_plot1 =comb_plot[!is.na(ms1),,drop=F]
  }else{
    comb_plot1 = comb_plot
  }
   ms1 = unique(comb_plot$maxsig); names(ms1) = ms1
  levs = sort(unlist(lapply(ms1, function(ms){
    inds1 = comb_plot$maxsig==ms
    min(comb_plot$cumulative[inds1])
  })),decreasing=F)
  labels = unlist(lapply(names(levs), function(lev){
      maxs = strsplit(lev,";")[[1]];
      paste(unlist(lapply(seq.int(1, length(maxs), by=step), function(st){
      
        mi = min(length(maxs), st+step)
        paste(maxs[st:mi], collapse=";")
      })), collapse="\n")
  }))
  
  comb_plot1$maxsig = factor(comb_plot1$maxsig, levels = names(levs), labels=labels)
  comb_plot1[[y]] = -1*comb_plot1[[y]]
  ggp=ggplot(comb_plot1, aes_string("nvar" ,y, color="maxsig"))+geom_point()+geom_line(alpha=.1)+facet_grid(facet)
  ggp=ggp+guides(color = "none")+theme(
    strip.text = element_text(size = txtsize, face = "bold", color = "black")
  )+geom_hline(yintercept = -1*maxl)+scale_y_log10()
  ggp
}
plot_ri=function(comb2_new, alpha =.5){
  ri = comb2_new[[1]]$pvs
  comb_ = comb2_new[[1]]$angles
  all_angles = attr(comb_,"all")
  typs = names(ri)
#  typ1 = typs[[1]]; ni1 =names(ri[[typ1]])[[1]] 
df1 = .merge1_new(lapply(typs, function(typ1){
    ni = names(ri[[typ1]]); 
    .merge1_new(lapply(ni, function(ni1){
      
      if(FALSE){
      pval = unlist(lapply(ri[[typ1]][[ni1]], function(x)  .sumChisq(x$pvs)))
      angle = comb_[[typ1]][[ni1]]$value
      comb22=cbind( comb_[[typ1]][[ni1]]$names, data.frame(cbind(angle,pval)),typ1, ni1)
      }else{
        pval_all = unlist(lapply(ri[[typ1]][[ni1]], function(x)  (x$pvs)))
        
        ang_all = all_angles[[typ1]][[ni1]]
        mat = t(as.matrix(ang_all))
        triples <- mat |> 
          as.data.frame() |> 
          tibble::rownames_to_column(var = "Row") |> 
          pivot_longer(
            cols = -Row, 
            names_to = "Column", 
            values_to = "angle"
          )
        comb11=cbind(triples, pval_all, typ1, ni1)
        names(comb11)=gsub("pval_all","pval", names(comb11))
        comb11
      }
      comb11
    }))
    
  }))
df2=pivot_wider(df1[,!(names(df1)=="angle")],names_from = c("Column"), values_from = pval)
df3=pivot_wider(df1[,!(names(df1)=="pval")],names_from = c("Column"), values_from = angle)
df22 = cbind(df2[,1:3],apply(df2[,-(1:3)],1,.sumChisq))
df33 = cbind(df3[,1:3],apply(df3[,-(1:3)],1,sum))
names(df33)[4] = "angle"
df4 = cbind(df33, df22[,4])
names(df4)[5] = "pval"
my_colors <- grDevices::colorRampPalette(brewer.pal(12, "Paired"))(100)
ggp1=ggplot(df1, aes(angle,pval, color=Column, shape=Row))+geom_point(alpha = alpha,size=2)+facet_grid("typ1 ~ni1")+ scale_shape_manual(values = rep(0:25, length.out = 100))+scale_color_manual(values = my_colors)
ggp2=ggplot(df1, aes(angle,pval, color=Row, shape=Column))+geom_point(alpha = alpha,size=2)+facet_grid("typ1 ~ni1")+ scale_shape_manual(values = rep(0:25, length.out = 100))+ scale_color_manual(values = my_colors)



ggp3 = ggplot(df4, aes(angle,pval, color=typ1, shape=ni1))+geom_point(alpha = alpha)
return(list(ggp1, ggp2, ggp3))

}

.extrAngles<-function(angleH,comb_all2, incl){
  angles1=angleH$angles;cols_incl1=angleH$cols_incl 
  nme_trans = names(angles1[[1]][[1]]); names(nme_trans) = nme_trans
  names(incl) = incl
  ang_extracted=lapply(nme_trans, function(nme_t1){
    nme_pow = names(angles1[[1]][[1]][[nme_t1]]); names(nme_pow)=nme_pow
    lapply(nme_pow, function(nme_p1){
      ea1 = lapply(incl, function(inc1){
      ang1 = angles1[[inc1]]
      if(is.null(ang1)) return(NULL)
      col_incl = cols_incl1[[inc1]]
      ang2=ang1[[1]][[nme_t1]][[nme_p1]]
      mi2 = match(comb_all2[[nme_t1]][[nme_p1]]$names,colnames(ang2))
      ang2[,mi2,drop=F]
      })
      ea_all = ea1[[1]]
      if(length(ea1)>1){
        for(kk2 in 2:length(ea1)){
          ea_all = ea1[[kk2]]+ea_all
        }
      }
      ea_all
    })
  })
  ang_extracted
}

.combineAngles1<-function(angleH, incl, flags,sumAngle, prev_signature, excl=list()){ 
  topn = .readFlag(flags,'topn', 20)
  onlyAll = .readFlag(flags,'only_all',F)
  angles1=angleH$angles;cols_incl1=angleH$cols_incl 
  nme_trans = names(angles1[[1]][[1]]); names(nme_trans) = nme_trans
  # nmes_angs1 = names(angles1); names(nmes_angs1)=nmes_angs1
  #nme_t1 = nme_trans[[1]]; nme_p1 = names(angles1[[1]][[1]][[nme_t1]])[[1]]; inc1 = incl[[1]]; jk=1
  names(incl) = incl
  comb_all2=lapply(nme_trans, function(nme_t1){
    nme_pow = names(angles1[[1]][[1]][[nme_t1]]); names(nme_pow)=nme_pow
    lapply(nme_pow, function(nme_p1){
      comb_all=lapply(incl, function(inc1){
        ang1 = angles1[[inc1]]
        if(is.null(ang1)) return(NULL)
        col_incl = cols_incl1[[inc1]]
        ang2=ang1[[1]][[nme_t1]][[nme_p1]]
        cs = Matrix::colSums(ang2)
        if(length(ang1)>1){
          for(jk in 1:length(ang1)){
            cs = cs+Matrix::colSums(ang1[[jk]][[nme_t1]][[nme_p1]])
          }
        }
        excl1 = excl[unlist(lapply(excl, function(ex) ex[3]==nme_t1 && ex[1] == inc1 && ex[4] ==nme_p1))]
        if(length(excl1)>0){
          col_incl[which(names(col_incl) %in% unlist(lapply(excl1, function(ex) ex[2])))]=F
        }
        cs[col_incl]
      })
      top_angles=whichpart1(comb_all, n=topn, return_scores=T)
      t1 = .merge1_new(lapply(top_angles, function(ta){
        data.frame(list(names = names(ta), value=ta))
      }),addName="data_type")
      t1 = t1[order(t1$value),]
      signature=if(prev_signature=="") t1$names else paste(prev_signature, t1$names, sep=";")
      
     t2 =   t1 |> tibble::add_column(sumAngle, signature);
     subset(t2, value<999)
    })
  })
  comb_all2
}

## this is a class which holds a data object and interacts with the coordination node
.getAllSparseMatrices<-function(data, hasNA=T, convertToBigMatrix=F, min_variance =0.001, max_na_proportion=0.99){
  rn = unlist(lapply(data, function(d1) rownames(d1)))
  rn = rn[!duplicated(rn)]
  lapply(data, function(mat){
   m1 =  .getSparseMatrices(mat, hasNA=hasNA, convertToBigMatrix = convertToBigMatrix,rn = rn)
   
     sv = sparse_variance(m1$matrix)
     na_cnt = colSums(m1$matrixNA)/nrow(m1$matrixNA)
     m2 = lapply(m1, function(m11) m11[,sv>min_variance & na_cnt<=max_na_proportion])
  m2
  })
  
}
##this function removes NAs
## if no NA matrixNA is just empty matrix
.getSparseMatrices<-function(mat, hasNA=T, convertToBigMatrix=FALSE,rn = rownames(mat)){
  mi1 =  match(rownames(mat), rn)
  mi0 =  match(rn,rownames(mat))
  newNA=T
  if(length(mi1)==length(rn)){
    if(max(abs(apply(cbind(mi1, 1:length(rn)),1,diff)))==0) newNA=F
  }
  newNA = length(which(is.na(mi0)))>0
  
  if(!hasNA& !newNA){
    if(convertToBigMatrix){
      stop(' not supported in this version')
#      m2=matrix(0, nrow = nrow(mat), ncol = ncol(mat))
#      res1 = list(matrix = as.big.matrix(mat),  matrixNA = as.big.matrix(m2)
                
      #)
    }else{
      mat1 = if(typeof(mat)=="S4") mat else Matrix(mat);
      res1 = list(matrix = mat1,
                  matrixNA = Matrix(0,nrow(mat) , ncol(mat), sparse = T)
      )
    }
  }else{
    if(newNA){
      m1=apply(mat,2,function(v){
        v1 = v[mi0]
        mv = mean(v, na.rm=T)
        if(is.na(mv)) mv =0
        v1[is.na(v1)]=mv
        v1
      })
      rownames(m1) = rn
      m2 = apply(mat,2,function(v){
        v1 = rep(1, length(rn))
        v1[mi1[!is.na(v)]]=0
        v1
      })
      rownames(m2) = rn
      
    }else{
      m1=apply(mat,2,function(v){
        mv = mean(v, na.rm=T)
        if(is.na(mv)) mv =0
        v[is.na(v)]=mv
        v
      })
      m2 = apply(mat,2,function(v){
        v1 = rep(0, length(v))
        v1[is.na(v)]=1
        v1
      })
    }
    if(convertToBigMatrix){
      stop("not supported")
#      res1 = list(matrix= as.big.matrix(m1), matrixNA = as.big.matrix(m2))
    }else{
      res1 = list( matrix = Matrix(m1, sparse=T),matrixNA = Matrix(m2, sparse=TRUE)) 
    }
  }
  res1
}

getFullModels<-function(all_models){
  .nonZero(lapply(all_models, function(all_model1){
    .nonZero(lapply(all_model1, function(all_model2){
      .nonZero(lapply(all_model2, function(all_model3){
        all_model3[names(all_model3) %in% "full"]
      }))
    }))
  }))
}



#' Infer the statistical families of the phenotypes
#' @param y_mat a matrix of phenotypes, with columns as variables and rows as samples
#' @param max_ordinal  maximum number of values before we consider an ordinal value as a continuous value
#' @returns a list with the statistical families detected
#' @export
getFamily<-function(y_mat, max_ordinal=getOption("max_ordinal",10)){
  types = attr(y_mat, "types")
  
  if(!is.null(types)){
    indst = 1:length(types)
    names(indst) = names(y_mat)
    family = (lapply(indst, function(i){
      typ = types[[i]]
      if(typ=="double") return("gaussian")
      if(typ=="boolean") return("binomial")
      if(typ=="integer") {
        tbl=if(is.list(y_mat)) table(y_mat[[i]]) else table(y_mat[,i])
        if(length(tbl)<=2) return("binomial")
        if(length(tbl)>max_ordinal) return("gaussian")
        return("ordinal")
      }
      if(typ=="character"){
        tbl=if(is.list(y_mat)) table(y_mat[[i]]) else table(y_mat[,i])
        if(length(tbl)<=2) return("binomial")
        return("multinomial")
      }
    }))
    subinds = unlist(lapply(family, length))==0
    return(family)
  }
  if(typeof(y_mat)=="list"){
    nmey = names(y_mat); names(nmey)=nmey
    family = lapply(nmey,function(ynme){
      y  = y_mat[[ynme]]
      if(!is.numeric(y)){
        y = as.factor(y);
      }
      if(is.factor(y)){
        return(if(length(levels(y))<=2 )"binomial" else "multinomial")
      }
      if(length(unique(y[!is.na(y)]))<=2) return("binomial")
      vals = unique(y)
      if(sum(abs(vals-round(vals)), na.rm=T)<1e-9) {
        if(length(table(vals))>max_ordinal) return("gaussian")
        return("ordinal")
      }
      return("gaussian")
    })
    return(family)
  }else{
    
    famsy = (apply(y_mat,2,function(y){
      if(!is.numeric(y)){
        y = as.factor(y);
      }
      if(is.factor(y)){
        return(if(length(levels(y))<=2 )"binomial" else "multinomial")
      }
      if(length(unique(y[!is.na(y)]))<=2) return("binomial")
      vals = unique(y)
      if(sum(abs(vals-round(vals)), na.rm=T)<1e-9){
        if(length(vals)>max_ordinal) return("gaussian")
        return("ordinal")
      }
      return("gaussian")
    }))
    famsy
  }
}




#' Data holder class
#'
#' @description
#' A class that encapsulates a dataset holder
#'
#' @details
#' Available methods:
#' \itemize{
#'   \item \code{new()} - Create a new instance
#'   \item \code{name()} - Get dataset name
#'   \item \code{clone()} - Clone the object
#' }
#'
#' @export
dataH<-R6::R6Class("dataH", 
                   inherit = analysisBase,
 private = list(
   data ="environment",
  # sigsdir="character",
  # sigs="environment",
  #nme="character",
   type="character",  
   flags="list",
   data_id="character",
   #transform_y="character",
   #var_t = "list",
   
   dbDir="character",
   makeModels=function(vars2, inds, phens,flags){
     checkRMSV = .readFlag(flags,"checkRMSV",FALSE)
     d = private$data
     logpthresh= log(.readFlag(flags,"pthresh",1e-3))
     project=.readFlag(flags,"project",TRUE)
     useoffset=.readFlag(flags,"useoffset",TRUE)
     #  train_nme = .readFlag(flags,'train', names(datas)[1])
     #  if(length(which(train_nme %in% names(private$datas)))==0)train_nme = names(private$datas)[[1]]
     verbose=.readFlag(flags,"verbose",FALSE)
     if(!is.null(flags[['useglm']])) stop("define useglmnet not useglm")
     useglm=.readFlag(flags,"useglmnet",TRUE)
     inds1 = inds#[[nmes_inds1]]
     phens1 = phens#[[nmes_inds1]]
     #k=inds1[[1]]; d = datas[[1]]
     mods1 = lapply(inds1, function(k){
       # print(k)
       # lapply(datas[names(datas) %in% train_nme], function(d){
       mods = d$makeModels(phens, vars2,k,logpthresh = logpthresh,project=project,
                           flags=flags,checkRMSV=checkRMSV,
                           useglm=useglm, useoffset=useoffset)
       mods
       #})
       #})
     })
     if(length(mods1)==0) stop("!!")
     models=mods1
     #})
     vars = names(models[[1]])
     names(vars) = vars
     models2 = lapply(vars, function(v){
       # lapply(models, function(models1){
       m3 = lapply(models, function(m){
         m[[v]]
         #m2[unlist(lapply(m2, length))>0]
       })
       m3[unlist(lapply(m3, length))>0]
       #})
     })
     #  print(names(models2))
     models2
   },
   select_k=function(analysis,phens,flags, k1,
                     vars_l_todo ,
                     verbose=F){
     expt_id = super$expt_id();
     get_plots=.readFlag(flags, "get_plots",F) 
     stop_y = .readFlag(flags, 'stop_y',"rand")
     logpvthresh = log(.readFlag(flags,"pthresh",0.1))
     beam= .readFlag(flags,"beam",1)
     comb2 = NULL;
     plot_results = list()
     # vars_l = analysis$nextVars(expt_id, flags)
     nvar=0;
     while(length(vars_l_todo$todo1)>0 ){
       comb2_new = try(self$multiAnglesAndPv(comb2, phens, k1,flags,expt_id, vars_l_todo))
       if(inherits(comb2_new,"try-error")) break;
       if(F && flags$plot){
         ggps=plot_ri(comb2_new,1)
         plot_grid(ggps[[1]], ggps[[2]], ggps[[3]])
       }
       data_nme=self$name();
       vars_l_todo_new=analysis$savePvalsAndNextVars(flags,phens,vars_l_todo,comb2,data_nme,  k1)
       if(length(vars_l_todo$todo1)==length(vars_l_todo_new$todo1)){ ## to account for continuing via shortening todo rather than adding variable
         comb2 = lapply(comb2_new, function(x) x$pvs)
         
       }
       vars_l_todo = vars_l_todo_new
       nvar = length(vars_l_todo$vars_l[[1]]$var)
       if(verbose) print(names(vars_l_todo$vars_l))
       if(length(vars_l_todo$vars_l[[1]]$var_names)>=flags$max) break;
       if(get_plots && length(vars_l_todo$todo1)>0){
         
         plot_results[[nvar]] = mergeAll(comb2_new, flags$beam)
       }
     }
     attr(vars_l_todo,"plots")=get_comb_plot(plot_results, nvar)
     
     vars_l_todo
   },
   var_thresh = function(qq_t){
     lapply(private$data$vars, function(v) quantile(v, qq_t))
   },
 
   updateTrain=function( phens, flags, transform_y,  verbose=FALSE, force=F){
     private$data$updateTrain( phens,flags,transform_y, verbose=verbose, force=force)
   },
   updateLOOC=function( phens, flags,varn=c(),verbose=FALSE, force=F){
     private$data$updateLOOC( phens,flags,varn=varn,force=force, verbose=verbose); ### update training object - updates all
   },
   findPrev=function(comb20, expt_id, prev_i, k){
     if(is.null(private$sigs)){
       nmes= unlist(lapply(prev_i$var_names, function(x) paste(x,collapse=".")))
       if(length(nmes)==0)nmes="empty"
       nme1 = strsplit(nmes[length(nmes)],"\\.")[[1]]
       str = paste(nmes[-length(nmes)], collapse=";")
       if(str=="") str="empty"
       prev_i2= comb20[[str]][[nme1[3]]][[nme1[4]]][[nme1[2]]]
     }else{
       
       prev_i2 =   private$sigs$loadPrev(expt_id, prev_i, k, data_nme = private$nme)
     }
     if(is.null(prev_i2)){ 
     #  warning("could not find")
      # print(prev_i)
       
       prev_i2 = prev_i
     }
     return(prev_i2);
   },
   res_inner=function(comb_,prev_i, flags,k, expt_id, phens){
     
     nme_comb = names(comb_); names(nme_comb) = nme_comb
     #nme_c1 = nme_comb[[1]]; nme_p1 = names(comb_[[nme_c1]])[[1]]; ik=1
     res_inner=lapply(nme_comb, function(nme_c1){
       nmesp1 = names(comb_[[nme_c1]]); names(nmesp1) = nmesp1
       lapply(nmesp1, function(nme_p1){
         comb = comb_[[nme_c1]][[nme_p1]]
         if(nrow(comb)==0) return(NULL)
         num_pvals1 = nrow(comb)
         inds1p = 1:num_pvals1; names(inds1p) = comb$names[1:length(inds1p)]
         nxt_vars = lapply(inds1p, function(ik){
           b_i_name = c(comb$data_type[[ik]], comb$names[[ik]], nme_c1,nme_p1)
           angle=comb$value[[ik]]
           if(!is.null(flags$angles_only) && flags$angles_only){
             b_i = private$data$convert(b_i_name)
             nv = list(angle = angle, var = c(prev_i$var, list(b_i)), 
                       angles = c(prev_i$angles, angle),
                       var_names = c(prev_i$var_names, list(b_i_name)),
                       varnames = c(prev_i$varnames, paste(b_i_name, collapse="."))
             )
            
             #  print(nv)
             nv$sumAngle = sum(nv$angles);
             
           }else{
             nv= private$getPvsAll(phens,prev_i, b_i_name,k,  prev_i$Wall,flags, angle=angle)
           }
           names(nv$var_names) = nv$varnames;
           if(inherits(nv,"try-error")) {
             print(paste(nme_c1, "error"))
             return(NULL)
           }
           nv
         })
         
       })
     })
     res_inner   
   },
   post_process_to_go=function(variables, flags, phens, useDB=FALSE){
     useAngles = !is.null(flags$angles_only) && flags$angles_only
     full_index = length(variables) 
     beams = 1:length(variables[[full_index]])
     names(beams)=beams
     names(variables) = 1:length(variables)
     
     vars_combined=lapply(beams, function(beam){
       print(beam)
       vars_all = list()
       vars_all1 = list()
       vars_all2 = list()
       #vars_all3 = list() #funcstr
       # func_inds = lapply(variables, function(vv) attr(vv,"func_ind"))
       
       for(repn in names(variables)){
         full = repn==full_index
         if(beam>length(variables[[repn]])) {
           print("skipping")
           next;
         }
         var1 = variables[[repn]][[beam]]   ### only taking the top1
         var2 = var1$var_names
         if(length(var2)>0){
           names(var2) = names(var1$var_names)
           cumpv = if(useAngles) var1$cum_angle else  var1$cum_pv#lapply(var1, function(vv) attr(vv,"cumpv"))
           varn = paste(names(var2),collapse=";") #paste(names(var2), collapse=";")
           if(is.null(vars_all[[varn]])){
             vars_all[[varn]] = list()
             vars_all1[[varn]] = var2
             vars_all2[[varn]] = list()
             # vars_all3[[varn]] =list(repn) 
             # names(vars_all3[[varn]]) = func_str1
           }else{
             #vars_all3[[varn]]=c( vars_all3[[varn]],repn)
           }
           repn1 = as.list(as.numeric(repn))
           repn2 = as.list(cumpv)
           names(repn1) = if(full) "full" else repn
           names(repn2) = if(full) "full" else repn
           if(is.null(vars_all[[varn]])){
             vars_all[[varn]] =repn1
             vars_all2[[varn]] =repn2
             # vars_all3[[varn]]
           }else{
             vars_all[[varn]] = c(vars_all[[varn]] , repn1)
             vars_all2[[varn]] = c(vars_all2[[varn]] , repn2)
           }
         }
       }
       list(variables = vars_all1, inds = vars_all,cumpv=vars_all2, beam=beam,flags = flags,phens = phens )# ,transf= vars_all3) 
     })
     if(useDB  && !is.null(private$sigs)){
       private$sigs$saveVars(vars_combined,replace=T)   #saving local
     }
     
     vars_combined
   },
   
   
   simplify = function(ri){
     ri_out=lapply(ri, function(ri1){
       lapply(ri1, function(ri2){
         lapply(ri2, function(ri3){
           ri3$simplify()
           
           #list(pvs = ri3$pvs, pvs_all = ri3$pvs_all, var_names = ri3$var_names)
           return(ri3)
         })
       })
     })
     ri_out
   },
   
   combinedAngles=function(phens, varnames, incl, k, g_incl, qq_t, flags, sumAngle, addPlot=F){ #phens, varnames, incl=incl, k=k, type=type
      type=private$type
     prev_signature =paste(unlist(lapply(varnames, function(vn)vn[2])),collapse=";")
     var_t = private$var_thresh(qq_t)
     
     angles=  private$data$getAngles1(phens,varnames,incl=incl,k=k, type=type)
     angles =angles[ unlist(lapply(angles, length))>0]
     angleH=list(angles=angles,
                 cols_incl = private$data$cols_incl(var_t,incl, g_incl,excl=varnames)) ### fix 
     comb_angle1 =  .combineAngles1(angleH, incl, flags, sumAngle, prev_signature, excl=varnames)
     if(addPlot){
       all_angles = .extrAngles(angleH,comb_angle1, incl)
        attr(comb_angle1,"all")=all_angles;
     }
     comb_angle1
   },
   getPvsAll=function(subphens, prev_i, b_i_name,k, #   prev_i = vars_l1[[nmed]]
                      Wall, # =lapply(subphens, function(f) matrix(nrow=0,ncol=0)),
                      flags, angle=0){
     #useglm=FALSE,,inv_transform=getOption("x_transform",T),
     #project=T, useoffset=T){
     inv_transform=T
     project=.readFlag(flags,"project",T)
     useoffset=.readFlag(flags,"useoffset",T)
     useglm = .readFlag(flags,'useglmnet',T)
     d = private$data
     family = strsplit(names(subphens)[[1]],"\\.")[[1]][1]
     if(family=="multinomial") useoffset=F
     
     prev_i1 = d$makeNextModel(prev_i,b_i_name,subphens,k, family, ypred=NULL, 
                               project=project, useglm=useglm, logpthresh =0, useoffset=useoffset)
     prev_i1$angle= angle
     prev_i1$angles = c(prev_i$angle, angle)
     prev_i1$sumAngle = sum(prev_i1$angles);
     prev_i1
     #   pvs_all = pvs_all[unlist(lapply(pvs_all, length))>0]
     # pvs_all
   }
 ),
 
 public = list(
   #' @description Create a new instance
   #' @param data a list of matrices (of different modalities), with columns as variables and rows as samples
   #' @param y phenotype matrix, with columns as outcomes and rows as samples
   #' @param nme The name of the data objetct. 
   #' @param flags list of options, described in the vignette
   #' @param transform_y a transformation object from the function getYTransformation
   #' @param family the statistical family of phenotype y, can be calculated by getFamily
   #' @param dbDir dir for database to store results to speed up re-reruns.  Can be NULL, if not required
     initialize=function(
    data,
    y,
    weights = rep(1, nrow(y)),
    nme,
    flags ,
    transform_y=getYTransform(pow = 1,  n_random=1, perm=F),
    family= getFamily(y),
      dbDir=tempdir()
   ){
       
    super$initialize(nme,dims =lapply(data, dim),  flags=flags, dbDir=dbDir);
       colnames(y) = gsub("\\.","_",colnames(y))  ## no . allowed
      for(k in 1:length(data)){
        colnames(data[[k]]) = gsub("\\.","_",colnames(data[[k]]))  ## no . allowed
        
      }
       .check_data(data, y)  ## checks rownames match
    memDir=NULL
    useDB=!is.null(dbDir);
    convertToBigMatrix=F #.readFlag(flags,"covertToBigMatrix", F)
    hasNA=.readFlag(flags,"hasNA", T)
    preprocessed=.readFlag(flags,"preprocessed", F)
    max_na_proportion =.readFlag(flags,"max_na_proportion",0.99)
    min_variance =.readFlag(flags,"max_na_proportion",0.001)
    
    mat = .getAllSparseMatrices(data,hasNA=hasNA, convertToBigMatrix=convertToBigMatrix,min_variance = min_variance, max_na_proportion=max_na_proportion)
    #print("HHHHH")
    #print(mat)
    private$dbDir = dbDir
    #nme=sub("/",".",nme)
  #  private$flags = flags
   # transform_y =.readFlag(flags, "transform_y",toJSON(list(x=list(unvfunc="function(y,param) y",func="function(y,param) y", param=1))))
    ### MAKE SIGNATURE DIRECTORY
   # private$nme=nme
  #  private$sigsdir=paste(dbDir,paste0("fspls_signatures__",private$nme,sep="/"))
    #private$sigs = list()
    #####
    private$data = 
      dataObj$new(mat, private$nme,dbDir,flags, 
                  incl_full=T,seed = getOption("seed",42), memDir=if(is.null(memDir)) NULL else paste(memDir, private$nme,sep="/"))
    private$type="slow1"
    types_all = getOption("types_all",names(private$data$data))
    names(types_all) = types_all
    batch=.readFlag(flags, "batchsize",0)
    all_v_all = .readFlag(flags,"all_v_all",F)
    one_v_rest = .readFlag(flags,"one_v_rest",F)
    nrep = .readFlag(flags,"nfold",if(batch>0) 0 else 1)
    if(nrep>0 && batch>0)warning("only one of nfold or batchsize should be non zero")
    pheno_balance=.readFlag(flags,"pheno_balance",NULL)
    varn = getOption("varn",c())
    if(is.null(weights)) stop("data weights need to be defined");
    
    #invisible(lapply(1:length(datas), function(ik) {
     # family = families[[ik]]
        private$data$updateY(y,weights, preprocessed=preprocessed, family=family, CHECK=T, all_v_all=all_v_all, one_v_rest = one_v_rest)
    
    #  private$sigsdir=paste(dbDir,"fspls_signatures1",sep="/")
      #dir.create(private$sigsdir, recursive=FALSE, showWarnings=F)
     # dims1 = list(private$data$dims()); names(dims1) = private$nme
      #private$sigs=   if(exists("dbConnect") && useDB) sigEnv$new(private$sigsdir,nme,flags, dims1, clear=F) else NULL;
      
    
  },
  #' @description Calculate the angles and pv across multiple phenotypes.  This is an internal function and should not need to be called by user
  #' @param comb20 values from previous iteration
  #' @param phens1 phenotypes being used
  #' @param k1 k1 is the current fold
  #' @param flags list of options
  #' @param expt_id and ID for the experiment (can be 0)
  #' @param vars_l_todo  the variables under consideration
    multiAnglesAndPv=function(comb20, phens1,  k1,flags, expt_id, vars_l_todo){
      get_plots=.readFlag(flags, "get_plots",F) 
      verbose=.readFlag(flags,'verbose',F)
    saveAngles=FALSE
    if(is.null(expt_id)) stop("expt_id is NULL")
    vars_l = vars_l_todo$vars_l
    todo1=vars_l_todo$todo1[[1]]
    incl=todo1$incl
    g_incl = todo1$g_incl
    qq_t = todo1$qq
    comb2_new=invisible( lapply(vars_l, function(prev_i){
      prev_i2 = private$findPrev(comb20, expt_id, prev_i, k1);
      varnames = prev_i2$var_names; 
      sumAngle =sum(prev_i2$angles)
      comb_=private$combinedAngles(phens1, varnames, incl, k1,  g_incl, qq_t, flags, sumAngle,addPlot=get_plots) ;
      if(length(comb_)==0) stop("length zero")
      # return(list(comb_angle1, all_angles));
      
      if(saveAngles) return(comb_)
      #comb_ = private$anglesAndPv(phens, prev_i, incl, k1, g_incl, qq_t, flags,expt_id, saveAngles=saveAngles, verbose=verbose)
      ri = private$res_inner( comb_,prev_i2,flags,k1, expt_id, phens1)
     super$savePvals(flags,phens1, ri, varnames,k1,useCurrVarnames=T)
      
      
      list(angles = comb_, pvs = ri) ;#private$simplify(ri))
    }))
    comb2_new
  },
  #' clear database. This is only used if you want to clear the database.  Internal function
  #'
  #' @param drop drop the tables instead of clear
  #' @param exclude which tables not to clear
  clear_db=function(drop=FALSE, exclude="vars"){
    if(drop){
      if(!is.null(private$sigs)) private$sigs$drop_all(exclude=exclude)
    }else{
      warning("need to set drop=T if you are sure, this will delete all saved signatures")
    }
    
  },

  #' get the data types
  #' @returns which data types are in the dataset
  data_types=function(){
    names(private$data$data)
  },
  #' @description split dataset into smaller datasets
  #' @param proportions what proportions to split into
  #' @returns a list of dataH objects with data partiioned according to proportions
 split=function(proportions = c(0.5,0.5)){
   datas = private$data$split(proportions);
   mats = datas$mats;
   ys = datas$ys;
   nme_d = names(mats); names(nme_d) = nme_d
   all_v_all = .readFlag(private$flags,"all_v_all",F)
   one_v_rest = .readFlag(private$flags,"one_v_rest",F)
   dbDir =private$dbDir;
   #nme = nme_d[[1]]
   lapply(nme_d, function(nme){
     mat = mats[[nme]]
     dh1 = dataH$new(NULL,nme=nme, y=NULL, y1=ys[[nme]], family = private$data$family, mat = mat,       dbDir = dbDir, flags=flags, useDB=!is.null(private$sigs))
         dh1
     
   })
 },
 
 #' @description get the phenotyeps
 #' @returns  phenotypes
  pheno=function(){
    maxpheno=1e9;sep=F; sep_group = F;exclude=NULL; code=NULL; memb=NULL
   res = private$data$pheno(maxpheno=maxpheno, sep=sep, sep_group = sep_group, code = code,memb=memb);
   if(!is.null(exclude)){
     lapply(res, function(res2){
       lapply(res2, function(res1) res1[-grep(exclude,res1)])
     })
   }
   
   res
 },
 #' @description get the dimensions
 #' @returns  dimensions
 dims=function(){
   private$data$dims()
 },
 #' get the nreps
 #' @returns  nreps
 nreps=function(){
   ncol(private$data$looc$incl)
 },
 
 #' @description update the phenotypes without remaking the entire object. You can provide many phenotypes in the initialisation stage, but only consider a subset in model fitting stage in this way.
 #' @param phens phenotypes
 #' @param flags list of options
 #' @param transform_y transformation object 
 
 update=function(phens, flags, transform_y=fromJSON(flags$transform_y)){
    verbose=.readFlag(flags,'verbose',F)
    force=.readFlag(flags,'force',F);
    flags$transform_y = transform_y;
    super$updateExpt(phens, flags)
#   flags[['data_types']] =toJSON(names(datasH[[1]]$data$data))

   if(is.null(flags[['data_types']]) || flags[['data_types']]=="{}")flags[['data_types']]=toJSON(names(private$data$data))
     private$updateLOOC(phens, flags, verbose=verbose,force=force)
     private$updateTrain(phens, flags,transform_y, verbose=verbose, force=force)
   ##updated after updateLOOC
   nreps1 =self$nreps()
   nreps = 1:nreps1
   names(nreps) = nreps
   nreps
 },
#' @description main function for variable selection
#' @param analysis an analysisEnv object
#' @param phens list of phenotyps
#' @param flags list of options
#' @param transform_y transformation object 
#' @param useDB boolean to indicate if results should be saved to database
 select=function(analysis, phens,flags,transform_y, useDB=F
                 ## expt_id specific to this database .. might be diff for global
               ){#c(y="function(y) y","function(y) y")
   verbose=.readFlag(flags,'verbose',F);
   super$updateExpt(phens, flags);
   flags$transform_y = toJSON(transform_y);
   if(is.null(flags[['data_types']]) || flags[['data_types']]=="{}")flags[['data_types']]=toJSON(self$data_types())
      if(flags$topn<flags$beam) stop("beam should be less than topn")
   if(is.null(flags[['data_types']])) flags[['data_types']] = names(private$data$data)
   nreps = self$update(phens, flags, transform_y);
   
     vars_all = private$loadVars(phens)
     if(!is.null(vars_all)) return(vars_all)
   
   vars_l_todo = analysis$getTodo(flags, phens)
  # expt_id=super$getExpt(flags, phens, add_new=T)
   variables=lapply(nreps, function(k1){
     if(verbose) print(paste("cv",k1,"of",length(nreps)))
    private$select_k(analysis, phens,flags, k1, vars_l_todo,verbose=verbose)
               })
   variables1 = lapply(variables, function(vars_l_todo){
     vars_l_todo$vars_l 
   })
  
   get_plots=.readFlag(flags, "get_plots",F) 
#   private$sigs$clearPvals(expt_id)
   vars_all=post_process(variables1,flags,phens)
   private$saveVars(vars_all);
   
   if(get_plots){
   attr(vars_all,"plots") = .merge1_new(lapply(variables, function(vars_l_todo){
     attr(vars_l_todo,"plots") 
   }), addName="nrep")
   }
  vars_all
  
 },
 #merge22=function(comb2, vars_l){
#  lapply(vars_l, function(vars_l1){
#      vars_l1$var_names
#  })
# },
 
#' @description provides access to internal storage of phenotype data
#' @returns list of matrices
y=function(){
  private$data$y
},

#' @description extract the predictions for the fitted models
#' @param all_modelsh fitted models from makeAllModels
#' @param phens list of phenotyps
#' @param flags list of options
#' @param CV cross validation predictions
#' @param liab return liability score, or probability (for binomial, ordinal multinomial)
#' @returns  a table with results
extractPredictions=function(all_modelsh,phens=all_modelsh$phens, flags=all_modelsh$flags){
  private$updateLOOC(phens, flags)
  
  all_models_y0 = all_modelsh$models#[[mod_nme]]
  # eval1 =  .merge1_new(lapply(nme_d2, function(nme1){
  #print(nme1)
  d = private$data
  #all_models_y = all_models_y0[[1]]
  predictions0 = lapply(all_models_y0, function(all_models_y){
    d$extractPredictions(all_models_y, phens, flags)
  ##  d$evaluateAllModels(all_models_y,phens,flags, verbose=verbose) |> tibble::add_column(data=private$nme, trainedOn=all_modelsh$trainedOn)#|> tibble::add_column(trainedOn=private$nam)
  })
 
#  predictions0=res3[unlist(lapply(res3, function(x) length(x[[1]])))>0]
  predictions0 # 
},

#' @description plot all data for selected variables
#' @param vars_all variables
#' @param phens list of phenotyps
#' @param all_types use all types?
#' @param transform_x which transform to use on x
#' @param violin violin plots
#' @param assoc use association
#' @returns  a table with results
plotData=function(vars_all, phens = vars_all$phens, all_types=FALSE, transform_x = NULL, violin=FALSE, assoc=F){
  vars_all = vars_all1
  df4= #.merge1_new( 
   # lapply(private$datas, function(d) 
      private$data$plotData(vars_all1, phens1 = phens, all_types=all_types, transform_x = transform_x, violin=violin, assoc=assoc)
               #     addName="dataset")
  df4 = df4|>tibble::add_column(dataset=private$nme);
  facet= if(!is.null(df4[['transform']]) ) "transform~pheno1" else "pheno1"
  df4$y = factor(df4$y)
  gene_levs = unlist(lapply(unlist(vars_all1$variables,recursive=FALSE), function(x) x[[2]]))
  gene_levs = gene_levs[!duplicated(gene_levs)]
  df4$gene = factor(df4$gene, levels = gene_levs)
  df5=df4|> separate(col="pheno",sep="\\.", into=c("family","pheno1"))
  df6=subset(df5, family=="gaussian")
  df7=(subset(df5, family!="gaussian"))
  # color= if(length(unique(df4$dataset))>1) "dataset" else #"pheno"
  ggp = NULL; ggp1 = NULL
  if(nrow(df7)>0){
    ggp<-ggplot(df7, aes(x=y, y=value, color=gene, shape=data, linetype=dataset))+facet_wrap(facet, scales="free");#+ggtitle(unlist(phens1))
    if(violin){
      ggp<-ggp+geom_violin()+geom_point()
    }else{
      ggp<-ggp+geom_boxplot()
    }
  }
  if(nrow(df6)>0){
    prbs=c(0.33,0.5,0.66)
    df6_1 =df6 |> unite("comb",gene,data,family,pheno1,dataset,sep="__")
    comb1=unique(df6_1$comb); names(comb1)=comb1
    quants=.merge1_new(lapply(comb1,function(c1){
      df6_2 = subset(df6_1, comb==c1)
      df6_2$y = as.numeric(as.character(df6_2$y))
      df6_2 = df6_2[order(df6_2$y),]
      yv = unique(df6_2$y)
      names(yv)=yv
      .merge1_new(lapply(yv, function(yv1){
        df6_3 = subset(df6_2, y<=yv1)
        q1=quantile(df6_3$value,na.rm=T,probs = prbs)
        df_=data.frame(t(data.frame(q1)))
        df_
      }),addName="y")
    }), addName="comb")
    quants1 = quants |> separate("comb", sep="__", into=c("gene","data","family","pheno1","dataset"))
    names(quants1) = gsub("\\.","",names(quants1))
    quants1$y = as.numeric(quants1$y)
    ggp1<-ggplot(quants1, aes(x=y, y=X50, ymin=X33, ymax = X66,color=gene, fill=gene,shape=data, linetype=dataset))+facet_wrap(facet, scales="free");#+ggtitle(unlist(phens1))
    ggp1<-ggp1+geom_line()+geom_ribbon(alpha=0.1)
  }
  
  list("binomial"=ggp,"gaussian"=ggp1)
},

#updateWeights=function(subphens = self$pheno()[[1]][1]){ ## upweights low count values
#  for(k in 1:length(private$datas)){
#    private$datas[[k]]$updateWeights(subphens)
#  } 
#},

#' @description get data after projection
#' @param varnames varnames
#' @returns projected data after projecting out varnames
getProjectedData=function(varnames){
#  lapply(private$datas, function(d){
  d = private$data
    d$getProjectedData(varnames = varnames);      
 # })
},

#' @description get variance of data
#' @param varnames varnames
#' @returns variance
getVariance=function(varnames){
  d = private$data
    d$getVariance();      
},

#' @description fit models based on variables
#' @param vars_all list of variables selected by select method
#' @param phens list of phenotypes
#' @param flags flags
#' @param transform_y an object to determine the y transformaton
#' @param useDB whether to use the DB to save the reuslt
#' @returns fitted models
makeAllModels=function(vars_all, 
                       phens=vars_all[[1]]$phens, flags=vars_all[[1]]$flags, 
                       transform_y =  fromJSON(flags$transform_y),
                       useDB=F){
  
  self$update(phens, flags,transform_y);
  verbose=.readFlag(flags,'verbose',F); max = .readFlag(flags,'max',1e6)
  sigDB = if(useDB) private$sigs else NULL
  if(!is.null(sigDB) ){
    all_models =try( sigDB$loadModels(flags,phens))
    if(inherits(all_models,"try-error")) {
      print(paste("problem reading from DB .. recalculating"))
    }else if(!is.null(all_models) && length(all_models$models)>0){
      return(all_models)
    }
  }
  private$updateLOOC(phens,flags)
  logpthresh= log(.readFlag(flags,"pthresh",1e-3))
  project=.readFlag(flags,"project",TRUE)
  beams = names(vars_all); names(beams)=beams
  all_models_full=lapply(beams, function(beam){
    vars_all0 = vars_all[[beam]]
  vars = vars_all0#[[nme_v_all]]
  all_models = list()
  variables = vars$variables
  var_inds = vars$inds
  rem_inds = 1:self$nreps()
  names(rem_inds) = as.character(rem_inds)
  names(rem_inds)[which(rem_inds==self$nreps())]="full"
  all_models = private$makeModels(list(),rem_inds , phens, flags)

  if(length(variables)==0)  return(all_models) ;#return(list(models=all_models, flags = flags, phens = phens, trainedOn=private$nme))
  ord = order(unlist(lapply(variables, length)),decreasing=T)
  variables = variables[ord]
  var_inds = var_inds[ord]
  #v_nme = names(vars_all$variables)[1]; max=10; verbose=T; k=1;variables =vars_all$variables; 
  
  for(v_nme in names(variables)){
    if(verbose)print(v_nme)
    vars2 = variables[[v_nme]]
    vars2 = vars2[1:min(length(vars2), max)]
    inds =var_inds[[v_nme]]
    nme_ = paste(names(vars2),collapse=";")
    models1 = all_models[[nme_]]
    if(is.null(models1)){
      models1 = private$makeModels( vars2, inds,phens,flags)
      for(k in 1:length(models1)){
        mod1 =   all_models[[names(models1)[[k]]]]
        if(is.null(mod1)){
          all_models[[names(models1)[[k]]]] = models1[[k]]
        }else{
          # for(p_nme in names(models1[[k]])){
          mod2 = all_models[[names(models1)[[k]]]]#[[p_nme]]
          if(is.null(mod2)){
            all_models[[names(models1)[[k]]]]= models1[[k]]#[[p_nme]]#[[p_nme]] 
            
          }else{
            for(r_nme in names(models1[[k]])){ #[[p_nme]])){
              mod3 = all_models[[names(models1)[[k]]]][[r_nme]] #[[p_nme]][[r_nme]]
              if(is.null(mod3)){
                all_models[[names(models1)[[k]]]][[r_nme]]= models1[[k]][[r_nme]] #[[p_nme]][[r_nme]] #[[p_nme]][[r_nme]] 
                
              }
            }
          }
          #}
        }
      }
    }else{ ## fill in gaps
      #for(p_i in 1:length(inds)){
      # p_nme = names(inds)[p_i]
      models2 = all_models[[nme_]]#[[p_nme]]  
      if(is.null(models2)){
        models2 = private$makeModels( vars2, inds#[p_i]
                                   ,phens#[which(names(phens) %in% p_nme)]
                                   ,flags)
        for(k in 1:length(models2)){
          all_models[[names(models2)[[k]]]] = models2[[k]]#[[p_nme]]#[[p_nme]]
          
        }
      }else{
        subinds= inds[which(is.na(match(names(inds), names(all_models[[nme_]]))))]
        if(length(subinds)>0){ ## missing inds
          models3 = private$makeModels( vars2, subinds,phens#[which(names(phens) %in% p_nme)]
                                     ,flags)
          for(nme1_ in names(models3)){
            for(r_i in names(models3[[nme1_]])){
              all_models[[nme1_]][[r_i]] = models3[[nme1_]][[r_i]]
            }
          }
        }
      }
      #}
      
    }
  }
  all_models

  })
  all_models_=list(models=all_models_full, flags = flags, phens = phens, trainedOn=private$nme)
  if(useDB  && !is.null(sigDB)){
    sigDB$saveModels (all_models_)
  }
  #combined_models
  all_models_
},
#' @description evaluate the fit of models
#' @param all_modelsh  models fitted from makeAllModels
#' @param phens list of phenotypes
#' @param flags flags
#' @param useDB whether to save results to DB
#' @returns evaluation of models
evaluateAllModels=function(all_modelsh, phens=all_modelsh$phens,flags=all_modelsh$flags,useDB=F){ ## different folds with same variables
  sigDB = if(useDB) private$sigs else NULL
  if(!is.null(sigDB) ){
    eval1 = sigDB$loadEval(flags,phens,)
    if(!is.null(eval1) && nrow(eval1)>0){
      return(eval1)
    }
  }
  verbose=.readFlag(flags,"verbose",T)
  inv_transform_y=F
  private$updateLOOC(phens, flags)
  if(length(all_modelsh$models)==0) return(NULL)
  #nme_d2 = .readFlag(flags,"test",names(private$datas))
  #names(nme_d2) = nme_d2
  all_models_y0 = all_modelsh$models#[[mod_nme]]
 # eval1 =  .merge1_new(lapply(nme_d2, function(nme1){
    #print(nme1)
    d = private$data
  #  all_models_y = all_models_y0[[1]]
    eval1 =   .merge1_new(lapply(all_models_y0, function(all_models_y){
   d$evaluateAllModels(all_models_y,phens,flags, verbose=verbose) |> tibble::add_column(data=private$nme, trainedOn=all_modelsh$trainedOn)#|> tibble::add_column(trainedOn=private$nam)
  }), addName="beam")  #if(inherits(resd,"try-error")) {
    #  print(resd)
    #  print(paste("problem", nme1))
    #  stop("!!")
    #  return(NULL)
    #}
#    resd
#  }),addName="data")
  #}),addName="transform_y")
  if(is.null(eval1)) return(NULL)
  #  eval1 = subset(eval1, model!="avg")
  eval2 = eval1|> pivot_wider(names_from="submeasure") #|> tibble::add_column(transform_y=strsplit(transform_y[[1]]," ")[[1]][2])
  
  #  isfull=eval2$model %in% full_model_nmes
  #  eval2|>tibble::add_column(isfull=isfull)
  if(!is.null(sigDB) ){
    sigDB$saveEval(eval2, flags,phens)
    eval1 = sigDB$loadEval(flags,phens)
    return(eval1)
  }
  #print("HH")
  #  if(T) return(eval2)
  eval3 = .calcEval1(eval2, rename=F)
  eval4 = eval3 |> tibble::add_column(variable= unlist(lapply(eval3$model, function(x){
    if(x=="cv" || x=="") return("");
    x1 =rev(strsplit(x,";")[[1]])[1]
    
    x2 = strsplit(x1,"\\.")[[1]]
    x2[min(length(x2), 2)]
  })))
  eval4$variable[is.na(eval4$variable)]=""
  eval4
}
    
)
)
  