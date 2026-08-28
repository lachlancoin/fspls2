#private = self[[".__enclos_env__"]]$private

mergeAll = function(comb2_new, beam){
  comb1 = .merge1_new(lapply(comb2_new, function(c){
    .merge1_new(lapply(c$angles, function(c1){
    
      .merge1_new(c1,addName="pow")
    }),addName="func")
  }), addName="prev") 
  if(is.null(comb1)) return(NULL)
   #comb1$prev 
   o =  order(comb1$value)
   comb1[o[1:beam],]
  
}
permuteLabel<-function(y,certainty,na_inds,
                       levs = levels(y[,1]),
                       offset = 2:length(levs)-1
                       ){
  y_orig = y[na_inds,,drop=FALSE]
  names(certainty) = rownames(y)
  alt_inds=lapply(offset, function(i){
    nrow(y)+(i-1) * length(na_inds) +  1:length(na_inds)
  })
  
  
  weights_new = rep(1, nrow(y))
  names(weights_new) = rownames(y)
  y_new = as.character(y[,1])
  y_new[na_inds] = levs[1]
  names(y_new)[1:nrow(y)] = rownames(y)
  
  for(i in 1:length(alt_inds)){
    ai = alt_inds[[i]]
    weights_new[ai] = rep(1, length(ai)) #(1-certainty[na_inds])/length(alt_inds)
    names(weights_new)[ai] = paste(rownames(y)[na_inds],"alt",i,sep=".")
   
    y_new[ai] = levs[i+1]
    names(y_new)[ai] = names(weights_new)[ai]
  }
  mi1 = match(y_orig[,1], levs)
  original_inds = rep(0, length(mi1));
  cert1 = certainty[na_inds]
  mat1 = cbind(mi1, cert1)
  mat2 = t(apply(mat1, 1,function(v){
    v1 = rep(0, length(levs))
    v1[v[1]] = v[2]
    v1[-v[1]] = (1-v[2])/(length(levs)-1)
    v1
  }))
  weights_new[na_inds] = mat2[,1]
  for(i in 1:length(alt_inds)){
    weights_new[alt_inds[[i]]] = mat2[,i+1] 
  }
  for(i in 1:length(mi1)){
    x = mi1[i];
    if(x==1){
      original_inds[i] =  (na_inds[i]) 
    }else{
      original_inds[i] =  alt_inds[[x-1]][i]
    }
  }
  
 #weights_new[original_inds]
 
  y_new2 = data.frame(list(y=factor(y_new, levels = levs)))
 
  list(y = y_new2,
#       y_orig = y_orig,
       original_inds = original_inds,
       weights = weights_new, offset = offset, alt_inds = alt_inds)
}

## expands to replace proportion with NA
expandData<-function(data, y, certainty,family,weights,  mult = 100 
){
  na_inds = which(certainty<1)
  non_na_inds = 1:nrow(y)
  non_na_inds = non_na_inds[-na_inds]
  if(length(na_inds)==0){
    return(list(dataset = data, y= y,  weights =weights, na_inds = na_inds, non_na_inds = non_na_inds, 
                original_inds = c(),
                alt_inds = c()))
  }
  .check_data(data, y)  ## checks rownames match
  if(max(certainty)>1.0) stop("!!")
  if(min(certainty<0.0)) stop("!!")
     y_alts = permuteLabel(y, certainty, na_inds)  
    y_new = y_alts$y
    weights = y_alts$weights * mult
  offset = y_alts$offset
  d_new =  lapply(data, function(d){
    d_out = d;
    for(i1 in offset){
      d1 = d[na_inds,,drop=FALSE];
      rownames(d1) = paste(rownames(d)[na_inds],"alt", i1, sep=".")
      d_out = rbind(d_out,d1)
    }
    d_out
  })
  if(family[[1]]=="binomial"){
    y_new[,1] = as.numeric(y_new[,1])
  }
  list(y = y_new, dataset=d_new, na_inds=na_inds,non_na_inds = non_na_inds, 
     
       original_inds = y_alts$original_inds,
      alt_inds = y_alts$alt_inds,
       weights =weights)
  
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

plot_angle_vs_pv<-function(comb2_new, alpha =.5, k1 = 1){
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
#my_colors <- grDevices::colorRampPalette(brewer.pal(12, "Paired"))(100)
#ggp1=ggplot(df1, aes(angle,pval, color=Column, shape=Row))+geom_point(alpha = alpha,size=2)+facet_grid("typ1 ~ni1")+ scale_shape_manual(values = rep(0:25, length.out = 100))+scale_color_manual(values = my_colors)
#ggp2=ggplot(df1, aes(angle,pval, color=Row, shape=Column))+geom_point(alpha = alpha,size=2)+facet_grid("typ1 ~ni1")+ scale_shape_manual(values = rep(0:25, length.out = 100))+ scale_color_manual(values = my_colors)



ggp3 = ggplot(df4, aes(angle,pval, color=typ1))+geom_text(alpha = alpha, aes(label=ni1))
ggp3+ggtitle(paste("fold=",k1))


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
      ang2[,mi2,drop=FALSE]
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


.combineAngles1<-function(angleH, incl, sumAngle, prev_signature, flags){ 
  #$types
  #flags = private$flags;
  topn = .readFlag(flags,'topn', 20)
  onlyAll = .readFlag(flags,'only_all',FALSE)
  angles1=angleH$angles;cols_incl1=angleH$cols_incl 
  nme_trans = names(angles1[[1]][[1]]); names(nme_trans) = nme_trans
  types = incl$types; names(types) = types
  excl = incl$excl
  # nmes_angs1 = names(angles1); names(nmes_angs1)=nmes_angs1
  #nme_t1 = nme_trans[[1]]; nme_p1 = names(angles1[[1]][[1]][[nme_t1]])[[1]]; inc1 = types[[1]]; jk=1
  comb_all2=lapply(nme_trans, function(nme_t1){
    nme_pow = names(angles1[[1]][[1]][[nme_t1]]); names(nme_pow)=nme_pow
    lapply(nme_pow, function(nme_p1){
      comb_all=lapply(types, function(inc1){
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
     #   excl1 = excl[unlist(lapply(excl, function(ex) ex[3]==nme_t1 && ex[1] == inc1 && ex[4] ==nme_p1))]
        if(length(excl)>0){
          mi2 = match(excl, names(col_incl))
          col_incl[mi2[!is.na(mi2)]]=FALSE
        }
        cs[col_incl]
      })
      top_angles=whichpart1(comb_all, n=topn, return_scores=TRUE)
      t1 = .merge1_new(lapply(top_angles, function(ta){
        data.frame(list(names = names(ta), value=ta))
      }),addName="data_type")
      if(is.null(t1)) return(t1)
      
      t1 = t1[order(t1$value),]
      signature=if(prev_signature=="") t1$names else paste(prev_signature, t1$names, sep=";")
      
     t2 =   t1 |> tibble::add_column(sumAngle, signature);
     subset(t2, value<999)
    })
  })
  comb_all2
}

## this is a class which holds a data object and interacts with the coordination node
.getAllSparseMatrices<-function(data, hasNA=TRUE, convertToBigMatrix=FALSE, min_variance =0.001, max_na_proportion=0.99){
  rn = unlist(lapply(data, function(d1) rownames(d1)))
  rn = rn[!duplicated(rn)]
  lapply(data, function(mat){
   m1 =  .getSparseMatrices(mat, hasNA=hasNA, convertToBigMatrix = convertToBigMatrix,rn = rn)
   
     sv = sparse_variance(m1$matrix)
     na_cnt = colSums(m1$matrixNA)/nrow(m1$matrixNA)
     m2 = lapply(m1, function(m11) m11[,!is.na(sv) & sv>min_variance & na_cnt<=max_na_proportion, drop=FALSE])
  m2
  })
  
}
##this function removes NAs
## if no NA matrixNA is just empty matrix
.getSparseMatrices<-function(mat, hasNA=TRUE, convertToBigMatrix=FALSE,rn = rownames(mat)){
  mi1 =  match(rownames(mat), rn)
  mi0 =  match(rn,rownames(mat))
  newNA=TRUE
  if(length(mi1)==length(rn)){
    if(max(abs(apply(cbind(mi1, 1:length(rn)),1,diff)))==0) newNA=FALSE
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
                  matrixNA = Matrix(0,nrow(mat) , ncol(mat), sparse = TRUE)
      )
    }
  }else{
    if(newNA){
      m1=apply(mat,2,function(v){
        v1 = v[mi0]
        mv = mean(v, na.rm=TRUE)
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
        mv = mean(v, na.rm=TRUE)
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
      res1 = list( matrix = Matrix(m1, sparse=TRUE),matrixNA = Matrix(m2, sparse=TRUE)) 
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
      if(sum(abs(vals-round(vals)), na.rm=TRUE)<1e-9) {
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
      if(sum(abs(vals-round(vals)), na.rm=TRUE)<1e-9){
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
    type="character",  
  plot_results="list",
  
  levs="list",
   data_id="character",
    dbDir="character",
  na_inds = "list",  ## these are for imputing uncertain values
  alt_inds = "list",
  mult="numeric",
  na_probs_ordered="list",
  original_inds ="list",  ## what is the indices in new matrix of the na_inds original values
  original_rows = "list", ## what is indices of original rows
  #returns the weights of the uncertain samples ordered
  sample_na_weights=function(){
    mult = private$mult
    if(length(private$na_inds)==0) return(c())
    w = private$data$weights[private$na_inds]
    for(alt in private$alt_inds){
      w = cbind(w,private$data$weights[alt])
    }
    sort(apply(w, 1, function(v) max(v)) )/mult
  },
  #auc for na samples
  sample_auc=function(pr){
    na_inds = private$na_inds
    y = self$y()[[1]]
    orig_inds = private$original_inds
  
    #y = y_orig ## already na_inds
    #df1 =cbind( pr[[1]][na_inds,,drop=FALSE],y[[1]][na_inds,,drop=FALSE])
    rocs =lapply(1:ncol(y),function(jk){
        roc1 = roc(  y[orig_inds,jk],
                     pr[[1]][na_inds,jk])
        roc1$auc
    })
    rocs
  },
   makeModels=function(vars2, inds){
     phens = private$phens; flags = private$flags
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
   select_k=function(analysis, k1,
                     vars_l_todo
                   ){
     verbose=getOption("verbose",FALSE)
     if(is.null(private$phens)){
       stop("run  dh$update(phens, flags, transform_x) before running select")
     }
     expt_id = super$expt_id();
     show_pvalue_plots=.readFlag(flags, "show_pvalue_plots",FALSE) 
     stop_y = .readFlag(flags, 'stop_y',"rand")
     logpvthresh = log(.readFlag(flags,"pthresh",0.1))
     beam= .readFlag(flags,"beam",1)
     comb20 = NULL;
     private$plot_results[[k1]] = list()
     # vars_l = analysis$nextVars(expt_id, flags)
     nvar=0;
     while(length(vars_l_todo$todo1)>0 ){
      # direct = direction[[min(nvar+1, length(direction))]]
       comb2_new1 = try(self$multiAnglesAndPv(comb20 , k1,expt_id, vars_l_todo))
       if(inherits(comb2_new1,"try-error")) break;
        comb21 =  lapply(comb2_new1, function(x) x$pvs)
        if(length(unlist(comb21))==0){
          vars_l_todo$todo1 = vars_l_todo$todo1[-1]
          if(length(vars_l_todo$todo1)==0) vars_l_todo$stop=TRUE;
          next;
        }
##         private$savePvals(flags,phens,k1, dh$name(), vars_l_todo$vars_l,comb21)# no need to save here, just keep
         private$savePvals(k1, dh$name(), vars_l_todo$vars_l,comb21)# no need to save here, just keep
         #try(self$multiAnglesAndPv(comb2, phens, k1,flags,expt_id, vars_l_todo))
      
       data_nme=self$name();
#       comb2_news1 = list(comb21); names(comb2_news1) = self$name()
       comb2_new = list(comb21);
       vars_l_todo_new=analysis$savePvalsAndNextVars(vars_l_todo,comb2_new,data_nme,  k1)
      
       if(vars_l_todo_new$moveNext){
       
         comb20 = comb21
       }
       vars_l_todo = vars_l_todo_new
       nvar = length(vars_l_todo$vars_l[[1]]$var)
       if(verbose) print(names(vars_l_todo$vars_l))
       if(length(vars_l_todo$vars_l[[1]]$var_names)>=flags$max) break;
       
      
     }
          
     vars_l_todo$vars_l
   },
   var_thresh = function(qq_t){
     lapply(private$data$vars, function(v) quantile(v, qq_t))
   },
 
   updateTrain=function(  verbose=FALSE, force=FALSE){
    
    
     private$data$updateTrain( private$phens,private$flags,private$transform_x, verbose=verbose, force=force)
   },
   updateLOOC=function( varn=c(),verbose=FALSE, force=FALSE){
     private$data$updateLOOC( private$phens,private$flags,varn=varn,force=force, verbose=verbose); ### update training object - updates all
   },
   findPrev=function(comb20, expt_id, prev_i, k){
     if(is.null(private$sigs)){
       nmes= unlist(lapply(prev_i$var_names, function(x) paste(x,collapse=".")))
       if(length(nmes)==0)nmes="empty"
       nme1 = strsplit(nmes[length(nmes)],"\\.")[[1]]
       if(length(nme1)>4){
         nme1[4] = paste(nme1[4:length(nme1)], collapse=".")
       }
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
   res_inner=function(comb_,prev_i2, k, expt_id){
     phens = private$phens;
     flags = private$flags;
     
     nme_comb = names(comb_); names(nme_comb) = nme_comb
     #nme_c1 = nme_comb[[2]]; nme_p1 = names(comb_[[nme_c1]])[[3]]; ik=1
     res_inner=lapply(nme_comb, function(nme_c1){
       nmesp1 = names(comb_[[nme_c1]]); names(nmesp1) = nmesp1
       lapply(nmesp1, function(nme_p1){
         #print(nme_p1)
         comb = comb_[[nme_c1]][[nme_p1]]
         if(length(comb)==0 || nrow(comb)==0) return(NULL)
         num_pvals1 = nrow(comb)
         inds1p = 1:num_pvals1; names(inds1p) = comb$names[1:length(inds1p)]
         nxt_vars = lapply(inds1p, function(ik){
          # print(ik)
           b_i_name = c(comb$data_type[[ik]], comb$names[[ik]], nme_c1,nme_p1)
           angle=comb$value[[ik]]
           if(!is.null(flags$angles_only) && flags$angles_only){
             b_i = private$data$convert(b_i_name)
             nv = list(angle = angle, var = c(prev_i2$var, list(b_i)), 
                       angles = c(prev_i2$angles, angle),
                       var_names = c(prev_i2$var_names, list(b_i_name)),
                       varnames = c(prev_i2$varnames, paste(b_i_name, collapse="."))
             )
            
             #  print(nv)
             nv$sumAngle = sum(nv$angles);
             
           }else{
             nv= private$getPvsAll(prev_i2, b_i_name,k,  prev_i2$Wall, angle=angle)
             if(is.na(nv$sumPv)) stop(paste(b_i_name, sep=","))
            }
           names(nv$var_names) = nv$varnames;
           if(inherits(nv,"try-error")) {
             warning(paste(nme_c1, "error"))
             return(NULL)
           }
           nv
         })
         
       })
     })
     res_inner   
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
##equivalent of combinedAngles1 but for predefined subsignature
predefined=function(incl1,prev_signature, sumAngle){
  dt1 = lapply(private$data$data, function(d) which(incl1 %in% colnames(d)))
  dt1 = dt1[unlist(lapply(dt1, length))>0]
  lapply(private$transform_x, function(tx){
    lapply(tx$params, function(tx1){
      .merge1_new(lapply(dt1, function(ind){
        signature=        apply(cbind(prev_signature, incl1[ind]),1,paste,collapse=";")
        
        data.frame(list(names = incl1[ind],value = 0,sumAngle=sumAngle, signature=signature ))
      }),addName="data_type")
    })
  })
  
},
   combinedAngles=function( varnames, incl, k, g_incl, qq_t, sumAngle,  addPlot=FALSE){ #phens, varnames, incl=incl, k=k, type=type
      type=private$type
      phens = private$phens;
      flags = private$flags; 
     prev_signature =paste(unlist(lapply(varnames, function(vn)vn[2])),collapse=";")
     var_t = private$var_thresh(qq_t)
     if(length(which(names(private$data$data) %in% incl$types))==0) stop("incl does not match data")
     v2 = unlist(lapply(varnames, function(x) x[[2]]))
     incl2 = incl$incl[!(incl$incl %in% v2)] 
     if(length(incl2>0)){  ## use pre-defined
        comb_angle1 = private$predefined(incl2, prev_signature, sumAngle)
     }else{
     
         angles=  private$data$getAngles1(phens,varnames,incl,k=k, type=type)
         angles =angles[ unlist(lapply(angles, length))>0]
         angleH=list(angles=angles,
                     cols_incl = private$data$cols_incl(var_t,incl$types, g_incl,excl=varnames)) ### fix 
         comb_angle1 =  .combineAngles1(angleH, incl,  sumAngle, prev_signature,flags)
         if(addPlot){
           all_angles = .extrAngles(angleH,comb_angle1, incl$types)
            attr(comb_angle1,"all")=all_angles;
         }
     }
     comb_angle1
   },
   getPvsAll=function( prev_i2, b_i_name,k, #   prev_i = vars_l1[[nmed]]
                      Wall, # =lapply(subphens, function(f) matrix(nrow=0,ncol=0)),
                      angle=0){
   flags = private$flags; phens = private$phens;
     project=.readFlag(flags,"project",TRUE)
     useoffset=.readFlag(flags,"useoffset",TRUE)
     useglm = .readFlag(flags,'useglmnet',TRUE)
     d = private$data
     family = strsplit(names(phens)[[1]],"\\.")[[1]][1]
     if(family=="multinomial") useoffset=FALSE
     
     prev_i1 = d$makeNextModel(prev_i2,b_i_name,phens,k, family, ypred=NULL, 
                               project=project, useglm=useglm, logpthresh =0, useoffset=useoffset)
     prev_i1$angle= angle
     prev_i1$angles = c(prev_i2$angle, angle)
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
   #' @param certainty the certainty of each observation, expressed as a probability of being correct. Default is to treat each observation as certain
      #' @param nme The name of the data objetct. 
   #' @param flags list of options, described in the vignette
   #' @param transform_x a transformation object from the function getTransformation
   #' @param family the statistical family of phenotype y, can be calculated by getFamily
   #' @param dbDir dir for database to store results to speed up re-reruns.  Can be NULL
     initialize=function(
    data,
    y,
    certainty = rep(1, nrow(y)),
    nme,
    flags ,
    family= getFamily(y),
         dbDir=NULL
   ){
       
    super$initialize(nme,dims =lapply(d$data, dim),  flags=flags, dbDir=dbDir);
       colnames(y) = gsub("\\.","_",colnames(y))  ## no . allowed
      for(k in 1:length(data)){
        colnames(data[[k]]) = gsub("\\.","_",colnames(data[[k]]))  ## no . allowed
        
      }
       weights = rep(1, nrow(y))
       private$plot_results=list();
    #   private$transform_x = transform_x;
    private$levs = lapply(y,function(yc) levels(yc) )       
    private$mult= .readFlag(flags,"mult",100)
    private$original_rows = 1:nrow(y)
    
    d = expandData(data, y, certainty,family,  weights, mult = private$mult)
    private$na_inds = d$na_inds; private$alt_inds = d$alt_inds;
   private$original_inds = d$original_inds
    memDir=NULL
    useDB=!is.null(dbDir);
    convertToBigMatrix=FALSE #.readFlag(flags,"covertToBigMatrix", FALSE)
    hasNA=.readFlag(flags,"hasNA", TRUE)
    preprocessed=.readFlag(flags,"preprocessed", FALSE)
    max_na_proportion =.readFlag(flags,"max_na_proportion",0.99)
    min_variance =.readFlag(flags,"max_na_proportion",0.001)
    
    mat = .getAllSparseMatrices(d$dataset,hasNA=hasNA, convertToBigMatrix=convertToBigMatrix,min_variance = min_variance, max_na_proportion=max_na_proportion)
    
    private$dbDir = dbDir
    
    private$data = 
      dataObj$new(mat, private$nme,dbDir,flags, 
                  incl_full=TRUE,seed = getOption("seed",42), memDir=if(is.null(memDir)) NULL else paste(memDir, private$nme,sep="/"))
    private$type="slow1"
    types_all = getOption("types_all",names(private$data$data))
    names(types_all) = types_all
    batch=.readFlag(flags, "batchsize",0)
    all_v_all = .readFlag(flags,"all_v_all",FALSE)
    one_v_rest = .readFlag(flags,"one_v_rest",FALSE)
    nrep = .readFlag(flags,"nfold",if(batch>0) 0 else 1)
    if(nrep>0 && batch>0)warning("only one of nfold or batchsize should be non zero")
    pheno_balance=.readFlag(flags,"pheno_balance",NULL)
    varn = getOption("varn",c())
    if(is.null(weights)) stop("data weights need to be defined");
    
    #invisible(lapply(1:length(datas), function(ik) {
     # family = families[[ik]]
        private$data$updateY(d$y,d$weights, preprocessed=preprocessed, family=family, CHECK=TRUE, all_v_all=all_v_all, one_v_rest = one_v_rest)
    private$na_probs_ordered =private$sample_na_weights()
    #  private$sigsdir=paste(dbDir,"fspls_signatures1",sep="/")
      #dir.create(private$sigsdir, recursive=FALSE, showWarnings=FALSE)
     # dims1 = list(private$data$dims()); names(dims1) = private$nme
      #private$sigs=   if(exists("dbConnect") && useDB) sigEnv$new(private$sigsdir,nme,flags, dims1, clear=FALSE) else NULL;
      
    
  },

  
  #' @description Calculate the angles and pv across multiple phenotypes.  This is an internal function and should not need to be called by user
  #' @param comb20 values from previous iteration
  #' @param phens1 phenotypes being used
  #' @param k1 k1 is the current fold
  #' @param flags list of options
  #' @param expt_id and ID for the experiment (can be 0)
  #' @param vars_l_todo  the variables under consideration
    multiAnglesAndPv=function(comb20,  k1, expt_id, vars_l_todo ){
      phens1 = private$phens;
      flags = private$flags; 
#      self$update(phens, flags, transform_x);
      show_pvalue_plots=.readFlag(flags, "show_pvalue_plots",FALSE) 
      verbose=.readFlag(flags,'verbose',FALSE)
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
      comb_=private$combinedAngles(varnames, incl, k1,  g_incl, qq_t, sumAngle,
                                   addPlot=show_pvalue_plots) ;
      if(length(comb_)==0) stop("length zero")
      # return(list(comb_angle1, all_angles));
      
      if(saveAngles) return(comb_)
      #comb_ = private$anglesAndPv(phens, prev_i, incl, k1, g_incl, qq_t, flags,expt_id, saveAngles=saveAngles, verbose=verbose)
      ri = private$res_inner( comb_,prev_i2,k1, expt_id)
     super$savePvals( ri, varnames,k1,useCurrVarnames=TRUE)
      
      
      list(angles = comb_, pvs = ri) ;#private$simplify(ri))
    }))
     if( show_pvalue_plots){ ## just prints the plot to screen, or to pdf
       
         ggps=try(plot_angle_vs_pv(comb2_new,1, k1))
        print(ggps)
       
      #  if(length(vars_l_todo$todo1)>0){
       
       # }
     }
    nvar = length(vars_l_todo$vars_l[[1]]$var)+1
    if(nvar==1){
      private$plot_results[[k1]] = list()
    }
    private$plot_results[[k1]][[nvar]] = mergeAll(comb2_new, flags$beam)
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
      warning("need to set drop=TRUE if you are sure, this will delete all saved signatures")
    }
    
  },
  
  #' @description returns data on trajectory
  #' @returns a dataframe with information to plot trajectory
  get_trajectory=function(){
    plot_results = private$plot_results
    if(length(plot_results)==0) return(NULL)
    
    names(plot_results) = 1:length(plot_results);
    .merge1_new(lapply(plot_results, function(pr){
    
        names(pr) = 1:length(pr)
        pr[[1]]$prev=""
        comb_plot = .merge1_new(pr, addName="nvar");
        comb_plot = comb_plot[comb_plot$func!="rand",,drop=FALSE]
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
    }), addName="repeat")
  },
  #' get the data types
  #' data_nmes a list of names of the data names to include, optional
  #' inds a list of positive and negative inds, optional
  #' direction a direction object
  #' incl a list of variables to include 
  #' excl a list of variables to exclude 
  #' max a list of maximum number of variables for each data type
  #' @returns an object which is used to specify both directionallity and data types
  data_types=function(
    data_nmes=list(all = names(private$data$data)),
    inds =lapply(data_nmes, function(x) list(pos_inds = c(), neg_inds = c(), max = 1000)),
    direction = lapply(inds, private$data$getDirection),
    incl = c(),
    excl = c()
  ){
#  match(data_nmes,     
    nmes = names(inds); names(nmes) = nmes;
    all = lapply(nmes, function(nme1){
      dt1 = data_nmes[[nme1]];
   
      dir1 = direction[[nme1]]
      max = inds[[nme1]]$max
      if(is.null(max)) max = 100
      if(is.null(dt1)){
        res =  lapply(data_nmes, function(dt2){
          list(direction = dir1, types = dt2, max = max, excl = excl, incl = incl, nvar=0)
        })
      }else if(is.null(dir1)){ ## use all directions
           res =  lapply(direction, function(dir2){
              list(direction = dir2, types = dt1, max = max, excl = excl, incl = incl, nvar=0)
            })
      }else{
        res = list(list(direction = dir1, types = dt1, max = max, excl = excl, incl = incl, nvar=0))
        #names(res) = nme1
      }
      res
    })
    all1 = unlist(all, recursive=F)
    lapply(all1, function(a1){
      if(length(a1$incl)>0){
       a1$incl=unlist(lapply(a1$types, function(t1) {
          cn =  colnames(private$data$data[[t1]])
          cn[cn%in% a1$incl]
        }))
      }
      if(length(a1$excl)>0){
        a1$excl=unlist(lapply(a1$types, function(t1) {
          cn =  colnames(private$data$data[[t1]])
          cn[cn%in% a1$excl]
        }))
      }
      a1
    })
    
  },
  #' @description split dataset into smaller datasets
  #' @param proportions what proportions to split into
  #' @returns a list of dataH objects with data partiioned according to proportions
 split=function(proportions = c(0.5,0.5)){
   datas = private$data$split(proportions);
   mats = datas$mats;
   ys = datas$ys;
   nme_d = names(mats); names(nme_d) = nme_d
   all_v_all = .readFlag(private$flags,"all_v_all",FALSE)
   one_v_rest = .readFlag(private$flags,"one_v_rest",FALSE)
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
    maxpheno=1e9;sep=FALSE; sep_group =FALSE;exclude=NULL; code=NULL; memb=NULL
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
   res = 1:ncol(private$data$looc$incl)
   names(res) = 1:length(res)
  
   names(res)[length(res)]="full";
   res
 },


 #' @description update the phenotypes without remaking the entire object. You can provide many phenotypes in the initialisation stage, but only consider a subset in model fitting stage in this way.Must run this before select
 #' @param phens phenotypes
 #' @param flags list of options
 #' @param transform_x transformation object 
 
 update=function(phens=self$pheno()$all, flags=private$flags, transform_x=fromJSON(flags$transform_x), 
                 data_types = self$data_types()){
   flags = super$updateExpt(phens, flags, transform_x, data_types);
    verbose=.readFlag(flags,'verbose',FALSE)
    force=.readFlag(flags,'force',FALSE);
    #flags$transform_x = transform_x;
  #  super$updateExpt(phens, flags)
     private$updateLOOC(verbose=verbose,force=force)
     private$updateTrain(verbose=verbose, force=force)
    
   ##updated after updateLOOC
  # nreps1 =self$nreps()
  # nreps = 1:nreps1
  # names(nreps) = nreps
  # nreps
 },
 
 
#' @description main function for variable selection
#' @param analysis an analysisEnv object

 select=function( 
                 analysis =analysisEnv$new(flags=private$flags, dbDir=NULL)
                              ){#c(y="function(y) y","function(y) y")

   if(is.null(private$phens)) stop("need to update first");
    nreps =self$nreps();
     variables = super$loadVars()
     if(!is.null(variables)) return(variables)
   
   vars_l_todo = analysis$getTodo(private$flags, private$phens)
   variables=lapply(nreps, function(k1){
     if(getOption("verbose",FALSE)) print(paste("cv",k1,"of",length(nreps)))
     private$select_k(analysis, k1, vars_l_todo)
               })
   attr(variables, "phens")=private$phens; attr(variables,"transform_x") = private$transform_x; attr(variables, "flags") = private$flags
   
   super$saveVars(variables);
  variables;
  
 },
 #merge22=function(comb2, vars_l){
#  lapply(vars_l, function(vars_l1){
#      vars_l1$var_names
#  })
# },
 
#' @description provides access to internal storage of phenotype data
#' @param phens list of phenotypes
#' @returns list of matrices
y=function(phens = self$pheno()[[1]]){
  y1 = private$data$y
  out1=.lapply_nme(y1, function(nme) {
   
    y2 = y1[[nme]]
    y2[,colnames(y2) %in% phens[[nme]],drop=FALSE]
  })
  out1
},


#' @description ggplot to visualise predictions
#' @param all_modelsh fitted models from makeAllModels
#' @param update whether to update
#' @param liab return liability score, or probability (for binomial, ordinal multinomial)
#' @returns  a ggplot
plotPredictions=function(all_modelsh,
                         update=FALSE,
                         liab=TRUE){
  if(update) self$update(all_modelsh$phens, all_modelsh$flags, transform_x =  fromJSON(all_modelsh$flags$transform_x))
  
 y = self$y(private$phens)
  preds= self$extractPredictions(all_modelsh, phens, flags, liab, transform_x = transform_x)
  
  #familys=names(preds[[beam]][[nv]][[cv]]);
  #family = familys[[1]];td=1
  ##beams = names(preds); names()  .merge_lapply_nme
toplot=.merge_lapply_nme(preds, "beam",function(beam){
 .merge_lapply_nme(preds[[beam]], "numvar",function(nv){
    cvs = preds[[beam]][[nv]]
   
    .merge_lapply_nme(cvs, "cv",function(cv){
      familys = preds[[beam]][[nv]][[cv]]
      .merge_lapply_nme(familys, "family",function(family){
        pred = preds[[beam]][[nv]][[cv]][[family]];
        y1 = y[[family]];
        todo =1:ncol(y1) ; names(todo) = colnames(y1);
        .merge1_new(lapply(todo, function(td){
          df = data.frame(pred[,td],y1[,td])
          names(df) = c("prediction","value")
          df;
        }),addName="subpheno")
        
      })
    })
  })
})
  
 
  
  toplot$numvar = factor(toplot$numvar, levels = sort(unique(as.numeric(toplot$numvar))))
  toplot1 <- toplot |> 
    tidyr::unite("cv_family_beam", cv,family,beam, remove = FALSE)
  ggplot(toplot1, aes(x=value, y=prediction, color=subpheno, shape=beam))+geom_point()+facet_grid("numvar~cv_family_beam")
  
},



#' @description extract the predictions for the fitted models
#' @param all_modelsh fitted models from makeAllModels
#' @param liab return liability score, or probability (for binomial, ordinal multinomial)
#' @returns  a table with results
extractPredictions=function(all_modelsh,
                            liab=TRUE){
  phens = private$phens; flags = private$flags; transform_x = private$transform_x;
  private$updateLOOC()
    all_models_y0 = all_modelsh$models#[[mod_nme]]
  # eval1 =  .merge1_new(lapply(nme_d2, function(nme1){
  #print(nme1)
  d = private$data
  #all_models_y = all_models_y0[[1]]
  predictions0 = lapply(all_models_y0, function(all_models_y){
    d$extractPredictions(all_models_y, phens, flags, transform_x = transform_x, liab= liab)
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
plotData=function(variables, all_types=FALSE, violin=FALSE, assoc=FALSE, update=FALSE){
  
  attrs = attributes(variables)
  if(update) self$update(attrs$phens, attrs$flags, transform_x =  attrs$transform_x)
  
  vars_all=super$integrate(variables)

  df4= #.merge1_new( 
   # lapply(private$datas, function(d) 
      private$data$plotData(vars_all, phens1 = private$phens, all_types=all_types, transform_x = private$transform_x, 
                            violin=violin, assoc=assoc)
               #     addName="dataset")
  df4 = df4|>tibble::add_column(dataset=private$nme);
  facet= if(!is.null(df4[['transform']]) ) "transform~pheno" else "pheno"
  df4$y = factor(df4$y)
  gene_levs = unlist(lapply(unlist(vars_all[[1]]$variables,recursive=FALSE), function(x) x[[2]]))
  gene_levs = gene_levs[!duplicated(gene_levs)]
  df4$gene = factor(df4$gene, levels = gene_levs)
  df5=df4|> separate(col="pheno",sep="\\.", into=c("family","pheno"),fill="right",extra="drop")
  df6=subset(df5, family=="gaussian")
  df7=(subset(df5, family!="gaussian"))
  # color= if(length(unique(df4$dataset))>1) "dataset" else #"pheno"
  ggp = NULL; ggp1 = NULL
  if(nrow(df7)>0){
    ggp<-ggplot(df7, aes(x=y, y=value, color=y, shape=data, linetype=dataset))+facet_wrap(facet, scales="free");#+ggtitle(unlist(phens1))
    if(violin){
      ggp<-ggp+ggplot2::geom_violin()+geom_point()
    }else{
      ggp<-ggp+ggplot2::geom_boxplot()
    }
    ggp<-ggp+facet_wrap("gene")
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
        q1=quantile(df6_3$value,na.rm=TRUE,probs = prbs)
        df_=data.frame(t(data.frame(q1)))
        df_
      }),addName="y")
    }), addName="comb")
    quants1 = quants |> separate("comb", sep="__", into=c("gene","data","family","pheno1","dataset"))
    names(quants1) = gsub("\\.","",names(quants1))
    quants1$y = as.numeric(quants1$y)
    ggp1<-ggplot(quants1, aes(x=y, y=X50, ymin=X33, ymax = X66,color=gene, fill=gene,shape=data, linetype=dataset))+facet_wrap(facet, scales="free");#+ggtitle(unlist(phens1))
    ggp1<-ggp1+geom_line()+geom_ribbon(alpha=0.1)
    ggp1<-ggp1+facet_wrap("gene")
  }
  
  list("binomial"=ggp,"gaussian"=ggp1)
},

#' @description update the weights based on predictions from the model
#' @param all_models the output of dataH function getAllModels
#' @returns list of auc beta and updated sumdiff of weights used as stopping criteria
updateWeights=function(all_models){
  predictions =     self$extractPredictions(all_models, liab=TRUE) ## extract predictions in liability space
  pr1 = predictions[[1]];
  fulls = lapply(pr1, function(x) x$full)
  fulls = fulls[unlist(lapply(fulls, length))>0]
  
  pr =fulls[[length(fulls)]]
  #calcAUC=TRUE;calcWeights=TRUE;
  auc = private$sample_auc(pr)
  am = all_models$models[[1]]
  betas = (am[[length(am)]]$full$betas[[1]])
  rownames(betas) =  names(am[[length(am)]]$full$var_names)
  mult = private$mult;
  y = self$y()[[1]]
  na_inds = private$na_inds; alt_inds = private$alt_inds; original_inds = private$original_inds
  #pr[[1]][na_inds,]
  #for(kk in 1:ncol(y)){
  #  zero_inds = which(y[,kk]==0)
  #  pr[[1]][zero_inds,kk] = 1-pr[[1]][zero_inds,kk]
  #}
  binom = names(pr)[[1]]=="binomial"
  pr2 = pr[[1]][na_inds,,drop=FALSE]
  if(binom){
    pr2 = cbind(1-pr2, pr2)
     
  }
  private$data$weights[na_inds] =round(pr2[,1]*mult)
    
  for(jk in 1:length(alt_inds)){
    private$data$weights[alt_inds[[jk]]] = round(mult*pr2[,jk+1])
  }

  probs_ordered = private$sample_na_weights()
  diff =(probs_ordered - private$na_probs_ordered)
  private$na_probs_ordered = probs_ordered
  
  list(auc = auc, betas = betas, sumdiff = sum(diff))
},


#' @description return updated y after iterative imputation
#' @returns updated y, updated certainty and indices which were updated
getYNew=function(){
  #predictions =     self$extractPredictions(all_models, liab=TRUE) ## extract predictions in liability space
  #pr =predictions[[1]][[length(predictions[[1]])]]$full
  na_inds = private$na_inds
  mult = private$mult
  w = private$data$weights[private$na_inds]
  for(alt in private$alt_inds){
    w = cbind(w,private$data$weights[alt])
  }
  w =  w/mult
  y = self$y()[[1]]
  if(names(self$y())[[1]]=="binomial"){  ## to make it like multinomial
    y = cbind(1-y,y)
  }
  y2 = apply(y,1,function(v){
    which.max(v)
  })
  y_orig = y2[private$original_inds]
  y2 = y2[private$original_rows]
  y_new=t(apply(w, 1, function(v){
   c(which.max(v), max(v))
  }))
  
  y2[na_inds] = y_new[,1]
   
  certainty = rep(1, length(private$original_rows))
  certainty[na_inds] = y_new[,2]
  error_rate = sum(abs(y2[na_inds] - y_orig))/ length(na_inds)
    levs = private$levs
    if(length(levs[[1]])>0){
     y2 =  factor(y2, levels = sort(unique(y2)), labels = levs[[1]])
    }
  list(y=y2, certainty =certainty, na_inds = na_inds, error_rate = error_rate)
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
#' @param variables list of variables selected by select method
#' @param update whether to automatically update phens, transform_x and flags , default TRUE
#' @returns fitted models
makeAllModels=function(variables,
                       update=FALSE
                       ){
  attrs = attributes(variables)
  if(update) self$update(attrs$phens, attrs$flags, transform_x =  attrs$transform_x)
  
  vars_all=super$integrate(variables)
#  phens=vars_all[[1]]$phens, flags=vars_all[[1]]$flags, 
  useDB = super$useDB;
  flags = private$flags;phens = private$phens;
  verbose=.readFlag(flags,'verbose',FALSE); max = .readFlag(flags,'max',1e6)
  all_models = super$loadModels();
  if(!is.null(all_models)) return(all_modelsh);
  
  private$updateLOOC()
  logpthresh= log(.readFlag(flags,"pthresh",1e-3))
  project=.readFlag(flags,"project",TRUE)
  beams = names(vars_all); names(beams)=beams
  all_models_full=lapply(beams, function(beam){
    vars_all0 = vars_all[[beam]]
    vars = vars_all0#[[nme_v_all]]
  all_models = list()
  variables = vars$variables
  var_inds = vars$inds
  rem_inds = self$nreps()
 
  all_models = private$makeModels(list(),rem_inds )

  if(length(variables)==0)  return(all_models) ;#return(list(models=all_models, flags = flags, phens = phens, trainedOn=private$nme))
  ord = order(unlist(lapply(variables, length)),decreasing=TRUE)
  variables = variables[ord]
  var_inds = var_inds[ord]
  #v_nme = names(variables)[1]; #max=10; verbose=TRUE; k=1;variables =vars_all$variables; 
  
  for(v_nme in names(variables)){
    if(verbose)print(v_nme)
    vars2 = variables[[v_nme]]
    vars2 = vars2[1:min(length(vars2), max)]
    inds =var_inds[[v_nme]]
    nme_ = paste(names(vars2),collapse=";")
    models1 = all_models[[nme_]]
    if(is.null(models1)){
      models1 = private$makeModels( vars2, inds)
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
        models2 = private$makeModels( vars2, inds)
                                  
        for(k in 1:length(models2)){
          all_models[[names(models2)[[k]]]] = models2[[k]]#[[p_nme]]#[[p_nme]]
          
        }
      }else{
        subinds= inds[which(is.na(match(names(inds), names(all_models[[nme_]]))))]
        if(length(subinds)>0){ ## missing inds
          models3 = private$makeModels( vars2, subinds)
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
  super$saveModels(all_models_);
  #combined_models
  all_models_
},
#' @description evaluate the fit of models
#' @param all_modelsh  models fitted from makeAllModels
#' @param update whether to automatically update phens, transform_x and flags , default TRUE
#' @returns evaluation of models
evaluateAllModels=function(all_modelsh, update=TRUE){ ## different folds with same variables
  if(update) self$update(all_modelsh$phens, all_modelsh$flags, transform_x =  fromJSON(all_modelsh$flags$transform_x))
  flags = private$flags; phens = private$phens;
  eval1 = super$loadEval();
  if(!is.null(eval1)) return(eval1);
  
  
    private$data$updateTransform(private$transform_x)
    verbose=.readFlag(flags,"verbose",TRUE)
  inv_transform_x=FALSE
  private$updateLOOC()
  if(length(all_modelsh$models)==0) return(NULL)
 
  all_models_y0 = all_modelsh$models#[[mod_nme]]
     d = private$data
    # all_models_y=all_models_y0[[1]]
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
  #}),addName="transform_x")
  if(is.null(eval1)) return(NULL)
  #  eval1 = subset(eval1, model!="avg")

  eval2 = eval1|> pivot_wider(names_from="submeasure") #|> tibble::add_column(transform_x=strsplit(transform_x[[1]]," ")[[1]][2])
  
  #  isfull=eval2$model %in% full_model_nmes
  #  eval2|>tibble::add_column(isfull=isfull)
  
  #print("HH")
  #  if(TRUE) return(eval2)
  eval3 = .calcEval1(eval2, rename=FALSE)
  eval4 = eval3 |> tibble::add_column(variable= unlist(lapply(eval3$model, function(x){
    if(x=="cv" || x=="") return("");
    x1 =rev(strsplit(x,";")[[1]])[1]
    x1
  })))
  eval4$variable[is.na(eval4$variable)]=""
  eval5 = eval4 |> tidyr::separate("variable", sep="\\.", into=c("type","variable","func","param"),remove=TRUE)
  super$saveEval(eval5);
  eval5
}
    
)
)
  