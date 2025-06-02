modelsEnv<-R6Class("modelsEnv", public = list(
  models = "list",
  type="character",
  initialize=function(str = NULL){
    self$models=if(is.null(str)) list() else fromJSON(str)
  },
  json=function(){
    toJSON(self$models)
  },
  getModels=function(index="full"){
    res=lapply(models, function(x){
      mi=match(index, names(x))
      if(is.na(mi)) return(NULL)
      x[mi]
    })
    res = res[lapply(res, length)>0]
  }
  
)
)
