setClass(
  "ts_inter_single", 
  slots = c(
    values = "matrix",
    time = "numeric",
    interventions = "matrix"
  )
)

setClass(
  "ts_inter",
  slots = c(series = "list")
)

setClass(
  "transfer_model",
  slots = c(
    parameters = "ANY",
    method = "character"
  )
)

single_subset <- function(x, i, j, ..., drop = FALSE) {
  values <- x@values[i, j, drop = drop]
  time <- x@time[j]
  interventions <- x@interventions[, j, drop = drop]
  new("ts_inter_single", values = values, time = time, interventions = interventions)
}

multi_subset <- function(x, i, j, ..., drop = FALSE) {
  result <- list()
  for (k in seq_along(x)) {
    result[[k]] <- x[[k]][i, j]
  }
  
  new("ts_inter", series = result)
}

setMethod("length", "ts_inter", function(x) length(x@series))
setMethod("nrow", "ts_inter_single", function(x) nrow(x@values))
setMethod("ncol", "ts_inter_single", function(x) ncol(x@values))
setMethod("dim", "ts_inter_single", function(x) dim(x@values))
setMethod("[", c("ts_inter_single", "numeric", "numeric"), single_subset)
setMethod("[", c("ts_inter_single", "numeric", "missing"), single_subset)
setMethod("[", c("ts_inter_single", "missing", "numeric"), single_subset)
setMethod("[", c("ts_inter", "numeric", "missing", "ANY"), function(x, i, j, ..., drop=TRUE) initialize(x, series=x@series[i]))
setMethod("[[", c("ts_inter", "numeric", "missing"), function(x, i, j, ...) x@series[[i]])
setMethod("[", c("ts_inter", "logical", "missing", "ANY"), function(x, i, j, ..., drop=TRUE) initialize(x, series=x@series[i]))
setMethod("[[", c("ts_inter", "logical", "missing"), function(x, i, j, ...) x@series[[i]])
#setMethod("predict",  c(object = "transfer_model"), transfer_predict)

setGeneric("values", function(x) standardGeneric("values"))
setMethod("values", "ts_inter_single", function(x) x@values)
setGeneric("values<-", function(x, values) standardGeneric("values<-"))
setMethod("values<-", "ts_inter_single", function(x, values) {
  x@values <- values
  x
})

setMethod("[[<-", "ts_inter", function(x, i, j, value) {
  x@series[[i]] <- value
  x
})

setGeneric("interventions", function(x) standardGeneric("interventions"))
setMethod("interventions", "ts_inter_single", function(x) x@interventions)
setGeneric("interventions<-", function(x, interventions) standardGeneric("interventions<-"))
setMethod("interventions<-", "ts_inter_single", function(x, interventions) {
  x@interventions <- interventions
  x
})

