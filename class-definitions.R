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
    parameters = "list",
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

transfer_predict <- function(object, newdata, n_ahead = 1) {
  result <- list()

  if (object@method == "zeros") {
    for (i in seq_along(newdata)) {
      result[[i]] <- initialize(
        newdata[[i]], 
        values = matrix(0, nrow(newdata[[i]]), ncol(newdata[[i]]))
      )
    }
  }
  
  result
}

setMethod("length", c("ts_inter"), function(x) length(x@series))
setMethod("nrow", c("ts_inter_single"), function(x) nrow(x@values))
setMethod("ncol", c("ts_inter_single"), function(x) ncol(x@values))
setMethod("dim", c("ts_inter_single"), function(x) dim(x@values))
setMethod("[", c("ts_inter_single", "numeric", "numeric"), single_subset)
setMethod("[", c("ts_inter_single", "numeric", "missing"), single_subset)
setMethod("[", c("ts_inter_single", "missing", "numeric"), single_subset)
setMethod("[", c("ts_inter", "numeric", "missing", "ANY"), function(x, i, j, ..., drop=TRUE) initialize(x, series=x@series[i]))
setMethod("[[", c("ts_inter", "numeric", "missing"), function(x, i, j, ...) x@series[[i]])
setMethod("[", c("ts_inter", "logical", "missing", "ANY"), function(x, i, j, ..., drop=TRUE) initialize(x, series=x@series[i]))
setMethod("[[", c("ts_inter", "logical", "missing"), function(x, i, j, ...) x@series[[i]])
setMethod("[[<-", c("ts_inter_single", "matrix"), function(x, i, j, ...) x@series[[i]])
setMethod("predict",  c(object = "transfer_model"), transfer_predict)



