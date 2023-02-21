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

transfer_predict <- function(object, newdata) {
  result <- list()
  if (object@method == "zeros") {
    result <- zeros_predict(newdata)
  } else if (object@method == "gaussian_latent") {
    result <- gaussian_latent_predict(object, newdata)
  }
  
  result
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

#' @export
setMethod("predict",  c(object = "transfer_model"), transfer_predict)

#' @export
setGeneric("values", function(x) standardGeneric("values"))
setMethod("values", "ts_inter_single", function(x) x@values)

#' @export
setGeneric("values<-", function(x, values) standardGeneric("values<-"))
setMethod("values<-", "ts_inter_single", function(x, values) {
  x@values <- values
  x
})

setMethod("[[<-", "ts_inter", function(x, i, j, value) {
  x@series[[i]] <- value
  x
})

#' @export
setGeneric("interventions", function(x) standardGeneric("interventions"))
setMethod("interventions", "ts_inter_single", function(x) x@interventions)

#' @export
setGeneric("interventions<-", function(x, value) standardGeneric("interventions<-"))
setMethod("interventions<-", "ts_inter_single", function(x, value) {
  x@interventions <- value
  x
})

