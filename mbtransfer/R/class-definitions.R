setClassUnion("data.frameOrNull", members=c("data.frame", "NULL"))

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
  slots = c(
    series = "list",
    subject_data = "data.frameOrNull"
  )
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

model_predict <- function(object, newdata) {
  fit <- object@parameters
  series <- list()
  new_interventions <- list()

  for (i in seq_along(newdata)) {
    split_data <- split_future(newdata[[i]])
    series[[i]] <- split_data$pre
    new_interventions[[i]] <- split_data$interventions
  }

  newdata <- new("ts_inter", series = series, subject_data = subject_data(newdata))
  model_predict_(fit, newdata, new_interventions)
}

train <- function(ts_inter, P = 1, Q = 1, ...) {
  train_data <- patchify_df(ts_inter, P, Q)
  fit <- list()
  for (j in seq_along(train_data$y)) {
    fit[[j]] <- xgboost(data = train_data$x, label = train_data$y[[j]], ...)
  }
  new("transfer_model", parameters = fit, method = method)
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
setMethod("predict",  c(object = "transfer_model"), model_predict)

#' @export
setMethod("train",  c(ts_inter = "ts_inter", P = "numeric", Q = "numeric"), train)

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

#' @export
setGeneric("times", function(x) standardGeneric("times"))
setMethod("times", "ts_inter_single", function(x) x@time)

#' @export
setGeneric("subject_data", function(x) standardGeneric("subject_data"))
setMethod("subject_data", "ts_inter", function(x) x@subject_data)

#' @export
setMethod("names", "ts_inter", function(x) names(x@series))
