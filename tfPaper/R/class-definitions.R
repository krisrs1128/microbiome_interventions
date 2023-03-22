
#' @importFrom mbtrasnfer split_future
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

setClass(
  "transfer_model",
  slots = c(
    parameters = "ANY",
    method = "character"
  )
)

#' @export
setMethod("predict",  c(object = "transfer_model"), model_predict)
