
#' @importFrom mbtransfer mbtransfer
#' @importFrom mdsine mdsine
#' @export
train <- function(ts_inter, method = "mbtransfer", hyper = list(), ...) {
  if (method == "mdsine") {
    
    fit <- mdsine(ts, hyper$taxonomy, ...)
    result <- new("transfer_model", parameters = fit, method = method)
    
  } else if(method == "mbtransfer") {
    
    fit <- mbtransfer(ts, hyper$P, hyper$Q, ...)
    result <- new("transfer_model", parameters = fit, method = method)

  }
  result
}
