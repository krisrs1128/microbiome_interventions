#' Obtain predictions from a fitted multinomial logit (ALR, CLR) Gaussian process model object and return a ts_inter object. 
#'
#' @param object a fitted object of class inheriting from "fido::basset".
#' @param newdata a ts_inter in which to look for variables with which to predict. 
#' @param design 
#' @importFrom fido predict alrInv
#' @export
fido_predict <- function(object, 
                         newdata, 
                         #design = paste0("~", "time", "+", "as.integer(as.factor(subject))","-1")){
                         design = "~ time") {
  
  data <- fido_data(newdata)
  X_predict <- t(model.matrix(formula(design), data=data$samples))

  predicted <- predict(object, X_predict)
  pred_median_alr <- apply(predicted, c(1,2), median)# taxa-1 x time
  pred_median_invalr <- alrInv(t(pred_median_alr)) |> 
    t()

  library_size <- colSums(data$Y)
  y_hat <- t(t(pred_median_invalr)*library_size) |> round(digits = 0)
  colnames(y_hat) <- colnames(data$Y)
  rownames(y_hat) <- rownames(data$Y)
  
  data$Y <- y_hat
  
  ts_from_dfs(
    reads = data$Y |> t(), 
    interventions = data$interventions |> column_to_rownames("sample"), 
    metadata = data$samples, 
    subject_data = data$samples[, colnames(newdata@subject_data)]
  )
}

#' @importFrom fido basset
#' @export
fido <- function(
    ts, 
    sigma, 
    rho, 
    design = "~ time") {
    #design = paste0("~", "time", "+", "as.integer(as.factor(subject))","-1")
  dat <- fido_data(ts)
  Y <- dat$Y
  D <- nrow(Y) # taxa
  N <- ncol(Y)
 
  X <-  t(
    model.matrix(formula(design), data=dat$samples
    )
  )
  
  Gamma <- function(X){
    Gamma_(X, sigma = sigma, rho = rho)
  }
  
  Theta <- function(X){
    Theta_(X, D)
  }
  
  upsilon <- D-1+3
  
  Xi <- matrix(.4, D-1, D-1)
  
  diag(Xi) <- 1
  fit <- basset(Y, X, upsilon, Theta, Gamma, Xi, n_samples = 0, verbose = TRUE, calcGradHess = FALSE)
  new("fido_model", parameters = fit, method = "fido", hyper = list(sigma = sigma, rho = rho))
}

setClass(
  "fido_model",
  slots = c(
    parameters = "ANY",
    method = "character",
    hyper = "list"
  )
)

#' @export
setMethod("predict",  c(object = "fido_model"), fido_predict)

#' @importFrom fido SE
#' @export
Gamma_ <- function(X, sigma, rho) {
  ## here we are assuming RBF kernel over time
  Gamma <- SE(X[1, ,drop = FALSE], sigma = sigma, rho = rho)
  ## next we add in the conditional independence assumption between subjects
  ## following code should zero out any entry where z_i != z_j,
  ## please check
  mask <- as.matrix(dist(factor(X[2, , drop = FALSE])))
  Gamma[mask >= 1] <- 0
  return(Gamma)
}

#' @export
Theta_ <- function(X, D){
  matrix(0, D-1, ncol(X))
}



#' @importFrom tibble as_tibble
#' @importFrom dplyr bind_rows
#' @importFrom purrr reduce
#' @export
fido_data <- function(
    ts_inter, 
    taxonomy=NULL, 
    subject_names=NULL
){
  
  if (is.null(ts_inter@subject_data)) {
    subject_names <- str_c("S", seq_along(ts_inter))
  }
  
  data <- list(
    Y = list(), # Y - taxa x samples - matrix
    samples = list(),
    interventions = list()
  )
  
  for (i in seq_along(ts_inter)) {
    data$Y[[i]] <- values(ts_inter[[i]]) %>%
      data.frame() %>%
      rownames_to_column("name") %>%
      as_tibble()
    
    data$samples[[i]] <- tibble(
      sample = colnames(values(ts_inter[[i]])),
      subject = ts_inter@subject_data$subject[i],
      time = ts_inter[[i]]@time)
    
    data$interventions[[i]] <- interventions(ts_inter[[i]])[,colnames(values(ts_inter[[i]])), drop = FALSE] |> 
      t() |>
      as_tibble() |>
      mutate(sample = colnames(values(ts_inter[[i]]))) |>
      column_to_rownames("sample")
  }
  
  data$Y <- reduce(data$Y, left_join, by = "name") |>
    column_to_rownames("name") |>
    as.matrix()
  
  data$samples <- bind_rows(data$samples)
  
  # intervention can be used as covariates in fido
  data$interventions <- bind_rows(data$interventions) |>
    rownames_to_column("sample")
  
  data$samples <- data$samples |> 
    left_join(ts_inter@subject_data) |>
    left_join(data$interventions)
  
  # # D - numer of taxa
  # D <- nrow(data$Y)
  # # N - number of samples
  # N <- ncol(data$Y)
  
  return(data)
}

