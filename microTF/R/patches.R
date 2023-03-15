
#' @importFrom slider slide
patchify_single <- function(ts_inter, p = 2, q = 3) {
  ix <- seq_len(ncol(ts_inter))
  x_indices <- slide(ix, ~ ., .before = p, .after = -1)
  z_indices <- slide(ix, ~ ., .before = q - 1, .after = 0)
  y_indices <- slide(ix, ~ ., .before = -1, .after = 1)
  
  # initialize result
  k <- 1
  data <- replicate(3, list())
  names(data) <- c("x", "y", "z")
  values_ <- values(ts_inter) 
  interventions_ <- interventions(ts_inter)
  
  # extract x (taxa), z (intervention), and y (future taxa) patches
  for (i in seq_along(x_indices)) {
    if (length(x_indices[[i]]) == p & 
        length(z_indices[[i]]) == q & 
        length(y_indices[[i]]) == 1) {
      data$x[[k]] <- values_[, x_indices[[i]], drop = FALSE]
      data$z[[k]] <- interventions_[, z_indices[[i]], drop = FALSE]
      data$y[[k]] <- values_[, y_indices[[i]], drop = FALSE]
      k <- k + 1
    }
  }
    
  data
}

#' @importFrom stringr str_c
patchify_single_df <- function(ts_inter) {
  data <- patchify_single(ts_inter)
  x <- data$x
  y <- data$y
  z <- data$z
  
  result <- list(
    x = matrix(nrow = length(x), ncol = length(x[[1]]) + length(z[[1]])),
    y = matrix(nrow = length(x), ncol = length(y[[1]]))
  )
  
  for (i in seq_along(x)) {
    xi <- matrix(x[[i]], nrow = 1)
    result$y[i, ] <- matrix(y[[i]], nrow = 1)
    zi <- matrix(z[[i]], nrow = 1)
    result$x[i, ] <- cbind(xi, zi)
  }
  
  colnames(result$x) <- predictor_names(dim(x[[1]]), dim(z[[1]]))
  colnames(result$y) <- str_c("taxon", seq_len(nrow(y[[1]])))
  result
}

#' @importFrom purrr map_dfr
#' @importFrom dplyr select bind_cols
#' @export
patchify_df <- function(ts_inter, p = 2, q = 3) {
  patches <- list()
  for (i in seq_along(ts_inter)) {
    patches[[i]] <- patchify_single_df(ts_inter[[i]])
    sdata <- subject_data(ts_inter)
    if (!is.null(sdata)) {
      patches[[i]]$x <- patches[[i]]$x |>
        bind_cols(select(sdata[i, ], -subject))
    }
  }
  
  x <- map_dfr(patches, ~ as_tibble(.$x))
  y <- map_dfr(patches, ~ as_tibble(.$y))
  list(x = x, y = y)
}

predictor_names <- function(x_dim, z_dim) {
  n1 <- rep(str_c("taxon", seq_len(x_dim[1])), x_dim[2])
  n2 <- rep(str_c("lag", seq(x_dim[2], 1)), each = x_dim[1])
  x_names <- str_c(n1, "_", n2)
  
  n1 <- rep(str_c("intervention", seq_len(z_dim[1])), z_dim[2])
  n2 <- rep(str_c("lag", seq(z_dim[2] - 1, 0)), each = z_dim[1])
  z_names <- str_c(n1, "_", n2)
  
  c(x_names, z_names)
}

lag_from_names <- function(names, group = "taxon") {
  names[str_detect(names, "taxon")] |>
    str_extract("lag[0-9]+") |>
    str_remove("lag") |>
    as.numeric() |>
    max()
}

time_lags <- function(fit) {
  inputs <- fit$var.names
  P <- lag_from_names(inputs, "taxon")
  Q <- lag_from_names(inputs, "intervention") + 1
  c(P, Q)
}

predictors <- function(ts_inter, z_next, lags) {
  x <- values(ts_inter)
  z <- interventions(ts_inter)
  x_prev <- x[, seq(ncol(x) - lags[1] + 1, ncol(x)), drop = FALSE]
  z_prev <- cbind(z[, seq(ncol(z) - lags[2] + 2, ncol(z)), drop = FALSE], z_next)
  
  cbind(matrix(x_prev, nrow = 1), matrix(z_prev, nrow = 1)) |>
    as.data.frame() |>
    set_names(predictor_names(dim(x_prev), dim(z_prev)))
}
