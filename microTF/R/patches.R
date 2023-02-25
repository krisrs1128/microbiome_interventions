
#' @importFrom slider slide
patchify_single <- function(ts_inter, p = 4) {
  ix <- seq_len(ncol(ts_inter))
  x_indices <- slide(ix, ~ ., .before = p, .after = -1)
  y_indices <- slide(ix, ~ ., .before = -1, .after = 1)
  
  k <- 1
  x <- list()
  y <- list()
  for (i in seq_along(x_indices)) {
    if (length(x_indices[[i]]) == p & length(y_indices[[i]]) == 1) {
      x[[k]] <- ts_inter[, x_indices[[i]]]
      y[[k]] <- ts_inter[, y_indices[[i]]]
      k <- k + 1
    }
  }
    
  list(x = x, y = y)
}

#' @importFrom purrr map2
#' @export
patchify <- function(ts_inter, p = 4) {
  x <- list()
  y <- list()
  for (i in seq_along(ts_inter)) {
    patches <- patchify_single(ts_inter[[i]], p)
    x[[i]] <- patches$x
    y[[i]] <- patches$y
  }
  x <- new("ts_inter", series = unlist(x))
  y <- new("ts_inter", series = unlist(y))
  
  source <- map2(seq_along(patches), patches, ~ rep(.x, length(.y$y)))
  list(x = x, y = y, source = unlist(source))
}

pivot_elem <- function(x, name = "taxon") {
  result <- matrix(x, nrow = 1)
  n1 <- rep(str_c(name, seq_len(nrow(x))), ncol(x))
  n2 <- rep(str_c("time", seq_len(ncol(x))), each = nrow(x))
  colnames(result) <- str_c(n1, "_", n2)
  result 
}

pivot_ts_<- function(ts_inter, all = TRUE) {
  x <- pivot_elem(values(ts_inter))
  if (all) {
    z <- pivot_elem(interventions(ts_inter), name = "intervention")
    x <- do.call(cbind, list(x, z))
  }
  
  as_tibble(x)
}

pivot_ts <- function(ts_inter, ...) {
  result <- list()
  for (i in seq_along(ts_inter)) {
    result[[i]] <- pivot_ts_(ts_inter[[i]], ...)
  }
  
  bind_rows(result, .id = "patch")
}
