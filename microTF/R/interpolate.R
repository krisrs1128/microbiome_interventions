
approx_mat <- function(time, y_mat, n, method) {
  new_y <- matrix(nrow = nrow(y_mat), ncol = n, dimnames = list(rownames(y_mat), seq_len(n)))
  
  for (i in seq_len(nrow(y_mat))) {
    new_y[i, ] <- approx(
      time, 
      y_mat[i, ], 
      n = n, 
      method = method
    )$y
  }
  
  new_y
}

interpolate_ <- function(ts_inter_single, n, method) {
  time <- ts_inter_single@time
  ts_inter_single@time <- approx(
    seq_along(time), 
    time, 
    n = n, 
    method = method
  )$y
  
  v <- values(ts_inter_single)
  inter <- interventions(ts_inter_single)
  
  values(ts_inter_single) <- approx_mat(time, v, n, method)
  interventions(ts_inter_single) <- approx_mat(time, inter, n, method)
  ts_inter_single
}

interpolate <- function(ts_inter, n = 20, method = "constant") {
  for (i in seq_along(ts_inter)) {
    ts_inter[[i]] <- interpolate_(ts_inter[[i]], n, method)
  }
  ts_inter
}