
pivot_slot <- function(ts_inter, f, name = "taxon") {
  result <- map(ts_inter, ~ f(.)) %>%
    map(
      ~ as_tibble(.) %>% 
          mutate(id = row_number()) %>%
          pivot_longer(-id, names_to = "time")
      )
  
  for (i in seq_along(result)) {
    result[[i]] <- result[[i]] %>%
      mutate(time = ts_inter[[i]]@time[as.integer(str_remove(time, "V"))])
  }
  
  result %>%
    bind_rows(.id = "subject") %>%
    set_names(c("subject", name, "time", "value"))
}

#' @export
tidy_ts <- function(ts_inter) {
  z <- pivot_slot(ts_inter, interventions, "intervention")
  x <- pivot_slot(ts_inter, values)
  list(z = z, x = x)
}