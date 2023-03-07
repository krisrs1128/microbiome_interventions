
#' @importFrom dplyr bind_rows
#' @importFrom purrr map
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

#' @export
ts_from_dfs <- function(reads, interventions, metadata) {
  subjects <- unique(metadata$subject)
  series <- list()

  for (i in seq_along(subjects)) {
    sample_ix <- metadata %>%
      filter(subject == subjects[i]) %>%
      pull(time, sample)
    
    x <- t(as.matrix(reads[names(sample_ix), ]))
    z <- t(as.matrix(interventions[names(sample_ix), ]))
    
    series[[i]] <- new(
      "ts_inter_single", 
      values = x[, order(sample_ix)],
      interventions = z[, order(sample_ix)],
      time = sort(sample_ix)
    )
  }
  
  new("ts_inter", series = series)
}