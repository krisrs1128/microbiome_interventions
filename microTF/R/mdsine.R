

endpoints <- function(z) {
  start_ix <- c()
  end_ix <- c()
  if (z[1] != 0) {
    start_ix <- c(start_ix, 1)
  }
  
  for (i in seq_len(length(z) - 1)) {
    if (z[i] == 0 & z[i + 1] != 0) {
      start_ix <- c(start_ix, i + 1)
    } else if (z[i] != 0 & z[i + 1] == 0) {
      end_ix <- c(end_ix, i + 1)
    }
  }
  
  list(start = start_ix, end = end_ix)
}

perturbation_windows <- function(z, times) {
  perturbations <- list()
  
  k <- 1
  for (i in seq_len(nrow(z))) {
    boundaries <- endpoints(z[i, ])
    perturbations[[i]] <- tibble(
      name = rownames(z)[i],
      start = times[boundaries$start],
      end = times[boundaries$end]
    )
    k <- k + 1

  }
  
  bind_rows(perturbations)
}

check_outputs <- function(data) {
  data  
}

#' @importFrom tibble as_tibble
#' @importFrom tibble bind_rows
#' @export
md_data <- function(ts_inter, taxonomy=NULL, qpcr=NULL, subject_names=NULL) {
  # outputs data.frames that can be written and input to md2.data.parse
  # taxonomy, reads, qpcr, metadata, perturbations
  
  if (is.null(subject_names)) {
    subject_names <- str_c("subject_", seq_along(ts_inter))
  }
  
  data <- list(
    reads = list(),
    perturbations = list()
  )
  for (i in seq_along(ts_inter)) {
    data$reads[[i]] <- values(ts_inter[[i]]) %>%
      data.frame() %>%
      rownames_to_column("name") %>%
      as_tibble()
    data$perturbations[[i]] <- interventions(ts_inter[[i]]) %>%
      perturbation_windows(ts_inter[[i]]@time) %>%
      mutate(subject = subject_names[i])
  }
  
  data <- map(data, bind_rows)
  check_outputs(data)
}

#for (i in seq_along(result)) {
#  rownames(interventions(result[[i]])) <- "1"
#}
#md <- md_data(result)
