
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
  if (is.null(subject_names)) {
    subject_names <- str_c("S", seq_along(ts_inter))
  }

  sample_names <- unlist(map(ts_inter, ~ colnames(values(.))))
  N <- length(sample_names)

  if (is.null(qpcr)) {
    qpcr <- tibble(
      sampleID = sample_names,
      measurement1 = rnorm(sample_names, 1e9, 100),
      measurement2 = rnorm(N, 1e9, 100),
      measurement3 = rnorm(N, 1e9, 100)
    )
  }

  data <- list(
    reads = list(),
    perturbations = list(),
    metadata = list(),
    taxonomy = taxonomy,
    qpcr = qpcr
  )

  for (i in seq_along(ts_inter)) {
    data$reads[[i]] <- values(ts_inter[[i]]) %>%
      data.frame() %>%
      rownames_to_column("name") %>%
      as_tibble()
    data$perturbations[[i]] <- interventions(ts_inter[[i]]) %>%
      perturbation_windows(ts_inter[[i]]@time) %>%
      mutate(subject = subject_names[i])
    data$metadata[[i]] <- tibble(
      subject = subject_names[i],
      sampleID = colnames(values(ts_inter[[i]])),
      time = ts_inter[[i]]@time
    )
  }
  
  data$perturbations <- bind_rows(data$perturbations)
  data$metadata <- bind_rows(data$metadata)
  data$reads <- reduce(data$reads, left_join, by = "name")
  check_outputs(data)
}

#' @importFrom glue glue
#' @export
add_names <- function(ts_inter) {
  for (i in seq_along(result)) {
    inter <- interventions(result[[i]])
   rownames(interventions(result[[i]])) <- str_c("I", seq_along(nrow(inter)))
   colnames(values(result[[i]])) <- glue("S{i}_T{seq_len(ncol(result[[i]]))}")
  }

  result
}
