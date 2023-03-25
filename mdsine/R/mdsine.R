
#' @export
mdsine <- function(ts_inter, taxonomy) {
  data <- md_data(ts_inter, taxonomy)
  do.call(mdsine_, data)
}

#' @importFrom reticulate conda_create conda_list
#' @importFrom dplyr filter pull
#' @importFrom glue glue
#' @export
install_mdsine <- function(envname = "mdsine") {
  conda_create(envname)
  path <- conda_list() |>
    filter(name == envname) |>
    pull(python)
  md_dir <- path(tempdir())
  download.file("https://github.com/gerberlab/MDSINE2/archive/fac880c8c096464cd9a80e72b930cfcf3bd88f21.zip", md_dir / "mdsine2.zip")
  unzip(md_dir / "mdsine2.zip")
  system(glue("{path} -m pip install MDSINE2-fac880c8c096464cd9a80e72b930cfcf3bd88f21/"))
}

#' @importFrom fs path
#' @importFrom reticulate py use_condaenv source_python
#' @export
mdsine_ <- function(taxonomy, reads, qpcr, metadata, perturbations, envname = "mdsine", ...) {
  use_condaenv(envname)
  source_python(system.file("mdsine.py", package = "mdsine"))
  f <- path(tempdir())
  vars <- c("taxonomy", "reads", "qpcr", "metadata", "perturbations")
  paths <- map(vars, ~ f / glue("{.}.tsv"))
  names(paths) <- vars

  map2(vars, paths, ~ write_tsv(get(.x), .y))
  dataset <- py$md_data(paths)
  py$mdsine(dataset, ...)
}

#' @importClassesFrom mbtransfer ts_inter
#' @export
forward_simulate <- function(object, newdata, dt=0.25) {
  fit <- object@parameters
  series <- list()

  for (i in seq_along(newdata)) {
    if (ncol(newdata[[i]]) == ncol(interventions(newdata[[i]]))) {
      series_i <- newdata[[i]]
    } else {

      # get the perturbation and initial conditions for simulation
      split_data <- split_future(newdata[[i]])
      pdata <- perturbation_intervals(newdata[[i]])
      val <- values(split_data$pre)
      x0 <- val[, ncol(val), drop = FALSE]

      # get the forward simulation values at the required times
      n_future <- ncol(split_data$interventions)
      sim <- py$forward_simulate(
        fit, x0, pdata$perturbations, pdata$starts, pdata$ends, dt, n_future
      )
      y_hat <- sim$X[, sim$times %in% seq(0, n_future - 1)]

      # input results to the original ts_inter object
      series_i <- newdata[[i]]
      values(series_i) <- cbind(values(series_i), y_hat)

    }
    series[[i]] <- series_i
  }
  new("ts_inter", series = series)
}

perturbation_intervals <- function(ts_inter) {
  inter <- interventions(ts_inter)
  perturbation_types <- rownames(inter)
  windows <- perturbation_windows(inter, ts_inter@time)

  perturbations <- list()
  for (i in seq_along(perturbation_types)) {
    if (perturbation_types[i] %in% windows$name) {
      perturbations[[i]] <- matrix(1, 1, nrow(ts_inter))
    } else {
      perturbations[[i]] <- matrix(0, 1, nrow(ts_inter))
      windows <- windows %>%
        bind_rows(tibble(name = perturbation_types[i], start = 0, end = 1))
    }
  }

  list(perturbations = perturbations, starts = windows$start, ends = windows$end)
}

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
  
  if (z[length(z)] != 0) {
    end_ix <- c(end_ix, length(z))
  }

  list(start = start_ix, end = end_ix)
}

#' @importFrom dplyr bind_rows
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

dummy_perturbations <- function(perturbations) {
  subjects <- unique(perturbations$subject)
  names <- unique(perturbations$name)
  expand.grid(subject = subjects, name = names) %>%
    mutate(start = 0, end = 1)
}

resolve_perturbations <- function(perturbations, dummies) {
  subjects <- perturbations$subject
  names <- perturbations$name

  for (i in seq_len(nrow(dummies))) {
    exists <- any(
      subjects == dummies$subject[i] &
      names == dummies$name[i]
    )

    if (!exists) {
      perturbations <- rbind(perturbations, dummies[i, ])
    }
  }

  perturbations
}

endpad_interventions <- function(perturbations, metadata) {
  subjects <- unique(perturbations$subject)
  for (s in seq_along(subjects)) {
    final_time <- metadata |>
      filter(subject == subjects[s]) |>
      slice_max(time) |>
      pull(time)
    perturbations <- perturbations |>
      mutate(
        current_subject = subject == subjects[s],
        end = ifelse(current_subject & (end == final_time), end - 2, end)
      ) |>
      select(-current_subject)
  }
  
  perturbations
}

check_outputs <- function(data) {
  dummies <- dummy_perturbations(data$perturbations)
  data$perturbations <- resolve_perturbations(data$perturbations, dummies) |>
    endpad_interventions(data$metadata)
  data
}

#' @importFrom tibble as_tibble
#' @importFrom dplyr bind_rows
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
      measurement1 = rnorm(N, 1e9, 100),
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
      mutate(subject = subject_names[i]) %>%
      select(subject, name, start, end)
    data$metadata[[i]] <- tibble(
      sampleID = colnames(values(ts_inter[[i]])),
      subject = subject_names[i],
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

