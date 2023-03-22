
#' @importFrom dplyr pull filter
#' @export
ts_from_dfs <- function(reads, interventions, metadata, subject_data = NULL) {
  subjects <- unique(metadata$subject)
  series <- list()

  # fill in each subject's intervention and abundance data
  for (i in seq_along(subjects)) {
    sample_ix <- metadata |>
      filter(subject == subjects[i]) |>
      pull(time, sample)
    
    x <- t(as.matrix(reads[names(sample_ix), ]))
    z <- t(as.matrix(interventions[names(sample_ix), ]))
    
    series[[i]] <- new(
      "ts_inter_single", 
      values = x[, order(sample_ix)],
      interventions = z[, order(sample_ix), drop=FALSE],
      time = sort(sample_ix)
    )
  }

  # ensure subject and metadata order agree
  if (!is.null(subject_data)) {
    names(series) <- subjects
    subject_data <- metadata |>
      distinct(subject) |>
      left_join(subject_data)
  }
  
  new("ts_inter", series = series, subject_data = subject_data)
}