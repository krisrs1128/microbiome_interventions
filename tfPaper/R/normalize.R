
#' @importFrom phyloseq phyloseq otu_table sample_data phyloseq_to_deseq2
#' @importFrom DESeq2 sizeFactors estimateSizeFactors
#' @importFrom mbImpute mbImpute
#' @importFrom parallel detectCores
#' @export
normalize <-  function(reads, method = "none", metadata = NULL, ...) {
  if (is.null(metadata)) {
    metadata <- data.frame(dummy = rep(1, nrow(reads)), sample = rownames(reads))
  }

  if (method == "none") {
    result <- reads
  } else if (method == "DESeq2") {

    # construct design for DESeq2 object
    if (is.null(metadata)) {
      fmla <- formula(~ 1)
    } else {
      fmla <- formula(~ condition)
    }

    # construct phyloseq and perform normalization
    metadata <- metadata |>
      column_to_rownames("sample") |>
      sample_data()
    ps <- phyloseq(
      otu_table(reads, taxa_are_rows = FALSE),
      metadata
    )
    
    dds <- phyloseq_to_deseq2(ps, fmla)
    size_factors <- sizeFactors(estimateSizeFactors(dds, "poscounts"))
    result <- reads / size_factors

  } else if (method == "relative_abundance") {
    result <- reads / rowSums(reads)
  } else if (method == "mbImpute") {
    condition <- pull(metadata, condition)
    metadata <- select(metadata, -condition:-sample)
    if (any(apply(metadata, 2, type) %in% c("factor", "character"))) {
      metadata <- model.matrix(~ -1 + . , metadata)
    }

    n_cores <- detectCores()
    if (ncol(metadata) > 0) {
      result <- mbImpute(condition, reads, metadata, ncores = n_cores, parallel = TRUE, ...)$imp_count_mat_lognorm
    } else {
      result <- mbImpute(condition, reads, ncores = n_cores, parallel = TRUE, ...)$imp_count_mat_lognorm
    }
  }

  result
}
