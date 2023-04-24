
#' @importFrom phyloseq phyloseq otu_table sample_data phyloseq_to_deseq2
#' @importFrom DESeq2 sizeFactors estimateSizeFactors
deseq_normalize <- function(reads, metadata) {
  # construct design for DESeq2 object
  if (!("condition" %in% colnames(metadata))) {
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
  reads / size_factors
}

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
    result <- deseq_normalize(reads, metadata)
  } else if (method == "relative_abundance") {
    result <- reads / rowSums(reads)
  } else if (method == "mbImpute") {
    condition <- pull(metadata, condition)
    metadata <- select(metadata, -condition:-sample)
    if (any(apply(metadata, 2, class) %in% c("factor", "character"))) {
      metadata <- model.matrix(~ . , metadata)[, -1]
    }

    if (ncol(metadata) > 0) {
      result <- mbImpute(condition, reads, metadata, ncores = detectCores(), parallel = TRUE, ...)$imp_count_mat_lognorm
    } else {
      result <- mbImpute(condition, reads, ncores = detectCores(), parallel = TRUE, ...)$imp_count_mat_lognorm
    }
  } else if (method == "DESeq2-asinh") {
    result <- asinh(deseq_normalize(reads, metadata))
  }

  result
}