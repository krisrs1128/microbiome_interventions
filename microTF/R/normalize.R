
#' @importFrom phyloseq phyloseq otu_table sample_data phyloseq_to_deseq2
#' @importFrom DESeq2 sizeFactors estimateSizeFactors
#' @importFrom mbImpute mbImpute
#' @export
normalize <-  function(reads, method = "none",  metadata = NULL, ...) {
  if (method == "none") {
    result <- x
  } else if (method == "DESeq2") {

    # construct design for DESeq2 object
    if (is.null(metadata)) {
      metadata <- data.frame(dummy = rep(1, nrow(reads)))
      rownames(metadata) <- rownames(reads)
      fmla <- formula(~ 1)
    } else {
      fmla <- formula(~ condition)
    }
    
    # construct phyloseq and perform normalization
    ps <- phyloseq(
      otu_table(reads, taxa_are_rows = FALSE),
      sample_data(metadata)
    )
    dds <- phyloseq_to_deseq2(ps, fmla)
    size_factors <- sizeFactors(estimateSizeFactors(dds, "poscounts"))
    result <- reads / size_factors

  } else if (method == "relative_abundance") {
    result <- reads / rowSums(result)
  } else if (method == "mbImpute") {
    condition <- pull(metadata, condition)
    metadata <- select(metadata, -condition)
    if (ncol(metadata) > 0) {
      result <- mbImpute(condition, reads, metadata, ...)$imp_count_mat_lognorm
    } else {
      result <- mbImpute(condition, reads, ...)$imp_count_mat_lognorm
    }
  }
  
  result
}


