
#' @importFrom phyloseq phyloseq otu_table sample_data phyloseq_to_deseq2
#' @importFrom DESeq2 sizeFactors estimateSizeFactors
normalize <-  function(reads, method = "none") {
  if (method == "none") {
    result <- x
  } else if (method == "DESeq2") {
    sdata <- data.frame(dummy = rep(1, nrow(reads)))
    rownames(sdata) <- rownames(reads)
    
    ps <- phyloseq(
      otu_table(reads, taxa_are_rows = FALSE),
      sample_data(sdata)
      )
    dds <- phyloseq_to_deseq2(ps, ~ 1)
    size_factors <- sizeFactors(estimateSizeFactors(dds, "poscounts"))
    result <- reads / size_factors
  }
  
  result
}