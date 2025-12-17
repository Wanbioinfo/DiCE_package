# Get low-dimensional representations of single-cell RNAseq data

#' Run ZINB-WaVE : Zero-Inflated Negative Binomial - Weighted Adaptation for Variance and Effects Generalized linear model (GLM)
#' Helper function (not for users)
#'
#' @param raw_geneExp Raw UMI count matrix
#'
#' @return Dataframe of denoised gene expression values
#' @noRd
zinbWave_model <- function(raw_geneExp,min_count=1,min_cells=1){
  raw_geneExp_tmp <- raw_geneExp[,-ncol(raw_geneExp)]
  class <- raw_geneExp[,ncol(raw_geneExp)]
  
  count_data <- t(raw_geneExp_tmp)
  count_data <- round(as.matrix(count_data))
  
  # remove genes that have zero counts in all cells 
  keep_genes <- rowSums(count_data >= min_count) >= min_cells
  counts_filtered <- count_data[keep_genes, ]
  
  # create singlecellexperiment object
  # exp_obj <- SingleCellExperiment(assays = list(counts = counts_filtered))
  
  # wrap into SE
  exp_obj <- SummarizedExperiment(assays = list(counts = counts_filtered))
  
  # run ZINB-WaVE
  BiocParallel::register(BiocParallel::MulticoreParam(workers = 4))
  
  zinb <- zinbwave(exp_obj,
                   K=0,
                   epsilon = 1e12,
                   normalizedValues = TRUE,
                   verbose = FALSE,
                   BPPARAM = BiocParallel::MulticoreParam(workers = 4))
  
  # extract denoised gene expression
  denoised <- assay(zinb,"normalizedValues") # gene x cells
  
  denoised_t <- t(denoised)
  
  # Convert to data frames (if not already)
  denoised_t_df <- as.data.frame(denoised_t)
  
  # Extract the class column with matching rownames
  class_col <- colnames(raw_geneExp)[ncol(raw_geneExp)]
  matched_class <- raw_geneExp[rownames(denoised_t_df), class_col]
  
  # Add class as a new column, ensuring alignment
  denoised_t_df$class <- matched_class
  
  return(denoised_t_df)
}

