# Set of correlation calculation methods


#' Run filtered_cor: Custom correlation function excluding shared-zero cells
#' Helper function (not for users)
#'
#' @param x Vector of gene expression of gene1 in all samples
#' @param y Vector of gene expression of gene2 in all samples
#' @param corr_method Correlation method ("pearson" OR "spearman")
#'
#' @return Absolute correlation between two genes 
#' @noRd
filtered_cor <- function(x, y, corr_method) {
  keep <- !(x == 0 & y == 0)
  if (sum(keep) < 3) return(NA)
  corr_matrix <- suppressWarnings(
    abs(cor(class_geneExp, method = corr_method))
  )
  
  # Check if any NA values (often caused by sd = 0 genes) exist
  if (any(is.na(corr_matrix))) {
    warning("Some correlations were undefined due to zero variance; these have been set to 0.",
            call. = FALSE)
  }
  corr_matrix[is.na(corr_matrix)] <- 0
  
  return(corr_matrix)
}

#' Run traditional_corr: Pearson OR spearman corrleation taking all cells
#' Helper function (not for users)
#'
#' @param class_geneExp Dataframe of gene expression of all genes in a speific group of samples
#' @param corr_method Correlation method ("pearson" OR "spearman")
#'
#' @return Absolute correlation matrix between genes 
#' @noRd
traditional_corr <- function(class_geneExp, corr_method){
  corr_matrix <- suppressWarnings(
    abs(cor(class_geneExp, method = corr_method))
  )
  
  # Check if any NA values (often caused by sd = 0 genes) exist
  if (any(is.na(corr_matrix))) {
    warning("Some correlations were undefined due to zero variance; these have been set to 0.",
            call. = FALSE)
  }
  corr_matrix[is.na(corr_matrix)] <- 0
  
  return(corr_matrix)
}

#' Run remove_Zerocells_corr: Pearson OR spearman corrleation removing cells where gene expression is zero for both genes in a gene-pair
#' Helper function (not for users)
#'
#' @param class_geneExp Dataframe of gene expression of all genes in a speific group of samples
#' @param corr_method Correlation method ("pearson" OR "spearman")
#'
#' @return Absolute correlation matrix between genes 
#' @noRd
remove_Zerocells_corr <- function(class_geneExp, corr_method){
  set.seed(42)
  n_cores = 4
  
  gene_names <- colnames(class_geneExp)
  n <- length(gene_names)
  cor_matrix <- matrix(NA, n, n)
  colnames(cor_matrix) <- rownames(cor_matrix) <- gene_names
  
  # All gene index pairs (upper triangle)
  gene_pairs <- combn(n, 2, simplify = FALSE)
  
  # Function to calculate for one gene pair
  compute_pair <- function(pair) {
    i <- pair[1]
    j <- pair[2]
    val <- filtered_cor(class_geneExp[, i], class_geneExp[, j], corr_method)
    list(i = i, j = j, value = val)
  }
  
  # Parallel compute
  results <- mclapply(gene_pairs, compute_pair, mc.cores = n_cores)
  
  # Fill matrix
  for (res in results) {
    cor_matrix[res$i, res$j] <- res$value
    cor_matrix[res$j, res$i] <- res$value
  }
  diag(cor_matrix) <- 1
  return(cor_matrix)
  
}

