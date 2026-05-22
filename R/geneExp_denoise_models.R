# Get low-dimensional representations of single-cell RNAseq data

#' Run ZINB-WaVE : Zero-Inflated Negative Binomial - Weighted Adaptation for Variance and Effects Generalized linear model (GLM)
#' Helper function (not for users)
#'
#' @param raw_geneExp A data.frame or matrix with cells as rows and
#'   genes as columns. The last column must contain class labels.raw UMI counts.
#' @param min_count Minimum count threshold to retain a gene.
#'   Genes must have counts >= \code{min_count} in at least
#'   \code{min_cells} cells.
#' @param min_cells Minimum number of cells required for a gene
#'   to be retained.
#' @param K Number of latent factors to estimate in ZinbWave
#'   Must be >= 1. Default is 5.
#' @param maxIters Maximum number of iterations for the optimization step (default 100).
#' @param workers Number of parallel workers to use in ZinbWave
#'
#' @return Dataframe of denoised gene expression values
#' @noRd
zinbWave_model <- function(raw_geneExp,
                           min_count            = 1,
                           min_cells            = 1,
                           K                    = 5,
                           maxIters             = 100,       # max optimization iterations
                           workers              = 4) {
  
  epsilon              = 1e6       # regularization (1e6–1e12)
  stop_epsilon         = 1e-4      # convergence criterion
  nb_repeat_initialize = 2        # initialization repeats

  # Split expression matrix and class label 
  class_col_name  <- colnames(raw_geneExp)[ncol(raw_geneExp)]
  raw_geneExp_tmp <- raw_geneExp[, -ncol(raw_geneExp), drop = FALSE]
  
  # genes x cells (ZINB-WaVE convention)
  count_data <- t(raw_geneExp_tmp)
  count_data <- round(as.matrix(count_data))
  
  # Filter low-expressed genes 
  keep_genes      <- rowSums(count_data >= min_count) >= min_cells
  counts_filtered <- count_data[keep_genes, , drop = FALSE]
  
  # Build SingleCellExperiment 
  # zinbwave requires SCE, not plain SummarizedExperiment
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts_filtered)
  )
  
  # Register parallel backend
  # Use SerialParam on Windows (MulticoreParam not supported)
  # Use MulticoreParam on Linux/Mac
  bp <- if (.Platform$OS.type == "windows") {
    BiocParallel::SnowParam(workers = workers)
  } else {
    BiocParallel::MulticoreParam(workers = workers)
  }
  BiocParallel::register(bp)
  
  # Fit ZINB-WaVE model via zinbFit() 
  # zinbFit() returns a zinbModel object — gives access to getLogMu()
  # epsilon controls regularization:
  #   higher = more regularization (1e12 recommended for DE analysis)
  #   lower  = less regularization (1e6 works well for denoising)
  # maxIter controls convergence iterations
  model_zw <- zinbwave::zinbFit(
    sce,
    K                    = K,
    epsilon              = epsilon,
    commondispersion     = FALSE,
    nb.repeat.initialize = nb_repeat_initialize,
    maxiter.optimize     = maxIters,
    stop.epsilon.optimize = stop_epsilon,
    BPPARAM              = bp,
    verbose              = FALSE
  )
  
  # Extract true denoised log(µ) 
  # getLogMu() returns log(µ) from the full fitted NB model
  # This is the TRUE denoised expression — NOT residuals
  # Shape: genes x cells (already on log scale — no further transform needed)
  log_mu <- zinbwave::getLogMu(model_zw)   # cells x genes

  # Assign correct dimnames (cells x genes)
  rownames(log_mu) <- colnames(counts_filtered)   # cell names
  colnames(log_mu) <- rownames(counts_filtered)   # gene names
  
  # Build output: cells x genes + class 
  out_df     <- as.data.frame(log_mu)    # cells x genes
  cell_names <- rownames(out_df)
  
  # Safe class label assignment by cell name
  out_df[[class_col_name]] <- raw_geneExp[cell_names, class_col_name]
  
  return(out_df)
}

#' Denoise scRNA-seq expression using NewWave and log-normalization
#' @param raw_geneExp A data.frame or matrix with cells as rows and
#'   genes as columns. The last column must contain class labels.
#' @param min_count Minimum count threshold to retain a gene.
#'   Genes must have counts >= \code{min_count} in at least
#'   \code{min_cells} cells.
#' @param min_cells Minimum number of cells required for a gene
#'   to be retained.
#' @param K Number of latent factors to estimate in NewWave.
#'   Must be >= 1. Default is 5.
#' @param maxIters Maximum number of iterations for the optimization step (default 100).
#' @param workers Number of parallel workers to use in NewWave. Default is 4.
#'
#' @return A data.frame with Log-normalized gene expression values (cells × genes) and class column
#' @noRd
newWave_model <- function(raw_geneExp,
                          min_count    = 1,
                          min_cells    = 1,
                          K            = 5,
                          maxIters     = 100,    # max optimization iterations
                          workers      = 4) {
  
  epsilon      = 1e6    # regularization (1e6–1e12)
  stop_epsilon = 1e-4   # convergence criterion (default)
  
  # split expression + class
  class_col_name  <- colnames(raw_geneExp)[ncol(raw_geneExp)]
  raw_geneExp_tmp <- raw_geneExp[, -ncol(raw_geneExp), drop = FALSE]
  
  
  # counts: genes x cells
  count_data <- t(raw_geneExp_tmp)
  count_data <- round(as.matrix(count_data))
  
  # filter genes
  keep_genes <- rowSums(count_data >= min_count) >= min_cells
  counts_filtered <- count_data[keep_genes, , drop = FALSE]
  
  # build SCE
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts_filtered)
  )
  
  # run NewWave
  # newFit() returns a newmodel object with all fitted parameters
  # epsilon          : regularization on W, alpha, beta — needs to be large
  #                    enough to compete with likelihood magnitude
  # commondispersion : FALSE = gene-wise dispersion θ_g, better µ per gene
  # maxiter_optimize : hard cap on iterations
  # stop_epsilon     : stop when relative likelihood gain < this (1e-4 default)
  model_fit <- NewWave::newFit(
    sce,
    K                = K,
    epsilon          = epsilon,
    commondispersion = FALSE,
    maxiter_optimize = maxIters,
    stop_epsilon     = stop_epsilon,
    children         = workers,
    verbose          = FALSE
  )
  
  # newLogMu() - directly gives log(mu), genes x cells 
  # newLogMu() = log(mu) = alpha + W * beta^T
  # This is the core denoised signal
  log_mu <- NewWave::newLogMu(model_fit) # cells x genes (n x J)
  
  # assign correct dimnames — rows are cells, cols are genes
  rownames(log_mu) <- colnames(counts_filtered)   # cell names
  colnames(log_mu) <- rownames(counts_filtered)   # gene names
  
  # Build output : cells x genes + class
  out_df     <- as.data.frame(log_mu)
  cell_names <- rownames(out_df)
  
  out_df[[class_col_name]] <- raw_geneExp[cell_names, class_col_name]
  
  return(out_df)

}


