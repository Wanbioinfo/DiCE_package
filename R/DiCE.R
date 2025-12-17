#' Run DiCE Pipeline
#'
#' This function runs the DiCE pipeline for gene biomarker identification from transcriptomic data using differential expression, feature selection, and network-based centrality metrics.
#'
#' @param data_type Sequencing data type ("bulkRNA-seq" OR "scRNA-seq"). Default "bulkRNA-seq".
#' @param species "human" or "mouse". Default "human".
#' @param dge_file_path File path to the differential gene expression '.Rds' file. The columns must include: "Gene.Name", "logFC", "P.Value", "adj.P.Val".
#' @param normGeneExp_file_path File path for the normalized gene expression '.Rds' file (cells/samples x genes + label column).
#' @param rawGeneExp_file_path File path for the raw UMI counts '.Rds' file (only needed when data_type = "scRNA-seq"). (cells x genes + label column)
#' @param treatment Label of treatment samples (e.g., "Tumor").
#' @param control Label of control samples (e.g., "Normal").
#' @param loose_criteria Statistical significance metric used for initial gene filtering ("P.Value" or "adj.P.Val"). Default "adj.P.Val".
#' @param loose_cutoff Numeric threshold for filtering based on 'loose_criteria'. Genes with values ≤ this cutoff are retained. Default 0.05.
#' @param logFC_cutoff Minimum absolute log2 fold change threshold for retaining genes. Default 0.
#' @param is_wIG_needed Whether to compute weighted IG. Options: "yes"/TRUE or "no"/FALSE. Default is "no"/FALSE.
#' @param B Number of bootstrap resamples used in weighted IG calculation. Default is 300.
#' @param ig_cutoff Method for selecting IG-filtered genes. Options include:
#'   \describe{
#'     \item{"all_mean"}{Retain all genes with IG greater than the mean IG computed across all genes (including zeros). This is the standard mean-based filter. (Default)}
#'     \item{"all_median"}{Retain all genes with IG greater than the median IG computed across all genes (including zeros). More robust than the mean for extremely skewed IG distributions.}
#'     \item{"nonzero_mean"}{Compute the mean IG only among genes with IG > 0, and retain genes whose IG exceeds this non-zero mean threshold. Useful when IG has many exact zeros.}
#'     \item{nonzero_median}{Compute the median IG only among genes with IG > 0, and retain genes whose IG exceeds this non-zero median threshold. Ideal when IG is sparse and strongly right-skewed.}
#'     \item{all_nonzero}{Retain all genes with IG > 0. This keeps every gene carrying information with respect to the class label and excludes only those with zero information gain.}
#'   }
#' @param norm_type Type of the gene expression normalization ("logNorm" OR "ZINBWaVE_denoised"). Default "logNorm".
#' @param corr_mode Mode for computing gene–gene correlation. Options are:
#' \describe{
#'   \item{"directCorr"}{Use raw or normalized expression values without dropping zero-expression cells. (Default)}
#'   \item{"remove_Zerocells"}{Exclude cell pairs where both genes have zero expression before computing correlation.}
#'   \item{"ZINB-WaVE"}{Apply ZINB-WaVE denoising before computing correlation.}
#' }
#' @param corr_method Correlation method ("pearson" OR "spearman"). Default "pearson".
#' @param centrality_list Character vector of centrality metrics to compute.  
#'   Valid options include: "betweenness", "eigen vector", "pagerank", "closeness", "harmonic", "authority", "strength". 
#'   (Names are normalized internally.)
#' @param min_passCount Minimum number of centrality metrics a gene must pass to be retained. Default is the length of \code{centrality_list}.
#' @param cutoff Method for selecting final candidate genes based on centrality metrics ("mean", "median", or "topK\%" such as "top10\%", "top25\%"). Default "mean".
#' 
#' @return A data frame containing the final DiCE-ranked genes with:
#'   \itemize{
#'     \item Centrality values across metrics
#'     \item Ensemble ranking score
#'     \item IG or weighted IG
#'     \item Differential expression statistics
#'   }
#'
#' @importFrom dplyr mutate filter as_tibble bind_rows select
#' @importFrom FSelectorRcpp information_gain
#' @importFrom igraph graph_from_data_frame set_edge_attr E E<- V V<-
#' @importFrom igraph eigen_centrality betweenness
#' @importFrom utils stack
#' @importFrom parallel mclapply
#' @importFrom NetWeaver ensemble_rank
#' @importFrom zinbwave zinbwave
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom stats cor quantile median
#' @importFrom BiocParallel register MulticoreParam
#' @importFrom SummarizedExperiment SummarizedExperiment assay
#' @importFrom stats as.formula setNames
#' @importFrom utils combn read.csv read.delim read.table
#' @importFrom readxl read_excel
#' @importFrom data.table fread
#' @importFrom igraph authority_score closeness harmonic_centrality
#' @importFrom igraph page_rank strength simplify cluster_louvain
#' @importFrom igraph membership sizes modularity induced_subgraph degree
#' @importFrom AnnotationDbi keys
#'
#' @examples
#' \dontrun{
#' dice_results_df <- perform_DiCE(
#'   data_type = "bulkRNA-seq",
#'   dge_file_path = "path/to/dge_results.rds",
#'   normGeneExp_file_path = "path/to/logNorm_geneExp.rds",
#'   treatment = "Tumor",
#'   control = "Normal",
#'   loose_criteria = "adj.P.Val",
#'   loose_cutoff = 0.05,
#'   logFC_cutoff = 1,
#'   species = "human",
#'   is_wIG_needed = "yes",
#'   B = 200,
#'   ig_cutoff = "nonzero_mean",
#'   norm_type = "logNorm",
#'   corr_mode = "directCorr",
#'   corr_method = "pearson",
#'   cutoff = "mean"
#' )
#' head(dice_results_df)
#' }
#'
#' @export
perform_DiCE <- function(
    data_type = "bulkRNA-seq",
    species = "human",
    dge_file_path,
    normGeneExp_file_path,
    rawGeneExp_file_path,
    treatment,
    control,
    loose_criteria = "adj.P.Val",
    loose_cutoff = 0.05,
    logFC_cutoff = 0,
    is_wIG_needed = "no",
    B = 300,
    ig_cutoff = "all_mean",
    norm_type = "logNorm",
    corr_mode = "directCorr",
    corr_method = "pearson",
    centrality_list = c("betweenness", "eigen vector"),
    min_passCount = length(centrality_list),
    cutoff = "mean"
    ) 
{
  
  if (data_type == "bulkRNA-seq"){
    norm_type <- "logNorm"
    corr_mode <- "directCorr"
    corr_method <- "pearson"
  }

  # String db downloaded files
  if(species == "human"){
    string_protInfo_file <- system.file("extdata/stringDB_v12/human/9606.protein.info.v12.0.txt", package = "DiCE")
    string_ppi_file <- system.file("extdata/stringDB_v12/human/9606.protein.links.v12.0.txt", package = "DiCE")
    taxonID <- 9606
  }else if(species == "mouse"){
    string_protInfo_file <- system.file("extdata/stringDB_v12/mouse/10090.protein.info.v12.0.txt", package = "DiCE")
    string_ppi_file <- system.file("extdata/stringDB_v12/mouse/10090.protein.links.v12.0.txt", package = "DiCE")
    taxonID = 10090
  }

  ############################### Load data #############################################
  
  message("Loading input data...")
  
  # Load DGE data and gene expression data
  dge_data <- read_any(dge_file_path)
  normGeneExp_data <- read_any(normGeneExp_file_path)
  
  if (data_type == "scRNA-seq") {
    rawGeneExp_data <- read_any(rawGeneExp_file_path)
  }
  
  
  # Change column names to DiCE col names
  dge_data <- normalize_dge_cols(dge_data)
  
  ###################### Keep only protein coding genes #################################
  
  message("Starting protein-coding gene filtering...")
  
  dge_data <- keep_protCoding_dgeData(taxonID, dge_data)
  
  prot_genes <- dge_data$Gene.Name
  
  normGeneExp_data <- keep_protCoding_expData(prot_genes, normGeneExp_data)
  
  if (data_type == "scRNA-seq") {
    rawGeneExp_data <- keep_protCoding_expData(prot_genes, rawGeneExp_data)
  }
  
  ############################### Phase 1 ###############################################
  
  # Construction of a candidate gene pool by Differential Gene Expression
  # analysis with a loose cutoff
  
  phase1_res <- run_phase1(dge_data, loose_criteria, loose_cutoff, logFC_cutoff)
  
  phase1_criteria <- paste0(loose_criteria,
                            "<",
                            loose_cutoff,
                            " & ",
                            "|logFC|>",
                            logFC_cutoff)
  
  message(paste0("#Genes in Phase 1 (", phase1_criteria, ") = ", nrow(phase1_res)))
  
  
  ############################### Phase 2 ###############################################
  
  # Selection of the top discriminative genes using the Information Gain
  # filter approach and construct the PPI interactions.
  
  phase2_res <- run_phase2(normGeneExp_data, phase1_res, is_wIG_needed, ig_cutoff, B)
  
  phase2_genes_df <- phase2_res$phase2_genes_df
  ig_df <- phase2_res$infoGain_df
  
  phase1_res <- merge(phase1_res, ig_df, by = "Gene.Name", all.x = TRUE) 
  
  if (isTRUE(is_wIG_needed) || tolower(is_wIG_needed) %in% c("yes", "true", "y", "1")){
    ig_term <- "WeightedIG"
  }else{
    ig_term <- "IG"
  }
  
  if (ig_cutoff == "all_mean"){
    phase2_criteria <- paste0("Mean(",ig_term,")")
    
  }else if (ig_cutoff == "all_median"){
    phase2_criteria <- paste0("Median(",ig_term,")")
    
  }else if (ig_cutoff == "nonzero_mean"){
    phase2_criteria <- paste0("Mean(",ig_term,">0)")
    
  }else if (ig_cutoff == "nonzero_median"){
    phase2_criteria <- paste0("Median(",ig_term,">0)")
  
  }else if (ig_cutoff == "all_nonzero"){
    phase2_criteria <- paste0("All non-zero IG")
  }else{
    stop("Invalid cutoff for IG!. Please select from 'all_mean', 'all_median', 'nonzero_mean', 'nonzero_median', and 'all_nonzero'.")
  }
  
  message(paste0("#Genes in Phase 2 (",phase2_criteria,") = " , nrow(phase2_genes_df)))
  
  
  
  ############################### Phase 3 ###############################################
  
  # Making weighted PPI using (1-|C.C|) for each phenotype
  
  # Extract interactions
  ppi_results_phase3 <- extract_PPI(string_protInfo_file, 
                                    string_ppi_file, 
                                    phase2_genes_df)
  pp_interactions <- ppi_results_phase3$interactions
  mapped_proteins <- ppi_results_phase3$mapped_proteins
  
  filtered_normGeneExp_data <- getNetworkGene_expression(normGeneExp_data,
                                                         mapped_proteins)
  
  if (data_type == "scRNA-seq"){
    filtered_raw_geneExp <- getNetworkGene_expression(rawGeneExp_data,
                                                      mapped_proteins)
  }
  
  if (corr_mode == "ZINB-WaVE"){
    print("Raw gene expression - ZINB-WaVE")
    filtered_denoised_geneExp <- zinbWave_model(filtered_raw_geneExp)
    geneExp_to_corr <- filtered_denoised_geneExp
  }else if ((corr_mode == "directCorr") | (corr_method == "remove_Zerocells")){
    print("Normalized gene expressions")
    geneExp_to_corr <- filtered_normGeneExp_data
  }else{
    stop("Invalid correlation mode. Use 'directCorr', 'remove_Zerocells', 'ZINB-WaVE'")
  }
  
  # Run Phase 3
  phase3_res <- run_phase3(pp_interactions, geneExp_to_corr, treatment, control, 
                           corr_method, corr_mode, centrality_list)
  
  phase3_centralities_df <- phase3_res$phase3_centralities_df
  interactions_withWeights_df <- phase3_res$interactions_withWeights_df
  
  phase3_centralities_df$Phase <- "III"
  
  message(paste("#Genes in Phase 3 (StringDB PPI) = ", nrow(phase3_centralities_df)))
  
  phase3_centralities_out_df <- merge(phase3_centralities_df, mapped_proteins, 
                               by.x = "Gene.Name", by.y = "Gene.Name")
  colnames(phase3_centralities_out_df)[1] <- "Gene.Name"
  
  # Rearrange the columns
  base_cols <- c("Gene.Name", "STRING_id", "logFC", "adj.P.Val", "P.Value", if ("IG" %in% colnames(phase3_centralities_out_df)) "IG")
  cent_cols <- colnames(phase3_centralities_out_df)[2:((length(centrality_list)*2)+1)]
  tail_cols <- "Phase"
  final_cols <- c(base_cols, cent_cols, tail_cols)
  
  phase3_centralities_out_df <- phase3_centralities_out_df[, final_cols, drop = FALSE]
  
  ############################### Phase 4 ###############################################
  
  # Ensemble ranking
  
  phase4_res <- run_phase4(phase3_centralities_out_df,cutoff, centrality_list, min_passCount)
  
  DiCE_genes_df <- phase4_res$DiCE_genes_df
  phase4_centralities_df <- phase4_res$phase4_centralities_df
  
  message(paste0("#Genes in Phase 4 (",cutoff,") = ", nrow(DiCE_genes_df)))
  
  
  ############################### Final ensemble ranking ################################
  
  
  # Create final ensembl ranking
  
  dice_results_df <- createFinalRanking(dge_data, phase1_res, 
                                      phase2_genes_df, 
                                      phase4_centralities_df,
                                      centrality_list)
  
  ############################### Return output #####################################
  
  return(dice_results_df)

}