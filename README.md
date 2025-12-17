# DiCE
## Version: V1.1.3

## Disease Biomarker Identification with scRNA-seq Data

### Required R Libraries

```r
dplyr,tibble,FSelectorRcpp,igraph,data.table,NetWeaver,praznik,zinbwave,SingleCellExperiment,
reticulate,stats,utils,openxlsx,BiocParallel,SummarizedExperiment,purrr,parallel,annotate,
org.Hs.eg.db,AnnotationDbi,org.Mm.eg.db,readxl
```
---


### Installation

```r
install.packages("path/to/DiCE_V1.1.3.tar.gz", repos = NULL, type = "source")

```
---

### Step 1: Prepare Input Files for DiCE

If you are running DiCE on **bulk RNA-seq** data, you need the following two input files:

1. **"dge_data"** – Rds/xlsx/csv/tsv file containing the differential gene expression analysis results of **protein coding genes**.
  - The columns must include: **"Gene.Symbol", "logFC", "P.Value", "adj.P.Val"**. 
  - Note: Column names in your data must exactly match the above.

2. **"logNorm_geneExp"** - Rds/xlsx/csv/tsv file containing the normalized (logCPM) gene expression data of **protein coding genes**.
  - Rows represent samples, and columns represent genes.
  - There should be one additional column indicating the sample group/class.
    - Example: If a sample is from a tumor tissue, the Group column should contain “Tumor”.


If you are running DiCE on **scRNA-seq** data, you need the following three input files:


1. **"dge_data"** – Rds/xlsx/csv/tsv file containing the differential gene expression analysis results of **protein coding genes**.
  - The columns must include: **"Gene.Symbol", "logFC", "P.Value", "adj.P.Val"**. 
  - Note: Column names in your data must exactly match the above.

2. **"logNorm_geneExp"** – Rds/xlsx/csv/tsv file containing the normalized (logCPM) gene expression data of **protein coding genes**.
  - Rows represent samples, and columns represent genes.
  - There should be one additional column indicating the sample group/class.
    - Example: If a sample is from a tumor tissue, the Group column should contain “Tumor”.

3. **"raw_gene_exp"** – Rds/xlsx/csv/tsv file containing the raw UMI counts of the **protein coding genes** .
  - Rows represent cells, and columns represent genes.
  - There should be one additional column indicating the cell group/class.
    - Example: If a sample is from a tumor tissue, the Group column should contain “Tumor”.

If you plan to use **ZINB-WaVE-denoised gene expression data** for Information Gain (IG) calculation, use:

- `zinbWaVE_denoising()` Function

Example:

```r
zinb_denoised <- zinbWaVE_denoising(
	dge_file = "path/to/dge_data.rds",
	rawGeneExp_file = "path/to/raw_gene_exp.rds",
	loose_criteria = "adj.P.Val",
	loose_cutoff = 0.05,
	logFC_cutoff = 1
)
```

Note: Use the same filtering option and cutoff value as in **Phase 1**.

---

### Step 2: Run `DiCE.R` to Identify DiCE Genes

Configure the parameters based on your application:

- **data_type** : Sequencing data type ("bulkRNA-seq" OR "scRNA-seq"). Default is "bulkRNA-seq".
- **species** : "human" or "mouse". Default is "human".
- **dge_file_path** : File path to the differential gene expression file. The columns must include: "Gene.Symbol", "logFC", "P.Value", "adj.P.Val".
- **normGeneExp_file_path** : File path for the differential gene expression analysis '.Rds' file. (cells/samples x genes + label column)
- **rawGeneExp_file_path** : File path for the raw UMI counts file (Needs only when the data_type = "scRNA-seq"). (cells/samples x genes + label column)
- **treatment** : Label of treatment samples (eg: Tumor)
- **control** : Label of control samples (eg: Normal)
- **loose_criteria** : Statistical significance metric used for initial gene filtering ("P.Value" or "adj.P.Val"). Default is "adj.P.Val"
- **loose_cutoff** : Numeric threshold for filtering based on 'loose_criteria'. Genes with values ≤ this cutoff are retained. Default is 0.05.
- **logFC_cutoff** : Minimum absolute log2 fold change threshold for retaining genes. Default is 0.
- **is_wIG_needed** : Whether to compute weighted IG. Options: "yes"/TRUE or "no"/FALSE. Default is "no"/FALSE.
- **B** : Number of bootstrap resamples used in weighted IG calculation. Default is 300.
- **ig_cutoff** - Method for selecting IG-filtered genes. Options include: 
	- "all_mean" : Retain all genes with IG greater than the mean IG computed across all genes (including zeros). Default.
	- "all_median" : Retain all genes with IG greater than the median IG computed across all genes (including zeros)
	- "nonzero_mean" : Compute the mean IG only among genes with IG > 0, and retain genes whose IG exceeds this non-zero mean threshold
	- "nonzero_median" : Compute the median IG only among genes with IG > 0, and retain genes whose IG exceeds this non-zero median threshold
	- "all_nonzero" : Retain all genes with IG > 0
- **norm_type** : Type of the gene expression normalization ("logNorm" OR "ZINBWaVE_denoised"). Default is "logNorm".
- **corr_mode** : Mode for computing gene–gene correlation. Determines how gene expression data are pre-processed and how correlation is computed. Options:
	- "directCorr" : Use raw or normalized expression values without dropping zero-expression cells. Computes correlation directly. (Default)
	- "remove_Zerocells" : Exclude cell pairs where both genes have zero expression before computing correlation. Useful for sparse data. (For scRNAseq data).
	- "ZINB-WaVE" : Apply ZINB-WaVE denoising to model zero inflation and overdispersion before computing correlation. (For scRNAseq data). 
- **corr_method** : "pearson" OR "spearman". Default is "pearson".
- **centrality_list** : Character vector of centrality metrics to compute. Valid options include: "betweenness", "eigen vector", "pagerank", "closeness", "harmonic", "authority", "strength". 
- **min_passCount** : Minimum number of centrality metrics a gene must pass to be retained. Default is the length of **centrality_list**.
- **cutoff** : Method for selecting final candidate genes based on centrality metrics. ("mean" OR "median" OR "topK%" ("top10%", "top25%", etc.)). Default is "mean".



Exmaple: Run DiCE pipeline on bulk RNA-seq data

```r
dice_results_df <- perform_DiCE(
  data_type = "bulkRNA-seq",
  dge_file_path = "path/to/dge_results.rds",
  normGeneExp_file_path = "path/to/logNorm_geneExp.rds",
  treatment = "Tumor",
  control = "Normal",
  loose_criteria = "adj.P.Val",
  loose_cutoff = 0.05,
  logFC_cutoff = 1,
  species = "human",
  norm_type = "logNorm",
  corr_mode = "directCorr",
  corr_method = "pearson",
  cutoff = "mean"
)

# View DiCE genes
head(dice_results_df)
```

Exmaple: Run DiCE pipeline on scRNA-seq data

```r
dice_results_df <- perform_DiCE(
  data_type = "bulkRNA-seq",
  dge_file_path = "path/to/dge_results.rds",
  normGeneExp_file_path = "path/to/logNorm_geneExp.rds",
  rawGeneExp_file_path = "path/to/raw_geneExp.rds",
  treatment = "Tumor",
  control = "Normal",
  loose_criteria = "adj.P.Val",
  loose_cutoff = 0.05,
  logFC_cutoff = 1,
  species = "human",
  norm_type = "logNorm",
  corr_mode = "ZINB-WaVE",
  corr_method = "spearman",
  cutoff = "mean"
)

# View DiCE genes
head(dice_results_df)
```

### Output

- **dice_results_df**: Data frame of all genes in Phase1,  Phase2, Phase3, and DiCE genes with their final ensemble ranking.


### Step 3 (Optional): Run `detect_DiCE_PPI_unweightedModules.R` to Identify PPI network communities among DiCE Genes

Create a dataframe with only the DiCE genes and input it to the `detect_DiCE_PPI_unweightedModules.R` function.

This function outputs:
- summary_df: Summary of number of modules and modularity
- membership_df: Module assignments and within-module degrees
- edges_by_module : List of intra-module interaction tables

Example for find communities among the DiCE genes"

```r
# Run module detection on DiCE genes
modules <- detect_DiCE_PPI_unweightedModules(
 dice_genes_df = dice_genes_df,
 species = "human",
 seed = 123
)

# View module membership
head(modules$membership_df)

```

### Last update - 12/05/2025

### Updates

- 12/05/2025​ - DiCE now supports input files in multiple formats, including .rds, .xlsx, .tsv, and .csv, allowing users to load data seamlessly regardless of file type. 
- 13/11/2025 - Added different IG cutoffs (all_mean, all_median, nonzero_mean, nonzero_median, all_nonzero)
- 13/11/2025 - Added a function at the beginning of the DiCE pipeline to retain only protein-coding genes before downstream analysis.
- 03/11/2025 - Introduced other network centralities (Authority, Strength, Closeness, pagerank, Harmonic). User can select the centralities they want.
- 21/10/2025 - Included function to find the network modules/communities among the DiCE genes with unweighted PPI and Louvain algorithm. 
- 15/10/2025 - Introduced class weighted Information Gain with Monte carlo approximation.
- 12/10/2025 - Handled errors/warnings due to column names mismatches in the input data.
- 08/10/2025 - Added IG information for Phase 1 genes, updated the script to return only the dataframe containing all Phase 1, Phase 2, Phase 3, and DiCE genes with their final ensemble ranking, and handled variations in column names provided by the user.
- 20/09/2025 - Corrected the final ensemble ranking calculation of all genes.
- 17/09/2025 - Changed the PPI construction method to use downloaded StringDB files instead of the StringDB R package and creates the final gene ranking table. 
- 10/09/2025 - Corrected eigen vector centrality calculation using correlation coefficient as the edge weights.
- 10/10/2025 - Forced correlations which were undefined due to zero variance to be zero.
- 08/10/2025 - Added IG information for Phase 1 genes, updated the script to return only the dataframe containing all Phase 1, Phase 2, Phase 3, and DiCE genes with their final ensemble ranking, and handled variations in column names provided by the user.
- 20/09/2025 - Corrected the final ensemble ranking calculation of all genes.
- 17/09/2025 - Changed the PPI construction method to use downloaded StringDB files instead of the StringDB R package and creates the final gene ranking table. 
- 10/09/2025 - Corrected eigen vector centrality calculation using correlation coefficient as the edge weights.
- 08/10/2025 - Added IG information of Phase1 genes output and changed to return only the dataframe with all genes in Phase1,  Phase2, Phase3, and DiCE genes with their final ensemble ranking.
- 20/09/2025 - Corrected the final ensemble ranking calculation of all genes.
- 17/09/2025 - Changed the PPI construction method to use downloaded StringDB files instead of the StringDB R package and creates the final gene ranking table. 
- 10/09/2025 - Corrected eigen vector centrality calculation using correlation coefficient as the edge weights.

			
