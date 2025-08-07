# Identification of hub genes and pathways in Cervical Cancer
This repository contains the code and analysis to reproduce the results from the study:  
*Bioinformatics analysis of differentially expressed genes and pathways in the development of cervical cancer* by Baojie Wu* and Shuyi Xi
 
#### Authors: Poojitha Kolli, Dillen Wischmeier, Manjinder Kaur  

## Description:
   This repository contains the full code, data processing steps, and analysis outputs used to reproduce the results of a published study on differential gene expression meta-analysis across three GEO datasets: GSE63514, GSE62417, and GSE138080. The workflow includes identification of differentially expressed genes (DEGs), functional enrichment, gene set enrichment analysis (GSEA), and protein-protein interaction (PPI) network construction and hub gene identification, following the procedures described in the original publication.

## Requirements:
#### Programming Languages:  
- R v4.3.3
- Python v.3.11.11  
#### R packages (DEG analysis): 
`Biobase_2.62.0`  
`GEOquery_2.70.0`  
`BiocGenerics_0.48.1`  
`limma_3.58.1`  
`umap_0.2.10.0`   
#### R packages (GSEA analysis):  
`dplyr 1.1.4`  
`janitor 2.2.1`  
#### Python packages:
`pandas==1.5.3`  
`venn==0.1.3`  
`matplotlib==3.7.1`  
 
#### External Tools:  
GEO2R tool - https://www.ncbi.nlm.nih.gov/geo/geo2r/  
DAVID v2024q4 (web tool) - https://david.ncifcrf.gov/  
GSEA v4.4.0 desktop app - http://www.broadinstitute.org/gsea/msigdb/downloads.jsp  
STRING v12.0 (web tool) - https://string-db.org/  
Cytoscape v3.10.3 (desktop app with CytoHubba and MCODE plugins) - https://cytoscape.org/download.html  
 
## Execution:
### 1. Create Output folder
```bash
   mkdir -p Output
```

### 2. DEG Analysis
**Scripts** (adapted from GEO2R workflow to automate volcano plot and DEG list generation)**:**  
`GSE63514_GEO2R.R`  
`GSE62417_GEO2R.R`  
`GSE138080_GEO2R.R`  

**Output:** Volcano plots, DEG lists (upregulated/downregulated genes)  
Stored in `Output/DEG/` (per dataset)  

### 3. **Comparative DEG Analysis**
**Jupyter notebooks:**  
`Upregulated_common_genes.ipynb`  
`Downregulated_common_genes.ipynb`  

**Output:** Venn diagrams for common DEGs across datasets, Common Upregulated genes list and Common Downregulated genes list  
Results stored in `Output/DEG`  

### 4. **DAVID Functional Enrichment**
Functional annotation via DAVID web tool  
**Input:** All common upregulated and downregulated genes

Analyzed for significant GO terms and KEGG pathways  
Filtering (Benjamini p < 0.05) done via `DAVID_analysis.ipynb`  
Results stored in `Output/DAVID`  

### 5. **Gene Set Enrichment Analysis (GSEA)**
Preprocessing with `GSE64217_GCT_file.R` for GSE64217 dataset and the gene expression data was taken for the other two datasets  
***Note:*** This script can also be adapted for preprocessing the other two datasets.

GSEA run using standalone tool  
**Input:** complete expression data from each dataset

**Output:** stored in `Output/GSEA/` (per dataset)  
Filtering common and significant gene sets using `common_gsea.Rmd`  

### 6. **STRING PPI Network + Cytoscape Analysis**
**STRING parameters:**  
Confidence ≥ 0.900  
Hide disconnected nodes  
K-means clustering (k = 3)  

Exported to Cytoscape app:  
**cytoHubbav2.0.3 :** hub gene analysis (12 scoring methods)  
**MCODEv0.1 :** subnetwork/module detection

Store results in `STRING/Cytoscape/`  

---

**Citations:**  
Original Paper: Wu, B., & Xi, S. (2021). Bioinformatics analysis of differentially expressed genes and pathways in the development of cervical cancer. BMC cancer, 21(1), 733. https://doi.org/10.1186/s12885-021-08412-4

> Completed: May 5th 2025
