# Deciphering Molecular Pathways in Myelodysplastic Syndromes (MDS)

# Overview
This project focuses on analyzing publicly available gene expression datasets to uncover dysregulated pathways and biomarkers associated with Myelodysplastic Syndromes (MDS). The datasets used are:
1. *GSE107490*: Mesenchymal stem cells (MSCs) derived from MDS patients and healthy controls.
2. *GSE61853*: Expression profiles from MDS patients and controls.

The project leverages R and Python for data preprocessing, differential expression analysis, pathway enrichment, and machine learning-based clustering.

# Data Preprocessing
- R:
  - Gene expression data normalized using `normalizeBetweenArrays` (quantile normalization).
  - Differentially expressed genes (DEGs) were identified using the `limma` package.
  - Top DEGs filtered based on log fold change (logFC ≥ 2) and adjusted p-value (p < 0.05).
  
- Python:
  - Applied K-means clustering on PCA-reduced data.
  - Hierarchical clustering with Spearman correlation.
  - Machine learning approaches using Random Forest for classification.

# Results
# Differential Expression Analysis
- *GSE107490*:
  - *Upregulated genes*: 1,040 (e.g., `IFI6`, `LDHA`).
  - *Downregulated genes*: 1,277 (e.g., `SERPINE2`, `CCN2`).

- *GSE61853*:
  - *Upregulated genes*: 1,852.
  - *Downregulated genes*: 1,519.

# Gene Ontology (GO) Enrichment
- Enriched biological processes include:
  - mRNA processing and protein synthesis.
  - Nuclear-transcribed mRNA catabolic processes.
  - SRP-dependent co-translational protein targeting to membranes.

# KEGG Pathway Enrichment
- Key enriched pathways:
  - Ribosome and protein synthesis pathways.
  - Cellular response to infection (COVID-19, Legionellosis, Salmonella).
  - Oxytocin signaling pathway, potentially involved in MDS pathogenesis.

# Clustering and Network Analysis
- K-means Clustering:
  - Identified three distinct clusters for both datasets, reflecting genetic subtypes.
  - Segments provide insights into variability in prognosis and treatment response.

- Hierarchical Clustering:
  - Generated dendrograms to visualize genetic similarities among samples.
  - Highlighted regulatory networks disrupted in MDS.

- Random Forest Classification:
  - GSE61853 model accuracy: 0.21%.
  - Indicates high complexity and heterogeneity in MDS datasets.

# Visualization
# R Analysis
1. Top DEGs Heatmaps:
   - Clustering of genes and samples based on expression profiles.
2. PCA and K-means Plots:
   - Visualized variance among clusters.

# Python Analysis
1. Spearman Correlation Heatmaps:
   - Highlighted gene-gene relationships.
2. Dendrograms:
   - Displayed hierarchical relationships among samples.

# Top 25 DEGs
The tables below summarize the top 25 upregulated and downregulated genes identified for both datasets.

# GSE107490 (R Results)
| **Gene Symbol** | **LogFC** | **P-value** | **Adj P-value** |
|------------------|-----------|-------------|------------------|
| SERPINE2        | -16,619   | 0.0139      | 0.6386          |
| IFI6            | 15,387    | 0.0022      | 0.6296          |
| LDHA            | -13,903   | 0.0037      | 0.6296          |

### GSE61853 (Python Results)
| **Gene Symbol** | **LogFC** | **P-value** | **Adj P-value** |
|------------------|-----------|-------------|------------------|
| EEF1A1          | -803,289  | 0.0257      | 0.3617          |
| COL1A1          | 773,312   | 0.0431      | 0.3978          |
| THBS1           | 463,324   | 0.0005      | 0.2361          |

(Complete tables available in the results CSV files.)

# Conclusion
The integration of statistical analysis and machine learning revealed critical insights into the molecular pathways associated with MDS, providing potential targets for therapeutic intervention. Despite challenges such as the low accuracy of Random Forest classification, the findings underscore the importance of multi-omics approaches in understanding complex diseases like MDS.




# How to Use
1. Clone the repository:
   ```bash
   git clone https://github.com/your-repo/MDS_analysis.git
   cd MDS_analysis
   ```
2. Run the R script for DEG analysis:
   ```bash
   Rscript scripts/MDS.R
   ```
3. Perform clustering using Python:
   ```bash
   python3 scripts/MDS_.ipynb
   ```

# References
1. GEO Accession Viewer: [GSE107490](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE107490), [GSE61853](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE61853).
2. Yin, C. et al., 2024. "Sequential gene expression analysis of MDS..." DOI: [10.1186/s12885-024-11859-w](https://doi.org/10.1186/s12885-024-11859-w).
