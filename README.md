# Single Nuclei RNA-seq Analysis of Atherosclerotic Plaque

## Description
snRNA-seq analysis of atherosclerotic plaque tissue using two samples (SC006, SC008).  
The pipeline covers the full workflow from raw CellRanger output to cell-cell communication analysis.

## Analysis Steps

**1. Data Loading**  
Loaded filtered feature-barcode matrices from CellRanger output for SC006 and SC008.

**2. Quality Control**  
Calculated mitochondrial gene percentage and filtered cells based on nFeature_RNA and percent.mt thresholds. SC008 used a stricter mitochondrial cutoff (< 1%) compared to SC006 (< 5%).

**3. Normalization**  
LogNormalize method — normalized each cell to 10,000 reads then log-transformed. Top 2000 highly variable genes selected using VST.

**4. CCA Integration**  
Integrated SC006 and SC008 using Canonical Correlation Analysis (CCA) to correct for batch effects while preserving biological variation.

**5. Dimensionality Reduction & Clustering**  
PCA followed by UMAP. KNN graph built on top 20 PCs, Louvain clustering at resolution 0.3.

**6. Cell Type Annotation**  
Automated annotation using SingleR against the Human Primary Cell Atlas reference, applied at the cluster level.

**7. DEG Analysis**  
Marker genes identified per cluster using `FindAllMarkers` (only positive markers, min.pct = 0.25, log2FC > 0.25). Top 50 significant markers per cluster extracted.

**8. Manual Annotation**  
Clusters manually labelled based on DEG results and SingleR output. Mitochondrial contaminated cells identified and removed.

**9. Composition Plot**  
Bar plot comparing cell type proportions across SC006 and SC008.

**10. Marker Gene Dot Plot**  
Dot plot showing expression of known plaque marker genes (macrophage, endothelial, VSMC markers) across cell types.

**11. CellChat Analysis**  
Cell-cell communication analysis run separately per sample using the human ligand-receptor database to infer intercellular signalling networks.

## Dependencies
- Seurat
- Harmony
- SingleR / celldex
- CellChat
- ggplot2 / dplyr / patchwork

## Author
harikrishnanm0731

## Last Updated
April 2026
