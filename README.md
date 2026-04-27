# Single Nuclei RNA-seq Analysis of Atherosclerotic Plaque

## Description
This repository contains an R script for single nuclei RNA sequencing (snRNA-seq) 
analysis of atherosclerotic plaque tissue. The analysis is actively ongoing and 
the script will be updated as the work progresses.

## Analysis Steps
1. **Data Loading** - Load raw count matrices from snRNA-seq data
2. **Quality Control** - Filter cells based on nFeature, nCount, and mitochondrial percentage
3. **Normalization** - Normalize and log-transform expression data
4. **Feature Selection** - Identify highly variable genes
5. **Dimensionality Reduction** - PCA, UMAP
6. **Clustering** - Identify cell populations
7. **Cell Type Annotation** - Annotate clusters based on marker genes

## Dependencies
- R (≥ 4.0)
- Seurat
- ggplot2

## Author
harikrishnanm0731

## Last Updated
<!-- update this when you make major changes -->
