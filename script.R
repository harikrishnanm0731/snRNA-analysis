# ============================================================
# Single Nuclei RNA-seq Analysis of Atherosclerotic Plaque
# Author  : harikrishnanm0731
# Project : snRNA-seq — SC006 & SC008
# Updated : April 2026
# ============================================================
# OVERVIEW:
# This script performs end-to-end snRNA-seq analysis of two
# atherosclerotic plaque samples (SC006, SC008) including:
#   1. Data loading & Seurat object creation
#   2. Quality control & filtering
#   3. Log normalization & variable feature selection
#   4. CCA integration across samples
#   5. Dimensionality reduction & clustering
#   6. Cell type annotation via SingleR
#   7. Differential gene expression (DEG) analysis
#   8. Cell type composition visualization
#   9. Marker gene dot plot
#  10. Cell-cell communication analysis via CellChat
# ============================================================


# ============================================================
# 0. SETUP — Working Directory & Libraries
# ============================================================

getwd()
setwd("~/snrna-project")

library(Seurat)
library(patchwork)
library(dplyr)
library(harmony)
library(ggplot2)
library(SingleR)
library(celldex)
library(metap)
library(SingleCellExperiment)

# Install CellChat from GitHub (run only once)
devtools::install_github("jinworks/CellChat")
library(CellChat)

# Load Human Primary Cell Atlas reference for SingleR annotation
ref <- celldex::HumanPrimaryCellAtlasData()


# ============================================================
# 1. DATA LOADING
# ============================================================
# CellRanger output path (base directory for all samples)
base_path <- "SO_14534/CellRanger_Output/"

# Read filtered feature-barcode matrices for each sample
data_SC006 <- Read10X(data.dir = file.path(base_path, "SC006", "sample_filtered_feature_bc_matrix"))
data_SC008 <- Read10X(data.dir = file.path(base_path, "SC008", "sample_filtered_feature_bc_matrix"))


# ============================================================
# 2. CREATE SEURAT OBJECTS
# ============================================================
# min.cells = 1  : keep genes detected in at least 1 cell
# min.features = 100 : keep cells with at least 100 detected genes

SC008 <- CreateSeuratObject(counts = data_SC008, project = "SC008", min.cells = 1, min.features = 100)
SC006 <- CreateSeuratObject(counts = data_SC006, project = "SC006", min.cells = 1, min.features = 100)


# ============================================================
# 3. QUALITY CONTROL
# ============================================================

# Add sample ID metadata
SC006$sample_id <- "SC006"
SC008$sample_id <- "SC008"

# Calculate mitochondrial gene percentage
# High % suggests damaged/dying cells
SC006[["percent.mt"]] <- PercentageFeatureSet(SC006, pattern = "^MT-|^mt-")
SC008[["percent.mt"]] <- PercentageFeatureSet(SC008, pattern = "^MT-|^mt-")

# --- QC Violin Plots BEFORE filtering ---
VlnPlot(SC006, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(SC008, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# --- Filter cells ---
# nFeature_RNA > 200     : remove empty droplets
# nFeature_RNA < 4500    : remove likely doublets
# percent.mt < 5 / < 1   : remove damaged cells
SC006 <- subset(SC006, subset = nFeature_RNA > 200 & nFeature_RNA < 4500 & percent.mt < 5)
SC008 <- subset(SC008, subset = nFeature_RNA > 200 & nFeature_RNA < 4500 & percent.mt < 1)

# --- QC Violin Plots AFTER filtering ---
VlnPlot(SC006, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(SC008, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


# ============================================================
# 4. NORMALIZATION & VARIABLE FEATURE SELECTION
# ============================================================
# LogNormalize: normalize to 10,000 reads per cell, then log-transform
# FindVariableFeatures: select top 2000 highly variable genes (HVGs)
# using variance-stabilizing transformation (VST)

SC006 <- NormalizeData(SC006, normalization.method = "LogNormalize", scale.factor = 10000)
SC006 <- FindVariableFeatures(SC006, selection.method = "vst", nfeatures = 2000)

SC008 <- NormalizeData(SC008, normalization.method = "LogNormalize", scale.factor = 10000)
SC008 <- FindVariableFeatures(SC008, selection.method = "vst", nfeatures = 2000)


# ============================================================
# 5. INTEGRATION USING CCA (Canonical Correlation Analysis)
# ============================================================
# CCA integration corrects for batch effects between SC006 & SC008
# while preserving biological variation

seurat_list <- list(SC006, SC008)

# Find shared anchors across samples
anchors <- FindIntegrationAnchors(
  object.list        = seurat_list,
  normalization.method = "LogNormalize",
  anchor.features    = 2000,
  reduction          = "cca"
)

# Integrate data into a single Seurat object
combined <- IntegrateData(
  anchorset            = anchors,
  normalization.method = "LogNormalize"
)


# ============================================================
# 6. SCALING & DIMENSIONALITY REDUCTION
# ============================================================
# ScaleData: zero-mean, unit-variance scaling (required before PCA)
# RunPCA: linear dimensionality reduction
# ElbowPlot: helps choose the number of PCs to use downstream

DefaultAssay(combined) <- "integrated"

combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, verbose = FALSE)
ElbowPlot(combined, ndims = 50)   # inspect to confirm dims = 1:20 is appropriate

# --- Clustering ---
# FindNeighbors: builds KNN graph using top 20 PCs
# FindClusters: Louvain algorithm, resolution = 0.3 (lower = fewer clusters)
combined <- FindNeighbors(combined, dims = 1:20)
combined <- FindClusters(combined, resolution = 0.3)

# --- UMAP Visualization ---
combined <- RunUMAP(combined, dims = 1:20)

# Plot by sample to check integration quality
DimPlot(combined, group.by = "sample_id")

# Plot clusters
DimPlot(combined, label = TRUE)

# Split by sample to compare cluster distribution
DimPlot(combined,
        split.by  = "sample_id",
        group.by  = "seurat_clusters",
        label     = TRUE,
        repel     = TRUE) +
  ggtitle("Clusters per Sample")


# ============================================================
# 7. CELL TYPE ANNOTATION USING SingleR
# ============================================================
# SingleR compares cluster expression profiles to a reference
# dataset (Human Primary Cell Atlas) to assign cell type labels

DefaultAssay(combined) <- "RNA"

# Join layers required for Seurat v5 before converting to SCE
combined[["RNA"]] <- JoinLayers(combined[["RNA"]])

# Convert to SingleCellExperiment format for SingleR
sce_com <- as.SingleCellExperiment(combined)

# Run SingleR — annotate at the cluster level (faster, more robust)
pred_com <- SingleR(
  test     = sce_com,
  ref      = ref,
  labels   = ref$label.fine,
  clusters = combined$seurat_clusters
)

# Map SingleR labels back to individual cells
combined$singler_labels <- pred_com$labels[match(
  combined$seurat_clusters,
  rownames(pred_com)
)]

# Visualize SingleR annotation alongside sample distribution
p1 <- DimPlot(combined, group.by = "singler_labels", label = TRUE, repel = TRUE) +
  ggtitle("SingleR Cell Type Annotation")
p2 <- DimPlot(combined, group.by = "sample_id")
p1 + p2


# ============================================================
# 8. DIFFERENTIAL GENE EXPRESSION (DEG) ANALYSIS
# ============================================================
# FindAllMarkers: finds marker genes for each cluster vs all others
# only.pos = TRUE       : only upregulated markers
# min.pct = 0.25        : gene must be detected in ≥25% of cells
# logfc.threshold = 0.25: minimum log2 fold-change

DefaultAssay(combined) <- "RNA"
combined <- JoinLayers(combined)   # merges SC006 + SC008 layers into one

markers_com <- FindAllMarkers(
  combined,
  only.pos       = TRUE,
  min.pct        = 0.25,
  logfc.threshold = 0.25
)

# Extract top 50 significant markers per cluster
top50 <- markers_com %>%
  filter(p_val_adj < 0.05, avg_log2FC > 0.5) %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 50)


# ============================================================
# 9. MANUAL CELL TYPE ANNOTATION
# ============================================================
# Clusters manually annotated based on marker genes and SingleR results

combined$cell_type <- dplyr::recode(
  as.character(combined$seurat_clusters),
  "0" = "Endothelial_General",
  "1" = "Macrophages",
  "2" = "Vascular_Smooth_Muscle_Cells",
  "3" = "Endothelial_EndoEMT",
  "4" = "Mitochondrial_Contaminated_Cells",
  "5" = "Macrophages_Foamy",
  "6" = "Unknown_Proliferating_Cells"
)

# UMAP split by sample with manual cell type labels
DimPlot(combined,
        split.by = "sample_id",
        group.by = "cell_type",
        repel    = TRUE) +
  ggtitle("Cell Types per Sample")

# Combined UMAP with cell type labels
DimPlot(combined,
        group.by = "cell_type",
        repel    = TRUE) +
  ggtitle("Cell Types — Combined")


# ============================================================
# 10. REMOVE LOW-QUALITY CLUSTERS
# ============================================================
# Remove mitochondrial contaminated cells before downstream analysis

combined_clean <- subset(combined,
                         subset = cell_type != "Mitochondrial_Contaminated_Cells")

# Verify cell counts per type after removal
table(combined_clean$cell_type)


# ============================================================
# 11. CELL TYPE COMPOSITION PLOT
# ============================================================
# Compare the number of cells per cell type across samples

library(ggplot2)
library(dplyr)

cell_counts <- combined_clean@meta.data %>%
  group_by(sample_id, cell_type) %>%
  summarise(n_cells = n(), .groups = "drop")

ggplot(cell_counts, aes(x = cell_type, y = n_cells, fill = sample_id)) +
  geom_bar(stat = "identity", position = "dodge",
           width = 0.7, color = "white", linewidth = 0.3) +
  scale_fill_manual(values = c(
    "SC006" = "#4E79A7",
    "SC008" = "#F28E2B"
  )) +
  labs(
    title = "Cell Type Composition per Sample",
    x     = NULL,
    y     = "Number of Cells",
    fill  = "Sample"
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title      = element_text(face = "bold", size = 11, hjust = 0.5, family = "Helvetica"),
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 8, family = "Helvetica"),
    axis.text.y     = element_text(size = 8, family = "Helvetica"),
    axis.title.y    = element_text(size = 9, family = "Helvetica"),
    legend.title    = element_text(size = 9, face = "bold", family = "Helvetica"),
    legend.text     = element_text(size = 8, family = "Helvetica"),
    legend.position = "top",
    panel.border    = element_rect(color = "grey80", fill = NA, linewidth = 0.5)
  )


# ============================================================
# 12. MARKER GENE DOT PLOT
# ============================================================
# Dot plot showing expression of known plaque marker genes
# Dot size = % cells expressing the gene
# Dot color = average expression level

genes_to_plot <- c(
  # Foamy Macrophages (Cluster 5)
  "CD68", "CHIT1", "CD36", "CD74", "LGALS3", "CTSD", "APOE", "MMP9",
  # Endothelial — General (Cluster 0)
  "VWF", "PECAM1", "ECSCR", "MGP",
  # Endothelial — EndoEMT (Cluster 3)
  "COL1A2", "COL1A1", "FN1",
  # Vascular Smooth Muscle Cells (Cluster 2)
  "ACTA2", "MYL9", "TAGLN"
)

genes_to_plot <- unique(genes_to_plot)   # remove any duplicates

Idents(combined_clean)        <- "cell_type"
DefaultAssay(combined_clean)  <- "RNA"

DotPlot(combined_clean,
        features = genes_to_plot,
        group.by = "cell_type",
        dot.scale = 6) +
  
  scale_color_gradient2(
    low      = "#4E79A7",
    mid      = "white",
    high     = "#E15759",
    midpoint = 0
  ) +
  
  labs(
    title = "Marker Gene Expression Across Cell Types",
    x     = NULL,
    y     = NULL,
    color = "Avg Expression",
    size  = "% Expressed"
  ) +
  
  theme_classic(base_size = 10) +
  theme(
    plot.title       = element_text(face = "bold", size = 11, hjust = 0.5, family = "Helvetica"),
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 8, family = "Helvetica", face = "italic"),
    axis.text.y      = element_text(size = 8, family = "Helvetica"),
    legend.title     = element_text(size = 8, face = "bold", family = "Helvetica"),
    legend.text      = element_text(size = 7, family = "Helvetica"),
    legend.position  = "right",
    panel.border     = element_rect(color = "grey80", fill = NA, linewidth = 0.5),
    panel.grid.major = element_line(color = "grey95", linewidth = 0.3)
  ) +
  
  guides(
    color = guide_colorbar(barwidth = 0.8, barheight = 4),
    size  = guide_legend(override.aes = list(color = "grey40"))
  )


# ============================================================
# 13. CELL-CELL COMMUNICATION ANALYSIS — CellChat
# ============================================================
# CellChat infers intercellular communication networks from
# ligand-receptor interaction databases
# Analysis run separately per sample for comparison

# --- Split cleaned object by sample ---
Idents(combined_clean) <- "sample_id"

seurat_SC006 <- subset(combined_clean, idents = "SC006")
seurat_SC008 <- subset(combined_clean, idents = "SC008")

seurat_list <- list(
  SC006 = seurat_SC006,
  SC008 = seurat_SC008
)

# --- CellChat wrapper function ---
# Runs full CellChat pipeline for a single Seurat object
run_cellchat <- function(seurat_obj) {
  
  library(CellChat)
  
  DefaultAssay(seurat_obj) <- "RNA"
  seurat_obj$cell_type     <- as.factor(seurat_obj$cell_type)
  
  # Create CellChat object grouped by cell type
  cellchat <- createCellChat(
    object   = seurat_obj,
    group.by = "cell_type"
  )
  
  # Use human ligand-receptor database
  CellChatDB        <- CellChatDB.human
  cellchat@DB       <- CellChatDB
  
  # Preprocessing: identify overexpressed genes and interactions
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat, do.fast = FALSE)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  # Compute communication probabilities
  cellchat <- computeCommunProb(cellchat)
  cellchat <- filterCommunication(cellchat, min.cells = 10)   # remove low-confidence interactions
  
  # Summarise at the signalling pathway level
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  return(cellchat)
}

# --- Run CellChat for each sample ---
cellchat_list <- lapply(seurat_list, run_cellchat)