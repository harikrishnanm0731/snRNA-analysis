
getwd()
setwd("~/snrna-project")  

library(Seurat)
library(patchwork)
library(dplyr)
library(harmony)

library(SingleR)
library(celldex)
library(SingleCellExperiment)
ref <- celldex::HumanPrimaryCellAtlasData()


#####set the path###############
base_path <- "SO_14534/CellRanger_Output/"

##########load datasets ###########
data_SC006 <- Read10X(data.dir = file.path(base_path, "SC006", "sample_filtered_feature_bc_matrix")) 
data_SC008 <- Read10X(data.dir = file.path(base_path, "SC008", "sample_filtered_feature_bc_matrix")) 



###########Create Seurat Objects####################################
SC008 <- CreateSeuratObject(counts = data_SC008, project = "SC008", min.cells = 1, min.features = 100) 
SC006 <- CreateSeuratObject(counts = data_SC006, project = "SC006", min.cells = 1, min.features = 100) 




# ================================
# SC006
# ================================
SC006[["percent.mt"]] <- PercentageFeatureSet(SC006, pattern = "^MT-|^mt-")
VlnPlot(SC006, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(SC006, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(SC006, feature1 = "nCount_RNA", feature2 = "percent.mt")

SC006<- subset(SC006, subset = nFeature_RNA > 500 & nFeature_RNA < 3000& percent.mt < 5)

SC006 <- NormalizeData(SC006)
SC006 <- FindVariableFeatures(SC006, selection.method = "vst", nfeatures = 2000)
VariableFeaturePlot(SC006)
SC006 <- ScaleData(SC006)
SC006 <- RunPCA(SC006, npcs = 50)
ElbowPlot(SC006,ndims = 40)

SC006 <- FindNeighbors(SC006, dims = 1:30)
SC006 <- FindClusters(SC006, resolution = 0.8)
SC006 <- RunUMAP(SC006, dims = 1:30)
DimPlot(SC006, reduction = "umap", label = TRUE)


################Sinlge R annotaion -use only just to check the clusters##############

sce006 <- as.SingleCellExperiment(SC006)
pred006 <- SingleR(test = sce006, ref = ref, labels = ref$label.main, clusters = SC006$seurat_clusters)
SC006$SingleR_labels <- pred006$labels[SC006$seurat_clusters]
DimPlot(SC006, group.by = "SingleR_labels", label = TRUE)


# ============================================================
#  FIND ALL MARKERS (DEGs per cluster)
# ============================================================
markers_SC006 <- FindAllMarkers(
  SC006,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

## ============================================================
# DEG ANALYSIS PER CLUSTER
# ============================================================

# Top 10 per cluster by log2FC
top20 <- markers_SC006 %>%
  filter(p_val_adj < 0.05, avg_log2FC > 0.5) %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 20)

print(top20 %>% select(cluster, gene, avg_log2FC, pct.1, pct.2, p_val_adj), n = 200)
write.csv(top20, "SC006_top20_markers_per_cluster.csv", row.names = FALSE)



# ================================
# SC008
# ================================
SC008[["percent.mt"]] <- PercentageFeatureSet(SC008, pattern = "^MT-|^mt-")
VlnPlot(SC008, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(SC008, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(SC008, feature1 = "nCount_RNA", feature2 = "percent.mt")

SC008 <- subset(SC008, subset = nFeature_RNA > 500 & nFeature_RNA < 4500 & percent.mt < 5)

SC008 <- NormalizeData(SC008)
SC008 <- FindVariableFeatures(SC008, selection.method = "vst", nfeatures = 2000)
VariableFeaturePlot(SC008)
SC008 <- ScaleData(SC008)
SC008 <- RunPCA(SC008, npcs = 50)
ElbowPlot(SC008, ndims=50)

SC008 <- FindNeighbors(SC008, dims = 1:30)
SC008 <- FindClusters(SC008, resolution = 0.5)
SC008 <- RunUMAP(SC008, dims = 1:30)
DimPlot(SC008, reduction = "umap", label = TRUE)


########single R

sce008 <- as.SingleCellExperiment(SC008)
pred08 <- SingleR(test = sce008, ref = ref, labels = ref$label.main, clusters = SC008$seurat_clusters)
SC008$SingleR_labels <- pred008$labels[SC008$seurat_clusters]
DimPlot(SC008, group.by = "SingleR_labels", label = TRUE)


############################checking the pca ####################################

# SC006
VizDimLoadings(SC006, dims = 1:4, reduction = "pca")

# SC008
VizDimLoadings(SC008, dims = 1:4, reduction = "pca")


load_SC006 <- SC006[["pca"]]@feature.loadings
load_SC008 <- SC008[["pca"]]@feature.loadings

common_genes <- intersect(rownames(load_SC006), rownames(load_SC008))


load_SC006 <- load_SC006[common_genes, ]
load_SC008 <- load_SC008[common_genes, ]

cor(load_SC006[,1], load_SC008[,1])


############merging###############just concatenating the datasets , not integrating

SC006$sample_id <- "SC006"
SC008$sample_id <- "SC008"

combined <- merge(SC006, SC008)

combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined, nfeatures = 2000)
combined <- ScaleData(combined)
combined <- RunPCA(combined)

ElbowPlot(combined,ndims = 50)


combined <- FindNeighbors(combined, dims = 1:20)
combined <- FindClusters(combined, resolution = 0.5)
combined <- RunUMAP(combined, dims = 1:20)

DimPlot(combined, group.by = "sample_id")

DimPlot(combined, group.by = "seurat_clusters")


################intergrating##########
library(harmony)

combined <- RunHarmony(combined, group.by.vars = "sample_id")


combined <- FindNeighbors(combined, reduction = "harmony", dims = 1:20)
combined <- FindClusters(combined, resolution = 0.5)
combined <- RunUMAP(combined, reduction = "harmony", dims = 1:20)

DimPlot(combined, group.by = "sample_id")



#####single R###################################

sce <- as.SingleCellExperiment(combined, assay = "RNA")
pred <- SingleR(test = sce, ref = ref, labels = ref$label.main, clusters = combined$seurat_clusters)
combined$SingleR_labels <- pred$labels[combined$seurat_clusters]
DimPlot(combined, group.by = "SingleR_labels", label = TRUE)



##########################################












library(dplyr)
library(tidyr)
library(tibble)
# below is the marker genes from integrated cell-atlas of atherosclerotic plaques

marker_list <- list(
  
  # ---------------- Structural ----------------
  Fibroblast = c("LUM","DCN","COL1A1","COL1A2","FBLN1","THY1","C3","C7"),
  Fibromyocyte = c("FN1","LUM","TNFRSF11B","ACTA2","TCF21"),
  Smooth_muscle_cell = c("ACTA2","MYH11","MYL9","TPM2","CALD1","TAGLN","TNFRSF11B","LUM","APOE","APOC1","AGT","NOTCH3","PDGFRB","MFAP4"),
  
  Endothelial_general = c("PECAM1","VWF","FABP4","CLDN5","IFI27","ECSCR","DYSF",
                          "CD34","COL4A1","COL4A2","SPARCL1","PLVAP","MPZL2","SULF1","EDN1"),
  
  Endothelial_proangiogenic = c("ACKR1","AQP1","FABP4","CXCL12"),
  Endothelial_EndoMT = c("COL1A2","FN1"),
  
  # ---------------- Myeloid ----------------
  Neutrophil = c("NAMPT","IFITM2","G0S2","CXCL8","NEAT1","SRGN",
                 "AQP9","SOD2","FCGR3B","IVNS1ABP"),
  
  Monocyte = c("FCN1","S100A8","S100A9","S100A12","VCAN","CD52","LYZ","CTSS"),
  
  Mast_cell = c("TPSAB1","TPSB2","KIT","HDC","CMA1"),
  
  DC_general = c("CLEC10A","FCER1A","CD1C","HLA-DRA","HLA-DRB1"),
  DC_cDC1 = c("CLEC9A","IRF8","SNX3"),
  DC_cDC2 = c("CD1C","CLEC10A","FCER1A"),
  DC_pDC = c("GTF2A","GZMB","TLR7","TLR9","NRP1","SCAMP5","CLEC4C","IRF7"),
  
  Macrophage_general = c("C1QA","C1QB","C1QC","CD74","CXCL8","AIF1","CD14",
                         "CD68","ITGAM","CSF1R","HLA-DRA","LGALS3"),
  
  Macrophage_foamy = c("TREM2","MARCO","FABP4","FABP5","CD36"),
  
  Macrophage_inflammatory = c("S100A8","IL1B","S100A9","IRF7","IFITM3","ISG15","IFIT2"),
  
  Macrophage_PLIN2_TREM1 = c("PLIN2","TREM1",
                             "CXCL1","CXCL2","CXCL3","CXCL8",
                             "CCL2","CCL7","CCL20",
                             "IL1B","TNF","CEBPB"),
  
  Macrophage_HMOX1 = c("HMOX1","ALOX5","GPNMB","LIPA","NPC2","PRDX1","SLC40A1",
                       "NUPR1","APOE","LAMP2","CTSB","SELENOP","LGMN","LRP1","CTSD","FTL"),
  
  # ---------------- Lymphoid ----------------
  B_cell = c("CD79A","CD79B","MS4A1","IGKC","CD22","FCER2"),
  
  Plasma_cell = c("IGKC","IGHM","IGHA1","IGLC2","IGLC3","JCHAIN"),
  
  NK_cell = c("NKG7","XCL1","CTSW","XCL2","CD160","FCGR3A","PRF1","GNLY"),
  
  T_cell_general = c("CD2","TRAC","CD69","CD3E","CD3D","CD4","CD8A","CD8B","EOMES","LAG3"),
  
  T_cell_CD4 = c("CD4","IFI1","IL7R","ANXA1","BATF","TNFRSF4","TNFRSF18"),
  
  T_cell_CD8 = c("CCL4L2","CRTAM","GZMK","CD8A","CCL4","GZMH")
)
############# markers from the paper-decoding....#########

marker_list2 <- list(
  Macrophages = c("CD14", "CD68", "AIF1", "LST1"),
  
  Endothelial_Cells = c("CLU", "VWF", "EDN1", "ECSCR", "SPARCL1", "PECAM1", "CALD1", "MGP"),
  
  ACKR1_Positive = c("ACKR1"),
  
  Natural_KT_Cells = c("NKG7", "XCL1", "CTSW", "CD69"),
  
  T_Cells = c("TRAC", "CD2"),
  
  Vascular_Smooth_Muscle_Cells = c("ACTA2", "TAGLN", "MYL9", "SPARCL1", "CALD1", "MGP", "DCN"),
  
  Neutrophils = c("S100A8")
)


######macrophage markers from-Macrophage subsets in atherosclerosis as defined by single‐cell technologies########

resident_genes <- c("LYVE1", "CX3CR1", "FOLR2", "MRC1", "F13A1", "CBR2", "SEPP1", "PF4", "GAS6")

inflammatory_genes <- c("TNF", "NLRP3", "IL1B", "EGR1", "TLR2", "IER3", "CEBPB",
                        "CXCL2", "CCL2", "CCL3", "CCL4", "CCL5", "NFKBIA")

trem2_genes <- c("TREM2", "CD9", "LGALS3", "CTSB", "SPP1")



resident_genes %in% rownames(SC006)
inflammatory_genes %in% rownames(SC006)
trem2_genes %in% rownames(SC006)


DotPlot(
  SC006,
  features = c(resident_genes, inflammatory_genes, trem2_genes)
) + RotatedAxis()
