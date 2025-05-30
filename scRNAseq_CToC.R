# Load the necessary packages
library(Matrix)
library(sctransform)
library(hdf5r)
library(harmony)
library(dplyr)
library(scran)
library(SingleCellExperiment)
library(glmGamPoi)
library(ggplot2)

# Read counts matrix files and preprocess them
bone4 <- Read10X_h5("BoneChipDay4/sample_filtered_feature_bc_matrix.h5")
bone7 <- Read10X_h5("BoneChipDay7/sample_filtered_feature_bc_matrix.h5")
bone_seurat4 <- CreateSeuratObject(counts = bone4)
bone_seurat7 <- CreateSeuratObject(counts = bone7)
bone_seurat4$orig.ident <- "Day4"
bone_seurat7$orig.ident <- "Day7"

#Based on the violin plots and the feature scatter plots choose the maximum number of genes to be included on downstream analysis
bone_seurat4 <- PercentageFeatureSet(bone_seurat4, pattern = "^MT-", col.name = "percent.mt")
plot1 <- FeatureScatter(bone_seurat4, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(bone_seurat4, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
plot_mt4 <- VlnPlot(bone_seurat4, features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3)
plot_mt4
bone_seurat4 <- subset(bone_seurat4, subset = nFeature_RNA > 100 & nFeature_RNA < 9500 & percent.mt < 5)


bone_seurat7 <- PercentageFeatureSet(bone_seurat7, pattern = "^MT-", col.name = "percent.mt")
plot3 <- FeatureScatter(bone_seurat7, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot4 <- FeatureScatter(bone_seurat7, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 + plot4
plot_mt7 <- VlnPlot(bone_seurat7, features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3)
plot_mt7
bone_seurat7 <- subset(bone_seurat7, subset = nFeature_RNA > 100 & nFeature_RNA < 9500 & percent.mt < 3)

## Merge, Normalize and Scale the data
#Normalize both datasets prior to combining them. 
# Merge datasets
combined <- merge(bone_seurat4, y = bone_seurat7, add.cell.ids = c("day4", "day7"), project = "BoneProject")


# Normalize cells from Day 4 & Day 7
combined <- NormalizeData(combined, normalization.method = "LogNormalize", scale.factor = 10000)
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 1000)
all.genes <- rownames(combined)
combined <- ScaleData(combined, features = all.genes)
combined <- RunPCA(combined, features = VariableFeatures(object = combined))

combined <- IntegrateLayers(object =combined, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = 'harmony', verbose = FALSE)
combined <- JoinLayers(combined)
combined <- FindNeighbors(combined, reduction = "harmony", dims = 1:20)
combined <- FindClusters(combined, resolution = c(0.15, 0.25, 0.35, 0.5, 0.7, 1))
combined <- RunUMAP(combined, reduction = "harmony", dims = 1:20)

DimPlot(combined, reduction = "harmony", group.by = c("orig.ident","RNA_snn_res.0.35"))
DimPlot(combined, reduction = "umap", split.by = "orig.ident", group.by = "RNA_snn_res.0.35")
    
# Find Markers for different clusters
# Gene expression exported and analyzed for annotation
FindMarkers(combined, ident.1 = 5, test.use = "roc", logfc.threshold = 0.25, only.pos = TRUE) 
FindMarkers(combined, ident.1 = 2, test.use = "roc", logfc.threshold = 0.25, only.pos = TRUE) 
FindMarkers(combined, ident.1 = 1, test.use = "roc", logfc.threshold = 0.25, only.pos = TRUE) 
FindMarkers(combined, ident.1 = 10, test.use = "roc", logfc.threshold = 0.25, only.pos = TRUE) 
FindMarkers(combined, ident.1 = 3, test.use = "roc", logfc.threshold = 0.25, only.pos = TRUE) 
FindMarkers(combined, ident.1 = 7, test.use = "roc", logfc.threshold = 0.25, only.pos = TRUE) 
FindMarkers(combined, ident.1 = 0, test.use = "roc", logfc.threshold = 0.25, only.pos = TRUE) 
FindMarkers(combined, ident.1 = 11, test.use = "roc", logfc.threshold = 0.25, only.pos = TRUE) 

new.cluster.ids <- c("Osteoblast", "Pre-adipocytes", "Endothelial", "Unknown", "Fibroblast", "Endothelial",
                     "MSC-Endothelial", "M2 macrophage", "Proliferative cells","Adipocytes")
names(new.cluster.ids) <- levels(combined)
comb.1 <- RenameIdents(combined, new.cluster.ids)
saveRDS(comb.1, "CToC_UMAP.rds")

# Analyzing "Unknown" subcluster of cells reveals a mixed population of cells
ctoc <- readRDS("CToC_UMAP.rds")
table(ctoc$cell_type)
subcluster <- subset(ctoc, subset = cell_type == "Unknown")
subcluster
subcluster <- NormalizeData(subcluster)
subcluster <- FindVariableFeatures(subcluster, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(subcluster), 25)
top10
plot1 <- VariableFeaturePlot(subcluster)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1+plot2
all.genes <- rownames(subcluster)
subcluster <- ScaleData(subcluster, features = all.genes)
subcluster <- RunPCA(subcluster, features = VariableFeatures(object = subcluster))
VizDimLoadings(subcluster, dims = 1:5, reduction = "pca")
DimHeatmap(subcluster, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(subcluster)
subcluster <- FindNeighbors(subcluster, dims = 1:15)
subcluster <- FindClusters(subcluster, resolution = 0.25)
subcluster <- RunUMAP(subcluster, dims = 1:20)
DimPlot(subcluster, reduction = "umap", split.by = "orig.ident")
DimPlot(subcluster, reduction = "umap", group.by = "orig.ident")
FindAllMarkers(subcluster, test.use = "roc", logfc.threshold = 0.25, pos = TRUE)




