# Load required libraries
library(Seurat)
library(tidyverse)

print("🚀 Starting Single-Cell Analysis Pipeline...")

# 1. Load the 10x Genomics dataset 
pbmc.data <- Read10X(data.dir = "data/filtered_feature_bc_matrix/")

# 2. Initialize the Seurat Object (creates the raw data matrix container)
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc10k", min.cells = 3, min.features = 200)
print(paste("Initial matrix size:", ncol(pbmc), "cells and", nrow(pbmc), "genes."))

# 3. Quality Control (QC): Calculate mitochondrial read percentages
# Dying cells spill cytoplasmic RNA, leaving mostly mitochondrial RNA behind
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# 4. Filter cells based on industry-standard thresholds
# Drop cells with fewer than 200 genes, more than 2500 genes, or > 5% mitochondrial data
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
print(paste("Matrix size after QC filtering:", ncol(pbmc), "cells."))

# 5. Global Normalization (LogNormalize scales total counts to 10,000 per cell)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# 6. Feature Selection (Find the top 2,000 highly variable genes driving variance)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# 7. Linear Dimension Reduction (Scale the data matrix and run PCA)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc), verbose = FALSE)

# 8. Graph-Based Clustering 
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# 9. Non-linear Dimension Reduction (Generate UMAP Coordinates)
pbmc <- RunUMAP(pbmc, dims = 1:10, verbose = FALSE)

# 10. Save the final cluster plot as a PNG image inside your figures folder
print(" Generate final UMAP cluster visualization...")
png("results/figures/umap_clusters.png", width = 800, height = 600)
DimPlot(pbmc, reduction = "umap", label = TRUE) + NoLegend() + ggtitle("UMAP Clustering of 10k Human PBMCs")
dev.off()

print("Pipeline complete! Your cluster plot is saved in results/figures/umap_clusters.png")
