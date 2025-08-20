# Load required libraries
library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(R.utils)
library(dplyr)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db) 
library(enrichplot)
library(patchwork)
# Extract the downloaded raw tar file into a directory
untar("C:/Users/3nany/Downloads/GSE180665_RAW.tar", exdir = "SE180665_RAW")
# Get current working directory
getwd()
# Define folder path containing extracted files
folder.path = "c:/Users/3nany/Documents/GitHub/myblog/SE180665_RAW/"
# List all .tar.gz files inside the folder
tar.gz.files = list.files(folder.path, pattern = "\\.tar\\.gz$", full.names = TRUE)
# Loop to untar each .tar.gz file into a separate folder
for(file in tar.gz.files){
  out = sub("\\.tar\\.gz$", "", file)
  untar(file, exdir = out)
}
# Get a list of directories created (each sample folder)
dirs = list.dirs(path = "C:/Users/3nany/Documents/GitHub/myblog/SE180665_RAW/", recursive = FALSE, full.names = TRUE)
# Loop through each directory, read matrix files, and create Seurat objects
for (x in dirs) {
  name = gsub("_filtered_feature_bc_matrix", "", basename(x))
  # Read count matrix, features, and barcodes
  cts = ReadMtx(
    mtx = file.path(x, "matrix.mtx.gz"),
    features = file.path(x, "features.tsv.gz"),
    cells = file.path(x, "barcodes.tsv.gz")
  )
  # Create Seurat object from counts
  seurat.obj = CreateSeuratObject(counts = cts)
  # Calculate mitochondrial gene percentage
  seurat.obj$mtpercent = PercentageFeatureSet(seurat.obj, pattern = "^MT-|^mt-")
  # Randomly sample up to 3000 cells per sample (for balancing)
  num_cells = min(3000, ncol(seurat.obj))
  cells_use = sample(colnames(seurat.obj), num_cells)
  seurat.obj = subset(seurat.obj, cells = cells_use)
  # Save the Seurat object with sample-specific name
  assign(name, seurat.obj)
  # Plot QC violin plots (nFeature_RNA, nCount_RNA, mtpercent) for each sample
  print(VlnPlot(seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "mtpercent"), 
                ncol = 3) + ggtitle(name))
}
# Apply QC filtering to each sample individually
HB17_background = subset(HB17_background, subset = nFeature_RNA >= 200 & nFeature_RNA < 4000
                         & nCount_RNA >= 500 & nCount_RNA < 15000
                         & mtpercent < 30)
HB17_PDX = subset(HB17_PDX, subset = nFeature_RNA >= 200 & nFeature_RNA < 8000 & 
                    nCount_RNA >= 500 & nCount_RNA < 25000)
HB17_tumor = subset(HB17_tumor, subset = nFeature_RNA >= 200 & nFeature_RNA < 8000 & 
                      nCount_RNA >= 500 & nCount_RNA < 24000 & mtpercent < 5)
HB30_PDX = subset(HB30_PDX, subset = nFeature_RNA >= 200 & nFeature_RNA < 9000 &
                    nCount_RNA >= 500 & nCount_RNA < 25000)
HB30_tumor = subset(HB30_tumor, subset = nFeature_RNA >= 200 & nFeature_RNA < 7500 &
                      nCount_RNA >= 500 & nCount_RNA < 25000)
HB53_background = subset(HB53_background, subset = nFeature_RNA >= 200 & nFeature_RNA < 7000 &
                           nCount_RNA >=500 & nCount_RNA < 45000)
HB53_tumor = subset(HB53_tumor, subset = nFeature_RNA >= 200 & nFeature_RNA < 9000 &
                      nCount_RNA >= 500 & nCount_RNA < 35000)
# Merge all samples into one combined Seurat object
merged.seurat = merge(HB17_background, y = c(HB17_tumor, HB17_PDX, 
                                             HB30_tumor, HB30_PDX,
                                             HB53_background, HB53_tumor),add.cell.ids = ls()[4:10])
# Save sample IDs into metadata
merged.seurat$sample = rownames(merged.seurat@meta.data)
# Split sample IDs into patient, type, and barcode metadata columns
merged.seurat@meta.data = separate(merged.seurat@meta.data, col =  "sample", into = c("patient", "Type", "Barcode"),
                                   
                                   sep = "_")
# View the final metadata table
view(merged.seurat@meta.data)
# Plot QC metrics for the merged dataset (violin plots)
VlnPlot(merged.seurat, features = c("nFeature_RNA", "nCount_RNA", "mtpercent"), ncol = 3)
# Feature scatter plot to visualize nCount_RNA vs nFeature_RNA
FeatureScatter(merged.seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = "lm")
# Apply additional QC filtering for the merged dataset
merged.seurat =  subset(merged.seurat,
                        subset = nFeature_RNA > 200 &
                          nFeature_RNA < 7000 &
                          nCount_RNA < 25000 &
                          mtpercent < 20)
# Normalize expression data
merged.seurat = NormalizeData(merged.seurat)
# Identify highly variable features
merged.seurat = FindVariableFeatures(merged.seurat)
# Scale the data
merged.seurat = ScaleData(merged.seurat)
# Perform PCA dimensionality reduction
merged.seurat = RunPCA(merged.seurat)
# Elbow plot to determine number of significant PCs
ElbowPlot(merged.seurat)
# Find nearest neighbors 
merged.seurat = FindNeighbors(merged.seurat, dims = 1:20)
# Perform clustering
merged.seurat = FindClusters(merged.seurat)
# Run UMAP for visualization
merged.seurat = RunUMAP(merged.seurat, dims = 1:20)
view(merged.seurat@meta.data)
# UMAP plots grouped by patient and by Type
p1 = DimPlot(merged.seurat, reduction = "umap", group.by = "patient")
p2 = DimPlot(merged.seurat, reduction = "umap", group.by = "Type", cols = c("red", "violet", "blue"))
# Arrange side by side
grid.arrange(p1, p2, ncol = 2, nrow = 2)
# Check RNA assay layers
names(merged.seurat@assays$RNA@layers)
# Set default assay to RNA
DefaultAssay(merged.seurat) = "RNA"
# Collapse layers into one 
merged.seurat = JoinLayers(merged.seurat, assay = "RNA")
# Split merged object into a list by patient for integration
obj.list = SplitObject(merged.seurat, split.by = "patient")
# Process each patient-specific object
for (i in seq_along(obj.list)) {
  # Subsample max 3000 cells per patient
  num_cells = min(3000, ncol(obj.list[[i]]))
  cells_use = sample(colnames(obj.list[[i]]), num_cells)
  obj.list[[i]] = subset(obj.list[[i]], cells = cells_use)
  # Normalize and find variable features
  obj.list[[i]] = NormalizeData(object = obj.list[[i]])
  obj.list[[i]] = FindVariableFeatures(object = obj.list[[i]], nfeatures = 1000)
  gc()
}
# Select integration features
features = SelectIntegrationFeatures(object.list = obj.list, nfeatures = 1000)
# Find integration anchors across patient datasets
anchors = FindIntegrationAnchors(object.list = obj.list,
                                 anchor.features = features,
                                 normalization.method = "LogNormalize")
# Integrate datasets into a single Seurat object
seurat.integrated = IntegrateData(anchorset = anchors)
# Set default assay to "integrated"
DefaultAssay(seurat.integrated) <- "integrated"
# Scale integrated data
seurat.integrated = ScaleData(object = seurat.integrated)
# Run PCA on integrated dataset
seurat.integrated = RunPCA(object = seurat.integrated)
# Run UMAP on integrated dataset
seurat.integrated = RunUMAP(object = seurat.integrated, dims = 1:20)
# UMAP plots grouped by patient and Type (after integration)
p3 = DimPlot(seurat.integrated, reduction = "umap", group.by = "patient")
p4 = DimPlot(seurat.integrated, reduction = "umap", group.by = "Type",
             cols = c("red", "violet", "blue"))
grid.arrange(p3, p4, ncol = 2)
# Graph-based clustering on integrated data
seurat.integrated = FindNeighbors(seurat.integrated, dims = 1:20)
seurat.integrated = FindClusters(seurat.integrated, resolution = 0.5)
# Visualize clusters on UMAP
DimPlot(seurat.integrated, reduction = "umap", label = T)
# Switch back to raw RNA assay for marker analysis
DefaultAssay(seurat.integrated) = "RNA"
seurat.integrated = JoinLayers(seurat.integrated)
# Find markers for all clusters 
markers = FindAllMarkers(seurat.integrated, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
# Select top 10 markers per cluster by log2FC
top10_markers = markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)
# Print and save markers
print(top10_markers)
write.csv(top10_markers, "top10_markers_per_clusters.csv", row.names = F)
view(top10_markers)