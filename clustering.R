# Define required packages
req_packages <- c("Seurat", "SeuratObject", "ggplot2", "dplyr", "tidyverse", "umap", "patchwork")
installed_packages <- find.packages(req_packages)

# Find missing packages
missing_packages <- setdiff(req_packages, installed_packages)

# Install missing packages
if (length(missing_packages) > 0) {
    message(paste("Installing missing packages:", paste(missing_packages, collapse = ", "), sep = " "))
    install.packages(missing_packages)
} else {
    message("All required packages are already installed.")
}
# Library setup
library(Seurat)
library(SeuratObject)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(umap)
library(patchwork)

# Take metadata file from command line
metadata <- readline("Enter the path to your metadata file: ")
if (!file.exists(metadata)) {
    stop(paste("Error:", metadata, "does not exist.", sep = " "))
}
metadata_df <- read.csv(metadata, sep = ",", header = TRUE)
message(paste("Metadata file loaded successfully from:", metadata, sep = " "))

# Generate Seurat object for each individual sample dataset
i <- 1
while (i <= length(metadata_df$file_path)) {
    cts <- Read10X_h5(filename = metadata_df$file_path[i])
    sample_name <- metadata_df$sample_name[i]

    # Create Seurat object
    seurat <- CreateSeuratObject(counts = cts, min.cells = 3, min.features = 200)
    seurat@meta.data$orig.ident <- plyr::mapvalues(seurat@meta.data$orig.ident,
                                    from = "SeuratProject", to = sample_name)
    assign(sample_name, seurat)
    i <- i + 1
} ## TESTED UNTIL THIS POINT
# Merge Seurat objects together
merged_seurat <- merge()

# Add metadata to the merged Seurat object
sample_labels <- sample(x = metadata_df$sample_name,
        size = ncol(x = merged_seurat), replace = TRUE)
merged_seurat$sample <- sample_labels
age_labels <- sample(x = metadata_df$age,
        size = ncol(merged_seurat), replace = TRUE)
merged_seurat$age <- age_labels
sex_labels <- sample(x = metadata_df$sex,
        size = ncol(merged_seurat), replace = TRUE)
merged_seurat$sex <- sex_labels
treamtent_labels <- sample(x = metadata_df$treatment,
        size = ncol(merged_seurat), replace = TRUE)
merged_seurat <- treamtent_labels

# remove unnecessary variables from the environment
remove("cts", "i", "seurat", "metadata", "metadata_df")

# Create output directory for QC
directory_path <- file.path("output", "QC_results")
if (!dir.exists(directory_path)) {
    dir.create(directory_path)
}
# Examine your data for mitochondrial and ribosomal RNA content
merged_seurat$sample <- rownames(merged_seurat@meta.data)
merged_seurat[["percent.mt"]] <- PercentageFeatureSet(merged_seurat, pattern = "^MT-")
merged_seurat[["percent.rb"]] <- PercentageFeatureSet(merged_seurat, pattern = "^RP[SL]")

# Visualize QC metrics
pdf(file = "output/QC_Results/QCmetrics.pdf")
par(mfrow = c(3,2))
VlnPlot(merged_seurat, group.by = "orig.ident", 
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 2)
FeatureScatter(merged_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(merged_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()

# Filter and normalize the data
merged_filtered <- subset(merged_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 5)
merged_filtered <- NormalizeData(merged_filtered, normalization.method = "LogNormalize", scale.factor = 10000)

# Create output directory for results
directory_path <- file.path("output", "Clustering_Results")
if (!dir.exists(directory_path)) {
    dir.create(directory_path)
}

# Identify highly variable features
merged_filtered <- FindVariableFeatures(merged_filtered, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(merged_filtered))
pdf(file = "output/Clustering_Results/variable_features.pdf")
par(mfrow = c(1,2))
plot1 <- VariableFeaturesPlot(merged_filtered)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1
plot2
dev.off()

# Scaling the data
all.genes <- rownames(merged_filtered)
merged_filtered <- ScaleData(merged_filtered, features = all.genes)

# Perform linear dimensional reduction
merged_filtered <- RunPCA(merged_filtered, 
    features = VariableFeatures(object = merged_filtered))

pdf(file = "output/Clustering_Results/pca.pdf")
print(merged_filtered[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(merged_filtered, dims = 1:15, reduction = "pca")
DimPlot(merged_filtered, reduction = "pca") + NoLegend()
DimHeatmap(merged_filtered, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(merged_filtered, dims = 1:15, cells = 500, balanced = TRUE)
Elbowplot(merged_filtered)
dev.off()

# Cluster the cells
merged_filtered <- FindNeighbors(merged_filtered, dims = 1:13)
merged_filtered <- FindClusters(merged_filtered, resolution = 0.5)

# Run non-linear dimensional reduction (UMAP)
merged_umap <- RunUMAP(merged_filtered, dims = 1:12)
joined_umap <- JoinLayers(object = merged_umap, assay = "RNA")
remove(merged_umap, merged_filtered, merged_seurat)

pdf()

dev.off()