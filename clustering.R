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

# Take path to CellRanger .h5 files from command line
# path2files <- readline("Enter the path to the directory containing your count files: ")
# files <- list.files(path = path2files)
# dataset <- c()

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
    seurat@meta.data$orig.ident <- plyr::mapvalues(seurat@meta.data$orig.ident, from = "SeuratProject", to = sample_name)
    assign(sample_name, seurat)
    i <- i + 1
} ## TESTED UNTIL THIS POINT
# Merge Seurat objects together
merged_seurat <- merge()

# Add metadata to the merged Seurat object
sample_labels <- sample(x = metadata_df$sample_name, size = ncol(x = merged_seurat), replace = TRUE)
merged_seurat$sample <- sample_labels
age_labels <- sample(x = metadata_df$age, size = ncol(merged_seurat), replace = TRUE)
merged_seurat$age <- age_labels
sex_labels <- sample(x = metadata_df$sex, size = ncol(merged_seurat), replace = TRUE)
merged_seurat$sex <- sex_labels
treamtent_labels <- sample(x = metadata_df$treatment, size = ncol(merged_seurat), replace = TRUE)
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