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

library(Seurat)
library(SeuratObject)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(umap)
library(patchwork)

# Take path to CellRanger .h5 files from command line
path2files <-
files <- list.files(path = )
Dataset <- c()

# Generate Seurat object for each individual sample dataset
for (x in files) {
    cts <- Read10X_h5(filename = paste0(path2files, x))
    x <- str_sub(x, start = 1, end = 3)
    Dataset <- append(Dataset, x)

    #create Seurat object
    assign(x, CreateSeuratObject(counts = cts, min.cells = 3, min.features = 200))
}