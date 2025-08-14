
# ====================================
# AERIAL Microbiome Manuscript – Figure S3
# ====================================
# Description: Procrustes analysis between raw and filtered microbiome datasets
# Author: Jose Caparros Martin, PhD
# Email: joseantonio.caparrosmartin@uwa.edu.au
# License: MIT
# Date: 2025-08-05

# -------------------------
# Load Required Libraries
# -------------------------
library(vegan)         
library(tidyverse)     
library(compositions)  
library(mixOmics)

# -------------------------
# Create Output Folder
# -------------------------
folder_path <- "figures"
if (!dir.exists(folder_path)) dir.create(folder_path, recursive = TRUE)

# -------------------------
# Load Custom Functions
# -------------------------
source("Rscripts/filter_and_normalize.R")

# -------------------------
# Load Input Data
# -------------------------
data_final <- read.csv("data/16S_data_final.csv", sep = ";")
raw_data   <- read.csv("data/Aerial_raw.csv", sep = ";")

# -------------------------
# Normalize and Transform Data
# -------------------------
filtered <- filter_and_normalize(t(data_final[, 1:71]), filter = FALSE, normalize = TRUE)$data.normalized
raw      <- filter_and_normalize(t(raw_data[, 1:71]),   filter = FALSE, normalize = TRUE)$data.normalized

filtered_clr <- logratio.transfo(as.matrix(filtered), logratio = "CLR", offset = 1)
raw_clr      <- logratio.transfo(as.matrix(raw),      logratio = "CLR", offset = 1)

# -------------------------
# PCA Analysis
# -------------------------
filtered_clr_pca <- prcomp(filtered_clr)
raw_clr_pca      <- prcomp(raw_clr)

# -------------------------
# Procrustes Analysis
# -------------------------
microbiota_pro <- procrustes(raw_clr_pca, filtered_clr_pca, symmetric = FALSE)

# Save Procrustes Plot
pdf("figures/FigS3.pdf", width = 4, height = 4)
plot(microbiota_pro)
dev.off()

# Perform PROTEST (permutation test of Procrustes correlation)
protest_result <- protest(X = raw_clr_pca, Y = filtered_clr_pca,
                          scores = "sites", permutations = 10000)

print(protest_result)


