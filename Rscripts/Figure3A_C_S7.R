# ========================================================
# AERIAL Microbiome Manuscript – Figure S7, Figure 3 A-C
# Dirichlet Multinomial Mixture Clustering
# ========================================================

# Author: Jose Caparros Martin, PhD
# Email: joseantonio.caparrosmartin@uwa.edu.au
# License: MIT
# Disclaimer: Provided "as is" for academic/research use.
# Date: 2025-08-05

# --------- Libraries ---------
library(phyloseq)
library(DirichletMultinomial)
library(reshape2)
library(ggplot2)
library(dplyr)
library(microbiome)
library(mixOmics)

# --------- Custom Functions ---------
source("Rscripts/filter_and_normalize.R")

# --------- Load Input Data ---------
data_final      <- read.csv("data/16S_data_final.csv", sep = ";")
taxonomy        <- read.csv("data/16S_taxonomy.csv", sep = ";")
metadata_cross  <- read.csv("data/metadata_sharing.csv", sep = ";")
metadata_cross2 <- read.csv("data/metadata_cross_shann_dirichlet.csv", sep = ";")


# --------- Set Output Folder ---------
folder_path <- "figures"
if (!dir.exists(folder_path)) dir.create(folder_path, recursive = TRUE)

# --------- Prepare Phyloseq Object ---------
# OTU Table
otumat <- as.matrix(data_final[, 1:71])
rownames(otumat) <- taxonomy$Species
OTU <- otu_table(otumat, taxa_are_rows = TRUE)

# Taxonomy Table
rownames(taxonomy) <- taxonomy$Species
TAX <- tax_table(as.matrix(taxonomy))

# Sample Metadata
metadata_dirichlet <- metadata_cross[metadata_cross$Seq_ID %in% colnames(data_final[, 1:71]), ]
rownames(metadata_dirichlet) <- metadata_dirichlet$Seq_ID
SAMPLE <- sample_data(metadata_dirichlet)

# Combine into phyloseq object
physeq <- phyloseq(OTU, TAX, SAMPLE)

# --------- Dirichlet Multinomial Clustering ---------
set.seed(1492)
count_matrix <- as.matrix(t(abundances(physeq)))

# Fit DirichletMultinomial models for 1–10 components
fit_list <- lapply(1:10, dmn, count = count_matrix, verbose = TRUE)

# Evaluate model fit
laplace_scores <- sapply(fit_list, laplace)
aic_scores <- sapply(fit_list, AIC)
bic_scores <- sapply(fit_list, BIC)

# Select best-fit model by Laplace
best_model <- fit_list[[which.min(laplace_scores)]]
cluster_assignments <- apply(mixture(best_model), 1, which.max)

# --------- Helper: Plot Cluster Drivers ---------
plot_cluster_drivers <- function(cluster_number, output_file) {
  df <- melt(fitted(best_model))
  colnames(df) <- c("OTU", "cluster", "value")
  
  df <- df %>%
    filter(cluster == cluster_number) %>%
    arrange(value) %>%
    mutate(OTU = factor(OTU, levels = unique(OTU))) %>%
    filter(abs(value) > quantile(abs(value), 0.9))  # Top 10% drivers
  
  p <- ggplot(df, aes(x = OTU, y = value)) +
    geom_bar(stat = "identity", fill = "red") +
    coord_flip() +
    labs(title = paste("Drivers Cluster", cluster_number)) +
    theme_classic()
  
  pdf(output_file, width = 6, height = 6)
  print(p)
  dev.off()
}

# --------- Generate Driver Plots ---------
plot_cluster_drivers(1, "figures/FigS7A.pdf")
plot_cluster_drivers(2, "figures/FigS7B.pdf")

# --------- Heatmap of Cluster Drivers ---------
# Optional: Normalize data (TSS or CLR as preferred)
filtered <- filter_and_normalize(count_matrix, filter = FALSE, normalize = TRUE)$data.normalized

pdf("figures/Fig3B_Heatmap.pdf", width = 10, height = 10)
heatmapdmn(filtered, fit_list[[1]], best_model, ntaxa = 10, col = hcl.colors(20, "Blues", rev = TRUE))
dev.off()


# --------- PCA vs Dirichlet Clusters (CLR-normalized data) ---------

# Use filtered + normalized data 
res <- filter_and_normalize(t(data_final[, 1:71]), filter = FALSE, normalize = TRUE)$data.normalized
x_for_pca <- res + 1  # Offset to avoid log(0)

# Run PCA
x_pca <- mixOmics::pca(x_for_pca, ncomp = 5, logratio = "CLR")

# Dirichlet cluster assignments from earlier
Y <- as.character(cluster_assignments)

# Plot PCA (Components 1 vs 3) colored by Dirichlet clusters
pdf("figures/Fig3A_PCA_DirichletClusters.pdf", width = 6, height = 6)
plotIndiv(
  x_pca,
  comp = c(1, 3),
  ind.names = FALSE,
  title = "PCA – CLR vs Dirichlet Clusters (Comp 1 & 3)",
  legend = TRUE,
  ellipse = TRUE,
  group = Y,
  col.per.group = color.mixo(1:2),
  legend.title = "Dirichlet Cluster"
)
dev.off()


# ------Shannon Diversity vs Dirichlet Clusters (Endotypes)-------
# ======================================

# Plot Shannon index by Dirichlet-inferred clusters (stored in Endotype column)
pdfname <- "figures/Fig3C_Shannon_DirichletClusters.pdf"
pdf(pdfname, width = 4, height = 4)
ggplot(metadata_cross2, aes(x = Endotype_label, y = shannon_index)) +
  geom_boxplot(notch = FALSE, width = 0.5, lwd = 1, outlier.shape = NA) +
  geom_violin(fill = 'blue', color = NA, alpha = 0.2) +
  geom_jitter(shape = 16, position = position_jitter(width = 0.2, height = 0), color = "red") +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "Shannon Diversity vs Dirichlet Cluster",
    x = "Dirichlet-Inferred Cluster",
    y = "Shannon Index"
  ) +
  theme_classic()
dev.off()

# ---------- Statistical Comparison ----------
# Wilcoxon rank-sum test
shannon_test <- metadata_cross2 %>%
  wilcox_test(shannon_index ~ Endotype_label, p.adjust.method = "bonferroni")

# Effect size
shannon_effsize <- metadata_cross2 %>%
  wilcox_effsize(shannon_index ~ Endotype_label)

# Summary stats per cluster
shannon_summary <- metadata_cross2 %>%
  group_by(Endotype_label) %>%
  get_summary_stats(shannon_index, type = "full")

# Print results to console
print(shannon_test)
print(shannon_effsize)
print(shannon_summary)