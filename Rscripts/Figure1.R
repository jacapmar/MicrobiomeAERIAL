# ======================================
# AERIAL Project – Figure 1 and ANCOMBC
# ======================================
# Author: Jose Caparros Martin, PhD
# Email: joseantonio.caparrosmartin@uwa.edu.au
# License: MIT
# Disclaimer: Provided "as is" for academic/research use.
# Date: 2025-08-05

# --------- Libraries ---------
library(ggplot2)
library(RColorBrewer)
library(vegan)
library(rstatix)
library(tidyverse)
library(reshape2)
library(phyloseq)
library(ANCOMBC)
library(compositions)

# --------- Custom Functions ---------
source("Rscripts/filter_and_normalize.R")

# --------- Create Output Folder ---------
folder_path <- "figures"
if (!dir.exists(folder_path)) dir.create(folder_path, recursive = TRUE)

# --------- Load Input Data ---------
data_final      <- read.csv("data/16S_data_final.csv", sep = ";")
taxonomy        <- read.csv("data/16S_taxonomy.csv", sep = ";")
metadata_cross  <- read.csv("data/metadata_sharing.csv", sep = ";")
metadata_longi  <- read.csv("data/metadata_sharing_longitudinal.csv", sep = ";")
qpcr            <- read.csv("data/qpcr_data.csv", sep = ";")

# --------- Compute Shannon Diversity ---------
metadata_cross$shannon_index <- vegan::diversity(data_final[,1:71], index = "shannon", MARGIN = 2)

# --------- Normalize Data ---------
rownames(data_final) <- taxonomy[,5]
res <- filter_and_normalize(t(data_final[,1:71]), filter = FALSE, normalize = TRUE)
long_compositional <- t(res$data.normalized)

# --------- Subset for Symptomatic, Virus-positive Samples ---------
metadata_viruses <- metadata_cross %>%
  filter(RHV == "Yes" | SARSCOV2 == "Yes") %>%
  filter(Typeofsample == "Symptomatic") %>%
  mutate(
    RHV_SARSCOV2 = ifelse(RHV == "Yes" & SARSCOV2 == "Yes", "Yes", "No"),
    virus_type = case_when(
      RHV_SARSCOV2 == "Yes" ~ "Both",
      SARSCOV2 == "Yes"     ~ "SARS-CoV-2",
      RHV == "Yes"          ~ "RHV",
      TRUE                  ~ NA_character_
    )
  ) %>%
  filter(virus_type != "Both")

# --------- Plot: Shannon Diversity ---------
pdf("figures/Fig1_Shannon.pdf", width = 4, height = 4)
ggplot(metadata_viruses, aes(x = virus_type, y = shannon_index)) +
  geom_boxplot(width = 0.5, lwd = 1, outlier.shape = NA) +
  geom_violin(fill = 'blue', alpha = 0.2, color = NA) +
  geom_jitter(color = "red", shape = 16, position = position_jitter(0.2)) +
  scale_color_brewer(palette = "Set1") +
  labs(title = "Shannon Diversity vs Virus", x = "Virus Detected", y = "Shannon Index") +
  theme_classic()
dev.off()

# Statistical test
wilcox.test(shannon_index ~ virus_type, data = metadata_viruses)  # p ~ 0.021

# --------- Add qPCR Data and Plot ---------
qpcr_subset <- qpcr[qpcr$Seq_ID %in% metadata_viruses$Seq_ID, ]
stopifnot(all.equal(qpcr_subset$Seq_ID, metadata_viruses$Seq_ID))

metadata_viruses$qpcr <- qpcr_subset$X1ct

# Plot qPCR
pdf("figures/Fig1_qPCR.pdf", width = 4, height = 4)
ggplot(metadata_viruses, aes(x = virus_type, y = qpcr)) +
  geom_boxplot(width = 0.5, lwd = 1, outlier.shape = NA) +
  geom_violin(fill = 'blue', alpha = 0.2, color = NA) +
  geom_jitter(color = "red", shape = 16, position = position_jitter(0.2)) +
  scale_color_brewer(palette = "Set1") +
  labs(title = "Bacterial Load vs Virus", x = "Virus Detected", y = "1/Ct") +
  theme_classic()
dev.off()

# Statistical test
wilcox.test(qpcr ~ virus_type, data = metadata_viruses)  # p ~ 0.6673



# ======================================
# Virus–Bacteria Interaction Analysis
# ======================================

# Select virus data (Yes/No values)
virus <- metadata_cross[, c("RHV", "SARSCOV2")]

# Transform microbiome data to CLR (Compositional Log-Ratio)
microbiome_clr <- compositions::clr(long_compositional + 1)

# Select specific taxa of interest
taxon <- c("Moraxella catarrhalis", "Moraxella nonliquefaciens", "Streptococcus")
subset_taxon <- as.data.frame(t(microbiome_clr[rownames(microbiome_clr) %in% taxon, ]))

# Convert virus Yes/No to 1/0
virus_table <- as.data.frame(ifelse(virus == "Yes", 1, 0))

# Initialize correlation and p-value matrices
cor_matrix <- matrix(NA, nrow = ncol(subset_taxon), ncol = ncol(virus_table),
                     dimnames = list(colnames(subset_taxon), colnames(virus_table)))
pval_matrix <- cor_matrix

# Compute Spearman correlations
for (i in 1:ncol(subset_taxon)) {
  for (j in 1:ncol(virus_table)) {
    test <- cor.test(subset_taxon[, i], virus_table[, j], method = "spearman")
    cor_matrix[i, j] <- test$estimate
    pval_matrix[i, j] <- test$p.value
  }
}

# FDR correction
pvals_vector <- as.vector(pval_matrix)
fdr_corrected <- p.adjust(pvals_vector, method = "fdr")
pval_matrix_fdr <- matrix(fdr_corrected,
                          nrow = nrow(pval_matrix),
                          ncol = ncol(pval_matrix),
                          dimnames = dimnames(pval_matrix))

# Prepare for heatmap
cor_df <- reshape2::melt(cor_matrix, varnames = c("ASV", "Virus"), value.name = "correlation")
pval_df <- reshape2::melt(pval_matrix_fdr, varnames = c("ASV", "Virus"), value.name = "p_value")

merged_df <- merge(cor_df, pval_df, by = c("ASV", "Virus"))
merged_df$p_label <- ifelse(merged_df$p_value < 0.001, "< 0.001", sprintf("%.3f", merged_df$p_value))

# Plot heatmap with correlations and FDR-corrected p-values
pdf("figures/Fig1_Correlation_Heatmap.pdf", width = 6, height = 6)
ggplot(merged_df, aes(x = Virus, y = ASV, fill = correlation)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0, limit = c(-1, 1), name = "Spearman\ncorrelation") +
  geom_text(aes(label = p_label), size = 3) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Virus–Bacteria Correlation Heatmap",
       x = "Virus", y = "ASV")
dev.off()



# ======================================
# Cumulative Abundance of Top Taxa
# ======================================

# Convert to relative abundance
microbiome_rel <- t(long_compositional)
microbiome_rel1 <- microbiome_rel[, order(colSums(microbiome_rel), decreasing = TRUE)]
# Use top 3 taxa
taxa_of_interest <- colnames(microbiome_rel1)[1:3]
# Total abundance per taxon across all samples (original absolute data)
total_abundance_all <- rowSums(data_final[, 1:71])
# Total abundance of selected top taxa
total_abundance_target <- sum(total_abundance_all[taxa_of_interest])
# Overall abundance across all taxa
total_abundance_overall <- sum(total_abundance_all)
# Calculate cumulative abundance percentage
cumulative_percentage <- (total_abundance_target / total_abundance_overall) * 100
message(paste("Cumulative abundance of top 3 taxa:", round(cumulative_percentage, 2), "%"))

# Barplot of Top 10ASVs


microbiota <- t(data_final[, 1:71])
n_asvs <- 10
top_asvs <- colnames(microbiome_rel1)[1:n_asvs]
microbiota_top <- microbiota[, top_asvs]
total_abundance_all <- sum(microbiota)
cumulative_percentages <- colSums(microbiota_top) / total_abundance_all * 100

sorted_indices <- order(cumulative_percentages, decreasing = TRUE)
cumulative_percentages <- cumulative_percentages[sorted_indices]
top_asvs <- top_asvs[sorted_indices]

pdf("figures/Fig1_Cumulative_TopASV.pdf", width = 6, height = 6)
barplot(cumulative_percentages,
        names.arg = top_asvs,
        las = 2,
        col = "steelblue",
        ylab = "Cumulative Abundance (%)",
        main = paste("Top", n_asvs, "ASVs by Cumulative Abundance"),
        cex.names = 0.7)
dev.off()





# ======================================
# Run ANCOM-BC2
# ======================================

# Define virus columns
virus_cols <- c("RHV", "SARSCOV2")  # replace with your actual column names
# Convert 'Yes'/'No' to TRUE/FALSE
metadata_cross[virus_cols] <- lapply(metadata_cross[virus_cols], function(x) x == "Yes") # convert into logit
# Count number of viruses per sample
metadata_cross$virus_count <- rowSums(metadata_cross[virus_cols])
# Keep only samples with 0 or 1 viruses
single_virus_metadata <- metadata_cross[metadata_cross$virus_count <= 1, ]


## create otu table for phyloseq OTU by Sample
community_virus <- data_final[,which(colnames(data_final) %in% single_virus_metadata$Seq_ID)]
x <- t(community_virus)
otumat <- as.matrix(t(x))
OTU = otu_table(otumat, taxa_are_rows = TRUE)

## ----------------------------- Check point
all.equal(rownames(otumat), taxonomy$Species) ## TRUE
## ------------------------------

rownames(taxonomy) <- taxonomy$Species
taxamat <- as.matrix(taxonomy)
TAX <- tax_table(taxamat)


## create sample data file
sampledata <- single_virus_metadata
rownames(sampledata) <- sampledata$Seq_ID

## ----------------------------- Check point
all.equal(colnames(otumat),rownames(sampledata)) ##TRUE
## ------------------------------

sampleinformation <- sample_data(sampledata)


## create phyloseq file and save

physeq <- phyloseq(OTU,TAX,sampleinformation)
microbiota_physeq <- saveRDS(physeq, file = "data/aerial_phyloseq.rds")

# Run ANCOM-BC2 for RHV DA
res_ancombc2 <- ancombc2(
  physeq,        
  tax_level= "Species",
  rand_formula = NULL, 
  fix_formula = "RHV",
  p_adj_method = "holm",
  pseudo = 0.01,
  pseudo_sens = TRUE,
  prv_cut = 0.10,
  lib_cut = 1000,
  s0_perc = 0.05,
  group = "RHV",                 
  struc_zero = TRUE,
  neg_lb = TRUE,
  alpha = 0.05,
  global = TRUE
)


# Access results
res_df <- res_ancombc2$res  

# Compute 95% CI for RHV effect

res_df$CI_lower_RHV <- res_df$lfc_RHVTRUE - 1.96 * res_df$se_RHVTRUE
res_df$CI_upper_RHV <- res_df$lfc_RHVTRUE + 1.96 * res_df$se_RHVTRUE

# Convert log fold changes and CI to fold change scale
res_df$fold_change_RHV <- exp(res_df$lfc_RHVTRUE)
res_df$CI_lower_RHV_FC <- exp(res_df$CI_lower_RHV)
res_df$CI_upper_RHV_FC <- exp(res_df$CI_upper_RHV)

# Select significant results
res_df %>%
  filter(diff_RHVTRUE == TRUE)





