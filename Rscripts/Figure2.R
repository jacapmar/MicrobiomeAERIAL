# ====================================
# AERIAL Project – Figure 2
# ====================================

# Description: Generates Figure 2
# Author: Jose Caparros Martin, PhD
# Email: joseantonio.caparrosmartin@uwa.edu.au
# License: MIT
# Disclaimer: Provided "as is" for academic/research use. 
# Date: 2025-08-05

# ---------
# Libraries
# ---------
library(ggplot2)
library(RColorBrewer)
library(vegan)
library(rstatix)
library(tidyverse)
library(PairedData)

# ---------
# Create Folder directory for figures (within current working directory)
# ---------

folder_path <- "figures" 

if (!dir.exists(folder_path)) {
  dir.create(folder_path, recursive = TRUE)
  message("Directory created at: ", folder_path)
} else {
  message("Directory already exists: ", folder_path)
}

# ====================================
# [1] Load Input Data
# ====================================
data_final      <- read.csv("data/16S_data_final.csv", sep = ";")
taxonomy        <- read.csv("data/16S_taxonomy.csv", sep = ";")
metadata        <- read.csv("data/metadata_sharing.csv", sep = ";")
metadata_longi  <- read.csv("data/metadata_sharing_longitudinal.csv", sep = ";")
qpcr            <- read.csv("data/qpcr_data.csv", sep = ";")

# ====================================
# [2] Shannon Diversity Analysis
# ====================================
# Subset to longitudinal samples
data_longi <- data_final[, colnames(data_final) %in% metadata_longi$Seq_ID]

# Add Shannon index (computed across taxa per sample)
metadata_longi$shannon_index <- vegan::diversity(data_longi, index = "shannon", MARGIN = 2)

# Prepare baseline and symptomatic vectors
symptomatic<- metadata_longi %>%
  filter(Typeofsample == "Symptomatic") %>%
  pull(shannon_index)

baseline <- metadata_longi %>%
  filter(Typeofsample == "Baseline") %>%
  pull(shannon_index)


# -----------------------------------
# Plot Paired Profile (Shannon Index)
# -----------------------------------
paired_data <- paired(baseline, symptomatic)

pdf("figures/Fig2C.pdf", width = 5, height = 5)
plot(paired_data, type = "profile") + theme_bw()
dev.off()

# -------------------------------
# Statistical Testing
# -------------------------------

# Difference vector for paired test
shannon_diff <- baseline - symptomatic

# Normality test
shapiro_res <- shapiro.test(shannon_diff)
print(shapiro_res)

# Paired wilcox-test
wilcox.test(baseline, symptomatic, paired = TRUE, alternative = "two.sided")


# ====================================
# [3] qPCR Biomass Analysis Fig2D
# ====================================

# Subset qPCR data to longitudinal samples
qpcr_longi <- qpcr %>% filter(Seq_ID %in% metadata_longi$Seq_ID)

# Confirm correct ordering
stopifnot(all.equal(metadata_longi$Seq_ID, qpcr_longi$Seq_ID))

# Add qPCR values to metadata
metadata_longi$qPCR <- qpcr_longi$X1ct

# Separate baseline and symptomatic values
qpcr_symptomatic<- metadata_longi %>%
  filter(Typeofsample == "Symptomatic") %>%
  pull(qPCR)

qpcr_baseline <- metadata_longi %>%
  filter(Typeofsample == "Baseline") %>%
  pull(qPCR)

# Paired plot
qpcr_paired <- paired(qpcr_baseline, qpcr_symptomatic)

pdf("figures/Fig2D.pdf", width = 5, height = 5)
plot(qpcr_paired, type = "profile") + theme_bw()
dev.off()

# Statistical tests
qpcr_diff <- qpcr_baseline - qpcr_symptomatic
shapiro.test(qpcr_diff)

t_test_res <- t.test(qpcr_baseline, qpcr_symptomatic, paired = TRUE, alternative = "two.sided")
print(t_test_res)


# ==========================================
# [4] Cross-sectional Biomass Analysis Fig 2B
# ==========================================

# Subset qPCR data to cross-sectional samples
qpcr_cross <- qpcr %>% filter(Seq_ID %in% metadata$Seq_ID)

# Add qPCR values to metadata
metadata$qPCR <- qpcr_cross$X1ct

# Plot: Biomass by Sample Type
pdf("figures/Fig2B.pdf", width = 4, height = 4)
ggplot(metadata, aes(x = Typeofsample, y = qPCR)) +
  geom_boxplot(notch = FALSE, width = 0.5, lwd = 1, outlier.shape = NA) +
  geom_violin(fill = "blue", color = NA, alpha = 0.2) +
  geom_jitter(shape = 16, position = position_jitter(0.2), color = "red") +
  scale_color_brewer(palette = "Set1") +
  labs(title = "Biomass vs Sample Type", x = "Sample Type", y = "1/Ct") +
  theme_classic()
dev.off()

# Statistical analysis: Biomass
metadata %>%
  wilcox_test(qPCR ~ Typeofsample, p.adjust.method = "bonferroni")

metadata %>%
  group_by(Typeofsample) %>%
  get_summary_stats(qPCR, type = "full")

metadata %>%
  wilcox_effsize(qPCR ~ Typeofsample)


# =============================================
# [5] Cross-sectional Shannon Diversity Fig 2A
# =============================================

# Subset abundance data to cross-sectional samples
data_cross <- data_final[, colnames(data_final) %in% metadata$Seq_ID]

# Compute Shannon diversity and add to metadata
metadata$shannon_index <- vegan::diversity(data_cross, index = "shannon", MARGIN = 2)

# Plot: Shannon Index by Sample Type
pdf("figures/Fig2A.pdf", width = 4, height = 4)
ggplot(metadata, aes(x = Typeofsample, y = shannon_index)) +
  geom_boxplot(notch = FALSE, width = 0.5, lwd = 1, outlier.shape = NA) +
  geom_violin(fill = "blue", color = NA, alpha = 0.2) +
  geom_jitter(shape = 16, position = position_jitter(0.2), color = "red") +
  scale_color_brewer(palette = "Set1") +
  labs(title = "Shannon Diversity vs Sample Type", x = "Sample Type", y = "Shannon Index") +
  theme_classic()
dev.off()

# Statistical analysis: Shannon Index
metadata %>%
  wilcox_test(shannon_index ~ Typeofsample, p.adjust.method = "bonferroni")

metadata %>%
  group_by(Typeofsample) %>%
  get_summary_stats(shannon_index, type = "full")

metadata %>%
  wilcox_effsize(shannon_index ~ Typeofsample)

