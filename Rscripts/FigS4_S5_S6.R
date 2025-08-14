# ====================================
# AERIAL Project – Figures S4, S5, S6
# ====================================

# Description: Generates Figures S4, S5, S6
# Author: Jose Caparros Martin, PhD
# Email: joseantonio.caparrosmartin@uwa.edu.au
# License: MIT
# Date: 2025-08-05

# ===============================
# Libraries
# ===============================
library(ggplot2)
library(RColorBrewer)
library(mixOmics)
library(vegan)
library(rstatix)
library(tidyverse)

# ===============================
# Create Folder directory for figures (within current working directory)
# ===============================

folder_path <- "figures" 

if (!dir.exists(folder_path)) {
  dir.create(folder_path, recursive = TRUE)
  message("Directory created at: ", folder_path)
} else {
  message("Directory already exists: ", folder_path)
}


# ===============================
# [1] Load custom function
# ===============================
source("Rscripts/filter_and_normalize.R")

# ===============================
# [2] Load input data
# ===============================
data_final <- read.csv("data/16S_data_final.csv", sep = ";")
taxonomy <- read.csv("data/16S_taxonomy.csv", sep = ";")
metadata <- read.csv("data/metadata_sharing.csv", sep = ";")

# Set row names from taxonomy 
rownames(data_final) <- taxonomy[,5]

# ===============================
# [3] Normalize data
# ===============================
res <- filter_and_normalize(t(data_final[,1:71]), filter = FALSE, normalize = TRUE)
long_compositional <- res$data.normalized

# ===============================
# [4] Reorder taxa and select top N
# ===============================
x <- long_compositional[, order(colSums(long_compositional), decreasing = TRUE)]

# Reorder to include Haemophilus influenzae manually
x <- x[, c(1:11, 14, 12:13, 15:133)]

# Select top N taxa
N <- 12
taxa_list <- colnames(x)[1:N]

# Create matrix with top taxa and "Others"
new_z <- data.frame(
  x[, colnames(x) %in% taxa_list],
  a_Others = rowSums(x[, !colnames(x) %in% taxa_list])
)

# ===============================
# [5] Match metadata to samples
# ===============================
vectortomatch <- match(rownames(new_z), metadata$Seq_ID)
treat.normal <- metadata$Typeofsample
treat.normal <- treat.normal[vectortomatch]

group.normal <- data.frame(row.names = rownames(new_z), treat.normal)

# ===============================
# [6] Long format dataframe
# ===============================
df <- NULL
for (i in 1:dim(new_z)[2]) {
  tmp <- data.frame(
    row.names = NULL,
    Sample = rownames(new_z),
    Taxa = rep(colnames(new_z)[i], dim(new_z)[1]),
    Value = new_z[, i],
    Type = group.normal[, 1]
  )
  if (i == 1) {
    df <- tmp
  } else {
    df <- rbind(df, tmp)
  }
}

# ===============================
# [7] Color palette
# ===============================
set3_colors <- brewer.pal(12, "Set3")
extra_color <- "#FF6347"  # Tomato red
extended_colors <- c(set3_colors, extra_color)

# ===============================
# [8] Barplot (Fig S4)
# ===============================
pdf("figures/FigS4.pdf", 14, 14)

p <- ggplot(df, aes(Sample, Value, fill = Taxa)) +
  geom_bar(stat = "identity") +
  facet_grid(. ~ Type, drop = TRUE, scale = "free", space = "free_x") +
  theme_bw() +
  ylab("Proportions") +
  scale_fill_manual(values = extended_colors) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(
    strip.background = element_rect(fill = "gray85"),
    panel.spacing = unit(0.3, "lines"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

print(p)
dev.off()


# ===============================
# [9] Boxplots FigS5A-B
# ===============================
## isolate moraxella data [2] and streptococcus pneummonia [9]
moraxella <- df[df$Taxa == unique(df$Taxa)[2], ]
moraxella %>% wilcox_effsize(Value~Type) ## effect size 0.01
moraxella %>% group_by(Type) %>% get_summary_stats(Value,type = "full") 
moraxella %>% wilcox_test(Value~Type, p.adjust.method="bonferroni") ## p value 0.373


pdf("figures/FigS5B.pdf",4,4)
ggplot(data=moraxella, aes(x=Type, y=Value)) +
  geom_boxplot(notch = F, width=0.5,lwd=1,outlier.shape = NA) +
  geom_violin(fill='blue', color=NA,alpha=0.2) +
  geom_jitter(shape=16, position=position_jitter(0.2),color="red") + 
  scale_color_brewer(palette="Set1") + 
  labs(title="Abundance Moraxella", x="Sample Type", y="% moraxella catharallis") + 
  theme_classic()
dev.off()


staph <- df[df$Taxa == unique(df$Taxa)[9], ]
staph %>% wilcox_effsize(Value~Type) ## effect size 0.09
staph %>% group_by(Type) %>% get_summary_stats(Value,type = "full") 
staph %>% wilcox_test(Value~Type, p.adjust.method="bonferroni") ## p value 0.024


pdf("figures/FigS5A.pdf",4,4)
ggplot(data=staph, aes(x=Type, y=Value)) +
  geom_boxplot(notch = F, width=0.5,lwd=1,outlier.shape = NA) +
  geom_violin(fill='blue', color=NA,alpha=0.2) +
  geom_jitter(shape=16, position=position_jitter(0.2),color="red") + 
  scale_color_brewer(palette="Set1")  + 
  labs(title="Abundance Strepto", x="Sample Type", y="% Streptococcus pneumoniae") +
  theme_classic()
dev.off()

pvalues <- c(0.319,0.0218)
correctedpvalues <- p.adjust(pvalues,method = "fdr")
correctedpvalues #0.3190 0.0436



# ===============================
# [10] PCA (Fig S6) and PERMANOVA
# ===============================

# ----- PCA analysis -----
all_samples_forpca <- long_compositional + 1
pca_forallsamples <- mixOmics::pca(all_samples_forpca, ncomp = 5, logratio = "CLR")

Y <- factor(metadata$Typeofsample, c("Baseline", "Symptomatic"))

pdf("figures/FigS6.pdf", 6, 6)
mixOmics::plotIndiv(
  pca_forallsamples,
  comp = c(1, 3),
  ind.names = FALSE,
  title = "PCA comp 1 - 3",
  legend = TRUE,
  ellipse = TRUE,
  group = Y,
  col.per.group = color.mixo(c(1:2)),
  legend.title = "Sample type"
)
dev.off()

# ----- PERMANOVA: Type of sample -----
Y <- metadata$Typeofsample
adonis2((long_compositional + 1) ~ Y, method = "aitchison", permutations = 9999)

# ----- PERMANOVA: Subject -----
Y <- metadata$Infant_number
adonis2((long_compositional + 1) ~ Y, method = "aitchison", permutations = 9999)

# ----- PERMANOVA: Gender -----
Y <- metadata$gender
adonis2((long_compositional + 1) ~ Y, method = "aitchison", permutations = 9999)
