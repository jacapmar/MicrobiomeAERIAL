# ====================================
# AERIAL Microbiome Manuscript – Figure S14
# ====================================
# Description: Alluvial plot and FigS10
# Author: Jose Caparros Martin, PhD
# Email: joseantonio.caparrosmartin@uwa.edu.au
# License: MIT
# Date: 2025-08-05

# -------------------------
# Load Required Libraries
# -------------------------
library(ggalluvial)
library(dplyr)
library(forcats)
library(ggplot2)

# --------- Load Input Data ---------
metadata_cross2 <- read.csv("data/metadata_cross_shann_dirichlet.csv", sep = ";")
metadata_longi  <- read.csv("data/metadata_sharing_longitudinal.csv", sep = ";")

# --------- Prepare Data ---------


metadata_alluvial <- metadata_cross2[which(metadata_cross2$Seq_ID %in%metadata_longi$Seq_ID),]

df_alluvial <- data.frame(
  ID = metadata_alluvial$Infant_number,
  Timepoint = metadata_alluvial$Typeofsample,
  Cluster = metadata_alluvial$Endotype_label,
  nwheezing = metadata_alluvial$nwheezing,
  nsymptom = metadata_alluvial$Total_symptomatic
)


df_alluvial$Timepoint <- factor(df_alluvial$Timepoint, levels = c("Baseline", "Symptomatic")) 

#reorder

df_alluvial$ID <- factor(df_alluvial$ID, levels = unique(df_alluvial$ID))

# --------- Alluvial plot ---------
p_alluvial <- ggplot(
  df_alluvial,
  aes(x = Timepoint, stratum = Cluster, alluvium = ID, fill = Cluster, label = Cluster)
) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback", color = "gray") +
  geom_stratum(alpha = 0.8) +
  geom_text(stat = "alluvium", aes(label = ID), size = 3, nudge_x = 0.05, check_overlap = TRUE) +
  theme_minimal() +
  labs(title = "Cluster Transitions", x = "Timepoint", y = "Participants") +
  scale_fill_brewer(palette = "Set3")

ggsave(
  filename = "figures/alluvial_plot.pdf", 
  plot = p_alluvial,
  width = 8,
  height = 8
)

## Extract the ordering actually used for labels at the Baseline axis
gb <- ggplot_build(p_alluvial)

# The geom_text(stat="alluvium") layer typically holds the ID labels and their y-positions
# Layer indices: [[1]] flow, [[2]] stratum, [[3]] text (ID labels)
labdat <- gb$data[[3]]

# Find the x-position corresponding to "Baseline"
# (If Timepoint is a factor, Baseline is usually x == 1. We compute it to be safe.)
x_baseline <- which(levels(df_alluvial$Timepoint) == "Baseline")

# Keep only Baseline labels and order them by their plotted y position (top-to-bottom)
baseline_ids_in_plot_order <- labdat %>%
  filter(round(x) == x_baseline) %>%
  arrange(y) %>%
  pull(label)

# Defensive: drop NAs/duplicates while preserving order
baseline_ids_in_plot_order <- baseline_ids_in_plot_order[!is.na(baseline_ids_in_plot_order)]
baseline_ids_in_plot_order <- unique(baseline_ids_in_plot_order)

# Apply that order to df_baseline$ID and plot bars

df_baseline <- df_alluvial %>% filter(Timepoint == "Baseline") %>% distinct(ID, .keep_all = TRUE)

df_baseline <- df_baseline %>%
  mutate(ID = factor(ID, levels = baseline_ids_in_plot_order))

# --------- Wheezing plot (same ID order) ---------
pdf("figures/alluvialplot_wheezing.pdf", 8, 8)
ggplot(df_baseline, aes(x = ID, y = nwheezing)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Wheezing Episodes", x = NULL, y = "Count") +
  theme(axis.text.y = element_text(size = 8),
        axis.ticks.y = element_blank())
dev.off()

# --------- Symptom plot (same ID order) ---------
pdf("figures/alluvialplot_sympto.pdf", 8, 8)
ggplot(df_baseline, aes(x = ID, y = nsymptom)) +
  geom_bar(stat = "identity", fill = "coral") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Symptom Episodes", x = NULL, y = "Count") +
  theme(axis.text.y =  element_text(size = 8),
        axis.ticks.y = element_blank())
dev.off()


# --------- Stats ---------
# Safety: ensure the indices exist given current data size
n_baseline <- nrow(df_baseline)
req_idx <- c(2,4,7,10,12,14,15,16,19,
             3,5,8,9,13,20,22,
             1,6,11,17,18,21)
if (any(req_idx > n_baseline)) {
  warning("Some hard-coded indices exceed df_baseline rows — adjust indices or group logic.")
}

g1_idx <- c(2,4,7,10,12,14,15,16,19)
g2_idx <- c(3,5,8,9,13,20,22)
g3_idx <- c(1,6,11,17,18,21)

g1 <- df_baseline[g1_idx[g1_idx <= n_baseline], ]
g2 <- df_baseline[g2_idx[g2_idx <= n_baseline], ]
g3 <- df_baseline[g3_idx[g3_idx <= n_baseline], ]

# Wilcoxon tests (report p-values)
w_p_12_wz <- tryCatch(wilcox.test(g1$nwheezing, g2$nwheezing)$p.value, error = function(e) NA_real_)
w_p_13_wz <- tryCatch(wilcox.test(g1$nwheezing, g3$nwheezing)$p.value, error = function(e) NA_real_)
w_p_12_sy <- tryCatch(wilcox.test(g1$nsymptom,  g2$nsymptom)$p.value,  error = function(e) NA_real_)
w_p_13_sy <- tryCatch(wilcox.test(g1$nsymptom,  g3$nsymptom)$p.value,  error = function(e) NA_real_)

cat("\nWilcoxon p-values:\n",
    "g1 vs g2 (nwheezing): ", w_p_12_wz, "\n",
    "g1 vs g3 (nwheezing): ", w_p_13_wz, "\n",
    "g1 vs g2 (nsymptom):  ", w_p_12_sy, "\n",
    "g1 vs g3 (nsymptom):  ", w_p_13_sy, "\n", sep = "")




