# ======================================
# AERIAL Project – Figure 3D
# ======================================
# Author: Jose Caparros Martin, PhD
# Email: joseantonio.caparrosmartin@uwa.edu.au
# License: MIT
# Disclaimer: Provided "as is" for academic/research use.
# Date: 2025-08-05

# --------- Libraries ---------
library(dplyr)
library(ggplot2)
library(emmeans)


# --------- Load Input Data ---------
metadata_cross2 <- read.csv("data/metadata_cross_shann_dirichlet.csv", sep = ";")

# Remove placeholder/missing values
df <- metadata_cross2 %>%
  filter(ever_wheeze_age != 99.0)

# ----------------------------------------
# Fit linear model
# ----------------------------------------

model <- lm(ever_wheeze_age ~ Endotype_label, data = df)
summary(model)

# ----------------------------------------
# Estimated Marginal Means (EMMs)
# ----------------------------------------

em <- emmeans(model, ~ Endotype_label)
em_df <- as.data.frame(em)

# ----------------------------------------
# Plotting Marginal Means
# ----------------------------------------

pdf("figures/Fig3D.pdf", width = 8, height = 8)

ggplot(em_df, aes(x = Endotype_label, y = emmean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
  labs(
    title = "Marginal Effects of Cluster on Wheeze Onset Age",
    x = "Cluster",
    y = "Predicted Age of First Wheeze"
  ) +
  theme_minimal()

dev.off()