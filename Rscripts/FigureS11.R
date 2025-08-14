# ======================================
# AERIAL Project – Figure S11
# ======================================
# Author: Jose Caparros Martin, PhD
# Email: joseantonio.caparrosmartin@uwa.edu.au
# License: MIT
# Date: 2025-08-05


# --------- Load Required Libraries ---------
library(MASS)
library(ggeffects)
library(AER)
library(ggplot2)
library(dplyr)
library(emmeans)

# --------- Load Input Data ---------
metadata_cross2 <- read.csv("data/metadata_cross_shann_dirichlet.csv", sep = ";")

metadata_cross2_sensi <- metadata_cross2[-13,]
metadata_cross2_sensi$Endotype_label <- as.factor(metadata_cross2_sensi$Endotype_label)
metadata_cross2_sensi$gender <- as.factor(metadata_cross2_sensi$gender)
# ====================================================
# Model 1: Negative Binomial Adjusted for Gender
# ====================================================
nb_model2 <- glm.nb(nwheezing ~ Endotype_label * gender, data = metadata_cross2_sensi)

summary(nb_model2)

# Check dispersion
dispersion_nb2 <- sum(residuals(nb_model2, type = "pearson")^2) / nb_model2$df.residual
cat("Dispersion ratio:", dispersion_nb2, "\n")  # Should be ~1

exp(coef(nb_model2))              # IRRs
exp(confint(nb_model2))          # 95% CIs

# ====================================================
# Generate Predicted vs Observed Plot
# ====================================================
# Extract model frame (actual data used)
model_data <- model.frame(nb_model2)

# Add predicted values
model_data$predicted <- predict(nb_model2, type = "response")

# Group-level means for plotting
summary_data <- model_data %>%
  group_by(Endotype_label, gender) %>%
  summarise(
    mean_obs = mean(nwheezing, na.rm = TRUE),
    se_obs = sd(nwheezing, na.rm = TRUE) / sqrt(n()),
    mean_pred = mean(predicted, na.rm = TRUE),
    .groups = "drop"
  )

# Plot: Observed vs Predicted Wheeze Episodes
pdf("figures/FigS11.pdf", 8, 8)
ggplot(model_data, aes(x = interaction(Endotype_label, gender), y = nwheezing)) +
  geom_jitter(
    position = position_jitter(width = 0.2, height = 0),
    alpha = 0.8,
    size = 3,
    color = "blue"
  ) +
  stat_summary(
    fun = mean,
    geom = "bar",
    fill = "#69b3a2",
    alpha = 0.6,
    width = 0.6
  ) +
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.2,
    color = "black"
  ) +
  geom_point(aes(y = predicted), shape = 18, size = 6, color = "red") +
  labs(
    x = "Cluster by Gender",
    y = "Number of Wheezing Episodes",
    title = "Observed vs Predicted Wheezing Episodes by Cluster and Gender"
  ) +
  theme_minimal()
dev.off()

# ====================================================
# Marginal Means and Contrast Test (Dunnett)
# ====================================================
emm <- emmeans(nb_model2, ~ Endotype_label | gender)

# Pairwise comparisons using Dunnett adjustment
contrast_results <- contrast(emm, method = "trt.vs.ctrl", adjust = "dunnett")
print(summary(contrast_results))