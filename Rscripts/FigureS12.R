# ====================================
# AERIAL Microbiome Manuscript – Figure S12
# ====================================
# Description: Marginal predicted values FigS12
# Author: Jose Caparros Martin, PhD
# Email: joseantonio.caparrosmartin@uwa.edu.au
# License: MIT
# Date: 2025-08-05

# -------------------------
# Load Required Libraries
# -------------------------
library(MASS)
library(ggeffects)
library(AER)
library(ggplot2)
library(dplyr)
library(emmeans)

# --------- Load Input Data ---------
metadata_cross2 <- read.csv("data/metadata_cross_shann_dirichlet.csv", sep = ";")

metadata_cross2$Endotype_label <- as.factor(metadata_cross2$Endotype_label)
metadata_cross2$gender <- as.factor(metadata_cross2$gender)
# ==================================================================
# Model 1: Poisson Regression (Unadjusted) symptomatic swabs
# ==================================================================
poisson_model1 <- glm(
  Total_symptomatic ~ Endotype_label,
  family = poisson(link = "log"),
  data = metadata_cross2
)

summary(poisson_model1)
exp(coef(poisson_model1))            # Incidence Rate Ratios (IRRs)
exp(confint(poisson_model1))        # 95% Confidence Intervals

# Check for overdispersion
disp_test <- AER::dispersiontest(poisson_model1) 
print(disp_test) 

# Marginal Effects Plot
marginal_poisson_model1 <- ggpredict(poisson_model1, terms = c("Endotype_label"))

pdf("figures/S12A.pdf", width = 8, height = 8)
plot(marginal_poisson_model1) +
  ggtitle("Marginal Effects: Total_symptomatic by Cluster") +
  ylab("Predicted Number of Symptomatic swabs") +
  xlab("Cluster")
dev.off()

# ===================================================================
# Model 2: Poisson Regression Adjusted for Gender (symptomatic swabs)
# ===================================================================
poisson_model2 <- glm(
  Total_symptomatic ~ Endotype_label * gender,
  family = poisson(link = "log"),
  data = metadata_cross2
)

summary(poisson_model2)
exp(coef(poisson_model2))            # Incidence Rate Ratios (IRRs)
exp(confint(poisson_model2))        # 95% Confidence Intervals

# Check for overdispersion
disp_test <- AER::dispersiontest(poisson_model2) 
print(disp_test) 

# Marginal Effects Plot
marginal_poisson_model2 <- ggpredict(poisson_model2, terms = c("Endotype_label", "gender"))

pdf("figures/S12C.pdf", width = 8, height = 8)
plot(marginal_poisson_model2) +
  ggtitle("Marginal Effects: Total_symptomatic by Cluster and gender") +
  ylab("Predicted Number of Symptomatic swabs adjusted by gender") +
  xlab("Cluster")
dev.off()



# ==================================================================
# Model 3: Poisson Regression (Unadjusted) Bronchilitis episodes
# ==================================================================


poisson_model3 <- glm(
  nbronchiolitis ~ Endotype_label,
  family = poisson(link = "log"),
  data = metadata_cross2
)

summary(poisson_model3)
exp(coef(poisson_model3))            # Incidence Rate Ratios (IRRs)
exp(confint(poisson_model3))        # 95% Confidence Intervals

# Check for overdispersion
disp_test <- AER::dispersiontest(poisson_model3) 
print(disp_test) 

# Marginal Effects Plot
marginal_poisson_model3 <- ggpredict(poisson_model3, terms = c("Endotype_label"))

pdf("figures/S12B.pdf", width = 8, height = 8)
plot(marginal_poisson_model3) +
  ggtitle("Marginal Effects: Bronchiolitis by Cluster") +
  ylab("Predicted Number of Bronchiolitis episodes") +
  xlab("Cluster")
dev.off()


# ======================================================================
# Model 4: Poisson Regression Adjusted for Gender (Bronchilitis episodes)
# ======================================================================
poisson_model4 <- glm(
  nbronchiolitis ~ Endotype_label * gender,
  family = poisson(link = "log"),
  data = metadata_cross2
)

summary(poisson_model4)
exp(coef(poisson_model4))            # Incidence Rate Ratios (IRRs)
exp(confint(poisson_model4))        # 95% Confidence Intervals

# Check for overdispersion
disp_test <- AER::dispersiontest(poisson_model4) 
print(disp_test) 

# Marginal Effects Plot
marginal_poisson_model4 <- ggpredict(poisson_model4, terms = c("Endotype_label", "gender"))

pdf("figures/S12D.pdf", width = 8, height = 8)
plot(marginal_poisson_model4) +
  ggtitle("Marginal Effects: Bronchiolitis by Cluster and gender") +
  ylab("Predicted Number of Bronchiolitis episodes adjusted by gender") +
  xlab("Cluster")
dev.off()