# ====================================
# AERIAL Microbiome Manuscript – Figure S13
# ====================================
# Description: Marginal predicted values FigS13
# Author: Jose Caparros Martin, PhD
# Email: joseantonio.caparrosmartin@uwa.edu.au
# License: MIT
# Date: 2025-08-05

# -------------------------
# Load Required Libraries
# -------------------------
library(broom)
library(ggeffects)
library(dplyr)
library(ggplot2)

# --------- Load Input Data ---------
metadata_cross2 <- read.csv("data/metadata_cross_shann_dirichlet.csv", sep = ";")

metadata_cross2 <- metadata_cross2 %>%
  mutate(
    Endotype_label = factor(Endotype_label, levels = c("cluster 1","cluster 2")),
    gender = factor(gender, levels = c("F","M")),
    ever_wheeze = as.factor(ever_wheeze),
    ever_recurrentwheeze = as.factor(ever_recurrentwheeze),
    ever_bronchiolitis = as.factor(ever_bronchiolitis)
  )

# ==================================================================
# Model 1: Binomial Regression (Unadjusted) - ever wheeze
# ==================================================================

# Fit logistic regression
binomial_model1 <- glm(
  ever_wheeze ~ Endotype_label,
  family = binomial,
  data = metadata_cross2
)

# Tidy model output and exponentiate
model_tidy <- broom::tidy(
  binomial_model1,
  conf.int = TRUE,
  conf.level = 0.95,
  exponentiate = TRUE
)

# Now rename safely if columns exist
if(all(c("term","estimate","conf.low","conf.high","p.value") %in% names(model_tidy))) {
  results_table <- model_tidy %>%
    dplyr::filter(.data$term != "(Intercept)") %>%     # note dplyr:: and .data$
    dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
    dplyr::rename(
      Term = term,
      OR = estimate,
      `CI Lower` = conf.low,
      `CI Upper` = conf.high,
      `p-value` = p.value
    )
  
  print(results_table)
} else {
  stop("Tidy output does not contain expected columns. Check `model_tidy` above.")
}

# Marginal Effects Plot
marginal_binomial_model1 <- ggpredict(binomial_model1, terms = "Endotype_label")

pdf("figures/S13A.pdf", width = 8, height = 8)
plot(marginal_binomial_model1) +
  ggtitle("Predicted probability of ever wheeze by cluster") +
  ylab("Predicted probability") +
  xlab("Cluster") +
  theme_classic()
dev.off()

# ==================================================================
# Model 2: Binomial Regression (Unadjusted) - ever recurrent wheeze
# ==================================================================

# Fit logistic regression
binomial_model2 <- glm(
  ever_recurrentwheeze ~ Endotype_label,
  family = binomial,
  data = metadata_cross2
)

# Tidy model output and exponentiate
model_tidy <- broom::tidy(
  binomial_model2,
  conf.int = TRUE,
  conf.level = 0.95,
  exponentiate = TRUE
)

# Now rename safely if columns exist
if(all(c("term","estimate","conf.low","conf.high","p.value") %in% names(model_tidy))) {
  results_table <- model_tidy %>%
    dplyr::filter(.data$term != "(Intercept)") %>%     # note dplyr:: and .data$
    dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
    dplyr::rename(
      Term = term,
      OR = estimate,
      `CI Lower` = conf.low,
      `CI Upper` = conf.high,
      `p-value` = p.value
    )
  
  print(results_table)
} else {
  stop("Tidy output does not contain expected columns. Check `model_tidy` above.")
}

# Marginal Effects Plot
marginal_binomial_model2 <- ggpredict(binomial_model2, terms = "Endotype_label")

pdf("figures/S13B.pdf", width = 8, height = 8)
plot(marginal_binomial_model2) +
  ggtitle("Predicted probability of ever recurrent wheeze by cluster") +
  ylab("Predicted probability") +
  xlab("Cluster") +
  theme_classic()
dev.off()


# ==================================================================
# Model 3: Binomial Regression (Unadjusted) - ever bronchiolitis
# ==================================================================

# Fit logistic regression
binomial_model3 <- glm(
  ever_bronchiolitis ~ Endotype_label,
  family = binomial,
  data = metadata_cross2
)

# Tidy model output and exponentiate
model_tidy <- broom::tidy(
  binomial_model3,
  conf.int = TRUE,
  conf.level = 0.95,
  exponentiate = TRUE
)

# Now rename safely if columns exist
if(all(c("term","estimate","conf.low","conf.high","p.value") %in% names(model_tidy))) {
  results_table <- model_tidy %>%
    dplyr::filter(.data$term != "(Intercept)") %>%     # note dplyr:: and .data$
    dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
    dplyr::rename(
      Term = term,
      OR = estimate,
      `CI Lower` = conf.low,
      `CI Upper` = conf.high,
      `p-value` = p.value
    )
  
  print(results_table)
} else {
  stop("Tidy output does not contain expected columns. Check `model_tidy` above.")
}

# Marginal Effects Plot
marginal_binomial_model3 <- ggpredict(binomial_model3, terms = "Endotype_label")

pdf("figures/S13C.pdf", width = 8, height = 8)
plot(marginal_binomial_model3) +
  ggtitle("Predicted probability of ever bronchiolitis by cluster") +
  ylab("Predicted probability") +
  xlab("Cluster") +
  theme_classic()
dev.off()



# ==================================================================
# Model 4: Binomial Regression (sex adjusted) - ever wheeze
# ==================================================================

# Fit logistic regression with interaction
binomial_model4 <- glm(
  ever_wheeze ~ Endotype_label * gender,
  family = binomial,
  data = metadata_cross2
)

# Tidy ORs with 95% CIs
model_tidy <- broom::tidy(binomial_model4, conf.int = TRUE, conf.level = 0.95, exponentiate = TRUE)

# Safe rename + drop intercept
results_table <- model_tidy %>%
  dplyr::filter(.data$term != "(Intercept)") %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  dplyr::rename(
    Term = term,
    OR = estimate,
    `CI Lower` = conf.low,
    `CI Upper` = conf.high,
    `p-value` = p.value
  )

print(results_table)

# Marginal effects for interaction
marginal_binomial_model4 <- ggpredict(binomial_model4, terms = c("Endotype_label","gender"))

pdf("figures/S13D.pdf", width = 6, height = 5)
plot(marginal_binomial_model4) +
  ggtitle("Predicted probability of ever wheeze by cluster and gender") +
  ylab("Predicted probability") +
  xlab("Cluster") +
  theme_classic()
dev.off()

# ==================================================================
# Model 5: Binomial Regression (sex adjusted) - ever recurrent wheeze
# ==================================================================

# Fit logistic regression with interaction
binomial_model5 <- glm(
  ever_recurrentwheeze ~ Endotype_label * gender,
  family = binomial,
  data = metadata_cross2
)

# Tidy ORs with 95% CIs
model_tidy <- broom::tidy(binomial_model5, conf.int = TRUE, conf.level = 0.95, exponentiate = TRUE)

# Safe rename + drop intercept
results_table <- model_tidy %>%
  dplyr::filter(.data$term != "(Intercept)") %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  dplyr::rename(
    Term = term,
    OR = estimate,
    `CI Lower` = conf.low,
    `CI Upper` = conf.high,
    `p-value` = p.value
  )

print(results_table)

# Marginal effects for interaction
marginal_binomial_model5 <- ggpredict(binomial_model5, terms = c("Endotype_label","gender"))

pdf("figures/S13E.pdf", width = 6, height = 5)
plot(marginal_binomial_model5) +
  ggtitle("Predicted probability of ever recurrent wheeze by cluster and gender") +
  ylab("Predicted probability") +
  xlab("Cluster") +
  theme_classic()
dev.off()


# ==================================================================
# Model 6: Binomial Regression (sex adjusted) - ever bronchiolitis
# ==================================================================

# Fit logistic regression with interaction
binomial_model6 <- glm(
  ever_bronchiolitis ~ Endotype_label * gender,
  family = binomial,
  data = metadata_cross2
)

# Tidy ORs with 95% CIs
model_tidy <- broom::tidy(binomial_model6, conf.int = TRUE, conf.level = 0.95, exponentiate = TRUE)

# Safe rename + drop intercept
results_table <- model_tidy %>%
  dplyr::filter(.data$term != "(Intercept)") %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  dplyr::rename(
    Term = term,
    OR = estimate,
    `CI Lower` = conf.low,
    `CI Upper` = conf.high,
    `p-value` = p.value
  )

print(results_table)

# Marginal effects for interaction
marginal_binomial_model6 <- ggpredict(binomial_model6, terms = c("Endotype_label","gender"))

pdf("figures/S13F.pdf", width = 6, height = 5)
plot(marginal_binomial_model6) +
  ggtitle("Predicted probability of ever bronchiolitis by cluster and gender") +
  ylab("Predicted probability") +
  xlab("Cluster") +
  theme_classic()
dev.off()



