
# ======================================
# AERIAL Project – Figure S8
# ======================================
# Author: Jose Caparros Martin, PhD
# Email: joseantonio.caparrosmartin@uwa.edu.au
# License: MIT
# Disclaimer: Provided "as is" for academic/research use.
# Date: 2025-08-05

# --------- Libraries ---------
library(ggplot2)
library(dplyr)
library(rstatix)
library(viridis)

# --------- Load Input Data ---------
metadata_cross2 <- read.csv("data/metadata_cross_shann_dirichlet.csv", sep = ";")

# --------- Helper: Plot stacked barchart ---------
plot_stacked_barchart <- function(df, variable, cluster_col = "Endotype_label", out_pdf, exclude_rows = NULL) {
  # Ensure input is a data.frame
  if (!is.data.frame(df)) {
    stop("`df` must be a data.frame.")
  }
  
  # Check if required columns exist
  if (!all(c(variable, cluster_col) %in% names(df))) {
    stop(paste("Missing column(s):", paste(setdiff(c(variable, cluster_col), names(df)), collapse = ", ")))
  }
  
  # Optional row exclusion
  if (!is.null(exclude_rows)) {
    df <- df[-exclude_rows, , drop = FALSE]
  }
  
  # Remove rows with NA
  df <- df %>% filter(!is.na(.data[[cluster_col]]) & !is.na(.data[[variable]]))
  
  # Convert the variable column to a factor to avoid viridis error
  df[[variable]] <- as.factor(df[[variable]])
  
  # Exit if data is empty
  if (nrow(df) == 0) {
    warning("No data left after filtering. Exiting function.")
    return(NULL)
  }
  
  # Counts
  cnt <- df %>%
    group_by(!!sym(cluster_col), !!sym(variable)) %>%
    summarise(n = n(), .groups = "drop")
  
  if (nrow(cnt) == 0) {
    warning("No data left to plot after grouping. Exiting function.")
    return(NULL)
  }
  
  # Proportions
  pcnt <- cnt %>%
    group_by(!!sym(cluster_col)) %>%
    mutate(pcnt = n / sum(n)) %>%
    ungroup()
  
  # Label positions
  ref_levels <- levels(df[[variable]])
  label_cond <- if (length(ref_levels) > 0) {
    pcnt[[variable]] == ref_levels[1]
  } else {
    rep(TRUE, nrow(pcnt))
  }
  
  KS <- pcnt %>%
    mutate(labelpos = ifelse(label_cond, pcnt / 2, 1 - pcnt / 2))
  
  # Plot
  pdf(out_pdf, 5, 5)
  gg <- ggplot(KS, aes(x = .data[[cluster_col]], y = pcnt, fill = .data[[variable]])) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_viridis(discrete = TRUE) +
    theme_classic()
  print(gg)
  dev.off()
  
  # Stats
  cat("\n--- Chi-square test for", variable, "---\n")
  tab <- table(df[[cluster_col]], df[[variable]])
  if (all(dim(tab) > 1)) {
    print(chisq.test(tab))
  } else {
    cat("Not enough levels to perform chi-square test.\n")
  }
}


# --------- Fig S8A ---------
# Run for gender
plot_stacked_barchart(
  df = metadata_cross2,
  variable = "gender",
  out_pdf = "figures/gender_S8A.pdf"
)


# --------- Fig S8B ---------
# Run for season
plot_stacked_barchart(
  df = metadata_cross2,
  variable = "Season2",
  out_pdf = "figures/gender_S8B.pdf"
)


# --------- Fig S8C ---------
# Run for gender
plot_stacked_barchart(
  df = metadata_cross2,
  variable = "ever_childcare",
  out_pdf = "figures/gender_S8C.pdf"
)