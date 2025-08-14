
# ====================================
# AERIAL microbiome manuscript – Figure S2
# ====================================

# Description: Generates Figure S2
# Author: Jose Caparros Martin, PhD
# Email: joseantonio.caparrosmartin@uwa.edu.au
# License: MIT
# Disclaimer: Provided "as is" for academic/research use. 
# Date: 2025-08-05


# Load required libraries
library(ggplot2)
library(RColorBrewer)
library(tidyverse)
library(rstatix)

plot_qpcr <- function(data_path = "data/qpcr_data.csv",
                      output_dir = "figures",
                      output_file = "paper_qpcr_FigS2.pdf") {
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  qpcr <- read.csv(data_path, sep = ";")
  
  pdfname <- file.path(output_dir, output_file)
  
  # Example: calculate y-intercept as mean + 3*SD of negative controls
  neg_controls <- qpcr %>% filter(Sampletype == "Negative")
  y_intercept <- mean(neg_controls$X1ct, na.rm = TRUE) + 
    3 * sd(neg_controls$X1ct, na.rm = TRUE)
  
  p <- ggplot(qpcr, aes(x = Sampletype, y = X1ct)) +
    geom_boxplot(notch = FALSE, width = 0.5, lwd = 1, outlier.shape = NA) +
    geom_violin(fill = 'blue', color = NA, alpha = 0.2) +
    geom_jitter(shape = 16, position = position_jitter(0.2), color = "red") +
    scale_color_brewer(palette = "Set1") +
    geom_hline(yintercept = y_intercept, linetype = "solid", color = "red", size = 1) +
    labs(
      title = "Bacterial DNA quantification",
      x = "Sample Type",
      y = "Bacterial load 1/Ct"
    ) +
    theme_classic()
  
  pdf(pdfname, width = 4, height = 4)
  print(p)
  dev.off()
  
  message("Plot saved to: ", pdfname)
  
  return(qpcr)
}

# Run the function and store the data for further stats
qpcr_data <- plot_qpcr(
  data_path = "data/qpcr_data.csv",
  output_dir = "figures",
  output_file = "paper_qpcr_FigS2.pdf"
)

# Statistical analysis
qpcr_data %>% wilcox_test(X1ct ~ Sampletype, p.adjust.method = "bonferroni")