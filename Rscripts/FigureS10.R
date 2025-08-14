# ====================================
# AERIAL Microbiome Manuscript – Figure S10
# ====================================
# Description: LDA model FigS10
# Author: Jose Caparros Martin, PhD
# Email: joseantonio.caparrosmartin@uwa.edu.au
# License: MIT
# Date: 2025-08-05

# -------------------------
# Load Required Libraries
# -------------------------
library(tidyverse)
library(tidytext)      # tidy() for topicmodels objects (matrix="beta"/"gamma")
library(topicmodels)
library(tm)
library(MASS)
library(ggeffects)
library(broom)
library(AER)
library(dplyr)

set.seed(1492)

# -------------------------
# Create Output Folder
# -------------------------
folder_path <- "figures"
if (!dir.exists(folder_path)) dir.create(folder_path, recursive = TRUE)


# -------------------------
# Load Input Data
# -------------------------
metadata_cross2 <- read.csv("data/metadata_cross_shann_dirichlet.csv", sep = ";")
data_final      <- read.csv("data/16S_data_final.csv", sep = ";")
taxonomy        <- read.csv("data/16S_taxonomy.csv", sep = ";")


# -------------------------
# Fit LDA Model (k = 2)
# -------------------------
data_micro <- data_final[,1:71]
rownames(data_micro) <- taxonomy$Species
x <- t(data_micro)
otu_table_matrix <- as.matrix(x)
otu_dtm <- as.DocumentTermMatrix(otu_table_matrix, weighting = weightTf)

# Fit the model k=2
lda_model <- LDA(otu_dtm, k = 2, control = list(seed = 1492))


# -------------------------
# Helper: Heatmap of top genera across topics
# -------------------------
generate_heatmap_lda <- function(lda_model, output_path) {
  # make sure needed namespaces are available
  stopifnot(requireNamespace("tidytext"),
            requireNamespace("dplyr"),
            requireNamespace("ggplot2"),
            requireNamespace("rlang"),
            requireNamespace("grid"))
  
  # 1) Get beta (per-topic term probabilities)
  beta_df <- tidytext::tidy(lda_model, matrix = "beta")
  
  # 2) Find the term/feature column name
  term_col <- intersect(names(beta_df), c("term", "token", "feature", "word"))[1]
  if (is.null(term_col) || is.na(term_col)) {
    stop("Could not find a term/feature column in `beta_df`. Columns are: ",
         paste(names(beta_df), collapse = ", "))
  }
  # Basic sanity: we also need 'topic' and 'beta'
  if (!all(c("topic", "beta") %in% names(beta_df))) {
    stop("`beta_df` must contain columns 'topic' and 'beta'. Columns are: ",
         paste(names(beta_df), collapse = ", "))
  }
  
  # 3) Standardize column names for downstream code
  beta_df <- beta_df |>
    dplyr::rename(.term  = !!rlang::sym(term_col),
                  .topic = topic,
                  .beta  = beta)
  
  # 4) Top 3 terms per topic
  top_genera <- beta_df |>
    dplyr::group_by(.topic) |>
    dplyr::slice_max(.beta, n = 3, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::mutate(home_topic = .topic)
  
  # 5) All topic probabilities for those top terms
  top_genus_across_topics <- beta_df |>
    dplyr::filter(.term %in% unique(top_genera$.term)) |>
    dplyr::left_join(dplyr::select(top_genera, .term, home_topic),
                     by = ".term") |>
    dplyr::mutate(term_grouped = paste0("T", home_topic, "_", .term),
                  label = sprintf("%.2f", .beta))
  
  # 6) Plot
  p <- ggplot2::ggplot(top_genus_across_topics,
                       ggplot2::aes(x = term_grouped,
                                    y = factor(.topic),
                                    fill = .beta)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = label), color = "black", size = 3.5) +
    ggplot2::scale_fill_gradientn(
      colours = c("#542788", "#998ec3", "#f7f7f7", "#f1a340", "#b35806"),
      name = "Probability",
      limits = c(0, max(top_genus_across_topics$.beta, na.rm = TRUE)),
      guide = ggplot2::guide_colorbar(
        direction = "horizontal",
        barwidth  = grid::unit(8, "cm"),
        barheight = grid::unit(0.5, "cm"),
        title.position = "top",
        title.hjust = 0.5
      )
    ) +
    ggplot2::labs(
      title = "Genus Distributions Across Microbial Assemblages",
      x = "Genus (Grouped by Home Topic)",
      y = "Assemblage"
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.title = ggplot2::element_text(size = 12),
      legend.text = ggplot2::element_text(size = 10),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid = ggplot2::element_blank()
    )
  
  ggplot2::ggsave(output_path, plot = p, width = 10, height = 8)
}

# Generate Heatmap
generate_heatmap_lda(lda_model, "figures/FigS14A.pdf")


# -------------------------
# Figure S10B
# -------------------------
topic_proportions <- posterior(lda_model)$topics

df <- as.data.frame(topic_proportions)

# Convert rownames to a column
df <- df %>% rownames_to_column("Sample")

colnames(df) <- c("Sample", "Topic1", "Topic2")

# Arrange by Topic 1 and factor the Sample to preserve order in the plot
df_ordered <- df %>%
  arrange(`Topic1`) %>%
  mutate(Sample = factor(Sample, levels = Sample)) %>%
  pivot_longer(cols = c(`Topic1`, `Topic2`), names_to = "Topic", values_to = "Probability")

# Plot
pdf("figures/FigS10B.pdf", 6,6)

ggplot(df_ordered, aes(x = Sample, y = Probability, color = Topic)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  theme_minimal() +
  coord_flip() +
  labs(title = "Topic Probabilities per Sample (Ordered by Topic 1)", x = "Sample", y = "Probability")
dev.off()


# -------------------------
# Figure S10C - Marginal effects
# -------------------------
df_work <- df %>%
  dplyr::mutate(
    Condition      = metadata_cross2$Typeofsample,
    numberWheezing = metadata_cross2$nwheezing,
    gender_raw     = metadata_cross2$gender,
    Cluster        = ifelse(Topic1 > Topic2, "Cluster1", "Cluster2")
  ) %>%
  dplyr::select(Topic1, Topic2, Condition, numberWheezing, gender_raw, Cluster) %>%
  dplyr::filter(stats::complete.cases(Topic1, numberWheezing, gender_raw, Condition, Cluster)) %>%
  dplyr::mutate(
    Condition = factor(Condition, levels = c("Baseline", "Symptomatic")),
    Cluster   = factor(Cluster,   levels = c("Cluster1", "Cluster2")),
    gender    = factor(gender_raw, levels = c("F", "M"))
  )


# ---- Poisson model (for dispersion check) ----
pois_mod <- glm(numberWheezing ~ Topic1 * gender,
                family = poisson(link = "log"),
                data   = df_work)

# Overdispersion test
disp <- AER::dispersiontest(pois_mod)
print(disp)
# If p < 0.05 → overdispersed → use Negative Binomial 

# ---- Negative Binomial model ----
nb_mod <- MASS::glm.nb(numberWheezing ~ Topic1 * gender, data = df_work)
summary(nb_mod)

# ---- Marginal effects ----
me <- ggeffects::ggpredict(nb_mod, terms = c("Topic1", "gender"))

# ---- Plot ----
outfile <- "figures/FigS10C.pdf"
pdf(outfile, width = 6, height = 6)

plot(me) +
  ggtitle("Predicted incidence rates by Topic1 and gender") +
  xlab("Topic1 probability") +
  ylab("Predicted incidence rate") +
  theme_minimal(base_size = 12)

dev.off()

message("Saved plot to: ", normalizePath(outfile))


