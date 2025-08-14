#' Filter and Normalize OTU or ASV Table
#'
#' Filters OTUs based on a minimum relative abundance percentage
#' and optionally normalizes samples using total sum scaling (TSS).
#'
#' @param data OTU count data frame (samples x OTUs)
#' @param percent Cutoff for filtering OTUs (default = 0.01%)
#' @param filter Logical, whether to filter low-abundance OTUs
#' @param normalize Logical, whether to normalize by total sum scaling
#'
#' @return A list with filtered data, optional normalized data, and OTU indices
filter_and_normalize <- function(data, percent = 0.01, filter = TRUE, normalize = TRUE) {
  data.filtered <- data
  keep.otu <- NULL
  
  # Step 1: Optional filtering
  if (filter) {
    keep.otu <- which(colSums(data) * 100 / sum(colSums(data)) > percent)
    data.filtered <- data[, keep.otu, drop = FALSE]
  }
  
  # Step 2: Optional normalization
  if (normalize) {
    data.normalized <- t(apply(data.filtered, 1, function(x) x / sum(x)))
    return(list(
      data.filtered = data.filtered,
      data.normalized = data.normalized,
      keep.otu = keep.otu
    ))
  } else {
    return(list(
      data.filtered = data.filtered,
      keep.otu = keep.otu
    ))
  }
}
