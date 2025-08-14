# Microbiome AERIAL 🧬 

📌 This repo contains the `R` code used to generate the microbiome analysis figures associated with the manuscript titled  *"Insights into the relationship between nasal bacterial composition and susceptibility to early-life respiratory disease: a pilot observational study"*, available through [MedRxiv](https://www.medrxiv.org/) and [Pubmed](https://pubmed.ncbi.nlm.nih.gov/) **_To be updated_**.


## Structure
- **data/** — input datasets used for analysis  
- **RScripts/** — R scripts for processing, modeling, and plotting  
- **AERIAL.Rproj** — RStudio project file  

## Reproducibility
We use [`renv`](https://rstudio.github.io/renv/) to manage package versions.

### To set up:
```r
install.packages("renv")
renv::restore()
