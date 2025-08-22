# Microbiome AERIAL 🧬 

📌 This repo contains the `R` code used to generate the microbiome analysis figures associated with the manuscript titled  *"Insights into the relationship between nasal bacterial composition and susceptibility to early-life respiratory disease: a pilot observational study"*, available through [MedRxiv](https://www.medrxiv.org/content/10.1101/2025.08.16.25333459v1).

## Structure
- **data/** — input datasets used for analysis (*not yet available*)  
- **RScripts/** — R scripts for processing, modeling, and plotting  
- **GITHUB.Rproj** — RStudio project file  

## Reproducibility
We use [`renv`](https://rstudio.github.io/renv/) to manage package versions.

### To set up:
```r
install.packages("renv")
renv::restore()
