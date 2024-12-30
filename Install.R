.onLoad <- function(libname, pkgname) {
  # List of CRAN packages that need to be checked and installed
  cran_packages <- c(
    "survival", "survminer", "kableExtra", "ggplot2", "data.table", "tidyverse",
    "parallel", "pbapply", "patchwork", "doParallel", "dplyr", "gridExtra"
  )

  # List of Bioconductor packages
  bioc_packages <- c(
    "TCGAimmunosurv", "TCGAbiolinks", "SummarizedExperiment", "DESeq2",
    "clusterProfiler", "msigdbr", "SingleR", "celldex", "Seurat", "SeuratObject",
    "BiocParallel", "SeuratWrappers"
  )

  # Install CRAN packages if not installed
  for (pkg in cran_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message("Installing CRAN package: ", pkg)
      install.packages(pkg)
    }
  }

  # Install Bioconductor packages if not installed
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    message("Installing BiocManager")
    install.packages("BiocManager")
  }

  for (pkg in bioc_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message("Installing Bioconductor package: ", pkg)
      BiocManager::install(pkg)
    }
  }

  # Install moncole3 from GitHub if not installed
  if (!requireNamespace("monocle3", quietly = TRUE)) {
    message("Installing monocle3 from GitHub")
    if (!requireNamespace("devtools", quietly = TRUE)) {
      message("Installing devtools package")
      install.packages("devtools")
    }
    devtools::install_github("monocle3/monocle3")
  }

  # Optional: You can include any other installation-related instructions here
}
