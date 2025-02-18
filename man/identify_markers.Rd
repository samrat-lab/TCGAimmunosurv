% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TCGAimmunosurv.R
\name{identify_markers}
\alias{identify_markers}
\title{Identify Cluster Markers}
\usage{
identify_markers(seurat_obj, min_pct = 0.25, logfc_thresh = 0.25)
}
\arguments{
\item{seurat_obj}{A Seurat object with clustered cells.}

\item{min_pct}{Minimum percentage of cells expressing the gene to be considered as a marker (default is 0.25).}

\item{logfc_thresh}{Minimum log fold change for a gene to be considered significant (default is 0.25).}
}
\value{
A data frame containing identified markers for each cluster.
}
\description{
This function identifies marker genes for each cluster.
}
