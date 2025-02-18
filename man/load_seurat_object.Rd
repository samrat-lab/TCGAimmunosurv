% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TCGAimmunosurv.R
\name{load_seurat_object}
\alias{load_seurat_object}
\title{Load Seurat Object for Single-Cell Analysis}
\usage{
load_seurat_object(
  directory,
  constant_prefix = "combined_seurat_",
  selected_cancer
)
}
\arguments{
\item{directory}{A character string specifying the directory where the Seurat object file is located.}

\item{constant_prefix}{A character string specifying the constant prefix for the Seurat object file name. The default is "combined_seurat_".}

\item{selected_cancer}{A character string specifying the cancer type for which the Seurat object file will be loaded. The cancer type will be appended to the constant prefix.}
}
\value{
A Seurat object that has been loaded from the specified file.
}
\description{
This function loads a Seurat object from a specified directory based on a given cancer type.
The Seurat object will be used to perform further single-cell RNA-seq analysis.
}
\details{
This function constructs the file name for the Seurat object based on the provided \code{constant_prefix} and \code{selected_cancer}. It checks if the file exists in the specified directory and loads the Seurat object if the file is found. If the file does not exist, an error message is displayed.
}
