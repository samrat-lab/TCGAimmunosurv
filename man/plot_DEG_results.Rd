% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TCGAimmunosurv.R
\name{plot_DEG_results}
\alias{plot_DEG_results}
\title{Plot Differentially Expressed Genes (DEGs) Results}
\usage{
plot_DEG_results(
  res_df,
  rna_data,
  condition_column,
  log2FC_threshold = 1,
  alpha = 0.05
)
}
\arguments{
\item{res_df}{A data frame containing the results of the differential expression analysis,
including columns for \code{log2FoldChange}, \code{padj} (adjusted p-values), and \code{significance}.
The \code{significance} column should indicate whether a gene is "Upregulated",
"Downregulated", or "Not Significant".}

\item{rna_data}{A \code{SummarizedExperiment} object containing RNA-Seq data with
sample metadata accessible via \code{SummarizedExperiment::colData}.}

\item{condition_column}{A character string specifying the name of the column in
\code{SummarizedExperiment::colData(rna_data)} that contains condition labels (e.g., "Control" and "Treatment").}

\item{log2FC_threshold}{A numeric value specifying the threshold for \code{log2FoldChange}.
Default is 1.}

\item{alpha}{A numeric value specifying the significance threshold for adjusted p-values.
Default is 0.05.}
}
\value{
None. This function outputs two plots:
\itemize{
\item A bar plot showing the distribution of samples across conditions.
\item A volcano plot visualizing the DEGs.
}
}
\description{
This function generates two plots to visualize the results of differential expression analysis:
a bar plot showing the sample distribution for different conditions and a volcano plot
for visualizing the differentially expressed genes (DEGs).
}
\details{
\itemize{
\item The bar plot displays the sample count for each condition in the dataset, with labels
for the number of samples above each bar.
\item The volcano plot uses \code{log2FoldChange} and \code{-log10(padj)} to visualize the significance
and magnitude of gene expression changes. Points are colored based on their significance
status.
\item Both plots use customizable font styles for better readability.
}
}
