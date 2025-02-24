% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TCGAimmunosurv.R
\name{plot_unique_mutation_types}
\alias{plot_unique_mutation_types}
\title{Plot Unique Mutation Types Across Genes}
\usage{
plot_unique_mutation_types(cancer_maf)
}
\arguments{
\item{cancer_maf}{A \code{MAF} object representing mutation data for a specific cancer type.
This object is typically created using the \code{maftools::read.maf} function.}
}
\value{
A ggplot object showing the bar plot of unique mutation types and their total number
of samples across all genes in the provided \code{MAF} object.
}
\description{
This function generates a bar plot to visualize the distribution of unique mutation types
across genes for a given cancer type based on the MAF (Mutation Annotation Format) data.
}
\details{
\itemize{
\item The function extracts the mutation type counts from the \code{MAF} object using the
\code{maftools::getGeneSummary} function.
\item Mutation types with zero counts are filtered out.
\item The mutation counts are aggregated across all genes for each mutation type.
\item The resulting plot displays mutation types on the x-axis and the total number of
samples on the y-axis.
}
}
