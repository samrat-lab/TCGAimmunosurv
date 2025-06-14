% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TCGAimmunosurv.R
\name{perform_survival_analysis_with_permutation}
\alias{perform_survival_analysis_with_permutation}
\title{Survival Analysis with Permutation Test for Robust Gene Selection}
\usage{
perform_survival_analysis_with_permutation(
  count_matrix_filtered,
  gene_metadata_dt,
  clinical_data_filtered,
  n_permutations = 1000
)
}
\arguments{
\item{count_matrix_filtered}{A numeric matrix of gene expression (genes x samples) with gene IDs as rows and sample IDs as columns.}

\item{gene_metadata_dt}{A data.table containing at least a column named \code{gene_name} corresponding to the genes in \code{count_matrix_filtered}.}

\item{clinical_data_filtered}{A data.table with clinical data, including columns \code{submitter_id}, \code{overall_survival}, and \code{deceased}.}

\item{n_permutations}{Number of permutations to perform for the permutation test (default: 1000).

Genes are reported as significant if their empirical permutation p-value (\code{perm_pval}) is less than or equal to 0.05.}
}
\value{
A data.frame with one row per gene, containing:
\describe{
\item{gene_name}{Gene name}
\item{hazard_ratio}{Hazard ratio (HIGH vs LOW expression)}
\item{pvalue}{Wald p-value from the Cox model}
\item{logrank_stat}{Observed log-rank (likelihood ratio) statistic}
\item{perm_pval}{Empirical p-value from permutation test}
\item{significant}{"Yes" if perm_pval <= 0.05, otherwise "No"}
}
}
\description{
This function performs survival analysis for each gene using a Cox proportional hazards model,
stratifies expression as HIGH/LOW per gene, and assesses robustness by performing a permutation test
(random shuffling of survival labels) to estimate empirical p-values for each gene's association with survival.
}
