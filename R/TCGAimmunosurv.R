#' Load Cancer Data for a Selected Cancer Type
#'
#' This function loads TCGA data for a specified cancer type, including RNA-Seq counts,
#' clinical data, and mutation data. If the data file exists, it is loaded directly.
#'
#' @param selected_cancer_type A character string specifying the TCGA cancer type
#'   (e.g., "BRCA" for breast cancer, "LUAD" for lung adenocarcinoma).
#'
#' @return None. This function loads the data into the global environment:
#'   - `tcga_count_data`: RNA-Seq count data in a `SummarizedExperiment` object.
#'   - `clinical_data`: Clinical data for survival analysis and metadata.
#'   - `mutation_data`: Mutation data from TCGA.
#'   - `sample_type`: A data frame of sample types.
#'
#' @details
#' If the data file exists locally, this function will:
#' - Load the RNA-Seq count data.
#' - Load the clinical data.
#' - Load the mutation data.
#'
#' @import TCGAbiolinks
#' @export
load_cancer_data <- function(selected_cancer_type) {
  rdata_filename <- paste0("tcga_", selected_cancer_type, "_data.RData")

  if (file.exists(rdata_filename)) {
    base::load(rdata_filename, envir = .GlobalEnv)
    base::cat("Data successfully loaded for analysis.\n")
  } else {
    stop("Error: Data file not found. Please download the data before loading.")
  }
}


#' Filter RNA Data
#'
#' This function filters RNA-Seq data stored in a `SummarizedExperiment` object.
#' It allows for filtering samples based on specific conditions, removing rows with
#' high NA values, replacing remaining NAs with the median, and removing low-count rows.
#'
#' @param rna_data A `SummarizedExperiment` object containing RNA-Seq count data.
#' @param condition_column A character string specifying the column name in `colData`
#'   that contains the conditions for filtering.
#' @param condition_levels A character vector of condition levels to retain in the data.
#' @param na_threshold A numeric value (default = 0.8) specifying the maximum allowed
#'   fraction of NAs in a row. Rows with a higher fraction of NAs will be removed.
#'
#' @return A filtered `SummarizedExperiment` object with updated counts, `rowRanges`,
#'   and `colData`.
#'
#'
#' @import SummarizedExperiment
#' @export
filter_rna_data <- function(rna_data, condition_column, condition_levels, na_threshold = 0.8) {
  # Check if condition_column exists
  if (!condition_column %in% colnames(SummarizedExperiment::colData(rna_data))) {
    stop(paste("Column", condition_column, "not found in colData of the provided SummarizedExperiment object."))
  }

  # Filter samples based on condition levels
  valid_samples <- SummarizedExperiment::colData(rna_data)[[condition_column]] %in% condition_levels
  rna_data <- rna_data[, valid_samples]

  # Update the condition column
  SummarizedExperiment::colData(rna_data)[[condition_column]] <- factor(SummarizedExperiment::colData(rna_data)[[condition_column]], levels = condition_levels)

  # Retrieve count data
  count_data <- SummarizedExperiment::assay(rna_data)

  # Remove rows with >80% NA
  na_fraction <- rowMeans(is.na(count_data))
  filtered_data <- count_data[na_fraction <= na_threshold, , drop = FALSE]

  # Replace remaining NAs with the median
  filtered_data <- t(apply(filtered_data, 1, function(x) {
    ifelse(is.na(x), median(x, na.rm = TRUE), x)
  }))

  # Remove low-count rows
  filtered_data <- filtered_data[rowSums(filtered_data) > 10, , drop = FALSE]

  # Update the metadata (rowRanges and colData)
  valid_rows <- rownames(filtered_data)
  valid_cols <- colnames(filtered_data)

  rna_data <- rna_data[valid_rows, valid_cols]

  # Ensure row and column names match
  rownames(filtered_data) <- rownames(rna_data)
  colnames(filtered_data) <- colnames(rna_data)

  # Assign the filtered data back to the assay
  SummarizedExperiment::assay(rna_data, withDimnames = FALSE) <- as.matrix(filtered_data)

  return(rna_data)
}


#' Perform Differential Expression Analysis
#'
#' This function calculates differentially expressed genes (DEGs) from RNA-Seq data
#' using the DESeq2 package. It applies thresholds for statistical significance (p-value)
#' and fold change to classify genes as upregulated, downregulated, or not significant.
#'
#' @param rna_data A `DESeqDataSet` object containing the RNA-Seq count matrix and metadata.
#' @param condition_column A string specifying the name of the column in the metadata
#'                         that defines the experimental conditions/groups.
#' @param alpha A numeric value specifying the significance threshold for the p-value (e.g., 0.05).
#' @param log2FC_threshold A numeric value specifying the log2 fold change threshold to
#'                         determine upregulated or downregulated genes.
#'
#' @return A data frame containing the DESeq2 results with an additional column, `significance`,
#'         indicating whether each gene is upregulated, downregulated, or not significant.
#'
#' @import DESeq2
#' @export
calculate_DEGs <- function(rna_data, condition_column, alpha, log2FC_threshold) {
  # Prepare DESeq2 dataset
  dds <- DESeq2::DESeqDataSet(rna_data, design = as.formula(paste0("~ ", condition_column)))
  dds <- DESeq2::DESeq(dds)

  # Extract results
  res <- DESeq2::results(dds, alpha = alpha,
                         contrast = c(condition_column,
                                      levels(SummarizedExperiment::colData(rna_data)[[condition_column]])[1],
                                      levels(SummarizedExperiment::colData(rna_data)[[condition_column]])[2]))

  # Add significance criteria
  res_df <- as.data.frame(res)
  res_df$significance <- with(res_df,
                              ifelse(pvalue <= alpha & log2FoldChange >= log2FC_threshold, "Upregulated",
                                     ifelse(pvalue <= alpha & log2FoldChange <= -log2FC_threshold, "Downregulated", "Not Significant")))

  return(res_df)
}


#' Plot Differentially Expressed Genes (DEGs) Results
#'
#' This function generates two plots to visualize the results of differential expression analysis:
#' a bar plot showing the sample distribution for different conditions and a volcano plot
#' for visualizing the differentially expressed genes (DEGs).
#'
#' @param res_df A data frame containing the results of the differential expression analysis,
#'   including columns for `log2FoldChange`, `padj` (adjusted p-values), and `significance`.
#'   The `significance` column should indicate whether a gene is "Upregulated",
#'   "Downregulated", or "Not Significant".
#' @param rna_data A `SummarizedExperiment` object containing RNA-Seq data with
#'   sample metadata accessible via `SummarizedExperiment::colData`.
#' @param condition_column A character string specifying the name of the column in
#'   `SummarizedExperiment::colData(rna_data)` that contains condition labels (e.g., "Control" and "Treatment").
#' @param log2FC_threshold A numeric value specifying the threshold for `log2FoldChange`.
#'   Default is 1.
#' @param alpha A numeric value specifying the significance threshold for adjusted p-values.
#'   Default is 0.05.
#'
#' @return None. This function outputs two plots:
#'   - A bar plot showing the distribution of samples across conditions.
#'   - A volcano plot visualizing the DEGs.
#'
#' @details
#' - The bar plot displays the sample count for each condition in the dataset, with labels
#'   for the number of samples above each bar.
#' - The volcano plot uses `log2FoldChange` and `-log10(padj)` to visualize the significance
#'   and magnitude of gene expression changes. Points are colored based on their significance
#'   status.
#' - Both plots use customizable font styles for better readability.
#'
#' @import ggplot2
#' @import SummarizedExperiment
#' @export
plot_DEG_results <- function(res_df, rna_data, condition_column, log2FC_threshold = 1, alpha = 0.05) {

  # Define font styles
  font_style <- ggplot2::theme(
    axis.text.x = ggplot2::element_text(size = 20, color = "black", face = "bold"),
    axis.text.y = ggplot2::element_text(size = 20, color = "black", face = "bold"),
    axis.title.x = ggplot2::element_text(size = 22, color = "black", face = "bold"),
    axis.title.y = ggplot2::element_text(size = 22, color = "black", face = "bold"),
    plot.title = ggplot2::element_text(size = 24, face = "bold", hjust = 0.5),
    legend.title = ggplot2::element_text(size = 20, face = "bold"),
    legend.text = ggplot2::element_text(size = 18)
  )

  # Bar plot for condition distribution
  sample_counts <- as.data.frame(table(SummarizedExperiment::colData(rna_data)[[condition_column]]))
  colnames(sample_counts) <- c("Condition", "Frequency")

  bar_plot <- ggplot2::ggplot(sample_counts, ggplot2::aes(x = Condition, y = Frequency)) +
    ggplot2::geom_bar(stat = "identity", fill = c("royalblue", "red"), width = 0.5) +
    ggplot2::geom_text(ggplot2::aes(label = Frequency), vjust = -0.3, size = 8, fontface = "bold") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = paste("Sample Distribution of", condition_column),
      x = "Condition",
      y = "Number of Samples"
    ) +
    font_style

  print(bar_plot)

  # Volcano plot
  volcano_plot <- ggplot2::ggplot(res_df, ggplot2::aes(x = log2FoldChange, y = -log10(padj))) +
    ggplot2::geom_point(ggplot2::aes(color = significance), alpha = 0.6) +
    ggplot2::scale_color_manual(
      values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "black"),
      name = "Direction"
    ) +
    ggplot2::guides(color = ggplot2::guide_legend(title = "Direction", override.aes = list(size = 3))) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Volcano Plot of Differentially Expressed Genes",
      x = "Log2(Fold Change)",
      y = "-Log10(P-Value)"
    ) +
    font_style

  print(volcano_plot)
}


#' Plot Unique Mutation Types Across Genes
#'
#' This function generates a bar plot to visualize the distribution of unique mutation types
#' across genes for a given cancer type based on the MAF (Mutation Annotation Format) data.
#'
#' @param cancer_maf A `MAF` object representing mutation data for a specific cancer type.
#'   This object is typically created using the `maftools::read.maf` function.
#'
#' @return A ggplot object showing the bar plot of unique mutation types and their total number
#'   of samples across all genes in the provided `MAF` object.
#'
#' @details
#' - The function extracts the mutation type counts from the `MAF` object using the
#'   `maftools::getGeneSummary` function.
#' - Mutation types with zero counts are filtered out.
#' - The mutation counts are aggregated across all genes for each mutation type.
#' - The resulting plot displays mutation types on the x-axis and the total number of
#'   samples on the y-axis.
#'
#' @import maftools
#' @import ggplot2
#' @import dplyr
#' @export
plot_unique_mutation_types <- function(cancer_maf) {
  # Extract gene summary
  gene_summary <- maftools::getGeneSummary(cancer_maf)

  # Extract mutation columns
  mutation_columns <- setdiff(
    colnames(gene_summary),
    c("Hugo_Symbol", "total", "MutatedSamples", "AlteredSamples")
  )

  # Reshape the data to focus only on mutation types and their counts across all genes
  mutation_df_list <- lapply(mutation_columns, function(mut_type) {
    data.frame(
      Mutation_Type = mut_type,
      NumSamples = gene_summary[[mut_type]]
    )
  })

  # Combine into a single DataFrame
  mutation_df <- do.call(rbind, mutation_df_list)

  # Filter out zero mutation counts and keep unique mutation types
  mutation_df <- mutation_df[mutation_df$NumSamples > 0, ]

  # Aggregate counts for each mutation type across all genes
  mutation_summary <- dplyr::group_by(mutation_df, Mutation_Type) %>%
    dplyr::summarise(NumSamples = sum(NumSamples), .groups = "drop")

  # Create the plot
  ggplot2::ggplot(mutation_summary, ggplot2::aes(x = Mutation_Type, y = NumSamples, fill = Mutation_Type)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = paste("All mutations present in", "selected_cancer"),
      x = "Mutation Type",
      y = "Total Number of Samples"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
}

#' Generate Mutation Summary for Genes of Interest
#'
#' This function generates a summary of mutation types and their counts for a given list of genes
#' from a MAF (Mutation Annotation Format) dataset.
#'
#' @param cancer_maf A MAF object created using the `maftools` package, containing mutation data.
#' @param genes_of_interest A character vector of gene names to summarize mutations for.
#'
#' @return A data table containing the mutation summary for the specified genes. The table includes:
#' \itemize{
#'   \item \code{Gene}: The gene name.
#'   \item \code{Mutation_Type}: The type of mutation.
#'   \item \code{NumSamples}: The number of samples with the specified mutation type.
#' }
#'
#' @details
#' - The function checks if each gene in `genes_of_interest` exists in the MAF data.
#' - If a gene is not found or has no non-synonymous variants, a warning is issued, and the gene is skipped.
#' - Mutation types and counts are extracted for each valid gene.
#'
#' @importFrom data.table data.table rbindlist
#' @importFrom maftools subsetMaf getGeneSummary
#' @export
generate_mutation_summary <- function(cancer_maf, genes_of_interest) {
  all_gene_mutation_counts <- list()  # Use a list to store individual gene data tables

  # Loop over each gene in the list of genes of interest
  for (gene in genes_of_interest) {
    # Check if the gene exists in the MAF data
    if (!gene %in% unique(cancer_maf@gene.summary$Hugo_Symbol)) {
      warning(paste("Gene", gene, "not found in the MAF data. Skipping."))
      next
    }

    # Filter the MAF data for the current gene's mutations
    gene_mutations <- maftools::subsetMaf(maf = cancer_maf, genes = gene)

    # Check if the gene has non-synonymous variants
    if (nrow(gene_mutations@data) == 0) {
      warning(paste("Gene", gene, "has no non-synonymous variants. Skipping."))
      next
    }

    # Get summary of mutations for the current gene
    gene_summary <- maftools::getGeneSummary(gene_mutations)

    # Extract mutation types and counts
    mutation_columns <- setdiff(colnames(gene_summary), c("Hugo_Symbol", "total", "MutatedSamples", "AlteredSamples"))
    gene_mutation_counts <- data.table::data.table(
      Gene = gene_summary$Hugo_Symbol[1],
      Mutation_Type = mutation_columns,
      NumSamples = as.numeric(gene_summary[1, ..mutation_columns])
    )

    # Store the data table in the list
    all_gene_mutation_counts[[gene]] <- gene_mutation_counts
  }

  # Combine all individual data tables into one
  final_table <- data.table::rbindlist(all_gene_mutation_counts, use.names = TRUE, fill = TRUE)

  return(final_table)
}


#' Get Unique Mutated Sample IDs for Selected Genes
#'
#' This function extracts unique tumor sample barcodes for selected genes with mutations from a MAF (Mutation Annotation Format) dataset.
#'
#' @param cancer_maf A MAF object created using the `maftools` package, containing mutation data.
#' @param genes_of_interest A character vector of gene names for which the mutated samples are to be extracted.
#'
#' @return A character vector containing the unique tumor sample barcodes of samples where mutations in the selected genes are present.
#'
#' @details
#' - The function subsets the MAF data to include only the mutations for the specified genes.
#' - It then extracts the unique tumor sample barcodes where mutations for those genes are present.
#'
#' @importFrom maftools subsetMaf
#' @export
get_mutated_samples <- function(cancer_maf, genes_of_interest) {
  # Subset MAF data for the genes of interest
  gene_mutations <- maftools::subsetMaf(maf = cancer_maf, genes = genes_of_interest)

  # Extract unique sample IDs
  return(unique(gene_mutations@data$Tumor_Sample_Barcode))
}


#' Plot Mutation Counts for Each Gene
#'
#' This function generates separate bar plots showing the mutation counts for each gene, based on mutation types and sample counts.
#'
#' @param mutation_counts A data.frame or data.table containing the mutation data, including columns for:
#' - `Gene`: The gene name.
#' - `Mutation_Type`: The type of mutation (e.g., missense, nonsense).
#' - `NumSamples`: The number of samples with that mutation type for the given gene.
#'
#' @return A series of bar plots, one for each gene, showing the count of mutations by mutation type.
#'
#' @details
#' The function iterates over each unique gene in the `mutation_counts` data and creates a bar plot for that gene. Each plot represents the count of mutations per mutation type for the gene, with the bars filled according to the mutation type.
#'
#' @importFrom ggplot2 ggplot geom_bar theme_minimal labs element_text unit aes
#' @export
plot_mutation_counts <- function(mutation_counts) {
  # Ensure the Gene column is treated as a factor to iterate over unique genes
  unique_genes <- unique(mutation_counts$Gene)

  # Iterate through each gene and create a separate plot
  for (gene in unique_genes) {
    # Filter data for the current gene
    gene_data <- subset(mutation_counts, Gene == gene)

    # Create the plot for the current gene
    p <- ggplot2::ggplot(gene_data, ggplot2::aes(x = Mutation_Type, y = NumSamples, fill = Mutation_Type)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = paste("Mutation Counts for Gene:", gene),
        x = "Mutation Type",
        y = "Number of Samples"
      ) +
      ggplot2::theme(
        # Customize axis text and labels
        axis.text.x = ggplot2::element_text(color = "black", face = "bold", size = 18, angle = 45, hjust = 1),
        axis.text.y = ggplot2::element_text(color = "black", face = "bold", size = 18),
        axis.title.x = ggplot2::element_text(color = "black", face = "bold", size = 20),
        axis.title.y = ggplot2::element_text(color = "black", face = "bold", size = 20),

        # Customize plot title
        plot.title = ggplot2::element_text(color = "black", face = "bold", size = 22, hjust = 0.5),

        # Customize legend
        legend.text = ggplot2::element_text(color = "black", size = 14, face = "bold"),
        legend.title = ggplot2::element_text(color = "black", size = 16, face = "bold"),
        legend.key.size = grid::unit(0.8, "cm")  # Adjust the size of the legend key
      )

    # Print the plot to the current graphics device
    print(p)
  }
}


#' Style the Mutation Counts Table
#'
#' This function applies styling to a mutation counts table, enhancing the readability and presentation of the data using the `kable` and `kableExtra` packages.
#'
#' @param mutation_counts A data frame or data table containing the mutation counts, with columns for:
#' - `Gene`: The gene name.
#' - `Mutation_Type`: The type of mutation.
#' - `NumSamples`: The number of samples with that mutation type for the given gene.
#'
#' @return A styled table with mutation counts for selected genes, formatted for better presentation.
#'
#' @details
#' The function uses `kableExtra::kbl()` to create a table and then applies various styling options:
#' - `striped`, `hover`, `condensed`, and `responsive` bootstrap options for table formatting.
#' - The first column (Gene) is made bold with a fixed width of 3 cm.
#' - The second column (Mutation Type) has a width of 4 cm.
#' - The third column (NumSamples) is colored blue.
#' - A footnote is added to summarize the content of the table.
#'
#' @importFrom kableExtra kbl kable_styling column_spec footnote
#' @importFrom dplyr %>%
#' @export
style_mutation_table <- function(mutation_counts) {
  mutation_counts %>%
    kableExtra::kbl(caption = "Mutation Counts Across Samples for Genes of Interest") %>%
    kableExtra::kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive")) %>%
    kableExtra::column_spec(1, bold = TRUE, width = "3cm") %>%
    kableExtra::column_spec(2, width = "4cm") %>%
    kableExtra::column_spec(3, width = "3cm", color = "blue") %>%
    kableExtra::footnote(general = "This table summarizes the types of mutations and the number of samples in which they occur for the selected genes.")
}


#' Standardize TCGA Sample IDs
#'
#' This function standardizes TCGA sample IDs by extracting the first three components of the ID before the first hyphen ("-"). If the ID format is invalid (does not contain "-"), it returns `NA`.
#'
#' @param sample_ids A character vector of TCGA sample IDs.
#'
#' @return A character vector of standardized TCGA sample IDs, where each ID is truncated to the first three components, or `NA` if the format is invalid.
#'
#' @details
#' The function ensures that all TCGA sample IDs are standardized to the format containing only the first three components separated by hyphens. If the input contains IDs in an incorrect format (i.e., not containing "-"), the function returns `NA` for those IDs.
#'
#' @importFrom stringr str_split str_detect
#' @export
standardize_tcga_ids <- function(sample_ids) {
  sample_ids <- as.character(sample_ids)  # Ensure input is character
  sapply(sample_ids, function(x) {
    if (stringr::str_detect(x, "-")) {  # Check if the string contains "-"
      paste(stringr::str_split(x, "-")[[1]][1:3], collapse = "-")
    } else {
      NA  # Return NA for invalid formats
    }
  })
}


#' Preprocess RNA-Seq and Clinical Data
#'
#' This function preprocesses RNA-Seq count data and clinical data, filtering both to retain only the samples that have mutation data and are present in both the RNA-Seq and clinical datasets. It also standardizes the sample IDs to a common format.
#'
#' @param rna_data A `SummarizedExperiment` or similar object containing RNA-Seq count data. It must have a count matrix and associated gene metadata.
#' @param clinical_data A data frame containing clinical data with a column `submitter_id` that matches sample IDs in `rna_data`.
#' @param extracted_mutated_sample_ids A character vector of sample IDs that have mutations of interest.
#'
#' @return A list containing two elements:
#'   - `count_matrix_filtered`: A matrix of RNA-Seq count data for samples that have mutation data and are present in both datasets.
#'   - `clinical_data_filtered`: A data frame of clinical data for samples that have mutation data and are present in both datasets.
#'
#' @details
#' The function filters the RNA-Seq count matrix and the clinical data to include only the samples that are present in both the mutation data and clinical data. The sample IDs are standardized to the first three components (using hyphen delimiters) for consistency.
#'
#' @importFrom SummarizedExperiment assay rowData
#' @importFrom dplyr filter
#' @importFrom data.table as.data.table
#' @export
preprocess_data <- function(rna_data, clinical_data, extracted_mutated_sample_ids) {

  count_matrix <- SummarizedExperiment::assay(rna_data, "unstranded")
  gene_metadata_dt <- as.data.frame(SummarizedExperiment::rowData(rna_data))

  clinical_submitter_ids <- data.frame(submitter_id = clinical_data$submitter_id)

  samples_with_mutation <- intersect(extracted_mutated_sample_ids, clinical_submitter_ids$submitter_id)

  clinical_data_filtered <- dplyr::filter(clinical_data, submitter_id %in% samples_with_mutation)

  clinical_data_filtered <- data.table::as.data.table(clinical_data_filtered)

  colnames <- data.frame(colnames(count_matrix))
  extracted_colnames <- sapply(colnames$colnames.count_matrix., function(x) {
    paste(strsplit(x, "-")[[1]][1:3], collapse = "-")
  })
  colnames$trimmed <- extracted_colnames

  colnames_filtered <- dplyr::filter(colnames, trimmed %in% samples_with_mutation)

  count_matrix_filtered <- data.table::as.data.table(count_matrix[, colnames(count_matrix) %in% colnames_filtered$colnames.count_matrix.])
  count_matrix_filtered <- as.matrix(count_matrix_filtered)

  return(list(count_matrix_filtered = count_matrix_filtered, clinical_data_filtered = clinical_data_filtered))
}


#' Perform Survival Analysis on RNA-Seq Data
#'
#' This function performs survival analysis on RNA-Seq data by grouping samples into two strata (HIGH/LOW) based on the median expression of each gene. It then performs Cox regression analysis for each gene and returns the hazard ratios and p-values.
#'
#' @param count_matrix_filtered A matrix of RNA-Seq count data, where rows represent genes and columns represent samples.
#' @param gene_metadata_dt A data frame containing metadata for the genes, including a `gene_name` column.
#' @param clinical_data_filtered A data frame containing clinical data, including `submitter_id`, `overall_survival`, and `deceased` columns.
#'
#' @return A data frame with the results of the survival analysis for each gene, including:
#'   - `gene_name`: The name of the gene.
#'   - `hazard_ratio`: The hazard ratio for the gene.
#'   - `pvalue`: The p-value from the Cox regression for the gene.
#'   - `significant`: Whether the gene is significantly associated with survival (p-value ≤ 0.05).
#'
#' @importFrom survival coxph Surv
#' @importFrom data.table as.data.table melt
#' @importFrom utils globalVariables
#' @export
perform_survival_analysis <- function(count_matrix_filtered, gene_metadata_dt, clinical_data_filtered) {

  # Identify duplicate gene names in gene_metadata_dt
  dup_genes <- gene_metadata_dt$gene_name[duplicated(gene_metadata_dt$gene_name)]
  unique_dup_genes <- unique(dup_genes)

  # Calculate mean for each set of duplicate rows in count_matrix_filtered
  rows_to_keep <- integer(0)
  for (gene in unique_dup_genes) {
    gene_rows <- which(gene_metadata_dt$gene_name == gene)
    mean_row <- colMeans(count_matrix_filtered[gene_rows, ], na.rm = TRUE)

    # Replace the first occurrence in count_matrix_filtered with the mean row
    count_matrix_filtered[gene_rows[1], ] <- mean_row
    rows_to_keep <- c(rows_to_keep, gene_rows[1])
  }

  # Add rows that are not duplicates
  non_duplicate_rows <- which(!gene_metadata_dt$gene_name %in% unique_dup_genes)
  final_rows_to_keep <- sort(c(rows_to_keep, non_duplicate_rows))

  # Filter count_matrix_filtered and gene_metadata_dt to keep only desired rows
  count_matrix_filtered <- count_matrix_filtered[final_rows_to_keep, ]
  gene_metadata_dt <- gene_metadata_dt[final_rows_to_keep, ]

  # Set gene names as rownames in count_matrix_filtered
  rownames(count_matrix_filtered) <- gene_metadata_dt$gene_name
  count_matrix_filtered_dt <- data.table::as.data.table(count_matrix_filtered, keep.rownames = "gene_name")

  # Melting to long format and calculate strata for each sample by gene name
  strata_data <- data.table::melt(count_matrix_filtered_dt, id.vars = "gene_name", variable.name = "case_id", value.name = "counts")[
    , strata := ifelse(counts >= median(counts, na.rm = TRUE), "HIGH", "LOW"), by = gene_name
  ]

  strata_data[, strata := factor(strata, levels = c("LOW", "HIGH"))]

  # Remove genes with only one strata level
  single_strata_genes <- strata_data[, .N, by = .(gene_name, strata)][, .N, by = gene_name][N == 1]
  strata_data <- strata_data[gene_name %in% single_strata_genes$gene_name == FALSE]

  # Merge with clinical data
  strata_data[, case_id := gsub('-01.*', '', case_id)]
  clinical_data_filtered[, submitter_id := gsub('-01.*', '', clinical_data_filtered$submitter_id)]

  assign("strata_data", strata_data, envir = .GlobalEnv)

  # Perform Cox regression for each gene
  results_list <- vector("list", length(unique(strata_data$gene_name)))

  gene_names <- unique(strata_data$gene_name)
  for (i in seq_along(gene_names)) {
    gene <- gene_names[i]
    gene_data_long <- strata_data[gene_name == gene]

    # Merging with clinical data on case_id
    gene_data_long$case_id <- gsub('-01.*', '', gene_data_long$case_id)
    gene_data_long <- merge(gene_data_long, clinical_data_filtered, by.x = 'case_id', by.y = 'submitter_id')

    # Performing Cox regression
    cox_fit <- tryCatch({
      survival::coxph(survival::Surv(overall_survival, deceased) ~ strata, data = gene_data_long)
    }, error = function(e) {
      cat("Cox model failed for gene:", gene, "-", e$message, "\n")
      NULL
    })

    if (is.null(cox_fit)) next

    hr <- exp(coef(cox_fit))
    pvalue <- summary(cox_fit)$coef["strataHIGH", "Pr(>|z|)"]

    results_list[[i]] <- data.frame(gene_name = gene, hazard_ratio = hr, pvalue = pvalue)
  }

  # Combining results into a single data frame
  results <- do.call(rbind, results_list)
  results$significant <- ifelse(results$pvalue <= 0.05, "Yes", "No")

  return(results)
}

<<<<<<< HEAD
<<<<<<< HEAD
#' Survival Analysis with Permutation Test for Robust Gene Selection
#'
#' This function performs survival analysis for each gene using a Cox proportional hazards model,
#' stratifies expression as HIGH/LOW per gene, and assesses robustness by performing a permutation test
#' (random shuffling of survival labels) to estimate empirical p-values for each gene's association with survival.
#'
#' @param count_matrix_filtered A numeric matrix of gene expression (genes x samples) with gene IDs as rows and sample IDs as columns.
#' @param gene_metadata_dt A data.table containing at least a column named `gene_name` corresponding to the genes in `count_matrix_filtered`.
#' @param clinical_data_filtered A data.table with clinical data, including columns `submitter_id`, `overall_survival`, and `deceased`.
#' @param n_permutations Number of permutations to perform for the permutation test (default: 1000).
#'
#'
#' Genes are reported as significant if their empirical permutation p-value (`perm_pval`) is less than or equal to 0.05.
#'
#' @return
#' A data.frame with one row per gene, containing:
#'   \describe{
#'     \item{gene_name}{Gene name}
#'     \item{hazard_ratio}{Hazard ratio (HIGH vs LOW expression)}
#'     \item{pvalue}{Wald p-value from the Cox model}
#'     \item{logrank_stat}{Observed log-rank (likelihood ratio) statistic}
#'     \item{perm_pval}{Empirical p-value from permutation test}
#'     \item{significant}{"Yes" if perm_pval <= 0.05, otherwise "No"}
#'   }
#'
#' @import data.table
#' @importFrom survival Surv coxph
#' @export
perform_survival_analysis_with_permutation <- function(
    count_matrix_filtered,
    gene_metadata_dt,
    clinical_data_filtered,
    n_permutations = 1000
) {
=======
perform_survival_analysis_with_permutation <- function(count_matrix_filtered, gene_metadata_dt, clinical_data_filtered, n_permutations = 1000) {
>>>>>>> 7633678ebafb876ed4b313c463b7f0caf11f23e4
  
  # Handle duplicate genes
  dup_genes <- gene_metadata_dt$gene_name[duplicated(gene_metadata_dt$gene_name)]
  unique_dup_genes <- unique(dup_genes)
  rows_to_keep <- integer(0)
  for (gene in unique_dup_genes) {
    gene_rows <- which(gene_metadata_dt$gene_name == gene)
    mean_row <- colMeans(count_matrix_filtered[gene_rows, ], na.rm = TRUE)
    count_matrix_filtered[gene_rows[1], ] <- mean_row
    rows_to_keep <- c(rows_to_keep, gene_rows[1])
  }
  non_duplicate_rows <- which(!gene_metadata_dt$gene_name %in% unique_dup_genes)
  final_rows_to_keep <- sort(c(rows_to_keep, non_duplicate_rows))
  count_matrix_filtered <- count_matrix_filtered[final_rows_to_keep, ]
  gene_metadata_dt <- gene_metadata_dt[final_rows_to_keep, ]
  rownames(count_matrix_filtered) <- gene_metadata_dt$gene_name
  
  count_matrix_filtered_dt <- data.table::as.data.table(count_matrix_filtered, keep.rownames = "gene_name")
  
  # Stratify expression
  strata_data <- melt(count_matrix_filtered_dt, id.vars = "gene_name", variable.name = "case_id", value.name = "counts")[
    , strata := ifelse(counts >= median(counts, na.rm = TRUE), "HIGH", "LOW"), by = gene_name]
  strata_data[, strata := factor(strata, levels = c("LOW", "HIGH"))]
  
  # Remove unstratifiable genes
  single_strata_genes <- strata_data[, .N, by = .(gene_name, strata)][, .N, by = gene_name][N == 1]
  strata_data <- strata_data[!gene_name %in% single_strata_genes$gene_name]
  
  # Merge with clinical
  strata_data[, case_id := gsub("-01.*", "", case_id)]
  clinical_data_filtered[, submitter_id := gsub("-01.*", "", submitter_id)]
  setnames(clinical_data_filtered, old = c("submitter_id"), new = c("case_id"))
  
  results_list <- list()
  gene_names <- unique(strata_data$gene_name)
  
  for (gene in gene_names) {
    gene_data <- strata_data[gene_name == gene]
    merged_data <- merge(gene_data, clinical_data_filtered, by = "case_id")
    
    if (length(unique(merged_data$strata)) < 2) next
    
    # Observed Cox model
    cox_obs <- tryCatch({
<<<<<<< HEAD
      survival::coxph(survival::Surv(overall_survival, deceased) ~ strata, data = merged_data)
=======
      coxph(Surv(overall_survival, deceased) ~ strata, data = merged_data)
>>>>>>> 7633678ebafb876ed4b313c463b7f0caf11f23e4
    }, error = function(e) NULL)
    
    if (is.null(cox_obs)) next
    
    obs_lr <- summary(cox_obs)$logtest[1]  # Chi-square log-rank stat
    hr <- exp(coef(cox_obs))
    pval <- summary(cox_obs)$coef["strataHIGH", "Pr(>|z|)"]
    
    # Permutation test
    perm_stats <- numeric(n_permutations)
    for (i in seq_len(n_permutations)) {
      perm_data <- merged_data
      perm_data$deceased <- sample(perm_data$deceased)
      perm_data$overall_survival <- sample(perm_data$overall_survival)
      
      cox_perm <- tryCatch({
<<<<<<< HEAD
        survival::coxph(survival::Surv(overall_survival, deceased) ~ strata, data = perm_data)
=======
        coxph(Surv(overall_survival, deceased) ~ strata, data = perm_data)
>>>>>>> 7633678ebafb876ed4b313c463b7f0caf11f23e4
      }, error = function(e) NULL)
      
      perm_stats[i] <- if (!is.null(cox_perm)) summary(cox_perm)$logtest[1] else NA
    }
    
    perm_stats <- perm_stats[!is.na(perm_stats)]
    perm_pval <- mean(perm_stats >= obs_lr)
    
    results_list[[gene]] <- data.frame(
      gene_name = gene,
      hazard_ratio = hr,
      pvalue = pval,
      logrank_stat = obs_lr,
      perm_pval = perm_pval
    )
  }
  
  final_results <- do.call(rbind, results_list)
  final_results$significant <- ifelse(final_results$perm_pval <= 0.05, "Yes", "No")
  
  return(final_results)
}
<<<<<<< HEAD
=======

>>>>>>> 7633678ebafb876ed4b313c463b7f0caf11f23e4

=======
>>>>>>> 7102884e2e691fd8c4486e47876a4be4f2d38b05
#' Plot Survival Analysis Results
#'
#' This function generates three types of plots to visualize the survival analysis results:
#' a bar plot showing the number of significant vs. non-significant genes,
#' a bar plot showing the count of significant genes categorized by hazard ratio,
#' and a volcano plot representing the relationship between hazard ratio and p-value.
#'
#' @param results A data frame containing the results of the survival analysis. The data frame should have the following columns:
#' \describe{
#'   \item{gene_name}{Character string, the name of the gene.}
#'   \item{hazard_ratio}{Numerical value representing the hazard ratio.}
#'   \item{pvalue}{Numerical value representing the p-value of the gene.}
#'   \item{significant}{Character string, "Yes" if the gene is statistically significant, otherwise "No".}
#' }
#'
#' @return A list containing three ggplot objects:
#' \describe{
#'   \item{SignificancePlot}{A bar plot showing the count of significant vs. non-significant genes.}
#'   \item{HRPlot}{A bar plot showing the count of significant genes categorized by hazard ratio.}
#'   \item{VolcanoPlot}{A volcano plot showing the relationship between log hazard ratio and -log10 p-value.}
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_bar geom_text geom_point scale_fill_manual scale_color_manual labs theme_minimal theme element_text element_blank element_line
#' @importFrom dplyr group_by summarise mutate filter arrange slice
#' @export
plot_results <- function(results) {

  # Data Preparation
  significance_counts <- results %>%
    dplyr::group_by(significant) %>%
    dplyr::summarise(count = dplyr::n()) %>%
    dplyr::mutate(significant = factor(significant, levels = c("Yes", "No")))

  significant_hr_counts <- results %>%
    dplyr::filter(significant == "Yes") %>%
    dplyr::mutate(Category = ifelse(hazard_ratio > 1, "HR > 1", "HR < 1")) %>%
    dplyr::group_by(Category) %>%
    dplyr::summarise(count = dplyr::n())

  # Plot 1: Significant vs Non-Significant
  plot1 <- ggplot2::ggplot(significance_counts, ggplot2::aes(x = significant, y = count, fill = significant)) +
    ggplot2::geom_bar(stat = "identity", width = 0.5) +
    ggplot2::scale_fill_manual(values = c("Yes" = "#FF6666", "No" = "#BFBFBF")) +
    ggplot2::labs(title = "Number of Significant vs Non-Significant Genes", x = "Significance", y = "Count") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = ggplot2::element_text(size = 20, color = "black", face = "bold"),
      axis.text.y = ggplot2::element_text(size = 20, color = "black", face = "bold"),
      axis.title.x = ggplot2::element_text(size = 22, color = "black", face = "bold"),
      axis.title.y = ggplot2::element_text(size = 22, color = "black", face = "bold"),
      plot.title = ggplot2::element_text(size = 24, face = "bold", hjust = 0.5)
    )

  # Plot 2: Significant Genes by Hazard Ratio
  plot2 <- ggplot2::ggplot(significant_hr_counts, ggplot2::aes(x = Category, y = count, fill = Category)) +
    ggplot2::geom_bar(stat = "identity", width = 0.5) +
    ggplot2::geom_text(ggplot2::aes(label = count), vjust = 1, size = 8, color = "black", fontface = "bold", position = ggplot2::position_stack(vjust = 0.5)) +
    ggplot2::scale_fill_manual(values = c("HR > 1" = "#66CCFF", "HR < 1" = "#FFCC66")) +
    ggplot2::labs(
      title = "Significant Genes by Hazard Ratio",
      x = "Hazard Ratio Category",
      y = "Count"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = ggplot2::element_text(size = 20, color = "black", face = "bold"),
      axis.text.y = ggplot2::element_text(size = 20, color = "black", face = "bold"),
      axis.title.x = ggplot2::element_text(size = 22, color = "black", face = "bold"),
      axis.title.y = ggplot2::element_text(size = 22, color = "black", face = "bold"),
      plot.title = ggplot2::element_text(size = 24, face = "bold", hjust = 0.5)
    )

  # Plot 3: Volcano Plot
  plot3 <- ggplot2::ggplot(results, ggplot2::aes(x = log(hazard_ratio), y = -log10(pvalue))) +
    ggplot2::geom_point(ggplot2::aes(color = significant)) +
    ggplot2::scale_color_manual(values = c("Yes" = "#FF6666", "No" = "#BFBFBF")) +
    ggplot2::labs(title = "Volcano Plot", x = "Log(Hazard Ratio)", y = "-Log10(P-value)") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "right",
      legend.title = ggplot2::element_text(size = 20, face = "bold"),
      legend.text = ggplot2::element_text(size = 18),
      axis.text.x = ggplot2::element_text(size = 20, color = "black", face = "bold"),
      axis.text.y = ggplot2::element_text(size = 20, color = "black", face = "bold"),
      axis.title.x = ggplot2::element_text(size = 22, color = "black", face = "bold"),
      axis.title.y = ggplot2::element_text(size = 22, color = "black", face = "bold"),
      plot.title = ggplot2::element_text(size = 24, face = "bold", hjust = 0.5)
    )

  # Print Separate Plots
  print(plot1)
  print(plot2)
  print(plot3)

  return(list(SignificancePlot = plot1, HRPlot = plot2, VolcanoPlot = plot3))
}


#' Plot Kaplan-Meier Curves for Significant Genes
#'
#' This function generates Kaplan-Meier survival curves for each significant gene based on its association with the clinical data.
#' The function filters the `strata_data` for the significant genes, merges it with clinical data, and fits the Kaplan-Meier model
#' for each gene. A survival plot is then created for each gene showing survival probability over time.
#'
#' @param strata_data A data frame containing gene expression data for multiple cases. This data frame should include at least the columns:
#' \describe{
#'   \item{gene_name}{The name of the gene.}
#'   \item{case_id}{The identifier for the case/sample.}
#' }
#' @param clinical_data_filtered A data frame containing clinical data filtered for the relevant cases. This data frame should include at least the columns:
#' \describe{
#'   \item{submitter_id}{The identifier for the case/sample.}
#'   \item{overall_survival}{The survival time of the patient (in days).}
#'   \item{deceased}{A binary indicator of whether the patient is deceased (1 = deceased, 0 = alive).}
#'   \item{strata}{A factor or categorical variable used to stratify the survival analysis (e.g., treatment groups).}
#' }
#' @param significant_genes A character vector containing the names of significant genes to plot. The function will generate Kaplan-Meier curves for each gene in this vector.
#'
#' @return A list of Kaplan-Meier plots (`ggsurvplot` objects) for each significant gene. The plots are printed individually for each gene.
#'
#' @importFrom survival Surv survfit
#' @importFrom survminer ggsurvplot
#' @importFrom dplyr filter mutate select
#' @export
plot_km_curves <- function(strata_data, clinical_data_filtered, significant_genes) {

  # Filter strata_data for significant genes
  filtered_strata_data <- dplyr::filter(strata_data, gene_name %in% significant_genes)

  for (gene in significant_genes) {
    cat("\nProcessing gene:", gene, "\n")

    # Subset data for the current gene and merge with clinical data
    gene_data_long <- merge(
      dplyr::filter(filtered_strata_data, gene_name == gene),
      clinical_data_filtered,
      by.x = "case_id",
      by.y = "submitter_id",
      all.x = TRUE
    )

    # Validate required columns and skip if issues are found
    if (!all(c("overall_survival", "deceased", "strata") %in% colnames(gene_data_long))) {
      cat("Missing required columns for gene:", gene, "- Skipping\n")
      next
    }

    # Convert 'deceased' to numeric (TRUE = 1, FALSE = 0)
    gene_data_long$deceased <- as.numeric(gene_data_long$deceased)

    # Check if 'strata' has at least two levels
    if (length(unique(gene_data_long$strata)) < 2) {
      cat("'strata' column has insufficient levels for gene:", gene, "- Skipping\n")
      next
    }

    # Fit the survival object and Kaplan-Meier model
    surv_object <- survival::Surv(gene_data_long$overall_survival, gene_data_long$deceased)
    assign("surv_object", surv_object, envir= .GlobalEnv)
    km_fit <- survival::survfit(surv_object ~ strata, data = gene_data_long)

    # Plot the Kaplan-Meier curve
    plot <- survminer::ggsurvplot(
      km_fit,
      data = gene_data_long,
      pval = TRUE,
      risk.table = TRUE,
      ggtheme = ggplot2::theme_minimal(),
      title = paste("Kaplan-Meier Curve for", gene),
      xlab = "Time (days)",
      ylab = "Survival Probability",
      risk.table.height = 0.25
    )
    print(plot)
  }
}


#' Map Differentially Expressed Genes (DEGs) with Gene Metadata
#'
#' This function takes in the results of differential gene expression analysis (`deg_results`), a gene metadata table (`gene_metadata_dt`),
#' and filters for upregulated and downregulated genes based on user-defined thresholds for log fold change and adjusted p-value.
#' It then returns a merged table with gene names and the associated differential expression data.
#'
#' @param deg_results A data frame containing differential expression results. The data frame should include at least the following columns:
#' \describe{
#'   \item{log2FoldChange}{The log2 fold change values for the genes.}
#'   \item{padj}{The adjusted p-values for the genes.}
#'   \item{ensemble_id}{A unique identifier for each gene (can be rownames if not present as a column).}
#' }
#' @param gene_metadata_dt A data frame containing gene metadata. This should include at least the following columns:
#' \describe{
#'   \item{ensemble_id}{The unique identifier for the gene, which matches the `ensemble_id` in `deg_results`.}
#'   \item{gene_name}{The gene names corresponding to the `ensemble_id`.}
#' }
#' @param logFC_threshold A numeric value specifying the threshold for log2 fold change. Genes with `log2FoldChange > logFC_threshold` are considered upregulated,
#' and genes with `log2FoldChange < -logFC_threshold` are considered downregulated.
#' @param padj_threshold A numeric value specifying the adjusted p-value threshold. Only genes with `padj < padj_threshold` will be selected.
#'
#' @return A data frame containing the upregulated and downregulated genes based on the given thresholds, with columns:
#' \describe{
#'   \item{ensemble_id}{The unique identifier for the gene.}
#'   \item{gene_name}{The gene name associated with the `ensemble_id`.}
#'   \item{log2FoldChange}{The log2 fold change values for the gene.}
#'   \item{padj}{The adjusted p-values for the gene.}
#' }
#' @importFrom dplyr left_join filter select mutate
#' @export
mapping_genes <- function(deg_results, gene_metadata_dt, logFC_threshold, padj_threshold) {

  # Step 1: Merge the DEG results with gene metadata
  # Ensure ensemble_id is present
  gene_metadata_dt <- dplyr::mutate(gene_metadata_dt, ensemble_id = rownames(gene_metadata_dt))
  deg_results <- dplyr::mutate(deg_results, ensemble_id = rownames(deg_results))

  # Merge the data frames by ensemble_id using dplyr::left_join for better performance
  merged_data <- dplyr::left_join(deg_results, gene_metadata_dt %>% dplyr::select(ensemble_id, gene_name), by = "ensemble_id")

  # Remove rows with NA values
  merged_data <- na.omit(merged_data)

  # Step 2: Filter for upregulated and downregulated genes based on the log fold change threshold
  upregulated_genes <- dplyr::filter(merged_data, log2FoldChange > logFC_threshold, padj < padj_threshold)
  downregulated_genes <- dplyr::filter(merged_data, log2FoldChange < -logFC_threshold, padj < padj_threshold)

  # Combine upregulated and downregulated genes
  selected_genes <- dplyr::bind_rows(upregulated_genes, downregulated_genes)

  return(selected_genes)
}


#' Extract Significant Genes with Hazard Ratio Greater Than 1 and Save to CSV
#'
#' This function filters significant genes with a hazard ratio (HR) greater than 1
#' from a given results DataFrame. It returns a DataFrame of the filtered genes
#' and their hazard ratios and saves the results to a `.csv` file.
#'
#' @param results A DataFrame containing gene analysis results with at least
#'                the following columns: \code{significant}, \code{hazard_ratio},
#'                and \code{gene_name}.
#' @param output_file Character. The name of the output `.csv` file to save
#'                    the results. Default is "significant_genes.csv".
#' @return A DataFrame with two columns:
#'         \itemize{
#'           \item \code{gene_name}: Names of significant genes with HR > 1.
#'           \item \code{hazard_ratio}: The hazard ratios of the significant genes.
#'         }
#' @importFrom dplyr filter select
#' @importFrom utils write.csv
#' @export
#'

extract_genes_hr_gt1 <- function(results, output_file = "significant_genes.csv") {
  if (!all(c("significant", "hazard_ratio", "gene_name") %in% colnames(results))) {
    stop("The input DataFrame must contain 'significant', 'hazard_ratio', and 'gene_name' columns.")
  }

  # Filter for significant genes with hazard ratio > 1
  significant_genes <- dplyr::filter(results, significant == "Yes", hazard_ratio > 1) %>%
    dplyr::select(gene_name, hazard_ratio)

  # Save the filtered results to a CSV file
  utils::write.csv(significant_genes, file = output_file, row.names = FALSE)

  return(significant_genes)
}

#' Load Seurat Object for Single-Cell Analysis
#'
#' This function loads a Seurat object from a specified directory based on a given cancer type.
#' The Seurat object will be used to perform further single-cell RNA-seq analysis.
#'
#' @param directory A character string specifying the directory where the Seurat object file is located.
#' @param constant_prefix A character string specifying the constant prefix for the Seurat object file name. The default is "combined_seurat_".
#' @param selected_cancer A character string specifying the cancer type for which the Seurat object file will be loaded. The cancer type will be appended to the constant prefix.
#'
#' @return A Seurat object that has been loaded from the specified file.
#'
#' @details This function constructs the file name for the Seurat object based on the provided `constant_prefix` and `selected_cancer`. It checks if the file exists in the specified directory and loads the Seurat object if the file is found. If the file does not exist, an error message is displayed.
#'@export

load_seurat_object <- function(directory, constant_prefix = "combined_seurat_", selected_cancer) {
  # Construct the full filename
  file_name <- paste0(constant_prefix, selected_cancer, ".rds")
  file_path <- file.path(directory, file_name)

  # Check if the file exists
  if (!file.exists(file_path)) {
    stop(paste("Error: File", file_name, "does not exist in the specified directory:", directory))
  }

  # Load the Seurat object from the .rds file
  seurat_object <- readRDS(file_path)
  message(paste("Seurat object loaded successfully from", file_name))

  return(seurat_object)
}


#' Filter Seurat Object Based on Quality Control Criteria
#'
#' This function filters a Seurat object based on the number of detected features (genes) and the percentage of mitochondrial genes.
#' It ensures that only high-quality cells are kept for downstream analysis.
#'
#' @param seurat_obj A Seurat object containing single-cell RNA-seq data.
#' @param min_features An integer specifying the minimum number of detected features (genes) required for a cell to be retained. Default is 200.
#' @param max_features An integer specifying the maximum number of detected features (genes) allowed for a cell. Default is 2500.
#' @param max_percent_mt A numeric value specifying the maximum percentage of mitochondrial genes allowed for a cell. Default is 5.
#'
#' @return A filtered Seurat object, where cells meeting the quality control criteria are retained.
#'
#' @details The function adds a new metadata column to the Seurat object, `percent.mt`, which represents the percentage of mitochondrial genes expressed in each cell. The cells are then filtered based on the following criteria:
#'
#' - The number of detected features (genes) must be between `min_features` and `max_features`.
#' - The percentage of mitochondrial gene expression must be below `max_percent_mt`.
#'
#' This helps to remove cells with poor quality (e.g., low gene count or high mitochondrial gene expression) that could skew downstream analyses.
#'
#' @importFrom Seurat PercentageFeatureSet
#' @export
#' 
filter_seurat_object <- function(seurat_obj, min_features = 200, max_features = 2500, max_percent_mt = 5) {
  # Calculate mitochondrial percentage and add to metadata
  seurat_obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(seurat_obj, pattern = "^MT-")

  # Filter cells based on the criteria
  seurat_obj <- subset(
    seurat_obj,
    subset = nFeature_RNA > min_features &
      nFeature_RNA < max_features &
      percent.mt < max_percent_mt
  )

  return(seurat_obj)
}

#' Perform Principal Component Analysis (PCA) on Seurat Object
#'
#' This function performs PCA on a Seurat object after normalizing the data, identifying variable features, and scaling the data.
#' The resulting PCA components can be used for dimensionality reduction and visualization in downstream analysis.
#'
#' @param seurat_obj A Seurat object containing single-cell RNA-seq data.
#'
#' @return A Seurat object with PCA results stored in the `pca` assay, ready for further analysis (e.g., visualization, clustering).
#'
#' @details The function performs the following steps:
#' \itemize{
#'   \item Normalizes the data using `Seurat::NormalizeData()`.
#'   \item Identifies highly variable features with `Seurat::FindVariableFeatures()`.
#'   \item Scales the data using `Seurat::ScaleData()`.
#'   \item Runs PCA using the identified variable features with `Seurat::RunPCA()`.
#' }
#' These steps are essential for reducing the dimensionality of the data and identifying principal components that explain most of the variance in the dataset.
#'
#' @importFrom Seurat NormalizeData FindVariableFeatures ScaleData RunPCA VariableFeatures
#'@export
#'
perform_pca <- function(seurat_obj) {
  seurat_obj <- Seurat::NormalizeData(seurat_obj)
  seurat_obj <- Seurat::FindVariableFeatures(seurat_obj)
  seurat_obj <- Seurat::ScaleData(seurat_obj)
  seurat_obj <- Seurat::RunPCA(seurat_obj, features = Seurat::VariableFeatures(object = seurat_obj))
  return(seurat_obj)
}

#' Visualize PCA and Elbow Plot
#'
#' This function generates an Elbow plot for PCA variance and visualizes the PCA results using a DimPlot.
#'
#' @param seurat_obj A Seurat object containing single-cell RNA-seq data.
#'
#' @return The Elbow plot, PCA plot, and heatmap of the first 5 principal components.
#'
#' @importFrom Seurat ElbowPlot DimPlot DimHeatmap
#' @importFrom ggplot2 theme_bw
#' @export
visualize_pca <- function(seurat_obj) {
  print(Seurat::ElbowPlot(seurat_obj))
  Seurat::DimPlot(seurat_obj, reduction = "pca") + ggplot2::theme_bw()
  Seurat::DimHeatmap(seurat_obj, dims = 1:5, cells = 500, balanced = TRUE)
}


#' Cluster and Perform UMAP
#'
#' This function clusters the cells and performs UMAP dimensionality reduction for visualization.
#'
#' @param seurat_obj A Seurat object containing single-cell RNA-seq data.
#' @param dims Number of dimensions to use for clustering (default is 20).
#' @param resolution Resolution parameter for clustering (default is 0.1).
#'
#' @return A Seurat object with clustering and UMAP results added.
#'
#' @importFrom Seurat FindNeighbors FindClusters RunUMAP
#' @export
cluster_and_umap <- function(seurat_obj, dims = 20, resolution = 0.1) {
  seurat_obj <- Seurat::FindNeighbors(seurat_obj, dims = 1:dims)
  seurat_obj <- Seurat::FindClusters(seurat_obj, resolution = resolution)
  seurat_obj <- Seurat::RunUMAP(seurat_obj, dims = 1:dims)
  return(seurat_obj)
}


#' Visualize UMAP
#'
#' This function visualizes the UMAP embedding of the Seurat object.
#'
#' @param seurat_obj A Seurat object containing UMAP results.
#'
#' @return UMAP plot with labels and custom themes.
#'
#' @importFrom Seurat DimPlot
#' @importFrom ggplot2 theme_bw theme element_text
#' @export
visualize_umap <- function(seurat_obj) {
  Seurat::DimPlot(seurat_obj, reduction = "umap", label = TRUE, pt.size = 0.5) +
    ggplot2::theme_bw() +
    Seurat::NoLegend() +
    ggplot2::theme(
      axis.title = ggplot2::element_text(face = "bold", size = 20),
      axis.text = ggplot2::element_text(size = 16, color = "black")
    )
}

#' Identify Cluster Markers
#'
#' This function identifies marker genes for each cluster.
#'
#' @param seurat_obj A Seurat object with clustered cells.
#' @param min_pct Minimum percentage of cells expressing the gene to be considered as a marker (default is 0.25).
#' @param logfc_thresh Minimum log fold change for a gene to be considered significant (default is 0.25).
#'
#' @return A data frame containing identified markers for each cluster.
#'
#' @importFrom Seurat FindAllMarkers
#' @export
identify_markers <- function(seurat_obj, min_pct = 0.25, logfc_thresh = 0.25) {
  Seurat::FindAllMarkers(
    seurat_obj,
    only.pos = TRUE,
    min.pct = min_pct,
    logfc.threshold = logfc_thresh
  )
}


#' Visualize CD8+ T cell Markers
#'
#' This function visualizes CD8+ T cell markers using a FeaturePlot and identifies clusters enriched for these markers.
#'
#' @param seurat_obj A Seurat object.
#' @param all_markers A data frame containing all identified markers.
#' @param cd8_markers A vector of CD8+ T cell marker genes.
#' @param reduction_type Type of dimensionality reduction to use for plotting (default is "umap").
#' @param color_scheme A vector of colors for plotting (default is c("red", "blue")).
#'
#' @return Feature plot of CD8+ T cell markers and printed list of enriched clusters.
#'
#' @importFrom Seurat FeaturePlot
#' @importFrom dplyr filter
#' @export
visualize_cd8_markers <- function(seurat_obj, all_markers, cd8_markers, reduction_type = "umap", color_scheme = c("red", "blue")) {

  cd8_clusters <- dplyr::filter(all_markers, gene %in% cd8_markers)

  print(cd8_clusters)
  Seurat::FeaturePlot(seurat_obj, features = cd8_markers, reduction = reduction_type, cols = color_scheme)
}


#' Assign CD8+ T cell Labels
#'
#' This function assigns labels to CD8+ T cells based on marker gene expression (e.g., naive, cytotoxic, exhausted).
#'
#' @param seurat_obj A Seurat object with CD8+ T cell expression data.
#' @param all_markers A data frame of all markers identified across clusters.
#' @param naive_markers A vector of naive CD8+ T cell markers.
#' @param cytotoxic_markers A vector of cytotoxic CD8+ T cell markers.
#' @param exhausted_markers A vector of exhausted CD8+ T cell markers.
#'
#' @return A list containing the updated Seurat object with CD8+ labels and a subset of CD8+ T cells.
#'
#' @importFrom Seurat Idents
#' @export
assign_cd8_labels <- function(seurat_obj, all_markers, naive_markers, cytotoxic_markers, exhausted_markers) {

  naive_clusters <- unique(all_markers$cluster[all_markers$gene %in% naive_markers])
  cytotoxic_clusters <- unique(all_markers$cluster[all_markers$gene %in% cytotoxic_markers])
  exhausted_clusters <- unique(all_markers$cluster[all_markers$gene %in% exhausted_markers])

  seurat_obj$cd8_cluster <- ifelse(
    Seurat::Idents(seurat_obj) %in% naive_clusters, "Naive CD8+ T",
    ifelse(Seurat::Idents(seurat_obj) %in% cytotoxic_clusters, "Cytotoxic CD8+ T",
           ifelse(Seurat::Idents(seurat_obj) %in% exhausted_clusters, "Exhausted CD8+ T", "Non-CD8"))
  )

  cd8_t_cells <- subset(seurat_obj, subset = cd8_cluster %in% c("Naive CD8+ T", "Cytotoxic CD8+ T", "Exhausted CD8+ T"))

  return(list(seurat_obj = seurat_obj, cd8_t_cells = cd8_t_cells))
}


#' Visualize CD8+ Subtypes on UMAP
#'
#' This function visualizes the CD8+ T cell subtypes on a UMAP plot.
#'
#' @param seurat_obj A Seurat object containing CD8+ T cells with assigned labels.
#'
#' @return UMAP plot of CD8+ subtypes with custom theme.
#'
#' @importFrom Seurat DimPlot
#' @importFrom ggplot2 theme_bw theme element_text
#' @export
visualize_cd8_subtypes <- function(seurat_obj) {
  Seurat::DimPlot(seurat_obj, reduction = "umap", group.by = "cd8_cluster", label = TRUE) +
    ggplot2::theme_bw() +
    Seurat::NoLegend() +
    ggplot2::theme(
      axis.title = ggplot2::element_text(face = "bold", size = 20),
      axis.text = ggplot2::element_text(size = 16, color = "black")
    )
}


#' Perform Pseudotime Analysis
#'
#' This function performs pseudotime analysis on CD8+ T cells using the Monocle3 package.
#'
#' @param seurat_obj A Seurat object with CD8+ T cells.
#' @param num_dim Number of dimensions for preprocessing (default is 100).
#' @param root_cluster The cluster to set as the root (default is "Naive CD8+ T").
#'
#' @return A Monocle3 CellDataSet with pseudotime trajectory learned.
#'
#' @importFrom SeuratWrappers as.cell_data_set
#' @importFrom monocle3 preprocess_cds reduce_dimension cluster_cells learn_graph order_cells
#' @export
perform_pseudotime <- function(seurat_obj, num_dim = 100, root_cluster = "Naive CD8+ T") {
  cds <- SeuratWrappers::as.cell_data_set(seurat_obj)
  cds <- monocle3::preprocess_cds(cds, num_dim = num_dim)
  cds <- monocle3::reduce_dimension(cds)
  cds <- monocle3::cluster_cells(cds)
  cds <- monocle3::learn_graph(cds)

  root_cells <- colnames(seurat_obj)[seurat_obj$cd8_cluster == root_cluster]
  cds <- monocle3::order_cells(cds, root_cells = root_cells)

  return(cds)
}


#' Visualize Pseudotime Trajectory
#'
#' This function visualizes the pseudotime trajectory for CD8+ T cells.
#'
#' @param cds A Monocle3 CellDataSet containing pseudotime data.
#' @param cd8_tcells A Seurat object containing CD8+ T cell clusters.
#' @param color_by_pseudotime A logical value to color by pseudotime (default is TRUE).
#' @param label_size Size of labels in the plot (default is 6).
#'
#' @return Pseudotime trajectory plot.
#'
#' @importFrom monocle3 plot_cells
#' @importFrom ggplot2 ggtitle
#' @export
plot_cd8_trajectory <- function(cds, cd8_tcells, color_by_pseudotime = TRUE, label_size = 6) {
  if (color_by_pseudotime) {
    x <- monocle3::plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster = TRUE, label_leaves = FALSE) +
      ggplot2::ggtitle("Trajectory colored by Pseudotime")
    print(x)
  }

  cds$cd8_cluster <- cd8_tcells$cd8_cluster

  monocle3::plot_cells(cds, color_cells_by = "cd8_cluster", label_groups_by_cluster = TRUE, group_label_size = label_size) +
    ggplot2::ggtitle("Trajectory colored by CD8+ T cell types")
}

#' Run Graph Test
#'
#' This function performs a graph test on the Monocle3 CellDataSet.
#'
#' @param cds A Monocle3 CellDataSet to run the graph test on.
#' @param neighbor_graph The neighbor graph type to use ("knn" or "principal_graph").
#' @param cores Number of CPU cores to use for parallel processing.
#'
#' @return Data frame with graph test results.
#'
#' @importFrom monocle3 graph_test
#' @export
run_graph_test <- function(cds, neighbor_graph = "knn", cores = 1) {
  result <- monocle3::graph_test(cds, neighbor_graph = neighbor_graph, cores = cores)
  return(result)
}


#' Get Top Combined Ranked Genes and Save to CSV
#'
#' This function identifies the top genes based on a combined rank of
#' q-values and Moran's I statistic from spatial transcriptomics analysis,
#' and saves the results to a `.csv` file.
#'
#' @param graph_test_res A data frame containing gene statistics,
#' including `q_value` (adjusted p-value) and `morans_I` (Moran's I statistic).
#' @param q_value_threshold Numeric. The threshold for significant q-values.
#' Default is 0.05.
#' @param top_n Integer. The number of top genes to return based on the
#' combined rank. Default is 10.
#' @param output_file Character. The name of the output `.csv` file to save
#' the results. Default is "top_genes.csv".
#'
#' @return A character vector containing the names of the top-ranked genes.
#'
#' @details
#' The function filters genes based on the `q_value_threshold`, then combines
#' the ranks of `q_value` (ascending order) and `morans_I` (descending order).
#' The combined rank is used to determine the top `top_n` genes.
#' The top genes and their ranks are saved to a `.csv` file.
#'
#' @importFrom utils write.csv
#' @export
get_top_combined_genes <- function(graph_test_res, q_value_threshold = 0.05, top_n = 10, output_file = "top_genes.csv") {
  # Subset significant genes based on q_value threshold
  significant_genes <- subset(graph_test_res, q_value < q_value_threshold)

  # Add a rank column combining ranks of q_value and Moran's I
  significant_genes$rank <- base::rank(significant_genes$q_value) + base::rank(-significant_genes$morans_I)

  assign("significant_genes", significant_genes, envir = .GlobalEnv)

  # Order by combined rank and extract the top genes
  top_combined_genes <- significant_genes[base::order(significant_genes$rank), ]
  top_combined_genes <- top_combined_genes[1:top_n, ]

  # Save the top genes to a CSV file
  utils::write.csv(top_combined_genes, file = output_file, row.names = TRUE)

  # Return the names of the top-ranked genes
  return(base::rownames(top_combined_genes))
}
