#' Load Cancer Data for a Selected Cancer Type
#'
#' This function loads or downloads TCGA data for a specified cancer type, including RNA-Seq counts,
#' clinical data, and mutation data. If the data file already exists, it is loaded directly;
#' otherwise, it downloads and processes the data before saving it locally.
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
#' If the data file does not exist locally, this function will:
#' - Query and download the RNA-Seq count data.
#' - Query and download the clinical data.
#' - Query and download the mutation data.
#' - Prepare the data for analysis.
#' - Save all datasets in a `.RData` file for future use.
#'
#' The function uses `GDCquery`, `GDCdownload`, and `GDCprepare` from the `TCGAbiolinks` package
#' for accessing and preparing TCGA data.
#'
#' @examples
#' \dontrun{
#' # Load or download data for breast cancer (BRCA)
#' load_cancer_data("BRCA")
#' }
#'
#' @import TCGAbiolinks
#' @export
load_cancer_data <- function(selected_cancer_type) {
  rdata_filename <- paste0("tcga_", selected_cancer_type, "_data.RData")

  if (file.exists(rdata_filename)) {
    load(rdata_filename, envir = .GlobalEnv)
    cat("Data successfully loaded for analysis.\n")
  } else {
    cat("Data file not found. Initiating download and preparation...\n")

    # Define save file name
    save_filename <- rdata_filename

    # ============================================
    # SAVING TRANSCRIPTOMICS DATA (RNA-SEQ COUNTS)
    # ============================================
    query <- GDCquery(project = paste0("TCGA-", selected_cancer_type),
                      data.category = "Transcriptome Profiling",
                      experimental.strategy = "RNA-Seq",
                      workflow.type = "STAR - Counts",
                      data.type = "Gene Expression Quantification",
                      sample.type = c("Primary Tumor", "Solid Tissue Normal"),
                      access = "open")

    GDCdownload(query)  # Download the data
    tcga_count_data <- GDCprepare(query, summarizedExperiment = TRUE)  # Prepare the count data
    sample_type <- data.frame(tcga_count_data$sample_type)

    # =======================
    # SAVING CLINICAL DATA
    # =======================
    clinical_data <- GDCquery_clinic(project = paste0("TCGA-", selected_cancer_type), type = "clinical")

    # Select relevant columns in clinical data
    cols_of_interest <- which(colnames(clinical_data) %in% c("vital_status", "days_to_last_follow_up", "days_to_death", "sample_type"))
    if (length(cols_of_interest) == 0) {
      stop("Error: Specified columns not found in clinical data.")
    }

    # Create new variables for survival analysis
    clinical_data$deceased <- ifelse(clinical_data$vital_status == "Alive", FALSE, TRUE)
    clinical_data$overall_survival <- ifelse(clinical_data$vital_status == "Alive",
                                             clinical_data$days_to_last_follow_up,
                                             clinical_data$days_to_death)

    # =======================
    # SAVING MUTATION DATA
    # =======================
    query_mutation <- GDCquery(
      project = paste0("TCGA-", selected_cancer_type),
      data.category = "Simple Nucleotide Variation",
      data.type = "Masked Somatic Mutation",
      workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
    )

    GDCdownload(query_mutation)  # Download mutation data
    mutation_data <- GDCprepare(query_mutation)  # Prepare the mutation data

    # Saving all datasets into a single RData file
    save(list = c("tcga_count_data", "clinical_data", "mutation_data", "sample_type"), file = save_filename)
    message(paste0("Data for ", selected_cancer_type, " saved as ", save_filename))

    # Load the data into the environment for immediate use
    load(save_filename, envir = .GlobalEnv)
    cat("Data successfully downloaded, prepared, and loaded for analysis.\n")
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
#' @examples
#' # Example usage with a SummarizedExperiment object (rna_data)
#' filtered_rna <- filter_rna_data(
#'   rna_data = rna_data,
#'   condition_column = "condition",
#'   condition_levels = c("treated", "control"),
#'   na_threshold = 0.8
#' )
#'
#' @import SummarizedExperiment
#' @export
filter_rna_data <- function(rna_data, condition_column, condition_levels, na_threshold = 0.8) {
  # Check if condition_column exists
  if (!condition_column %in% colnames(colData(rna_data))) {
    stop(paste("Column", condition_column, "not found in colData of the provided SummarizedExperiment object."))
  }

  # Filter samples based on condition levels
  valid_samples <- colData(rna_data)[[condition_column]] %in% condition_levels
  rna_data <- rna_data[, valid_samples]

  # Update the condition column
  colData(rna_data)[[condition_column]] <- factor(colData(rna_data)[[condition_column]], levels = condition_levels)

  # Retrieve count data
  count_data <- assay(rna_data)

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
  assay(rna_data, withDimnames = FALSE) <- as.matrix(filtered_data)

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
#' @examples
#' \dontrun{
#' # Example usage:
#' library(DESeq2)
#' # Assuming `rna_counts` is a count matrix and `metadata` contains sample info
#' rna_data <- DESeqDataSetFromMatrix(countData = rna_counts,
#'                                    colData = metadata,
#'                                    design = ~ condition)
#' # Perform differential expression analysis
#' results_df <- calculate_DEGs(rna_data, condition_column = "condition", alpha = 0.05, log2FC_threshold = 1)
#' }
#'
#' @import DESeq2
#' @export
calculate_DEGs <- function(rna_data, condition_column, alpha, log2FC_threshold) {
  # Prepare DESeq2 dataset
  dds <- DESeqDataSet(rna_data, design = as.formula(paste0("~ ", condition_column)))
  dds <- DESeq(dds)

  # Extract results
  res <- results(dds, alpha = alpha,
                 contrast = c(condition_column,
                              levels(colData(rna_data)[[condition_column]])[1],
                              levels(colData(rna_data)[[condition_column]])[2]))

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
#'   sample metadata accessible via `colData`.
#' @param condition_column A character string specifying the name of the column in
#'   `colData(rna_data)` that contains condition labels (e.g., "Control" and "Treatment").
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
#' @examples
#' \dontrun{
#' # Assuming `res_df` is the result data frame and `rna_data` is a SummarizedExperiment object
#' plot_DEG_results(res_df, rna_data, condition_column = "Condition")
#' }
#'
#' @import ggplot2
#' @import SummarizedExperiment
#' @export
plot_DEG_results <- function(res_df, rna_data, condition_column, log2FC_threshold = 1, alpha = 0.05) {

  # Define font styles
  font_style <- theme(
    axis.text.x = element_text(size = 20, color = "black", face = "bold"),  # Adjusted x-axis text
    axis.text.y = element_text(size = 20, color = "black", face = "bold"),  # Adjusted y-axis text
    axis.title.x = element_text(size = 22, color = "black", face = "bold"),
    axis.title.y = element_text(size = 22, color = "black", face = "bold"),
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 18)
  )

  # Bar plot for condition distribution
  sample_counts <- as.data.frame(table(colData(rna_data)[[condition_column]]))
  colnames(sample_counts) <- c("Condition", "Frequency")  # Rename columns for better readability

  bar_plot <- ggplot(sample_counts, aes(x = Condition, y = Frequency)) +
    geom_bar(stat = "identity", fill = c("royalblue", "red"), width = 0.5) +  # Adjusted bar width
    geom_text(aes(label = Frequency), vjust = -0.3, size = 8, fontface = "bold") +  # Add numbers above bars
    theme_minimal() +
    labs(
      title = paste("Sample Distribution of", condition_column),
      x = "Condition",
      y = "Number of Samples"
    ) +
    font_style

  print(bar_plot)

  # Volcano plot
  volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = significance), alpha = 0.6) +
    scale_color_manual(
      values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "black"),
      name = "Direction"
    ) +
    guides(color = guide_legend(title = "Direction", override.aes = list(size = 3))) +
    theme_minimal() +
    labs(
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
#'   This object is typically created using the `read.maf` function from the `maftools` package.
#'
#' @return A ggplot object showing the bar plot of unique mutation types and their total number
#'   of samples across all genes in the provided `MAF` object.
#'
#' @details
#' - The function extracts the mutation type counts from the `MAF` object using the
#'   `getGeneSummary` function from the `maftools` package.
#' - Mutation types with zero counts are filtered out.
#' - The mutation counts are aggregated across all genes for each mutation type.
#' - The resulting plot displays mutation types on the x-axis and the total number of
#'   samples on the y-axis.
#'
#' @examples
#' \dontrun{
#' # Assuming `cancer_maf` is a MAF object for a selected cancer type
#' plot_unique_mutation_types(cancer_maf)
#' }
#'
#' @import maftools
#' @import ggplot2
#' @import dplyr
#' @export
plot_unique_mutation_types <- function(cancer_maf) {
  # Extract gene summary
  gene_summary <- getGeneSummary(cancer_maf)

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
  mutation_summary <- mutation_df %>%
    group_by(Mutation_Type) %>%
    summarise(NumSamples = sum(NumSamples), .groups = "drop")

  # Create the plot
  ggplot(mutation_summary, aes(x = Mutation_Type, y = NumSamples, fill = Mutation_Type)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    labs(
      title = paste("All mutations present in", "selected_cancer"),
      x = "Mutation Type",
      y = "Total Number of Samples"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
}

#' Get Genes of Interest
#'
#' This function captures a list of gene names provided as arguments, validates the input,
#' and returns the gene names as a character vector.
#'
#' @param ... A list of gene names to be captured as a character vector. Each gene should
#'   be provided as a separate argument.
#'
#' @return A character vector containing the names of the genes provided as input.
#'
#' @details
#' - The function ensures that at least one valid gene name is provided.
#' - If no valid gene names are supplied, it stops with an error message.
#' - This function can be used to capture user-specified gene names dynamically.
#'
#' @examples
#' \dontrun{
#' # Provide a list of genes
#' genes <- get_genes_of_interest("TP53", "BRCA1", "KRAS")
#' print(genes)
#'
#' # Error case: No genes provided
#' genes <- get_genes_of_interest()
#' }
#'
#' @export
get_genes_of_interest <- function(...) {
  # Capture all arguments as a character vector
  genes_input <- c(...)

  # Validate input: Check if genes_input is provided and not empty
  if (is.null(genes_input) || length(genes_input) == 0 || all(genes_input == "")) {
    stop("Error: No valid genes provided. Please provide at least one gene.")
  }

  # Return the genes as a character vector
  return(genes_input)
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
#' @examples
#' \dontrun{
#' # Load a MAF file using maftools
#' cancer_maf <- read.maf("example.maf")
#'
#' # Specify genes of interest
#' genes <- c("TP53", "BRCA1", "KRAS")
#'
#' # Generate mutation summary
#' mutation_summary <- generate_mutation_summary(cancer_maf, genes)
#' print(mutation_summary)
#' }
#'
#' @importFrom data.table data.table
#' @importFrom maftools subsetMaf getGeneSummary
#' @export
generate_mutation_summary <- function(cancer_maf, genes_of_interest) {
  all_gene_mutation_counts <- data.table()

  # Loop over each gene in the list of genes of interest
  for (gene in genes_of_interest) {
    # Check if the gene exists in the MAF data
    if (!gene %in% unique(cancer_maf@gene.summary$Hugo_Symbol)) {
      warning(paste("Gene", gene, "not found in the MAF data. Skipping."))
      next
    }

    # Filter the MAF data for the current gene's mutations
    gene_mutations <- subsetMaf(maf = cancer_maf, genes = gene)

    # Check if the gene has non-synonymous variants
    if (nrow(gene_mutations@data) == 0) {
      warning(paste("Gene", gene, "has no non-synonymous variants. Skipping."))
      next
    }

    # Get summary of mutations for the current gene
    gene_summary <- getGeneSummary(gene_mutations)

    # Extract mutation types and counts
    mutation_columns <- setdiff(colnames(gene_summary), c("Hugo_Symbol", "total", "MutatedSamples", "AlteredSamples"))
    gene_mutation_counts <- data.table(
      Gene = gene_summary$Hugo_Symbol[1],
      Mutation_Type = mutation_columns,
      NumSamples = as.numeric(gene_summary[1, ..mutation_columns])
    )

    # Append to the final table
    all_gene_mutation_counts <- rbind(all_gene_mutation_counts, gene_mutation_counts)
  }

  return(all_gene_mutation_counts)
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
#' @examples
#' \dontrun{
#' # Load a MAF file using maftools
#' cancer_maf <- read.maf("example.maf")
#'
#' # Specify genes of interest
#' genes <- c("TP53", "BRCA1", "KRAS")
#'
#' # Get mutated sample IDs
#' mutated_samples <- get_mutated_samples(cancer_maf, genes)
#' print(mutated_samples)
#' }
#'
#' @importFrom maftools subsetMaf
#' @export
get_mutated_samples <- function(cancer_maf, genes_of_interest) {
  # Subset MAF data for the genes of interest
  gene_mutations <- subsetMaf(maf = cancer_maf, genes = genes_of_interest)

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
#' @examples
#' \dontrun{
#' # Example data frame with mutation counts
#' mutation_counts <- data.frame(
#'   Gene = c("TP53", "TP53", "BRCA1", "BRCA1"),
#'   Mutation_Type = c("Missense", "Nonsense", "Missense", "Frameshift"),
#'   NumSamples = c(50, 20, 40, 10)
#' )
#'
#' # Generate mutation counts plots
#' plot_mutation_counts(mutation_counts)
#' }
#'
#' @importFrom ggplot2 ggplot geom_bar theme_minimal labs element_text unit
#' @export
plot_mutation_counts <- function(mutation_counts) {
  # Ensure the Gene column is treated as a factor to iterate over unique genes
  unique_genes <- unique(mutation_counts$Gene)

  # Iterate through each gene and create a separate plot
  for (gene in unique_genes) {
    # Filter data for the current gene
    gene_data <- subset(mutation_counts, Gene == gene)

    # Create the plot for the current gene
    p <- ggplot(gene_data, aes(x = Mutation_Type, y = NumSamples, fill = Mutation_Type)) +
      geom_bar(stat = "identity") +
      theme_minimal() +
      labs(
        title = paste("Mutation Counts for Gene:", gene),
        x = "Mutation Type",
        y = "Number of Samples"
      ) +
      theme(
        # Customize axis text and labels
        axis.text.x = element_text(color = "black", face = "bold", size = 18, angle = 45, hjust = 1),
        axis.text.y = element_text(color = "black", face = "bold", size = 18),
        axis.title.x = element_text(color = "black", face = "bold", size = 20),
        axis.title.y = element_text(color = "black", face = "bold", size = 20),

        # Customize plot title
        plot.title = element_text(color = "black", face = "bold", size = 22, hjust = 0.5),

        # Customize legend
        legend.text = element_text(color = "black", size = 14, face = "bold"),
        legend.title = element_text(color = "black", size = 16, face = "bold"),
        legend.key.size = unit(0.8, "cm")  # Adjust the size of the legend key
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
#' The function uses `kable` to create a table and then applies various styling options:
#' - `striped`, `hover`, `condensed`, and `responsive` bootstrap options for table formatting.
#' - The first column (Gene) is made bold with a fixed width of 3 cm.
#' - The second column (Mutation Type) has a width of 4 cm.
#' - The third column (NumSamples) is colored blue.
#' - A footnote is added to summarize the content of the table.
#'
#' @examples
#' \dontrun{
#' # Example mutation counts data
#' mutation_counts <- data.frame(
#'   Gene = c("TP53", "BRCA1", "EGFR"),
#'   Mutation_Type = c("Missense", "Frameshift", "Nonsense"),
#'   NumSamples = c(50, 20, 10)
#' )
#'
#' # Generate styled mutation counts table
#' style_mutation_table(mutation_counts)
#' }
#'
#' @importFrom kableExtra kbl kable_styling column_spec footnote
#' @importFrom dplyr %>%
#' @export
style_mutation_table <- function(mutation_counts) {
  mutation_counts %>%
    kbl(caption = "Mutation Counts Across Samples for Genes of Interest") %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive")) %>%
    column_spec(1, bold = TRUE, width = "3cm") %>%
    column_spec(2, width = "4cm") %>%
    column_spec(3, width = "3cm", color = "blue") %>%
    footnote(general = "This table summarizes the types of mutations and the number of samples in which they occur for the selected genes.")
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
#' @examples
#' \dontrun{
#' # Example TCGA sample IDs
#' sample_ids <- c("TCGA-01-1234", "TCGA-02-5678", "XYZ-12-3456")
#'
#' # Standardize sample IDs
#' standardized_ids <- standardize_tcga_ids(sample_ids)
#' print(standardized_ids)
#' }
#'
#' @export
standardize_tcga_ids <- function(sample_ids) {
  sample_ids <- as.character(sample_ids)  # Ensure input is character
  sapply(sample_ids, function(x) {
    if (grepl("-", x)) {  # Check if the string contains "-"
      paste(strsplit(x, "-")[[1]][1:3], collapse = "-")
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
#' @examples
#' \dontrun{
#' # Example RNA-Seq and clinical data
#' rna_data <- readRDS("rna_data.rds")
#' clinical_data <- read.csv("clinical_data.csv")
#' mutated_samples <- c("TCGA-01-1234", "TCGA-02-5678")
#'
#' # Preprocess the data
#' processed_data <- preprocess_data(rna_data, clinical_data, mutated_samples)
#'
#' # View processed count matrix and clinical data
#' processed_data$count_matrix_filtered
#' processed_data$clinical_data_filtered
#' }
#'
#' @export
preprocess_data <- function(rna_data, clinical_data, extracted_mutated_sample_ids) {

  count_matrix <- assay(rna_data, "unstranded")
  gene_metadata_dt <- as.data.frame(rowData(rna_data))

  clinical_submitter_ids <- data.frame(submitter_id = clinical_data$submitter_id)

  samples_with_mutation <- intersect(extracted_mutated_sample_ids, clinical_submitter_ids$submitter_id)

  clinical_data_filtered <- clinical_data %>%
    filter(submitter_id %in% samples_with_mutation)

  clinical_data_filtered <- as.data.table(clinical_data_filtered)

  colnames <- data.frame(colnames(count_matrix))
  extracted_colnames <- sapply(colnames$colnames.count_matrix., function(x) {
    paste(strsplit(x, "-")[[1]][1:3], collapse = "-")
  })
  colnames$trimmed <- extracted_colnames

  colnames_filtered <- colnames %>%
    filter(trimmed %in% samples_with_mutation)

  count_matrix_filtered <- as.data.table(count_matrix[, colnames(count_matrix) %in% colnames_filtered$colnames.count_matrix.])
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
#' @details
#' The function calculates the median expression of each gene across samples, stratifies the samples into `HIGH` and `LOW` based on this median, and then performs Cox regression on the association between the gene's expression and overall survival. Genes with only one strata level are excluded from the analysis.
#'
#' @examples
#' \dontrun{
#' # Example RNA-Seq and clinical data
#' rna_data <- readRDS("rna_data.rds")
#' clinical_data <- read.csv("clinical_data.csv")
#'
#' # Preprocess data to obtain filtered RNA-Seq count matrix and clinical data
#' processed_data <- preprocess_data(rna_data, clinical_data, mutated_samples)
#' count_matrix_filtered <- processed_data$count_matrix_filtered
#' clinical_data_filtered <- processed_data$clinical_data_filtered
#' gene_metadata_dt <- data.frame(gene_name = rownames(count_matrix_filtered))
#'
#' # Perform survival analysis
#' survival_results <- perform_survival_analysis(count_matrix_filtered, gene_metadata_dt, clinical_data_filtered)
#'
#' # View results
#' head(survival_results)
#' }
#'
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
  count_matrix_filtered_dt <- as.data.table(count_matrix_filtered, keep.rownames = "gene_name")

  # Melting to long format and calculate strata for each sample by gene name
  strata_data <- melt(count_matrix_filtered_dt, id.vars = "gene_name", variable.name = "case_id", value.name = "counts")[
    , strata := ifelse(counts >= median(counts, na.rm = TRUE), "HIGH", "LOW"), by = gene_name
  ]

  strata_data[, strata := factor(strata, levels = c("LOW", "HIGH"))]

  # Remove genes with only one strata level
  single_strata_genes <- strata_data[, .N, by = .(gene_name, strata)][, .N, by = gene_name][N == 1]
  strata_data <- strata_data[gene_name %in% single_strata_genes$gene_name == FALSE]

  # Merge with clinical data
  strata_data[, case_id := gsub('-01.*', '', case_id)]
  clinical_data_filtered[, submitter_id := gsub('-01.*', '', submitter_id)]

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
      coxph(Surv(overall_survival, deceased) ~ strata, data = gene_data_long)
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
#' @import ggplot2
#' @import dplyr
#'
#' @examples
#' # Assuming 'results' is a data frame containing the survival analysis results
#' plot_results(results)
#'
plot_results <- function(results) {
  library(ggplot2)
  library(dplyr)

  # Data Preparation
  significance_counts <- results %>%
    group_by(significant) %>%
    summarise(count = n()) %>%
    mutate(significant = factor(significant, levels = c("Yes", "No")))

  significant_hr_counts <- results %>%
    filter(significant == "Yes") %>%
    mutate(Category = ifelse(hazard_ratio > 1, "HR > 1", "HR < 1")) %>%
    group_by(Category) %>%
    summarise(count = n())

  top_genes <- results %>%
    arrange(desc(hazard_ratio)) %>%
    dplyr::slice(1:5)

  # Plot 1: Significant vs Non-Significant
  plot1 <- ggplot(significance_counts, aes(x = significant, y = count, fill = significant)) +
    geom_bar(stat = "identity", width = 0.5) +
    scale_fill_manual(values = c("Yes" = "#FF6666", "No" = "#BFBFBF")) +
    labs(title = "Number of Significant vs Non-Significant Genes", x = "Significance", y = "Count") +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 20, color = "black", face = "bold"),
      axis.text.y = element_text(size = 20, color = "black", face = "bold"),
      axis.title.x = element_text(size = 22, color = "black", face = "bold"),
      axis.title.y = element_text(size = 22, color = "black", face = "bold"),
      plot.title = element_text(size = 24, face = "bold", hjust = 0.5)
    )

  # Plot 2: Significant Genes by Hazard Ratio
  plot2 <- ggplot(significant_hr_counts, aes(x = Category, y = count, fill = Category)) +
    geom_bar(stat = "identity", width = 0.5) +
    geom_text(aes(label = count), vjust = 1, size = 8, color = "black", fontface = "bold", position = position_stack(vjust = 0.5)) +  # Add frequency labels
    scale_fill_manual(values = c("HR > 1" = "#66CCFF", "HR < 1" = "#FFCC66")) +
    labs(
      title = "Significant Genes by Hazard Ratio",
      x = "Hazard Ratio Category",
      y = "Count"
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 20, color = "black", face = "bold"),
      axis.text.y = element_text(size = 20, color = "black", face = "bold"),
      axis.title.x = element_text(size = 22, color = "black", face = "bold"),
      axis.title.y = element_text(size = 22, color = "black", face = "bold"),
      plot.title = element_text(size = 24, face = "bold", hjust = 0.5)
    )

  # Plot 3: Volcano Plot
  plot3 <- ggplot(results, aes(x = log(hazard_ratio), y = -log10(pvalue))) +
    geom_point(aes(color = significant)) +
    scale_color_manual(values = c("Yes" = "#FF6666", "No" = "#BFBFBF")) +
    labs(title = "Volcano Plot", x = "Log(Hazard Ratio)", y = "-Log10(P-value)") +
    theme_minimal() +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 20, face = "bold"),
      legend.text = element_text(size = 18),
      axis.text.x = element_text(size = 20, color = "black", face = "bold"),
      axis.text.y = element_text(size = 20, color = "black", face = "bold"),
      axis.title.x = element_text(size = 22, color = "black", face = "bold"),
      axis.title.y = element_text(size = 22, color = "black", face = "bold"),
      plot.title = element_text(size = 24, face = "bold", hjust = 0.5)
    )

  # Print Separate Plots
  print(plot1)
  print(plot2)
  print(plot3)

  return(list(SignificancePlot = plot1, HRPlot = plot2, VolcanoPlot = plot3))
}

#' Plot Kaplan-Meier Curves for Significant Genes
#'
#' This function generates Kaplan-Meier survival curves for each significant gene based on its association with the clinical data. The function filters the `strata_data` for the significant genes, merges it with clinical data, and fits the Kaplan-Meier model for each gene. A survival plot is then created for each gene showing survival probability over time.
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
#' @return A list of Kaplan-Meier plots (ggsurvplot objects) for each significant gene. The plots are printed individually for each gene.
#'
#' @import survival
#' @import survminer
#' @import dplyr
#'
#' @examples
#' # Assuming 'strata_data' contains gene expression data, 'clinical_data_filtered' contains clinical data,
#' # and 'significant_genes' is a vector of significant gene names
#' plot_km_curves(strata_data, clinical_data_filtered, significant_genes)
#'
plot_km_curves <- function(strata_data, clinical_data_filtered, significant_genes) {

  # Filter strata_data for significant genes
  filtered_strata_data <- strata_data[gene_name %in% significant_genes]

  for (gene in significant_genes) {
    cat("\nProcessing gene:", gene, "\n")

    # Subset data for the current gene and merge with clinical data
    gene_data_long <- merge(
      filtered_strata_data[gene_name == gene],
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
    surv_object <- Surv(gene_data_long$overall_survival, gene_data_long$deceased)
    assign("surv_object", surv_object, envir= .GlobalEnv)
    km_fit <- survfit(surv_object ~ strata, data = gene_data_long)

    # Plot the Kaplan-Meier curve
    plot <- ggsurvplot(
      km_fit,
      data = gene_data_long,
      pval = TRUE,
      risk.table = TRUE,
      ggtheme = theme_minimal(),
      title = paste("Kaplan-Meier Curve for", gene),
      xlab = "Time (days)",
      ylab = "Survival Probability",
      risk.table.height = 0.25
    )
    print(plot)
  }
}

#' Plot Kaplan-Meier Curve for Combined Significant Genes
#'
#' This function generates a Kaplan-Meier survival curve to evaluate the combined
#' effect of significant genes (e.g., genes with HR > 1) on survival probability.
#'
#' @param strata_data A data frame containing gene-wise survival data with columns:
#'                    \code{gene_name}, \code{case_id}, and survival-related metrics
#'                    like \code{strata}.
#' @param clinical_data_filtered A data frame with clinical information, including:
#'                                \code{submitter_id}, \code{overall_survival}, and
#'                                \code{deceased}.
#' @param significant_genes A data frame of significant genes, including at least
#'                          a \code{gene_name} column. This should be the output
#'                          of a gene selection process like hazard ratio filtering.
#'
#' @return A Kaplan-Meier plot showing the survival probability of patients based on
#'         the combined effect of significant genes.
#' @importFrom dplyr filter
#' @importFrom survival Surv survfit
#' @importFrom survminer ggsurvplot
#' @importFrom ggplot2 theme_minimal
#' @export
#'
#' @examples
#' # Example significant_genes
#' significant_genes <- data.frame(
#'   gene_name = c("GeneA", "GeneB", "GeneC")
#' )
#'
#' # Example strata_data
#' strata_data <- data.frame(
#'   gene_name = rep(c("GeneA", "GeneB", "GeneC"), each = 10),
#'   case_id = paste0("Sample", 1:30),
#'   strata = sample(c("High", "Low"), 30, replace = TRUE)
#' )
#'
#' # Example clinical_data_filtered
#' clinical_data_filtered <- data.frame(
#'   submitter_id = paste0("Sample", 1:30),
#'   overall_survival = sample(200:1000, 30),
#'   deceased = sample(c(TRUE, FALSE), 30, replace = TRUE)
#' )
#'
#' # Plot Kaplan-Meier Curve
#' plot_km_combined(strata_data, clinical_data_filtered, significant_genes)
plot_km_combined <- function(strata_data, clinical_data_filtered, significant_genes) {

  # Extract gene names
  gene_names <- significant_genes$gene_name

  # Filter strata_data for significant genes
  filtered_strata_data <- strata_data %>%
    filter(gene_name %in% gene_names)

  # Merge filtered strata data with clinical data
  combined_data <- merge(
    filtered_strata_data,
    clinical_data_filtered,
    by.x = "case_id",
    by.y = "submitter_id",
    all.x = TRUE
  )

  # Validate required columns
  required_columns <- c("overall_survival", "deceased", "strata")
  if (!all(required_columns %in% colnames(combined_data))) {
    stop("Missing required columns in the combined data.")
  }

  # Convert 'deceased' to numeric (TRUE = 1, FALSE = 0)
  combined_data$deceased <- as.numeric(combined_data$deceased)
  assign("combined_data", combined_data, envir = .GlobalEnv)

  # Fit the survival object and Kaplan-Meier model
  surv_object_new <- Surv(combined_data$overall_survival, combined_data$deceased)
  assign("surv_object_new", surv_object_new, envir = .GlobalEnv)
  km_fit <- survfit(surv_object_new ~ strata, data = combined_data)

  # Generate the Kaplan-Meier plot
  km_plot <- ggsurvplot(
    km_fit,
    data = combined_data,
    pval = TRUE,
    risk.table = TRUE,
    ggtheme = theme_minimal(),
    title = "Kaplan-Meier Curve for Combined Significant Genes (HR > 1)",
    xlab = "Time (days)",
    ylab = "Survival Probability",
    risk.table.height = 0.25
  )

  print(km_plot)
}



#' Map Differentially Expressed Genes (DEGs) with Gene Metadata
#'
#' This function takes in the results of differential gene expression analysis (`deg_results`), a gene metadata table (`gene_metadata_dt`), and filters for upregulated and downregulated genes based on user-defined thresholds for log fold change and adjusted p-value. It then returns a merged table with gene names and the associated differential expression data.
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
#' @param logFC_threshold A numeric value specifying the threshold for log2 fold change. Genes with `log2FoldChange > logFC_threshold` are considered upregulated, and genes with `log2FoldChange < -logFC_threshold` are considered downregulated.
#' @param padj_threshold A numeric value specifying the adjusted p-value threshold. Only genes with `padj < padj_threshold` will be selected.
#'
#' @return A data frame containing the upregulated and downregulated genes based on the given thresholds, with columns:
#' \describe{
#'   \item{ensemble_id}{The unique identifier for the gene.}
#'   \item{gene_name}{The gene name associated with the `ensemble_id`.}
#'   \item{log2FoldChange}{The log2 fold change values for the gene.}
#'   \item{padj}{The adjusted p-values for the gene.}
#' }
#' @import dplyr
#'
#' @examples
#' # Example usage:
#' # Assuming 'deg_results' contains DEG results and 'gene_metadata_dt' contains metadata
#' selected_genes <- mapping_genes(deg_results, gene_metadata_dt, logFC_threshold = 1, padj_threshold = 0.05)
#' head(selected_genes)
#'
mapping_genes <- function(deg_results, gene_metadata_dt, logFC_threshold, padj_threshold) {

  # Step 1: Merge the DEG results with gene metadata
  # We use rownames of gene_metadata_dt for merging instead of 'ensemble_id' column
  gene_metadata_dt$ensemble_id <- rownames(gene_metadata_dt)  # Create a column for ensemble_id
  deg_results$ensemble_id <- rownames(deg_results)  # Add ensemble_id as a column in deg_results

  # Merge the data frames by ensemble_id
  merged_data <- merge(deg_results, gene_metadata_dt[, c("ensemble_id", "gene_name")], by = "ensemble_id", all.x = TRUE)

  # Remove rows with NA values
  merged_data <- na.omit(merged_data)

  # Step 2: Filter for upregulated and downregulated genes based on the log fold change threshold
  # Upregulated: logFC > logFC_threshold
  # Downregulated: logFC < -logFC_threshold
  upregulated_genes <- merged_data[merged_data$log2FoldChange > logFC_threshold & merged_data$padj < padj_threshold, ]
  downregulated_genes <- merged_data[merged_data$log2FoldChange < -logFC_threshold & merged_data$padj < padj_threshold, ]

  # Combine upregulated and downregulated genes
  selected_genes <- rbind(upregulated_genes, downregulated_genes)

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
#' @export
#'
#' @examples
#' # Example input DataFrame
#' results <- data.frame(
#'   gene_name = c("GeneA", "GeneB", "GeneC", "GeneD"),
#'   significant = c("Yes", "No", "Yes", "Yes"),
#'   hazard_ratio = c(1.5, 0.8, 2.3, 0.6)
#' )
#'
#' # Extract significant genes with HR > 1 and save to a CSV file
#' significant_genes <- extract_genes_hr_gt1(results, output_file = "filtered_genes.csv")
#' print(significant_genes)
extract_genes_hr_gt1 <- function(results, output_file = "significant_genes.csv") {
  if (!all(c("significant", "hazard_ratio", "gene_name") %in% colnames(results))) {
    stop("The input DataFrame must contain 'significant', 'hazard_ratio', and 'gene_name' columns.")
  }

  # Filter for significant genes with hazard ratio > 1
  significant_genes <- results %>%
    filter(significant == "Yes" & hazard_ratio > 1) %>%
    select(gene_name, hazard_ratio)

  # Save the filtered results to a CSV file
  write.csv(significant_genes, file = output_file, row.names = FALSE)

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
#'
#' @import Seurat
#'
#' @examples
#' # Example usage:
#' # Load Seurat object for "PAAD" (pancreatic adenocarcinoma) cancer type from the specified directory
#' seurat_object <- load_seurat_object(directory = "/path/to/directory", selected_cancer = "PAAD")
#'
#' # Start single-cell analysis from this Seurat object
#' # Perform clustering, dimensionality reduction, and other single-cell analysis steps
#'
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
#' @import Seurat
#'
#' @examples
#' # Example usage:
#' # Filter Seurat object based on quality control criteria
#' filtered_seurat_obj <- filter_seurat_object(seurat_obj = seurat_object, min_features = 300, max_features = 5000, max_percent_mt = 10)
#'
#' # Continue with downstream analysis on the filtered Seurat object
#'
filter_seurat_object <- function(seurat_obj, min_features = 200, max_features = 2500, max_percent_mt = 5) {
  # Calculate mitochondrial percentage and add to metadata
  seurat_obj <- JoinLayers(seurat_obj)
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

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
#'   \item Normalizes the data using `NormalizeData()`.
#'   \item Identifies highly variable features with `FindVariableFeatures()`.
#'   \item Scales the data using `ScaleData()`.
#'   \item Runs PCA using the identified variable features with `RunPCA()`.
#' }
#' These steps are essential for reducing the dimensionality of the data and identifying principal components that explain most of the variance in the dataset.
#'
#' @import Seurat
#'
#' @examples
#' # Example usage:
#' # Perform PCA on a Seurat object
#' seurat_obj <- perform_pca(seurat_obj)
#'
#' # Visualize PCA results
#' DimPlot(seurat_obj, reduction = "pca")
#'
perform_pca <- function(seurat_obj) {
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
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
#' @import Seurat
#' @export
visualize_pca <- function(seurat_obj) {
  print(ElbowPlot(seurat_obj))
  DimPlot(seurat_obj, reduction = "pca") + theme_bw()
  DimHeatmap(seurat_obj, dims = 1:5, cells = 500, balanced = TRUE)
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
#' @import Seurat
#' @export
cluster_and_umap <- function(seurat_obj, dims = 20, resolution = 0.1) {
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:dims)
  seurat_obj <- FindClusters(seurat_obj, resolution = resolution)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:dims)
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
#' @import Seurat
#' @export
visualize_umap <- function(seurat_obj) {
  DimPlot(seurat_obj, reduction = "umap", label = TRUE, pt.size = 0.5) +
    theme_bw() +
    NoLegend() +
    theme(
      axis.title = element_text(face = "bold", size = 20),
      axis.text = element_text(size = 16, color = "black")
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
#' @import Seurat
#' @export
identify_markers <- function(seurat_obj, min_pct = 0.25, logfc_thresh = 0.25) {
  FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = min_pct, logfc.threshold = logfc_thresh)
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
#' @import Seurat dplyr
#' @export
visualize_cd8_markers <- function(seurat_obj, all_markers, cd8_markers, reduction_type = "umap", color_scheme = c("red", "blue")) {

  cd8_clusters <- all_markers %>%
    filter(gene %in% cd8_markers)

  print(cd8_clusters)
  FeaturePlot(seurat_obj, features = cd8_markers, reduction = reduction_type, cols = color_scheme)
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
#' @import Seurat dplyr
#' @export
assign_cd8_labels <- function(seurat_obj, all_markers, naive_markers, cytotoxic_markers, exhausted_markers) {

  naive_clusters <- unique(all_markers$cluster[all_markers$gene %in% naive_markers])
  cytotoxic_clusters <- unique(all_markers$cluster[all_markers$gene %in% cytotoxic_markers])
  exhausted_clusters <- unique(all_markers$cluster[all_markers$gene %in% exhausted_markers])

  seurat_obj$cd8_cluster <- ifelse(
    Idents(seurat_obj) %in% naive_clusters, "Naive CD8+ T",
    ifelse(Idents(seurat_obj) %in% cytotoxic_clusters, "Cytotoxic CD8+ T",
           ifelse(Idents(seurat_obj) %in% exhausted_clusters, "Exhausted CD8+ T", "Non-CD8"))
  )

  cd8_t_cells <- subset(seurat_obj, cd8_cluster %in% c("Naive CD8+ T", "Cytotoxic CD8+ T", "Exhausted CD8+ T"))

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
#' @import Seurat
#' @export
visualize_cd8_subtypes <- function(seurat_obj) {
  DimPlot(seurat_obj, reduction = "umap", group.by = "cd8_cluster", label = TRUE) +
    theme_bw() +
    NoLegend() +
    theme(
      axis.title = element_text(face = "bold", size = 20),
      axis.text = element_text(size = 16, color = "black")
    )
}

#' Perform Pseudotime Analysis
#'
#' This function performs pseudotime analysis on CD8+ T cells using the Monocle package.
#'
#' @param seurat_obj A Seurat object with CD8+ T cells.
#' @param num_dim Number of dimensions for preprocessing (default is 100).
#' @param root_cluster The cluster to set as the root (default is "Naive CD8+ T").
#'
#' @return A Monocle CellDataSet with pseudotime trajectory learned.
#'
#' @import Seurat monocle3
#' @export
perform_pseudotime <- function(seurat_obj, num_dim = 100, root_cluster = "Naive CD8+ T") {
  cds <- as.cell_data_set(seurat_obj)
  cds <- preprocess_cds(cds, num_dim = num_dim)
  cds <- reduce_dimension(cds)
  cds <- cluster_cells(cds)
  cds <- learn_graph(cds)

  root_cells <- colnames(seurat_obj)[seurat_obj$cd8_cluster == root_cluster]
  cds <- order_cells(cds, root_cells = root_cells)

  return(cds)
}

#' Visualize Pseudotime Trajectory
#'
#' This function visualizes the pseudotime trajectory for CD8+ T cells.
#'
#' @param cds A Monocle CellDataSet containing pseudotime data.
#' @param cd8_tcells A Seurat object containing CD8+ T cell clusters.
#' @param color_by_pseudotime A logical value to color by pseudotime (default is TRUE).
#' @param label_size Size of labels in the plot (default is 6).
#'
#' @return Pseudotime trajectory plot.
#'
#' @import monocle3 ggplot2
#' @export
plot_cd8_trajectory <- function(cds, cd8_tcells, color_by_pseudotime = TRUE, label_size = 6) {
  if (color_by_pseudotime) {
    x <- plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster = TRUE, label_leaves = FALSE) +
      ggtitle("Trajectory colored by Pseudotime")
    print(x)
  }

  cds$cd8_cluster <- cd8_tcells$cd8_cluster

  plot_cells(cds, color_cells_by = "cd8_cluster", label_groups_by_cluster = TRUE, group_label_size = label_size) +
    ggtitle("Trajectory colored by CD8+ T cell types")
}

#' Run Graph Test
#'
#' This function performs a graph test on the Monocle CellDataSet.
#'
#' @param cds A Monocle CellDataSet to run the graph test on.
#'
#' @return Data frame with graph test results.
#'
#' @import monocle3
#' @export
run_graph_test <- function(cds, neighbor_graph, cores) {
  result <- graph_test(cdscds, neighbor_graph = neighbor_graph, cores = cores)
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
#' @examples
#' # Example usage:
#' # Assuming `graph_test_res` is a data frame with columns `q_value` and `morans_I`
#' top_genes <- get_top_combined_genes(graph_test_res, q_value_threshold = 0.05, 
#'                                     top_n = 10, output_file = "top_genes.csv")
#' print(top_genes)
#'
#' @export
get_top_combined_genes <- function(graph_test_res, q_value_threshold = 0.05, top_n = 10, output_file = "top_genes.csv") {
  # Subset significant genes based on q_value threshold
  significant_genes <- subset(graph_test_res, q_value < q_value_threshold)
  
  # Add a rank column combining ranks of q_value and Moran's I
  significant_genes$rank <- rank(significant_genes$q_value) + rank(-significant_genes$morans_I)
  
  assign("significant_genes", significant_genes, envir = .GlobalEnv)
  
  # Order by combined rank and extract the top genes
  top_combined_genes <- significant_genes[order(significant_genes$rank), ]
  top_combined_genes <- top_combined_genes[1:top_n, ]
  
  # Save the top genes to a CSV file
  write.csv(top_combined_genes, file = output_file, row.names = TRUE)
  
  # Return the names of the top-ranked genes
  return(rownames(top_combined_genes))
}

