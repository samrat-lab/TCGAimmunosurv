setwd("set path where you store dowoaded data from zenodo")
####################### Load the Packages ###############################

#### First step Load the cancer list ############################
# You have a list of 33 cancer types in TCGA
cancer_list <- list(
  uterine_carcino_sarcoma = "UCS",
  ovarian_serous_cystadeno_carcinoma = "OV",
  esophageal_carcinoma = "ESCA",
  lung_squamous_cell_carcinoma = "LUSC",
  rectum_adenocarcinoma = "READ",
  head_and_neck_squamous_cell_carcinoma = "HNSC",
  pancreatic_adenocarcinoma = "PAAD", 
  colon_adenocarcinoma = "COAD",
  lung_adenocarcinoma = "LUAD",
  brain_lower_grade_glioma = "LGG",
  bladder_urothelial_carcinoma = "BLCA",
  stomach_adenocarcinoma = "STAD",
  sarcoma = "SARC",
  kidney_chromophobe = "KICH",
  liver_hepato_cellular_carcinoma = "LIHC",
  breast_invasive_carcinoma = "BRCA",
  uterine_corpus_endometrial_carcinoma = "UCEC",
  glioblastoma_multiforme = "GBM",
  adrenocortical_carcinoma = "ACC",
  skincutaneous_melanoma = "SKCM",
  prostate_adenocarcinoma = "PRAD",
  cholangio_carcinoma = "CHOL",
  lymphoid_neoplasm_diffuse_large_Bcell_lymphoma = "DLBC",
  acute_myeloid_leukemia = "LAML",
  cervical_squamous_cell_carcinoma_and_endocervical_adenocarcinoma = "CESC",
  thymoma = "THYM",
  kidney_renal_papillary_cell_carcinoma = "KIRP",
  kidney_renal_clear_cell_carcinoma = "KIRC",
  testicular_germ_cell_tumors = "TGCT",
  thyroid_carcinoma = "THCA",
  pheochromocytoma_and_paraganglioma = "PCPG",
  uveal_melanoma = "UVM",
  Mesothelioma = "MESO"
)

################################################################################################

##### Select the cancer you want to study #####################
selected_cancer <- cancer_list$breast_invasive_carcinoma

##############################################################

#### Load the data for the Bulk RNA-seq analysis #########################
load_cancer_data(selected_cancer)

#########################################################################

##### Analysis Portion #################################################

#data filtering
rna_data_filtered <- filter_rna_data(rna_data = tcga_count_data, 
                                     condition_column = "sample_type", 
                                     condition_levels = c("Primary Tumor", "Solid Tissue Normal"))

# DEGs Calculation
deg_results <- calculate_DEGs(rna_data_filtered, condition_column = "sample_type", 
                              alpha = 0.05, log2FC_threshold = 2)

# Plot the results

plot_DEG_results(
  res_df = deg_results, 
  rna_data = rna_data_filtered, 
  condition_column = "sample_type", 
  log2FC_threshold = 2,
  alpha = 0.05
)

#Read cancer data

cancer_maf <- read.maf(mutation_data)

#plot summary of the mutations 
plotmafSummary(maf = cancer_maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)

#select the gene of interest
genes_of_interest <- c("KRAS", "SMAD4")

#oncoplots
oncoplot(maf = cancer_maf, fontSize = 1, titleFontSize = 1.5, legendFontSize = 1.5, top = 10)

#Transition and Transversion plot
titv = titv(maf = cancer_maf, plot = TRUE)

#Lollipop Plot of the gene of interest
lollipopPlot(maf = cancer_maf, gene = "KRAS", AACol = "HGVSp_Short", showMutationRate = TRUE)

# somatic interactions
somaticInteractions(maf = cancer_maf, top = 10, pvalue = c(0.05, 0.1))


plot_unique_mutation_types(cancer_maf)

# Generate mutation summary
mutation_counts <- generate_mutation_summary(cancer_maf, genes_of_interest)

# Plot mutation counts
plot_mutation_counts(mutation_counts)

# Style and display the mutation table
style_mutation_table(mutation_counts) # see in viewer

# Get mutated sample IDs
mutated_samples <- get_mutated_samples(cancer_maf, genes_of_interest)

# Standardize TCGA sample IDs
standardized_ids <- standardize_tcga_ids(mutated_samples)

processed_data <- preprocess_data(tcga_count_data, clinical_data, standardized_ids)

gene_metadata_dt <- as.data.frame(rowData(tcga_count_data))

count_matrix_filtered<-processed_data$count_matrix_filtered
clinical_data_filtered<-processed_data$clinical_data_filtered

#Perform Survival analysis
survival_results <- perform_survival_analysis(count_matrix_filtered, gene_metadata_dt, clinical_data_filtered)

survival_results_permutation <- perform_survival_analysis_with_permutation(count_matrix_filtered, gene_metadata_dt, clinical_data_filtered, n_permutations = 1000)

#plot the Results
plot_results(survival_results)

#extract the genes for km plot of combined genes
significant_genes_all <- extract_genes_hr_gt1(survival_results)
#or for single km plot of multiple genes
significant_genes <- ("PRSS21")

#single_single plots
plot_km_curves(strata_data, clinical_data_filtered, significant_genes)


#### for mapping of the genes and extract upregulated and downregulated genes only ###########
selected_genes <- mapping_genes(
  deg_results = deg_results,
  gene_metadata_dt = gene_metadata_dt,
  logFC_threshold = 2,
  padj_threshold = 0.05
)

#for exporting the genes survival related
export_hr_counts_to_csv(results)



####################################### Single-cell Analysis function usage #################################################

directory = ("set path where you store Single-cell dowoaded data from zenodo")
seurat_obj <- load_seurat_object(directory, selected_cancer = selected_cancer)

cancer <- filter_seurat_object(seurat_obj, min_features = 200, max_features = 2500, max_percent_mt = 5)

cancer <- perform_pca(cancer)

visualize_pca(cancer)

cancer <- cluster_and_umap(cancer, dims = 20, resolution = 0.5)

visualize_umap(cancer)

markers<- identify_markers(cancer, min_pct = 0.25, logfc_thresh = 0.25)

naive_markers <- c("CD45RA", "CD62L", "CCR7", "TCF7", "LEF1")
cytotoxic_markers <- c("GZMB", "PRF1", "IFNG", "TNF", "KLRG1")
exhausted_markers <- c("PDCD1", "HAVCR2", "LAG3", "CTLA4", "TOX")


cd8_markers <- c(naive_markers, cytotoxic_markers, exhausted_markers)

visualize_cd8_markers(cancer, markers, cd8_markers, reduction_type = "umap", color_scheme = c("red", "blue"))

results <- assign_cd8_labels(cancer, markers, naive_markers, cytotoxic_markers, exhausted_markers)

cd8_tcells <- results$cd8_t_cells

visualize_cd8_subtypes(cd8_tcells)

cds <- perform_pseudotime(cd8_tcells, num_dim = 3, root_cluster = "Naive CD8+ T")

plot_cd8_trajectory(cds, cd8_tcells, label_size = 4)

graph_test_res <- run_graph_test(cds, neighbor_graph = "principal_graph", cores = 18)

genes <- get_top_combined_genes(graph_test_res, q_value_threshold = 0.05, top_n = 300)

######################################################################################################################
