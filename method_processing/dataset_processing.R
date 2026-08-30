#!/usr/bin/env Rscript
# Unified Preprocessing R Script for Disease Datasets
# ===================================================

# -------------------------------------------------------------------------
# 1. Setup paths and directories
# -------------------------------------------------------------------------
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
if (length(script.name) == 0) {
  script_dir <- getwd()
} else {
  script_dir <- dirname(normalizePath(script.name))
}

DATA_DIR <- file.path(script_dir, "data")
WORKSPACE_ROOT <- dirname(dirname(script_dir))

# Ensure output directories exist
rest_methods_out_dir <- file.path(DATA_DIR, "test_datasets", "rest_methods")
gexmolgen_out_dir <- file.path(DATA_DIR, "test_datasets", "gexmolgen")
paccmann_out_dir  <- file.path(DATA_DIR, "test_datasets", "paccmann_rl")
triomphe_out_dir   <- file.path(DATA_DIR, "test_datasets", "triomphe")
# smilesgen_model_dir <- file.path(WORKSPACE_ROOT, "SmilesGEN", "datasets", "test_protein")

dir.create(rest_methods_out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(gexmolgen_out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(paccmann_out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(triomphe_out_dir, showWarnings = FALSE, recursive = TRUE)
# dir.create(smilesgen_model_dir, showWarnings = FALSE, recursive = TRUE)

# -------------------------------------------------------------------------
# 2. Check and load libraries
# -------------------------------------------------------------------------
required_packages <- c("dplyr", "tidyr", "tibble", "readr", "edgeR", 
                       "limma", "SummarizedExperiment", "reticulate", "GEOquery", "org.Hs.eg.db")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Installing missing package:", pkg))
    if (pkg %in% c("edgeR", "limma", "SummarizedExperiment", "GEOquery", "org.Hs.eg.db")) {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager", repos = "https://cloud.r-project.org")
      BiocManager::install(pkg, update = FALSE)
    } else {
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
  }
  library(pkg, character.only = TRUE)
}
select <- dplyr::select

# -------------------------------------------------------------------------
# 3. Global Configurations & Reference Templates
# -------------------------------------------------------------------------
paccmann_pkl_path <- file.path(DATA_DIR, "2128_genes.pkl") # can you copy this into the TheraFlow/method_processing/data directory 
paccmann_genes <- unlist(py_load_object(paccmann_pkl_path))
paccmann_template <- tibble(gene_symbol = paccmann_genes)

# Load template from rest_methods
lincs_m_path <- file.path(rest_methods_out_dir, "lincs1000_1.csv")
if (!file.exists(lincs_m_path)) {
  stop("lincs1000_1.csv template not found in rest_methods folder. Please run LINCS preprocessing first.")
}
MCF7_AKT1_mean <- read_csv(lincs_m_path, show_col_types = FALSE)
target_order <- colnames(MCF7_AKT1_mean)[-1]
template <- tibble(gene_id = target_order)

# Load GExMolGen templates
gex_ctl_path <- file.path(gexmolgen_out_dir, "lincs1000_1_ctl.csv")
if (!file.exists(gex_ctl_path)) {
  stop("lincs1000_1_ctl.csv template not found in gexmolgen folder. Please run LINCS preprocessing first.")
}
server_test_ctl_AKT2 <- read_csv(gex_ctl_path, show_col_types = FALSE)
target_symbols <- colnames(server_test_ctl_AKT2)[-1]

# Triomphe reference columns
reference_df <- read_csv(file.path(DATA_DIR, "L1000_ideal_cp_example.csv"), show_col_types = FALSE) # can you copy this into the TheraFlow/method_processing/data directory 
target_cols <- colnames(reference_df)

# Landmark gene mapping info (Consolidated LINCS folder)
gene_info_path <- file.path(DATA_DIR, "LINCS", "GSE70138", "GSE70138_Broad_LINCS_gene_info_2017-03-06.txt.gz")
gene_info <- read_tsv(gene_info_path, show_col_types = FALSE)
gene_info_lm <- gene_info[gene_info$pr_is_lm == 1, ]

# Function to format matrix into "server_test" style (Level 3 pairs)
format_lincs_level3 <- function(sce_obj, samples, target_cols) {
  sce_subset <- sce_obj[, colnames(sce_obj) %in% samples]
  mat <- edgeR::cpm(assay(sce_subset), log = TRUE, prior.count = 1)
  df_wide <- as.data.frame(t(mat))
  missing <- setdiff(target_cols, colnames(df_wide))
  if (length(missing) > 0) {
    if (nrow(df_wide) == 0) {
      df_wide <- as.data.frame(matrix(0, nrow = length(samples), ncol = length(target_cols)))
      colnames(df_wide) <- target_cols
      rownames(df_wide) <- samples
    } else {
      df_wide[missing] <- 0
    }
  }
  df_final <- df_wide %>%
    mutate(`...1` = row_number() - 1) %>%
    dplyr::select(`...1`, all_of(target_cols))
  return(df_final)
}

# Helper to format Triomphe datasets
format_triomphe <- function(sig_final, name) {
  n_samples <- nrow(sig_final)
  metadata_dummy <- data.frame(
    `Unnamed: 0` = 0:(n_samples - 1),
    sig_id = sig_final$sig_id,
    InChIKey = "-666",
    drug_id = "-",
    pert_iname = paste0(name, "_Signature"),
    pert_type = "trt_cp",
    cell_line = "Primary_Glia",
    dose = "1 uL",
    time = "96 h",
    check.names = FALSE
  )
  expression_data <- sig_final %>% rename_with(~ paste0("hsa:", .), -sig_id)
  lincs_formatted_df <- metadata_dummy %>% left_join(expression_data, by = "sig_id")
  missing_cols <- setdiff(target_cols, colnames(lincs_formatted_df))
  if (length(missing_cols) > 0) {
    zero_mat <- matrix(0, nrow = nrow(lincs_formatted_df), ncol = length(missing_cols))
    colnames(zero_mat) <- missing_cols
    lincs_formatted_df <- cbind(lincs_formatted_df, as.data.frame(zero_mat))
  }
  lincs_formatted_df <- lincs_formatted_df %>% dplyr::select(all_of(target_cols))
  return(lincs_formatted_df)
}

# Helper to format paccmann datasets
format_paccmann <- function(level5_sig, paccmann_genes, paccmann_template) {
  paccmann_sig_df <- level5_sig %>%
    enframe(name = "gene_symbol", value = "expression") %>%
    group_by(gene_symbol) %>%
    summarise(expression = expression[which.max(abs(expression))], .groups = "drop")
  paccmann_final <- paccmann_template %>%
    left_join(paccmann_sig_df, by = "gene_symbol") %>%
    mutate(expression = replace_na(expression, 0)) %>%
    pivot_wider(names_from = gene_symbol, values_from = expression)
  return(paccmann_final)
}

# Helper to format rest_methods datasets
format_rest_methods <- function(level5_symbols, name) {
  my_sig_df <- level5_symbols %>%
    enframe(name = "gene_id", value = "expression") %>%
    mutate(sig_id = name) %>%
    group_by(gene_id) %>%
    summarise(expression = expression[which.max(abs(expression))], .groups = "drop")
  my_sig_complete <- template %>%
    left_join(my_sig_df, by = "gene_id") %>%
    mutate(expression = replace_na(expression, 0), sig_id = name)
  my_sig_final <- my_sig_complete %>%
    pivot_wider(names_from = gene_id, values_from = expression) %>%
    dplyr::select(sig_id, all_of(target_order))
  return(my_sig_final)
}

# -------------------------------------------------------------------------
# 4. Processing dilated cardiomyopathy (GSE55296)
# -------------------------------------------------------------------------
message("Processing GSE55296 (dilated cardiomyopathy)...")
bulk_counts <- read_tsv(file.path(DATA_DIR, "bulk_transcriptomics_data", "GSE55296", "GSE55296_count_data.txt.gz"), show_col_types = FALSE)
count_matrix <- as.matrix(bulk_counts[, 3:38])
count_matrix[is.na(count_matrix)] <- 0
cpm_matrix <- edgeR::cpm(count_matrix, log = TRUE)
gene_names <- bulk_counts[, 2]
colnames(gene_names) <- "Gene"

sample_info <- read_tsv(file.path(DATA_DIR, "bulk_transcriptomics_data", "GSE55296", "GSE55296_processed_data_readme.txt"), show_col_types = FALSE)
sample_info$CONDITION <- sub(" .*", "", sample_info$`sample title`)
sample_info$CONDITION <- case_when(sample_info$CONDITION == "Ischemic" ~ "Cardiomyopathy",
                                   sample_info$CONDITION == "Dilated" ~ "Cardiomyopathy",
                                   .default = "Control")
sample_info$CONDITION <- factor(sample_info$CONDITION, levels = c("Control", "Cardiomyopathy"))

design <- model.matrix(~ CONDITION, data = sample_info)
cpm_matrix <- cpm_matrix[, sample_info$`column headers (in the count_data.txt)`]
level3_lincs <- cbind(gene_names, cpm_matrix)

fit <- lmFit(level3_lincs, design)
fit <- eBayes(fit)

level5_signature <- fit$coefficients[, "CONDITIONCardiomyopathy"]
names(level5_signature) <- gene_names$Gene
level5_symbols_cardio <- level5_signature[names(level5_signature) %in% gene_info_lm$pr_gene_symbol]
map_idx_cardio <- match(names(level5_symbols_cardio), gene_info_lm$pr_gene_symbol)
names(level5_symbols_cardio) <- gene_info_lm$pr_gene_id[map_idx_cardio]

# Save paccmann
paccmann_final_cardio <- format_paccmann(level5_signature, paccmann_genes, paccmann_template)
write_csv(paccmann_final_cardio, file.path(paccmann_out_dir, "cardiomyopathy.csv"))

# Save rest_methods
my_sig_final_cardio <- format_rest_methods(level5_symbols_cardio, "cardiomyopathy")
## write_csv(my_sig_final_cardio, file.path(smilesgen_model_dir, "MCF7_cardiomyopathy.csv")) # can be ignored
write_csv(my_sig_final_cardio, file.path(rest_methods_out_dir, "cardiomyopathy.csv"))

# Save gexmolgen
counts_mat_cardio <- as.matrix(bulk_counts[, 3:ncol(bulk_counts)])
rownames(counts_mat_cardio) <- bulk_counts$...2
counts_mat_cardio <- counts_mat_cardio[, colnames(counts_mat_cardio) %in% sample_info$`column headers (in the count_data.txt)`]
counts_mat_cardio[is.na(counts_mat_cardio)] <- 0
sce_bulk_cardio <- SummarizedExperiment(assays = list(counts = counts_mat_cardio), colData = sample_info)
ctrl_samples_cardio <- colData(sce_bulk_cardio)$`column headers (in the count_data.txt)`[colData(sce_bulk_cardio)$CONDITION == "Control"]
disease_samples_cardio <- colData(sce_bulk_cardio)$`column headers (in the count_data.txt)`[colData(sce_bulk_cardio)$CONDITION != "Control"]
server_test_ctl_cardio <- format_lincs_level3(sce_bulk_cardio, ctrl_samples_cardio, target_symbols)
server_test_poag_cardio <- format_lincs_level3(sce_bulk_cardio, disease_samples_cardio, target_symbols)
write_csv(server_test_ctl_cardio, file.path(gexmolgen_out_dir, "cardiomyopathy_ctl.csv"))
write_csv(server_test_poag_cardio, file.path(gexmolgen_out_dir, "cardiomyopathy_disease.csv"))

# Save Triomphe
triomphe_cardio <- format_triomphe(my_sig_final_cardio, "cardiomyopathy")
write_csv(triomphe_cardio, file.path(triomphe_out_dir, "cardiomyopathy_pert.csv"))

# -------------------------------------------------------------------------
# 5. Processing Parkinson's Disease (GSE110716)
# -------------------------------------------------------------------------
message("Processing GSE110716 (Parkinson's Disease)...")
untar_dir <- file.path(DATA_DIR, "bulk_transcriptomics_data", "GSE110716", "extracted_raw")
dir.create(untar_dir, showWarnings = FALSE, recursive = TRUE)
untar(file.path(DATA_DIR, "bulk_transcriptomics_data", "GSE110716", "GSE110716_RAW.tar"), exdir = untar_dir)

pd_metadata <- read_csv(file.path(DATA_DIR, "bulk_transcriptomics_data", "GSE110716_clinical_metadata.csv"), show_col_types = FALSE)
pd_metadata$CONDITION <- ifelse(grepl("control", pd_metadata$`patient status/id:ch1`, ignore.case = TRUE), "Control", "PD")
pd_metadata$CONDITION <- factor(pd_metadata$CONDITION, levels = c("Control", "PD"))

pd_files <- list.files(untar_dir, pattern = "\\.txt\\.gz$", full.names = TRUE)
pd_counts_list <- list()
for (f in pd_files) {
  gsm_id <- sub("_.*", "", basename(f))
  df_f <- read_tsv(f, show_col_types = FALSE)
  df_f <- df_f %>% 
    dplyr::select(Gene_Symbol, Count) %>%
    group_by(Gene_Symbol) %>%
    summarise(Count = sum(Count, na.rm = TRUE), .groups = "drop") %>%
    filter(!is.na(Gene_Symbol) & Gene_Symbol != "")
  colnames(df_f) <- c("Gene", gsm_id)
  pd_counts_list[[gsm_id]] <- df_f
}

pd_counts <- Reduce(function(x, y) full_join(x, y, by = "Gene"), pd_counts_list)
pd_counts[is.na(pd_counts)] <- 0
gsm_cols <- pd_metadata$geo_accession
counts_mat_pd <- as.matrix(pd_counts[, gsm_cols])
rownames(counts_mat_pd) <- pd_counts$Gene

cpm_mat_pd <- edgeR::cpm(counts_mat_pd, log = TRUE)
design_pd <- model.matrix(~ CONDITION, data = pd_metadata)

fit_pd <- lmFit(cpm_mat_pd, design_pd)
fit_pd <- eBayes(fit_pd)

level5_sig_pd <- fit_pd$coefficients[, "CONDITIONPD"]
names(level5_sig_pd) <- rownames(cpm_mat_pd)
level5_symbols_pd <- level5_sig_pd[names(level5_sig_pd) %in% gene_info_lm$pr_gene_symbol]
map_idx_pd <- match(names(level5_symbols_pd), gene_info_lm$pr_gene_symbol)
names(level5_symbols_pd) <- gene_info_lm$pr_gene_id[map_idx_pd]

# Save paccmann
paccmann_final_pd <- format_paccmann(level5_sig_pd, paccmann_genes, paccmann_template)
write_csv(paccmann_final_pd, file.path(paccmann_out_dir, "parkinsons.csv"))

# Save rest_methods
my_sig_final_pd <- format_rest_methods(level5_symbols_pd, "parkinsons")
# # write_csv(my_sig_final_pd, file.path(smilesgen_model_dir, "MCF7_parkinsons.csv"))
write_csv(my_sig_final_pd, file.path(rest_methods_out_dir, "parkinsons.csv"))

# Save gexmolgen
sce_pd <- SummarizedExperiment(assays = list(counts = counts_mat_pd), colData = pd_metadata)
ctrl_samples_pd <- colData(sce_pd)$geo_accession[colData(sce_pd)$CONDITION == "Control"]
disease_samples_pd <- colData(sce_pd)$geo_accession[colData(sce_pd)$CONDITION == "PD"]
server_test_ctl_pd <- format_lincs_level3(sce_pd, ctrl_samples_pd, target_symbols)
server_test_poag_pd <- format_lincs_level3(sce_pd, disease_samples_pd, target_symbols)
write_csv(server_test_ctl_pd, file.path(gexmolgen_out_dir, "parkinsons_ctl.csv"))
write_csv(server_test_poag_pd, file.path(gexmolgen_out_dir, "parkinsons_disease.csv"))

# Save Triomphe
triomphe_pd <- format_triomphe(my_sig_final_pd, "parkinsons")
write_csv(triomphe_pd, file.path(triomphe_out_dir, "parkinsons_pert.csv"))

# -------------------------------------------------------------------------
# 6. Processing HCC (GSE143233)
# -------------------------------------------------------------------------
message("Processing GSE143233 (HCC)...")
hcc_data <- read_tsv(file.path(DATA_DIR, "bulk_transcriptomics_data", "GSE143233", "GSE143233_patients_normalHCC_and_sorafenib_resistantHCC_genes_fpkm_tracking.txt.gz"), show_col_types = FALSE)
level5_sig_hcc <- log2(hcc_data$`sorafenib-resistant_FPKM` + 1) - log2(hcc_data$normal_FPKM + 1)
names(level5_sig_hcc) <- hcc_data$gene_short_name
level5_symbols_hcc <- level5_sig_hcc[names(level5_sig_hcc) %in% gene_info_lm$pr_gene_symbol]
map_idx_hcc <- match(names(level5_symbols_hcc), gene_info_lm$pr_gene_symbol)
names(level5_symbols_hcc) <- gene_info_lm$pr_gene_id[map_idx_hcc]

# Save paccmann
paccmann_final_hcc <- format_paccmann(level5_sig_hcc, paccmann_genes, paccmann_template)
write_csv(paccmann_final_hcc, file.path(paccmann_out_dir, "hcc.csv"))

# Save rest_methods
my_sig_final_hcc <- format_rest_methods(level5_symbols_hcc, "hcc")
# write_csv(my_sig_final_hcc, file.path(smilesgen_model_dir, "MCF7_hcc.csv"))
write_csv(my_sig_final_hcc, file.path(rest_methods_out_dir, "hcc.csv"))

# Save gexmolgen (using FPKM mock data control/perturbed)
hcc_ctl_mock <- data.frame(
  `...1` = 0,
  t(log2(hcc_data$normal_FPKM[match(target_symbols, hcc_data$gene_short_name)] + 1))
)
colnames(hcc_ctl_mock) <- c("...1", target_symbols)
hcc_ctl_mock[is.na(hcc_ctl_mock)] <- 0

hcc_disease_mock <- data.frame(
  `...1` = 0,
  t(log2(hcc_data$`sorafenib-resistant_FPKM`[match(target_symbols, hcc_data$gene_short_name)] + 1))
)
colnames(hcc_disease_mock) <- c("...1", target_symbols)
hcc_disease_mock[is.na(hcc_disease_mock)] <- 0

write_csv(hcc_ctl_mock, file.path(gexmolgen_out_dir, "hcc_ctl.csv"))
write_csv(hcc_disease_mock, file.path(gexmolgen_out_dir, "hcc_disease.csv"))

# Save Triomphe
triomphe_hcc <- format_triomphe(my_sig_final_hcc, "hcc")
write_csv(triomphe_hcc, file.path(triomphe_out_dir, "hcc_pert.csv"))

# -------------------------------------------------------------------------
# 7. Processing Glaucoma Bulk (GSE236302)
# -------------------------------------------------------------------------
message("Processing GSE236302 (Glaucoma Bulk)...")
raw_counts_path <- file.path(DATA_DIR, "bulk_transcriptomics_data", "GSE236302", "GSE236302_raw_counts.csv.gz")
metadata_path <- file.path(DATA_DIR, "bulk_transcriptomics_data", "GSE236302", "GSE236302_sample_metadata.csv.gz")

glaucoma_counts <- read_delim(raw_counts_path, delim = ";", show_col_types = FALSE)
glaucoma_metadata <- read_delim(metadata_path, delim = ";", show_col_types = FALSE)
colnames(glaucoma_counts)[1] <- "Gene"
counts_mat_glaucoma <- as.matrix(glaucoma_counts[, -1])
rownames(counts_mat_glaucoma) <- glaucoma_counts$Gene

# Map Ensembl IDs to Gene Symbols
ensembl_ids_all <- rownames(counts_mat_glaucoma)
mapping_all <- AnnotationDbi::select(org.Hs.eg.db, keys = ensembl_ids_all, columns = "SYMBOL", keytype = "ENSEMBL")
mapping_all <- mapping_all %>% filter(!is.na(SYMBOL)) %>% distinct(SYMBOL, .keep_all = TRUE)

counts_mat_glaucoma <- counts_mat_glaucoma[mapping_all$ENSEMBL, ]
rownames(counts_mat_glaucoma) <- mapping_all$SYMBOL

glaucoma_metadata$CONDITION <- ifelse(glaucoma_metadata$condition == "CTL", "Control", "Glaucoma")
glaucoma_metadata$CONDITION <- factor(glaucoma_metadata$CONDITION, levels = c("Control", "Glaucoma"))
glaucoma_metadata$geo_accession <- glaucoma_metadata$sample_id

sample_ids <- intersect(colnames(counts_mat_glaucoma), glaucoma_metadata$sample_id)
counts_mat_glaucoma <- counts_mat_glaucoma[, sample_ids]
glaucoma_metadata <- glaucoma_metadata[match(sample_ids, glaucoma_metadata$sample_id), ]

cpm_mat_glaucoma <- edgeR::cpm(counts_mat_glaucoma, log = TRUE, prior.count = 1)
design_glaucoma <- model.matrix(~ CONDITION, data = glaucoma_metadata)

fit_glaucoma <- lmFit(cpm_mat_glaucoma, design_glaucoma)
fit_glaucoma <- eBayes(fit_glaucoma)

level5_sig_glaucoma <- fit_glaucoma$coefficients[, "CONDITIONGlaucoma"]
names(level5_sig_glaucoma) <- rownames(cpm_mat_glaucoma)
level5_symbols_glaucoma <- level5_sig_glaucoma[names(level5_sig_glaucoma) %in% gene_info_lm$pr_gene_symbol]
map_idx_glaucoma <- match(names(level5_symbols_glaucoma), gene_info_lm$pr_gene_symbol)
names(level5_symbols_glaucoma) <- gene_info_lm$pr_gene_id[map_idx_glaucoma]

# Save paccmann
paccmann_final_glaucoma <- format_paccmann(level5_sig_glaucoma, paccmann_genes, paccmann_template)
write_csv(paccmann_final_glaucoma, file.path(paccmann_out_dir, "glaucoma_bulk.csv"))

# Save rest_methods
my_sig_final_glaucoma <- format_rest_methods(level5_symbols_glaucoma, "glaucoma_bulk")
# write_csv(my_sig_final_glaucoma, file.path(smilesgen_model_dir, "MCF7_glaucoma_bulk.csv"))
write_csv(my_sig_final_glaucoma, file.path(rest_methods_out_dir, "glaucoma_bulk.csv"))

# Save gexmolgen
sce_glaucoma <- SummarizedExperiment(assays = list(counts = counts_mat_glaucoma), colData = glaucoma_metadata)
ctrl_samples_glaucoma <- colData(sce_glaucoma)$geo_accession[colData(sce_glaucoma)$CONDITION == "Control"]
disease_samples_glaucoma <- colData(sce_glaucoma)$geo_accession[colData(sce_glaucoma)$CONDITION == "Glaucoma"]
server_test_ctl_glaucoma <- format_lincs_level3(sce_glaucoma, ctrl_samples_glaucoma, target_symbols)
server_test_poag_glaucoma <- format_lincs_level3(sce_glaucoma, disease_samples_glaucoma, target_symbols)
write_csv(server_test_ctl_glaucoma, file.path(gexmolgen_out_dir, "glaucoma_bulk_ctl.csv"))
write_csv(server_test_poag_glaucoma, file.path(gexmolgen_out_dir, "glaucoma_bulk_disease.csv"))

# Save Triomphe
triomphe_glaucoma <- format_triomphe(my_sig_final_glaucoma, "glaucoma_bulk")
write_csv(triomphe_glaucoma, file.path(triomphe_out_dir, "glaucoma_bulk_pert.csv"))

# -------------------------------------------------------------------------
# 8. Processing Glaucoma Single Cell (Maciej Daniszewski 2022)
# -------------------------------------------------------------------------
message("Processing Maciej Daniszewski 2022 (Glaucoma Single Cell)...")
sc_rds_path <- file.path(DATA_DIR, "single_cell", "Maciej_Daniszewski_2022.rds")

if (file.exists(sc_rds_path)) {
  sc_obj <- readRDS(sc_rds_path)
  
  donor_factor <- factor(colData(sc_obj)$DONOR)
  donor_matrix <- model.matrix(~ 0 + donor_factor)
  colnames(donor_matrix) <- levels(donor_factor)
  
  # Pseudobulk by summing counts per donor
  pb_counts <- assay(sc_obj, "counts") %*% donor_matrix
  
  donor_meta <- as.data.frame(colData(sc_obj)) %>%
    dplyr::select(DONOR, CONDITION) %>%
    distinct(DONOR, .keep_all = TRUE)
  donor_meta <- donor_meta[match(colnames(pb_counts), donor_meta$DONOR), ]
  
  donor_meta$CONDITION <- ifelse(donor_meta$CONDITION == "control", "Control", "Glaucoma")
  donor_meta$CONDITION <- factor(donor_meta$CONDITION, levels = c("Control", "Glaucoma"))
  donor_meta$geo_accession <- donor_meta$DONOR
  rownames(donor_meta) <- donor_meta$DONOR
  
  cpm_mat_sc <- edgeR::cpm(pb_counts, log = TRUE, prior.count = 1)
  design_sc <- model.matrix(~ CONDITION, data = donor_meta)
  
  fit_sc <- lmFit(cpm_mat_sc, design_sc)
  fit_sc <- eBayes(fit_sc)
  
  level5_sig_sc <- fit_sc$coefficients[, "CONDITIONGlaucoma"]
  names(level5_sig_sc) <- rownames(cpm_mat_sc)
  
  level5_symbols_sc <- level5_sig_sc[names(level5_sig_sc) %in% gene_info_lm$pr_gene_symbol]
  map_idx_sc <- match(names(level5_symbols_sc), gene_info_lm$pr_gene_symbol)
  names(level5_symbols_sc) <- gene_info_lm$pr_gene_id[map_idx_sc]
  
  # Save paccmann
  paccmann_final_sc <- format_paccmann(level5_sig_sc, paccmann_genes, paccmann_template)
  write_csv(paccmann_final_sc, file.path(paccmann_out_dir, "glaucoma.csv"))
  
  # Save rest_methods
  my_sig_final_sc <- format_rest_methods(level5_symbols_sc, "glaucoma")
  # write_csv(my_sig_final_sc, file.path(smilesgen_model_dir, "MCF7_glaucoma.csv"))
  write_csv(my_sig_final_sc, file.path(rest_methods_out_dir, "glaucoma.csv"))
  
  # Save gexmolgen
  sce_sc <- SummarizedExperiment(assays = list(counts = pb_counts), colData = donor_meta)
  ctrl_samples_sc <- colData(sce_sc)$geo_accession[colData(sce_sc)$CONDITION == "Control"]
  disease_samples_sc <- colData(sce_sc)$geo_accession[colData(sce_sc)$CONDITION == "Glaucoma"]
  server_test_ctl_sc <- format_lincs_level3(sce_sc, ctrl_samples_sc, target_symbols)
  server_test_poag_sc <- format_lincs_level3(sce_sc, disease_samples_sc, target_symbols)
  write_csv(server_test_ctl_sc, file.path(gexmolgen_out_dir, "glacoma_ctl.csv"))
  write_csv(server_test_poag_sc, file.path(gexmolgen_out_dir, "glacoma_disease.csv"))
  
  # Save Triomphe
  triomphe_sc <- format_triomphe(my_sig_final_sc, "glaucoma")
  write_csv(triomphe_sc, file.path(triomphe_out_dir, "glaucoma_pert.csv"))
} else {
  message("Warning: Maciej_Daniszewski_2022.rds not found. Skipping single cell glaucoma processing.")
}

message("R Disease Preprocessing Completed Successfully!")
