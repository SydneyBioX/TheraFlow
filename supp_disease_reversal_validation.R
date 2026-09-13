#!/usr/bin/env Rscript
# ==============================================================================
# TheraFlow Supplementary Disease Reversal Validation Script
# ==============================================================================
# Mimics TheraFlow/pipeline/disease_reversal.R for the experimental validation set
# evaluating known drug-target mechanisms across disease transcriptomics signatures
# ==============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(DTSEA)
  library(igraph)
  library(parallel)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

# ------------------------------------------------------------------------------
# 1. CLI Arguments & Path Resolution
# ------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

parse_arg <- function(flag, default_val = NULL) {
  idx <- which(args == flag)
  if (length(idx) > 0 && idx < length(args)) {
    return(args[idx + 1])
  }
  return(default_val)
}

# Resolve script directory for relative fallbacks
get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  return(getwd())
}

script_dir <- get_script_dir()

find_file <- function(candidates) {
  for (cand in candidates) {
    if (file.exists(cand)) return(cand)
  }
  return(candidates[1])
}

val_dir_rel <- "method_processing/data/disease_reversal_validation"

opt <- list(
  targets_csv  = parse_arg("--targets_csv", find_file(c(
    file.path(script_dir, val_dir_rel, "target_genes.csv"),
    file.path("TheraFlow", val_dir_rel, "target_genes.csv"),
    file.path(val_dir_rel, "target_genes.csv"),
    "target_genes.csv"
  ))),
  de_csv       = parse_arg("--de_csv", find_file(c(
    file.path(script_dir, val_dir_rel, "DE_embl_download.csv"),
    file.path("TheraFlow", val_dir_rel, "DE_embl_download.csv"),
    file.path(val_dir_rel, "DE_embl_download.csv"),
    "DE_embl_download.csv"
  ))),
  glaucoma_rds = parse_arg("--glaucoma_rds", find_file(c(
    file.path(script_dir, val_dir_rel, "glaucoma_deg.rds"),
    file.path("TheraFlow", val_dir_rel, "glaucoma_deg.rds"),
    file.path(val_dir_rel, "glaucoma_deg.rds"),
    "glaucoma_deg.rds"
  ))),
  download_dir = parse_arg("--download_dir", find_file(c(
    file.path(script_dir, val_dir_rel, "downloaded"),
    file.path("TheraFlow", val_dir_rel, "downloaded"),
    file.path(val_dir_rel, "downloaded"),
    "downloaded"
  ))),
  ppi_rds      = parse_arg("--ppi_rds", find_file(c(
    file.path(script_dir, "method_processing/data/full_human_ppi.rds"),
    "TheraFlow/method_processing/data/full_human_ppi.rds",
    "method_processing/data/full_human_ppi.rds",
    "../method_processing/data/full_human_ppi.rds"
  ))),
  output_csv   = parse_arg("--output_csv", find_file(c(
    file.path(script_dir, val_dir_rel, "disease_reversal_results.csv"),
    file.path("TheraFlow", val_dir_rel, "disease_reversal_results.csv"),
    file.path(val_dir_rel, "disease_reversal_results.csv"),
    "disease_reversal_results.csv"
  ))),
  workers      = as.numeric(parse_arg("--workers", "4")),
  nperm        = as.numeric(parse_arg("--nperm", "1000"))
)

message("==================================================")
message("TheraFlow Disease Reversal Validation Pipeline")
message("==================================================")
message("Targets CSV:        ", opt$targets_csv)
message("DE EMBL CSV:        ", opt$de_csv)
message("Glaucoma RDS:       ", opt$glaucoma_rds)
message("PPI Network:        ", opt$ppi_rds)
message("Output CSV:         ", opt$output_csv)
message("Workers (Parallel): ", opt$workers)
message("DTSEA Permutations: ", opt$nperm)
message("==================================================")

# ------------------------------------------------------------------------------
# 2. Network & Metadata Ingestion
# ------------------------------------------------------------------------------
ensure_ppi_network <- function(ppi_path) {
  if (file.exists(ppi_path)) {
    message("Loading Human PPI Network from: ", ppi_path)
    return(readRDS(ppi_path))
  }
  stop(paste("PPI network file not found at:", ppi_path))
}

full_human_ppi <- ensure_ppi_network(opt$ppi_rds)
ppi_nodes <- igraph::V(full_human_ppi)$name
message("PPI network loaded with ", length(ppi_nodes), " nodes and ", igraph::ecount(full_human_ppi), " edges.")

# ------------------------------------------------------------------------------
# 3. Target Parsing & Standardization
# ------------------------------------------------------------------------------
parse_target_symbols <- function(target_str) {
  if (is.na(target_str) || target_str == "") return(character(0))
  
  # Remove extraneous descriptions while preserving gene symbols
  cleaned <- target_str
  # Split by / , or whitespace
  tokens <- unlist(strsplit(cleaned, "[/,\\s()]+"))
  # Filter to typical uppercase HGNC gene symbol tokens
  tokens <- tokens[grepl("^[A-Z0-9_-]{2,}$", tokens)]
  # Remove known non-gene descriptor keywords
  tokens <- tokens[!tokens %in% c("NMDA", "HER2", "PARP", "AMPK", "INDIRECT", "RECEPTOR", "AGGREGATES")]
  
  # Special explicit fallbacks for known validation rows if needed
  if (grepl("AMPK", target_str, ignore.case = TRUE)) {
    tokens <- unique(c(tokens, "PRKAA1", "PRKAA2"))
  }
  if (grepl("APP", target_str, ignore.case = TRUE)) {
    tokens <- unique(c(tokens, "APP"))
  }
  if (grepl("GRIN2B", target_str, ignore.case = TRUE)) {
    tokens <- unique(c(tokens, "GRIN2B"))
  }
  if (grepl("ERBB2", target_str, ignore.case = TRUE)) {
    tokens <- unique(c(tokens, "ERBB2"))
  }
  if (grepl("ROCK", target_str, ignore.case = TRUE)) {
    tokens <- unique(c(tokens, "ROCK1", "ROCK2"))
  }
  if (grepl("CDK4", target_str, ignore.case = TRUE)) {
    tokens <- unique(c(tokens, "CDK4", "CDK6"))
  }
  if (grepl("ACHE", target_str, ignore.case = TRUE) && grepl("BCHE", target_str, ignore.case = TRUE)) {
    tokens <- unique(c(tokens, "ACHE", "BCHE"))
  }
  
  return(unique(tokens[tokens != ""]))
}

message("\nLoading and parsing target genes from: ", opt$targets_csv)
targets_raw <- read_csv(opt$targets_csv, show_col_types = FALSE)

# Build expanded drug-target dataframe
targets_expanded_list <- list()
for (i in seq_len(nrow(targets_raw))) {
  row <- targets_raw[i, ]
  symbols <- parse_target_symbols(row$target)
  # Filter symbols present in PPI
  valid_symbols <- intersect(symbols, ppi_nodes)
  if (length(valid_symbols) == 0) {
    # If none found directly, retain first token as fallback
    valid_symbols <- intersect(row$target, ppi_nodes)
  }
  for (sym in valid_symbols) {
    targets_expanded_list[[length(targets_expanded_list) + 1]] <- data.frame(
      disease = row$disease,
      drug = row$drug,
      target_gene = sym,
      action = row$action,
      modality = row$modality,
      stringsAsFactors = FALSE
    )
  }
}

targets_df <- bind_rows(targets_expanded_list) %>% distinct()
message("Parsed ", nrow(targets_df), " valid drug-target mappings across ", n_distinct(targets_df$drug), " drugs.")

# ------------------------------------------------------------------------------
# 4. Disease Signature Extraction
# ------------------------------------------------------------------------------
get_disease_signatures <- function(disease_name, opt, ppi_nodes) {
  # 1. Check Glaucoma RDS
  if (disease_name == "glaucoma" && file.exists(opt$glaucoma_rds)) {
    g_deg <- readRDS(opt$glaucoma_rds)
    if (is.list(g_deg) && all(c("up", "down") %in% names(g_deg))) {
      up_genes <- if (is.null(names(g_deg$up))) as.character(g_deg$up) else names(g_deg$up)
      down_genes <- if (is.null(names(g_deg$down))) as.character(g_deg$down) else names(g_deg$down)
      return(list(
        up = intersect(unique(up_genes), ppi_nodes),
        down = intersect(unique(down_genes), ppi_nodes)
      ))
    }
  }
  
  # 2. Check individual files in downloaded directory
  dl_files_map <- list(
    "breast_cancer"     = "differentialResults.tsv",
    "alzheimer_disease" = "differentialResults_ad.tsv",
    "type_2_diabetes"   = "differentialResults_dm.tsv",
    "ovarian_cancer"    = "differentialResults-ovarian.tsv"
  )
  
  deg_data <- NULL
  if (disease_name %in% names(dl_files_map) && file.exists(file.path(opt$download_dir, dl_files_map[[disease_name]]))) {
    file_path <- file.path(opt$download_dir, dl_files_map[[disease_name]])
    deg_data <- read_tsv(file_path, show_col_types = FALSE)
  } else if (file.exists(opt$de_csv)) {
    # Fallback to consolidated DE_embl_download.csv
    de_all <- read_csv(opt$de_csv, show_col_types = FALSE)
    # Match comparison patterns
    comp_patterns <- list(
      "breast_cancer"     = "breast",
      "alzheimer_disease" = "Alzheimer",
      "type_2_diabetes"   = "diabetes",
      "ovarian_cancer"    = "ovarian|HE4|SKOV3|TIGAR|CB-5083"
    )
    if (disease_name %in% names(comp_patterns)) {
      pat <- comp_patterns[[disease_name]]
      deg_data <- de_all %>% filter(grepl(pat, Comparison, ignore.case = TRUE))
    }
  }
  
  if (is.null(deg_data) || nrow(deg_data) == 0) {
    warning(paste("No differential expression data found for disease:", disease_name))
    return(list(up = character(0), down = character(0)))
  }
  
  # Map Ensembl Gene IDs to HGNC Gene Symbols
  if ("Gene" %in% colnames(deg_data)) {
    ensembl_ids <- as.character(deg_data$Gene)
    symbols <- suppressMessages(
      mapIds(org.Hs.eg.db, keys = ensembl_ids, keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first")
    )
    deg_data$hgnc_symbol <- unname(symbols)
  } else if ("hgnc_symbol" %in% colnames(deg_data)) {
    deg_data$hgnc_symbol <- as.character(deg_data$hgnc_symbol)
  }
  
  # Filter valid symbols in PPI network
  logfc_col <- grep("log.*fold.*change", colnames(deg_data), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(logfc_col)) logfc_col <- "log_2 fold change"
  
  valid_degs <- deg_data %>%
    filter(!is.na(hgnc_symbol) & hgnc_symbol %in% ppi_nodes)
  
  up_genes <- valid_degs %>% filter(.data[[logfc_col]] > 0) %>% pull(hgnc_symbol) %>% unique()
  down_genes <- valid_degs %>% filter(.data[[logfc_col]] < 0) %>% pull(hgnc_symbol) %>% unique()
  
  return(list(
    up = up_genes,
    down = down_genes
  ))
}

# ------------------------------------------------------------------------------
# 5. DTSEA Disease Reversal Execution
# ------------------------------------------------------------------------------
unique_diseases <- unique(targets_df$disease)
num_cores <- min(length(unique_diseases), as.numeric(opt$workers))
message("\nProcessing DTSEA across ", length(unique_diseases), " diseases using ", num_cores, " parallel cores...")

process_single_disease <- function(dis) {
  message("  [+] Running disease: ", dis)
  sig <- get_disease_signatures(dis, opt, ppi_nodes)
  
  up_valid <- intersect(sig$up, ppi_nodes)
  down_valid <- intersect(sig$down, ppi_nodes)
  
  dis_targets <- targets_df %>% filter(disease == dis)
  if (nrow(dis_targets) == 0) return(NULL)
  
  # Form DTSEA drug target input: drug_id, gene_target
  ds_drugs <- dis_targets %>%
    select(drug_id = drug, gene_target = target_gene) %>%
    distinct()
  
  if (nrow(ds_drugs) == 0) return(NULL)
  
  # Up-regulated DTSEA
  up_result <- if (length(up_valid) > 0) {
    tryCatch({
      DTSEA(network = full_human_ppi, disease = up_valid, drugs = ds_drugs, nPermSimple = opt$nperm, verbose = FALSE)
    }, error = function(e) {
      warning(paste("DTSEA Up Error for", dis, ":", e$message))
      data.frame(drug_id = character(), NES = numeric(), padj = numeric())
    })
  } else {
    data.frame(drug_id = character(), NES = numeric(), padj = numeric())
  }
  
  # Down-regulated DTSEA
  down_result <- if (length(down_valid) > 0) {
    tryCatch({
      DTSEA(network = full_human_ppi, disease = down_valid, drugs = ds_drugs, nPermSimple = opt$nperm, verbose = FALSE)
    }, error = function(e) {
      warning(paste("DTSEA Down Error for", dis, ":", e$message))
      data.frame(drug_id = character(), NES = numeric(), padj = numeric())
    })
  } else {
    data.frame(drug_id = character(), NES = numeric(), padj = numeric())
  }
  
  all_drugs <- unique(dis_targets$drug)
  res_df <- data.frame(drug_id = all_drugs, stringsAsFactors = FALSE)
  
  if (nrow(up_result) > 0 && "NES" %in% colnames(up_result)) {
    up_sub <- as.data.frame(up_result)[, c("drug_id", "NES", "padj")]
    colnames(up_sub) <- c("drug_id", "nes_up", "up_padj")
    res_df <- merge(res_df, up_sub, by = "drug_id", all.x = TRUE)
  } else {
    res_df$nes_up <- 0
    res_df$up_padj <- 1.0
  }
  
  if (nrow(down_result) > 0 && "NES" %in% colnames(down_result)) {
    down_sub <- as.data.frame(down_result)[, c("drug_id", "NES", "padj")]
    colnames(down_sub) <- c("drug_id", "nes_down", "down_padj")
    res_df <- merge(res_df, down_sub, by = "drug_id", all.x = TRUE)
  } else {
    res_df$nes_down <- 0
    res_df$down_padj <- 1.0
  }
  
  res_df$nes_up[is.na(res_df$nes_up)] <- 0
  res_df$nes_down[is.na(res_df$nes_down)] <- 0
  res_df$up_padj[is.na(res_df$up_padj)] <- 1.0
  res_df$down_padj[is.na(res_df$down_padj)] <- 1.0
  
  # Disease Reversal Score: DRS = NES(down) - NES(up)
  res_df$reversal_score <- res_df$nes_down - res_df$nes_up
  res_df$disease <- dis
  
  # Join modality & metadata
  meta_sub <- targets_raw %>%
    filter(disease == dis) %>%
    select(drug_id = drug, modality, action, target) %>%
    distinct()
  
  res_df <- res_df %>% left_join(meta_sub, by = "drug_id")
  return(res_df)
}

# Run disease reversal calculation
disease_results <- mclapply(unique_diseases, process_single_disease, mc.cores = num_cores)
disease_results <- disease_results[!sapply(disease_results, is.null)]

if (length(disease_results) > 0) {
  final_df <- bind_rows(disease_results)
  
  # Select and order columns matching validation specification
  final_out <- final_df %>%
    select(drug_id, nes_up, nes_down, modality, reversal_score, disease)
  
  # Write output
  dir.create(dirname(opt$output_csv), showWarnings = FALSE, recursive = TRUE)
  fwrite(final_out, opt$output_csv)
  message("\n>>> Successfully wrote ", nrow(final_out), " validation reversal results to: ", opt$output_csv)
  print(head(final_out, 15))
} else {
  warning("No disease reversal results generated.")
}
