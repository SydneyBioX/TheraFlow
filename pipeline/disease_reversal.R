#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(DTSEA)
  library(igraph)
  library(parallel)
})

args <- commandArgs(trailingOnly = TRUE)

parse_arg <- function(flag, default_val = NULL) {
  idx <- which(args == flag)
  if (length(idx) > 0 && idx < length(args)) {
    return(args[idx + 1])
  }
  return(default_val)
}

opt <- list(
  targets_csv       = parse_arg("--targets_csv", "TheraFlow/outputs/chembl_targets.csv"),
  potency_csv       = parse_arg("--potency_csv", "TheraFlow/outputs/potency_scores.csv"),
  output_csv        = parse_arg("--output_csv", "TheraFlow/outputs/reversal_res.csv"),
  ppi_rds           = parse_arg("--ppi_rds", "TheraFlow/method_processing/data/full_human_ppi.rds"),
  lincs_meta        = parse_arg("--lincs_meta", "TheraFlow/method_processing/data/LINCS/lincs_gene_meta.csv"),
  test_datasets_dir = parse_arg("--test_datasets_dir", "TheraFlow/method_processing/data/test_datasets/rest_methods"),
  workers           = as.numeric(parse_arg("--workers", "7"))
)

ensure_ppi_network <- function(ppi_path) {
  if (file.exists(ppi_path)) {
    message("Loading existing PPI network from: ", ppi_path)
    return(readRDS(ppi_path))
  }
  stop(paste("PPI network not found at:", ppi_path))
}

process_lincs_signatures <- function(path, id2sym) {
  dt <- read.csv(path, check.names = FALSE)
  if (ncol(dt) < 2) return(list(up = character(), down = character()))
  
  g_ids <- colnames(dt)[2:ncol(dt)]
  z_vals <- as.numeric(dt[1, 2:ncol(dt)])
  g_syms <- unname(id2sym[as.character(g_ids)])
  
  valid <- !is.na(g_syms) & g_syms != ""
  df_sig <- data.frame(gene = g_syms[valid], z = z_vals[valid], stringsAsFactors = FALSE)
  
  up_genes <- unique(df_sig$gene[df_sig$z > 0])
  down_genes <- unique(df_sig$gene[df_sig$z < 0])
  
  return(list(up = up_genes, down = down_genes))
}

get_disease_signature <- function(test_dataset, test_datasets_dir, id2sym) {
  target_file <- paste0(test_dataset, ".csv")
  possible_paths <- c(
    file.path(test_datasets_dir, target_file),
    file.path(test_datasets_dir, "rest_methods", target_file),
    file.path(test_datasets_dir, "transgem", target_file),
    file.path(test_datasets_dir, "paccmann_rl", target_file),
    file.path(test_datasets_dir, "smilesgen", target_file),
    file.path("TheraFlow/method_processing/data/test_datasets/rest_methods", target_file)
  )
  
  matched_file <- NULL
  for (p in possible_paths) {
    if (file.exists(p)) {
      matched_file <- p
      break
    }
  }
  if (is.null(matched_file)) {
    all_files <- list.files("TheraFlow/method_processing/data/test_datasets", recursive = TRUE, full.names = TRUE, pattern = "\\.csv$")
    base_names <- tolower(basename(all_files))
    idx <- which(base_names == tolower(target_file))
    if (length(idx) > 0) matched_file <- all_files[idx[1]]
  }
  if (is.null(matched_file)) {
    stop(paste("Disease signature file not found for dataset:", test_dataset))
  }
  return(process_lincs_signatures(matched_file, id2sym))
}

full_human_ppi <- ensure_ppi_network(opt$ppi_rds)

message("Loading LINCS gene metadata from: ", opt$lincs_meta)
meta <- read.csv(opt$lincs_meta, check.names = FALSE)
id2sym <- setNames(as.character(meta[['pr_gene_symbol']]), as.character(meta[['pr_gene_id']]))

message("Loading predicted targets: ", opt$targets_csv)
targets_df <- read_csv(opt$targets_csv, show_col_types = FALSE)

props_candidates <- c(
  "TheraFlow/outputs/molecular_properties_metrics.csv",
  "TheraFlow/method_processing/data/processed_data.csv"
)
for (pf in props_candidates) {
  if (file.exists(pf)) {
    p_df <- read_csv(pf, show_col_types = FALSE) %>% 
      dplyr::select(any_of(c("SMILES", "method", "test_dataset"))) %>% 
      distinct()
    targets_df <- targets_df %>% left_join(p_df, by = "SMILES")
    break
  }
}

if (!("method_dataset" %in% colnames(targets_df))) {
  if ("method" %in% colnames(targets_df) && "test_dataset" %in% colnames(targets_df)) {
    targets_df$method_dataset <- paste0(targets_df$method, "_", targets_df$test_dataset)
  } else if ("test_dataset" %in% colnames(targets_df)) {
    targets_df$method_dataset <- targets_df$test_dataset
  } else {
    targets_df$method_dataset <- "unknown"
  }
}

df_clean <- targets_df %>%
  filter(!is.na(Predicted_Targets) & Predicted_Targets != "") %>%
  tidyr::separate_longer_delim(Predicted_Targets, delim = ", ") %>%
  dplyr::rename(hgnc_symbol = Predicted_Targets) %>%
  dplyr::mutate(hgnc_symbol = trimws(hgnc_symbol)) %>%
  dplyr::filter(hgnc_symbol %in% igraph::V(full_human_ppi)$name) %>%
  dplyr::distinct(SMILES, method, test_dataset, method_dataset, hgnc_symbol)

if (nrow(df_clean) == 0) {
  message("Warning: No valid HGNC targets found in PPI network. Writing empty reversal table.")
  dir.create(dirname(opt$output_csv), showWarnings = FALSE, recursive = TRUE)
  empty_df <- data.frame(method_dataset=character(), SMILES=character(), reversal_score=numeric(), reversal_score_scaled=numeric(), confidence=character(), significant=logical())
  fwrite(empty_df, opt$output_csv)
  quit(save = "no", status = 0)
}

unique_datasets <- unique(df_clean$test_dataset)
num_cores <- min(length(unique_datasets), as.numeric(opt$workers))
message(paste("Processing DTSEA in parallel across", length(unique_datasets), "unique disease datasets using", num_cores, "cores..."))

process_single_dataset <- function(ds) {
  if (is.na(ds) || is.null(ds) || ds == "") return(NULL)
  
  sig <- tryCatch({
    get_disease_signature(ds, opt$test_datasets_dir, id2sym)
  }, error = function(e) {
    NULL
  })
  
  if (is.null(sig)) return(NULL)
  
  ds_drugs <- df_clean %>%
    filter(test_dataset == ds) %>%
    select(SMILES, hgnc_symbol) %>%
    distinct()
  colnames(ds_drugs) <- c("drug_id", "gene_target")
  
  if (nrow(ds_drugs) == 0) return(NULL)
  
  up_valid <- intersect(sig$up, igraph::V(full_human_ppi)$name)
  down_valid <- intersect(sig$down, igraph::V(full_human_ppi)$name)
  
  up_result <- if (length(up_valid) > 0) {
    tryCatch({
      DTSEA(network = full_human_ppi, disease = up_valid, drugs = ds_drugs, nPermSimple = 1000, verbose = FALSE)
    }, error = function(e) data.frame(drug_id = character(), NES = numeric(), padj = numeric()))
  } else data.frame(drug_id = character(), NES = numeric(), padj = numeric())
  
  down_result <- if (length(down_valid) > 0) {
    tryCatch({
      DTSEA(network = full_human_ppi, disease = down_valid, drugs = ds_drugs, nPermSimple = 1000, verbose = FALSE)
    }, error = function(e) data.frame(drug_id = character(), NES = numeric(), padj = numeric()))
  } else data.frame(drug_id = character(), NES = numeric(), padj = numeric())
  
  all_d_ids <- unique(ds_drugs$drug_id)
  res_df <- data.frame(SMILES = all_d_ids, stringsAsFactors = FALSE)
  
  if (nrow(up_result) > 0 && "NES" %in% colnames(up_result)) {
    up_sub <- as.data.frame(up_result)[, c("drug_id", "NES", "padj")]
    colnames(up_sub) <- c("SMILES", "nes_up", "up_padj")
    res_df <- merge(res_df, up_sub, by = "SMILES", all.x = TRUE)
  } else {
    res_df$nes_up <- 0
    res_df$up_padj <- 1.0
  }
  
  if (nrow(down_result) > 0 && "NES" %in% colnames(down_result)) {
    down_sub <- as.data.frame(down_result)[, c("drug_id", "NES", "padj")]
    colnames(down_sub) <- c("SMILES", "nes_down", "down_padj")
    res_df <- merge(res_df, down_sub, by = "SMILES", all.x = TRUE)
  } else {
    res_df$nes_down <- 0
    res_df$down_padj <- 1.0
  }
  
  res_df$nes_up[is.na(res_df$nes_up)] <- 0
  res_df$nes_down[is.na(res_df$nes_down)] <- 0
  res_df$up_padj[is.na(res_df$up_padj)] <- 1.0
  res_df$down_padj[is.na(res_df$down_padj)] <- 1.0
  res_df$reversal_score <- res_df$nes_down - res_df$nes_up
  res_df$test_dataset <- ds
  
  return(res_df)
}

all_dataset_results <- mclapply(unique_datasets, process_single_dataset, mc.cores = num_cores)
all_dataset_results <- all_dataset_results[!sapply(all_dataset_results, is.null)]

if (length(all_dataset_results) > 0) {
  reversal_merged <- bind_rows(all_dataset_results)
  targets_info <- targets_df %>% 
    select(any_of(c("SMILES", "Predicted_Targets", "ChEMBL_Neighbors"))) %>% 
    distinct()
  
  mapping_df <- df_clean %>% 
    select(SMILES, method, test_dataset, method_dataset) %>% 
    distinct() %>% 
    left_join(targets_info, by = "SMILES")
    
  final_df <- mapping_df %>% 
    inner_join(reversal_merged, by = c("SMILES", "test_dataset"))
} else {
  final_df <- data.frame(method_dataset = character(), method = character(), test_dataset = character(), SMILES = character(), Predicted_Targets = character(), reversal_score = numeric(), nes_up = numeric(), nes_down = numeric(), up_padj = numeric(), down_padj = numeric())
}

if (file.exists(opt$potency_csv)) {
  potency <- read_csv(opt$potency_csv, show_col_types = FALSE)
  if ("query_smiles" %in% colnames(potency)) {
    potency <- potency %>% dplyr::rename(SMILES = query_smiles)
  }
  
  final_df <- final_df %>% left_join(potency, by = "SMILES")
  
  action_col <- if ("action" %in% colnames(final_df)) final_df$action else rep("Unknown", nrow(final_df))
  ic50_col <- if ("IC50" %in% colnames(final_df)) as.numeric(final_df$IC50) else rep(NA, nrow(final_df))
  
  final_df$reversal_score_scaled <- case_when(
    action_col %in% c("INHIBITOR", "ANTAGONIST") ~ final_df$reversal_score,
    action_col %in% c("AGONIST") ~ -final_df$reversal_score,
    !is.na(ic50_col) & ic50_col < 10000 ~ final_df$reversal_score,
    TRUE ~ final_df$reversal_score
  )
  final_df$confidence <- case_when(
    action_col %in% c("INHIBITOR", "AGONIST", "ANTAGONIST") ~ "high",
    !is.na(ic50_col) & ic50_col < 10000 ~ "high",
    TRUE ~ "low"
  )
  final_df$significant <- ifelse(
    final_df$up_padj <= 0.05 | final_df$down_padj <= 0.05, TRUE, FALSE
  )
} else {
  final_df$reversal_score_scaled <- final_df$reversal_score
  final_df$confidence <- "low"
  final_df$significant <- FALSE
}

dir.create(dirname(opt$output_csv), showWarnings = FALSE, recursive = TRUE)
fwrite(final_df, opt$output_csv)
message("\n>>> Successfully wrote ", nrow(final_df), " disease reversal results to: ", opt$output_csv)
