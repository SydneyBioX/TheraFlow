# R script for step 5: DTSEA disease reversal calculations
args <- commandArgs(trailingOnly = TRUE)
opt <- list(
  targets_csv = "TheraFlow/outputs/chembl_targets.csv",
  potency_csv = "TheraFlow/outputs/potency_scores.csv",
  output_reversal_csv = "TheraFlow/outputs/reversal_res.csv",
  ppi_rds = "TheraFlow/method_processing/data/full_human_ppi.rds",
  lincs_meta = "TheraFlow/method_processing/data/LINCS/lincs_gene_meta.csv",
  test_datasets_dir = "TheraFlow/method_processing/data/test_datasets/rest_methods",
  workers = 2
)

i <- 1
while (i <= length(args)) {
  arg <- args[i]
  if (arg == "--targets_csv" && i < length(args)) {
    opt$targets_csv <- args[i + 1]
    i <- i + 2
  } else if (arg == "--potency_csv" && i < length(args)) {
    opt$potency_csv <- args[i + 1]
    i <- i + 2
  } else if (arg == "--output_csv" && i < length(args)) {
    opt$output_reversal_csv <- args[i + 1]
    i <- i + 2
  } else if (arg == "--ppi_rds" && i < length(args)) {
    opt$ppi_rds <- args[i + 1]
    i <- i + 2
  } else if (arg == "--lincs_meta" && i < length(args)) {
    opt$lincs_meta <- args[i + 1]
    i <- i + 2
  } else if (arg == "--test_datasets_dir" && i < length(args)) {
    opt$test_datasets_dir <- args[i + 1]
    i <- i + 2
  } else if (arg == "--workers" && i < length(args)) {
    opt$workers <- as.integer(args[i + 1])
    i <- i + 2
  } else {
    i <- i + 1
  }
}

suppressPackageStartupMessages({
  library(furrr)
  library(tidyr)
  library(dplyr)
  library(data.table)
  library(tidyverse)
  library(readr)
  library(DTSEA)
  library(igraph)
})

# -------------------------------------------------------------------------
# 1. Helper function: ensure PPI network exists, build from STRING DB if missing
# -------------------------------------------------------------------------
ensure_ppi_network <- function(ppi_path = "TheraFlow/method_processing/data/full_human_ppi.rds", min_score = 700) {
  if (file.exists(ppi_path) && file.info(ppi_path)$size > 1000) {
    message("Loading existing PPI network from: ", ppi_path)
    return(readRDS(ppi_path))
  }
  
  message("PPI network not found at: ", ppi_path)
  message("Downloading and constructing human PPI network from STRING DB v12.0 (combined_score >= ", min_score, ")...")
  
  dir.create(dirname(ppi_path), showWarnings = FALSE, recursive = TRUE)
  temp_dir <- tempdir()
  links_gz <- file.path(temp_dir, "9606.protein.links.v12.0.txt.gz")
  info_gz  <- file.path(temp_dir, "9606.protein.info.v12.0.txt.gz")
  
  links_url <- "https://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz"
  info_url  <- "https://stringdb-downloads.org/download/protein.info.v12.0/9606.protein.info.v12.0.txt.gz"
  
  if (!file.exists(links_gz)) {
    message("  Fetching STRING links from ", links_url, " ...")
    download.file(links_url, links_gz, mode = "wb", quiet = TRUE)
  }
  if (!file.exists(info_gz)) {
    message("  Fetching STRING protein info from ", info_url, " ...")
    download.file(info_url, info_gz, mode = "wb", quiet = TRUE)
  }
  
  message("  Parsing and filtering STRING interactions...")
  raw_string <- fread(links_gz)
  mapping <- fread(info_gz)[, c("#string_protein_id", "preferred_name")]
  
  ppi_edges <- raw_string[combined_score >= min_score] %>%
    left_join(mapping, by = c("protein1" = "#string_protein_id")) %>%
    rename(from = preferred_name) %>%
    left_join(mapping, by = c("protein2" = "#string_protein_id")) %>%
    rename(to = preferred_name) %>%
    dplyr::select(from, to, combined_score) %>%
    filter(!is.na(from) & !is.na(to))
    
  message("  Building and simplifying igraph network...")
  full_human_ppi <- graph_from_data_frame(d = ppi_edges, directed = FALSE)
  full_human_ppi <- simplify(full_human_ppi, remove.multiple = TRUE, remove.loops = TRUE)
  
  saveRDS(full_human_ppi, ppi_path)
  message("  Successfully generated and saved PPI network to: ", ppi_path)
  
  if (file.exists(links_gz)) unlink(links_gz)
  if (file.exists(info_gz)) unlink(info_gz)
  
  return(full_human_ppi)
}

# -------------------------------------------------------------------------
# 2. Helper function: process LINCS signature files
# -------------------------------------------------------------------------
process_lincs_signatures <- function(path, lincs_meta) {
  data <- fread(path)
  colnames(data) <- as.character(data[1, ])
  
  # Take row 2 (z-scores) and columns 2 to end (skip sig_id column)
  data <- data[2, 2:ncol(data)] %>% t() %>% as.data.frame()
  colnames(data) <- "z_score"
  data$pr_gene_id <- rownames(data)
  
  data <- data %>%
    mutate(
      pr_gene_id = as.character(pr_gene_id),
      z_score = as.numeric(as.character(z_score))
    ) %>%
    inner_join(
      lincs_meta %>% mutate(pr_gene_id = as.character(pr_gene_id)), 
      by = "pr_gene_id"
    ) %>%
    dplyr::select(pr_gene_symbol, z_score) %>%
    filter(!is.na(pr_gene_symbol))
    
  signatures <- list(
    up = data %>% filter(z_score > 0) %>% deframe(),
    down = data %>% filter(z_score < 0) %>% deframe()
  )
  
  return(signatures)
}

# -------------------------------------------------------------------------
# 3. Helper function: load disease signature for the target diseases
# -------------------------------------------------------------------------
get_disease_signature <- function(test_dataset, test_datasets_dir, lincs_meta) {
  target_file <- paste0(test_dataset, ".csv")
  
  possible_paths <- c(
    file.path(test_datasets_dir, target_file),
    file.path(test_datasets_dir, "rest_methods", target_file),
    file.path(test_datasets_dir, "transgem", target_file),
    file.path(test_datasets_dir, "paccmann_rl", target_file),
    file.path(test_datasets_dir, "smilesgen", target_file),
    file.path("TheraFlow/method_processing/data/test_datasets/rest_methods", target_file),
    file.path("TheraFlow/method_processing/data/test_datasets/transgem", target_file)
  )
  
  matched_file <- NULL
  for (p in possible_paths) {
    if (file.exists(p)) {
      matched_file <- p
      break
    }
  }
  
  if (is.null(matched_file)) {
    all_files <- list.files(test_datasets_dir, recursive = TRUE, full.names = TRUE, pattern = "\\.csv$")
    base_names <- tolower(basename(all_files))
    idx <- which(base_names == tolower(target_file))
    if (length(idx) > 0) {
      matched_file <- all_files[idx[1]]
    }
  }
  
  if (is.null(matched_file)) {
    stop(paste("Disease signature file not found for dataset:", test_dataset))
  }
  
  message("Loading disease signature from: ", matched_file)
  return(process_lincs_signatures(matched_file, lincs_meta))
}

# -------------------------------------------------------------------------
# 4. Main Data Processing Workflow
# -------------------------------------------------------------------------
full_human_ppi <- ensure_ppi_network(opt$ppi_rds)

message("Loading LINCS gene metadata from: ", opt$lincs_meta)
if (!file.exists(opt$lincs_meta)) {
  stop(paste("LINCS metadata not found at:", opt$lincs_meta))
}
lincs_meta <- fread(opt$lincs_meta)

message("Loading predicted targets: ", opt$targets_csv)
metrics <- read_csv(opt$targets_csv, show_col_types = FALSE)

if (!"method_dataset" %in% colnames(metrics)) {
  if ("method" %in% colnames(metrics) && "test_dataset" %in% colnames(metrics)) {
    metrics$method_dataset <- paste0(metrics$method, "_", metrics$test_dataset)
  } else if ("test_dataset" %in% colnames(metrics)) {
    metrics$method_dataset <- metrics$test_dataset
  } else {
    metrics$method_dataset <- "unknown"
  }
}

df_clean <- metrics %>%
  dplyr::select(SMILES, Predicted_Targets, method_dataset) %>%
  filter(!is.na(Predicted_Targets) & Predicted_Targets != "None" & Predicted_Targets != "") %>%
  separate_rows(Predicted_Targets, sep = ",\\s*") %>%
  dplyr::rename(hgnc_symbol = Predicted_Targets) %>%
  distinct(SMILES, method_dataset, hgnc_symbol)

if (nrow(df_clean) == 0) {
  message("Warning: No valid HGNC targets found in Predicted_Targets column. Writing empty reversal table.")
  dir.create(dirname(opt$output_reversal_csv), showWarnings = FALSE, recursive = TRUE)
  empty_df <- data.frame(method_dataset=character(), SMILES=character(), reversal_score=numeric(), reversal_score_scaled=numeric(), confidence=character(), significant=logical())
  fwrite(empty_df, opt$output_reversal_csv)
  quit(save = "no", status = 0)
}

# -------------------------------------------------------------------------
# 5. DTSEA Disease Reversal Execution
# -------------------------------------------------------------------------
unique_groups <- unique(df_clean$method_dataset)
message(paste("Processing DTSEA for", length(unique_groups), "unique method_dataset groups..."))

dtsea_function <- function(drugs_long, m_ds, up, down) {
  message(paste("Processing group:", m_ds))
  drugs <- drugs_long %>%
    filter(method_dataset == m_ds) %>%
    dplyr::select(SMILES, hgnc_symbol)
  colnames(drugs) <- c("drug_id", "gene_target")
  
  up_result <- tryCatch({
    DTSEA(
      network = full_human_ppi,
      disease = names(up),
      drugs = drugs,
      verbose = FALSE
    )
  }, error = function(e) {
    data.frame(drug_id = character(), NES = numeric(), padj = numeric())
  })
  
  down_result <- tryCatch({
    DTSEA(
      network = full_human_ppi,
      disease = names(down),
      drugs = drugs,
      verbose = FALSE
    )
  }, error = function(e) {
    data.frame(drug_id = character(), NES = numeric(), padj = numeric())
  })
  
  if (nrow(up_result) == 0 || nrow(down_result) == 0) {
    return(data.frame(smiles_id = character(), reversal_score = numeric(), up_padj = numeric(), down_padj = numeric()))
  }
  
  up_df <- data.frame(smiles_id = up_result$drug_id, nes_up = up_result$NES, up_padj = up_result$padj)
  down_df <- data.frame(smiles_id = down_result$drug_id, nes_down = down_result$NES, down_padj = down_result$padj)
  
  combined_df <- merge(up_df, down_df, by = "smiles_id")
  combined_df$reversal_score <- combined_df$nes_up - combined_df$nes_down
  
  return(combined_df[, c("smiles_id", "reversal_score", "up_padj", "down_padj")])
}

all_results <- list()
for (m_ds in unique_groups) {
  parts <- strsplit(m_ds, "_")[[1]]
  ds_name <- paste(parts[-1], collapse = "_")
  if (ds_name == "") ds_name <- m_ds
  
  sig <- tryCatch({
    get_disease_signature(ds_name, opt$test_datasets_dir, lincs_meta)
  }, error = function(e) {
    tryCatch({
      get_disease_signature(parts[length(parts)], opt$test_datasets_dir, lincs_meta)
    }, error = function(e2) NULL)
  })
  
  if (!is.null(sig)) {
    res <- dtsea_function(df_clean, m_ds, sig$up, sig$down)
    if (nrow(res) > 0) {
      all_results[[m_ds]] <- res
    }
  }
}

if (length(all_results) > 0) {
  reversal_df <- bind_rows(all_results, .id = "method_dataset")
} else {
  reversal_df <- data.frame(method_dataset = character(), smiles_id = character(), reversal_score = numeric(), up_padj = numeric(), down_padj = numeric())
}

# -------------------------------------------------------------------------
# 6. Potency & Mechanism-Scaled Reversal Score
# -------------------------------------------------------------------------
if (file.exists(opt$potency_csv)) {
  potency <- read_csv(opt$potency_csv, show_col_types = FALSE)
  if ("query_smiles" %in% colnames(potency)) {
    potency <- potency %>% dplyr::rename(smiles_id = query_smiles)
  } else if ("SMILES" %in% colnames(potency)) {
    potency <- potency %>% dplyr::rename(smiles_id = SMILES)
  }
  
  final_reversal_df <- reversal_df %>% left_join(potency, by = "smiles_id")
  
  if (!"action" %in% colnames(final_reversal_df)) final_reversal_df$action <- "Unknown"
  if (!"IC50" %in% colnames(final_reversal_df)) final_reversal_df$IC50 <- 99999
  
  final_reversal_df$reversal_score_scaled <- case_when(
    final_reversal_df$action %in% c("INHIBITOR", "AGONIST") ~ final_reversal_df$reversal_score,
    final_reversal_df$action %in% c("ANTAGONIST") ~ -final_reversal_df$reversal_score,
    final_reversal_df$action == "Unknown" & final_reversal_df$IC50 < 10 ~ final_reversal_df$reversal_score,
    TRUE ~ final_reversal_df$reversal_score
  )
  
  final_reversal_df$confidence <- case_when(
    final_reversal_df$action %in% c("INHIBITOR", "AGONIST", "ANTAGONIST") ~ "high",
    final_reversal_df$action == "Unknown" & final_reversal_df$IC50 < 10 ~ "high",
    TRUE ~ "low"
  )
  
  final_reversal_df$significant <- ifelse(
    final_reversal_df$down_padj <= 0.05 & final_reversal_df$up_padj <= 0.05, 
    TRUE, 
    FALSE
  )
} else {
  final_reversal_df <- reversal_df
  final_reversal_df$reversal_score_scaled <- final_reversal_df$reversal_score
  final_reversal_df$confidence <- "low"
  final_reversal_df$significant <- FALSE
}

if ("smiles_id" %in% colnames(final_reversal_df)) {
  final_reversal_df <- final_reversal_df %>% dplyr::rename(SMILES = smiles_id)
}

dir.create(dirname(opt$output_reversal_csv), showWarnings = FALSE, recursive = TRUE)
fwrite(final_reversal_df, opt$output_reversal_csv)
message("Successfully wrote disease reversal results to: ", opt$output_reversal_csv)
