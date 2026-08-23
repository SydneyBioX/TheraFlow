# Setup command line arguments using standard commandArgs
args <- commandArgs(trailingOnly = TRUE)
opt <- list(
  targets_csv = NULL,
  potency_csv = NULL,
  output_reversal_csv = NULL,
  ppi_rds = "data/ppi.rds",
  mapping_file = "data/chembl_uniprot_mapping.txt",
  lincs_meta = "Omics2Drug/data/lincs_gene_meta.csv",
  test_datasets_dir = "Omics2Drug/data/test_datasets",
  workers = 4
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
  } else if (arg == "--output_reversal_csv" && i < length(args)) {
    opt$output_reversal_csv <- args[i + 1]
    i <- i + 2
  } else if (arg == "--ppi_rds" && i < length(args)) {
    opt$ppi_rds <- args[i + 1]
    i <- i + 2
  } else if (arg == "--mapping_file" && i < length(args)) {
    opt$mapping_file <- args[i + 1]
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

if (is.null(opt$targets_csv) || is.null(opt$potency_csv) || is.null(opt$output_reversal_csv)) {
  stop("Missing required arguments: --targets_csv, --potency_csv, --output_reversal_csv. Usage: Rscript disease_reversal.R --targets_csv <path> --potency_csv <path> --output_reversal_csv <path>", call. = FALSE)
}

library(furrr)
library(tidyr)
library(dplyr)
library(data.table)
library(tidyverse)
library(biomaRt)
library(readr)
library(DTSEA)
library(igraph)

# -------------------------------------------------------------------------
# 1. Helper function: process LINCS signature files
# -------------------------------------------------------------------------
process_lincs_signatures <- function(path, lincs_meta) {
  # fread by default reads the glaucoma.csv as header=FALSE due to mixed types
  data <- fread(path)
  colnames(data) <- as.character(data[1, ])
  
  # Take row 2 (first row of actual z-scores) and columns 2 to end (skip sig_id column)
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
      by = 'pr_gene_id'
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
# 2. Helper function: load disease signature dynamically
# -------------------------------------------------------------------------
get_disease_signature <- function(test_dataset, test_datasets_dir, lincs_meta) {
  # 1. Check for Excape / control-perturbed limma pipeline
  if (test_dataset == "excape_1" || test_dataset == "Excape_001") {
    ctl_file <- file.path(test_datasets_dir, "gexmolgen", "Excape_001_ctl_AKT2.csv")
    pert_file <- file.path(test_datasets_dir, "gexmolgen", "Excape_001_pert_AKT2.csv")
    
    if (file.exists(ctl_file) && file.exists(pert_file)) {
      message("Found control/perturbed files for Excape. Running limma differential expression...")
      library(limma)
      library(tibble)
      AKT2_ctl <- read_csv(ctl_file, show_col_types = FALSE)
      AKT2_pert <- read_csv(pert_file, show_col_types = FALSE)
      
      exp_matrix <- rbind(
        as.matrix(dplyr::select(AKT2_ctl, -1)), 
        as.matrix(dplyr::select(AKT2_pert, -1))
      ) %>% t()
      
      group <- factor(c(rep("Control", 100), rep("Perturbed", 100)))
      design <- model.matrix(~ group)
      fit <- lmFit(exp_matrix, design)
      fit <- eBayes(fit)
      
      de_genes <- topTable(fit, coef = 2, number = Inf, sort.by = "P") %>% 
        rownames_to_column('gene') %>%
        filter(P.Value <= 0.05) %>%
        dplyr::select(gene, logFC)
        
      return(list(
        up = de_genes %>% filter(logFC > 0) %>% deframe(),
        down = de_genes %>% filter(logFC < 0) %>% deframe()
      ))
    }
  }
  
  # 2. Check for standard dataset CSV files
  matched_file <- NULL
  subdirs <- c("smilesgen", "gexmolgen", "transgem", "smilesgen_old")
  for (sd in subdirs) {
    p <- file.path(test_datasets_dir, sd, paste0(test_dataset, ".csv"))
    if (file.exists(p)) {
      matched_file <- p
      break
    }
    # Casing variant (e.g. Lincs1000_1_m.csv vs lincs1000_1_m.csv)
    cap_name <- paste0(toupper(substr(test_dataset, 1, 1)), substr(test_dataset, 2, nchar(test_dataset)))
    p_cap <- file.path(test_datasets_dir, sd, paste0(cap_name, ".csv"))
    if (file.exists(p_cap)) {
      matched_file <- p_cap
      break
    }
  }
  
  if (is.null(matched_file)) {
    # Recursive search fallback
    all_files <- list.files(test_datasets_dir, recursive = TRUE, full.names = TRUE, pattern = "\\.csv$")
    base_names <- basename(all_files)
    idx <- which(tolower(base_names) == tolower(paste0(test_dataset, ".csv")))
    if (length(idx) > 0) {
      matched_file <- all_files[idx[1]]
    }
  }
  
  if (is.null(matched_file)) {
    stop(paste("Disease signature file not found for dataset:", test_dataset))
  }
  
  message("Loading disease signature: ", matched_file)
  return(process_lincs_signatures(matched_file, lincs_meta))
}

# -------------------------------------------------------------------------
# 3. Main Data Processing Workflow
# -------------------------------------------------------------------------
message("Loading PPI network: ", opt$ppi_rds)
full_human_ppi <- readRDS(opt$ppi_rds)

message("Loading LINCS gene metadata: ", opt$lincs_meta)
lincs_meta <- fread(opt$lincs_meta)

message("Loading targets and properties: ", opt$targets_csv)
metrics <- read_csv(opt$targets_csv, show_col_types = FALSE)

# Force method_dataset column if missing
if (!"method_dataset" %in% colnames(metrics)) {
  metrics$method_dataset <- paste0(metrics$method, "_", metrics$test_dataset)
}

message("Loading ChEMBL to UniProt mapping: ", opt$mapping_file)
mapping_uniprot <- fread(opt$mapping_file, header = FALSE, sep = "\t")
colnames(mapping_uniprot) <- c("uniprot_id", "chembl_id", "description", "protein_type")
lookup_vec <- setNames(mapping_uniprot$uniprot_id, mapping_uniprot$description)

# Canonicalize function
canonicalize <- function(smiles) {
  tryCatch({
    library(rdkit) # if available
    # fallback to R/Python or assume pre-canonicalized
    smiles
  }, error = function(e) {
    smiles
  })
}

# Extract targets list and map to UniProt ID
df_clean <- metrics %>%
  dplyr::select(SMILES, Predicted_Targets) %>%
  distinct(SMILES, .keep_all = TRUE) %>%
  mutate(description = gsub("\\[|\\]|'", "", Predicted_Targets)) %>% 
  separate_rows(description, sep = ", ") %>%
  filter(!description %in% c("", "Unchecked", "[]")) %>%
  mutate(uniprot_id = lookup_vec[description]) %>%
  filter(!is.na(uniprot_id)) %>%
  distinct(SMILES, uniprot_id, .keep_all = TRUE)

if (nrow(df_clean) == 0) {
  stop("No valid targets could be mapped to UniProt IDs.")
}

# Map UniProt SwissProt ID to HGNC symbol (with caching)
cache_file <- "data/biomart_uniprot_hgnc_cache.rds"
dir.create(dirname(cache_file), showWarnings = FALSE, recursive = TRUE)

if (file.exists(cache_file)) {
  message("Loading biomaRt conversion from cache...")
  conversion <- readRDS(cache_file)
} else {
  message("Querying Ensembl biomaRt for UniProt to HGNC mapping...")
  conversion <- tryCatch({
    mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    res <- getBM(
      attributes = c("uniprotswissprot", "hgnc_symbol"),
      filters = "uniprotswissprot", 
      values = unique(df_clean$uniprot_id), 
      mart = mart
    )
    saveRDS(res, cache_file)
    res
  }, error = function(e) {
    message("biomaRt query failed: ", e$message)
    message("Trying AnnotationDbi/org.Hs.eg.db fallback...")
    tryCatch({
      library(org.Hs.eg.db)
      res_fallback <- select(
        org.Hs.eg.db, 
        keys = unique(df_clean$uniprot_id), 
        columns = c("SYMBOL"), 
        keytype = "UNIPROT"
      )
      colnames(res_fallback) <- c("uniprotswissprot", "hgnc_symbol")
      res_fallback
    }, error = function(e2) {
      message("Fallback failed: ", e2$message)
      data.frame(uniprotswissprot = character(), hgnc_symbol = character(), stringsAsFactors = FALSE)
    })
  })
}

df_clean <- df_clean %>% 
  inner_join(conversion, by = join_by(uniprot_id == uniprotswissprot))

df_new <- df_clean %>%
  na.omit() %>% 
  group_by(SMILES) %>%
  summarise(hgnc_symbol = paste(hgnc_symbol, collapse = ","), .groups = "drop")

disease_reversal <- metrics %>% 
  left_join(df_new, by = 'SMILES') %>%
  filter(!is.na(hgnc_symbol))

if (nrow(disease_reversal) == 0) {
  stop("No molecules had HGNC symbols mapped. Cannot run DTSEA.")
}

disease_reversal$drug_id <- 1:nrow(disease_reversal)

drugs <- disease_reversal %>% dplyr::select(drug_id, SMILES, method_dataset, hgnc_symbol)
drugs_long <- drugs %>%
  separate_rows(hgnc_symbol, sep = ",\\s*") %>%
  dplyr::rename(target = hgnc_symbol)

# -------------------------------------------------------------------------
# 4. DTSEA Execution Loop
# -------------------------------------------------------------------------
dtsea_function <- function(drugs_long, m_ds, up, down) {
  message("Processing dataset: ", m_ds)
  drugs_subset <- drugs_long %>% 
    filter(method_dataset == m_ds) %>% 
    dplyr::select(SMILES, target)
  colnames(drugs_subset) <- c("drug_id", "gene_target")
  
  if (nrow(drugs_subset) == 0) {
    return(data.frame(smiles_id = character(), reversal_score = numeric(), up_padj = numeric(), down_padj = numeric()))
  }
  
  # enrichment for up-regulated signatures
  up_result <- DTSEA(
    network = full_human_ppi,
    disease = names(up),
    drugs = drugs_subset,
    verbose = FALSE
  )
  
  # enrichment for down-regulated signatures
  down_result <- DTSEA(
    network = full_human_ppi,
    disease = names(down),
    drugs = drugs_subset,
    verbose = FALSE
  )
  
  up_df <- data.frame(smiles_id = up_result$drug_id, nes_up = up_result$NES, up_padj = up_result$padj)
  down_df <- data.frame(smiles_id = down_result$drug_id, nes_down = down_result$NES, down_padj = down_result$padj)
  
  combined_df <- merge(up_df, down_df, by = "smiles_id", all = TRUE)
  combined_df$nes_up[is.na(combined_df$nes_up)] <- 0
  combined_df$nes_down[is.na(combined_df$nes_down)] <- 0
  combined_df$reversal_score <- combined_df$nes_up - combined_df$nes_down
  
  final_df <- combined_df[, c("smiles_id", "reversal_score", "up_padj", "down_padj")]
  return(final_df)
}

methods_datasets <- unique(drugs_long$method_dataset)

message("Setting up parallel multisession with ", opt$workers, " workers...")
plan(multisession, workers = opt$workers)

all_results <- seq_along(methods_datasets) %>%
  future_map(~ {
    m_ds <- methods_datasets[.x]
    test_dataset <- sub("^[^_]+_", "", m_ds)
    sig <- get_disease_signature(test_dataset, opt$test_datasets_dir, lincs_meta)
    
    dtsea_function(
      drugs_long = drugs_long, 
      m_ds = m_ds, 
      up = sig$up, 
      down = sig$down
    )
  }, .options = furrr_options(seed = TRUE)) %>%
  set_names(methods_datasets)

reversal_df <- bind_rows(all_results, .id = "method_dataset")

# -------------------------------------------------------------------------
# 5. Potency Scoring and Scaling
# -------------------------------------------------------------------------
message("Reading potency scores: ", opt$potency_csv)
potency <- fread(opt$potency_csv)

if (nrow(potency) > 0) {
  # Standardize colnames
  if (!"query_smiles" %in% colnames(potency) && "SMILES" %in% colnames(potency)) {
    setnames(potency, "SMILES", "query_smiles")
  }
  
  # Join DTSEA with potency data
  final_reversal_df <- reversal_df %>% left_join(potency, by = join_by(smiles_id == query_smiles))
  
  # Handle missing cols gracefully
  if (!"action" %in% colnames(final_reversal_df)) {
    final_reversal_df$action <- "Unknown"
  }
  if (!"IC50" %in% colnames(final_reversal_df)) {
    final_reversal_df$IC50 <- 99999.0
  }
  
  final_reversal_df <- final_reversal_df %>%
    mutate(
      reversal_score_scaled = case_when(
        action %in% c("INHIBITOR", "ANTAGONIST", "BLOCKER") ~ reversal_score,
        action %in% c("AGONIST") ~ -reversal_score,
        action == "Unknown" & IC50 < 10 ~ reversal_score,
        .default = reversal_score
      ),
      confidence = case_when(
        action %in% c("INHIBITOR", "AGONIST", "ANTAGONIST") ~ "high",
        action == "Unknown" & IC50 < 10 ~ "high",
        .default = "low"
      ),
      significant = ifelse(down_padj <= 0.05 & up_padj <= 0.05, TRUE, FALSE)
    )
} else {
  message("Potency file was empty. Skipping potency scaling.")
  final_reversal_df <- reversal_df %>%
    mutate(
      reversal_score_scaled = reversal_score,
      confidence = "low",
      significant = ifelse(down_padj <= 0.05 & up_padj <= 0.05, TRUE, FALSE),
      action = "Unknown",
      IC50 = NA,
      EC50 = NA
    )
}

setnames(final_reversal_df, "smiles_id", "SMILES")

# Write output file
dir.create(dirname(opt$output_reversal_csv), showWarnings = FALSE, recursive = TRUE)
fwrite(final_reversal_df, opt$output_reversal_csv)
message("Successfully wrote disease reversal results to ", opt$output_reversal_csv)
