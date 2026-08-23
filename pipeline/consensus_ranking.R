# R script for step 8: consensus ranking using Pareto front tiers and Robust Rank Aggregation (RRA)
args <- commandArgs(trailingOnly = TRUE)
opt <- list(
  input_csv = "results/pipeline_runs/step7_merged_master_dataframe.csv",
  output_csv = "results/pipeline_runs/step8_consensus_top10.csv"
)

i <- 1
while (i <= length(args)) {
  arg <- args[i]
  if (arg == "--input_csv" && i < length(args)) {
    opt$input_csv <- args[i + 1]
    i <- i + 2
  } else if (arg == "--output_csv" && i < length(args)) {
    opt$output_csv <- args[i + 1]
    i <- i + 2
  } else {
    i <- i + 1
  }
}

library(dplyr)
library(readr)
library(data.table)
library(rPref)
library(RobustRankAggreg)

message("Loading merged master dataframe: ", opt$input_csv)
df <- read_csv(opt$input_csv, show_col_types = FALSE)

if (nrow(df) == 0) {
  stop("Input master dataframe is empty. Cannot perform consensus ranking.")
}

# Define the 15 metrics to optimize
selected_metrics <- c(
  "QED", "SA_Score", "MW", "LogP", "TPSA", "HBA", "HBD", "F_Csp3",
  "structural_validity", "bond_angle_consistency", "steric_feasibility", 
  "double_bond_geometry", "weighted_sanitization", "Sim_Morg_Tani", 
  "energy_ratio", "ensemble_avg_energy"
)

# Calculate weighted_sanitization from Chemical Validity columns
validity_cols <- c("sanitization", "passes_valence_checks", "passes_kekulization", "no_radicals_before_sanitization")
present_val_cols <- intersect(validity_cols, colnames(df))

if (length(present_val_cols) > 0) {
  message("Calculating weighted_sanitization from: ", paste(present_val_cols, collapse = ", "))
  df <- df %>%
    mutate(weighted_sanitization = rowMeans(across(all_of(present_val_cols)), na.rm = TRUE))
} else {
  message("Warning: No chemical validity columns found. Setting weighted_sanitization to 1.0")
  df$weighted_sanitization <- 1.0
}

# Set default values for missing metrics to prevent psel/RRA failures
for (metric in selected_metrics) {
  if (!metric %in% colnames(df)) {
    if (metric %in% c("QED", "F_Csp3", "structural_validity", "bond_angle_consistency", "steric_feasibility", "double_bond_geometry", "weighted_sanitization", "Sim_Morg_Tani")) {
      df[[metric]] <- 0.0 # higher is better, default to 0
    } else {
      df[[metric]] <- 99999.0 # lower is better, default to high value
    }
  }
}

# Force boolean / numeric casts
df <- df %>%
  mutate(
    Passes_Safety = as.logical(Passes_Safety),
    number_clashes = as.numeric(number_clashes)
  )

# Default Passes_Safety and number_clashes if missing
if (!"Passes_Safety" %in% colnames(df)) df$Passes_Safety <- TRUE
if (!"number_clashes" %in% colnames(df)) df$number_clashes <- 0

df$Passes_Safety[is.na(df$Passes_Safety)] <- TRUE
df$number_clashes[is.na(df$number_clashes)] <- 0

# Filter for initial safety
message("Applying safety filters (number_clashes == 0 and Passes_Safety == TRUE)...")
data_ready <- df %>%
  filter(number_clashes == 0 & Passes_Safety == TRUE)

if (nrow(data_ready) == 0) {
  message("Warning: No molecules passed the safety filters. Falling back to entire dataset for ranking.")
  data_ready <- df
}

# -------------------------------------------------------------------------
# 1. Pareto Front Ranking
# -------------------------------------------------------------------------
message("Running Pareto front optimization...")
p_11d <- high(QED) * low(SA_Score) * low(MW) * low(LogP) * low(TPSA) * low(HBA) * 
         low(HBD) * high(bond_angle_consistency) * high(steric_feasibility) * 
         high(double_bond_geometry) * high(weighted_sanitization) * 
         high(Sim_Morg_Tani) * low(energy_ratio) * low(ensemble_avg_energy) * 
         high(F_Csp3)

pareto_results <- psel(data_ready, p_11d, top = nrow(data_ready))

# -------------------------------------------------------------------------
# 2. Robust Rank Aggregation (RRA)
# -------------------------------------------------------------------------
message("Running Robust Rank Aggregation (RRA)...")
polarity <- c(
  "QED"                     = 1,   # high
  "SA_Score"                = -1,  # low
  "MW"                      = -1,  # low
  "LogP"                    = -1,  # low
  "TPSA"                    = -1,  # low
  "HBA"                     = -1,  # low
  "HBD"                     = -1,  # low
  "F_Csp3"                  = 1,   # high
  "structural_validity"     = 1,   # high
  "bond_angle_consistency"  = 1,   # high
  "steric_feasibility"      = 1,   # high
  "double_bond_geometry"    = 1,   # high
  "weighted_sanitization"   = 1,   # high
  "Sim_Morg_Tani"           = 1,   # high
  "energy_ratio"            = -1,  # low
  "ensemble_avg_energy"     = -1   # low
)

perform_rra_ranking <- function(data, metric_map) {
  N_total <- length(unique(data$SMILES))
  all_rank_lists <- list()
  
  for (metric in names(metric_map)) {
    temp_df <- data %>%
      dplyr::select(SMILES, !!sym(metric)) %>%
      filter(!is.na(!!sym(metric)))
    
    if (nrow(temp_df) == 0) next
    
    if (metric_map[metric] == 1) {
      sorted_smiles <- temp_df$SMILES[order(-temp_df[[metric]], temp_df$SMILES)]
    } else {
      sorted_smiles <- temp_df$SMILES[order(temp_df[[metric]], temp_df$SMILES)]
    }
    all_rank_lists[[metric]] <- sorted_smiles
  }
  
  if (length(all_rank_lists) == 0) return(tibble(SMILES = character(), RRA_Score = numeric()))
  
  result <- RobustRankAggreg::aggregateRanks(glist = all_rank_lists, method = "RRA", N = N_total)
  colnames(result) <- c("SMILES", "RRA_Score")
  return(as_tibble(result))
}

# Run RRA grouped by method and test_dataset
final_rankings <- data_ready %>%
  group_by(method, test_dataset) %>%
  group_modify(~ perform_rra_ranking(.x, polarity)) %>%
  ungroup()

# -------------------------------------------------------------------------
# 3. Consensus Sorting (Rank by Pareto Tier then RRA Score)
# -------------------------------------------------------------------------
message("Merging Pareto and RRA rankings to select top 10 consensus molecules per group...")
consensus <- pareto_results %>%
  dplyr::select(method, test_dataset, SMILES, .level) %>%
  inner_join(final_rankings, by = c("method", "test_dataset", "SMILES")) %>%
  group_by(method, test_dataset) %>%
  arrange(.level, RRA_Score, .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  ungroup()

# Write consensus results
dir.create(dirname(opt$output_csv), showWarnings = FALSE, recursive = TRUE)
write.csv(consensus, opt$output_csv, row.names = FALSE)
message("Successfully wrote top 10 consensus rankings to ", opt$output_csv)
