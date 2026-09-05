#!/usr/bin/env Rscript
# Step 6: Multi-objective Pareto optimization & RRA Consensus Ranking + AF3 JSON Prep
# Logic directly sourced from metrics_calculation/7.2. result5.qmd (lines 231 to 484)

args <- commandArgs(trailingOnly = TRUE)
opt <- list(
  master_csv = "TheraFlow/outputs/merged_master_dataframe.csv",
  output_consensus = "TheraFlow/outputs/consensus_top.csv",
  output_af3_json = "TheraFlow/outputs/af3_batch_input.json",
  top_n = 10
)

i <- 1
while (i <= length(args)) {
  arg <- args[i]
  if (arg == "--master_csv" && i < length(args)) {
    opt$master_csv <- args[i + 1]
    i <- i + 2
  } else if (arg == "--output_consensus" && i < length(args)) {
    opt$output_consensus <- args[i + 1]
    i <- i + 2
  } else if (arg == "--output_af3_json" && i < length(args)) {
    opt$output_af3_json <- args[i + 1]
    i <- i + 2
  } else if (arg == "--top_n" && i < length(args)) {
    opt$top_n <- as.integer(args[i + 1])
    i <- i + 2
  } else {
    i <- i + 1
  }
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(data.table)
  library(tidyr)
  library(rPref)
  library(RobustRankAggreg)
  library(jsonlite)
})

message("Loading merged master dataframe from: ", opt$master_csv)
if (!file.exists(opt$master_csv)) {
  stop(paste("Master dataframe not found at:", opt$master_csv))
}
df <- read_csv(opt$master_csv, show_col_types = FALSE)

# -------------------------------------------------------------------------
# Helper functions for geometric metrics calculation
# -------------------------------------------------------------------------
safe_div <- function(a, b) {
  ifelse(is.na(b) | b == 0, 0, a / b)
}

clip01 <- function(x) {
  pmin(pmax(x, 0), 1)
}

normalize_count <- function(x, cap) {
  clip01(x / cap)
}

# -------------------------------------------------------------------------
# Combine geometry metrics
# -------------------------------------------------------------------------
# Ensure required columns exist with defaults if PoseBusters columns were omitted
ensure_col <- function(data, colname, default_val = 0) {
  if (!colname %in% colnames(data)) {
    data[[colname]] <- default_val
  }
  return(data)
}

pb_cols <- c(
  "inchi_convertible", "all_atoms_connected", "no_radicals",
  "number_valid_bonds", "number_bonds", "number_valid_angles", "number_angles",
  "number_valid_noncov_pairs", "number_noncov_pairs",
  "number_short_outlier_bonds", "number_long_outlier_bonds", "number_outlier_angles",
  "number_clashes", "number_aromatic_rings_pass", "number_aromatic_rings_checked",
  "number_non-aromatic_rings_pass", "number_non-aromatic_rings_checked",
  "aromatic_ring_maximum_distance_from_plane", "non-aromatic_ring_maximum_distance_from_plane",
  "number_double_bonds_pass", "number_double_bonds_checked", "double_bond_maximum_distance_from_plane",
  "energy_ratio", "ensemble_avg_energy", "Passes_Safety"
)
for (c in pb_cols) {
  df <- ensure_col(df, c, default_val = 0)
}

df <- df %>%
  mutate(
    # -------------------------
    # structural validity
    # -------------------------
    structural_validity = rowMeans(cbind(
      as.numeric(inchi_convertible),
      as.numeric(all_atoms_connected),
      as.numeric(no_radicals),
      safe_div(number_valid_bonds, number_bonds),
      safe_div(number_valid_angles, number_angles),
      safe_div(number_valid_noncov_pairs, number_noncov_pairs)
    ), na.rm = TRUE),

    # -------------------------
    # bond + angle consistency
    # -------------------------
    bond_angle_consistency = 1 - clip01(
      (normalize_count(number_short_outlier_bonds, 10) +
       normalize_count(number_long_outlier_bonds, 10) +
       normalize_count(number_outlier_angles, 10)) / 3
    ),

    # -------------------------
    # steric feasibility
    # -------------------------
    steric_feasibility = 1 - clip01(
      (normalize_count(number_clashes, 10) +
       normalize_count(number_noncov_pairs, 20)) / 2
    ),

    # -------------------------
    # ring planarity
    # -------------------------
    ring_planarity = rowMeans(cbind(
      safe_div(number_aromatic_rings_pass, number_aromatic_rings_checked),
      safe_div(`number_non-aromatic_rings_pass`, `number_non-aromatic_rings_checked`),
      1 - clip01(aromatic_ring_maximum_distance_from_plane / 1),
      1 - clip01(`non-aromatic_ring_maximum_distance_from_plane` / 1)
    ), na.rm = TRUE),

    # -------------------------
    # double bond geometry
    # -------------------------
    double_bond_geometry = rowMeans(cbind(
      safe_div(number_double_bonds_pass, number_double_bonds_checked),
      1 - clip01(double_bond_maximum_distance_from_plane / 1)
    ), na.rm = TRUE)
  )

# Chemical validity weighted sanitization
if (file.exists("TheraFlow/outputs/metrics_metadata.csv")) {
  meta <- fread("TheraFlow/outputs/metrics_metadata.csv")
  selected_validity <- meta[Classification == "Chemical validity", ]$Metric 
  validity_present <- intersect(selected_validity, colnames(df))
  if (length(validity_present) > 0) {
    df <- df %>% mutate(weighted_sanitization = rowMeans(across(all_of(validity_present)), na.rm = TRUE))
  } else {
    df$weighted_sanitization <- 1.0
  }
} else {
  validity_cols <- intersect(c("sanitization", "passes_valence_checks", "passes_kekulization", "no_radicals_before_sanitization", "no_radicals"), colnames(df))
  if (length(validity_cols) > 0) {
    df <- df %>% mutate(weighted_sanitization = rowMeans(across(all_of(validity_cols)), na.rm = TRUE))
  } else {
    df$weighted_sanitization <- 1.0
  }
}

# -------------------------------------------------------------------------
# Pareto Optimization with rPref
# -------------------------------------------------------------------------
selected_metrics <- c(
  "QED", "SA_Score", "MW", "LogP", "TPSA", "HBA", "HBD", "F_Csp3",
  "structural_validity", "bond_angle_consistency", "steric_feasibility", 
  "double_bond_geometry", "weighted_sanitization", "Sim_Morg_Tani", 
  "energy_ratio", "ensemble_avg_energy"
)

# Fill missing metrics with neutral defaults
for (m in selected_metrics) {
  if (!m %in% colnames(df)) {
    df[[m]] <- 0.0
  } else {
    df[[m]][is.na(df[[m]])] <- 0.0
  }
}

data_ready <- df %>%
  filter(
    (number_clashes == 0 | is.na(number_clashes)) & (Passes_Safety == TRUE | is.na(Passes_Safety))
  )
if (nrow(data_ready) == 0) {
  data_ready <- df
}

p_11d <- high(QED) * low(SA_Score) * low(MW) * low(LogP) * low(TPSA) * low(HBA) * low(HBD) * 
         high(bond_angle_consistency) * high(steric_feasibility) * high(double_bond_geometry) * 
         high(weighted_sanitization) * high(Sim_Morg_Tani) * low(energy_ratio) * 
         low(ensemble_avg_energy) * high(F_Csp3)

pareto_11d_results <- tryCatch({
  psel(data_ready, p_11d)
}, error = function(e) {
  message("rPref psel warning: ", e$message)
  data_ready
})

# -------------------------------------------------------------------------
# Robust Rank Aggregation (RRA)
# -------------------------------------------------------------------------
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
    if (!metric %in% colnames(data)) next
    
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

if ("method" %in% colnames(df) && "test_dataset" %in% colnames(df)) {
  final_rankings <- df %>%
    group_by(method, test_dataset) %>%
    group_modify(~ perform_rra_ranking(.x, polarity)) %>%
    ungroup()

  consensus <- pareto_11d_results %>%
    dplyr::select(any_of(c("method", "test_dataset", "SMILES"))) %>%
    inner_join(final_rankings, by = intersect(c("method", "test_dataset", "SMILES"), colnames(pareto_11d_results))) %>%
    group_by(method, test_dataset) %>%
    arrange(RRA_Score, .by_group = TRUE) %>%
    slice_head(n = opt$top_n) %>%
    ungroup()
} else {
  final_rankings <- perform_rra_ranking(df, polarity)
  consensus <- pareto_11d_results %>%
    inner_join(final_rankings, by = "SMILES") %>%
    arrange(RRA_Score) %>%
    slice_head(n = opt$top_n)
}

dir.create(dirname(opt$output_consensus), showWarnings = FALSE, recursive = TRUE)
write.csv(consensus, opt$output_consensus, row.names = FALSE)
message("Saved top consensus candidates to: ", opt$output_consensus)

# -------------------------------------------------------------------------
# AF3 Target Sequences & JSON Generation
# -------------------------------------------------------------------------
glaucoma_target <- "MSMNNSKQPVSPAAGLLSNTTCQTENRLSVFFSIIFMTVGILSNSLAIAILMKAYQRFRQKSKASFLLLASALVITDFFGHLINGAIAVFVYASDKDWIRFDQSNVLCSIFGICMVFSGLCPLLLGSAMAIERCIGVTKPIFHSTKITSKHVKIMLSGVCLFAVFIALLPILGNRDYKIQASRTWCFYNTEDIKDWEDRFYLLLFSFLGLLALVVSFLCNAITGITLLRVKFKNQQHRQGRSHHFEMVIQLLAIMCVSCICWSPFLVTMANIGINGNHSLETCETTLFALRMATWNQILDPWVYILLRKAVLKNLYKLARRCCGVHVISLHVWELSSIKNSLKVAAISESPVAEKINQQAPSLIGQ"
cardiomyopathy_target <- "MGDSEMAVFGAAAPYLRKSEKERLEAQTRPFDLKKDVFVPDDKQEFVKAKIVSREGGKVTAETEYGKTVTVKEDQVMQQNPPKFDKIEDMAMLTFLHEPAVLYNLKDRYGSWMIYTYSGLFCVTVNPYKWLPVYTPEVVAAYRGKKRSEAPPHIFSISDNAYQYMLTDRENQSILITGESGAGKTVNTKRVIQYFAVIAAIGDRSKKDQSPGKGTLEDQIIQANPALEAFGNAKTVRNDNSSRFGKFIRIHFGATGKLASADIETYLLEKSRVIFQLKAERDYHIFYQILSNKKPELLDMLLITNNPYDYAFISQGETTVASIDDAEELMATDNAFDVLGFTSEEKNSMYKLTGAIMHFGNMKFKLKQREEQAEPDGTEEADKSAYLMGLNSADLLKGLCHPRVKVGNEYVTKGQNVQQVIYATGALAKAVYERMFNWMVTRINATLETKQPRQYFIGVLDIAGFEIFDFNSFEQLCINFTNEKLQQFFNHHMFVLEQEEYKKEGIEWTFIDFGMDLQACIDLIEKPMGIMSILEEECMFPKATDMTFKAKLFDNHLGKSANFQKPRNIKGKPEAHFSLIHYAGIVDYNIIGWLQKNKDPLNETVVGLYQKSSLKLLSTLFANYAGADAPIEKGKGKAKKGSSFQTVSALHRENLNKLMTNLRSTHPHFVRCIIPNETKSPGVMDNPLVMHQLRCNGVLEGIRICRKGFPNRILYGDFRQRYRILNPAAIPEGQFIDSRKGAEKLLSSLDIDHNQYKFGHTKVFFKAGLLGLLEEMRDERLSRIITRIQAQSRGVLARMEYKKLLERRDSLLVIQWNIRAFMGVKNWPWMKLYFKIKPLLKSAEREKEMASMKEEFTRLKEALEKSEARRKELEEKMVSLLQEKNDLQLQVQAEQDNLADAEERCDQLIKNKIQLEAKVKEMNERLEDEEEMNAELTAKKRKLEDECSELKRDIDDLELTLAKVEKEKHATENKVKNLTEEMAGLDEIIAKLTKEKKALQEAHQQALDDLQAEEDKVNTLTKAKVKLEQQVDDLEGSLEQEKKVRMDLERAKRKLEGDLKLTQESIMDLENDKQQLDERLKKKDFELNALNARIEDEQALGSQLQKKLKELQARIEELEEELEAERTARAKVEKLRSDLSRELEEISERLEEAGGATSVQIEMNKKREAEFQKMRRDLEEATLQHEATAAALRKKHADSVAELGEQIDNLQRVKQKLEKEKSEFKLELDDVTSNMEQIIKAKANLEKMCRTLEDQMNEHRSKAEETQRSVNDLTSQRAKLQTENGELSRQLDEKEALISQLTRGKLTYTQQLEDLKRQLEEEVKAKNALAHALQSARHDCDLLREQYEEETEAKAELQRVLSKANSEVAQWRTKYETDAIQRTEELEEAKKKLAQRLQEAEEAVEAVNAKCSSLEKTKHRLQNEIEDLMVDVERSNAAAAALDKKQRNFDKILAEWKQKYEESQSELESSQKEARSLSTELFKLKNAYEESLEHLETFKRENKNLQEEISDLTEQLGSSGKTIHELEKVRKQLEAEKMELQSALEEAEASLEHEEGKILRAQLEFNQIKAEIERKLAEKDEEMEQAKRNHLRVVDSLQTSLDAETRSRNEALRVKKKMEGDLNEMEIQLSHANRMAAEAQKQVKSLQSLLKDTQIQLDDAVRANDDLKENIAIVERRNNLLQAELEELRAVVEQTERSRKLAEQELIETSERVQLLHSQNTSLINQKKKMDADLSQLQTEVEEAVQECRNAEEKAKKAITDAAMMAEELKKEQDTSAHLERMKKNMEQTIKDLQHRLDEAEQIALKGGKKQLQKLEARVRELENELEAEQKRNAESVKGMRKSERRIKELTYQTEEDRKNLLRLQDLVDKLQLKVKAYKRQAEEAEEQANTNLSKFRKVQHELDEAEERADIAESQVNKLRAKSRDIGTKGLNEE"
parkinsons_target <- "MDVFMKGLSKAKEGVVAAAEKTKQGVAEAAGKTKEGVLYVGSKTKEGVVHGVATVAEKTKEQVTNVGGAVVTGVTAVAQKTVEGAGSIAAATGFVKKDQLGKNEEGAPQEGILEDMPVDPDNEAYEMPSEEGYQDYEPEA"
hcc_target <- "MSPSAKKRPKNSRVSKMQDEKLRDETEQPVSKVIERNRLRTVLKNLSLLKLLKSSNRRIQELHKLAKRCWHSLLSVPKILRISSGENSACNKTKQNNEEFQEIGCSEKELKSKKLESTGDPKKKEYKEWKSQVQSGMRNKEKTSLAAMPRKEKHIEPEVPRTSRDDSLNPGVQGRQPLTEGPRVIFIKPYRNRTPMGHMKQLDVADQWIWFEGLPTRIHLPAPRVMCRSSTLRWVKRRCTRFCSASLEMPMWHPYKVDVTWTRARGASRGWRSRHQLKGRNGWRNSRVYK"

if ("test_dataset" %in% colnames(consensus)) {
  consensus$target_sequence <- case_when(
    consensus$test_dataset == "glaucoma" ~ glaucoma_target,
    consensus$test_dataset == "cardiomyopathy" ~ cardiomyopathy_target,
    consensus$test_dataset == "glaucoma_bulk" ~ glaucoma_target,
    consensus$test_dataset == "hcc" ~ hcc_target,
    consensus$test_dataset == "parkinsons" ~ parkinsons_target,
    TRUE ~ glaucoma_target
  )
} else {
  consensus$target_sequence <- glaucoma_target
}

if (!"method" %in% colnames(consensus)) consensus$method <- "method"
if (!"test_dataset" %in% colnames(consensus)) consensus$test_dataset <- "dataset"

af3_jobs <- lapply(1:nrow(consensus), function(i) {
  list(
    name = paste0("docking_", consensus$method[i], "_", consensus$test_dataset[i], "_", i),
    modelSeeds = list(1),
    sequences = list(
      list(
        protein = list(
          id = "A",
          sequence = consensus$target_sequence[i],
          count = 1
        )
      ),
      list(
        ligand = list(
          id = "B",
          smiles = consensus$SMILES[i],
          count = 1
        )
      )
    )
  )
})

dir.create(dirname(opt$output_af3_json), showWarnings = FALSE, recursive = TRUE)
write_json(af3_jobs, opt$output_af3_json, auto_unbox = TRUE, pretty = TRUE)
message("Generated AF3 batch JSON with ", length(af3_jobs), " docking jobs to: ", opt$output_af3_json)
