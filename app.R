# ==============================================================================
# TheraFlow: Interactive Multi-Stage Generative Molecule Evaluation Dashboard
# ==============================================================================
# An R Shiny web application featuring rich interactive Plotly visualizations,
# recreating all figures and analyses from TheraFlow/Report.qmd.
# ==============================================================================

library(shiny)
library(bslib)
library(plotly)
library(DT)
library(funkyheatmap)
library(data.table)
library(dplyr)
library(tidyr)
library(tibble)
library(scales)
library(ggplot2)
library(ggcorrplot)
library(viridis)
library(purrr)
library(stringr)

# ------------------------------------------------------------------------------
# 1. DATA INGESTION & GLOBAL PREPARATION
# ------------------------------------------------------------------------------

find_data_path <- function() {
  candidates <- c(
    "outputs_new/merged_master_dataframe.csv",
    "TheraFlow/outputs_new/merged_master_dataframe.csv",
    "outputs/merged_master_dataframe.csv",
    "../outputs_new/merged_master_dataframe.csv",
    "/users/stgrad/rojashrj/Projects/omics2mol/TheraFlow/outputs_new/merged_master_dataframe.csv"
  )
  for (cand in candidates) {
    if (file.exists(cand)) return(cand)
  }
  return(candidates[1])
}


# Asset encoding for Figure 2a
fig2a_path <- file.path(normalizePath(getwd(), winslash = "/"), "TheraFlow/Figures/Figure2a.png")
if (!file.exists(fig2a_path)) fig2a_path <- file.path(normalizePath(getwd(), winslash = "/"), "Figures/Figure2a.png")
if (!file.exists(fig2a_path)) fig2a_path <- "Figures/Figure2a.png"

fig2a_uri <- if (file.exists(fig2a_path)) {
  paste0("data:image/png;base64,", base64enc::base64encode(fig2a_path))
} else ""

raw_data_path <- find_data_path()
if (!file.exists(raw_data_path)) {
  stop(paste("Data file not found at:", raw_data_path))
}

# Read raw dataset
master_df <- fread(raw_data_path)

# Method name standardization
method_map <- c(
  "biaee"       = "BiAAE",
  "biaae"       = "BiAAE",
  "bicev"       = "BiCEV",
  "fame"        = "FAME",
  "gex2sgen"    = "Gex2SGen",
  "gexmolgen"   = "GexMolGen",
  "gxrnn"       = "GxRNN",
  "paccmann_rl" = "Paccmann RL",
  "smilesgen"   = "SmilesGEN",
  "transgem"    = "TransGEM",
  "triomphe"    = "Triomphe"
)

master_df[, method := ifelse(tolower(method) %in% names(method_map), 
                             method_map[tolower(method)], method)]

# Dataset display name mapping
dataset_display_map <- c(
  "cardiomyopathy" = "Cardiomyopathy",
  "parkinsons"     = "Parkinson\x27s Disease",
  "hcc"            = "Hepatocellular Carcinoma",
  "glaucoma"       = "Glaucoma (Single-cell)",
  "glaucoma_bulk"  = "Glaucoma (Bulk)",
  "lincs1000_1"    = "Lincs dataset 1",
  "lincs1000_2"    = "Lincs dataset 2"
)

master_df[, dataset_display := ifelse(tolower(test_dataset) %in% names(dataset_display_map),
                                      dataset_display_map[tolower(test_dataset)], test_dataset)]

# Clean boolean columns
bool_cols <- c(
  "Passes_Safety", "mol_pred_loaded", "sanitization", "inchi_convertible",
  "all_atoms_connected", "no_radicals", "bond_lengths", "bond_angles",
  "internal_steric_clash", "aromatic_ring_flatness",
  "non-aromatic_ring_non-flatness", "double_bond_flatness", "internal_energy",
  "passes_valence_checks", "passes_kekulization", "is_top_consensus", "significant"
)

for (col in bool_cols) {
  if (col %in% names(master_df)) {
    master_df[[col]] <- as.logical(master_df[[col]] == TRUE | tolower(as.character(master_df[[col]])) %in% c("true", "1", "yes"))
  }
}

safe_div <- function(a, b) ifelse(is.na(b) | b == 0, 0, a / b)
clip01 <- function(x) pmin(pmax(x, 0), 1)
normalize_count <- function(x, cap) clip01(x / cap)

# Calculate engineered PoseBusters & Validity metrics if not present
if (!("structural_validity" %in% names(master_df))) {
  master_df[, structural_validity := rowMeans(cbind(
    as.numeric(inchi_convertible),
    as.numeric(all_atoms_connected),
    as.numeric(no_radicals),
    safe_div(number_valid_bonds, number_bonds),
    safe_div(number_valid_angles, number_angles),
    safe_div(number_valid_noncov_pairs, number_noncov_pairs)
  ), na.rm = TRUE)]
}
if (!("bond_angle_consistency" %in% names(master_df))) {
  master_df[, bond_angle_consistency := 1 - clip01(
    (normalize_count(number_short_outlier_bonds, 10) +
     normalize_count(number_long_outlier_bonds, 10) +
     normalize_count(number_outlier_angles, 10)) / 3
  )]
}
if (!("steric_feasibility" %in% names(master_df))) {
  master_df[, steric_feasibility := 1 - clip01(
    (normalize_count(number_clashes, 10) +
     normalize_count(number_noncov_pairs, 20)) / 2
  )]
}
if (!("double_bond_geometry" %in% names(master_df))) {
  master_df[, double_bond_geometry := rowMeans(cbind(
    safe_div(number_double_bonds_pass, number_double_bonds_checked),
    safe_div(double_bond_maximum_distance_from_plane, 1)
  ), na.rm = TRUE)]
}
if (!("weighted_sanitization" %in% names(master_df))) {
  master_df[, weighted_sanitization := rowMeans(cbind(
    as.numeric(sanitization),
    as.numeric(passes_valence_checks),
    as.numeric(passes_kekulization),
    as.numeric(no_radicals_before_sanitization)
  ), na.rm = TRUE)]
}

# Precalculate stage flags
master_df[, stage1_valid := !is.na(QED)]
master_df[, stage2_druglike := (!is.na(QED) & QED >= 0.5 & SA_Score <= 4.5 & MW <= 500 & LogP <= 5.0 & HBA <= 10 & HBD <= 5 & Passes_Safety == TRUE)]
master_df[, stage2_posebusters := (internal_steric_clash == TRUE & aromatic_ring_flatness == TRUE & 
                                   `non-aromatic_ring_non-flatness` == TRUE & double_bond_flatness == TRUE & internal_energy == TRUE)]
master_df[, stage3_reversal := !is.na(reversal_score_scaled)]
master_df[, stage4_consensus := (is_top_consensus == TRUE)]
master_df[, in_success_zone := (Ligand_RMSD_Pocket_Aligned < 2.0 & HOMO_LUMO_Gap_eV >= 0.05 & HOMO_LUMO_Gap_eV <= 0.20)]

all_num_cols <- names(master_df)[sapply(master_df, is.numeric)]
default_funky_metrics <- intersect(
  c("Validity", "Uniqueness", "Novelty", "IntDiv", "Passes_Safety", "QED", "SA_Score",
    "MW", "LogP", "TPSA", "HBA", "HBD", "Sim_Morg_Tani", "internal_energy",
    "aromatic_ring_flatness", "reversal_score_scaled", "Ranking_Score",
    "Ligand_RMSD_Pocket_Aligned", "HOMO_LUMO_Gap_eV", "Total_SCF_Energy_Eh"),
  c(all_num_cols, "Passes_Safety", "internal_energy", "aromatic_ring_flatness")
)
if (length(default_funky_metrics) < 5) default_funky_metrics <- head(all_num_cols, 15)

all_methods <- sort(unique(master_df$method))
all_datasets <- sort(unique(master_df$dataset_display))

# Palette for Generative Methods
method_colors <- c(
  "BiAAE"       = "#ff9999",
  "BiCEV"       = "#ffbb99",
  "FAME"        = "#ffecb3",
  "Gex2SGen"    = "#86efac",
  "GexMolGen"   = "#93c5fd",
  "GxRNN"       = "#a5b4fc",
  "Paccmann RL" = "#f472b6",
  "SmilesGEN"   = "#fbcfe8",
  "TransGEM"    = "#cbd5e1",
  "Triomphe"    = "#67e8f9"
)

# ------------------------------------------------------------------------------
# 2. UI DESIGN (bslib theme with modern dark scientific styling)
# ------------------------------------------------------------------------------

app_theme <- bs_theme(
  version = 5,
  bootswatch = "darkly",
  bg = "#0d131f",
  fg = "#e2e8f0",
  primary = "#38bdf8",
  secondary = "#64748b",
  success = "#10b981",
  info = "#818cf8",
  warning = "#f59e0b",
  danger = "#f43f5e",
  base_font = font_google("Inter"),
  code_font = font_google("JetBrains Mono")
)

ui <- page_navbar(
  id = "main_navbar",
  theme = app_theme,
  title = div(
    class = "d-flex align-items-center gap-2",
    tags$span(style = "font-size: 1.4rem; font-weight: 800; color: #38bdf8;", "TheraFlow"),
    tags$span(style = "font-size: 0.85rem; background: #1e293b; color: #94a3b8; padding: 2px 8px; border-radius: 6px; border: 1px solid #334155;", "Benchmark v1.0"),
    tags$span(style = "font-size: 0.85rem; color: #64748b; margin-left: 8px;", "A Multi-Stage Framework for evaluating generative drug models from transcriptomics and prioritizing candidates")
  ),
  
  # Inject Custom CSS
  header = tags$head(
    tags$style(HTML("
      body { background-color: #0b101b; font-family: \x27Inter\x27, sans-serif; }
      .navbar { border-bottom: 1px solid #1e293b; background-color: #0f172a !important; padding: 0.75rem 1.5rem; }
      .nav-link { font-weight: 500; font-size: 0.95rem; color: #94a3b8 !important; transition: all 0.2s ease; }
      .nav-link.active { color: #38bdf8 !important; border-bottom: 2px solid #38bdf8; font-weight: 600; }
      .card { background-color: #111827; border: 1px solid #1f2937; border-radius: 10px; box-shadow: 0 4px 6px -1px rgba(0,0,0,0.3); }
      .card-header { background-color: #1a2234; border-bottom: 1px solid #1f2937; font-weight: 600; color: #f8fafc; }
      .stat-card { background: linear-gradient(135deg, #131d31 0%, #1e293b 100%); border: 1px solid #283548; border-radius: 10px; padding: 15px; }
      .stat-val { font-size: 1.8rem; font-weight: 800; }
      .stat-lbl { font-size: 0.8rem; text-transform: uppercase; letter-spacing: 0.05em; color: #94a3b8; }
      .stat-sub { font-size: 0.75rem; color: #64748b; }
      .sidebar { background-color: #0f172a; border-right: 1px solid #1e293b; padding: 15px; }
      .btn-primary { background-color: #0284c7; border-color: #0284c7; }
      .btn-primary:hover { background-color: #0369a1; }
      .badge-custom { padding: 4px 8px; border-radius: 4px; font-size: 0.75rem; font-weight: 600; }
      .table-dark { background-color: #111827; }
    "))
  ),
  
  # ----------------------------------------------------------------------------
  # SIDEBAR (Global Filters)
  # ----------------------------------------------------------------------------
  sidebar = sidebar(
    id = "main_sidebar",
    title = div(icon("sliders-h"), " Benchmark Filters"),
    width = 300,
    open = "open",
    
    # Dataset Filter
    selectizeInput(
      "selected_datasets",
      label = "Disease / Context Datasets:",
      choices = all_datasets,
      selected = all_datasets,
      multiple = TRUE,
      options = list(plugins = list("remove_button"))
    ),
    
    # Method Selection
    div(
      class = "d-flex justify-content-between align-items-center mb-1",
      tags$label(class = "control-label mb-0", "Generative Methods:"),
      actionLink("toggle_all_methods", "Select All", style = "font-size: 0.75rem; color: #38bdf8;")
    ),
    checkboxGroupInput(
      "selected_methods",
      label = NULL,
      choices = all_methods,
      selected = all_methods
    ),
    
    hr(style = "border-color: #334155;"),
    
    # Stage Gate Filter
    radioButtons(
      "stage_filter",
      label = "Pipeline Progression Filter:",
      choices = c(
        "All Generated Molecules" = "all",
        "Stage 1: Syntax Valid Only" = "stage1",
        "Stage 2: Drug-Like & 3D Valid" = "stage2",
        "Stage 3: Reversal Prioritized" = "stage3",
        "Stage 4: Top Consensus Candidates" = "stage4"
      ),
      selected = "all"
    ),
    
    # Numeric Sliders
    sliderInput(
      "qed_filter",
      label = "Min QED Score:",
      min = 0.0, max = 1.0, value = 0.0, step = 0.05
    ),
    sliderInput(
      "sa_filter",
      label = "Max SA Score (Synthesizability):",
      min = 1.0, max = 10.0, value = 10.0, step = 0.5
    ),
    
    hr(style = "border-color: #334155;"),
    
    # Reset & Export
    div(
      class = "d-grid gap-2",
      actionButton("reset_filters", "Reset All Filters", icon = icon("undo"), class = "btn-outline-secondary btn-sm"),
      downloadButton("download_filtered_csv", "Export Filtered CSV", class = "btn-primary btn-sm")
    )
  ),
  
    # ----------------------------------------------------------------------------
  # LANDING PAGE: Overview & Framework Staging (Full-Width, No Sidebar)
  # ----------------------------------------------------------------------------
  nav_panel(
    title = div(icon("house"), " Overview"),
    value = "Overview",
    
    div(
      class = "container-fluid py-3",
      
      # KPI Stat Row
      div(
        class = "row g-3 mb-4",
        div(class = "col-md", div(class = "stat-card",
          div(class = "stat-lbl", "Total Molecules"),
          div(class = "stat-val text-info", textOutput("kpi_total_mols", inline = TRUE)),
          div(class = "stat-sub", "Across all generative models")
        )),
        div(class = "col-md", div(class = "stat-card",
          div(class = "stat-lbl", "Validity Rate"),
          div(class = "stat-val text-primary", textOutput("kpi_validity_rate", inline = TRUE)),
          div(class = "stat-sub", "Syntactically parsable SMILES")
        )),
        div(class = "col-md", div(class = "stat-card",
          div(class = "stat-lbl", "Drug-Like & 3D"),
          div(class = "stat-val text-success", textOutput("kpi_druglike_rate", inline = TRUE)),
          div(class = "stat-sub", "Ro5 + PoseBusters pass")
        )),
        div(class = "col-md", div(class = "stat-card",
          div(class = "stat-lbl", "Disease Reversal"),
          div(class = "stat-val text-warning", textOutput("kpi_reversal_count", inline = TRUE)),
          div(class = "stat-sub", "Transcriptomics evaluated")
        )),
        div(class = "col-md", div(class = "stat-card",
          div(class = "stat-lbl", "Docking and Molecular Dynamics"),
          div(class = "stat-val text-danger", textOutput("kpi_consensus_count", inline = TRUE)),
          div(class = "stat-sub", "AlphaFold3 & QM validated")
        ))
      ),
      
      # Staged Evaluation Framework Hierarchy Diagram (Figure 2a)
      card(
        card_header(
          div(
            class = "d-flex justify-content-between align-items-center",
            tags$span("TheraFlow Staged Evaluation Framework"),
            tags$span(class = "badge bg-primary", "5-Tier Pipeline Architecture")
          )
        ),
        card_body(
          class = "text-center p-3",
          style = "background-color: #0f172a; border-radius: 8px;",
          if (fig2a_uri != "") {
            tags$img(
              src = fig2a_uri,
              alt = "TheraFlow Staged Evaluation Framework",
              style = "max-width: 100%; max-height: 560px; object-fit: contain; border-radius: 8px; box-shadow: 0 4px 10px rgba(0,0,0,0.5);"
            )
          } else {
            tags$iframe(
              src = "Figure2a.pdf",
              width = "100%",
              height = "560px",
              style = "border: none; border-radius: 8px;"
            )
          }
        )
      )
    )
  ),
  
  # ----------------------------------------------------------------------------
  # TAB 1: BENCHMARK MATRIX & FUNKY HEATMAP
  # ----------------------------------------------------------------------------
  nav_panel(
    title = div(icon("table-cells-large"), " 1. Benchmark Matrix & Funky Heatmap"),
    value = "Funky_Heatmap",
    
    div(
      class = "container-fluid py-2",
      
      # Main Funky Heatmap Card with Dynamic Metric Selector
      card(
        card_header(
          div(
            class = "d-flex justify-content-between align-items-center flex-wrap gap-2",
            tags$span("Comprehensive Multi-Stage Benchmark Funky Heatmap"),
            div(
              class = "d-flex align-items-center gap-2",
              tags$label(class = "text-muted small mb-0", "Sort Rows:"),
              selectInput(
                "funky_sort_by",
                label = NULL,
                choices = c("Method Name" = "method", "Sort by QED" = "QED", "Sort by Validity" = "Validity", "Sort by DRS" = "reversal_score_scaled", "Sort by Ranking Score" = "Ranking_Score"),
                selected = "method",
                width = "180px"
              )
            )
          )
        ),
        card_body(
          div(
            class = "mb-3",
            div(
              class = "d-flex justify-content-between align-items-center mb-1",
              tags$label(class = "control-label fw-bold text-light", "Select Metrics to Benchmark in Funky Heatmap (80+ available):"),
              div(
                class = "btn-group btn-group-sm",
                actionButton("funky_select_default", "Default Set", class = "btn btn-outline-info btn-sm py-0 px-2", style = "font-size: 0.75rem;"),
                actionButton("funky_select_all", "Select All", class = "btn btn-outline-primary btn-sm py-0 px-2", style = "font-size: 0.75rem;"),
                actionButton("funky_select_clear", "Clear All", class = "btn btn-outline-secondary btn-sm py-0 px-2", style = "font-size: 0.75rem;")
              )
            ),
            selectizeInput(
              "funky_selected_metrics",
              label = NULL,
              choices = all_num_cols,
              selected = default_funky_metrics,
              multiple = TRUE,
              options = list(plugins = list("remove_button"), placeholder = "Choose metrics to display...")
            )
          ),
          p(class = "text-muted small", "Marker circle sizes and Plasma colors correspond to min-max normalized performance (0 to 1) per metric across models. Hover over any point to inspect exact raw means."),
          plotlyOutput("plot_funky_plotly", height = "580px")
        )
      )
    )
  ),
# ----------------------------------------------------------------------------
  # TAB 2: STAGE ATTRITION & QUALITY GATES (Figures 2, 3, 4, 5, 6, 7)
  # ----------------------------------------------------------------------------
  nav_panel(
    title = div(icon("filter"), " 2. Stage Attrition & Quality"),
    
    div(
      class = "row g-3",
      
      # Fig 2: Cumulative Attrition Bar Chart
      div(
        class = "col-lg-6",
        card(
          card_header("Pipeline Cumulative Retention Rate Across Stages (%)"),
          card_body(plotlyOutput("plot_fig2_attrition", height = "400px"))
        )
      ),
      
      # Fig 4 & 5: Multi-Baseline Physicochemical Retention (Figure 6 in Report.qmd)
      div(
        class = "col-lg-6",
        card(
          card_header(
            div(
              class = "d-flex justify-content-between align-items-center flex-wrap gap-2",
              tags$span(tags$b("Structural 3D Quality vs. 2D Drug-Likeness Pass Rates Across Stages")),
              tags$small(style = "color: #94a3b8;", "Total, Syntax-Valid & Med-Chem Baselines")
            )
          ),
          card_body(
            div(
              style = "background: #111827; border-radius: 8px; padding: 4px;",
              plotlyOutput("plot_fig4_scatter", height = "420px")
            )
          )
        )
      ),
      
      # Fig 6: Drug-like Funnel
      div(
        class = "col-lg-6",
        card(
          card_header("Multi-Gate Drug-Like Property Retention Funnel"),
          card_body(plotlyOutput("plot_fig6_funnel", height = "420px"))
        )
      ),
      
      # Fig 7: PoseBusters 3D Checks
      div(
        class = "col-lg-6",
        card(
          card_header("3D PoseBusters Conformation Geometric Checks"),
          card_body(plotlyOutput("plot_fig7_posebusters", height = "420px"))
        )
      )
    )
  ),
  
  # ----------------------------------------------------------------------------
  # TAB 3: GENERATIVE ARCHITECTURE & LATENT SPACE (Figures 8, 9, 10)
  # ----------------------------------------------------------------------------
  nav_panel(
    title = div(icon("dna"), " 3. Architecture Profiles"),
    
    div(
      class = "row g-3 mb-3",
      
      # Fig 8: Multi-Axis Radar Plot
      div(
        class = "col-12",
        card(
          card_header(
            div(
              class = "d-flex justify-content-between align-items-center flex-wrap gap-2",
              tags$span(tags$b("Figure 8: 5-Axis Multi-Dimensional Architecture Radar Profile")),
              div(
                class = "d-flex align-items-center gap-2",
                tags$span(style = "font-size: 0.85rem; color: #94a3b8; font-weight: 500;", "Display Mode:"),
                radioButtons(
                  "radar_mode_select",
                  label = NULL,
                  choices = c(
                    "Top Molecules (Pareto)" = "top",
                    "All Molecules" = "all"
                  ),
                  selected = "top",
                  inline = TRUE
                )
              )
            )
          ),
          card_body(plotlyOutput("plot_fig8_radar", height = "520px"))
        )
      )
    ),
    
    div(
      class = "row g-3",
      
      # Fig 9: Covariance Heatmap with selector for top molecules (Fig 11), all molecules (Fig 20), diff (Fig 21)
      div(
        class = "col-12",
        card(
          card_header(
            div(
              class = "d-flex justify-content-between align-items-center flex-wrap gap-2",
              tags$span(tags$b("Figure 9: Covariance & Property Trade-offs")),
              div(
                class = "d-flex align-items-center gap-2",
                tags$span(style = "font-size: 0.85rem; color: #94a3b8; font-weight: 500;", "Display Mode:"),
                radioButtons(
                  "cov_mode_select",
                  label = NULL,
                  choices = c(
                    "Top Molecules" = "top",
                    "All Molecules" = "all",
                    "Covariance Diff" = "diff"
                  ),
                  selected = "top",
                  inline = TRUE
                )
              )
            )
          ),
          card_body(
            div(
              style = "background: #111827; border-radius: 8px; padding: 6px;",
              plotOutput("plot_fig9_covariance", height = "750px")
            )
          )
        )
      )
    )
  ),
  
  # ----------------------------------------------------------------------------
  # TAB 4: DISEASE REVERSAL & CHEMBL TARGETS (Figures 15, 16, 17)
  # ----------------------------------------------------------------------------
  nav_panel(
    title = div(icon("bullseye"), " 4. Disease Reversal & Targets"),
    
    div(
      class = "row g-3",
      
      # Fig 15: DRS Distribution
      div(
        class = "col-lg-7",
        card(
          card_header("Transcriptomics Disease Reversal Score (DRS) Distribution"),
          card_body(plotlyOutput("plot_fig15_drs", height = "450px"))
        )
      ),
      
      # Targets & Mechanism Exploration
      div(
        class = "col-lg-5",
        card(
          card_header("Predicted ChEMBL Targets & Mechanism of Action"),
          card_body(plotlyOutput("plot_targets_breakdown", height = "450px"))
        )
      )
    )
  ),
  
  # ----------------------------------------------------------------------------
  # TAB 5: DOWNSTREAM BIOPHYSICS & QUANTUM STABILITY (Figures 11, 12, 13, 14)
  # ----------------------------------------------------------------------------
  nav_panel(
    title = div(icon("atom"), " 5. Downstream MD & Quantum"),
    
    div(
      class = "row g-3",
      
      # Fig 11: MD Distributions
      div(
        class = "col-lg-6",
        card(
          card_header(
            div(
              class = "d-flex justify-content-between align-items-center",
              tags$span("AlphaFold3 & Quantum Property Distributions"),
              selectInput(
                "md_metric_select",
                label = NULL,
                choices = c(
                  "Ranking Score" = "Ranking_Score",
                  "Ligand RMSD (Pocket Aligned)" = "Ligand_RMSD_Pocket_Aligned",
                  "HOMO-LUMO Gap (eV)" = "HOMO_LUMO_Gap_eV",
                  "iPTM Confidence" = "iPTM",
                  "Pocket Backbone RMSD" = "Pocket_Backbone_RMSD",
                  "Total SCF Energy (Eh)" = "Total_SCF_Energy_Eh"
                ),
                selected = "Ligand_RMSD_Pocket_Aligned",
                width = "240px"
              )
            )
          ),
          card_body(plotlyOutput("plot_fig11_distributions", height = "420px"))
        )
      ),
      
      # Fig 13: HOMO-LUMO vs Ligand RMSD (Success Zone)
      div(
        class = "col-lg-6",
        card(
          card_header("HOMO-LUMO Gap vs Ligand RMSD (Pocket Aligned) Success Zone"),
          card_body(plotlyOutput("plot_fig13_success_zone", height = "420px"))
        )
      ),
      
      # Fig 12: Success Proportions
      div(
        class = "col-lg-6",
        card(
          card_header("Downstream Biophysical Success Zone Proportions"),
          card_body(plotlyOutput("plot_fig12_success_rate", height = "380px"))
        )
      ),
      
      # Fig 14: Biophysical Correlation Matrix
      div(
        class = "col-lg-6",
        card(
          card_header("Quantum & MD Property Correlation Matrix"),
          card_body(plotlyOutput("plot_fig14_md_correlation", height = "380px"))
        )
      )
    )
  ),
  
  # ----------------------------------------------------------------------------
  # TAB 6: MOLECULE EXPLORER & CANDIDATE DATA TABLE
  # ----------------------------------------------------------------------------
  nav_panel(
    title = div(icon("flask"), " 6. Molecule Explorer"),
    
    card(
      card_header(
        div(
          class = "d-flex justify-content-between align-items-center flex-wrap gap-2",
          tags$span(tags$b("Complete Master Benchmark Dataset (All 109 Columns & 17,931 Rows)")),
          div(
            class = "d-flex align-items-center gap-2",
            tags$span(class = "badge bg-primary", paste0(format(nrow(master_df), big.mark = ","), " Molecules")),
            tags$span(class = "badge bg-info", paste0(ncol(master_df), " Features")),
            downloadButton("download_master_csv", "Export Full Dataset (CSV)", class = "btn-sm btn-outline-light")
          )
        )
      ),
      card_body(
        div(
          style = "overflow-x: auto;",
          DTOutput("molecule_table")
        )
      )
    )
  )
)

# ------------------------------------------------------------------------------
# 3. SERVER LOGIC
# ------------------------------------------------------------------------------

server <- function(input, output, session) {
  # Automatically toggle sidebar: hidden on Landing Page, open on all benchmark tabs
  observe({
    req(input$main_navbar)
    if (identical(input$main_navbar, "Overview")) {
      sidebar_toggle(id = "main_sidebar", open = FALSE)
    } else {
      sidebar_toggle(id = "main_sidebar", open = TRUE)
    }
  })

  
  # Reactive toggle for selecting all methods
  observeEvent(input$toggle_all_methods, {
    if (length(input$selected_methods) == length(all_methods)) {
      updateCheckboxGroupInput(session, "selected_methods", selected = character(0))
    } else {
      updateCheckboxGroupInput(session, "selected_methods", selected = all_methods)
    }
  })
  
  # Reset filters
  observeEvent(input$reset_filters, {
    updateSelectizeInput(session, "selected_datasets", selected = all_datasets)
    updateCheckboxGroupInput(session, "selected_methods", selected = all_methods)
    updateRadioButtons(session, "stage_filter", selected = "all")
    updateSliderInput(session, "qed_filter", value = 0.0)
    updateSliderInput(session, "sa_filter", value = 10.0)
  })
  
  # Filtered Data Reactive
  filtered_data <- reactive({
    df <- master_df
    
    # Dataset filter
    if (length(input$selected_datasets) > 0) {
      df <- df[dataset_display %in% input$selected_datasets]
    }
    
    # Method filter
    if (length(input$selected_methods) > 0) {
      df <- df[method %in% input$selected_methods]
    }
    
    # Stage filter
    if (input$stage_filter == "stage1") {
      df <- df[stage1_valid == TRUE]
    } else if (input$stage_filter == "stage2") {
      df <- df[stage2_druglike == TRUE]
    } else if (input$stage_filter == "stage3") {
      df <- df[stage3_reversal == TRUE]
    } else if (input$stage_filter == "stage4") {
      df <- df[stage4_consensus == TRUE]
    }
    
    # Numeric sliders
    if (!is.null(input$qed_filter) && input$qed_filter > 0) {
      df <- df[!is.na(QED) & QED >= input$qed_filter]
    }
    if (!is.null(input$sa_filter) && input$sa_filter < 10.0) {
      df <- df[!is.na(SA_Score) & SA_Score <= input$sa_filter]
    }
    
    return(df)
  })
  
  # ----------------------------------------------------------------------------
  # KPI Output Handlers
  # ----------------------------------------------------------------------------
  output$kpi_total_mols <- renderText({
    format(nrow(filtered_data()), big.mark = ",")
  })
  
  output$kpi_validity_rate <- renderText({
    df <- filtered_data()
    if (nrow(df) == 0) return("0.0%")
    rate <- mean(df$stage1_valid, na.rm = TRUE) * 100
    sprintf("%.1f%%", rate)
  })
  
  output$kpi_druglike_rate <- renderText({
    df <- filtered_data()
    if (nrow(df) == 0) return("0.0%")
    rate <- mean(df$stage2_druglike, na.rm = TRUE) * 100
    sprintf("%.1f%%", rate)
  })
  
  output$kpi_reversal_count <- renderText({
    df <- filtered_data()
    format(sum(df$stage3_reversal, na.rm = TRUE), big.mark = ",")
  })
  
  output$kpi_consensus_count <- renderText({
    df <- filtered_data()
    format(sum(df$stage4_consensus, na.rm = TRUE), big.mark = ",")
  })
  
  # ----------------------------------------------------------------------------
  # TAB 1: FUNKY HEATMAP RENDERING
  # ----------------------------------------------------------------------------
  
  # Compute Aggregated Summary Table for Funky Heatmap
    # Funky Heatmap Metric Selection Helpers
  observeEvent(input$funky_select_default, {
    updateSelectizeInput(session, "funky_selected_metrics", selected = default_funky_metrics)
  })
  observeEvent(input$funky_select_all, {
    updateSelectizeInput(session, "funky_selected_metrics", selected = all_num_cols)
  })
  observeEvent(input$funky_select_clear, {
    updateSelectizeInput(session, "funky_selected_metrics", selected = character(0))
  })

  # Dynamic Funky Heatmap Aggregation
  funky_summary_data <- reactive({
    df <- filtered_data()
    req(nrow(df) > 0, length(input$funky_selected_metrics) > 0)
    
    sel_metrics <- intersect(input$funky_selected_metrics, names(df))
    req(length(sel_metrics) > 0)
    
    # Compute mean for each selected metric per method
    summary_list <- list()
    for (m in sel_metrics) {
      vals <- df[, .(Metric = m, Mean_Raw = mean(as.numeric(get(m)), na.rm = TRUE)), by = method]
      summary_list[[m]] <- vals
    }
    agg_df <- rbindlist(summary_list)
    agg_df[is.nan(Mean_Raw) | is.na(Mean_Raw), Mean_Raw := 0]
    
    # Compute Min-Max Normalized Score per metric [0, 1]
    agg_df[, Norm_Score := {
      mn <- min(Mean_Raw, na.rm = TRUE)
      mx <- max(Mean_Raw, na.rm = TRUE)
      if (is.na(mn) || is.na(mx) || mn == mx) rep(0.5, .N) else (Mean_Raw - mn) / (mx - mn)
    }, by = Metric]
    
    # Row Sorting
    sort_var <- input$funky_sort_by
    if (!is.null(sort_var) && sort_var != "method" && sort_var %in% agg_df$Metric) {
      method_order <- agg_df[Metric == sort_var][order(-Mean_Raw)]$method
    } else {
      method_order <- sort(unique(agg_df$method))
    }
    
    list(data = agg_df, methods = method_order, metrics = sel_metrics)
  })
  
  # Plotly Funky Heatmap
  output$plot_funky_plotly <- renderPlotly({
    res <- funky_summary_data()
    if (is.null(res) || nrow(res$data) == 0) return(plotly_empty())
    
    plot_df <- as.data.frame(res$data)
    plot_df$Method <- factor(plot_df$method, levels = rev(res$methods))
    plot_df$Metric <- factor(plot_df$Metric, levels = res$metrics)
    
    plot_ly(
      data = plot_df,
      x = ~Metric,
      y = ~Method,
      type = "scatter",
      mode = "markers",
      marker = list(
        size = ~scales::rescale(Norm_Score, to = c(12, 28)),
        color = ~Norm_Score,
        colorscale = "Plasma",
        showscale = TRUE,
        colorbar = list(title = "Score", len = 0.6, thickness = 14, tickfont = list(color = "#cbd5e1")),
        line = list(color = "#ffffff", width = 1)
      ),
      hoverinfo = "text",
      text = ~paste0(
        "<b>Method:</b> ", Method, "<br>",
        "<b>Metric:</b> ", Metric, "<br>",
        "<b>Raw Mean:</b> ", ifelse(abs(Mean_Raw) >= 100, sprintf("%.1f", Mean_Raw), sprintf("%.4f", Mean_Raw)), "<br>",
        "<b>Normalized Score:</b> ", sprintf("%.2f", Norm_Score)
      )
    ) %>%
      layout(
        paper_bgcolor = "#111827",
        plot_bgcolor = "#111827",
        xaxis = list(
          title = "",
          tickangle = -40,
          tickfont = list(color = "#cbd5e1", size = 10),
          gridcolor = "#1f2937"
        ),
        yaxis = list(
          title = "",
          tickfont = list(color = "#f8fafc", size = 12),
          gridcolor = "#1f2937"
        ),
        margin = list(l = 100, r = 30, t = 20, b = 100)
      )
  })
  
# ----------------------------------------------------------------------------
  # TAB 2: STAGE ATTRITION & QUALITY GATES (Figs 2, 4, 5, 6, 7)
  # ----------------------------------------------------------------------------
  
  # Figure 2: Cumulative Pipeline Attrition
  output$plot_fig2_attrition <- renderPlotly({
    df <- filtered_data()
    if (nrow(df) == 0) return(plotly_empty())
    
    attr_df <- df %>%
      group_by(method) %>%
      summarise(
        Total = n(),
        Stage1 = sum(stage1_valid) / n() * 100,
        Stage2 = sum(stage2_druglike) / n() * 100,
        Stage3 = sum(stage3_reversal) / n() * 100,
        Stage4 = sum(stage4_consensus) / n() * 100,
        .groups = "drop"
      ) %>%
      pivot_longer(cols = c(Stage1, Stage2, Stage3, Stage4), names_to = "Stage", values_to = "Pass_Rate")
    
    stage_labels <- c(
      "Stage1" = "Stage 1: Syntax Valid",
      "Stage2" = "Stage 2: Drug-Like & 3D",
      "Stage3" = "Stage 3: Disease Reversal",
      "Stage4" = "Stage 4: MD Consensus"
    )
    attr_df$Stage_Label <- stage_labels[attr_df$Stage]
    
    p <- plot_ly(
      data = attr_df,
      x = ~method,
      y = ~Pass_Rate,
      color = ~Stage_Label,
      colors = c("#38bdf8", "#ec4899", "#f59e0b", "#10b981"),
      type = "bar"
    ) %>%
      layout(
        barmode = "group",
        paper_bgcolor = "#111827",
        plot_bgcolor = "#111827",
        font = list(color = "#f8fafc"),
        xaxis = list(title = "Generative Method", gridcolor = "#1f2937"),
        yaxis = list(title = "Cumulative Molecule Retention (%)", gridcolor = "#1f2937", range = c(0, 105)),
        legend = list(orientation = "h", y = 1.15, x = 0.1)
      )
    return(p)
  })
  
  # Figure 4 & 5: Structural 3D Quality vs. 2D Drug-Likeness Pass Rates (Report.qmd lines 710-745)
  output$plot_fig4_scatter <- renderPlotly({
    raw <- filtered_data()
    if (nrow(raw) == 0) return(plotly_empty())
    
    external_totals <- raw %>%
      group_by(method) %>%
      summarise(Total_Gen = n(), .groups = "drop")
    
    drug_like_analysis <- raw %>%
      group_by(method) %>%
      summarise(
        Stage1_Valid = sum(!is.na(QED)),
        Stage2_Valid = sum((MW < 500) & (LogP < 5) & (HBA <= 10) & (HBD <= 5), na.rm = TRUE),
        Count_druglike = sum((QED > 0.5) & (SA_Score < 4.5) & (MW < 500) & (LogP < 5) & (HBA <= 10) & (HBD <= 5), na.rm = TRUE),
        .groups = "drop"
      ) %>%
      pivot_longer(
        cols = starts_with("Count_"),
        names_to = "Metric",
        values_to = "Passed_Count"
      ) %>%
      left_join(external_totals, by = "method") %>%
      mutate(
        Metric = gsub("Count_", "", Metric),
        Pass_Rate_Total  = (Passed_Count / Total_Gen) * 100,
        Pass_Rate_Stage1 = (Passed_Count / Stage1_Valid) * 100,
        Pass_Rate_Stage2 = (Passed_Count / Stage2_Valid) * 100
      )
    
    posebusters_analysis <- raw %>%
      mutate(
        pass_Basic_Chem = sanitization & all_atoms_connected & passes_valence_checks & passes_kekulization,
        pass_Geometry   = bond_lengths & bond_angles & (number_clashes == 0),
        pass_Planarity  = aromatic_ring_flatness & double_bond_flatness,
        pass_PoseBusters_Total = pass_Basic_Chem & pass_Geometry & pass_Planarity
      ) %>%
      group_by(method) %>%
      summarise(
        Stage1_Valid = sum(!is.na(QED)),
        Stage2_Valid = sum((MW < 500) & (LogP < 5) & (HBA <= 10) & (HBD <= 5), na.rm = TRUE),
        Passed_Count   = sum(pass_PoseBusters_Total == TRUE, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      left_join(external_totals, by = "method") %>%
      mutate(
        Pass_Rate_Total  = (Passed_Count / Total_Gen) * 100,
        Pass_Rate_Stage1 = (Passed_Count / Stage1_Valid) * 100,
        Pass_Rate_Stage2 = (Passed_Count / Stage2_Valid) * 100
      )
    
    combined <- drug_like_analysis %>%
      dplyr::select(method, 
             DrugLike_Total  = Pass_Rate_Total, 
             DrugLike_Stage1 = Pass_Rate_Stage1, 
             DrugLike_Stage2 = Pass_Rate_Stage2,
             DrugLike_Count  = Passed_Count) %>%
      inner_join(
        posebusters_analysis %>%
          dplyr::select(method, 
                 Structural_Total  = Pass_Rate_Total, 
                 Structural_Stage1 = Pass_Rate_Stage1, 
                 Structural_Stage2 = Pass_Rate_Stage2,
                 Structural_Count  = Passed_Count),
        by = "method"
      )
    
    combined_long <- combined %>%
      pivot_longer(
        cols = c(DrugLike_Total, DrugLike_Stage1, DrugLike_Stage2,
                 Structural_Total, Structural_Stage1, Structural_Stage2),
        names_to = c("Metric", "Stage"),
        names_pattern = "(DrugLike|Structural)_(.*)",
        values_to = "Pass_Rate"
      ) %>%
      pivot_wider(
        names_from = Metric,
        values_from = Pass_Rate
      ) %>%
      mutate(Stage = factor(Stage, levels = c("Total", "Stage1", "Stage2"),
                             labels = c("vs Total Generated", "vs Stage 1 Syntax-Valid", "vs Stage 2 Med-Chem-Valid"))) %>%
      mutate(
        hover_text = paste0(
          "<b>Method:</b> ", method,
          "<br><b>Stage:</b> ", Stage,
          "<br><b>Drug-Like Pass Rate:</b> ", round(DrugLike, 1), "%",
          "<br><b>Structural (PB) Pass Rate:</b> ", round(Structural, 1), "%",
          "<br><b>Drug-Like Count:</b> ", DrugLike_Count
        )
      )
    
    p <- ggplot(combined_long, aes(x = DrugLike, y = Structural, fill = method, size = DrugLike_Count, text = hover_text)) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#f43f5e", linewidth = 0.8) +
      geom_point(shape = 21, color = "black", stroke = 0.8, alpha = 0.95) +
      facet_wrap(~ Stage, nrow = 1) +
      scale_fill_manual(values = method_colors) +
      scale_size_continuous(range = c(3, 8), name = "Count") +
      xlim(0, 100) + ylim(0, 100) +
      theme_classic(base_size = 12) +
      labs(
        x = "Drug-like Pass Rate (%)",
        y = "Structural (PoseBusters) Pass Rate (%)",
        fill = "Generative Method"
      ) +
      theme(
        strip.text = element_text(face = "bold", color = "#f8fafc", size = 11),
        strip.background = element_rect(fill = "#1e293b", color = "#334155"),
        axis.text = element_text(size = 9, color = "#cbd5e1"),
        axis.title = element_text(face = "bold", color = "#f8fafc", size = 11),
        axis.line = element_line(color = "#475569"),
        axis.ticks = element_line(color = "#475569"),
        panel.grid.major = element_line(color = "#1e293b", linewidth = 0.4),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        plot.background = element_rect(fill = "#111827", color = NA),
        panel.background = element_rect(fill = "#0f172a", color = NA)
      )
    
    ply <- ggplotly(p, tooltip = "text") %>%
      layout(
        paper_bgcolor = "#111827",
        plot_bgcolor = "#0f172a",
        font = list(color = "#f8fafc"),
        showlegend = FALSE
      )
    
    for (i in seq_along(ply$x$layout$annotations)) {
      ply$x$layout$annotations[[i]]$font$color <- "#f8fafc"
    }
    
    return(ply)
  })
  
  # Figure 6: Drug-like Multi-Gate Retention Funnel
  output$plot_fig6_funnel <- renderPlotly({
    df <- filtered_data()
    if (nrow(df) == 0) return(plotly_empty())
    
    funnel_agg <- df %>%
      group_by(method) %>%
      summarise(
        `Total Generated` = n(),
        `Syntax Valid`    = sum(stage1_valid),
        `MW [200-500]`    = sum(MW >= 200 & MW <= 500, na.rm = TRUE),
        `LogP <= 5.0`     = sum(LogP <= 5.0, na.rm = TRUE),
        `HBA <= 10`       = sum(HBA <= 10, na.rm = TRUE),
        `HBD <= 5`        = sum(HBD <= 5, na.rm = TRUE),
        `QED >= 0.5`      = sum(QED >= 0.5, na.rm = TRUE),
        `SA Score <= 4.5` = sum(SA_Score <= 4.5, na.rm = TRUE),
        `Safety Passed`   = sum(Passes_Safety == TRUE, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      pivot_longer(cols = -method, names_to = "Gate", values_to = "Count")
    
    funnel_agg$Gate <- factor(funnel_agg$Gate, levels = c(
      "Total Generated", "Syntax Valid", "MW [200-500]", "LogP <= 5.0",
      "HBA <= 10", "HBD <= 5", "QED >= 0.5", "SA Score <= 4.5", "Safety Passed"
    ))
    
    p <- plot_ly(
      data = funnel_agg,
      x = ~Gate,
      y = ~Count,
      color = ~method,
      colors = method_colors,
      type = "scatter",
      mode = "lines+markers"
    ) %>%
      layout(
        paper_bgcolor = "#111827",
        plot_bgcolor = "#111827",
        font = list(color = "#f8fafc"),
        xaxis = list(title = "Chemical Filtration Gate", tickangle = -30, gridcolor = "#1f2937"),
        yaxis = list(title = "Molecule Count Remaining", gridcolor = "#1f2937"),
        legend = list(orientation = "h", y = 1.15)
      )
    return(p)
  })
  
  # Figure 7: PoseBusters 3D Checks Matrix
  output$plot_fig7_posebusters <- renderPlotly({
    df <- filtered_data()
    if (nrow(df) == 0) return(plotly_empty())
    
    pose_agg <- df %>%
      group_by(method) %>%
      summarise(
        `Bond Lengths` = mean(as.numeric(bond_lengths), na.rm = TRUE) * 100,
        `Bond Angles` = mean(as.numeric(bond_angles), na.rm = TRUE) * 100,
        `No Steric Clash` = mean(as.numeric(internal_steric_clash), na.rm = TRUE) * 100,
        `Aromatic Flatness` = mean(as.numeric(aromatic_ring_flatness), na.rm = TRUE) * 100,
        `Non-Arom Non-Flatness` = mean(as.numeric(`non-aromatic_ring_non-flatness`), na.rm = TRUE) * 100,
        `Double Bond Flatness` = mean(as.numeric(double_bond_flatness), na.rm = TRUE) * 100,
        `Internal Energy` = mean(as.numeric(internal_energy), na.rm = TRUE) * 100,
        .groups = "drop"
      ) %>%
      pivot_longer(cols = -method, names_to = "Check", values_to = "Pass_Rate")
    
    p <- plot_ly(
      data = pose_agg,
      x = ~Check,
      y = ~Pass_Rate,
      color = ~method,
      colors = method_colors,
      type = "bar"
    ) %>%
      layout(
        barmode = "group",
        paper_bgcolor = "#111827",
        plot_bgcolor = "#111827",
        font = list(color = "#f8fafc"),
        xaxis = list(title = "PoseBusters 3D Quality Check", tickangle = -30, gridcolor = "#1f2937"),
        yaxis = list(title = "Conformation Pass Rate (%)", gridcolor = "#1f2937", range = c(0, 105)),
        legend = list(orientation = "h", y = 1.15)
      )
    return(p)
  })
  
  # ----------------------------------------------------------------------------
  # TAB 3: GENERATIVE ARCHITECTURE PROFILES (Figs 8, 9, 10)
  # ----------------------------------------------------------------------------
  
  # Figure 8: 5-Axis Multi-Dimensional Architecture Radar Profile (Report.qmd p1 & p2)
  output$plot_fig8_radar <- renderPlotly({
    df <- filtered_data()
    if (nrow(df) == 0) return(plotly_empty())
    
    mode <- input$radar_mode_select
    if (is.null(mode)) mode <- "top"
    
    target_data <- if (mode == "top") {
      sub <- df[is_top_consensus == TRUE]
      if (nrow(sub) < 5) df else sub
    } else {
      df
    }
    
    safe_rescale <- function(x) {
      if (all(is.na(x))) return(rep(0, length(x)))
      r <- range(x, na.rm = TRUE)
      if (r[1] == r[2]) return(rep(0.5, length(x)))
      scales::rescale(x, to = c(0, 1))
    }
    
    # Construct 5 macro-pillars as engineered in Report.qmd lines 1090-1230
    radar_df <- target_data %>%
      mutate(
        norm_QED = safe_rescale(QED),
        norm_SA = safe_rescale(10 - SA_Score),
        norm_MW = safe_rescale(MW),
        norm_LogP = safe_rescale(LogP),
        norm_TPSA = safe_rescale(TPSA),
        norm_HBA = safe_rescale(HBA),
        norm_HBD = safe_rescale(HBD),
        norm_F_Csp3 = safe_rescale(F_Csp3),
        norm_structural_validity = safe_rescale(structural_validity),
        norm_bond_consistency = safe_rescale(bond_angle_consistency),
        norm_steric_feasibility = safe_rescale(steric_feasibility),
        norm_double_bond = safe_rescale(double_bond_geometry),
        norm_sanitization = safe_rescale(weighted_sanitization),
        norm_similarity = if ("Sim_Morg_Tani" %in% names(target_data)) Sim_Morg_Tani else 0,
        norm_energy_ratio = safe_rescale(energy_ratio),
        norm_energy_stability = safe_rescale(-1 * ensemble_avg_energy)
      ) %>%
      mutate(
        `Drug-likeness` = rowMeans(cbind(
          norm_QED, norm_SA, norm_MW, norm_LogP, norm_TPSA, norm_HBA, norm_HBD, norm_F_Csp3
        ), na.rm = TRUE),
        `Structural Properties` = rowMeans(cbind(
          norm_structural_validity, norm_bond_consistency, norm_steric_feasibility,
          norm_double_bond, norm_sanitization
        ), na.rm = TRUE),
        `Energy (MD)` = rowMeans(cbind(
          norm_energy_ratio, norm_energy_stability
        ), na.rm = TRUE)
      ) %>%
      group_by(method) %>%
      summarise(
        `Drug-likeness` = mean(`Drug-likeness`, na.rm = TRUE),
        `Structural Properties` = mean(`Structural Properties`, na.rm = TRUE),
        `Energy (MD)` = mean(`Energy (MD)`, na.rm = TRUE),
        `Disease Reversal (DRS)` = mean(reversal_score_scaled, na.rm = TRUE),
        `Similarity` = mean(norm_similarity, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(across(-method, safe_rescale))
    
    categories <- c("Drug-likeness", "Structural Properties", "Energy (MD)", "Disease Reversal (DRS)", "Similarity")
    
    p <- plot_ly(type = "scatterpolar", fill = "toself", mode = "lines+markers")
    
    for (i in 1:nrow(radar_df)) {
      m <- radar_df$method[i]
      vals <- c(
        radar_df$`Drug-likeness`[i],
        radar_df$`Structural Properties`[i],
        radar_df$`Energy (MD)`[i],
        radar_df$`Disease Reversal (DRS)`[i],
        radar_df$`Similarity`[i]
      )
      vals <- c(vals, vals[1]) # Close loop
      cats <- c(categories, categories[1])
      
      col <- if (m %in% names(method_colors)) method_colors[m] else "#38bdf8"
      
      p <- p %>% add_trace(
        r = vals,
        theta = cats,
        name = m,
        line = list(color = col, width = 2),
        marker = list(size = 6),
        fillcolor = paste0(col, "33")
      )
    }
    
    subtitle_txt <- if (mode == "top") "Top Pareto-Optimal Molecules" else "All Generated Molecules"
    
    p <- p %>% layout(
      polar = list(
        bgcolor = "#111827",
        radialaxis = list(visible = TRUE, range = c(0, 1), gridcolor = "#1f2937", linecolor = "#334155", tickfont = list(color = "#94a3b8", size = 10)),
        angularaxis = list(gridcolor = "#1f2937", linecolor = "#334155", tickfont = list(color = "#f8fafc", size = 11))
      ),
      paper_bgcolor = "#111827",
      font = list(color = "#f8fafc"),
      legend = list(orientation = "h", y = -0.15, x = 0.5, xanchor = "center")
    )
    return(p)
  })
  
  # Figure 9: Covariance & Property Trade-offs (Report Figures 11, 20, 21)
  output$plot_fig9_covariance <- renderPlot({
    df <- filtered_data()
    if (nrow(df) == 0) {
      return(ggplot() + theme_void() + annotate("text", x = 0.5, y = 0.5, label = "No data available matching filters", color = "#94a3b8", size = 6))
    }
    
    mode <- input$cov_mode_select
    if (is.null(mode)) mode <- "top"
    
    safe_rescale <- function(x) {
      if (all(is.na(x))) return(rep(0, length(x)))
      r <- range(x, na.rm = TRUE)
      if (r[1] == r[2]) return(rep(0.5, length(x)))
      scales::rescale(x, to = c(0, 1))
    }
    
    # 1. Normalize individual metrics across the dataset (Report.qmd lines 1255-1275)
    cov_prep <- df %>%
      mutate(
        norm_QED = safe_rescale(QED),
        norm_SA = safe_rescale(10 - SA_Score),
        norm_MW = safe_rescale(MW),
        norm_LogP = safe_rescale(LogP),
        norm_TPSA = safe_rescale(TPSA),
        norm_HBA = safe_rescale(HBA),
        norm_HBD = safe_rescale(HBD),
        norm_F_Csp3 = safe_rescale(F_Csp3),
        norm_structural_validity = safe_rescale(structural_validity),
        norm_bond_consistency = safe_rescale(bond_angle_consistency),
        norm_steric_feasibility = safe_rescale(steric_feasibility),
        norm_double_bond = safe_rescale(double_bond_geometry),
        norm_sanitization = safe_rescale(weighted_sanitization),
        norm_similarity = if ("Sim_Morg_Tani" %in% names(df)) Sim_Morg_Tani else 0,
        norm_energy_ratio = safe_rescale(energy_ratio),
        norm_energy_stability = safe_rescale(-1 * ensemble_avg_energy)
      )
    
    individual_metrics <- c(
      "norm_QED", "norm_SA", "norm_MW", "norm_LogP", "norm_TPSA", "norm_HBA", "norm_HBD", "norm_F_Csp3",
      "norm_structural_validity", "norm_bond_consistency", "norm_steric_feasibility", "norm_double_bond", "norm_sanitization",
      "norm_energy_ratio", "norm_energy_stability", "reversal_score_scaled", "norm_similarity"
    )
    
    # Helper to compute lower-triangle covariance per method
    calc_cov_lower <- function(data_sub) {
      cov_list <- data_sub %>%
        split(by = "method") %>%
        map(function(d) {
          cov(d[, ..individual_metrics], use = "pairwise.complete.obs")
        })
      
      purrr::imap_dfr(cov_list, function(cov_mat, arch_name) {
        as.data.frame(as.table(cov_mat)) %>%
          rename(class1 = Var1, class2 = Var2, abs_cov = Freq) %>%
          mutate(method = arch_name)
      }) %>%
        mutate(across(c(class1, class2), ~ . %>%
          stringr::str_remove_all("norm_") %>%
          stringr::str_replace_all("_", " ") %>%
          tools::toTitleCase())) %>%
        filter(as.character(class1) >= as.character(class2))
    }
    
    if (mode == "top") {
      cov_top_dat <- cov_prep[is_top_consensus == TRUE]
      if (nrow(cov_top_dat) < 5) cov_top_dat <- cov_prep
      plot_df <- calc_cov_lower(cov_top_dat)
      plot_df$fill_val <- plot_df$abs_cov
      fill_title <- "Absolute\nCovariance"
      plot_title <- "Covariance of Metrics on Top Molecules"
    } else if (mode == "all") {
      plot_df <- calc_cov_lower(cov_prep)
      plot_df$fill_val <- plot_df$abs_cov
      fill_title <- "Covariance"
      plot_title <- "Covariance of Metrics on All Molecules"
    } else {
      cov_top_dat <- cov_prep[is_top_consensus == TRUE]
      if (nrow(cov_top_dat) < 5) cov_top_dat <- cov_prep
      plot_top <- calc_cov_lower(cov_top_dat)
      plot_all <- calc_cov_lower(cov_prep)
      plot_df <- plot_top
      plot_df$fill_val <- plot_top$abs_cov - plot_all$abs_cov
      fill_title <- "Covariance\nDifference"
      plot_title <- "Covariance Difference: Top vs All Molecules"
    }
    
    finite_vals <- plot_df$fill_val[is.finite(plot_df$fill_val)]
    if (length(finite_vals) == 0) finite_vals <- c(-0.1, 0.1)
    cov_range <- quantile(finite_vals, probs = c(0.01, 0.99), na.rm = TRUE)
    if (cov_range[1] == cov_range[2]) {
      cov_range[1] <- cov_range[1] - 0.01
      cov_range[2] <- cov_range[2] + 0.01
    }
    dynamic_breaks <- seq(from = cov_range[1], to = cov_range[2], length.out = 6)
    
    p <- ggplot(plot_df, aes(x = class1, y = class2, fill = fill_val)) +
      geom_tile(color = "#1e293b", linewidth = 0.25) +
      facet_wrap(~method, nrow = 2) +
      scale_fill_stepsn(
        colors = c("#e0f3f8", "#abd9e9", "#bfd3e6", "#9ebcda", "#8c96c6", "#8c6bb1", "#810f7c"),
        breaks = dynamic_breaks,
        labels = round(dynamic_breaks, 3),
        limits = cov_range,
        na.value = "#1e293b"
      ) +
      theme_minimal(base_family = "sans") +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8, color = "#cbd5e1", face = "bold"),
        axis.text.y = element_text(size = 8, color = "#cbd5e1", face = "bold"),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "#1e293b", color = "#334155"),
        strip.text = element_text(face = "bold", size = 11, color = "#f8fafc"),
        legend.position = "right",
        legend.title = element_text(size = 10, face = "bold", color = "#f8fafc"),
        legend.text = element_text(size = 9, color = "#cbd5e1"),
        plot.title = element_text(face = "bold", size = 13, color = "#38bdf8", hjust = 0, margin = margin(b = 10)),
        plot.background = element_rect(fill = "#111827", color = NA),
        panel.background = element_rect(fill = "#0f172a", color = NA),
        panel.border = element_rect(color = "#334155", fill = NA, linewidth = 0.5)
      ) +
      labs(
        title = plot_title,
        fill = fill_title,
        x = NULL, y = NULL
      ) +
      coord_fixed()
    
    return(p)
  }, bg = "transparent", res = 110)
  
  # ----------------------------------------------------------------------------
  # TAB 4: DISEASE REVERSAL & CHEMBL TARGETS (Figs 15, 16, 17)
  # ----------------------------------------------------------------------------
  
  # Figure 15: DRS Distribution Plot
  output$plot_fig15_drs <- renderPlotly({
    df <- filtered_data()
    drs_df <- df[!is.na(reversal_score_scaled)]
    if (nrow(drs_df) == 0) return(plotly_empty())
    
    p <- plot_ly(
      data = drs_df,
      x = ~method,
      y = ~reversal_score_scaled,
      color = ~confidence,
      colors = c("high" = "#f43f5e", "low" = "#38bdf8"),
      type = "box",
      boxpoints = "all",
      jitter = 0.3,
      pointpos = -1.8
    ) %>%
      layout(
        paper_bgcolor = "#111827",
        plot_bgcolor = "#111827",
        font = list(color = "#f8fafc"),
        xaxis = list(title = "Generative Method", gridcolor = "#1f2937"),
        yaxis = list(title = "Scaled Disease Reversal Score (DRS)", gridcolor = "#1f2937"),
        legend = list(title = list(text = "Confidence"), orientation = "h", y = 1.15)
      )
    return(p)
  })
  
  # Target Breakdown Plot
  output$plot_targets_breakdown <- renderPlotly({
    df <- filtered_data()
    targ_df <- df[!is.na(action) & action != ""]
    if (nrow(targ_df) == 0) {
      targ_df <- df[!is.na(Predicted_Targets)]
      if (nrow(targ_df) == 0) return(plotly_empty())
    }
    
    agg <- targ_df %>%
      group_by(method, action) %>%
      summarise(Count = n(), .groups = "drop")
    
    p <- plot_ly(
      data = agg,
      x = ~method,
      y = ~Count,
      color = ~action,
      type = "bar"
    ) %>%
      layout(
        barmode = "stack",
        paper_bgcolor = "#111827",
        plot_bgcolor = "#111827",
        font = list(color = "#f8fafc"),
        xaxis = list(title = "Generative Method", gridcolor = "#1f2937"),
        yaxis = list(title = "Targeted Molecule Count", gridcolor = "#1f2937"),
        legend = list(title = list(text = "Action Mechanism"), orientation = "h", y = 1.15)
      )
    return(p)
  })
  
  # ----------------------------------------------------------------------------
  # TAB 5: DOWNSTREAM BIOPHYSICS & QUANTUM STABILITY (Figs 11, 12, 13, 14)
  # ----------------------------------------------------------------------------
  
  # Figure 11: MD Distributions
  output$plot_fig11_distributions <- renderPlotly({
    df <- filtered_data()
    metric <- input$md_metric_select
    if (nrow(df) == 0 || !(metric %in% names(df))) return(plotly_empty())
    
    plot_sub <- df[!is.na(df[[metric]])]
    if (nrow(plot_sub) == 0) return(plotly_empty())
    
    p <- plot_ly(
      data = plot_sub,
      x = ~method,
      y = as.formula(paste0("~", metric)),
      color = ~method,
      colors = method_colors,
      type = "box",
      boxpoints = "all",
      jitter = 0.3,
      pointpos = -1.8
    ) %>%
      layout(
        paper_bgcolor = "#111827",
        plot_bgcolor = "#111827",
        font = list(color = "#f8fafc"),
        xaxis = list(title = "Generative Method", gridcolor = "#1f2937"),
        yaxis = list(title = metric, gridcolor = "#1f2937")
      )
    return(p)
  })
  
  # Figure 13: HOMO-LUMO Gap vs Ligand RMSD (Success Zone)
  output$plot_fig13_success_zone <- renderPlotly({
    df <- filtered_data()
    plot_sub <- df[!is.na(Ligand_RMSD_Pocket_Aligned) & !is.na(HOMO_LUMO_Gap_eV)]
    if (nrow(plot_sub) == 0) return(plotly_empty())
    
    p <- plot_ly() %>%
      # Success Zone Shading
      add_polygons(
        x = c(0, 2.0, 2.0, 0),
        y = c(0.05, 0.05, 0.20, 0.20),
        fillcolor = "rgba(16, 185, 129, 0.15)",
        line = list(color = "#10b981", dash = "dash", width = 1.5),
        name = "Biophysical Success Zone",
        hoverinfo = "none"
      ) %>%
      # Data Points
      add_trace(
        data = plot_sub,
        x = ~Ligand_RMSD_Pocket_Aligned,
        y = ~HOMO_LUMO_Gap_eV,
        color = ~Ranking_Score,
        colorscale = "Plasma",
        type = "scatter",
        mode = "markers",
        marker = list(size = 9, line = list(color = "#ffffff", width = 1)),
        text = ~paste0(
          "<b>Method:</b> ", method, "<br>",
          "<b>Dataset:</b> ", dataset_display, "<br>",
          "<b>Pocket RMSD:</b> ", sprintf("%.2f Å", Ligand_RMSD_Pocket_Aligned), "<br>",
          "<b>HOMO-LUMO Gap:</b> ", sprintf("%.3f eV", HOMO_LUMO_Gap_eV), "<br>",
          "<b>AlphaFold3 Rank:</b> ", sprintf("%.3f", Ranking_Score)
        ),
        hoverinfo = "text"
      ) %>%
      layout(
        paper_bgcolor = "#111827",
        plot_bgcolor = "#111827",
        font = list(color = "#f8fafc"),
        xaxis = list(title = "Ligand RMSD (Pocket Aligned) (Å)", range = c(0, 6), gridcolor = "#1f2937"),
        yaxis = list(title = "HOMO-LUMO Gap (eV)", range = c(0, 0.5), gridcolor = "#1f2937")
      )
    return(p)
  })
  
  # Figure 12: Success Proportions
  output$plot_fig12_success_rate <- renderPlotly({
    df <- filtered_data()
    plot_sub <- df[!is.na(in_success_zone)]
    if (nrow(plot_sub) == 0) return(plotly_empty())
    
    succ_agg <- plot_sub %>%
      group_by(method, in_success_zone) %>%
      summarise(Count = n(), .groups = "drop") %>%
      group_by(method) %>%
      mutate(Proportion = Count / sum(Count) * 100)
    
    succ_agg$Status <- ifelse(succ_agg$in_success_zone, "Success Zone", "Outside Zone")
    
    p <- plot_ly(
      data = succ_agg,
      x = ~method,
      y = ~Proportion,
      color = ~Status,
      colors = c("Success Zone" = "#10b981", "Outside Zone" = "#ef4444"),
      type = "bar"
    ) %>%
      layout(
        barmode = "stack",
        paper_bgcolor = "#111827",
        plot_bgcolor = "#111827",
        font = list(color = "#f8fafc"),
        xaxis = list(title = "Generative Method", gridcolor = "#1f2937"),
        yaxis = list(title = "Proportion (%)", range = c(0, 100), gridcolor = "#1f2937"),
        legend = list(orientation = "h", y = 1.15)
      )
    return(p)
  })
  
  # Figure 14: Downstream Biophysical Correlation Matrix
  output$plot_fig14_md_correlation <- renderPlotly({
    df <- filtered_data()
    md_cols <- c("Ranking_Score", "Ligand_RMSD_Pocket_Aligned", "HOMO_LUMO_Gap_eV", "Total_SCF_Energy_Eh", "iPTM", "Pocket_Backbone_RMSD", "MW")
    md_cols <- intersect(md_cols, names(df))
    
    sub_df <- na.omit(df[, ..md_cols])
    if (nrow(sub_df) < 5) return(plotly_empty())
    
    cor_mat <- cor(sub_df, use = "pairwise.complete.obs")
    
    p <- plot_ly(
      x = colnames(cor_mat),
      y = rownames(cor_mat),
      z = cor_mat,
      type = "heatmap",
      colorscale = "Viridis",
      zmin = -1, zmax = 1
    ) %>%
      layout(
        paper_bgcolor = "#111827",
        plot_bgcolor = "#111827",
        font = list(color = "#f8fafc"),
        xaxis = list(tickangle = -30, gridcolor = "#1f2937"),
        yaxis = list(gridcolor = "#1f2937")
      )
    return(p)
  })
  
  # ----------------------------------------------------------------------------
  # TAB 6: MOLECULE DATA TABLE
  # ----------------------------------------------------------------------------
  output$molecule_table <- renderDT({
    # Full dataset with all columns and all rows
    df <- master_df
    num_cols <- names(df)[sapply(df, is.numeric)]
    
    dt <- datatable(
      df,
      options = list(
        pageLength = 25,
        lengthMenu = list(c(10, 25, 50, 100, 500, -1), c("10", "25", "50", "100", "500", "All")),
        scrollX = TRUE,
        dom = "Bfrtip",
        autoWidth = TRUE
      ),
      class = "cell-border stripe hover table-dark",
      selection = "single",
      rownames = FALSE
    )
    
    if (length(num_cols) > 0) {
      dt <- dt %>% formatRound(columns = num_cols, digits = 3)
    }
    return(dt)
  })
  
  # Download handler for master CSV
  output$download_master_csv <- downloadHandler(
    filename = function() {
      paste0("theraflow_master_molecules_", Sys.Date(), ".csv")
    },
    content = function(file) {
      fwrite(master_df, file)
    }
  )
  
  # Download handler for filtered CSV
  output$download_filtered_csv <- downloadHandler(
    filename = function() {
      paste0("theraflow_filtered_molecules_", Sys.Date(), ".csv")
    },
    content = function(file) {
      fwrite(filtered_data(), file)
    }
  )
}

# ----------------------------------------------------------------------------
# 4. LAUNCH APPLICATION
# ----------------------------------------------------------------------------
shinyApp(ui = ui, server = server)
