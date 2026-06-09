## =============================================================================
## DeepARG_2.R
## CRASH — Cancer-associated gut micRobiome Analysis of Strains in
##         Haematological patients
##
## Purpose:
##   Cross-cohort overlap analysis of deepARG differential abundance results.
##   Identifies ARG types/subtypes that replicate directionally across multiple
##   cohorts. Set LEVEL to match DeepARG_1.R and DeepARG_3.R.
##
## Usage:
##   Set LEVEL to "type" or "subtype" (must match Scripts 1 and 3).
##
## Input:
##   cross_method_deeparg_<LEVEL>_<DATASET>.tsv — per cohort, from DeepARG_1.R
##
## Outputs (written to cross_dataset_deeparg/):
##   cross_dataset_overlap_deeparg_<LEVEL>.tsv
##   summary_table_deeparg_<LEVEL>.tsv
##
## Dependencies:
##   tidyverse
##
## Notes:
##   - Run DeepARG_1.R on all four cohorts before this script
##   - High-confidence: MaAsLin2 sig + ≥1 method directional consensus
## =============================================================================
library(tidyverse)

# ---- Configuration -------------------------------------------------------

Q_THRESH <- 0.25
LEVEL    <- "subtype"   # "type" or "subtype"

datasets <- list(
  prjna813705      = file.path("D:/MasterThesis/Vik/DeepARG", paste0("prjna813705_",      LEVEL), "tables", paste0("cross_method_deeparg_", LEVEL, "_prjna813705.tsv")),
  kulecka_leukemia = file.path("D:/MasterThesis/Vik/DeepARG", paste0("kulecka_leukemia_", LEVEL), "tables", paste0("cross_method_deeparg_", LEVEL, "_kulecka_leukemia.tsv")),
  kulecka_lymphoma = file.path("D:/MasterThesis/Vik/DeepARG", paste0("kulecka_lymphoma_", LEVEL), "tables", paste0("cross_method_deeparg_", LEVEL, "_kulecka_lymphoma.tsv")),
  cra              = file.path("D:/MasterThesis/Vik/DeepARG", paste0("cra_",              LEVEL), "tables", paste0("cross_method_deeparg_", LEVEL, "_cra.tsv"))
)

cancer_labels <- c(
  prjna813705      = "AML (PRJNA813705)",
  kulecka_leukemia = "AML (Kulecka)",
  kulecka_lymphoma = "LN (Kulecka)",
  cra              = "NKTCL (CRA)"
)

out_dir <- file.path("D:/MasterThesis/Vik/DeepARG", paste0("cross_dataset_", LEVEL))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Load high-confidence hits per cohort --------------------------------

sig_list <- lapply(names(datasets), function(ds) {
  df <- read_tsv(datasets[[ds]], show_col_types = FALSE)
  
  df |>
    filter(confidence == "high_confidence") |>
    select(feature, direction_consensus, maaslin_q, maaslin_coef, confidence) |>
    rename(
      !!paste0(ds, "_dir")  := direction_consensus,
      !!paste0(ds, "_q")    := maaslin_q,
      !!paste0(ds, "_coef") := maaslin_coef,
      !!paste0(ds, "_conf") := confidence
    )
})
names(sig_list) <- names(datasets)

# ---- Build full overlap table --------------------------------------------

all_features <- sig_list |>
  lapply(\(x) x$feature) |>
  unlist() |>
  unique() |>
  sort()

overlap <- tibble(feature = all_features)
for (ds in names(datasets)) {
  overlap <- overlap |> left_join(sig_list[[ds]], by = "feature")
}

dir_cols <- paste0(names(datasets), "_dir")
q_cols   <- paste0(names(datasets), "_q")

overlap <- overlap |>
  mutate(
    n_datasets = rowSums(!is.na(across(all_of(dir_cols)))),
    direction_consistent = apply(
      across(all_of(dir_cols)), 1,
      function(x) { vals <- unique(na.omit(x)); length(vals) == 1 }
    ),
    consensus_direction = apply(
      across(all_of(dir_cols)), 1,
      function(x) {
        vals <- unique(na.omit(x))
        if (length(vals) == 1) vals else NA_character_
      }
    )
  ) |>
  arrange(desc(n_datasets), feature)

write_tsv(overlap, file.path(out_dir, paste0("cross_dataset_overlap_deeparg_", LEVEL, ".tsv")))
cat("Full overlap table saved:", nrow(overlap), "features\n")

# ---- Summary stats -------------------------------------------------------

cat("\nFeatures high-confidence in >=1 dataset:", nrow(overlap), "\n")
cat("Features high-confidence in >=2 datasets:", sum(overlap$n_datasets >= 2), "\n")
cat("Features high-confidence in >=3 datasets:", sum(overlap$n_datasets >= 3), "\n")
cat("Features high-confidence in all 4 datasets:", sum(overlap$n_datasets == 4), "\n")
cat("\n=== Present in all 4 datasets ===\n")
overlap |>
  filter(n_datasets == 4) |>
  select(feature, consensus_direction, all_of(dir_cols), all_of(q_cols)) |>
  print(n = Inf)

# ---- Summary table for thesis -------------------------------------------

summary_table <- overlap |>
  filter(n_datasets >= 2, direction_consistent == TRUE) |>
  select(
    Feature                   = feature,
    `PRJNA813705\n(AML)`      = prjna813705_dir,
    `Kulecka\nLeukemia (AML)` = kulecka_leukemia_dir,
    `Kulecka\nLymphoma (LN)`  = kulecka_lymphoma_dir,
    `CRA\n(NKTCL)`            = cra_dir,
    n_datasets,
    consensus_direction
  ) |>
  mutate(
    across(
      c(`PRJNA813705\n(AML)`, `Kulecka\nLeukemia (AML)`,
        `Kulecka\nLymphoma (LN)`, `CRA\n(NKTCL)`),
      ~ case_when(. == "depleted" ~ "\u2b07",
                  . == "enriched" ~ "\u2b06",
                  is.na(.)        ~ "\u2014")
    )
  ) |>
  arrange(desc(n_datasets), consensus_direction, Feature)

write_tsv(summary_table,
          file.path(out_dir, paste0("summary_table_deeparg_", LEVEL, ".tsv")))
cat("\nSummary table saved:", nrow(summary_table), "features\n")

# ---- ARG class checks ----------------------------------------------------

# Clinically important ARG categories to highlight
broad_spectrum <- c("tetracycline", "multidrug", "beta-lactam", "aminoglycoside",
                    "fluoroquinolone", "glycopeptide")

mls_class <- c("MLS")  # macrolide-lincosamide-streptogramin — dominant in gut

last_resort <- c("colistin", "polymyxin", "carbapenem", "oxazolidinone",
                 "glycopeptide")

# Key tetracycline subtypes (dominant signal from single-sample QC)
tet_subtypes <- c("TETW", "TETM", "TETQ", "TETX", "TETO", "TETP",
                  "TET32", "TET40", "TETA", "TETB")

if (LEVEL == "type") {
  cat("\n=== Broad-spectrum ARG categories (>=2 datasets) ===\n")
  overlap |>
    filter(feature %in% broad_spectrum, n_datasets >= 2) |>
    select(feature, n_datasets, consensus_direction, all_of(dir_cols)) |>
    print(n = Inf)
  
  cat("\n=== Last-resort ARG categories (any dataset) ===\n")
  overlap |>
    filter(feature %in% last_resort) |>
    select(feature, n_datasets, consensus_direction, all_of(dir_cols)) |>
    print(n = Inf)
}

if (LEVEL == "subtype") {
  cat("\n=== Tetracycline subtypes (>=2 datasets) ===\n")
  overlap |>
    filter(feature %in% tet_subtypes, n_datasets >= 2) |>
    select(feature, n_datasets, consensus_direction, all_of(dir_cols)) |>
    print(n = Inf)
}
