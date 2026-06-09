## =============================================================================
## MetaAnalysis_2.R
## CRASH — Cancer-associated gut micRobiome Analysis of Strains in
##         Haematological patients
##
## Purpose:
##   Cross-cohort overlap analysis of MetaPhlAn4 differential abundance results.
##   Reads per-cohort high-confidence hits from MetaAnalysis_1.R and identifies
##   taxa that replicate directionally across multiple cohorts.
##
##   High-confidence is defined as: MaAsLin2 significant (q < 0.25) AND at
##   least one supporting method (Wilcoxon or ALDEx2) agreeing on direction.
##   This guards against single-method false positives before cross-cohort
##   comparison.
##
## Input:
##   cross_method_metaphlan_<DATASET>.tsv — one per cohort, produced by
##   MetaAnalysis_1.R (tables_metaphlan_<DATASET>/ subdirectories)
##
## Outputs (written to cross_dataset/):
##   cross_dataset_overlap_metaphlan.tsv  Full overlap table with n_datasets,
##                                        direction consistency, and per-cohort
##                                        q-values and coefficients
##   summary_table_metaphlan.tsv          Arrow-formatted thesis table (≥2
##                                        cohorts, direction consistent)
##
## Also prints to console:
##   - Feature counts at ≥1/2/3/4 cohort thresholds
##   - Pan-cohort hits (n=4)
##   - Butyrate producer depletion and opportunistic enrichment summaries
##
## Dependencies:
##   tidyverse
##
## Notes:
##   - Run MetaAnalysis_1.R on all four cohorts before this script
##   - Direction consensus requires all non-NA cohorts to agree
##   - Summary table uses Unicode arrows (↑/↓/—) for formatting
## =============================================================================

library(tidyverse)

# ---- Configuration -------------------------------------------------------

Q_THRESH <- 0.25

datasets <- list(
  prjna813705      = "D:/MasterThesis/Vik/PRJNA813705/tables_metaphlan_prjna813705/cross_method_metaphlan_prjna813705.tsv",
  kulecka_leukemia = "D:/MasterThesis/Vik/Kulecka/Leukemia/tables_metaphlan_kulecka_leukemia/cross_method_metaphlan_kulecka_leukemia.tsv",
  kulecka_lymphoma = "D:/MasterThesis/Vik/Kulecka/Lymphoma/tables_metaphlan_kulecka_lymphoma/cross_method_metaphlan_kulecka_lymphoma.tsv",
  cra              = "D:/MasterThesis/Vik/CRA007433/tables_metaphlan_cra/cross_method_metaphlan_cra.tsv"
)

cancer_labels <- c(
  prjna813705      = "AML (PRJNA813705)",
  kulecka_leukemia = "AML (Kulecka)",
  kulecka_lymphoma = "LN (Kulecka)",
  cra              = "NKTCL (CRA)"
)

out_dir <- "D:/MasterThesis/Vik/cross_dataset"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Load and extract high-confidence hits per cohort --------------------
#
# Change vs previous version: filter(maaslin_sig == TRUE) -> filter(confidence == "high_confidence")
# Rationale: per-cohort hits must already be method-robust (MaAsLin2 sig + at
# least one other method agreeing on direction) before we ask whether they
# replicate across cohorts. Reduces single-method false positives.

sig_list <- lapply(names(datasets), function(ds) {
  df <- read_tsv(datasets[[ds]], show_col_types = FALSE)
  df |>
    filter(confidence == "high_confidence") |>
    select(feature, maaslin_direction, maaslin_q, maaslin_coef, confidence) |>
    rename(
      !!paste0(ds, "_dir")  := maaslin_direction,
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

# Count datasets and check direction consistency
dir_cols <- paste0(names(datasets), "_dir")

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

write_tsv(overlap, file.path(out_dir, "cross_dataset_overlap_metaphlan.tsv"))
cat("Full overlap table saved:", nrow(overlap), "features\n")

# ---- Summary stats -------------------------------------------------------

cat("\nFeatures high-confidence in >=1 dataset:", nrow(overlap), "\n")
cat("Features high-confidence in >=2 datasets:", sum(overlap$n_datasets >= 2), "\n")
cat("Features high-confidence in >=3 datasets:", sum(overlap$n_datasets >= 3), "\n")
cat("Features high-confidence in all 4 datasets:", sum(overlap$n_datasets == 4), "\n")

cat("\n=== Present in all 4 datasets ===\n")
overlap |>
  filter(n_datasets == 4) |>
  select(feature, all_of(dir_cols), direction_consistent) |>
  print(n = Inf)

# ---- Formatted summary table for thesis ----------------------------------

summary_table <- overlap |>
  filter(n_datasets >= 2, direction_consistent == TRUE) |>
  select(
    Feature = feature,
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

write_tsv(summary_table, file.path(out_dir, "summary_table_metaphlan.tsv"))
cat("\nSummary table saved:", nrow(summary_table), "features\n")

# ---- Functional group checks ---------------------------------------------

butyrate_producers <- c(
  "Roseburia_hominis", "Roseburia_intestinalis", "Roseburia_faecis", "Roseburia_inulinivorans",
  "Eubacterium_rectale", "Faecalibacterium_prausnitzii",
  "Coprococcus_catus", "Coprococcus_comes", "Coprococcus_eutactus",
  "Agathobaculum_butyriciproducens", "Ruminococcus_callidus", "Lachnospira_eligens",
  "Lachnospira_pectinoschiza", "Faecalicatena_fissicatena", "Anaerostipes_hadrus"
)

opportunistic <- c(
  "Clostridium_innocuum", "Enterocloster_bolteae", "Eggerthella_lenta",
  "Hungatella_hathewayi", "Clostridium_symbiosum", "Escherichia_coli",
  "Ruminococcus_gnavus", "Erysipelatoclostridium_ramosum"
)

cat("\n=== Butyrate producers depleted (>=2 datasets) ===\n")
overlap |>
  filter(feature %in% butyrate_producers, n_datasets >= 2) |>
  select(feature, n_datasets, consensus_direction, all_of(dir_cols)) |>
  arrange(desc(n_datasets), feature) |>
  print(n = Inf)

cat("\n=== Opportunistic taxa enriched (>=2 datasets) ===\n")
overlap |>
  filter(feature %in% opportunistic, n_datasets >= 2) |>
  select(feature, n_datasets, consensus_direction, all_of(dir_cols)) |>
  arrange(desc(n_datasets), feature) |>
  print(n = Inf)
