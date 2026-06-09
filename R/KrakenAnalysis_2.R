## =============================================================================
## KrakenAnalysis_2.R
## CRASH — Cancer-associated gut micRobiome Analysis of Strains in
##         Haematological patients
##
## Purpose:
##   Cross-cohort overlap analysis of Kraken2/Bracken differential abundance
##   results. Mirrors MetaAnalysis_2.R exactly but reads from Kraken outputs.
##   Identifies taxa that replicate directionally across multiple cohorts as
##   independent validation of MetaPhlAn4 findings.
##
## Input:
##   cross_method_kraken_<DATASET>.tsv — one per cohort, produced by
##   KrakenAnalysis_1.R (tables_kraken_<DATASET>/ subdirectories)
##
## Outputs (written to cross_dataset_kraken/):
##   cross_dataset_overlap_kraken.tsv   Full overlap table
##   summary_table_kraken.tsv           Arrow-formatted thesis table
##
## Dependencies:
##   tidyverse
##
## Notes:
##   - Run KrakenAnalysis_1.R on all four cohorts before this script
##   - High-confidence filter: MaAsLin2 sig + ≥1 method directional consensus
## =============================================================================

library(tidyverse)

# ---- Configuration -------------------------------------------------------

Q_THRESH <- 0.25

datasets <- list(
  prjna813705      = "D:/MasterThesis/Vik/PRJNA813705/tables_kraken_prjna813705/cross_method_kraken_prjna813705.tsv",
  kulecka_leukemia = "D:/MasterThesis/Vik/Kulecka/Leukemia/tables_kraken_kulecka_leukemia/cross_method_kraken_kulecka_leukemia.tsv",
  kulecka_lymphoma = "D:/MasterThesis/Vik/Kulecka/Lymphoma/tables_kraken_kulecka_lymphoma/cross_method_kraken_kulecka_lymphoma.tsv",
  cra              = "D:/MasterThesis/Vik/CRA007433/tables_kraken_cra/cross_method_kraken_cra.tsv"
)

cancer_labels <- c(
  prjna813705      = "AML (PRJNA813705)",
  kulecka_leukemia = "AML (Kulecka)",
  kulecka_lymphoma = "LN (Kulecka)",
  cra              = "NKTCL (CRA)"
)

out_dir <- "D:/MasterThesis/Vik/cross_dataset_kraken"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Normalizer must match the one in the per-cohort pipeline, because Kraken
# species names vary slightly across databases (brackets, renamed genera).
# Cross-cohort matching uses the normalized key, then a canonical display
# name is chosen from whichever cohort first reported the feature.
normalize_name <- function(x) {
  x <- tolower(x)
  x <- gsub("^x\\.+", "", x)
  x <- gsub("[^a-z0-9]", "", x)
  x
}

# ---- Load per-cohort high-confidence hits -------------------------------

sig_list <- lapply(names(datasets), function(ds) {
  df <- read_tsv(datasets[[ds]], show_col_types = FALSE) |>
    filter(confidence == "high_confidence") |>
    mutate(norm_key = normalize_name(feature)) |>
    select(feature, norm_key, maaslin_direction, maaslin_q, maaslin_coef)
  
  df |>
    rename(
      !!paste0(ds, "_name") := feature,
      !!paste0(ds, "_dir")  := maaslin_direction,
      !!paste0(ds, "_q")    := maaslin_q,
      !!paste0(ds, "_coef") := maaslin_coef
    )
})
names(sig_list) <- names(datasets)

# ---- Build overlap keyed on normalized name ------------------------------

all_keys <- sig_list |>
  lapply(\(x) x$norm_key) |>
  unlist() |>
  unique() |>
  sort()

overlap <- tibble(norm_key = all_keys)
for (ds in names(datasets)) {
  overlap <- overlap |> left_join(sig_list[[ds]], by = "norm_key")
}

# Choose a canonical display name per feature: first non-NA cohort name.
name_cols <- paste0(names(datasets), "_name")
overlap <- overlap |>
  rowwise() |>
  mutate(
    feature = {
      vals <- c_across(all_of(name_cols))
      vals <- vals[!is.na(vals)]
      if (length(vals) == 0) NA_character_ else vals[1]
    }
  ) |>
  ungroup() |>
  select(feature, norm_key, everything(), -all_of(name_cols))

# ---- Direction consistency and cohort counts -----------------------------

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

write_tsv(overlap, file.path(out_dir, "cross_dataset_overlap_kraken.tsv"))
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

# ---- Formatted summary for thesis ---------------------------------------

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

write_tsv(summary_table, file.path(out_dir, "summary_table_kraken.tsv"))
cat("\nSummary table saved:", nrow(summary_table), "features\n")

# ---- Functional group checks (Kraken naming) ----------------------------

# Kraken species names use spaces, and some taxa appear with "[Genus]"
# brackets (disputed genus assignments) or under renamed genera.
# Any list entries not present in the overlap are ignored silently.
butyrate_producers <- c(
  "Roseburia hominis", "Roseburia intestinalis", "Roseburia faecis",
  "Roseburia inulinivorans",
  "Eubacterium rectale", "[Eubacterium] rectale",
  "Faecalibacterium prausnitzii",
  "Coprococcus catus", "Coprococcus comes", "Coprococcus eutactus",
  "Agathobaculum butyriciproducens",
  "Lachnospira eligens", "[Eubacterium] eligens",
  "Ruminococcus callidus", "Anaerostipes hadrus"
)

opportunistic <- c(
  "[Clostridium] innocuum", "Clostridium innocuum",
  "Enterocloster bolteae", "[Clostridium] bolteae",
  "Eggerthella lenta",
  "Hungatella hathewayi", "[Clostridium] hathewayi",
  "[Clostridium] symbiosum", "Clostridium symbiosum",
  "Escherichia coli",
  "Mediterraneibacter gnavus", "Ruminococcus gnavus",
  "Erysipelatoclostridium ramosum", "[Clostridium] ramosum"
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
