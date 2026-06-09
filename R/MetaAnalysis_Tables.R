## =============================================================================
## MetaAnalysis_Tables.R
## CRASH — Cancer-associated gut micRobiome Analysis of Strains in
##         Haematological patients
##
## Purpose:
##   Compiles all per-cohort MetaPhlAn4 analysis outputs into thesis-ready
##   tables. Reads from MetaAnalysis_1.R and MetaAnalysis_2.R outputs.
##
##   Table 1  — Cohort summary: sample counts (case/ctrl), cancer type,
##               accession, country
##   Table 2  — Cross-cohort consensus species (≥3 cohorts, consistent
##               direction); formatted with Unicode arrows
##   Table 3a — Alpha diversity per cohort × metric (median case/ctrl,
##               Wilcoxon U, BH-adjusted p)
##   Table 3b — Beta diversity per cohort × metric (PERMANOVA R²/F/p,
##               betadisper p)
##   Table 4  — Per-cohort DA method triangulation: MaAsLin2 / Wilcoxon /
##               ALDEx2 hit counts and pairwise overlaps
##   Table 5  — Curated functional group results (butyrate producers +
##               opportunistic taxa); enriched/depleted by cohort
##
## Inputs:
##   - Per-cohort cross_method_metaphlan_<DATASET>.tsv
##   - Per-cohort alpha_diversity_<DATASET>.tsv
##   - Per-cohort beta_panel_<metric>_<DATASET>.rds
##   - cross_dataset/cross_dataset_overlap_metaphlan.tsv
##   - Per-cohort metadata Excel files
##
## Outputs (written to thesis_tables/):
##   table1_cohort_summary.tsv / .xlsx
##   table2_consensus_species.tsv / .xlsx
##   table3a_alpha_diversity.tsv / .xlsx
##   table3b_beta_diversity.tsv / .xlsx
##   table4_da_method_counts.tsv / .xlsx
##   table5_functional_groups.tsv / .xlsx
##   thesis_tables_combined.xlsx   (all tables, one sheet each)
##
## Dependencies:
##   tidyverse, readxl, janitor, writexl
##
## Notes:
##   - Run MetaAnalysis_1.R (all 4 cohorts) and MetaAnalysis_2.R first
##   - writexl is used (no Java dependency); auto-installed if missing
## =============================================================================

library(tidyverse)
library(readxl)
library(janitor)

# writexl is lightweight (no Java dep). Install once if missing.
if (!requireNamespace("writexl", quietly = TRUE)) {
  install.packages("writexl")
}
library(writexl)

# ---- Paths & constants ---------------------------------------------------

out_dir <- "D:/MasterThesis/Vik/thesis_tables"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

Q_THRESH <- 0.25

cohorts <- list(
  prjna813705 = list(
    meta_path = "D:/MasterThesis/Vik/PRJNA813705/PRJNA813705_metadata.xlsx",
    base_dir  = "D:/MasterThesis/Vik/PRJNA813705",
    maaslin   = "D:/MasterThesis/Vik/PRJNA813705/maaslin2_metaphlan_prjna813705/all_results.tsv",
    case = "AML",  ctrl = "Control"
  ),
  kulecka_leukemia = list(
    meta_path = "D:/MasterThesis/Vik/Kulecka/Leukemia/Kuleckaleukemia_metadata.xlsx",
    base_dir  = "D:/MasterThesis/Vik/Kulecka/Leukemia",
    maaslin   = "D:/MasterThesis/Vik/Kulecka/Leukemia/maaslin2_metaphlan_kulecka_leukemia/all_results.tsv",
    case = "AML",  ctrl = "Control"
  ),
  kulecka_lymphoma = list(
    meta_path = "D:/MasterThesis/Vik/Kulecka/Lymphoma/Kuleckalymphoma_metadata.xlsx",
    base_dir  = "D:/MasterThesis/Vik/Kulecka/Lymphoma",
    maaslin   = "D:/MasterThesis/Vik/Kulecka/Lymphoma/maaslin2_metaphlan_kulecka_lymphoma/all_results.tsv",
    case = "LN",   ctrl = "Control"
  ),
  cra = list(
    meta_path = "D:/MasterThesis/Vik/CRA007433/CRA007433_metadata.xlsx",
    base_dir  = "D:/MasterThesis/Vik/CRA007433",
    maaslin   = "D:/MasterThesis/Vik/CRA007433/maaslin2_metaphlan_cra/all_results.tsv",
    case = "NKTCL", ctrl = "Control"
  )
)

cohort_name_map <- c(
  prjna813705      = "AML (PRJNA813705)",
  kulecka_leukemia = "AML (Kulecka)",
  kulecka_lymphoma = "LN (Kulecka)",
  cra              = "NKTCL (CRA)"
)

case_map <- c(
  prjna813705 = "AML", kulecka_leukemia = "AML",
  kulecka_lymphoma = "LN", cra = "NKTCL"
)

# Study-level info not present in metadata files. EDIT to match canonical
# citations + platforms; cross-check against the project proposal.
fixed_info <- tribble(
  ~cohort_key,        ~dataset_id,    ~cancer_type,               ~geography, ~platform,          ~reference,
  "prjna813705",      "PRJNA813705",  "Acute Myeloid Leukaemia",  "Belgium",  "Illumina NovaSeq", "Potgens et al. 2022",
  "kulecka_leukemia", "PRJNA1116523", "Acute Myeloid Leukaemia",  "Poland",   "Illumina NovaSeq", "Kulecka et al. 2024",
  "kulecka_lymphoma", "PRJNA1116523", "Lymphoma",                 "Poland",   "Illumina NovaSeq", "Kulecka et al. 2024",
  "cra",              "CRA007433",    "NK/T-cell lymphoma",       "China",    "Illumina (GSA)",   "Shi et al. 2023"
)

# Curated functional group lists (must match the lists used in Script 3)
butyrate_producers <- c(
  "Roseburia_hominis", "Roseburia_intestinalis", "Roseburia_faecis",
  "Roseburia_inulinivorans", "Eubacterium_rectale", "Faecalibacterium_prausnitzii",
  "Coprococcus_catus", "Coprococcus_comes", "Coprococcus_eutactus",
  "Agathobaculum_butyriciproducens", "Lachnospira_eligens", "Ruminococcus_callidus"
)

opportunistic <- c(
  "Clostridium_innocuum", "Enterocloster_bolteae", "Eggerthella_lenta",
  "Hungatella_hathewayi", "Clostridium_symbiosum", "Escherichia_coli",
  "Ruminococcus_gnavus", "Erysipelatoclostridium_ramosum"
)

# ---- Helpers -------------------------------------------------------------

fmt_median_range <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (all(is.na(x))) return(NA_character_)
  sprintf("%.0f (%.0f\u2013%.0f)",
          median(x, na.rm = TRUE),
          min(x, na.rm = TRUE),
          max(x, na.rm = TRUE))
}

find_col <- function(df, candidates) {
  hit <- intersect(candidates, names(df))
  if (length(hit) == 0) NA_character_ else hit[1]
}

# ---- Table 1: Cohort summary ---------------------------------------------

extract_demographics <- function(info) {
  m <- read_excel(info$meta_path) |>
    clean_names() |>
    filter(disease %in% c(info$case, info$ctrl))
  
  n_case <- sum(m$disease == info$case, na.rm = TRUE)
  n_ctrl <- sum(m$disease == info$ctrl, na.rm = TRUE)
  
  age_col <- find_col(m, c("age", "age_years", "age_at_sampling", "age_y"))
  sex_col <- find_col(m, c("sex", "gender"))
  
  # Age split by group if available, so confounding is visible
  if (!is.na(age_col)) {
    age_case <- fmt_median_range(m[[age_col]][m$disease == info$case])
    age_ctrl <- fmt_median_range(m[[age_col]][m$disease == info$ctrl])
  } else {
    age_case <- NA_character_
    age_ctrl <- NA_character_
  }
  
  if (!is.na(sex_col)) {
    sx <- toupper(substr(as.character(m[[sex_col]]), 1, 1))
    n_m <- sum(sx == "M", na.rm = TRUE)
    n_f <- sum(sx == "F", na.rm = TRUE)
    sex_str <- if (n_m + n_f == 0) NA_character_ else sprintf("%d / %d", n_m, n_f)
  } else {
    sex_str <- NA_character_
  }
  
  tibble(
    n_total = nrow(m), n_case = n_case, n_ctrl = n_ctrl,
    age_case = age_case, age_ctrl = age_ctrl, sex_mf = sex_str
  )
}

demographics <- imap_dfr(cohorts, function(info, key) {
  extract_demographics(info) |> mutate(cohort_key = key)
})

table1 <- fixed_info |>
  left_join(demographics, by = "cohort_key") |>
  mutate(
    n_combined = sprintf("%d (%d / %d)", n_total, n_case, n_ctrl)
  ) |>
  transmute(
    `Cohort`               = cohort_name_map[cohort_key],
    `Dataset accession`    = dataset_id,
    `Cancer type`          = cancer_type,
    `Geography`            = geography,
    `N (case / control)`   = n_combined,
    `Age case (median, range)`    = age_case,
    `Age control (median, range)` = age_ctrl,
    `Sex (M / F)`          = sex_mf,
    `Sequencing platform`  = platform,
    `Reference`            = reference
  )

write_tsv(table1, file.path(out_dir, "Table1_cohort_summary.tsv"))

# ---- Load cross-dataset overlap (for Table 2) ---------------------------

overlap <- read_tsv(
  "D:/MasterThesis/Vik/cross_dataset/cross_dataset_overlap_metaphlan.tsv",
  show_col_types = FALSE
) |>
  mutate(
    functional_group = case_when(
      feature %in% butyrate_producers ~ "Butyrate producer",
      feature %in% opportunistic      ~ "Opportunistic",
      TRUE                             ~ "Other"
    )
  )

# ---- Table 2: Cross-cohort consensus species ----------------------------

table2 <- overlap |>
  filter(n_datasets >= 3, direction_consistent == TRUE) |>
  arrange(desc(n_datasets), consensus_direction, feature) |>
  mutate(
    across(ends_with("_coef"), ~ round(.x, 2)),
    across(ends_with("_q"),    ~ signif(.x, 3)),
    feature = str_replace_all(feature, "_", " ")
  ) |>
  transmute(
    Species                = feature,
    `Cohorts sig`          = n_datasets,
    Direction              = consensus_direction,
    `Functional group`     = functional_group,
    `PRJNA813705 coef`     = prjna813705_coef,
    `PRJNA813705 q`        = prjna813705_q,
    `Kulecka AML coef`     = kulecka_leukemia_coef,
    `Kulecka AML q`        = kulecka_leukemia_q,
    `Kulecka LN coef`      = kulecka_lymphoma_coef,
    `Kulecka LN q`         = kulecka_lymphoma_q,
    `CRA NKTCL coef`       = cra_coef,
    `CRA NKTCL q`          = cra_q
  )

write_tsv(table2, file.path(out_dir, "Table2_cross_cohort_consensus.tsv"))

# ---- Table 3a: Alpha diversity -----------------------------------------

alpha_combined <- imap_dfr(cohorts, function(info, key) {
  p <- file.path(info$base_dir, sprintf("tables_metaphlan_%s", key),
                 sprintf("alpha_diversity_%s.tsv", key))
  if (!file.exists(p)) return(NULL)
  read_tsv(p, show_col_types = FALSE) |> mutate(cohort_key = key)
})

table3a <- alpha_combined |>
  mutate(
    Cohort = cohort_name_map[cohort_key],
    ratio  = median_case / median_ctrl,
    across(c(median_ctrl, median_case, ratio), ~ round(.x, 3)),
    p       = signif(p, 3),
    p.adj   = signif(p.adj, 3)
  ) |>
  transmute(
    Cohort, Metric = metric,
    `Median (control)`  = median_ctrl,
    `Median (case)`     = median_case,
    `Case/ctrl ratio`   = ratio,
    `Wilcoxon U`        = U,
    `p`                 = p,
    `q (BH)`            = p.adj
  ) |>
  arrange(Cohort, factor(Metric, levels = c("Observed", "Shannon", "Simpson")))

write_tsv(table3a, file.path(out_dir, "Table3a_alpha_diversity.tsv"))

# ---- Table 3b: Beta diversity ------------------------------------------

read_permanova <- function(path, key, metric_lbl) {
  if (!file.exists(path)) return(NULL)
  df <- read_tsv(path, show_col_types = FALSE)
  row <- df |> filter(term == "disease")
  if (nrow(row) == 0) return(NULL)
  
  # Column name can be `Pr(>F)` or similar after read_tsv
  p_col <- intersect(c("Pr(>F)", "Pr(>=F)", "pval", "p"), names(row))[1]
  
  tibble(
    cohort_key    = key,
    Metric        = metric_lbl,
    Df            = row$Df,
    R2            = round(row$R2, 4),
    `F`           = round(row$`F`, 3),
    `PERMANOVA p` = signif(row[[p_col]], 3)
  )
}

read_betadisper <- function(path, key, metric_lbl) {
  if (!file.exists(path)) return(NULL)
  df <- read_tsv(path, show_col_types = FALSE)
  row <- df |> filter(term == "Groups")
  if (nrow(row) == 0) return(NULL)
  
  tibble(
    cohort_key     = key,
    Metric         = metric_lbl,
    `Betadisper F` = round(row$statistic, 3),
    `Betadisper p` = signif(row$p.value, 3)
  )
}

perm_rows <- map_dfr(names(cohorts), function(key) {
  base <- file.path(cohorts[[key]]$base_dir, sprintf("tables_metaphlan_%s", key))
  bind_rows(
    read_permanova(file.path(base, sprintf("permanova_bray_%s.tsv", key)),    key, "Bray-Curtis"),
    read_permanova(file.path(base, sprintf("permanova_jaccard_%s.tsv", key)), key, "Jaccard")
  )
})

disp_rows <- map_dfr(names(cohorts), function(key) {
  base <- file.path(cohorts[[key]]$base_dir, sprintf("tables_metaphlan_%s", key))
  bind_rows(
    read_betadisper(file.path(base, sprintf("betadisper_anova_bray_%s.tsv", key)),    key, "Bray-Curtis"),
    read_betadisper(file.path(base, sprintf("betadisper_anova_jaccard_%s.tsv", key)), key, "Jaccard")
  )
})

table3b <- perm_rows |>
  left_join(disp_rows, by = c("cohort_key", "Metric")) |>
  mutate(Cohort = cohort_name_map[cohort_key]) |>
  select(Cohort, Metric, Df, R2, `F`, `PERMANOVA p`, `Betadisper F`, `Betadisper p`) |>
  arrange(Cohort, Metric)

write_tsv(table3b, file.path(out_dir, "Table3b_beta_diversity.tsv"))

# ---- Table 4: Per-cohort DA method triangulation -----------------------

table4 <- imap_dfr(cohorts, function(info, key) {
  p <- file.path(info$base_dir, sprintf("tables_metaphlan_%s", key),
                 sprintf("cross_method_metaphlan_%s.tsv", key))
  if (!file.exists(p)) return(NULL)
  
  cm <- read_tsv(p, show_col_types = FALSE)
  
  m_sig <- replace_na(cm$maaslin_sig, FALSE)
  w_sig <- replace_na(cm$wilcox_sig,  FALSE)
  a_sig <- replace_na(cm$aldex_sig,   FALSE)
  
  tibble(
    cohort_key            = key,
    `N tested`            = nrow(cm),
    `MaAsLin2 sig`        = sum(m_sig),
    `Wilcoxon sig`        = sum(w_sig),
    `ALDEx2 sig`          = sum(a_sig),
    `MaAsLin2 & Wilcoxon` = sum(m_sig & w_sig),
    `MaAsLin2 & ALDEx2`   = sum(m_sig & a_sig),
    `All three`           = sum(m_sig & w_sig & a_sig),
    `High-confidence`     = sum(cm$confidence == "high_confidence", na.rm = TRUE)
  )
}) |>
  mutate(Cohort = cohort_name_map[cohort_key]) |>
  select(Cohort, everything(), -cohort_key)

write_tsv(table4, file.path(out_dir, "Table4_DA_triangulation.tsv"))

# ---- Table 5: Curated functional groups (from full MaAsLin2 results) ---
# Pulls coef + q from every cohort's all_results.tsv so taxa that were
# subthreshold in a cohort still appear with their estimated effect.
# This matches Figures 3 and 4 in the lab-meeting output.

maaslin_all <- imap_dfr(cohorts, function(info, key) {
  if (!file.exists(info$maaslin)) return(NULL)
  read_tsv(info$maaslin, show_col_types = FALSE) |>
    filter(metadata == "disease") |>
    mutate(
      cohort_key = key,
      coef_case  = if_else(value == case_map[[key]], coef, -coef)
    ) |>
    select(feature, cohort_key, coef_case, qval)
})

curated_long <- maaslin_all |>
  filter(feature %in% c(butyrate_producers, opportunistic)) |>
  mutate(
    functional_group = if_else(feature %in% butyrate_producers,
                               "Butyrate producer", "Opportunistic"),
    coef_case = round(coef_case, 2),
    qval      = signif(qval, 3)
  )

# Pivot to wide: coef_<cohort> and q_<cohort> columns
table5 <- curated_long |>
  pivot_wider(
    id_cols     = c(feature, functional_group),
    names_from  = cohort_key,
    values_from = c(coef_case, qval),
    names_glue  = "{cohort_key}_{.value}"
  ) |>
  mutate(feature = str_replace_all(feature, "_", " ")) |>
  # Build a sort helper: n cohorts with q<0.25
  rowwise() |>
  mutate(
    n_sig = sum(c_across(ends_with("_qval")) < Q_THRESH, na.rm = TRUE)
  ) |>
  ungroup() |>
  arrange(functional_group, desc(n_sig), feature) |>
  transmute(
    `Functional group`  = functional_group,
    Species             = feature,
    `Cohorts sig (q<0.25)` = n_sig,
    `PRJNA813705 coef`  = prjna813705_coef_case,
    `PRJNA813705 q`     = prjna813705_qval,
    `Kulecka AML coef`  = kulecka_leukemia_coef_case,
    `Kulecka AML q`     = kulecka_leukemia_qval,
    `Kulecka LN coef`   = kulecka_lymphoma_coef_case,
    `Kulecka LN q`      = kulecka_lymphoma_qval,
    `CRA NKTCL coef`    = cra_coef_case,
    `CRA NKTCL q`       = cra_qval
  )

write_tsv(table5, file.path(out_dir, "Table5_functional_groups.tsv"))

# ---- Combined Excel workbook -------------------------------------------

write_xlsx(
  list(
    `Table 1 - Cohorts`           = table1,
    `Table 2 - Consensus species` = table2,
    `Table 3a - Alpha diversity`  = table3a,
    `Table 3b - Beta diversity`   = table3b,
    `Table 4 - DA triangulation`  = table4,
    `Table 5 - Functional groups` = table5
  ),
  path = file.path(out_dir, "Thesis_Tables_MetaPhlAn.xlsx")
)

cat("\n=== Thesis tables complete ===\n")
cat("Directory: ", out_dir, "\n", sep = "")
cat("  - Table1_cohort_summary.tsv       (", nrow(table1), " rows)\n", sep = "")
cat("  - Table2_cross_cohort_consensus.tsv (", nrow(table2), " rows)\n", sep = "")
cat("  - Table3a_alpha_diversity.tsv     (", nrow(table3a), " rows)\n", sep = "")
cat("  - Table3b_beta_diversity.tsv      (", nrow(table3b), " rows)\n", sep = "")
cat("  - Table4_DA_triangulation.tsv     (", nrow(table4), " rows)\n", sep = "")
cat("  - Table5_functional_groups.tsv    (", nrow(table5), " rows)\n", sep = "")
cat("  - Thesis_Tables_MetaPhlAn.xlsx    (all sheets)\n")
