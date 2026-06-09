## =============================================================================
## DeepARG_1.R
## CRASH — Cancer-associated gut micRobiome Analysis of Strains in
##         Haematological patients
##
## Purpose:
##   Per-cohort differential abundance analysis of antimicrobial resistance
##   genes (ARGs) from deepARG v1.0.4 short_reads_pipeline output.
##   Runs at two resolution levels (ARG category "type" and "subtype") and
##   incorporates a Fisher exact prevalence test alongside the standard
##   continuous DA framework.
##
##   Methods applied:
##     - MaAsLin2 (primary; LOG transform on normalised ARG abundances)
##     - Wilcoxon rank-sum (rank-biserial effect size)
##     - ALDEx2 (sensitivity; pseudo-count scaled)
##     - Fisher exact test (presence/absence; prevalence shift)
##     - High-confidence: MaAsLin2 sig + ≥1 method + directional consensus
##
## Usage:
##   Set DATASET and LEVEL (line 46), then source the script.
##     DATASET: "kulecka_leukemia" | "kulecka_lymphoma" | "cra" | "prjna813705"
##     LEVEL:   "type" | "subtype"
##
## Outputs (per cohort/level):
##   plots/    PNG figures
##   tables/   TSV result tables (cross_method, wilcoxon, aldex2, fisher)
##   maaslin/  MaAsLin2 full output folder
##
## Dependencies:
##   readxl, janitor, tidyverse, phyloseq, Maaslin2, ALDEx2,
##   ggpubr, vegan, viridis, pairwiseAdonis, broom
##
## Notes:
##   - deepARG normalised by 16S rRNA (`--bowtie_16s_identity 100`);
##     stringent threshold — state in Methods
##   - Location covariate included in MaAsLin2 and alpha LM for Kulecka cohorts
##   - ALDEx2 direction derived from log2FC_mean, not effect sign
## =============================================================================
############################################################################
# CONFIGURATION
############################################################################

DATASET <- "kulecka_lymphoma" # Options: "kulecka_leukemia", "kulecka_lymphoma", "cra", "prjna813705"
LEVEL   <- "type" # Options: "type", "subtype"

cfg <- list(
  kulecka_leukemia = list(
    wdir         = "D:/MasterThesis/Vik/Kulecka/Leukemia",
    meta_file    = "D:/MasterThesis/Vik/Kulecka/Leukemia/Kuleckaleukemia_metadata.xlsx",
    deeparg_dir  = "D:/MasterThesis/Vik/DeepARG_ALL",
    case_label   = "AML",   ctrl_label = "Control",
    comparison   = "AML",   location_col = "location"
  ),
  kulecka_lymphoma = list(
    wdir         = "D:/MasterThesis/Vik/Kulecka/Lymphoma",
    meta_file    = "D:/MasterThesis/Vik/Kulecka/Lymphoma/Kuleckalymphoma_metadata.xlsx",
    deeparg_dir  = "D:/MasterThesis/Vik/DeepARG_ALL",
    case_label   = "LN",    ctrl_label = "Control",
    comparison   = "LN",    location_col = "location"
  ),
  cra = list(
    wdir         = "D:/MasterThesis/Vik/CRA007433",
    meta_file    = "D:/MasterThesis/Vik/CRA007433/CRA007433_metadata.xlsx",
    deeparg_dir  = "D:/MasterThesis/Vik/DeepARG_ALL",
    case_label   = "NKTCL", ctrl_label = "Control",
    comparison   = "NKTCL", location_col = "location"
  ),
  prjna813705 = list(
    wdir         = "D:/MasterThesis/Vik/PRJNA813705",
    meta_file    = "D:/MasterThesis/Vik/PRJNA813705/PRJNA813705_metadata.xlsx",
    deeparg_dir  = "D:/MasterThesis/Vik/DeepARG_ALL",
    case_label   = "AML",   ctrl_label = "Control",
    comparison   = "AML",   location_col = "location"
  )
)

stopifnot(DATASET %in% names(cfg))
stopifnot(LEVEL %in% c("type", "subtype"))

C          <- cfg[[DATASET]]
wdir       <- C$wdir
CASE       <- C$case_label
CTRL       <- C$ctrl_label
COMPARISON <- C$comparison
LOC_COL    <- C$location_col

Q_THRESH       <- 0.25
MIN_PREVALENCE <- 0.10
ALDEX_SCALE    <- 1e6
SEED           <- 42

############################################################################
# Packages
############################################################################

library(readxl)
library(janitor)
library(tidyverse)
library(Maaslin2)
library(ALDEx2)
library(vegan)
library(viridis)
library(ggpubr)
library(forcats)
library(broom)
library(pheatmap)

select <- dplyr::select

############################################################################
# Output directories
############################################################################

out_suffix  <- paste0(DATASET, "_", LEVEL)
plot_dir    <- file.path("D:/MasterThesis/Vik/DeepARG", out_suffix, "plots")
tsv_dir     <- file.path("D:/MasterThesis/Vik/DeepARG", out_suffix, "tables")
maaslin_dir <- file.path("D:/MasterThesis/Vik/DeepARG", out_suffix, "maaslin2")

dir.create(plot_dir,    showWarnings = FALSE, recursive = TRUE)
dir.create(tsv_dir,     showWarnings = FALSE, recursive = TRUE)
dir.create(maaslin_dir, showWarnings = FALSE, recursive = TRUE)

############################################################################
# Helpers
############################################################################

tidy_ids <- function(x) gsub("\\s+", "_", trimws(as.character(x)))

safe_median <- function(x) {
  if (length(x) == 0 || all(is.na(x))) return(NA_real_)
  median(x, na.rm = TRUE)
}

safe_wilcox <- function(x, y) {
  if (length(x) < 2 || length(y) < 2) return(NULL)
  if (all(is.na(x)) || all(is.na(y))) return(NULL)
  if (sd(c(x, y), na.rm = TRUE) == 0) return(NULL)
  wilcox.test(x, y, exact = FALSE)
}

calc_prevalence <- function(mat) rowMeans(mat > 0, na.rm = TRUE)

############################################################################
# Load metadata
############################################################################

meta <- read_excel(C$meta_file) |>
  clean_names() |>
  filter(disease %in% c(CASE, CTRL)) |>
  mutate(
    disease   = factor(disease, levels = c(CTRL, CASE)),
    sample_id = tidy_ids(sample)
  ) |>
  distinct(sample_id, .keep_all = TRUE)

cat("Metadata loaded:", nrow(meta), "samples\n")
print(table(meta$disease, useNA = "ifany"))

############################################################################
# Load deepARG quant files
############################################################################

ext_map      <- list(type    = ".clean.deeparg.mapping.ARG.merged.quant.type",
                     subtype = ".clean.deeparg.mapping.ARG.merged.quant.subtype")
feat_col_map <- list(type    = "#ARG-category",
                     subtype = "#ARG-group")

ext      <- ext_map[[LEVEL]]
feat_col <- feat_col_map[[LEVEL]]

cat("Loading deepARG quant files (", LEVEL, ")...\n", sep = "")

files <- list.files(C$deeparg_dir,
                    pattern    = paste0(gsub("\\.", "\\\\.", ext), "$"),
                    recursive  = TRUE,
                    full.names = TRUE)
cat("Files found:", length(files), "\n")

df_raw <- lapply(files, function(f) {
  sample_id <- basename(dirname(f))
  read_tsv(f, show_col_types = FALSE) |>
    mutate(sample = tidy_ids(sample_id))
}) |> bind_rows()

valid_samples <- intersect(unique(df_raw$sample), meta$sample_id)
df_raw        <- df_raw |> filter(sample %in% valid_samples)
cat("Samples with both metadata and deepARG output:", length(valid_samples), "\n")

# ---- Safe matrix construction (sum duplicates) --------------------------
mat <- df_raw |>
  select(sample, feature = !!sym(feat_col), `16s-NormalizedReadCount`) |>
  group_by(sample, feature) |>
  summarise(value = sum(`16s-NormalizedReadCount`, na.rm = TRUE), .groups = "drop") |>
  pivot_wider(names_from  = sample,
              values_from = value,
              values_fill = 0) |>
  column_to_rownames("feature") |>
  as.matrix()

cat("Raw matrix:", nrow(mat), "features x", ncol(mat), "samples\n")

prev     <- calc_prevalence(mat)
mat_filt <- mat[prev >= MIN_PREVALENCE, , drop = FALSE]
cat("After prevalence filter (>=", MIN_PREVALENCE, "):", nrow(mat_filt), "features\n")

# ---- Align metadata -----------------------------------------------------
meta_df <- meta |>
  filter(sample_id %in% colnames(mat_filt)) |>
  as.data.frame()
rownames(meta_df) <- meta_df$sample_id

shared   <- intersect(colnames(mat_filt), rownames(meta_df))
mat_filt <- mat_filt[, shared, drop = FALSE]
meta_df  <- meta_df[shared, , drop = FALSE]

# ---- Zero-sample check --------------------------------------------------
keep_samples <- colSums(mat_filt, na.rm = TRUE) > 0
if (any(!keep_samples)) {
  cat("Dropping samples with zero ARG signal after filtering:",
      sum(!keep_samples), "\n")
}
mat_filt <- mat_filt[, keep_samples, drop = FALSE]
meta_df  <- meta_df[colnames(mat_filt), , drop = FALSE]

grp <- as.character(meta_df$disease)
names(grp) <- rownames(meta_df)

cat("Final:", nrow(mat_filt), "features x", ncol(mat_filt), "samples\n")
print(table(grp))

############################################################################
# 1. MaAsLin2 (location-adjusted where applicable)
############################################################################

cat("\n=== MaAsLin2 ===\n")

input_maaslin <- as.data.frame(t(mat_filt))

meta_maaslin <- data.frame(
  disease   = factor(grp[rownames(input_maaslin)], levels = c(CTRL, CASE)),
  row.names = rownames(input_maaslin)
)

fixed_effects_use <- "disease"
use_location <- !is.null(LOC_COL) &&
  LOC_COL %in% colnames(meta_df) &&
  length(unique(na.omit(meta_df[[LOC_COL]]))) > 1

if (use_location) {
  meta_maaslin[[LOC_COL]] <- factor(meta_df[rownames(input_maaslin), LOC_COL])
  fixed_effects_use        <- c("disease", LOC_COL)
  cat("MaAsLin2: adjusting for", LOC_COL, "\n")
}

fit_data <- Maaslin2(
  input_data       = input_maaslin,
  input_metadata   = meta_maaslin,
  output           = maaslin_dir,
  fixed_effects    = fixed_effects_use,
  normalization    = "NONE",
  transform        = "LOG",
  min_prevalence   = 0.0,
  analysis_method  = "LM",
  correction       = "BH",
  standardize      = FALSE,
  max_significance = Q_THRESH
)

# Keep ALL disease rows before filtering so non-sig coefs are not lost
results <- read_tsv(file.path(maaslin_dir, "all_results.tsv"), show_col_types = FALSE)

res_case_all <- results |>
  filter(metadata == "disease") |>
  mutate(
    case_oriented_coef = case_when(
      value == CASE ~ coef,
      value == CTRL ~ -coef,
      TRUE          ~ NA_real_
    ),
    direction_case = case_when(
      case_oriented_coef > 0 ~ "enriched",
      case_oriented_coef < 0 ~ "depleted",
      TRUE                   ~ NA_character_
    )
  )

res_case    <- res_case_all |> filter(qval < Q_THRESH)
maaslin_sig <- res_case$feature
cat("MaAsLin2 significant (q <", Q_THRESH, "):", nrow(res_case), "\n")

############################################################################
# 2. Wilcoxon rank-sum (direction from rank_biserial; mean-based metrics)
############################################################################

cat("\n=== Wilcoxon rank-sum ===\n")

wilcox_results <- lapply(rownames(mat_filt), function(feat) {
  vals <- mat_filt[feat, ]
  ctrl <- vals[grp == CTRL]
  case <- vals[grp == CASE]
  
  wt <- safe_wilcox(case, ctrl)
  if (is.null(wt)) return(NULL)
  
  n_case <- length(case)
  n_ctrl <- length(ctrl)
  U      <- as.numeric(wt$statistic)
  
  data.frame(
    feature       = feat,
    prevalence    = prev[feat],
    prev_ctrl     = mean(ctrl > 0, na.rm = TRUE),
    prev_case     = mean(case > 0, na.rm = TRUE),
    U             = U,
    p             = wt$p.value,
    rank_biserial = (2 * U) / (n_case * n_ctrl) - 1,
    median_ctrl   = safe_median(ctrl),
    median_case   = safe_median(case),
    mean_ctrl     = mean(ctrl, na.rm = TRUE),
    mean_case     = mean(case, na.rm = TRUE),
    log2FC        = log2((safe_median(case) + 1e-10) /
                           (safe_median(ctrl) + 1e-10)),
    log2FC_mean   = log2((mean(case, na.rm = TRUE) + 1e-10) /
                           (mean(ctrl, na.rm = TRUE) + 1e-10)),
    stringsAsFactors = FALSE
  )
}) |>
  bind_rows() |>
  mutate(
    p.adj = p.adjust(p, method = "BH"),
    wilcox_direction = case_when(
      rank_biserial > 0 ~ "enriched",
      rank_biserial < 0 ~ "depleted",
      TRUE              ~ NA_character_
    )
  ) |>
  arrange(p)

wilcox_sig <- wilcox_results |> filter(p.adj < Q_THRESH) |> pull(feature)
cat("Wilcoxon significant (q <", Q_THRESH, "):", length(wilcox_sig), "\n")
cat("MaAsLin2-Wilcoxon overlap:", length(intersect(maaslin_sig, wilcox_sig)), "\n")

write_tsv(wilcox_results,
          file.path(tsv_dir, paste0("wilcoxon_deeparg_", LEVEL, "_", DATASET, ".tsv")))

############################################################################
# 3. ALDEx2 (explicit condition ordering; direction from log2FC_mean)
############################################################################

cat("\n=== ALDEx2 ===\n")
cat("Note: scaled deepARG 16S-normalised counts as pseudo-counts. Sensitivity only.\n")

pseudo_counts <- round(mat_filt * ALDEX_SCALE)
mode(pseudo_counts) <- "integer"
pseudo_counts <- pseudo_counts[rowSums(pseudo_counts) > 0, , drop = FALSE]

aldex_conditions <- as.character(factor(grp[colnames(pseudo_counts)], levels = c(CTRL, CASE)))

set.seed(SEED)
aldex_res <- aldex(
  reads      = pseudo_counts,
  conditions = aldex_conditions,
  mc.samples = 128,
  test       = "t",
  effect     = TRUE,
  denom      = "all",
  verbose    = FALSE
) |>
  rownames_to_column("feature") |>
  left_join(
    wilcox_results |> select(feature, log2FC_mean),
    by = "feature"
  ) |>
  mutate(
    aldex_direction = case_when(
      log2FC_mean > 0 ~ "enriched",
      log2FC_mean < 0 ~ "depleted",
      TRUE            ~ NA_character_
    )
  )

aldex_sig <- aldex_res |> filter(wi.eBH < Q_THRESH) |> pull(feature)
cat("ALDEx2 significant (wi.eBH <", Q_THRESH, "):", length(aldex_sig), "\n")

write_tsv(aldex_res,
          file.path(tsv_dir, paste0("aldex2_deeparg_", LEVEL, "_", DATASET, ".tsv")))

############################################################################
# 4. Fisher prevalence test
############################################################################

cat("\n=== Fisher prevalence test ===\n")

fisher_results <- lapply(rownames(mat_filt), function(feat) {
  vals <- mat_filt[feat, ]
  tab  <- table(
    carrier = vals > 0,
    disease = factor(grp, levels = c(CTRL, CASE))
  )
  if (!all(dim(tab) == c(2, 2))) return(NULL)
  ft <- fisher.test(tab)
  data.frame(
    feature    = feat,
    fisher_p   = ft$p.value,
    odds_ratio = unname(ft$estimate),
    prev_ctrl  = mean(vals[grp == CTRL] > 0),
    prev_case  = mean(vals[grp == CASE] > 0),
    stringsAsFactors = FALSE
  )
}) |>
  bind_rows() |>
  mutate(
    fisher_q         = p.adjust(fisher_p, method = "BH"),
    fisher_direction = case_when(
      odds_ratio > 1 ~ "enriched",
      odds_ratio < 1 ~ "depleted",
      TRUE           ~ NA_character_
    )
  ) |>
  arrange(fisher_p)

fisher_sig <- fisher_results |> filter(fisher_q < Q_THRESH) |> pull(feature)
cat("Fisher significant (q <", Q_THRESH, "):", length(fisher_sig), "\n")

write_tsv(fisher_results,
          file.path(tsv_dir, paste0("fisher_deeparg_", LEVEL, "_", DATASET, ".tsv")))

############################################################################
# 5. Cross-method consensus (direction from significant methods only)
############################################################################

full_comparison <- wilcox_results |>
  select(feature, prevalence, prev_ctrl, prev_case,
         wilcox_q = p.adj, rank_biserial, log2FC, log2FC_mean,
         wilcox_direction) |>
  full_join(
    res_case_all |>
      select(feature,
             maaslin_coef      = case_oriented_coef,
             maaslin_q         = qval,
             maaslin_direction = direction_case),
    by = "feature"
  ) |>
  full_join(
    aldex_res |> select(feature, aldex_q = wi.eBH,
                        aldex_effect = effect, aldex_direction),
    by = "feature"
  ) |>
  full_join(
    fisher_results |> select(feature, fisher_q, odds_ratio, fisher_direction),
    by = "feature"
  ) |>
  mutate(
    wilcox_sig  = !is.na(wilcox_q)  & wilcox_q  < Q_THRESH,
    maaslin_sig = !is.na(maaslin_q) & maaslin_q < Q_THRESH,
    aldex_sig   = !is.na(aldex_q)   & aldex_q   < Q_THRESH,
    n_sig       = rowSums(cbind(wilcox_sig, maaslin_sig, aldex_sig), na.rm = TRUE)
  ) |>
  rowwise() |>
  mutate(
    sig_dirs = list(na.omit(c(
      if (wilcox_sig)  wilcox_direction  else NA_character_,
      if (maaslin_sig) maaslin_direction else NA_character_,
      if (aldex_sig)   aldex_direction   else NA_character_
    ))),
    direction_consensus = if (length(sig_dirs) >= 2 && length(unique(sig_dirs)) == 1) {
      sig_dirs[[1]]
    } else {
      NA_character_
    },
    confidence = case_when(
      maaslin_sig & n_sig >= 2 & !is.na(direction_consensus) ~ "high_confidence",
      n_sig >= 2  & !is.na(direction_consensus)              ~ "moderate_confidence",
      n_sig == 1                                             ~ "single_method",
      TRUE                                                   ~ "not_significant"
    ),
    fisher_supports_consensus = case_when(
      is.na(direction_consensus)                                         ~ NA,
      !is.na(fisher_direction) & fisher_direction == direction_consensus ~ TRUE,
      !is.na(fisher_direction) & fisher_direction != direction_consensus ~ FALSE,
      TRUE                                                               ~ NA
    )
  ) |>
  ungroup() |>
  select(-sig_dirs) |>
  arrange(maaslin_q, wilcox_q, aldex_q)

write_tsv(full_comparison,
          file.path(tsv_dir, paste0("cross_method_deeparg_", LEVEL, "_", DATASET, ".tsv")))

cat("\nConfidence summary:\n")
print(table(full_comparison$confidence, useNA = "ifany"))

############################################################################
# 6. Alpha diversity
############################################################################

cat("\n=== Alpha diversity ===\n")

alpha <- data.frame(
  sample   = colnames(mat_filt),
  Richness = colSums(mat_filt > 0),
  Shannon  = apply(mat_filt, 2, function(x) vegan::diversity(x, index = "shannon")),
  stringsAsFactors = FALSE
)
alpha$disease <- factor(grp[alpha$sample], levels = c(CTRL, CASE))

# Unadjusted Wilcoxon
alpha_stats <- lapply(c("Richness", "Shannon"), function(metric) {
  ctrl <- alpha[[metric]][alpha$disease == CTRL]
  case <- alpha[[metric]][alpha$disease == CASE]
  wt   <- safe_wilcox(case, ctrl)
  if (is.null(wt)) return(data.frame(metric = metric, U = NA, p = NA,
                                     median_ctrl = safe_median(ctrl),
                                     median_case = safe_median(case)))
  data.frame(metric      = metric,
             U           = as.numeric(wt$statistic),
             p           = wt$p.value,
             median_ctrl = safe_median(ctrl),
             median_case = safe_median(case))
}) |>
  bind_rows() |>
  mutate(p.adj = p.adjust(p, method = "BH"))

print(alpha_stats)
write_tsv(alpha_stats,
          file.path(tsv_dir, paste0("alpha_diversity_deeparg_", DATASET, ".tsv")))

# Location-adjusted LM
alpha_lm_stats <- lapply(c("Richness", "Shannon"), function(metric) {
  df <- alpha
  
  if (!is.null(LOC_COL) && LOC_COL %in% colnames(meta_df) &&
      length(unique(na.omit(meta_df[[LOC_COL]]))) > 1) {
    df[[LOC_COL]] <- meta_df[df$sample, LOC_COL]
  }
  
  form <- if (!is.null(LOC_COL) && LOC_COL %in% colnames(df) &&
              length(unique(na.omit(df[[LOC_COL]]))) > 1) {
    as.formula(paste(metric, "~ disease +", LOC_COL))
  } else {
    as.formula(paste(metric, "~ disease"))
  }
  
  fit <- lm(form, data = df)
  broom::tidy(fit) |>
    filter(term == paste0("disease", CASE)) |>
    mutate(metric = metric)
}) |>
  bind_rows() |>
  mutate(q = p.adjust(p.value, method = "BH"))

cat("\nAlpha diversity LM (location-adjusted where applicable):\n")
print(alpha_lm_stats)
write_tsv(alpha_lm_stats,
          file.path(tsv_dir, paste0("alpha_diversity_lm_deeparg_", DATASET, ".tsv")))

p_alpha <- alpha |>
  pivot_longer(c(Richness, Shannon), names_to = "metric", values_to = "value") |>
  mutate(disease = factor(disease, levels = c(CTRL, CASE))) |>
  ggplot(aes(x = disease, y = value)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(aes(color = disease), width = 0.15, size = 2, alpha = 0.7) +
  facet_wrap(~metric, scales = "free_y") +
  scale_color_viridis_d(option = "plasma", end = 0.8) +
  labs(x = NULL, y = "Value",
       title = paste("ARG", LEVEL, "diversity —", DATASET)) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

ggsave(file.path(plot_dir, paste0("alpha_diversity_deeparg_", DATASET, ".png")),
       p_alpha, width = 8, height = 5, dpi = 300)

############################################################################
# 7. Beta diversity (location modelled as explanatory variable)
############################################################################

cat("\n=== Beta diversity ===\n")

dist_obj <- vegdist(t(mat_filt), method = "bray")
pcoa     <- cmdscale(dist_obj, k = 2, eig = TRUE)
eig      <- pcoa$eig
rel_eig  <- eig / sum(eig[eig > 0]) * 100

coords           <- as.data.frame(pcoa$points)
colnames(coords) <- c("Axis1", "Axis2")
coords$sample    <- rownames(coords)
coords$disease   <- factor(grp[coords$sample], levels = c(CTRL, CASE))

use_location_beta <- !is.null(LOC_COL) &&
  LOC_COL %in% colnames(meta_df) &&
  length(unique(na.omit(meta_df[[LOC_COL]]))) > 1

perm <- if (use_location_beta) {
  adonis2(
    dist_obj ~ location + disease,
    data         = meta_df |> mutate(location = factor(.data[[LOC_COL]])),
    permutations = 999,
    by           = "margin"
  )
} else {
  adonis2(dist_obj ~ disease, data = meta_df, permutations = 999, by = "margin")
}

cat("\nPERMANOVA (Bray-Curtis):\n"); print(perm)

bd       <- betadisper(dist_obj, meta_df$disease)
bd_anova <- anova(bd)

write_tsv(as.data.frame(perm) |> rownames_to_column("term"),
          file.path(tsv_dir, paste0("permanova_bray_deeparg_", DATASET, ".tsv")))
write_tsv(broom::tidy(bd_anova),
          file.path(tsv_dir, paste0("betadisper_deeparg_", DATASET, ".tsv")))

p_pcoa <- ggplot(coords, aes(x = Axis1, y = Axis2, color = disease)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_viridis_d(option = "plasma", end = 0.8) +
  labs(title = paste0("PCoA Bray-Curtis (ARG ", LEVEL, ") — ", DATASET),
       x     = paste0("PCo1 (", round(rel_eig[1], 1), "%)"),
       y     = paste0("PCo2 (", round(rel_eig[2], 1), "%)")) +
  theme_minimal(base_size = 14)

if (min(table(grp)) >= 3)
  p_pcoa <- p_pcoa + stat_ellipse(linetype = 2, linewidth = 1)

ggsave(file.path(plot_dir, paste0("pcoa_bray_deeparg_", LEVEL, "_", DATASET, ".png")),
       p_pcoa, width = 8, height = 6, dpi = 300)

saveRDS(
  list(dataset      = DATASET,
       level        = LEVEL,
       coords       = coords,
       rel_eig      = rel_eig,
       permanova    = as.data.frame(perm) |> rownames_to_column("term"),
       betadisper_p = bd_anova[["Pr(>F)"]][1],
       betadisper_F = bd_anova[["F value"]][1]),
  file.path(tsv_dir, paste0("beta_panel_deeparg_", LEVEL, "_", DATASET, ".rds"))
)

############################################################################
# 8. MaAsLin2 bar plot — top 40 by |coef|
############################################################################

if (length(maaslin_sig) > 0) {
  p_bar <- res_case |>
    slice_max(abs(case_oriented_coef), n = 40) |>
    mutate(feature = fct_reorder(feature, case_oriented_coef)) |>
    ggplot(aes(x = feature, y = case_oriented_coef, fill = direction_case)) +
    geom_col(width = 0.8) +
    coord_flip() +
    scale_fill_manual(values = c(enriched = "#d62728", depleted = "#1f77b4")) +
    labs(title    = paste("MaAsLin2 deepARG", LEVEL, "—", DATASET,
                          "(", COMPARISON, "vs", CTRL, ") — top 40 by |coef|"),
         subtitle = paste("Total significant (q<", Q_THRESH, "):", nrow(res_case)),
         x = NULL, y = "Case-oriented coefficient") +
    theme_minimal(base_size = 14) +
    theme(legend.position = "none")
  
  ggsave(file.path(plot_dir,
                   paste0("maaslin2_barplot_deeparg_", LEVEL, "_", DATASET, ".png")),
         p_bar, width = 10, height = 14, dpi = 300)
}

############################################################################
# 9. Consensus heatmap — high-confidence features
############################################################################

hc_features <- full_comparison |>
  filter(confidence == "high_confidence") |>
  pull(feature)

if (length(hc_features) > 1) {
  heatmap_z <- t(scale(t(log1p(mat_filt[hc_features, , drop = FALSE]))))
  
  ann_col <- data.frame(
    Condition = factor(grp[colnames(heatmap_z)], levels = c(CTRL, CASE)),
    row.names = colnames(heatmap_z)
  )
  
  pheatmap(
    heatmap_z,
    annotation_col = ann_col,
    show_colnames  = FALSE,
    cluster_rows   = TRUE,
    cluster_cols   = TRUE,
    fontsize_row   = 9,
    main     = paste("High-confidence DA ARG", LEVEL, "—", DATASET),
    filename = file.path(plot_dir,
                         paste0("heatmap_hc_deeparg_", LEVEL, "_", DATASET, ".png")),
    width  = 12,
    height = max(4, length(hc_features) * 0.4)
  )
}

############################################################################
# 10. Summary
############################################################################

cat("\n=== Pipeline complete:", DATASET, "(", LEVEL, ") ===\n")
cat("MaAsLin2 significant (q<", Q_THRESH, "): ", length(maaslin_sig), "\n", sep = "")
cat("Wilcoxon significant (q<", Q_THRESH, "): ", length(wilcox_sig),  "\n", sep = "")
cat("ALDEx2 significant   (q<", Q_THRESH, "): ", length(aldex_sig),   "\n", sep = "")
cat("Fisher significant   (q<", Q_THRESH, "): ", length(fisher_sig),  "\n", sep = "")
cat("High-confidence ARGs:",
    sum(full_comparison$confidence == "high_confidence", na.rm = TRUE), "\n")
cat("Plots  :", plot_dir,    "\n")
cat("Tables :", tsv_dir,     "\n")
cat("MaAsLin2:", maaslin_dir, "\n")
