## =============================================================================
## KrakenAnalysis_1.R
## CRASH — Cancer-associated gut micRobiome Analysis of Strains in
##         Haematological patients
##
## Purpose:
##   Per-cohort differential abundance and diversity analysis using
##   Kraken2/Bracken taxonomic profiles. Runs in parallel to MetaAnalysis_1.R
##   as an independent classifier validation layer.
##
##   Key differences from the MetaPhlAn pipeline:
##     - Input: raw Bracken integer counts (not relative abundances)
##     - MaAsLin2: TSS + LOG normalisation (appropriate for count data)
##     - ALDEx2: receives raw counts directly (primary method, not sensitivity)
##     - Prevalence filter: 30% (Kraken detects more low-abundance taxa)
##     - Feature name sanitisation: normalize_name() key handles bracket taxa
##       (e.g. "[Clostridium]") and MaAsLin2's internal make.names() conversion
##
## Usage:
##   Set DATASET (line 43) to one of:
##     "kulecka_leukemia" | "kulecka_lymphoma" | "cra" | "prjna813705"
##   then source the script.
##
## Outputs (per cohort, dataset-namespaced subdirs):
##   plots_kraken_<DATASET>/         PNG figures
##   tables_kraken_<DATASET>/        TSV result tables + RDS beta payloads
##   maaslin2_kraken_<DATASET>/      MaAsLin2 full output folder
##
## Dependencies:
##   readxl, janitor, tidyverse, phyloseq, Maaslin2, ALDEx2,
##   ggpubr, forcats, vegan, viridis, pairwiseAdonis, broom
##
## Notes:
##   - select <- dplyr::select set explicitly to prevent namespace masking
##   - Beta diversity RDS payloads feed KrakenAnalysis_3.R cross-cohort panel
##   - Kraken2 v2.1.3 + Bracken v2.8 (PlusPF database)
## =============================================================================

############################################################################
# CONFIGURATION
############################################################################

DATASET <- "kulecka_lymphoma"   # "kulecka_leukemia", "kulecka_lymphoma", "cra", "prjna813705"

cfg <- list(
  kulecka_leukemia = list(
    wdir         = "D:/MasterThesis/Vik/Kulecka/Leukemia",
    meta_file    = "D:/MasterThesis/Vik/Kulecka/Leukemia/Kuleckaleukemia_metadata.xlsx",
    bracken_file = "D:/MasterThesis/Vik/Kulecka/Leukemia/KRAKEN2/bracken_species_counts.tsv",
    case_label   = "AML",    ctrl_label = "Control",
    comparison   = "AML",    location_col = "location"
  ),
  kulecka_lymphoma = list(
    wdir         = "D:/MasterThesis/Vik/Kulecka/Lymphoma",
    meta_file    = "D:/MasterThesis/Vik/Kulecka/Lymphoma/Kuleckalymphoma_metadata.xlsx",
    bracken_file = "D:/MasterThesis/Vik/Kulecka/Lymphoma/KRAKEN2/bracken_species_counts.tsv",
    case_label   = "LN",     ctrl_label = "Control",
    comparison   = "LN",     location_col = "location"
  ),
  cra = list(
    wdir         = "D:/MasterThesis/Vik/CRA007433",
    meta_file    = "D:/MasterThesis/Vik/CRA007433/CRA007433_metadata.xlsx",
    bracken_file = "D:/MasterThesis/Vik/CRA007433/KRAKEN2/bracken_species_counts.tsv",
    case_label   = "NKTCL",  ctrl_label = "Control",
    comparison   = "NKTCL",  location_col = NULL
  ),
  prjna813705 = list(
    wdir         = "D:/MasterThesis/Vik/PRJNA813705",
    meta_file    = "D:/MasterThesis/Vik/PRJNA813705/PRJNA813705_metadata.xlsx",
    bracken_file = "D:/MasterThesis/Vik/PRJNA813705/KRAKEN2/bracken_species_counts.tsv",
    case_label   = "AML",    ctrl_label = "Control",
    comparison   = "AML",    location_col = NULL
  )
)

stopifnot(DATASET %in% names(cfg))

C          <- cfg[[DATASET]]
wdir       <- C$wdir
CASE       <- C$case_label
CTRL       <- C$ctrl_label
COMPARISON <- C$comparison
LOC_COL    <- C$location_col

Q_THRESH        <- 0.25
TOP_N           <- 5
MIN_PREVALENCE  <- 0.30
MIN_MEAN_AB     <- 1e-4
SEED            <- 42

############################################################################
# Packages
############################################################################

library(readxl)
library(janitor)
library(tidyverse)
library(phyloseq)
library(Maaslin2)
library(ALDEx2)
library(ggpubr)
library(forcats)
library(vegan)
library(viridis)
library(pairwiseAdonis)
library(broom)

select <- dplyr::select

############################################################################
# Output directories
############################################################################

plot_dir    <- file.path(wdir, paste0("plots_kraken_", DATASET))
tsv_dir     <- file.path(wdir, paste0("tables_kraken_", DATASET))
maaslin_dir <- file.path(wdir, paste0("maaslin2_kraken_", DATASET))

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

# Language-invariant key used to join MaAsLin2 results (where make.names()
# has mangled species names) with Wilcoxon / ALDEx2 results (which keep
# the original names).
normalize_name <- function(x) {
  x <- tolower(x)
  x <- gsub("^x\\.+", "", x)
  x <- gsub("[^a-z0-9]", "", x)
  x
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
# Load Bracken species counts
############################################################################

bracken <- read.delim(C$bracken_file, row.names = 1, check.names = FALSE)
colnames(bracken) <- tidy_ids(colnames(bracken))

cat("\nBracken table:", nrow(bracken), "species x", ncol(bracken), "samples\n")

valid_samples <- intersect(colnames(bracken), meta$sample_id)
bracken <- bracken[, valid_samples, drop = FALSE]
cat("After metadata intersection:", ncol(bracken), "samples\n")

meta <- meta |> filter(sample_id %in% valid_samples)
meta <- meta[match(colnames(bracken), meta$sample_id), ]
stopifnot(identical(colnames(bracken), meta$sample_id))

############################################################################
# Prevalence / abundance filter
############################################################################

relab_full <- sweep(as.matrix(bracken), 2, colSums(bracken), "/")
prev       <- calc_prevalence(relab_full)
mean_ab    <- rowMeans(relab_full, na.rm = TRUE)
keep       <- prev >= MIN_PREVALENCE & mean_ab >= MIN_MEAN_AB

cat("\nFeatures before filter:", nrow(bracken), "\n")
cat("Features after filter:", sum(keep), "\n")

bracken_filt <- bracken[keep, , drop = FALSE]
relab_filt   <- relab_full[keep, , drop = FALSE]

grp <- as.character(meta$disease)
names(grp) <- meta$sample_id

############################################################################
# Build phyloseq (for diversity, not DA)
############################################################################

tax_mat <- matrix(rownames(bracken_filt), ncol = 1,
                  dimnames = list(rownames(bracken_filt), "species"))

meta_ps <- as.data.frame(meta)
rownames(meta_ps) <- meta$sample_id

ps <- phyloseq(
  otu_table(as.matrix(bracken_filt), taxa_are_rows = TRUE),
  tax_table(tax_mat),
  sample_data(meta_ps)
)

cat("\nPhyloseq:", nsamples(ps), "samples,", ntaxa(ps), "taxa\n")
print(table(sample_data(ps)$disease, useNA = "ifany"))

############################################################################
# 1. MaAsLin2 — TSS + LOG on raw counts
############################################################################

cat("\n=== MaAsLin2 ===\n")

input_maaslin <- as.data.frame(t(bracken_filt))
metadata_maaslin <- data.frame(
  disease   = factor(as.character(meta$disease), levels = c(CTRL, CASE)),
  row.names = meta$sample_id
)
metadata_maaslin <- metadata_maaslin[rownames(input_maaslin), , drop = FALSE]

cat("MaAsLin2 input:", nrow(input_maaslin), "samples x", ncol(input_maaslin), "features\n")

fit_data <- Maaslin2(
  input_data       = input_maaslin,
  input_metadata   = metadata_maaslin,
  output           = maaslin_dir,
  fixed_effects    = "disease",
  normalization    = "TSS",
  transform        = "LOG",
  min_prevalence   = 0.0,       # already filtered above
  min_abundance    = 0.0,
  analysis_method  = "LM",
  correction       = "BH",
  standardize      = FALSE,
  max_significance = Q_THRESH
)

results <- read_tsv(file.path(maaslin_dir, "all_results.tsv"),
                    show_col_types = FALSE)

res_case <- results |>
  filter(metadata == "disease") |>
  mutate(
    case_oriented_coef = case_when(
      value == CASE ~ coef,
      value == CTRL ~ -coef,
      TRUE ~ NA_real_
    ),
    direction_case = case_when(
      case_oriented_coef > 0 ~ "enriched",
      case_oriented_coef < 0 ~ "depleted",
      TRUE ~ NA_character_
    ),
    norm_key = normalize_name(feature)
  ) |>
  filter(qval < Q_THRESH)

maaslin_sig_norm <- res_case$norm_key
cat("MaAsLin2 significant hits (q <", Q_THRESH, "):", nrow(res_case), "\n")

############################################################################
# 2. Wilcoxon rank-sum on relative abundances
############################################################################

cat("\n=== Wilcoxon rank-sum ===\n")

wilcox_results <- lapply(rownames(relab_filt), function(feat) {
  vals <- relab_filt[feat, ]
  ctrl <- vals[grp == CTRL]
  case <- vals[grp == CASE]
  
  wt <- safe_wilcox(case, ctrl)   # case first -> large U = case-enriched
  if (is.null(wt)) return(NULL)
  
  n_case <- length(case)
  n_ctrl <- length(ctrl)
  U      <- as.numeric(wt$statistic)
  
  data.frame(
    feature       = feat,
    prevalence    = prev[feat],
    U             = U,
    p             = wt$p.value,
    rank_biserial = (2 * U) / (n_case * n_ctrl) - 1,
    median_ctrl   = safe_median(ctrl),
    median_case   = safe_median(case),
    log2FC        = log2((safe_median(case) + 1e-10) / (safe_median(ctrl) + 1e-10)),
    stringsAsFactors = FALSE
  )
}) |>
  bind_rows() |>
  mutate(
    p.adj = p.adjust(p, method = "BH"),
    wilcox_direction = case_when(
      log2FC > 0 ~ "enriched",
      log2FC < 0 ~ "depleted",
      TRUE ~ NA_character_
    ),
    norm_key = normalize_name(feature)
  ) |>
  arrange(p)

wilcox_sig_norm <- wilcox_results |> filter(p.adj < Q_THRESH) |> pull(norm_key)

cat("Features tested:", nrow(wilcox_results), "\n")
cat("Wilcoxon significant (q <", Q_THRESH, "):", length(wilcox_sig_norm), "\n")
cat("MaAsLin2-Wilcoxon overlap:", length(intersect(maaslin_sig_norm, wilcox_sig_norm)), "\n")

write_tsv(
  wilcox_results,
  file.path(tsv_dir, paste0("wilcoxon_kraken_", DATASET, ".tsv"))
)

############################################################################
# 3. ALDEx2 — primary method on raw counts
############################################################################

cat("\n=== ALDEx2 ===\n")

counts_mat <- as.matrix(bracken_filt)
mode(counts_mat) <- "integer"
conditions <- grp[colnames(counts_mat)]

cat("Input:", nrow(counts_mat), "features x", ncol(counts_mat), "samples\n")

set.seed(SEED)
aldex_res <- aldex(
  reads      = counts_mat,
  conditions = conditions,
  mc.samples = 128,
  test       = "t",
  effect     = TRUE,
  denom      = "all",
  verbose    = FALSE
) |>
  rownames_to_column("feature") |>
  mutate(
    aldex_direction = case_when(
      effect > 0 & CASE == levels(meta$disease)[2] ~ "enriched",
      effect < 0 & CASE == levels(meta$disease)[2] ~ "depleted",
      TRUE ~ NA_character_
    ),
    norm_key = normalize_name(feature)
  )

aldex_sig_norm <- aldex_res |> filter(wi.eBH < Q_THRESH) |> pull(norm_key)
cat("ALDEx2 significant (wi.eBH <", Q_THRESH, "):", length(aldex_sig_norm), "\n")

write_tsv(
  aldex_res,
  file.path(tsv_dir, paste0("aldex2_kraken_", DATASET, ".tsv"))
)

############################################################################
# 4. Cross-method comparison and confidence tiers
############################################################################

full_comparison <- wilcox_results |>
  select(feature, norm_key, prevalence,
         wilcox_q = p.adj, rank_biserial, log2FC, wilcox_direction) |>
  full_join(
    res_case |>
      select(norm_key,
             maaslin_value     = value,
             maaslin_coef_raw  = coef,
             maaslin_coef      = case_oriented_coef,
             maaslin_q         = qval,
             maaslin_direction = direction_case),
    by = "norm_key"
  ) |>
  full_join(
    aldex_res |>
      select(norm_key,
             aldex_q      = wi.eBH,
             aldex_effect = effect,
             aldex_direction),
    by = "norm_key"
  ) |>
  mutate(
    wilcox_sig  = !is.na(wilcox_q)  & wilcox_q  < Q_THRESH,
    maaslin_sig = !is.na(maaslin_q) & maaslin_q < Q_THRESH,
    aldex_sig   = !is.na(aldex_q)   & aldex_q   < Q_THRESH,
    n_sig = rowSums(cbind(wilcox_sig, maaslin_sig, aldex_sig), na.rm = TRUE),
    
    direction_consensus = case_when(
      !is.na(maaslin_direction) & !is.na(wilcox_direction) &
        maaslin_direction == wilcox_direction ~ maaslin_direction,
      !is.na(maaslin_direction) & !is.na(aldex_direction) &
        maaslin_direction == aldex_direction ~ maaslin_direction,
      !is.na(wilcox_direction) & !is.na(aldex_direction) &
        wilcox_direction == aldex_direction ~ wilcox_direction,
      TRUE ~ NA_character_
    ),
    
    confidence = case_when(
      maaslin_sig & n_sig >= 2 & !is.na(direction_consensus) ~ "high_confidence",
      n_sig >= 2 & !is.na(direction_consensus)               ~ "moderate_confidence",
      n_sig == 1                                             ~ "single_method",
      TRUE                                                   ~ "not_significant"
    )
  ) |>
  arrange(maaslin_q, wilcox_q, aldex_q)

# Drop the join_key from the final output
full_comparison <- full_comparison |> select(-norm_key)

write_tsv(
  full_comparison,
  file.path(tsv_dir, paste0("cross_method_kraken_", DATASET, ".tsv"))
)

############################################################################
# 5. Alpha diversity (on counts — native for Kraken)
############################################################################

cat("\n=== Alpha diversity ===\n")

alpha <- estimate_richness(ps, measures = c("Observed", "Shannon", "Simpson"))
alpha$sample <- rownames(alpha)
alpha$disease <- as.character(sample_data(ps)[alpha$sample, "disease"]$disease)

alpha_stats <- lapply(c("Observed", "Shannon", "Simpson"), function(metric) {
  ctrl <- alpha[[metric]][alpha$disease == CTRL]
  case <- alpha[[metric]][alpha$disease == CASE]
  
  wt <- safe_wilcox(case, ctrl)
  if (is.null(wt)) {
    return(data.frame(
      metric = metric, U = NA_real_, p = NA_real_,
      median_ctrl = safe_median(ctrl), median_case = safe_median(case)
    ))
  }
  data.frame(
    metric      = metric,
    U           = as.numeric(wt$statistic),
    p           = wt$p.value,
    median_ctrl = safe_median(ctrl),
    median_case = safe_median(case)
  )
}) |>
  bind_rows() |>
  mutate(p.adj = p.adjust(p, method = "BH"))

print(alpha_stats)

write_tsv(
  alpha_stats,
  file.path(tsv_dir, paste0("alpha_diversity_kraken_", DATASET, ".tsv"))
)

p_alpha <- alpha |>
  pivot_longer(c(Observed, Shannon, Simpson), names_to = "metric", values_to = "value") |>
  mutate(disease = factor(disease, levels = c(CTRL, CASE))) |>
  ggplot(aes(x = disease, y = value)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(aes(color = disease), width = 0.15, size = 2, alpha = 0.7) +
  facet_wrap(~ metric, scales = "free_y") +
  scale_color_viridis_d(option = "plasma", end = 0.8) +
  labs(x = NULL, y = "Value",
       title = paste("Alpha diversity (Kraken) —", DATASET)) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

ggsave(
  file.path(plot_dir, paste0("alpha_diversity_kraken_", DATASET, ".png")),
  p_alpha, width = 12, height = 5, dpi = 300
)

############################################################################
# 6. Beta diversity — Bray-Curtis & Jaccard (on relative abundances)
############################################################################

cat("\n=== Beta diversity ===\n")

ps_rel <- transform_sample_counts(ps, function(x) {
  s <- sum(x); if (s == 0) return(x); x / s
})

for (metric in c("bray", "jaccard")) {
  binary <- metric == "jaccard"
  dist_obj <- phyloseq::distance(ps_rel, method = metric, binary = binary)
  ord <- ordinate(ps_rel, method = "PCoA", distance = dist_obj)
  
  md <- data.frame(sample_data(ps_rel))[labels(dist_obj), , drop = FALSE]
  md$disease <- factor(md$disease, levels = c(CTRL, CASE))
  
  use_strata <- !is.null(LOC_COL) &&
    LOC_COL %in% colnames(md) &&
    length(unique(na.omit(md[[LOC_COL]]))) > 1
  
  if (use_strata) {
    perm <- adonis2(dist_obj ~ disease, data = md, permutations = 999,
                    by = "margin", strata = md[[LOC_COL]])
  } else {
    perm <- adonis2(dist_obj ~ disease, data = md, permutations = 999, by = "margin")
  }
  
  cat("\nPERMANOVA (", metric, "):\n", sep = ""); print(perm)
  
  bd <- betadisper(dist_obj, md$disease)
  bd_anova <- anova(bd)
  
  write_tsv(
    as.data.frame(perm) |> rownames_to_column("term"),
    file.path(tsv_dir, paste0("permanova_kraken_", metric, "_", DATASET, ".tsv"))
  )
  write_tsv(
    broom::tidy(bd_anova),
    file.path(tsv_dir, paste0("betadisper_kraken_anova_", metric, "_", DATASET, ".tsv"))
  )
  
  p_pcoa <- plot_ordination(ps_rel, ord, color = "disease") +
    geom_point(size = 3, alpha = 0.8) +
    scale_color_viridis_d(option = "plasma", end = 0.8) +
    labs(title = paste0("PCoA (", metric, ") Kraken — ", DATASET)) +
    theme_minimal(base_size = 14)
  
  if (min(table(md$disease)) >= 3) {
    p_pcoa <- p_pcoa + stat_ellipse(linetype = 2, linewidth = 1)
  }
  
  ggsave(
    file.path(plot_dir, paste0("pcoa_kraken_", metric, "_", DATASET, ".png")),
    p_pcoa, width = 8, height = 6, dpi = 300
  )
  
  png(
    file.path(plot_dir, paste0("betadisper_kraken_", metric, "_", DATASET, ".png")),
    width = 8, height = 6, units = "in", res = 300
  )
  boxplot(bd, col = viridis(2, end = 0.8),
          ylab = "Distance to centroid",
          main = paste0("Betadisper (", metric, ") Kraken — ", DATASET))
  dev.off()
}

############################################################################
# 7. MaAsLin2 bar plot
############################################################################

if (nrow(res_case) > 0) {
  top_hits <- res_case |>
    arrange(desc(abs(case_oriented_coef))) |>
    slice_head(n = TOP_N * 2)
  
  p_bar <- top_hits |>
    mutate(feature = fct_reorder(feature, case_oriented_coef)) |>
    ggplot(aes(x = feature, y = case_oriented_coef, fill = direction_case)) +
    geom_col(width = 0.8) +
    coord_flip() +
    scale_fill_manual(values = c(enriched = "#d62728", depleted = "#1f77b4")) +
    labs(
      title = paste("MaAsLin2 (Kraken) —", DATASET, "(", COMPARISON, "vs", CTRL, ")"),
      x = NULL, y = "Case-oriented effect size (coef)"
    ) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "none")
  
  ggsave(
    file.path(plot_dir, paste0("maaslin2_kraken_barplot_", DATASET, ".png")),
    p_bar, width = 10, height = 6, dpi = 300
  )
}

############################################################################
# 8. Summary
############################################################################

cat("\n=== Kraken pipeline complete:", DATASET, "===\n")
cat("MaAsLin2 significant (q<", Q_THRESH, "): ", length(maaslin_sig_norm), "\n", sep = "")
cat("Wilcoxon significant (q<", Q_THRESH, "): ", length(wilcox_sig_norm), "\n", sep = "")
cat("ALDEx2 significant  (q<", Q_THRESH, "): ", length(aldex_sig_norm), "\n", sep = "")
cat("MaAsLin2-Wilcoxon overlap:", length(intersect(maaslin_sig_norm, wilcox_sig_norm)), "\n")
cat("MaAsLin2-ALDEx2 overlap:",   length(intersect(maaslin_sig_norm, aldex_sig_norm)),   "\n")
cat(
  "High-confidence taxa:",
  sum(full_comparison$confidence == "high_confidence", na.rm = TRUE), "\n"
)
cat("Plots:", plot_dir, "\n")
cat("Tables:", tsv_dir, "\n")
