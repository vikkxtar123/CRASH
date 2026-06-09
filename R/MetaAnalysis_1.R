## =============================================================================
## MetaAnalysis_1.R
## CRASH — Cancer-associated gut micRobiome Analysis of Strains in
##         Haematological patients
##
## Purpose:
##   Per-cohort differential abundance and diversity analysis of MetaPhlAn4
##   taxonomic profiles across four haematological cancer cohorts:
##     - Belgian AML      (PRJNA813705)
##     - Polish AML       (Kulecka_leukemia / PRJNA1116523)
##     - Polish Lymphoma  (Kulecka_lymphoma)
##     - Chinese NKTCL    (CRA007433)
##
##   For each cohort the script runs:
##     1. MaAsLin2      — primary linear model DA (q < 0.25, BH)
##     2. Wilcoxon      — non-parametric support; rank-biserial effect size
##     3. ALDEx2        — compositional sensitivity analysis (pseudo-counts)
##     4. Cross-method  — confidence tiers (high / moderate / single_method)
##     5. Alpha diversity — Observed, Shannon, Simpson; Wilcoxon + BH
##     6. Beta diversity  — Bray-Curtis & Jaccard PCoA; PERMANOVA (location-
##                          stratified for Kulecka cohorts); betadisper
##     7. MaAsLin2 bar plot of top hits
##
## Usage:
##   Set DATASET (line 51) to one of:
##     "kulecka_leukemia" | "kulecka_lymphoma" | "cra" | "prjna813705"
##   then source the script. Outputs are written to dataset-namespaced
##   subdirectories under the cohort working directory defined in cfg[].
##
## Outputs (per cohort):
##   plots_metaphlan_<DATASET>/        PNG figures
##   tables_metaphlan_<DATASET>/       TSV result tables + RDS beta payloads
##   maaslin2_metaphlan_<DATASET>/     MaAsLin2 full output folder
##
## Dependencies:
##   readxl, janitor, tidyverse, phyloseq, Maaslin2, ALDEx2,
##   ggpubr, forcats, vegan, viridis, pairwiseAdonis, broom
##
## Notes:
##   - All samples are pre-treatment / drug-naive
##   - select <- dplyr::select is set explicitly to prevent namespace masking
##   - Wilcoxon called case-first so positive rank-biserial = case-enriched
##   - ALDEx2 uses scaled pseudo-counts (×1e6); treat as sensitivity only
##   - Beta diversity RDS payloads feed MetaAnalysis_3.R cross-cohort panel
## =============================================================================

############################################################################
# CONFIGURATION
############################################################################

DATASET <- "kulecka_leukemia"   # "kulecka_leukemia", "kulecka_lymphoma", "cra", "prjna813705"

cfg <- list(
  kulecka_leukemia = list(
    wdir         = "D:/MasterThesis/Vik/Kulecka/Leukemia",
    meta_file    = "D:/MasterThesis/Vik/Kulecka/Leukemia/Kuleckaleukemia_metadata.xlsx",
    profile_file = "D:/MasterThesis/Vik/Kulecka/Leukemia/leukemia_merged_profiles.tsv",
    case_label   = "AML",  ctrl_label = "Control",
    comparison   = "AML",  location_col = "location"
  ),
  kulecka_lymphoma = list(
    wdir         = "D:/MasterThesis/Vik/Kulecka/Lymphoma",
    meta_file    = "D:/MasterThesis/Vik/Kulecka/Lymphoma/Kuleckalymphoma_metadata.xlsx",
    profile_file = "D:/MasterThesis/Vik/Kulecka/Lymphoma/lymphoma_merged_profiles.tsv",
    case_label   = "LN",   ctrl_label = "Control",
    comparison   = "LN",   location_col = "location"
  ),
  cra = list(
    wdir         = "D:/MasterThesis/Vik/CRA007433",
    meta_file    = "D:/MasterThesis/Vik/CRA007433/CRA007433_metadata.xlsx",
    profile_file = "D:/MasterThesis/Vik/CRA007433/cra_merged_profiles.tsv",
    case_label   = "NKTCL", ctrl_label = "Control",
    comparison   = "NKTCL", location_col = NULL
  ),
  prjna813705 = list(
    wdir         = "D:/MasterThesis/Vik/PRJNA813705",
    meta_file    = "D:/MasterThesis/Vik/PRJNA813705/PRJNA813705_metadata.xlsx",
    profile_file = "D:/MasterThesis/Vik/PRJNA813705/merged_abundance.txt",
    case_label   = "AML",  ctrl_label = "Control",
    comparison   = "AML",  location_col = NULL
  )
)

stopifnot(DATASET %in% names(cfg))

C          <- cfg[[DATASET]]
wdir       <- C$wdir
CASE       <- C$case_label
CTRL       <- C$ctrl_label
COMPARISON <- C$comparison
LOC_COL    <- C$location_col

Q_THRESH         <- 0.25
TOP_N            <- 5
RICHNESS_THRESH  <- 1e-5
MIN_PREVALENCE   <- 0.10
ALDEX_SCALE      <- 1e6
SEED             <- 42

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

plot_dir    <- file.path(wdir, paste0("plots_metaphlan_", DATASET))
tsv_dir     <- file.path(wdir, paste0("tables_metaphlan_", DATASET))
maaslin_dir <- file.path(wdir, paste0("maaslin2_metaphlan_", DATASET))

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

read_metaphlan_merged <- function(path) {
  read_tsv(path, comment = "#", show_col_types = FALSE) |>
    rename(clade_name = 1) |>
    rename_with(~ make.names(., unique = TRUE))
}

keep_nonredundant_leafs <- function(df) {
  clades <- df$clade_name
  is_unclass <- toupper(clades) == "UNCLASSIFIED"
  is_t       <- grepl("(^|\\|)t__", clades)
  is_s       <- grepl("(^|\\|)s__", clades)
  
  species_prefix_from_t <- sub("(.*\\|s__[^|]+)\\|t__.*", "\\1", clades[is_t])
  species_with_t        <- unique(species_prefix_from_t)
  
  species_prefix_all <- ifelse(
    is_s,
    sub("(.*\\|s__[^|]+).*", "\\1", clades),
    NA_character_
  )
  species_has_t_child <- !is.na(species_prefix_all) & species_prefix_all %in% species_with_t
  
  df[is_unclass | is_t | (is_s & !species_has_t_child), , drop = FALSE]
}

normalize_to_100 <- function(df) {
  nums <- df |>
    select(-clade_name) |>
    mutate(across(everything(), ~ suppressWarnings(as.numeric(.))))
  nums[is.na(nums)] <- 0
  s <- colSums(nums); s[s == 0] <- 1
  scaled <- sweep(as.matrix(nums), 2, s, "/") * 100
  bind_cols(df["clade_name"], as.data.frame(scaled, check.names = FALSE))
}

parse_tax <- function(clades, ref_rownames) {
  toks    <- strsplit(clades, "\\|")
  max_len <- max(lengths(toks))
  pad <- lapply(toks, function(x) { length(x) <- max_len; x })
  tax_df <- as.data.frame(do.call(rbind, pad), stringsAsFactors = FALSE)
  ranks <- c("kingdom", "phylum", "class", "order", "family", "genus", "species", "subspecies")
  names(tax_df) <- head(ranks, ncol(tax_df))
  tax_df[] <- lapply(tax_df, function(col) {
    col[is.na(col)] <- ""
    col <- sub("^[a-z]__", "", col)
    ifelse(col == "", NA, col)
  })
  rownames(tax_df) <- ref_rownames
  as.matrix(tax_df)
}

make_safe_feature_names <- function(ps_obj, fallback_prefix = "taxon") {
  tax_df <- as.data.frame(tax_table(ps_obj), stringsAsFactors = FALSE)
  feat_names <- tax_df$species
  bad <- is.na(feat_names) | feat_names == ""
  feat_names[bad] <- taxa_names(ps_obj)[bad]
  feat_names[is.na(feat_names) | feat_names == ""] <-
    paste0(fallback_prefix, "_", seq_len(sum(is.na(feat_names) | feat_names == "")))
  make.unique(feat_names)
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
# Load and process MetaPhlAn profiles
############################################################################

prof <- read_metaphlan_merged(C$profile_file) |>
  keep_nonredundant_leafs()

colnames(prof)[-1] <- tidy_ids(colnames(prof)[-1])

valid_samples <- intersect(colnames(prof)[-1], meta$sample_id)
prof <- prof |> select(clade_name, all_of(valid_samples))

cat("Profiles:", ncol(prof) - 1, "samples after metadata intersection\n")

############################################################################
# Normalize and build phyloseq
############################################################################

prof_noU <- prof |> filter(toupper(clade_name) != "UNCLASSIFIED")

prof_norm <- normalize_to_100(prof_noU) |>
  group_by(clade_name) |>
  summarise(
    across(everything(), ~ sum(suppressWarnings(as.numeric(.)), na.rm = TRUE)),
    .groups = "drop"
  )

otu_mat <- prof_norm |> column_to_rownames("clade_name") |> as.matrix()
mode(otu_mat) <- "numeric"
otu_mat <- otu_mat / 100

tax <- tax_table(parse_tax(rownames(otu_mat), rownames(otu_mat)))

meta_ps <- meta |>
  filter(sample_id %in% colnames(otu_mat)) |>
  as.data.frame()
rownames(meta_ps) <- meta_ps$sample_id

shared <- intersect(colnames(otu_mat), rownames(meta_ps))
otu_mat <- otu_mat[, shared, drop = FALSE]
meta_ps <- meta_ps[shared, , drop = FALSE]

ps <- phyloseq(
  otu_table(otu_mat, taxa_are_rows = TRUE),
  tax_table(tax),
  sample_data(meta_ps)
)

cat("\nPhyloseq:", nsamples(ps), "samples,", ntaxa(ps), "taxa\n")
print(table(sample_data(ps)$disease, useNA = "ifany"))

############################################################################
# Species-level object
############################################################################

ps_spp <- tax_glom(ps, taxrank = "species", NArm = TRUE)
ps_spp_rel <- transform_sample_counts(ps_spp, function(x) {
  s <- sum(x); if (s == 0) return(x); x / s
})

feature_names <- make_safe_feature_names(ps_spp_rel)

relab_mat <- as(otu_table(ps_spp_rel), "matrix")
if (!taxa_are_rows(ps_spp_rel)) relab_mat <- t(relab_mat)
rownames(relab_mat) <- feature_names

grp <- as.character(sample_data(ps_spp_rel)$disease)
names(grp) <- sample_names(ps_spp_rel)

prev <- calc_prevalence(relab_mat)
keep_features <- prev >= MIN_PREVALENCE
relab_mat_filt <- relab_mat[keep_features, , drop = FALSE]
cat("\nSpecies tested after prevalence filter:", nrow(relab_mat_filt), "\n")

############################################################################
# 1. MaAsLin2
############################################################################

cat("\n=== MaAsLin2 ===\n")

input_maaslin <- as.data.frame(t(relab_mat_filt))

md_full <- data.frame(sample_data(ps_spp_rel), stringsAsFactors = FALSE)
md_full <- md_full[rownames(input_maaslin), , drop = FALSE]

has_loc <- !is.null(LOC_COL) && LOC_COL %in% colnames(md_full) &&
  length(unique(na.omit(md_full[[LOC_COL]]))) > 1

metadata_maaslin <- data.frame(
  disease   = factor(as.character(md_full$disease), levels = c(CTRL, CASE)),
  row.names = rownames(md_full)
)
fixed_fx <- "disease"
if (has_loc) {
  metadata_maaslin[[LOC_COL]] <- factor(md_full[[LOC_COL]])
  fixed_fx <- c("disease", LOC_COL)
  cat("MaAsLin2: adjusting for", LOC_COL, "\n")
}

cat("MaAsLin2 input:", nrow(input_maaslin), "samples x", ncol(input_maaslin), "features\n")

fit_data <- Maaslin2(
  input_data       = input_maaslin,
  input_metadata   = metadata_maaslin,
  output           = maaslin_dir,
  fixed_effects    = "disease",
  normalization    = "NONE",
  transform        = "LOG",
  min_prevalence   = 0.0,   # already filtered above
  analysis_method  = "LM",
  correction       = "BH",
  standardize      = FALSE, # no continuous predictors -> no-op, but explicit
  max_significance = Q_THRESH
)

results <- read_tsv(file.path(maaslin_dir, "all_results.tsv"), show_col_types = FALSE)

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
    )
  ) |>
  filter(qval < Q_THRESH)

maaslin_sig <- res_case$feature
cat("MaAsLin2 significant hits (q <", Q_THRESH, "):", nrow(res_case), "\n")

############################################################################
# 2. Wilcoxon rank-sum (case-first, so positive U = case-enriched)
############################################################################

cat("\n=== Wilcoxon rank-sum ===\n")

wilcox_results <- lapply(rownames(relab_mat_filt), function(feat) {
  vals <- relab_mat_filt[feat, ]
  ctrl <- vals[grp == CTRL]
  case <- vals[grp == CASE]
  
  # Call with case first so the returned statistic is U(case vs ctrl):
  # large U -> case stochastically greater than ctrl -> case-enriched
  wt <- safe_wilcox(case, ctrl)
  if (is.null(wt)) return(NULL)
  
  n_case <- length(case)
  n_ctrl <- length(ctrl)
  
  # R's wilcox.test already returns U; do NOT subtract n*(n+1)/2
  U <- as.numeric(wt$statistic)
  
  data.frame(
    feature       = feat,
    prevalence    = prev[feat],
    U             = U,
    p             = wt$p.value,
    rank_biserial = (2 * U) / (n_case * n_ctrl) - 1,   # in [-1, 1]; +ve = case-enriched
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
      rank_biserial > 0 ~ "enriched",
      rank_biserial < 0 ~ "depleted",
      TRUE              ~ NA_character_
    )
  ) |>
  arrange(p)

wilcox_sig <- wilcox_results |> filter(p.adj < Q_THRESH) |> pull(feature)

cat("Features tested:", nrow(wilcox_results), "\n")
cat("Wilcoxon significant (q <", Q_THRESH, "):", length(wilcox_sig), "\n")
cat("MaAsLin2-Wilcoxon overlap:", length(intersect(maaslin_sig, wilcox_sig)), "\n")

write_tsv(
  wilcox_results,
  file.path(tsv_dir, paste0("wilcoxon_metaphlan_", DATASET, ".tsv"))
)

############################################################################
# 3. ALDEx2 (supportive sensitivity analysis)
############################################################################

cat("\n=== ALDEx2 ===\n")
cat("Note: ALDEx2 uses scaled MetaPhlAn relative abundances, not raw counts.\n")
cat("Treat as sensitivity analysis, not primary.\n")

pseudo_counts <- round(relab_mat_filt * ALDEX_SCALE)
mode(pseudo_counts) <- "integer"
pseudo_counts <- pseudo_counts[rowSums(pseudo_counts) > 0, , drop = FALSE]

set.seed(SEED)
aldex_res <- aldex(
  reads      = pseudo_counts,
  conditions = grp[colnames(pseudo_counts)],
  mc.samples = 128,
  test       = "t",
  effect     = TRUE,
  denom      = "all",
  verbose    = FALSE
) |>
  rownames_to_column("feature") |>
  mutate(
    aldex_direction = case_when(
      effect > 0 ~ "enriched",
      effect < 0 ~ "depleted",
      TRUE ~ NA_character_
    )
  )

aldex_sig <- aldex_res |> filter(wi.eBH < Q_THRESH) |> pull(feature)
cat("ALDEx2 significant (wi.eBH <", Q_THRESH, "):", length(aldex_sig), "\n")

write_tsv(
  aldex_res,
  file.path(tsv_dir, paste0("aldex2_metaphlan_", DATASET, ".tsv"))
)

############################################################################
# 4. Cross-method comparison and confidence tiers
############################################################################

full_comparison <- wilcox_results |>
  select(feature, prevalence, wilcox_q = p.adj, rank_biserial, log2FC, wilcox_direction) |>
  full_join(
    res_case |>
      select(
        feature,
        maaslin_value     = value,
        maaslin_coef_raw  = coef,
        maaslin_coef      = case_oriented_coef,
        maaslin_q         = qval,
        maaslin_direction = direction_case
      ),
    by = "feature"
  ) |>
  full_join(
    aldex_res |> select(feature, aldex_q = wi.eBH, aldex_effect = effect, aldex_direction),
    by = "feature"
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

write_tsv(
  full_comparison,
  file.path(tsv_dir, paste0("cross_method_metaphlan_", DATASET, ".tsv"))
)

############################################################################
# 5. Alpha diversity
############################################################################

cat("\n=== Alpha diversity ===\n")

alpha <- data.frame(
  sample   = colnames(relab_mat),
  Observed = apply(relab_mat > RICHNESS_THRESH, 2, sum),
  Shannon  = apply(relab_mat, 2, function(x) vegan::diversity(x, index = "shannon")),
  Simpson  = apply(relab_mat, 2, function(x) vegan::diversity(x, index = "simpson")),
  stringsAsFactors = FALSE
)
alpha$disease <- as.character(sample_data(ps_spp_rel)[alpha$sample, "disease"]$disease)

alpha_stats <- lapply(c("Observed", "Shannon", "Simpson"), function(metric) {
  ctrl <- alpha[[metric]][alpha$disease == CTRL]
  case <- alpha[[metric]][alpha$disease == CASE]
  
  wt <- safe_wilcox(case, ctrl)   # case-first for consistency
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
  file.path(tsv_dir, paste0("alpha_diversity_", DATASET, ".tsv"))
)

p_alpha <- alpha |>
  pivot_longer(c(Observed, Shannon, Simpson), names_to = "metric", values_to = "value") |>
  mutate(disease = factor(disease, levels = c(CTRL, CASE))) |>
  ggplot(aes(x = disease, y = value)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(aes(color = disease), width = 0.15, size = 2, alpha = 0.7) +
  facet_wrap(~ metric, scales = "free_y") +
  scale_color_viridis_d(option = "plasma", end = 0.8) +
  labs(x = NULL, y = "Value", title = paste("Alpha diversity —", DATASET)) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

ggsave(
  file.path(plot_dir, paste0("alpha_diversity_", DATASET, ".png")),
  p_alpha, width = 12, height = 5, dpi = 300
)

############################################################################
# 6. Beta diversity
############################################################################

cat("\n=== Beta diversity ===\n")

for (metric in c("bray", "jaccard")) {
  binary <- metric == "jaccard"
  dist_obj <- phyloseq::distance(ps_spp_rel, method = metric, binary = binary)
  ord <- ordinate(ps_spp_rel, method = "PCoA", distance = dist_obj)
  
  md <- data.frame(sample_data(ps_spp_rel))[labels(dist_obj), , drop = FALSE]
  md$disease <- factor(md$disease, levels = c(CTRL, CASE))
  
  use_strata <- !is.null(LOC_COL) &&
    LOC_COL %in% colnames(md) &&
    length(unique(na.omit(md[[LOC_COL]]))) > 1
  
  if (use_strata) {
    perm <- adonis2(dist_obj ~ location + disease,
                    data = md |> mutate(location = factor(.data[[LOC_COL]])),
                    permutations = 999, by = "margin")
  } else {
    perm <- adonis2(dist_obj ~ disease, data = md, permutations = 999, by = "margin")
  }
  
  cat("\nPERMANOVA (", metric, "):\n", sep = ""); print(perm)
  
  pw <- NULL
  if (length(unique(md$disease)) > 2) {
    pw <- pairwiseAdonis::pairwise.adonis(
      as.matrix(dist_obj), factors = md$disease, perm = 999, p.adjust.m = "BH"
    )
  }
  
  bd <- betadisper(dist_obj, md$disease)
  bd_anova <- anova(bd)
  bd_tukey <- tryCatch(TukeyHSD(bd), error = function(e) NULL)
  
  write_tsv(
    as.data.frame(perm) |> rownames_to_column("term"),
    file.path(tsv_dir, paste0("permanova_", metric, "_", DATASET, ".tsv"))
  )
  if (!is.null(pw)) {
    write_tsv(pw, file.path(tsv_dir, paste0("pairwise_permanova_", metric, "_", DATASET, ".tsv")))
  }
  write_tsv(
    broom::tidy(bd_anova),
    file.path(tsv_dir, paste0("betadisper_anova_", metric, "_", DATASET, ".tsv"))
  )
  if (!is.null(bd_tukey)) {
    capture.output(
      bd_tukey,
      file = file.path(tsv_dir, paste0("betadisper_tukey_", metric, "_", DATASET, ".txt"))
    )
  }
  
  p_pcoa <- plot_ordination(ps_spp_rel, ord, color = "disease") +
    geom_point(size = 3, alpha = 0.8) +
    scale_color_viridis_d(option = "plasma", end = 0.8) +
    labs(title = paste0("PCoA (", metric, ") — ", DATASET)) +
    theme_minimal(base_size = 14)
  
  if (min(table(md$disease)) >= 3) {
    p_pcoa <- p_pcoa + stat_ellipse(linetype = 2, linewidth = 1)
  }
  
  ggsave(
    file.path(plot_dir, paste0("pcoa_", metric, "_", DATASET, ".png")),
    p_pcoa, width = 8, height = 6, dpi = 300
  )
  
  # Save ordination payload for the cross-cohort beta panel in Script 3.
  # We keep the eigenvalue vector (for axis % of variance), the sample
  # coordinates on the first two axes, and the disease label per sample.
  eig     <- ord$values$Eigenvalues
  rel_eig <- eig / sum(eig[eig > 0]) * 100
  coords  <- as.data.frame(ord$vectors[, 1:2, drop = FALSE])
  colnames(coords) <- c("Axis1", "Axis2")
  coords$sample <- rownames(coords)
  coords$disease <- md[coords$sample, "disease"]
  
  saveRDS(
    list(
      dataset      = DATASET,
      metric       = metric,
      coords       = coords,
      rel_eig      = rel_eig,
      permanova    = as.data.frame(perm) |> rownames_to_column("term"),
      betadisper_p = bd_anova[["Pr(>F)"]][1],
      betadisper_F = bd_anova[["F value"]][1]
    ),
    file.path(tsv_dir, paste0("beta_panel_", metric, "_", DATASET, ".rds"))
  )
  
  png(
    file.path(plot_dir, paste0("betadisper_", metric, "_", DATASET, ".png")),
    width = 8, height = 6, units = "in", res = 300
  )
  boxplot(
    bd, col = viridis(2, end = 0.8),
    ylab = "Distance to centroid",
    main = paste0("Betadisper (", metric, ") — ", DATASET)
  )
  dev.off()
}

############################################################################
# 7. MaAsLin2 bar plot
############################################################################

if (length(maaslin_sig) > 0) {
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
      title = paste("MaAsLin2 —", DATASET, "(", COMPARISON, "vs", CTRL, ")"),
      x = NULL, y = "Case-oriented effect size (coef)"
    ) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "none")
  
  ggsave(
    file.path(plot_dir, paste0("maaslin2_barplot_", DATASET, ".png")),
    p_bar, width = 10, height = 6, dpi = 300
  )
}

############################################################################
# 8. Summary
############################################################################

cat("\n=== Pipeline complete:", DATASET, "===\n")
cat("MaAsLin2 significant (q<", Q_THRESH, "): ", length(maaslin_sig), "\n", sep = "")
cat("Wilcoxon significant (q<", Q_THRESH, "): ", length(wilcox_sig), "\n", sep = "")
cat("ALDEx2 significant  (q<", Q_THRESH, "): ", length(aldex_sig), "\n", sep = "")
cat("MaAsLin2-Wilcoxon overlap:", length(intersect(maaslin_sig, wilcox_sig)), "\n")
cat("MaAsLin2-ALDEx2 overlap:", length(intersect(maaslin_sig, aldex_sig)), "\n")
cat(
  "High-confidence taxa:",
  sum(full_comparison$confidence == "high_confidence", na.rm = TRUE), "\n"
)
cat("Plots:", plot_dir, "\n")
cat("Tables:", tsv_dir, "\n")
