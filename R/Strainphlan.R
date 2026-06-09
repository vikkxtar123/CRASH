## =============================================================================
## Strainphlan.R
## CRASH — Cancer-associated gut micRobiome Analysis of Strains in
##         Haematological patients
##
## Purpose:
##   Strain-level phylogenetic signal analysis using StrainPhlAn4 consensus
##   marker trees and Fritz & Purvis D-statistic (caper::phylo.d). Tests
##   whether disease status (case vs control) shows non-random clustering on
##   per-species phylogenetic trees, as a proxy for disease-associated strain
##   enrichment.
##
## MODE-based cohort separation (critical design decision):
##   "AML"      → PRJNA813705 (Belgium) + Kulecka leukemia (Poland)
##   "Lymphoma" → Kulecka lymphoma (Poland) + CRA007433 (China/NKTCL)
##   "all"      → all four cohorts (QC/exploratory only; do not report)
##   Without separation, samples cluster by geography (Truong et al. 2017),
##   confounding disease signal with country of origin.
##
## Run order:
##   1. MODE = "AML"      — primary analysis, report these results
##   2. MODE = "Lymphoma" — secondary analysis
##   3. MODE = "all"      — QC only
##
## Key analytical decisions vs reference implementation:
##   - Dataset confounding D-test: flags species where clustering tracks
##     cohort origin rather than disease (upper-triangle pairwise distances)
##   - PGLS removed: invalid on binary response variable
##   - Upper-triangle-only distance matrices: avoids double-counting pairs
##   - Per-tree seeding (seed 42 + tree index) for reproducibility
##   - Minimum sample filters (MIN_MATCHED_SAMPLES, MIN_GROUP_N) with
##     logged skip reasons per tree
##   - Single CSV summary output: all species, D values, p-values, sample
##     counts, confounding flags — ready for thesis reporting
##
## Results (after BH correction):
##   AML mode:      min q = 0.296 — null result
##   Lymphoma mode: min q = 0.527 — null result
##
## Input:
##   StrainPhlAn4 Newick trees (.nwk) per species, one per MODE directory
##   Per-cohort metadata Excel files
##
## Outputs (per MODE, written to strainphlan_results/):
##   d_test_results_<MODE>.csv     All species; D, p-values, BH q, flags
##   tree_plots/                   Per-species phylogenetic tree PNGs
##
## Dependencies:
##   ape, ggtree, patchwork, phytools, tidyverse, vegan, reshape2,
##   caper, readxl, janitor
##
## Notes:
##   - select <- dplyr::select set to prevent namespace masking
##   - options(bitmapType = "cairo") required on LUNARC compute nodes
##   - StrainPhlAn4 run with --phylophlan_mode fast
## =============================================================================

library(ape)
library(ggtree)
library(patchwork)
library(phytools)
library(tidyr)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(vegan)
library(reshape2)
library(caper)
library(readxl)
library(janitor)

# Prevent MASS/ALDEx2 masking dplyr::select
select <- dplyr::select

# Reproducibility and Cairo graphics (required on compute nodes)
options(bitmapType = "cairo")
set.seed(42)

# ============================================================
# [CHANGE 1] CONFIGURATION — set MODE before running
# ============================================================

MODE <- "Lymphoma"   # CHANGE THIS: "AML", "Lymphoma", or "all"

MIN_TREE_TIPS       <- 10   # trees with fewer tips are skipped
MIN_MATCHED_SAMPLES <- 5    # trees where <5 samples match metadata are skipped
MIN_GROUP_N         <- 5    # trees where either group has <5 samples are skipped

# ============================================================
# PATHS — update if directory structure changes
# ============================================================

treepath   <- "D:/MasterThesis/Vik/STRAINPHLAN_ALL/trees_only"
resultpath <- paste0("D:/MasterThesis/Vik/STRAINPHLAN_ALL/Vikresults_", MODE, "/")

SGB2GTDB_PATH <- "D:/MasterThesis/Vik/STRAINPHLAN_ALL/mpa_vJun23_CHOCOPhlAnSGB_202403_SGB2GTDB.tsv"

# Metadata Excel files — one per cohort
META_PATHS <- list(
  prjna813705 = "D:/MasterThesis/Vik/PRJNA813705/PRJNA813705_metadata.xlsx",
  leukemia    = "D:/MasterThesis/Vik/Kulecka/Leukemia/Kuleckaleukemia_metadata.xlsx",
  lymphoma    = "D:/MasterThesis/Vik/Kulecka/Lymphoma/Kuleckalymphoma_metadata.xlsx",
  cra007433   = "D:/MasterThesis/Vik/CRA007433/CRA007433_metadata.xlsx"
)

dir.create(resultpath,                         recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(resultpath, "Small"),     recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(resultpath, "Big"),       recursive = TRUE, showWarnings = FALSE)

# ============================================================
# HELPER FUNCTIONS
# ============================================================

move_to_small <- function(filename) {
  src <- file.path(resultpath, filename)
  dst <- file.path(resultpath, "Small", filename)
  if (file.exists(src)) file.rename(src, dst)
}

move_to_big <- function(filename) {
  src <- file.path(resultpath, filename)
  dst <- file.path(resultpath, "Big", filename)
  if (file.exists(src)) file.rename(src, dst)
}

# Safe extractors for D-test result lists — handle NULLs gracefully
safe_num <- function(lst, field) {
  vapply(lst, function(x) {
    if (is.null(x)) return(NA_real_)
    val <- x[[field]]
    if (is.null(val) || length(val) == 0) return(NA_real_)
    as.numeric(val[1])
  }, FUN.VALUE = numeric(1))
}

safe_chr <- function(lst) {
  vapply(lst, function(x) {
    if (is.null(x) || length(x) == 0) return(NA_character_)
    as.character(x[1])
  }, FUN.VALUE = character(1))
}

make_safe_filename <- function(x) {
  x <- gsub("[/\\\\:*?\"<>|]", "_", x)
  x <- gsub("\\s+", "_", x)
  x <- gsub("__+", "_", x)
  x
}

# ============================================================
# COMBINED METADATA
# [CHANGE 1] All four cohorts loaded and harmonised here.
# group1 = "cancer" or "control" — the D-test binary variable.
# dataset = cohort identifier — used for the confounding D-test.
# ============================================================

meta_raw <- bind_rows(
  read_excel(META_PATHS$prjna813705) |>
    clean_names() |>
    mutate(dataset = "prjna813705", cancer_type = "AML") |>
    select(sample, disease, dataset, cancer_type),
  
  read_excel(META_PATHS$leukemia) |>
    clean_names() |>
    mutate(dataset = "leukemia", cancer_type = "AML") |>
    select(sample, disease, dataset, cancer_type),
  
  read_excel(META_PATHS$lymphoma) |>
    clean_names() |>
    mutate(dataset = "lymphoma", cancer_type = "Lymphoma") |>
    select(sample, disease, dataset, cancer_type),
  
  read_excel(META_PATHS$cra007433) |>
    clean_names() |>
    mutate(dataset = "cra007433", cancer_type = "Lymphoma") |>
    select(sample, disease, dataset, cancer_type)
) |>
  mutate(
    sample_id     = gsub("\\s+", "_", trimws(as.character(sample))),
    disease_clean = tolower(trimws(as.character(disease))),
    # Map all control-related labels to "control", everything else to "cancer"
    group1 = case_when(
      disease_clean %in% c(
        "control", "healthy", "hc", "ctl", "normal",
        "non cancer", "non-cancer", "healthy control",
        "non-cancer control", "non cancer control"
      ) ~ "control",
      TRUE ~ "cancer"
    )
  )

# Print mapping table so you can verify labels are correct before running
cat("=== Disease label mapping (verify this before proceeding) ===\n")
print(
  meta_raw |>
    distinct(dataset, disease, disease_clean, group1) |>
    arrange(dataset, disease_clean)
)

# [CHANGE 1] Filter to cohorts relevant to the current MODE
metadata_pre_dedup <- switch(
  MODE,
  all      = meta_raw,
  AML      = meta_raw |> filter(dataset %in% c("prjna813705", "leukemia")),
  Lymphoma = meta_raw |> filter(dataset %in% c("lymphoma", "cra007433")),
  stop("Unknown MODE: '", MODE, "'. Use 'all', 'AML', or 'Lymphoma'.")
)

# Deduplicate — Kulecka datasets share some control samples
metadata_df <- metadata_pre_dedup |>
  distinct(sample_id, .keep_all = TRUE) |>
  select(Sample = sample_id, group1, dataset, cancer_type, disease)

cat("\n=== Metadata after MODE filter and dedup ===\n")
cat("MODE:", MODE, "| Rows:", nrow(metadata_df), "\n")
print(table(metadata_df$group1, metadata_df$dataset, useNA = "ifany"))

# ============================================================
# SGB → SPECIES NAME LOOKUP
# [NO CHANGE] Same logic as Ramsha's, just using read.table for .tsv
# ============================================================

SGB <- read.table(
  SGB2GTDB_PATH,
  sep = "\t", header = FALSE, stringsAsFactors = FALSE
) |>
  setNames(c("SGB", "Taxonomy")) |>
  separate(
    Taxonomy,
    into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
    sep = ";", fill = "right", extra = "merge"
  ) |>
  mutate(across(everything(), ~ gsub("^.*__", "", .)))

SGB$SGB <- gsub("_group", "", SGB$SGB)

# ============================================================
# TREE FILE LIST
# ============================================================

fnFs_full <- sort(list.files(
  treepath,
  pattern   = "^RAxML_bestTree.*\\.StrainPhlAn4\\.tre$",
  full.names = TRUE,
  recursive  = TRUE
))
fnFs <- basename(fnFs_full)
cat("\nTrees found:", length(fnFs), "\n")

extracted_SGB <- str_extract(fnFs, "SGB\\d+")
species <- sapply(extracted_SGB, function(sgb) {
  m <- SGB |> filter(SGB == sgb) |> pull(Species)
  if (length(m) > 0 && !is.na(m[1]) && nzchar(m[1])) m[1] else NA_character_
})

cfnFs <- data.frame(
  filename = fnFs,
  fullpath = fnFs_full,
  SGB      = extracted_SGB,
  species  = species,
  stringsAsFactors = FALSE
)

# ============================================================
# LOOP STORAGE VECTORS
# ============================================================

DI          <- vector("list", length(fnFs))   # disease D-test results
DI_dataset  <- vector("list", length(fnFs))   # dataset confounding D-test results
SP          <- vector("list", length(fnFs))   # species labels
COMP        <- vector("list", length(fnFs))   # sample composition strings
N_TIPS      <- rep(NA_integer_, length(fnFs))
N_MATCHED   <- rep(NA_integer_, length(fnFs))
N_CANCER    <- rep(NA_integer_, length(fnFs))
N_CONTROL   <- rep(NA_integer_, length(fnFs))
SKIP_REASON <- rep(NA_character_, length(fnFs))

# ============================================================
# TREE PLOTTING FUNCTION
# Shows disease status (fill ring) and dataset origin (colour dots)
# as two independent aesthetics — avoids new_scale_fill() issues
# in newer ggtree/ggplot2 versions.
# ============================================================

plot_tree <- function(tree, cc, MR, title, branch.length = TRUE) {
  
  bl <- if (branch.length) "branch.length" else "none"
  p  <- ggtree::ggtree(tree, layout = "circular", branch.length = bl)
  
  if (branch.length) p <- p + geom_treescale()
  if (length(tree$tip.label) <= 30) p <- p + geom_tiplab(size = 2)
  
  if (!is.null(MR) && length(MR) > 0) {
    p <- p + geom_hilight(
      data = data.frame(node = MR),
      aes(node = node),
      type = "roundrect", fill = "grey80", alpha = 0.25,
      colour = NA, show.legend = FALSE
    )
  }
  
  tip_data <- p$data |>
    dplyr::filter(isTip) |>
    dplyr::select(label, x, y) |>
    dplyr::left_join(
      cc |> tibble::rownames_to_column("label") |>
        dplyr::select(label, group1, dataset),
      by = "label"
    )
  
  x_range    <- diff(range(p$data$x, na.rm = TRUE))
  if (is.na(x_range) || x_range == 0) x_range <- 1
  max_x      <- max(p$data$x, na.rm = TRUE)
  disease_x  <- max_x + 0.03 * x_range
  dataset_x  <- max_x + 0.08 * x_range
  tile_width <- 0.04 * x_range
  point_size <- if (length(tree$tip.label) <= 60) 3.5 else 2.5
  
  if (branch.length) {
    p <- p + geom_segment(
      data = tip_data,
      aes(x = x, xend = disease_x - 0.5 * tile_width, y = y, yend = y),
      color = "grey85", linewidth = 0.2, inherit.aes = FALSE
    )
  }
  
  p +
    geom_tile(
      data = tip_data,
      aes(x = disease_x, y = y, fill = group1),
      width = tile_width, height = 1, inherit.aes = FALSE, na.rm = TRUE
    ) +
    scale_fill_manual(
      values  = c(cancer = "#D55E00", control = "#F0E442"),
      breaks  = c("cancer", "control"),
      drop    = FALSE, name = "Disease", na.value = "grey90"
    ) +
    geom_point(
      data = tip_data,
      aes(x = dataset_x, y = y, color = dataset),
      size = point_size, shape = 15, inherit.aes = FALSE, na.rm = TRUE
    ) +
    scale_color_manual(
      values = c(
        cra007433   = "#66C2A5",
        lymphoma    = "#8DA0CB",
        leukemia    = "#FC8D62",
        prjna813705 = "#E78AC3"
      ),
      name = "Dataset", na.value = "grey70"
    ) +
    ggtitle(title) +
    theme(plot.margin = margin(5.5, 30, 5.5, 5.5))
}

# ============================================================
# MAIN LOOP
# ============================================================

for (j in seq_along(fnFs)) {
  
  sp_label <- if (!is.na(cfnFs$species[j]) && nzchar(cfnFs$species[j])) {
    cfnFs$species[j]
  } else {
    cfnFs$filename[j]
  }
  
  SP[[j]] <- sp_label
  
  nana <- paste0(
    MODE, "_",
    make_safe_filename(sp_label), "_",
    gsub("\\.StrainPhlAn4\\.tre$", "", cfnFs$filename[j]),
    "_tree.pdf"
  )
  
  cat("\n---", j, "/", length(fnFs), ":", sp_label, "---\n")
  
  pdf(file = file.path(resultpath, nana), width = 20, height = 20)
  
  # ----------------------------------------------------------
  # Load and root tree
  # [NO CHANGE] Midpoint rooting — standard for RAxML output
  # without an outgroup (Truong et al. 2017)
  # ----------------------------------------------------------
  
  tree        <- read.tree(cfnFs$fullpath[j])
  rooted_tree <- midpoint.root(tree)
  N_TIPS[j]  <- length(tree$tip.label)
  
  p1 <- ggtree(tree)        + ggtitle("Unrooted")         + geom_tiplab(size = 2.5)
  p2 <- ggtree(rooted_tree) + ggtitle("Midpoint rooted")  + geom_tiplab(size = 2.5)
  print(p1 | p2)
  
  # ----------------------------------------------------------
  # [CHANGE 7] Skip filters with logged reasons
  # ----------------------------------------------------------
  
  if (N_TIPS[j] < MIN_TREE_TIPS) {
    cat("  <", MIN_TREE_TIPS, "strains — skipping\n")
    SKIP_REASON[j] <- "too_few_tree_tips"
    dev.off(); move_to_small(nana); next
  }
  
  matched_metadata <- metadata_df |>
    filter(Sample %in% rooted_tree$tip.label)
  
  N_MATCHED[j] <- nrow(matched_metadata)
  
  if (N_MATCHED[j] < MIN_MATCHED_SAMPLES) {
    cat("  <", MIN_MATCHED_SAMPLES, "samples matched metadata — skipping\n")
    SKIP_REASON[j] <- "too_few_matched_samples"
    dev.off(); move_to_small(nana); next
  }
  
  if (length(unique(matched_metadata$group1)) < 2) {
    cat("  Only one disease group present — skipping\n")
    SKIP_REASON[j] <- "only_one_disease_group"
    dev.off(); move_to_small(nana); next
  }
  
  group_counts <- table(matched_metadata$group1)
  if (any(group_counts < MIN_GROUP_N)) {
    cat("  <", MIN_GROUP_N, "samples in one group — skipping\n")
    print(group_counts)
    SKIP_REASON[j] <- "too_few_samples_per_group"
    dev.off(); move_to_small(nana); next
  }
  
  subset_tree <- keep.tip(rooted_tree, matched_metadata$Sample)
  N_CANCER[j]  <- sum(matched_metadata$group1 == "cancer")
  N_CONTROL[j] <- sum(matched_metadata$group1 == "control")
  
  cat("  Samples:", nrow(matched_metadata),
      "| Cancer:", N_CANCER[j], "| Control:", N_CONTROL[j], "\n")
  cat("  Dataset breakdown:\n")
  print(table(matched_metadata$group1, matched_metadata$dataset))
  
  composition <- matched_metadata |>
    count(dataset, group1) |>
    mutate(txt = paste(dataset, group1, n, sep = ":")) |>
    pull(txt) |> paste(collapse = "; ")
  COMP[[j]] <- composition
  
  # ----------------------------------------------------------
  # [CHANGE 4] Pairwise distances — upper triangle only
  # Ramsha's script used the full symmetric matrix (each pair
  # counted twice + self-distances = 0). Upper triangle only
  # gives one unbiased distance per pair.
  # ----------------------------------------------------------
  
  dista_matrix <- cophenetic.phylo(subset_tree)
  
  gg2unif <- melt(
    dista_matrix,
    varnames   = c("Var1", "Var2"),
    value.name = "value"
  ) |>
    mutate(Var1 = as.character(Var1), Var2 = as.character(Var2)) |>
    filter(Var1 < Var2, value > 0) |>
    drop_na()
  
  if (nrow(gg2unif) == 0) {
    cat("  No valid pairwise distances — skipping\n")
    SKIP_REASON[j] <- "no_valid_pairwise_distances"
    dev.off(); move_to_small(nana); next
  }
  
  # Quantile thresholds
  quantile_5pc <- quantile(gg2unif$value, 0.05, na.rm = TRUE)
  quantile_3pc <- quantile(gg2unif$value, 0.03, na.rm = TRUE)
  quantile_1pc <- quantile(gg2unif$value, 0.01, na.rm = TRUE)
  
  # KDE thresholds
  density_est <- density(as.numeric(gg2unif$value), adjust = 1)
  deriv  <- diff(density_est$y)
  deriv2 <- diff(deriv)
  
  derivX  <- data.frame(V1 = density_est$x[-1],      deriv  = deriv)
  deriv2X <- data.frame(V1 = density_est$x[-c(1,2)], deriv2 = deriv2)
  
  min_pts  <- density_est$x[which(diff(sign(deriv))  > 0)]
  min_pts2 <- density_est$x[which(diff(sign(deriv2)) > 0)]
  
  threshold_kernel  <- if (length(min_pts)  > 0) min_pts[which.min(min_pts)]   else NA_real_
  threshold_kernel2 <- if (length(min_pts2) > 0) min_pts2[which.min(min_pts2)] else NA_real_
  
  densplot <- ggplot(gg2unif, aes(x = value)) +
    geom_density(fill = "steelblue", alpha = 0.4) +
    geom_vline(aes(xintercept = quantile_5pc,      color = "5% quantile"),   linetype = "dotted", lwd = 1) +
    geom_vline(aes(xintercept = quantile_3pc,      color = "3% quantile"),   linetype = "dotted", lwd = 1) +
    geom_vline(aes(xintercept = quantile_1pc,      color = "1% quantile"),   linetype = "dotted", lwd = 1) +
    geom_vline(aes(xintercept = threshold_kernel,  color = "KDE 1st deriv"), linetype = "dotted", lwd = 1) +
    geom_vline(aes(xintercept = threshold_kernel2, color = "KDE 2nd deriv"), linetype = "dotted", lwd = 1) +
    geom_line(data = derivX,  aes(x = V1, y = deriv  * 10), color = "lightgray", lwd = 1) +
    geom_line(data = deriv2X, aes(x = V1, y = deriv2 * 20), color = "gray40",    lwd = 1) +
    scale_color_manual(values = c(
      "5% quantile"  = "red", "3% quantile" = "orange", "1% quantile" = "blue",
      "KDE 1st deriv" = "purple", "KDE 2nd deriv" = "green"
    )) +
    theme_minimal() +
    ggtitle(sp_label) +
    xlab("Pairwise phylogenetic distance") + ylab("Density") +
    theme(legend.title = element_blank(), legend.position = "bottom")
  
  print(densplot)
  
  # ----------------------------------------------------------
  # Tight clade detection (exploratory visualisation only)
  # Based on KDE 2nd derivative threshold. Not used for
  # statistical inference — only for highlighting on tree plots.
  # Note: Do not call these "same strains" in the thesis — call
  # them "closely related strains below the KDE threshold".
  # Truong et al. 2017 calibrate same-strain using intra-
  # individual SNV rates; we lack longitudinal data for this.
  # ----------------------------------------------------------
  
  MR <- NULL
  
  if (!is.na(threshold_kernel2)) {
    gg2unif_k2 <- gg2unif |> filter(value < threshold_kernel2)
    
    if (nrow(gg2unif_k2) > 0) {
      for (aa in unique(gg2unif_k2$Var1)) {
        partners <- as.character(gg2unif_k2$Var2[gg2unif_k2$Var1 == aa])
        branch   <- unique(c(aa, partners))
        if (length(branch) >= 2) {
          mrca <- tryCatch(findMRCA(subset_tree, tips = branch), error = function(e) NULL)
          if (!is.null(mrca)) MR <- c(MR, mrca)
        }
      }
      MR <- unique(MR)
    }
  }
  
  # ----------------------------------------------------------
  # Tree plots
  # ----------------------------------------------------------
  
  cc <- matched_metadata |>
    select(Sample, group1, dataset) |>
    distinct() |>
    tibble::column_to_rownames("Sample")
  cc <- cc[subset_tree$tip.label, , drop = FALSE]
  cc$group1  <- factor(cc$group1,  levels = c("cancer", "control"))
  cc$dataset <- factor(cc$dataset)
  
  print(plot_tree(subset_tree, cc, MR,
                  paste0(sp_label, " (no branch length)"),
                  branch.length = FALSE))
  
  print(plot_tree(subset_tree, cc, MR,
                  sp_label,
                  branch.length = TRUE))
  
  # ----------------------------------------------------------
  # [NO CHANGE] Disease D-test
  # phylo.d tests whether cancer/control status shows non-random
  # phylogenetic clustering (Pval1) and whether that clustering
  # is compatible with Brownian motion (Pval0).
  # Significant result = Pval1 < 0.05 (non-random).
  # D < 0 = stronger than Brownian (highly conserved trait).
  # D = 0 = Brownian. D = 1 = random. D > 1 = overdispersed.
  # ----------------------------------------------------------
  
  D_result <- NULL
  
  trait <- data.frame(Sample = subset_tree$tip.label) |>
    left_join(metadata_df, by = "Sample") |>
    mutate(binarygroup = if_else(group1 == "cancer", 1, 0))
  rownames(trait) <- trait$Sample
  
  if (length(unique(trait$binarygroup)) > 1) {
    tryCatch({
      set.seed(42 + j)
      D_result <- phylo.d(
        trait, subset_tree,
        names.col = Sample, binvar = binarygroup
      )
      print(D_result)
    }, error = function(e) {
      cat("  Disease D-test failed:", conditionMessage(e), "\n")
    })
  } else {
    cat("  Disease D-test skipped: only one group\n")
  }
  
  DI[j] <- list(D_result)
  
  # ----------------------------------------------------------
  # [CHANGE 2] Dataset confounding D-test
  # Only runs when exactly two datasets are present in the tree.
  # Tests whether strain clustering tracks dataset origin rather
  # than disease. If significant alongside the disease D-test,
  # flag the species as potentially confounded.
  # ----------------------------------------------------------
  
  cat("\n=== Length check before summary ===\n")
  print(c(
    fnFs       = length(fnFs),
    SP         = length(SP),
    DI         = length(DI),
    DI_dataset = length(DI_dataset),
    COMP       = length(COMP),
    N_TIPS     = length(N_TIPS),
    SKIP       = length(SKIP_REASON)
  ))
  
  D_dataset_result <- NULL
  
  dataset_trait <- data.frame(Sample = subset_tree$tip.label) |>
    left_join(metadata_df, by = "Sample")
  
  if (length(unique(dataset_trait$dataset)) == 2) {
    dataset_trait <- dataset_trait |>
      mutate(binarydataset = as.numeric(factor(dataset)) - 1)
    rownames(dataset_trait) <- dataset_trait$Sample
    
    tryCatch({
      set.seed(4200 + j)
      D_dataset_result <- phylo.d(
        dataset_trait, subset_tree,
        names.col = Sample, binvar = binarydataset
      )
      cat("  Dataset confounding D-test:\n")
      print(D_dataset_result)
    }, error = function(e) {
      cat("  Dataset D-test failed:", conditionMessage(e), "\n")
    })
  } else {
    cat("  Dataset D-test skipped: not exactly two datasets in this tree\n")
  }
  
  DI_dataset[j] <- list(D_dataset_result)
  SKIP_REASON[j]  <- "processed"
  
  dev.off()
  move_to_big(nana)
}

# ============================================================
# SAVE RDS OBJECTS
# ============================================================

saveRDS(DI,         file.path(resultpath, paste0("DI_disease_",  MODE, ".rds")))
saveRDS(DI_dataset, file.path(resultpath, paste0("DI_dataset_",  MODE, ".rds")))
saveRDS(SP,         file.path(resultpath, paste0("SP_",          MODE, ".rds")))
saveRDS(COMP,       file.path(resultpath, paste0("composition_", MODE, ".rds")))

# ============================================================
# [CHANGE 6] CONSOLIDATED SUMMARY CSV
# Ramsha's script had no summary table — results were buried in
# per-tree RDS files. This table is what you'll filter and
# report in your thesis results section.
# ============================================================

d_summary <- data.frame(
  Species           = safe_chr(SP),
  D                 = safe_num(DI, "DEstimate"),
  pval_rand         = safe_num(DI, "Pval1"),
  pval_bm           = safe_num(DI, "Pval0"),
  D_dataset         = safe_num(DI_dataset, "DEstimate"),
  pval_dataset_rand = safe_num(DI_dataset, "Pval1"),
  pval_dataset_bm   = safe_num(DI_dataset, "Pval0"),
  n_tree_tips       = N_TIPS,
  n_matched         = N_MATCHED,
  n_cancer          = N_CANCER,
  n_control         = N_CONTROL,
  composition       = safe_chr(COMP),
  skip_reason       = SKIP_REASON,
  mode              = MODE,
  stringsAsFactors  = FALSE
)

# FDR across all trees actually tested in this MODE (skipped trees stay NA).
# Subset before p.adjust so the NA trees don't inflate n.
ok    <- !is.na(d_summary$pval_rand)
ok_ds <- !is.na(d_summary$pval_dataset_rand)
d_summary$qval_rand         <- NA_real_
d_summary$qval_dataset_rand <- NA_real_
d_summary$qval_rand[ok]            <- p.adjust(d_summary$pval_rand[ok], "BH")
d_summary$qval_dataset_rand[ok_ds] <- p.adjust(d_summary$pval_dataset_rand[ok_ds], "BH")

d_summary <- d_summary |>
  mutate(
    # PRIMARY criterion: clustering vs random, FDR-corrected
    significant        = !is.na(qval_rand) & qval_rand < 0.05,
    # DESCRIPTOR only — no longer gates significance
    compatible_with_bm = !is.na(pval_bm) & pval_bm > 0.05,
    dataset_signal     = !is.na(qval_dataset_rand) & qval_dataset_rand < 0.05,
    
    interpretation = case_when(
      is.na(D)           ~ "not tested",
      !significant       ~ "indistinguishable from random (q >= 0.05)",
      D < 0              ~ "stronger than Brownian — highly conserved trait",
      compatible_with_bm ~ "Brownian-compatible clustering",
      D < 1              ~ "phylogenetic clustering, weaker than Brownian",
      D > 1              ~ "overdispersed — cancer/control interleaved",
      TRUE               ~ "clustered"
    ),
    
    confounding_flag = case_when(
      significant  &  dataset_signal ~ "disease_and_dataset_signal — interpret with caution",
      significant  & !dataset_signal ~ "disease_signal_only — reliable",
      !significant &  dataset_signal ~ "dataset_signal_only — geography confound",
      TRUE                           ~ "no_signal"
    )
  )

write.csv(d_summary,
          file.path(resultpath, paste0("D_test_summary_", MODE, ".csv")),
          row.names = FALSE)

# Significant, named species only — this is your thesis results table
sig_named <- d_summary |>
  filter(significant == TRUE, !is.na(Species), nzchar(Species),
         !grepl("^RAxML_bestTree", Species)) |>
  arrange(qval_rand)

write.csv(sig_named,
          file.path(resultpath, paste0("D_test_significant_", MODE, ".csv")),
          row.names = FALSE)

cat("\n=== COMPLETE ===\n")
cat("Full summary:", nrow(d_summary), "species\n")
cat("Strict significant named species:", nrow(sig_named), "\n\n")
print(sig_named[, c("Species", "D", "qval_rand", "pval_rand", "pval_bm",
                    "n_cancer", "n_control", "confounding_flag")])

cat("\nSkip reason summary:\n")
print(table(d_summary$skip_reason, useNA = "ifany"))
