## =============================================================================
## DeepARG_3.R
## CRASH — Cancer-associated gut micRobiome Analysis of Strains in
##         Haematological patients
##
## Purpose:
##   Cross-cohort publication figures for deepARG AMR results. Set LEVEL to
##   match DeepARG_1.R and DeepARG_2.R.
##
##     Fig 1. Alpha diversity panel (ARG richness / Shannon; case/ctrl ratio)
##     Fig 2. Cross-dataset heatmap (ARGs consistent across ≥3 cohorts)
##     Fig 3. ARG category coefficient heatmap (type-level)
##     Fig 4. ARG subtype coefficient heatmap
##     Fig 5. Top-hits bar panel (top enriched + depleted per cohort)
##     Fig 6. Volcano panel (2×2 with top hits labelled)
##
## Usage:
##   Set LEVEL ("type" or "subtype") at top; must match Scripts 1 and 2.
##
## Inputs:
##   - Per-cohort deepARG MaAsLin2 all_results.tsv
##   - Per-cohort tables/alpha_diversity_deeparg_<LEVEL>_<DATASET>.tsv
##   - cross_dataset_deeparg/cross_dataset_overlap_deeparg_<LEVEL>.tsv
##
## Outputs (written to DeepARG/figures_<LEVEL>):
##   01_alpha_panel_clean.png  through  06_volcano_panel_clean.png
##
## Dependencies:
##   tidyverse, viridis, scales, ggrepel, patchwork, readxl, janitor
##
## Notes:
##   - Alpha significance labels from location-adjusted LM q-values
##   - File existence checks in place; missing inputs skip gracefully
## =============================================================================

library(tidyverse)
library(viridis)
library(scales)
library(ggrepel)
library(patchwork)
library(readxl)
library(janitor)

LEVEL    <- "type"   # must match Script 1 and 2
Q_THRESH <- 0.25

out_dir <- file.path("D:/MasterThesis/Vik/DeepARG", paste0("figures_", LEVEL))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

theme_set(theme_minimal(base_size = 15))

############################################################################
# Common labels / helpers
############################################################################

cohort_labels <- c(
  prjna813705      = "AML\n(PRJNA813705)",
  kulecka_leukemia = "AML\n(Kulecka)",
  kulecka_lymphoma = "LN\n(Kulecka)",
  cra              = "NKTCL\n(CRA)"
)

case_map <- c(
  prjna813705      = "AML",
  kulecka_leukemia = "AML",
  kulecka_lymphoma = "LN",
  cra              = "NKTCL"
)

tables_path <- function(ds) {
  file.path("D:/MasterThesis/Vik/DeepARG", paste0(ds, "_", LEVEL), "tables")
}

maaslin_files <- setNames(
  sapply(names(cohort_labels), function(ds)
    file.path("D:/MasterThesis/Vik/DeepARG",
              paste0(ds, "_", LEVEL), "maaslin2", "all_results.tsv")),
  names(cohort_labels)
)

alpha_files <- setNames(
  sapply(names(cohort_labels), function(ds)
    file.path(tables_path(ds),
              paste0("alpha_diversity_deeparg_", ds, ".tsv"))),
  unname(cohort_labels)
)

alpha_lm_files <- setNames(
  sapply(names(cohort_labels), function(ds)
    file.path(tables_path(ds),
              paste0("alpha_diversity_lm_deeparg_", ds, ".tsv"))),
  unname(cohort_labels)
)

dir_cols <- paste0(names(cohort_labels), "_dir")
q_cols   <- paste0(names(cohort_labels), "_q")

############################################################################
# File existence checks
############################################################################

required_files <- c(
  unname(maaslin_files),
  unlist(alpha_files),
  unlist(alpha_lm_files),
  file.path("D:/MasterThesis/Vik/DeepARG",
            paste0("cross_dataset_", LEVEL),
            paste0("cross_dataset_overlap_deeparg_", LEVEL, ".tsv"))
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required files — run Script 1 for all cohorts and Script 2 first:\n",
       paste(" -", missing_files, collapse = "\n"))
}
cat("All required files present.\n")

############################################################################
# Helpers
############################################################################

sig_label_fun <- function(q) {
  case_when(
    q < 0.001 ~ "***",
    q < 0.01  ~ "**",
    q < 0.05  ~ "*",
    q < 0.25  ~ "\u2020",
    TRUE      ~ "ns"
  )
}

fmt_p <- function(p) {
  if (is.na(p)) return("= NA")
  if (p < 0.001) "< 0.001" else sprintf("= %.3f", p)
}

# Load all MaAsLin2 disease associations, case-oriented
maaslin_all <- imap_dfr(maaslin_files, function(path, ds) {
  read_tsv(path, show_col_types = FALSE) |>
    filter(metadata == "disease") |>
    mutate(
      cohort_key = ds,
      coef_case  = if_else(value == case_map[[ds]], coef, -coef),
      is_sig     = qval < Q_THRESH
    ) |>
    select(feature, cohort_key, coef_case, qval, is_sig)
})

# Cohort sample counts from metadata
meta_info <- list(
  prjna813705      = list(path = "D:/MasterThesis/Vik/PRJNA813705/PRJNA813705_metadata.xlsx",
                          case = "AML",   ctrl = "Control"),
  kulecka_leukemia = list(path = "D:/MasterThesis/Vik/Kulecka/Leukemia/Kuleckaleukemia_metadata.xlsx",
                          case = "AML",   ctrl = "Control"),
  kulecka_lymphoma = list(path = "D:/MasterThesis/Vik/Kulecka/Lymphoma/Kuleckalymphoma_metadata.xlsx",
                          case = "LN",    ctrl = "Control"),
  cra              = list(path = "D:/MasterThesis/Vik/CRA007433/CRA007433_metadata.xlsx",
                          case = "NKTCL", ctrl = "Control")
)

cohort_counts <- imap_dfr(meta_info, function(info, key) {
  m <- read_excel(info$path) |>
    clean_names() |>
    filter(disease %in% c(info$case, info$ctrl))
  tibble(cohort_key = key,
         n_case  = sum(m$disease == info$case,  na.rm = TRUE),
         n_ctrl  = sum(m$disease == info$ctrl, na.rm = TRUE),
         n_total = nrow(m))
})

cohort_labels_with_n <- cohort_counts |>
  mutate(label = sprintf("%s\nN=%d (%d/%d) [metadata]",
                         cohort_labels[cohort_key], n_total, n_case, n_ctrl)) |>
  select(cohort_key, label) |>
  deframe()

cat("Cohort sample counts (metadata N — DeepARG N may differ after zero-signal filter):\n")
print(cohort_counts)

############################################################################
# 1. Alpha diversity panel (significance from location-adjusted LM)
############################################################################

alpha_combined <- imap_dfr(alpha_files, function(path, label) {
  read_tsv(path, show_col_types = FALSE) |> mutate(cohort = label)
}) |>
  mutate(
    cohort = factor(cohort, levels = names(alpha_files)),
    metric = factor(metric, levels = c("Richness", "Shannon")),
    ratio  = if_else(median_ctrl > 0, median_case / median_ctrl, NA_real_)
  )

# Load LM q-values and join for significance labels
alpha_lm_combined <- imap_dfr(alpha_lm_files, function(path, label) {
  read_tsv(path, show_col_types = FALSE) |>
    select(metric, q) |>
    rename(q_lm = q) |>
    mutate(cohort = label)
})

alpha_combined <- alpha_combined |>
  left_join(alpha_lm_combined, by = c("cohort", "metric")) |>
  mutate(sig_label = sig_label_fun(q_lm))

p_alpha <- ggplot(alpha_combined, aes(x = cohort, y = ratio, fill = metric)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.72) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40", linewidth = 0.8) +
  geom_text(aes(label = sig_label, y = ratio + 0.02),
            position = position_dodge(width = 0.8), size = 5, vjust = 0) +
  scale_fill_viridis_d(option = "plasma", end = 0.82) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(title    = paste("ARG", LEVEL, "diversity: median case/control ratio across cohorts"),
       subtitle = "Significance from location-adjusted LM where applicable. \u2020 q<0.25   * q<0.05   ** q<0.01   *** q<0.001",
       x = NULL, y = "Median case/control ratio", fill = NULL) +
  theme(legend.position  = "top",
        panel.grid.minor = element_blank(),
        axis.text.x      = element_text(size = 14),
        plot.title       = element_text(face = "bold"))

ggsave(file.path(out_dir, "01_alpha_diversity_deeparg_panel.png"),
       p_alpha, width = 12, height = 6, dpi = 300)
cat("Figure 1 saved\n")

############################################################################
# 2. Cross-dataset overlap heatmap
############################################################################

overlap <- read_tsv(
  file.path("D:/MasterThesis/Vik/DeepARG",
            paste0("cross_dataset_", LEVEL),
            paste0("cross_dataset_overlap_deeparg_", LEVEL, ".tsv")),
  show_col_types = FALSE
)

heatmap_features <- overlap |>
  filter(n_datasets >= 2, direction_consistent == TRUE) |>
  rowwise() |>
  mutate(min_q = suppressWarnings(min(c_across(all_of(q_cols)), na.rm = TRUE))) |>
  ungroup() |>
  mutate(min_q = if_else(is.infinite(min_q), NA_real_, min_q))

heatmap_features <- bind_rows(
  heatmap_features |> filter(consensus_direction == "depleted") |>
    arrange(desc(n_datasets), min_q) |> slice_head(n = 10),
  heatmap_features |> filter(consensus_direction == "enriched") |>
    arrange(desc(n_datasets), min_q) |> slice_head(n = 10)
)

feat_order <- heatmap_features |>
  arrange(consensus_direction, desc(n_datasets), min_q) |>
  pull(feature)

if (length(feat_order) > 0) {
  heatmap_data <- heatmap_features |>
    select(feature, all_of(dir_cols), all_of(q_cols)) |>
    pivot_longer(cols = -feature,
                 names_to  = c("cohort_key", ".value"),
                 names_pattern = "(.+)_(dir|q)") |>
    mutate(
      cohort = factor(cohort_labels_with_n[cohort_key],
                      levels = unname(cohort_labels_with_n)),
      feature = factor(feature, levels = feat_order),
      direction_plot = case_when(
        dir == "depleted" ~ "Depleted in cancer",
        dir == "enriched" ~ "Enriched in cancer",
        TRUE              ~ "Not a hit"
      ),
      direction_plot = factor(direction_plot,
                              levels = c("Depleted in cancer",
                                         "Enriched in cancer", "Not a hit")),
      sig_label = case_when(
        is.na(q)  ~ "",
        q < 0.001 ~ "***", q < 0.01 ~ "**",
        q < 0.05  ~ "*",   q < 0.25 ~ "\u2020",
        TRUE ~ ""
      )
    )
  
  p_heat <- ggplot(heatmap_data,
                   aes(x = cohort, y = feature, fill = direction_plot)) +
    geom_tile(colour = "white", linewidth = 0.35) +
    geom_text(aes(label = sig_label), colour = "white",
              fontface = "bold", size = 3.2, vjust = 0.55) +
    scale_fill_manual(values = c("Depleted in cancer" = "#2166ac",
                                 "Enriched in cancer" = "#d6604d",
                                 "Not a hit"          = "grey88")) +
    labs(title    = paste0("Top 10 enriched + top 10 depleted ARG ",
                           LEVEL, "s across \u22652 cohorts"),
         subtitle = "Ranked by n cohorts then minimum q. \u2020 q<0.25  * q<0.05  ** q<0.01  *** q<0.001",
         x = NULL, y = NULL, fill = NULL) +
    theme(legend.position = "top", panel.grid = element_blank(),
          axis.text.y     = element_text(size = 9),
          axis.text.x     = element_text(size = 11, lineheight = 0.9),
          plot.title      = element_text(face = "bold"),
          plot.subtitle   = element_text(size = 11, colour = "grey30"))
  
  ggsave(file.path(out_dir, "02_cross_dataset_heatmap_deeparg.png"),
         p_heat, width = 10,
         height = max(6, length(feat_order) * 0.4 + 3), dpi = 300)
  cat("Figure 2 saved\n")
} else {
  message("No features in >=2 datasets with consistent direction — Figure 2 skipped.")
}

############################################################################
# 3/4/5. Coefficient heatmaps (level-specific)
############################################################################

build_coef_heatmap <- function(feature_list, title_text, subtitle_text,
                               out_file, fig_height = 7, fig_width = 10) {
  
  dat <- maaslin_all |>
    filter(feature %in% feature_list) |>
    complete(feature = feature_list, cohort_key = names(maaslin_files))
  
  if (nrow(dat) == 0 || all(is.na(dat$coef_case))) {
    message("No data for: ", title_text, " — skipping.")
    return(FALSE)
  }
  
  rank_df <- dat |>
    group_by(feature) |>
    summarise(
      n_sig         = sum(is_sig, na.rm = TRUE),
      min_q         = suppressWarnings(min(qval, na.rm = TRUE)),
      mean_abs_coef = mean(abs(coef_case), na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(min_q = if_else(is.infinite(min_q), NA_real_, min_q)) |>
    arrange(desc(n_sig), min_q, desc(mean_abs_coef))
  
  feat_order <- rank_df$feature
  max_abs    <- max(abs(dat$coef_case), na.rm = TRUE)
  
  if (!is.finite(max_abs) || max_abs == 0) max_abs <- 1
  
  dat <- dat |>
    mutate(
      cohort      = factor(cohort_labels[cohort_key],
                           levels = unname(cohort_labels)),
      feature     = factor(feature, levels = rev(feat_order)),
      is_sig_safe = coalesce(is_sig, FALSE),
      sig_marker  = case_when(
        !is_sig_safe | is.na(coef_case) ~ "",
        qval < 0.001 ~ "***",
        qval < 0.01  ~ "**",
        qval < 0.05  ~ "*",
        qval < 0.25  ~ "\u2020",
        TRUE ~ ""
      ),
      label     = if_else(is.na(coef_case), "\u2014",
                          paste0(sprintf("%.2f", coef_case), sig_marker)),
      label_col = if_else(is.na(coef_case) | abs(coef_case) < 0.6 * max_abs,
                          "black", "white")
    )
  
  p <- ggplot(dat, aes(x = cohort, y = feature, fill = coef_case)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = label, colour = label_col),
              size = 3.6, show.legend = FALSE) +
    scale_colour_identity() +
    scale_fill_gradient2(
      low = "#2166ac", mid = "white", high = "#d6604d",
      midpoint = 0, na.value = "grey88",
      limits = c(-max_abs, max_abs),
      name = "MaAsLin2\ncoefficient\n(case-oriented)"
    ) +
    labs(title = title_text, subtitle = subtitle_text, x = NULL, y = NULL) +
    theme(
      axis.text.y = element_text(size = 10),
      panel.grid  = element_blank(),
      legend.position     = "right",
      plot.title          = element_text(face = "bold"),
      plot.subtitle       = element_text(size = 10, colour = "grey30"),
      plot.title.position = "plot"
    )
  
  ggsave(out_file, p, width = fig_width, height = fig_height, dpi = 300)
  return(TRUE)
}
# Feature lists
broad_spectrum_cats <- c("tetracycline", "multidrug",
                         "aminoglycoside", "fluoroquinolone", "MLS",
                         "glycopeptide", "bacitracin", "rifamycin",
                         "phenicol", "diaminopyrimidine")

tet_subtypes <- c("TETW", "TETM", "TETQ", "TETX", "TETO",
                  "TETP", "TET32", "TET40", "TETA", "TETB")

mls_subtypes <- c("ERMB", "ERMC", "ERMF", "ERMG", "MEFA",
                  "MSRA", "LNUC", "CFXA2")

# Level-specific calls
if (LEVEL == "type") {
  build_coef_heatmap(
    broad_spectrum_cats,
    title_text    = "ARG categories: MaAsLin2 effect sizes across cohorts",
    subtitle_text = "\u2020 q<0.25  * q<0.05  ** q<0.01  *** q<0.001",
    out_file      = file.path(out_dir, "03_ARG_categories_coef.png"),
    fig_height    = 7
  )
  cat("Figure 3 saved\n")
}

if (LEVEL == "subtype") {
  build_coef_heatmap(
    tet_subtypes,
    title_text    = "Tetracycline resistance subtypes: MaAsLin2 effect sizes",
    subtitle_text = "\u2020 q<0.25  * q<0.05  ** q<0.01  *** q<0.001",
    out_file      = file.path(out_dir, "04_tetracycline_subtypes_coef.png"),
    fig_height    = 6
  )
  cat("Figure 4 saved\n")
  
  build_coef_heatmap(
    mls_subtypes,
    title_text    = "MLS resistance subtypes: MaAsLin2 effect sizes",
    subtitle_text = "\u2020 q<0.25  * q<0.05  ** q<0.01  *** q<0.001",
    out_file      = file.path(out_dir, "05_MLS_subtypes_coef.png"),
    fig_height    = 6
  )
  cat("Figure 5 saved\n")
}

############################################################################
# 6. Top-hits bar panel — top 5 enriched + top 5 depleted per cohort (2x2)
############################################################################

TOP_N_BARS <- 5

top_bars_df <- maaslin_all |>
  filter(!is.na(coef_case), !is.na(qval), qval < Q_THRESH) |>
  group_by(cohort_key) |>
  group_modify(~ {
    ups   <- .x |> arrange(desc(coef_case))  |> slice_head(n = TOP_N_BARS)
    downs <- .x |> arrange(coef_case)         |> slice_head(n = TOP_N_BARS)
    bind_rows(ups, downs) |> distinct(feature, .keep_all = TRUE)
  }) |>
  ungroup() |>
  mutate(
    cohort = factor(cohort_labels_with_n[cohort_key],
                    levels = unname(cohort_labels_with_n)),
    direction_plot = if_else(coef_case > 0, "Enriched in cancer",
                             "Depleted in cancer"),
    direction_plot = factor(direction_plot,
                            levels = c("Depleted in cancer",
                                       "Enriched in cancer")),
    sig_label = case_when(
      qval < 0.001 ~ "***", qval < 0.01 ~ "**",
      qval < 0.05  ~ "*",   qval < 0.25 ~ "\u2020",
      TRUE ~ ""
    )
  )

build_top_bars <- function(df, cohort_name) {
  d <- df |>
    filter(cohort == cohort_name) |>
    mutate(feature = fct_reorder(feature, coef_case))
  
  if (nrow(d) == 0) {
    return(ggplot() +
             annotate("text", x = 0, y = 0,
                      label = "No q<0.25 hits", colour = "grey50") +
             labs(title = cohort_name) + theme_void() +
             theme(plot.title = element_text(face = "bold", size = 13)))
  }
  
  x_max <- max(abs(d$coef_case), na.rm = TRUE) * 1.18
  
  ggplot(d, aes(x = coef_case, y = feature, fill = direction_plot)) +
    geom_col(width = 0.75) +
    geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3) +
    geom_text(aes(label = sig_label,
                  hjust = if_else(coef_case > 0, -0.25, 1.25)),
              size = 3.8, fontface = "bold", colour = "grey20") +
    scale_fill_manual(values = c("Depleted in cancer" = "#2166ac",
                                 "Enriched in cancer" = "#d6604d"),
                      drop = FALSE) +
    scale_x_continuous(limits = c(-x_max, x_max)) +
    labs(title = cohort_name,
         x = "MaAsLin2 coefficient (case-oriented)", y = NULL) +
    theme(plot.title       = element_text(face = "bold", size = 13),
          legend.position  = "none",
          axis.text.y      = element_text(size = 9),
          panel.grid.minor = element_blank())
}

bar_plots <- lapply(levels(top_bars_df$cohort),
                    function(ch) build_top_bars(top_bars_df, ch))

p_top_bars <- wrap_plots(bar_plots, ncol = 2) +
  plot_annotation(
    title    = sprintf(
      "Top %d enriched + top %d depleted ARG %ss per cohort (MaAsLin2, q<%.2f)",
      TOP_N_BARS, TOP_N_BARS, LEVEL, Q_THRESH),
    subtitle = "Bar = case-oriented coefficient. \u2020 q<0.25  * q<0.05  ** q<0.01  *** q<0.001",
    theme    = theme(plot.title    = element_text(face = "bold", size = 14),
                     plot.subtitle = element_text(size = 11, colour = "grey30"))
  ) &
  theme(plot.margin = margin(4, 6, 4, 6))

ggsave(file.path(out_dir, "06_top_hits_bars_deeparg.png"),
       p_top_bars, width = 14, height = 10, dpi = 300)
cat("Figure 6 saved\n")

############################################################################
# 7. Volcano panel — full MaAsLin2 landscape per cohort (2x2)
############################################################################

VOLCANO_LABEL_N <- 5

volcano_df <- maaslin_all |>
  filter(!is.na(qval)) |>
  mutate(
    qval_plot = pmax(qval, 1e-300),
    cohort    = factor(cohort_labels[cohort_key], levels = unname(cohort_labels)),
    neglog10q = -log10(qval_plot),
    direction_plot = case_when(
      qval >= Q_THRESH ~ "ns",
      coef_case > 0    ~ "Enriched in cancer",
      coef_case < 0    ~ "Depleted in cancer",
      TRUE             ~ "ns"
    ),
    direction_plot = factor(direction_plot,
                            levels = c("Depleted in cancer",
                                       "Enriched in cancer", "ns"))
  )

volcano_labels <- volcano_df |>
  filter(qval < Q_THRESH, !is.na(coef_case)) |>
  group_by(cohort) |>
  group_modify(~ {
    ups   <- .x |> filter(coef_case > 0) |>
      arrange(desc(abs(coef_case) * neglog10q)) |>
      slice_head(n = VOLCANO_LABEL_N)
    downs <- .x |> filter(coef_case < 0) |>
      arrange(desc(abs(coef_case) * neglog10q)) |>
      slice_head(n = VOLCANO_LABEL_N)
    bind_rows(ups, downs)
  }) |>
  ungroup()

x_max <- ceiling(max(abs(volcano_df$coef_case), na.rm = TRUE) * 10) / 10
y_max <- max(volcano_df$neglog10q, na.rm = TRUE) * 1.05

build_volcano <- function(df, labels_df, cohort_name) {
  d  <- df        |> filter(cohort == cohort_name)
  lb <- labels_df |> filter(cohort == cohort_name)
  
  ggplot(d, aes(x = coef_case, y = neglog10q)) +
    geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.4) +
    geom_hline(yintercept = -log10(Q_THRESH), linetype = "dashed",
               colour = "grey40", linewidth = 0.5) +
    geom_point(aes(colour = direction_plot), alpha = 0.75, size = 2.5) +
    geom_text_repel(data = lb, aes(label = feature),
                    size = 3, min.segment.length = 0,
                    segment.size = 0.3, segment.colour = "grey60",
                    max.overlaps = Inf, box.padding = 0.35, force = 2) +
    scale_colour_manual(
      values = c("Depleted in cancer" = "#2166ac",
                 "Enriched in cancer" = "#d6604d",
                 "ns"                 = "grey75"),
      drop = FALSE
    ) +
    scale_x_continuous(limits = c(-x_max, x_max)) +
    scale_y_continuous(limits = c(0, y_max),
                       expand = expansion(mult = c(0, 0.02))) +
    labs(title = cohort_name,
         x = "MaAsLin2 coefficient (case-oriented)",
         y = expression(-log[10](q))) +
    theme(plot.title       = element_text(face = "bold", size = 13),
          legend.position  = "none",
          panel.grid.minor = element_blank())
}

volcano_plots <- lapply(levels(volcano_df$cohort),
                        function(ch) build_volcano(volcano_df, volcano_labels, ch))

p_volcano <- wrap_plots(volcano_plots, ncol = 2) +
  plot_annotation(
    title    = paste0("ARG ", LEVEL,
                      " differential abundance: MaAsLin2 across cohorts"),
    subtitle = paste0("Dashed line: q = ", Q_THRESH, ". Top ",
                      VOLCANO_LABEL_N, " depleted + top ", VOLCANO_LABEL_N,
                      " enriched per cohort labelled (by |coef| x -log10 q)."),
    theme    = theme(plot.title    = element_text(face = "bold", size = 15),
                     plot.subtitle = element_text(size = 11, colour = "grey30"))
  ) &
  theme(plot.margin = margin(4, 6, 4, 6))

ggsave(file.path(out_dir, "07_volcano_panel_deeparg.png"),
       p_volcano, width = 13, height = 11, dpi = 300)
cat("Figure 7 saved\n")

############################################################################
# 8. Beta diversity panel — 4 cohorts (Bray-Curtis PCoA)
############################################################################

beta_rds_paths <- tibble(
  cohort_key = names(cohort_labels),
  path = map_chr(cohort_key, function(ck)
    file.path(tables_path(ck),
              paste0("beta_panel_deeparg_", LEVEL, "_", ck, ".rds")))
)

missing_rds <- beta_rds_paths |> filter(!file.exists(path))
if (nrow(missing_rds) > 0) {
  warning("Beta RDS files missing for: ",
          paste(missing_rds$cohort_key, collapse = ", "),
          "\nRerun Script 1 for those datasets first. Skipping Figure 8.")
} else {
  
  build_beta_panel <- function(rds_path, cohort_key) {
    payload <- readRDS(rds_path)
    coords  <- payload$coords
    lv      <- unique(as.character(coords$disease))
    case_lv <- setdiff(lv, "Control")
    coords$disease <- factor(coords$disease, levels = c("Control", case_lv))
    coords$Group   <- factor(
      if_else(coords$disease == "Control", "Control", "Case"),
      levels = c("Control", "Case")
    )
    
    perm     <- payload$permanova |> filter(term == "disease")
    subtitle <- sprintf(
      "PERMANOVA: R\u00B2 = %.3f, F = %.2f, p %s   \u00B7   Betadisper p %s",
      perm$R2[1], perm$`F`[1],
      fmt_p(perm[["Pr(>F)"]][1]),
      fmt_p(payload$betadisper_p)
    )
    
    rel_eig <- payload$rel_eig
    
    p_pcoa <- ggplot(coords, aes(x = Axis1, y = Axis2, colour = Group)) +
      geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
      geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
      geom_point(size = 2.2, alpha = 0.8) +
      scale_colour_manual(values = c(Control = "#1b2a6b", Case = "#e6a23c")) +
      guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
      labs(title    = cohort_labels_with_n[cohort_key],
           subtitle = subtitle,
           x = sprintf("PCo1 (%.1f%%)", rel_eig[1]),
           y = sprintf("PCo2 (%.1f%%)", rel_eig[2]),
           colour = NULL) +
      theme(plot.title    = element_text(face = "bold", size = 10),
            plot.subtitle = element_text(size = 8, colour = "grey30"),
            axis.title    = element_text(size = 9),
            panel.grid.minor = element_blank())
    
    if (min(table(coords$Group)) >= 3) {
      p_pcoa <- p_pcoa +
        stat_ellipse(aes(fill = Group), type = "norm", level = 0.68,
                     geom = "polygon", alpha = 0.12, linewidth = 0.6) +
        scale_fill_manual(values = c(Control = "#1b2a6b", Case = "#e6a23c"),
                          guide = "none")
    }
    p_pcoa
  }
  
  beta_plots <- pmap(beta_rds_paths, function(cohort_key, path) {
    build_beta_panel(path, cohort_key)
  })
  
  p_beta <- wrap_plots(beta_plots, ncol = 2) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title    = paste0("Beta diversity: PCoA (Bray-Curtis) across cohorts — ARG ",
                        LEVEL, " profiles"),
      subtitle = "Ellipses: 68% (1 SD) contours. PERMANOVA and betadisper stats in subtitles.",
      theme    = theme(plot.title    = element_text(face = "bold", size = 15),
                       plot.subtitle = element_text(size = 11, colour = "grey30"))
    ) &
    theme(plot.margin     = margin(3, 5, 3, 5),
          legend.position = "bottom",
          legend.title    = element_blank(),
          legend.text     = element_text(size = 11))
  
  ggsave(file.path(out_dir, "08_beta_diversity_deeparg_panel.png"),
         p_beta, width = 14, height = 10, dpi = 300)
  cat("Figure 8 saved\n")
}

library(magick)
a <- image_read(alpha_path); b <- image_read(beta_path)
w <- max(image_info(a)$width, image_info(b)$width)     # match widths so edges align
a <- image_resize(a, paste0(w, "x")); b <- image_resize(b, paste0(w, "x"))
combined <- image_append(c(a, b), stack = TRUE)
image_write(combined, file.path(out_dir, "09_alpha_beta_panel_deeparg.png"))


cat("\nAll deepARG figures saved to:", out_dir, "\n")
