## =============================================================================
## KrakenAnalysis_3.R
## CRASH — Cancer-associated gut micRobiome Analysis of Strains in
##         Haematological patients
##
## Purpose:
##   Cross-cohort publication figures for Kraken2/Bracken results. Parallel
##   to MetaAnalysis_3.R; uses the same design language so MetaPhlAn and
##   Kraken panels sit side-by-side for classifier comparison.
##
##     Fig 1. Alpha diversity panel (case/ctrl ratio per cohort)
##     Fig 2. Cross-dataset heatmap (species consistent across ≥3 cohorts)
##     Fig 3. Butyrate producer coefficient heatmap
##     Fig 4. Opportunistic taxon coefficient heatmap
##     Fig 5. Top-hits bar panel (top 5 enriched + depleted per cohort)
##     Fig 6. Volcano panel (2×2, top 5 labelled per cohort)
##
## Inputs:
##   - Per-cohort maaslin2_kraken_<DATASET>/all_results.tsv
##   - Per-cohort tables_kraken_<DATASET>/alpha_diversity_<DATASET>.tsv
##   - Per-cohort tables_kraken_<DATASET>/beta_panel_<metric>_<DATASET>.rds
##   - cross_dataset_kraken/cross_dataset_overlap_kraken.tsv
##
## Outputs (written to lab_meeting_figs_kraken/):
##   01_alpha_diversity_panel_clean.png
##   02_cross_dataset_heatmap_clean.png
##   03_butyrate_producers_coef_clean.png
##   04_opportunistic_taxa_coef_clean.png
##   05_top_hits_bar_panel_clean.png
##   06_volcano_panel_clean.png
##
## Dependencies:
##   tidyverse, viridis, scales, ggrepel, patchwork, readxl, janitor
##
## Notes:
##   - Run KrakenAnalysis_1.R (all cohorts) and KrakenAnalysis_2.R first
## =============================================================================

library(tidyverse)
library(viridis)
library(scales)
library(ggrepel)
library(patchwork)
library(readxl)
library(janitor)

out_dir <- "D:/MasterThesis/Vik/CRASH_figs_kraken"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

theme_set(theme_minimal(base_size = 15))

Q_THRESH <- 0.25

# -------------------------------------------------------------------------
# Common labels / helpers
# -------------------------------------------------------------------------

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

maaslin_files <- c(
  prjna813705      = "D:/MasterThesis/Vik/PRJNA813705/maaslin2_kraken_prjna813705/all_results.tsv",
  kulecka_leukemia = "D:/MasterThesis/Vik/Kulecka/Leukemia/maaslin2_kraken_kulecka_leukemia/all_results.tsv",
  kulecka_lymphoma = "D:/MasterThesis/Vik/Kulecka/Lymphoma/maaslin2_kraken_kulecka_lymphoma/all_results.tsv",
  cra              = "D:/MasterThesis/Vik/CRA007433/maaslin2_kraken_cra/all_results.tsv"
)

dir_cols <- c(
  prjna813705      = "prjna813705_dir",
  kulecka_leukemia = "kulecka_leukemia_dir",
  kulecka_lymphoma = "kulecka_lymphoma_dir",
  cra              = "cra_dir"
)

q_cols <- c(
  prjna813705      = "prjna813705_q",
  kulecka_leukemia = "kulecka_leukemia_q",
  kulecka_lymphoma = "kulecka_lymphoma_q",
  cra              = "cra_q"
)

clean_taxon_name <- function(x) {
  # Kraken names already use spaces. Strip bracket-disputed-genus notation
  # for display (e.g. "[Clostridium] innocuum" -> "Clostridium innocuum")
  # so figures read cleanly. Underlying analysis keeps the bracketed form.
  x |>
    str_replace_all("\\[|\\]", "") |>
    str_squish()
}

sig_label_fun <- function(q) {
  case_when(
    q < 0.001 ~ "***",
    q < 0.01  ~ "**",
    q < 0.05  ~ "*",
    q < 0.25  ~ "\u2020",
    TRUE      ~ "ns"
  )
}

normalize_name <- function(x) {
  x <- tolower(x)
  x <- gsub("^x\\.+", "", x)
  x <- gsub("[^a-z0-9]", "", x)
  x
}

# -------------------------------------------------------------------------
# Load full MaAsLin2 results (all features, all cohorts)
# -------------------------------------------------------------------------
# MaAsLin2 mangles feature names via make.names(). We keep the mangled
# form alongside a normalized key for downstream joins and a cleaned
# display name for figures.

maaslin_all <- imap_dfr(maaslin_files, function(path, ds) {
  read_tsv(path, show_col_types = FALSE) |>
    filter(metadata == "disease") |>
    mutate(
      cohort_key = ds,
      coef_case  = if_else(value == case_map[[ds]], coef, -coef),
      is_sig     = qval < Q_THRESH,
      norm_key   = normalize_name(feature)
    ) |>
    select(feature, norm_key, cohort_key, coef_case, qval, is_sig)
})

# -------------------------------------------------------------------------
# Cohort sample counts from metadata files
# -------------------------------------------------------------------------

meta_info <- list(
  prjna813705      = list(
    path = "D:/MasterThesis/Vik/PRJNA813705/PRJNA813705_metadata.xlsx",
    case = "AML",   ctrl = "Control"
  ),
  kulecka_leukemia = list(
    path = "D:/MasterThesis/Vik/Kulecka/Leukemia/Kuleckaleukemia_metadata.xlsx",
    case = "AML",   ctrl = "Control"
  ),
  kulecka_lymphoma = list(
    path = "D:/MasterThesis/Vik/Kulecka/Lymphoma/Kuleckalymphoma_metadata.xlsx",
    case = "LN",    ctrl = "Control"
  ),
  cra              = list(
    path = "D:/MasterThesis/Vik/CRA007433/CRA007433_metadata.xlsx",
    case = "NKTCL", ctrl = "Control"
  )
)

cohort_counts <- imap_dfr(meta_info, function(info, key) {
  m <- read_excel(info$path) |>
    clean_names() |>
    filter(disease %in% c(info$case, info$ctrl))
  tibble(
    cohort_key = key,
    n_case  = sum(m$disease == info$case, na.rm = TRUE),
    n_ctrl  = sum(m$disease == info$ctrl, na.rm = TRUE),
    n_total = nrow(m)
  )
})

cohort_labels_with_n <- cohort_counts |>
  mutate(
    label = sprintf("%s\nN=%d (%d/%d)",
                    cohort_labels[cohort_key], n_total, n_case, n_ctrl)
  ) |>
  select(cohort_key, label) |>
  deframe()

# -------------------------------------------------------------------------
# 1. Alpha diversity across cohorts
# -------------------------------------------------------------------------

alpha_files <- c(
  "AML\n(PRJNA813705)" = "D:/MasterThesis/Vik/PRJNA813705/tables_kraken_prjna813705/alpha_diversity_kraken_prjna813705.tsv",
  "AML\n(Kulecka)"     = "D:/MasterThesis/Vik/Kulecka/Leukemia/tables_kraken_kulecka_leukemia/alpha_diversity_kraken_kulecka_leukemia.tsv",
  "LN\n(Kulecka)"      = "D:/MasterThesis/Vik/Kulecka/Lymphoma/tables_kraken_kulecka_lymphoma/alpha_diversity_kraken_kulecka_lymphoma.tsv",
  "NKTCL\n(CRA)"       = "D:/MasterThesis/Vik/CRA007433/tables_kraken_cra/alpha_diversity_kraken_cra.tsv"
)

alpha_combined <- imap_dfr(alpha_files, function(path, label) {
  read_tsv(path, show_col_types = FALSE) |> mutate(cohort = label)
}) |>
  mutate(
    cohort    = factor(cohort, levels = names(alpha_files)),
    metric    = factor(metric, levels = c("Observed", "Shannon", "Simpson")),
    ratio     = if_else(median_ctrl > 0, median_case / median_ctrl, NA_real_),
    sig_label = sig_label_fun(p.adj)
  )

p_alpha <- ggplot(alpha_combined, aes(x = cohort, y = ratio, fill = metric)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.72) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40", linewidth = 0.8) +
  geom_text(
    aes(label = sig_label, y = ratio + 0.02),
    position = position_dodge(width = 0.8),
    size = 5, vjust = 0
  ) +
  scale_fill_viridis_d(option = "plasma", end = 0.82) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    title    = "Alpha diversity (Kraken): median case/control ratio across cohorts",
    subtitle = "\u2020 q<0.25   * q<0.05   ** q<0.01   *** q<0.001",
    x = NULL, y = "Median case/control ratio", fill = NULL
  ) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 14),
    plot.title = element_text(face = "bold")
  )

ggsave(file.path(out_dir, "01_alpha_diversity_kraken.png"),
       p_alpha, width = 12, height = 6, dpi = 300)

# -------------------------------------------------------------------------
# 2. Cross-dataset heatmap
# -------------------------------------------------------------------------

overlap <- read_tsv(
  "D:/MasterThesis/Vik/cross_dataset_kraken/cross_dataset_overlap_kraken.tsv",
  show_col_types = FALSE
)

heatmap_features <- overlap |>
  filter(n_datasets >= 3, direction_consistent == TRUE) |>
  rowwise() |>
  mutate(min_q = suppressWarnings(min(c_across(all_of(unname(q_cols))), na.rm = TRUE))) |>
  ungroup() |>
  mutate(min_q = if_else(is.infinite(min_q), NA_real_, min_q))

feat_order <- heatmap_features |>
  arrange(consensus_direction, desc(n_datasets), min_q) |>
  pull(feature)

heatmap_data <- heatmap_features |>
  select(feature, all_of(unname(dir_cols)), all_of(unname(q_cols))) |>
  pivot_longer(
    cols = -feature,
    names_to = c("cohort_key", ".value"),
    names_pattern = "(.+)_(dir|q)"
  ) |>
  mutate(
    cohort = factor(cohort_labels_with_n[cohort_key],
                    levels = unname(cohort_labels_with_n)),
    feature_clean = clean_taxon_name(feature),
    feature_clean = factor(feature_clean, levels = clean_taxon_name(feat_order)),
    direction_plot = case_when(
      dir == "depleted" ~ "Depleted in cancer",
      dir == "enriched" ~ "Enriched in cancer",
      TRUE ~ "Not a hit"
    ),
    direction_plot = factor(direction_plot,
                            levels = c("Depleted in cancer",
                                       "Enriched in cancer", "Not a hit")),
    sig_label = case_when(
      is.na(q)  ~ "",
      q < 0.001 ~ "***",
      q < 0.01  ~ "**",
      q < 0.05  ~ "*",
      q < 0.25  ~ "\u2020",
      TRUE      ~ ""
    )
  )

p_heat <- ggplot(heatmap_data, aes(x = cohort, y = feature_clean, fill = direction_plot)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  geom_text(aes(label = sig_label), colour = "white",
            fontface = "bold", size = 3.2, vjust = 0.55) +
  scale_fill_manual(
    values = c(
      "Depleted in cancer" = "#2166ac",
      "Enriched in cancer" = "#d6604d",
      "Not a hit" = "grey88"
    )
  ) +
  labs(
    title = "Species consistently altered across \u22653 cohorts (Kraken)",
    subtitle = paste(
      "Rows ordered by n cohorts, then by minimum q. \u00A0",
      "Tile markers:\u2020 q<0.25  * q<0.05  ** q<0.01  *** q<0.001"
    ),
    x = NULL, y = NULL, fill = NULL
  ) +
  theme(
    legend.position = "top",
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 8.5, face = "italic"),
    axis.text.x = element_text(size = 11, lineheight = 0.9),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 11, colour = "grey30")
  )

# Height scales with number of rows
fig2_height <- max(8, 0.22 * length(feat_order) + 3)
ggsave(file.path(out_dir, "02_cross_dataset_heatmap_kraken.png"),
       p_heat, width = 10, height = fig2_height, dpi = 300)

# -------------------------------------------------------------------------
# 3/4. Curated functional group heatmaps
# -------------------------------------------------------------------------
# Kraken species-name variants (brackets, renamed genera) are included so
# the match finds whichever form your Bracken DB uses. We then dedupe by
# normalized key so only one row per biological taxon appears.

butyrate_key <- c(
  "Roseburia hominis", "Roseburia intestinalis", "Roseburia faecis",
  "Roseburia inulinivorans",
  "Eubacterium rectale", "[Eubacterium] rectale",
  "Faecalibacterium prausnitzii",
  "Coprococcus catus", "Coprococcus comes", "Coprococcus eutactus",
  "Agathobaculum butyriciproducens",
  "Lachnospira eligens", "[Eubacterium] eligens",
  "Ruminococcus callidus"
)

opportunistic_key <- c(
  "[Clostridium] innocuum", "Clostridium innocuum",
  "Enterocloster bolteae", "[Clostridium] bolteae",
  "Eggerthella lenta",
  "Hungatella hathewayi", "[Clostridium] hathewayi",
  "[Clostridium] symbiosum", "Clostridium symbiosum",
  "Escherichia coli",
  "Mediterraneibacter gnavus", "Ruminococcus gnavus",
  "Erysipelatoclostridium ramosum", "[Clostridium] ramosum"
)

build_coef_heatmap <- function(feature_list, title_text, subtitle_text,
                               out_file, fig_height = 7) {
  
  keys <- normalize_name(feature_list)
  
  # Pull every cohort x feature combination. Use a map of key -> display
  # name (first matched variant) to collapse synonyms.
  display_name <- tibble(
    norm_key = normalize_name(feature_list),
    display  = feature_list
  ) |>
    group_by(norm_key) |>
    slice(1) |>
    ungroup()
  
  dat <- maaslin_all |>
    filter(norm_key %in% keys) |>
    group_by(norm_key, cohort_key) |>
    slice(1) |>   # dedupe in case MaAsLin2 mangling produced duplicates
    ungroup() |>
    left_join(display_name, by = "norm_key") |>
    complete(norm_key = keys, cohort_key = names(maaslin_files)) |>
    left_join(display_name, by = "norm_key", suffix = c("", "_disp")) |>
    mutate(display = coalesce(display, display_disp)) |>
    select(-display_disp)
  
  if (all(is.na(dat$coef_case))) {
    message("No matches found for: ", title_text)
    return(invisible(NULL))
  }
  
  rank_df <- dat |>
    group_by(display) |>
    summarise(
      n_sig         = sum(is_sig, na.rm = TRUE),
      min_q         = suppressWarnings(min(qval, na.rm = TRUE)),
      mean_abs_coef = mean(abs(coef_case), na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(min_q = if_else(is.infinite(min_q), NA_real_, min_q)) |>
    arrange(desc(n_sig), min_q, desc(mean_abs_coef))
  
  feat_order <- rank_df$display
  
  max_abs <- max(abs(dat$coef_case), na.rm = TRUE)
  
  dat <- dat |>
    mutate(
      cohort = factor(cohort_labels_with_n[cohort_key],
                      levels = unname(cohort_labels_with_n)),
      feature_clean = clean_taxon_name(display),
      feature_clean = factor(feature_clean,
                             levels = clean_taxon_name(rev(feat_order))),
      is_sig_safe = coalesce(is_sig, FALSE),
      sig_marker  = if_else(is_sig_safe & !is.na(coef_case), "*", ""),
      label = if_else(
        is.na(coef_case),
        "\u2014",
        paste0(sprintf("%.1f", coef_case), sig_marker)
      ),
      label_col = if_else(
        is.na(coef_case) | abs(coef_case) < 0.6 * max_abs,
        "black", "white"
      )
    )
  
  p <- ggplot(dat, aes(x = cohort, y = feature_clean, fill = coef_case)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = label, colour = label_col), size = 3.6,
              show.legend = FALSE) +
    scale_colour_identity() +
    scale_fill_gradient2(
      low = "#2166ac", mid = "white", high = "#d6604d",
      midpoint = 0, na.value = "grey88",
      limits   = c(-max_abs, max_abs),
      name     = "MaAsLin2\ncoefficient\n(case-oriented)"
    ) +
    labs(title = title_text, subtitle = subtitle_text, x = NULL, y = NULL) +
    theme(
      axis.text.y = element_text(face = "italic"),
      panel.grid  = element_blank(),
      legend.position = "right",
      plot.title = element_text(face = "bold")
    )
  
  ggsave(out_file, p, width = 10, height = fig_height, dpi = 300)
  invisible(p)
}

build_coef_heatmap(
  butyrate_key,
  title_text    = "Butyrate producers (Kraken): MaAsLin2 effect sizes",
  subtitle_text = "Rows ordered by n cohorts sig then min q; * marks q < 0.25",
  out_file      = file.path(out_dir, "03_butyrate_producers_kraken.png"),
  fig_height    = 7
)

build_coef_heatmap(
  opportunistic_key,
  title_text    = "Opportunistic taxa (Kraken): MaAsLin2 effect sizes",
  subtitle_text = "Rows ordered by n cohorts sig then min q; * marks q < 0.25",
  out_file      = file.path(out_dir, "04_opportunistic_taxa_kraken.png"),
  fig_height    = 6
)

# -------------------------------------------------------------------------
# 5. Top-hits bar panel per cohort
# -------------------------------------------------------------------------

TOP_N_BARS <- 5

top_bars_df <- maaslin_all |>
  filter(!is.na(coef_case), !is.na(qval), qval < Q_THRESH) |>
  group_by(cohort_key) |>
  group_modify(~ {
    ups   <- .x |> arrange(desc(coef_case)) |> slice_head(n = TOP_N_BARS)
    downs <- .x |> arrange(coef_case)       |> slice_head(n = TOP_N_BARS)
    bind_rows(ups, downs) |> distinct(feature, .keep_all = TRUE)
  }) |>
  ungroup() |>
  mutate(
    cohort = factor(cohort_labels_with_n[cohort_key],
                    levels = unname(cohort_labels_with_n)),
    feature_clean = clean_taxon_name(feature),
    direction_plot = if_else(coef_case > 0, "Enriched in cancer", "Depleted in cancer"),
    direction_plot = factor(direction_plot,
                            levels = c("Depleted in cancer", "Enriched in cancer")),
    sig_label = case_when(
      qval < 0.001 ~ "***",
      qval < 0.01  ~ "**",
      qval < 0.05  ~ "*",
      qval < 0.25  ~ "\u2020",
      TRUE         ~ ""
    )
  )

build_top_bars <- function(df, cohort_name) {
  d <- df |>
    filter(cohort == cohort_name) |>
    mutate(feature_clean = fct_reorder(feature_clean, coef_case))
  
  if (nrow(d) == 0) {
    return(
      ggplot() +
        annotate("text", x = 0, y = 0, label = "No q<0.25 hits", colour = "grey50") +
        labs(title = cohort_name) + theme_void() +
        theme(plot.title = element_text(face = "bold", size = 13))
    )
  }
  
  x_max <- max(abs(d$coef_case), na.rm = TRUE) * 1.18
  
  ggplot(d, aes(x = coef_case, y = feature_clean, fill = direction_plot)) +
    geom_col(width = 0.75) +
    geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3) +
    geom_text(
      aes(label = sig_label,
          hjust = if_else(coef_case > 0, -0.25, 1.25)),
      size = 3.8, fontface = "bold", colour = "grey20"
    ) +
    scale_fill_manual(
      values = c("Depleted in cancer" = "#2166ac",
                 "Enriched in cancer" = "#d6604d"),
      drop = FALSE
    ) +
    scale_x_continuous(limits = c(-x_max, x_max)) +
    labs(title = cohort_name,
         x = "MaAsLin2 coefficient (case-oriented)",
         y = NULL) +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      legend.position = "none",
      axis.text.y = element_text(face = "italic", size = 9),
      panel.grid.minor = element_blank()
    )
}

bar_plots <- lapply(levels(top_bars_df$cohort),
                    function(ch) build_top_bars(top_bars_df, ch))

p_top_bars <- wrap_plots(bar_plots, ncol = 2) +
  plot_annotation(
    title = sprintf(
      "Top %d enriched and top %d depleted species per cohort (Kraken, MaAsLin2 q<%.2f)",
      TOP_N_BARS, TOP_N_BARS, Q_THRESH
    ),
    subtitle = "Bar = case-oriented coefficient. Tip marker: \u2020 q<0.25  * q<0.05  ** q<0.01  *** q<0.001",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11, colour = "grey30")
    )
  ) &
  theme(plot.margin = margin(4, 6, 4, 6))

ggsave(file.path(out_dir, "05_top_hits_bars_kraken.png"),
       p_top_bars, width = 14, height = 10, dpi = 300)

# -------------------------------------------------------------------------
# 6. Volcano panel (2 x 2)
# -------------------------------------------------------------------------

VOLCANO_LABEL_N <- 5

volcano_df <- maaslin_all |>
  filter(!is.na(qval), qval > 0) |>
  mutate(
    cohort    = factor(cohort_labels[cohort_key], levels = unname(cohort_labels)),
    neglog10q = -log10(qval),
    direction_plot = case_when(
      qval >= Q_THRESH ~ "ns",
      coef_case > 0    ~ "Enriched in cancer",
      coef_case < 0    ~ "Depleted in cancer",
      TRUE             ~ "ns"
    ),
    direction_plot = factor(direction_plot,
                            levels = c("Depleted in cancer",
                                       "Enriched in cancer", "ns")),
    feature_clean = clean_taxon_name(feature)
  )

volcano_labels <- volcano_df |>
  filter(qval < Q_THRESH, !is.na(coef_case)) |>
  group_by(cohort) |>
  arrange(desc(abs(coef_case) * neglog10q)) |>
  slice_head(n = VOLCANO_LABEL_N) |>
  ungroup()

x_max <- max(abs(volcano_df$coef_case), na.rm = TRUE)
x_max <- ceiling(x_max * 10) / 10
y_max <- max(volcano_df$neglog10q, na.rm = TRUE) * 1.05

build_volcano <- function(df, labels_df, cohort_name) {
  d  <- df        |> filter(cohort == cohort_name)
  lb <- labels_df |> filter(cohort == cohort_name)
  
  ggplot(d, aes(x = coef_case, y = neglog10q)) +
    geom_vline(xintercept = 0, linetype = "solid", colour = "grey70", linewidth = 0.4) +
    geom_hline(yintercept = -log10(Q_THRESH), linetype = "dashed",
               colour = "grey40", linewidth = 0.5) +
    geom_point(aes(colour = direction_plot), alpha = 0.75, size = 2) +
    geom_text_repel(
      data = lb,
      aes(label = feature_clean),
      size = 3, fontface = "italic",
      min.segment.length = 0, segment.size = 0.3, segment.colour = "grey60",
      max.overlaps = Inf, box.padding = 0.35, force = 2
    ) +
    scale_colour_manual(
      values = c("Depleted in cancer" = "#2166ac",
                 "Enriched in cancer" = "#d6604d",
                 "ns"                 = "grey75"),
      drop = FALSE
    ) +
    scale_x_continuous(limits = c(-x_max, x_max)) +
    scale_y_continuous(limits = c(0, y_max), expand = expansion(mult = c(0, 0.02))) +
    labs(title = cohort_name,
         x = "MaAsLin2 coefficient (case-oriented)",
         y = expression(-log[10](q))) +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      legend.position = "none",
      panel.grid.minor = element_blank()
    )
}

plots <- lapply(levels(volcano_df$cohort),
                function(ch) build_volcano(volcano_df, volcano_labels, ch))

p_volcano <- wrap_plots(plots, ncol = 2) +
  plot_annotation(
    title = "Species-level differential abundance: Kraken2/Bracken MaAsLin2 across cohorts",
    subtitle = paste0(
      "Dashed line: q = ", Q_THRESH,
      ". Top ", VOLCANO_LABEL_N, " hits per cohort labelled (by |coef| x -log10 q)."
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(size = 11, colour = "grey30")
    )
  ) &
  theme(plot.margin = margin(4, 6, 4, 6))

ggsave(file.path(out_dir, "06_volcano_panel_kraken.png"),
       p_volcano, width = 13, height = 11, dpi = 300)

cat("Kraken lab meeting figures saved to:", out_dir, "\n")
