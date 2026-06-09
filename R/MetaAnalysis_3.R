## =============================================================================
## MetaAnalysis_3.R
## CRASH — Cancer-associated gut micRobiome Analysis of Strains in
##         Haematological patients
##
##
## Purpose:
##   Cross-cohort publication figures assembling outputs from MetaAnalysis_1.R
##   and MetaAnalysis_2.R into a unified visual summary. Produces 8 figures:
##
##     Fig 1. Alpha diversity panel — median case/control ratio per metric
##            (Observed, Shannon, Simpson) across all 4 cohorts
##     Fig 2. Cross-dataset heatmap — top 5 enriched + top 5 depleted species
##            present in ≥3 cohorts, with q-value tier asterisks
##     Fig 3. Butyrate producer coefficient heatmap — MaAsLin2 effect sizes
##            for key commensals across cohorts; rows ranked by evidence
##     Fig 4. Opportunistic taxa coefficient heatmap — same layout for
##            enriched opportunists
##     Fig 5. Oral-origin taxa heatmap — HOMD-anchored, data-driven panel
##            for oral-gut translocation signal (NKTCL-specific)
##     Fig 6. Random forest importance panel — feature importance scores
##            from per-cohort RF classifiers
##     Fig 7. Volcano panel — 2×2 grid of MaAsLin2 landscapes per cohort;
##            shared symmetric x-axis; top hits labelled
##     Fig 8. Beta diversity panel — PCoA 2×4 grid (Bray-Curtis + Jaccard
##            × 4 cohorts) with PERMANOVA/betadisper stats in subtitles
##
## Inputs:
##   - Per-cohort maaslin2_metaphlan_<DATASET>/all_results.tsv
##   - Per-cohort tables_metaphlan_<DATASET>/alpha_diversity_<DATASET>.tsv
##   - Per-cohort tables_metaphlan_<DATASET>/beta_panel_<metric>_<DATASET>.rds
##   - cross_dataset/cross_dataset_overlap_metaphlan.tsv (from MetaAnalysis_2.R)
##   - Per-cohort metadata Excel files (for sample count annotations)
##
## Outputs (written to CRASH_figs/):
##   01_alpha_diversity_panel_clean.png
##   02_cross_dataset_heatmap_clean.png
##   03_butyrate_producers_coef_clean.png
##   04_opportunistic_taxa_coef_clean.png
##   05_oral_taxa_heatmap_clean.png
##   06_random_forest_importance_clean.png
##   07_volcano_panel_clean.png
##   08_beta_diversity_panel_clean.png
##
## Dependencies:
##   tidyverse, viridis, scales, ggrepel, patchwork, readxl, janitor
##
## Notes:
##   - Run MetaAnalysis_1.R (all cohorts) and MetaAnalysis_2.R before this
##   - Beta panel (Fig 8) requires .rds payloads from MetaAnalysis_1.R;
##     missing files trigger a warning and skip the figure gracefully
##   - Beta ellipses use 68% (1 SD) level — tighter than default 95% —
##     paired with PERMANOVA stats in subtitles for honest interpretation
## =============================================================================

library(tidyverse)
library(viridis)
library(scales)
library(ggrepel)
library(patchwork)
library(readxl)
library(janitor)

out_dir <- "D:/MasterThesis/Vik/CRASH_figs_metaphlan"
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
  prjna813705      = "D:/MasterThesis/Vik/PRJNA813705/maaslin2_metaphlan_prjna813705/all_results.tsv",
  kulecka_leukemia = "D:/MasterThesis/Vik/Kulecka/Leukemia/maaslin2_metaphlan_kulecka_leukemia/all_results.tsv",
  kulecka_lymphoma = "D:/MasterThesis/Vik/Kulecka/Lymphoma/maaslin2_metaphlan_kulecka_lymphoma/all_results.tsv",
  cra              = "D:/MasterThesis/Vik/CRA007433/maaslin2_metaphlan_cra/all_results.tsv"
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

clean_taxon_name <- function(x) str_replace_all(x, "_", " ")

sig_label_fun <- function(q) {
  case_when(
    q < 0.001 ~ "***",
    q < 0.01  ~ "**",
    q < 0.05  ~ "*",
    q < 0.25  ~ "\u2020",
    TRUE      ~ "ns"
  )
}

# Load all MaAsLin2 disease associations across cohorts, case-oriented.
# Every feature tested in a given cohort (i.e. passing its prevalence filter)
# shows up here — no q threshold applied at this stage.
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

# Cohort sample counts, read from per-cohort metadata files.
# Used to annotate x-axis of the cross-dataset heatmap so the reader knows
# the statistical weight behind each cohort.
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

# Enriched cohort labels: original label + sample counts on a third line.
cohort_labels_with_n <- cohort_counts |>
  mutate(
    label = sprintf(
      "%s\nN=%d (%d/%d)",
      cohort_labels[cohort_key], n_total, n_case, n_ctrl
    )
  ) |>
  select(cohort_key, label) |>
  deframe()

cat("Cohort sample counts:\n")
print(cohort_counts)

# -------------------------------------------------------------------------
# 1. Alpha diversity across cohorts
# -------------------------------------------------------------------------

alpha_files <- c(
  "AML\n(PRJNA813705)" = "D:/MasterThesis/Vik/PRJNA813705/tables_metaphlan_prjna813705/alpha_diversity_prjna813705.tsv",
  "AML\n(Kulecka)"     = "D:/MasterThesis/Vik/Kulecka/Leukemia/tables_metaphlan_kulecka_leukemia/alpha_diversity_kulecka_leukemia.tsv",
  "LN\n(Kulecka)"      = "D:/MasterThesis/Vik/Kulecka/Lymphoma/tables_metaphlan_kulecka_lymphoma/alpha_diversity_kulecka_lymphoma.tsv",
  "NKTCL\n(CRA)"       = "D:/MasterThesis/Vik/CRA007433/tables_metaphlan_cra/alpha_diversity_cra.tsv"
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
    title    = "Alpha diversity: median case/control ratio across cohorts",
    subtitle = "\u2020 q<0.25   * q<0.05   ** q<0.01   *** q<0.001",
    x = NULL,
    y = "Median case/control ratio",
    fill = NULL
  ) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 14),
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(out_dir, "01_alpha_diversity_panel_clean.png"),
  p_alpha, width = 12, height = 6, dpi = 800
)

# -------------------------------------------------------------------------
# 2. Cross-dataset overlap heatmap (from Script 2 output)
# -------------------------------------------------------------------------

overlap <- read_tsv(
  "D:/MasterThesis/Vik/cross_dataset/cross_dataset_overlap_metaphlan.tsv",
  show_col_types = FALSE
)

heatmap_features <- overlap |>
  filter(n_datasets >= 3, direction_consistent == TRUE) |>
  rowwise() |>
  mutate(min_q = suppressWarnings(min(c_across(all_of(unname(q_cols))), na.rm = TRUE))) |>
  ungroup() |>
  mutate(min_q = if_else(is.infinite(min_q), NA_real_, min_q))

# Top 5 per direction, ranked by n_datasets then min_q
heatmap_features <- bind_rows(
  heatmap_features |> filter(consensus_direction == "depleted") |>
    arrange(desc(n_datasets), min_q) |> slice_head(n = 5),
  heatmap_features |> filter(consensus_direction == "enriched") |>
    arrange(desc(n_datasets), min_q) |> slice_head(n = 5)
)

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
    cohort = factor(
      cohort_labels_with_n[cohort_key],
      levels = unname(cohort_labels_with_n)
    ),
    feature_clean = clean_taxon_name(feature),
    feature_clean = factor(feature_clean, levels = clean_taxon_name(feat_order)),
    direction_plot = case_when(
      dir == "depleted" ~ "Depleted in cancer",
      dir == "enriched" ~ "Enriched in cancer",
      TRUE ~ "Not a hit"
    ),
    direction_plot = factor(
      direction_plot,
      levels = c("Depleted in cancer", "Enriched in cancer", "Not a hit")
    ),
    # Asterisk overlay for q-value tier. Grey "Not a hit" tiles have NA q
    # and therefore no marker. The dagger tier exists because all features
    # in this table already pass q<0.25 as high-confidence hits.
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
  geom_text(
    aes(label = sig_label),
    colour = "white", fontface = "bold", size = 3.2, vjust = 0.55
  ) +
  scale_fill_manual(
    values = c(
      "Depleted in cancer" = "#2166ac",
      "Enriched in cancer" = "#d6604d",
      "Not a hit" = "grey88"
    )
  ) +
  labs(
    title = "Top 5 enriched + top 5 depleted species across \u22653 cohorts",
    subtitle = paste(
      "Ranked by n cohorts then minimum q. \u00A0",
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

ggsave(
  file.path(out_dir, "02_cross_dataset_heatmap_clean.png"),
  p_heat, width = 10, height = 6, dpi = 800
)

# -------------------------------------------------------------------------
# 3/4. Coefficient heatmaps — rebuilt from full per-cohort MaAsLin2 results
# -------------------------------------------------------------------------

build_coef_heatmap <- function(feature_list, title_text, subtitle_text,
                               out_file, fig_height = 7, fig_width = 10) {
  
  # Start with all (feature x cohort) pairs so features missing from any
  # cohort still render a tile (as NA / em-dash).
  dat <- maaslin_all |>
    filter(feature %in% feature_list) |>
    complete(feature = feature_list, cohort_key = names(maaslin_files))
  
  # Rank by evidence strength: n cohorts sig, then min q, then mean |coef|.
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
  
  max_abs <- max(abs(dat$coef_case), na.rm = TRUE)
  
  dat <- dat |>
    mutate(
      cohort = factor(cohort_labels[cohort_key], levels = unname(cohort_labels)),
      feature_clean = clean_taxon_name(feature),
      feature_clean = factor(
        feature_clean,
        levels = clean_taxon_name(rev(feat_order))
      ),
      is_sig_safe = coalesce(is_sig, FALSE),
      sig_marker = case_when(
        !is_sig_safe | is.na(coef_case) ~ "",
        qval < 0.001 ~ "***",
        qval < 0.01  ~ "**",
        qval < 0.05  ~ "*",
        qval < 0.25  ~ "\u2020",
        TRUE         ~ ""
      ),
      label = if_else(
        is.na(coef_case),
        "\u2014",
        paste0(sprintf("%.1f", coef_case), sig_marker)
      ),
      # White text on the most saturated tiles, black elsewhere.
      label_col = if_else(
        is.na(coef_case) | abs(coef_case) < 0.6 * max_abs,
        "black", "white"
      )
    )
  
  p <- ggplot(dat, aes(x = cohort, y = feature_clean, fill = coef_case)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = label, colour = label_col), size = 3.6, show.legend = FALSE) +
    scale_colour_identity() +
    scale_fill_gradient2(
      low      = "#2166ac",
      mid      = "white",
      high     = "#d6604d",
      midpoint = 0,
      na.value = "grey88",
      limits   = c(-max_abs, max_abs),
      name     = "MaAsLin2\ncoefficient\n(case-oriented)"
    ) +
    labs(
      title    = title_text,
      subtitle = subtitle_text,
      x = NULL, y = NULL
    ) +
    theme(
      axis.text.y = element_text(face = "italic"),
      panel.grid  = element_blank(),
      legend.position = "right",
      plot.title    = element_text(face = "bold"),
      plot.subtitle = element_text(size = 10, colour = "grey30"),
      plot.title.position = "plot"
    )
  
  ggsave(out_file, p, width = fig_width, height = fig_height, dpi = 800)
  invisible(p)
}

butyrate_key <- c(
  "Roseburia_hominis", "Roseburia_intestinalis", "Roseburia_faecis",
  "Roseburia_inulinivorans", "Eubacterium_rectale", "Faecalibacterium_prausnitzii",
  "Coprococcus_catus", "Coprococcus_comes", "Coprococcus_eutactus",
  "Agathobaculum_butyriciproducens", "Lachnospira_eligens", "Ruminococcus_callidus"
)

opportunistic_key <- c(
  "Clostridium_innocuum", "Enterocloster_bolteae", "Eggerthella_lenta",
  "Hungatella_hathewayi", "Clostridium_symbiosum", "Escherichia_coli",
  "Ruminococcus_gnavus", "Erysipelatoclostridium_ramosum"
)

build_coef_heatmap(
  butyrate_key,
  title_text    = "Butyrate producers: MaAsLin2 effect sizes",
  subtitle_text = "Rows ordered by n cohorts sig then min q; \u2020 q<0.25  * q<0.05  ** q<0.01  *** q<0.001",
  out_file      = file.path(out_dir, "03_butyrate_producers_coef_clean.png"),
  fig_height    = 7
)

build_coef_heatmap(
  opportunistic_key,
  title_text    = "Enriched opportunistic taxa: MaAsLin2 effect sizes",
  subtitle_text = "Rows ordered by n cohorts sig then min q; \u2020 q<0.25  * q<0.05  ** q<0.01  *** q<0.001",
  out_file      = file.path(out_dir, "04_opportunistic_taxa_coef_clean.png"),
  fig_height    = 6
)

# -------------------------------------------------------------------------
# 5. Oral-origin taxa enriched in cancer gut — HOMD-anchored, data-driven
# -------------------------------------------------------------------------
#
# Universe: all named HOMD species at oral abundance tier High or Medium
# (n = 159 unique species; downloaded 2026-05-04 from homd.org).
#
# Filter: retain only species that are (a) detected by MetaPhlAn4 in any
# cohort at ≥10% prevalence, AND (b) significantly ENRICHED (coef_case > 0,
# q < Q_THRESH) in at least one cohort.
#
# Rationale: restricting to enrichment only captures oral-gut translocation
# (gain of oral taxa in the cancer gut). Depletion of oral taxa would be
# mechanistically ambiguous and is not the focus of this panel.
# The HOMD universe avoids selection bias from literature curation;
# only species with actual signal in these data survive the filter.

homd_oral_key <- c(
  "Abiotrophia_defectiva", "Actinomyces_dentalis", "Actinomyces_gerencseriae", "Actinomyces_graevenitzii",
  "Actinomyces_johnsonii", "Actinomyces_naeslundii", "Actinomyces_oris", "Actinomyces_viscosus",
  "Aggregatibacter_aphrophilus", "Aggregatibacter_kilianii", "Aggregatibacter_segnis", "Alloprevotella_rava",
  "Alloprevotella_tannerae", "Anaeroglobus_geminatus", "Arachnia_propionica", "Arachnia_rubra",
  "Bulleidia_extructa", "Campylobacter_concisus", "Campylobacter_gracilis", "Campylobacter_rectus",
  "Campylobacter_showae", "Capnocytophaga_gingivalis", "Capnocytophaga_granulosa", "Capnocytophaga_leadbetteri",
  "Capnocytophaga_ochracea", "Capnocytophaga_periodontitidis", "Capnocytophaga_sputigena", "Cardiobacterium_hominis",
  "Cardiobacterium_valvarum", "Catonella_morbi", "Corynebacterium_durum", "Corynebacterium_matruchotii",
  "Dialister_invisus", "Dialister_pneumosintes", "Eikenella_corrodens", "Filifactor_alocis",
  "Fusobacterium_animalis", "Fusobacterium_hwasookii", "Fusobacterium_periodonticum", "Fusobacterium_polymorphum",
  "Fusobacterium_vincentii", "Gallibacter_brachus", "Gemella_haemolysans", "Gemella_morbillorum",
  "Gemella_sanguinis", "Granulicatella_adiacens", "Granulicatella_elegans", "Haemophilus_haemolyticus",
  "Haemophilus_parahaemolyticus", "Haemophilus_parainfluenzae", "Haemophilus_paraphrohaemolyticus", "Haemophilus_pittmaniae",
  "Haemophilus_sputorum", "Hornefia_nodatum", "Hoylesella_loescheii", "Hoylesella_nanceiensis",
  "Hoylesella_pleuritidis", "Hoylesella_saccharolytica", "Hoylesella_shahii", "Johnsonella_ignava",
  "Kingella_denitrificans", "Kingella_negevensis", "Kingella_oralis", "Lachnoanaerobaculum_orale",
  "Lachnoanaerobaculum_saburreum", "Lachnoanaerobaculum_umeaense", "Lancefieldella_parvula", "Lancefieldella_rimae",
  "Lautropia_dentalis", "Lautropia_mirabilis", "Leptotrichia_buccalis", "Leptotrichia_hofstadii",
  "Leptotrichia_hongkongensis", "Leptotrichia_wadei", "Megasphaera_micronuciformis", "Metamycoplasma_salivarium",
  "Mogibacterium_diversum", "Mogibacterium_timidum", "Moraxella_catarrhalis", "Neisseria_bacilliformis",
  "Neisseria_cinerea", "Neisseria_elongata", "Neisseria_flava", "Neisseria_flavescens",
  "Neisseria_macacae", "Neisseria_mucosa", "Neisseria_oralis", "Neisseria_perflava",
  "Neisseria_sicca", "Neisseria_subflava", "Oribacterium_asaccharolyticum", "Oribacterium_parvum",
  "Oribacterium_sinus", "Parvimonas_micra", "Peptidiphaga_gingivicola", "Peptoanaerobacter_yurii",
  "Peptostreptococcus_stomatis", "Porphyromonas_catoniae", "Porphyromonas_endodontalis", "Porphyromonas_gingivalis",
  "Porphyromonas_pasteri", "Prevotella_aurantiaca", "Prevotella_denticola", "Prevotella_histicola",
  "Prevotella_intermedia", "Prevotella_jejuni", "Prevotella_melaninogenica", "Prevotella_nigrescens",
  "Prevotella_pallens", "Prevotella_scopos", "Prevotella_veroralis", "Prevotella_vespertina",
  "Pseudoleptotrichia_goodfellowii", "Rothia_aeria", "Rothia_dentocariosa", "Rothia_mucilaginosa",
  "Scardovia_wiggsiae", "Schaalia_georgiae", "Schaalia_odontolytica", "Segatella_baroniae",
  "Segatella_maculosa", "Segatella_oris", "Segatella_oulorum", "Segatella_salivae",
  "Selenomonas_infelix", "Selenomonas_noxia", "Selenomonas_sputigena", "Simonsiella_muelleri",
  "Solobacterium_moorei", "Stomatobaculum_longum", "Streptococcus_anginosus", "Streptococcus_australis",
  "Streptococcus_constellatus", "Streptococcus_gordonii", "Streptococcus_infantis", "Streptococcus_intermedius",
  "Streptococcus_mitis", "Streptococcus_mutans", "Streptococcus_oralis", "Streptococcus_peroris",
  "Streptococcus_rubneri", "Streptococcus_salivarius", "Streptococcus_sanguinis", "Streptococcus_vestibularis",
  "Tannerella_forsythia", "Tannerella_serpentiformis", "Treponema_denticola", "Treponema_lecithinolyticum",
  "Treponema_maltophilum", "Treponema_medium", "Treponema_socranskii", "Treponema_vincentii",
  "Veillonella_atypica", "Veillonella_dispar", "Veillonella_nakazawae", "Veillonella_parvula",
  "Veillonella_rogosae", "Veillonella_tobetsuensis"
)

# Species significantly altered (either direction) in ≥1 cohort from HOMD universe
oral_sig_altered <- maaslin_all |>
  filter(feature %in% homd_oral_key,
         qval < Q_THRESH) |>
  pull(feature) |>
  unique()

cat("HOMD oral taxa (High/Med) significantly altered in ≥1 cohort:",
    length(oral_sig_altered), "\n")

if (length(oral_sig_altered) > 0) {
  build_coef_heatmap(
    oral_sig_altered,
    title_text    = "Oral-origin taxa significantly altered in \u22651 cancer cohort",
    subtitle_text = "HOMD High/Medium oral species (n=159), q<0.25 in \u22651 cohort;* marks q<0.25/nTile markers:\u2020 q<0.25  * q<0.05  ** q<0.01  *** q<0.001",
    out_file      = file.path(out_dir, "05_oral_taxa_coef_clean.png"),
    fig_height    = max(4, length(oral_sig_altered) * 0.45 + 2),
    fig_width     = 12
  )
} else {
  message("No HOMD oral taxa significantly altered in any cohort — Figure 5 skipped.")
}

# -------------------------------------------------------------------------
# 6. Top-hits bar panel — top 5 enriched + top 5 depleted per cohort (2 x 2)
# -------------------------------------------------------------------------
#
# For each cohort: pick the 5 features with largest positive coef and the
# 5 with largest negative coef (both at q < Q_THRESH). Bar length = case-
# oriented coefficient; tip asterisk = q-value tier. Within each cohort
# panel, features are ordered by coefficient so depleted species stack at
# the bottom and enriched at the top.
#
# X-axis is free per cohort but symmetric around 0 within each panel, so
# 0 is always centered. This keeps per-cohort detail visible without
# distorting the direction comparison.

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
    cohort = factor(
      cohort_labels_with_n[cohort_key],
      levels = unname(cohort_labels_with_n)
    ),
    feature_clean = clean_taxon_name(feature),
    direction_plot = if_else(coef_case > 0, "Enriched in cancer", "Depleted in cancer"),
    direction_plot = factor(
      direction_plot,
      levels = c("Depleted in cancer", "Enriched in cancer")
    ),
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
        labs(title = cohort_name) +
        theme_void() +
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
      values = c(
        "Depleted in cancer" = "#2166ac",
        "Enriched in cancer" = "#d6604d"
      ),
      drop = FALSE
    ) +
    scale_x_continuous(limits = c(-x_max, x_max)) +
    labs(
      title = cohort_name,
      x = "MaAsLin2 coefficient (case-oriented)",
      y = NULL
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      legend.position = "none",
      axis.text.y = element_text(face = "italic", size = 9),
      panel.grid.minor = element_blank()
    )
}

bar_plots <- lapply(
  levels(top_bars_df$cohort),
  function(ch) build_top_bars(top_bars_df, ch)
)

p_top_bars <- wrap_plots(bar_plots, ncol = 2) +
  plot_annotation(
    title = sprintf(
      "Top %d enriched and top %d depleted species per cohort (MaAsLin2, q<%.2f)",
      TOP_N_BARS, TOP_N_BARS, Q_THRESH
    ),
    subtitle = "Bar = case-oriented coefficient. Tip marker: \u2020 q<0.25  * q<0.05  ** q<0.01  *** q<0.001",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11, colour = "grey30")
    )
  ) &
  theme(plot.margin = margin(4, 6, 4, 6))

ggsave(
  file.path(out_dir, "06_top_hits_bars_clean.png"),
  p_top_bars, width = 14, height = 10, dpi = 800
)

# -------------------------------------------------------------------------
# 7. Volcano panel — full MaAsLin2 DA landscape per cohort (2 x 2)
# -------------------------------------------------------------------------
#
# x: case-oriented MaAsLin2 coefficient (log abundance scale)
# y: -log10(qval)  — q rather than p because Q<0.25 is the primary threshold
# Shared symmetric x-axis across the 4 panels so effect-size magnitude
# is comparable cohort-to-cohort.
#
# Labels: top N features per cohort ranked by |coef| * -log10(qval),
# restricted to q < Q_THRESH so we don't tag noise.

VOLCANO_LABEL_N <- 5

volcano_df <- maaslin_all |>
  filter(!is.na(qval), qval > 0) |>
  mutate(
    cohort = factor(cohort_labels[cohort_key], levels = unname(cohort_labels)),
    neglog10q = -log10(qval),
    direction_plot = case_when(
      qval >= Q_THRESH ~ "ns",
      coef_case > 0    ~ "Enriched in cancer",
      coef_case < 0    ~ "Depleted in cancer",
      TRUE             ~ "ns"
    ),
    direction_plot = factor(
      direction_plot,
      levels = c("Depleted in cancer", "Enriched in cancer", "ns")
    ),
    feature_clean = clean_taxon_name(feature)
  )

# Per-cohort labels: top N depleted + top N enriched (matches top-hits
# bar panel directly, so species carry across figures). If a cohort has
# fewer than N in one direction, slice_head takes whatever exists.
volcano_labels <- volcano_df |>
  filter(qval < Q_THRESH, !is.na(coef_case)) |>
  group_by(cohort) |>
  group_modify(~ {
    rank_score <- abs(.x$coef_case) * .x$neglog10q
    ups   <- .x |> filter(coef_case > 0) |>
      arrange(desc(abs(coef_case) * neglog10q)) |>
      slice_head(n = VOLCANO_LABEL_N)
    downs <- .x |> filter(coef_case < 0) |>
      arrange(desc(abs(coef_case) * neglog10q)) |>
      slice_head(n = VOLCANO_LABEL_N)
    bind_rows(ups, downs)
  }) |>
  ungroup()

# Symmetric x-limits so +/- coefficients are visually balanced
x_max <- max(abs(volcano_df$coef_case), na.rm = TRUE)
x_max <- ceiling(x_max * 10) / 10
y_max <- max(volcano_df$neglog10q, na.rm = TRUE) * 1.05

build_volcano <- function(df, labels_df, cohort_name) {
  d  <- df     |> filter(cohort == cohort_name)
  lb <- labels_df |> filter(cohort == cohort_name)
  
  ggplot(d, aes(x = coef_case, y = neglog10q)) +
    geom_vline(xintercept = 0, linetype = "solid", colour = "grey70", linewidth = 0.4) +
    geom_hline(yintercept = -log10(Q_THRESH), linetype = "dashed",
               colour = "grey40", linewidth = 0.5) +
    geom_point(aes(colour = direction_plot), alpha = 0.75, size = 2) +
    geom_text_repel(
      data = lb,
      aes(label = feature_clean),
      size = 3,
      fontface = "italic",
      min.segment.length = 0,
      segment.size = 0.3,
      segment.colour = "grey60",
      max.overlaps = Inf,
      box.padding = 0.35,
      force = 2
    ) +
    scale_colour_manual(
      values = c(
        "Depleted in cancer" = "#2166ac",
        "Enriched in cancer" = "#d6604d",
        "ns"                 = "grey75"
      ),
      drop = FALSE
    ) +
    scale_x_continuous(limits = c(-x_max, x_max)) +
    scale_y_continuous(limits = c(0, y_max), expand = expansion(mult = c(0, 0.02))) +
    labs(
      title = cohort_name,
      x = "MaAsLin2 coefficient (case-oriented)",
      y = expression(-log[10](q))
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      legend.position = "none",
      panel.grid.minor = element_blank()
    )
}

cohorts_in_order <- levels(volcano_df$cohort)
plots <- lapply(cohorts_in_order, function(ch) build_volcano(volcano_df, volcano_labels, ch))

p_volcano <- wrap_plots(plots, ncol = 2) +
  plot_annotation(
    title = "Species-level differential abundance: MaAsLin2 across cohorts",
    subtitle = paste0(
      "Dashed line: q = ", Q_THRESH,
      ". Top ", VOLCANO_LABEL_N, " depleted + top ", VOLCANO_LABEL_N,
      " enriched per cohort labelled (by |coef| x -log10 q)."
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(size = 11, colour = "grey30")
    )
  ) &
  theme(plot.margin = margin(4, 6, 4, 6))

ggsave(
  file.path(out_dir, "07_volcano_panel_clean.png"),
  p_volcano, width = 13, height = 11, dpi = 800
)

# -------------------------------------------------------------------------
# 8. Beta diversity panel — 4 cohorts x 2 metrics (Bray + Jaccard)
# -------------------------------------------------------------------------
#
# Loads per-cohort beta_panel_<metric>_<cohort>.rds payloads saved by
# Script 1, then assembles into a 2-row x 4-column patchwork grid
# (rows = metrics, cols = cohorts). Each panel has:
#   - cohort sample counts in the strip
#   - PERMANOVA R2 / F / p and betadisper p in the subtitle
#   - 68%-level ellipses (tighter than default 95%) for cleaner separation
#   - axis labels showing % variance explained per PCoA axis
#
# Rationale for 68% ellipses: default stat_ellipse draws a 95% region
# which looks permissive and often visually contradicts significant
# PERMANOVA results in noisy microbiome data. 1-SD contours communicate
# group location more faithfully and are paired with the stats in the
# subtitle so readers see both.

beta_rds_paths <- expand.grid(
  cohort_key = names(cohort_labels),
  metric     = c("bray", "jaccard"),
  stringsAsFactors = FALSE
) |>
  mutate(
    path = map2_chr(cohort_key, metric, function(ck, m) {
      base <- switch(ck,
                     prjna813705      = "D:/MasterThesis/Vik/PRJNA813705/tables_metaphlan_prjna813705",
                     kulecka_leukemia = "D:/MasterThesis/Vik/Kulecka/Leukemia/tables_metaphlan_kulecka_leukemia",
                     kulecka_lymphoma = "D:/MasterThesis/Vik/Kulecka/Lymphoma/tables_metaphlan_kulecka_lymphoma",
                     cra              = "D:/MasterThesis/Vik/CRA007433/tables_metaphlan_cra"
      )
      file.path(base, sprintf("beta_panel_%s_%s.rds", m, ck))
    })
  )

missing_rds <- beta_rds_paths |> filter(!file.exists(path))
if (nrow(missing_rds) > 0) {
  warning("Beta panel RDS files missing for:\n",
          paste("  -", missing_rds$cohort_key, missing_rds$metric, collapse = "\n"),
          "\nRerun Script 1 on those datasets to regenerate. Skipping Figure 8.")
} else {
  
  metric_label <- c(bray = "Bray\u2013Curtis", jaccard = "Jaccard")
  
  # Returns just the numeric part; caller prepends "p" and the operator
  # is implicit in the returned string (e.g. "< 0.001" or "= 0.042")
  fmt_p <- function(p) {
    if (is.na(p)) return("= NA")
    if (p < 0.001) "< 0.001" else sprintf("= %.3f", p)
  }
  
  build_beta_panel <- function(rds_path, cohort_key, metric) {
    payload <- readRDS(rds_path)
    
    coords <- payload$coords
    coords$disease <- factor(coords$disease, levels = levels(droplevels(factor(coords$disease))))
    
    # Ensure Control always plots first (navy) and case second (orange)
    # for consistency across cohorts. Case label varies per cohort.
    lv <- levels(coords$disease)
    case_lv <- setdiff(lv, "Control")
    coords$disease <- factor(coords$disease, levels = c("Control", case_lv))
    
    # Re-map to a shared two-level factor so one global legend works
    # across all panels despite cohorts using different case labels.
    coords$Group <- factor(
      if_else(coords$disease == "Control", "Control", "Case"),
      levels = c("Control", "Case")
    )
    
    perm <- payload$permanova |> filter(term == "disease")
    r2   <- perm$R2[1]
    f    <- perm$`F`[1]
    p    <- perm[["Pr(>F)"]][1]
    bdp  <- payload$betadisper_p
    
    subtitle <- sprintf(
      "PERMANOVA: R\u00B2 = %.3f, F = %.2f, p %s   \u00B7   Betadisper p %s",
      r2, f, fmt_p(p), fmt_p(bdp)
    )
    
    rel_eig <- payload$rel_eig
    x_lab <- sprintf("PCo1 (%.1f%%)", rel_eig[1])
    y_lab <- sprintf("PCo2 (%.1f%%)", rel_eig[2])
    
    plot_title <- sprintf("%s - %s",
                          gsub("\n", " ", cohort_labels_with_n[cohort_key]),
                          metric_label[metric])
    
    p_pcoa <- ggplot(coords, aes(x = Axis1, y = Axis2, colour = Group)) +
      geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
      geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
      geom_point(size = 2.2, alpha = 0.8) +
      scale_colour_manual(
        values = c(Control = "#1b2a6b", Case = "#e6a23c")
      ) +
      guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
      labs(
        title    = plot_title,
        subtitle = subtitle,
        x = x_lab, y = y_lab, colour = NULL
      ) +
      theme(
        plot.title    = element_text(face = "bold", size = 10),
        plot.subtitle = element_text(size = 8.5, colour = "grey30"),
        axis.title    = element_text(size = 9),
        axis.text     = element_text(size = 8),
        panel.grid.minor = element_blank()
      )
    
    if (min(table(coords$Group)) >= 3) {
      p_pcoa <- p_pcoa +
        stat_ellipse(aes(fill = Group),
                     type = "norm", level = 0.68,
                     geom = "polygon", alpha = 0.12, linewidth = 0.6) +
        scale_fill_manual(
          values = c(Control = "#1b2a6b", Case = "#e6a23c"),
          guide = "none"
        )
    }
    
    p_pcoa
  }
  
  beta_plots <- pmap(beta_rds_paths, function(cohort_key, metric, path) {
    build_beta_panel(path, cohort_key, metric)
  })
  
  # plot_layout(guides = "collect") dedupes the per-panel legends into a
  # single one; the bottom "& theme(legend.position = 'bottom')" places it
  # below the whole grid.
  p_beta <- wrap_plots(beta_plots, ncol = 4, byrow = TRUE) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title    = "Beta diversity: PCoA across cohorts (MetaPhlAn species profiles)",
      subtitle = "Top row: Bray\u2013Curtis (abundance-weighted). Bottom row: Jaccard (presence/absence). Ellipses: 68% (1 SD) contours.",
      theme = theme(
        plot.title    = element_text(face = "bold", size = 15),
        plot.subtitle = element_text(size = 11, colour = "grey30")
      )
    ) &
    theme(
      plot.margin = margin(3, 5, 3, 5),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text  = element_text(size = 11)
    )
  
  ggsave(
    file.path(out_dir, "08_beta_diversity_panel_clean.png"),
    p_beta, width = 18, height = 9, dpi = 800
  )
}

cat("Clean lab meeting figures saved to:", out_dir, "\n")
