## =============================================================================
## humann_figures.R
## CRASH — Cancer-associated gut micRobiome Analysis of Strains in
##         Haematological patients
##
## Purpose:
##   KO heatmap with live KEGG enzyme name annotations. Reads per-disease-group
##   HUMAnN4 MaAsLin2 results (AML pooled, LN, NKTCL) and visualises the
##   top cross-cohort KO hits with human-readable enzyme names fetched via
##   KEGGREST and cached to disk to avoid redundant API calls.
##
## Input:
##   Combined/MaAsLin2/<GROUP>/community_ko/all_results.tsv
##     where GROUP ∈ {AML, LN, NKTCL}
##
## Outputs (written to CRASH_figs_humann/):
##   ko_heatmap_annotated.png   Top KOs (≥MIN_COHORTS groups, q<Q_THRESH)
##                              with KEGG enzyme names; rows capped at TOP_N_FEAT
##   kegg_name_cache.rds        Disk cache of fetched KEGG names
##
## Dependencies:
##   tidyverse, KEGGREST
##
## Notes:
##   - KEGG API rate-limited; fetch_one_ko() sleeps 0.35 s between calls
##   - Cache (kegg_name_cache.rds) persists across runs; delete to refresh
##   - HUMAnN4 KO tables normalised with `--special n` (UNGROUPED dropped
##     before RPK renormalisation); omitting this flag inflates all DA results
##   - AML group pools Belgian (PRJNA813705) + Polish (Kulecka) AML cohorts
## =============================================================================

library(tidyverse)
library(KEGGREST)

out_dir <- "D:/MasterThesis/Vik/CRASH_figs_humann"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

theme_set(theme_minimal(base_size = 13))

Q_THRESH    <- 0.25
MIN_COHORTS <- 2
TOP_N_FEAT  <- 30
MAX_LABEL   <- 45

MAASLIN_ROOT <- "D:/MasterThesis/Vik/Combined/MaAsLin2"

maaslin_paths <- list(
  AML   = file.path(MAASLIN_ROOT, "AML/community_ko/all_results.tsv"),
  LN    = file.path(MAASLIN_ROOT, "LN/community_ko/all_results.tsv"),
  NKTCL = file.path(MAASLIN_ROOT, "NKTCL/community_ko/all_results.tsv")
)

cohort_labels <- c(
  AML   = "AML\n(n=100 + 70 ctrl)",
  LN    = "LN\n(n=60 + 59 ctrl)",
  NKTCL = "NKTCL\n(n=42 + 33 ctrl)"
)

# ---- KEGG fetch + cache ---------------------------------------------

fetch_one_ko <- function(kid) {
  Sys.sleep(0.35)
  tryCatch({
    entry <- keggGet(kid)[[1]]
    name  <- entry[["NAME"]]
    if (is.null(name) || length(name) == 0 || is.na(name[1])) return(kid)
    name <- strsplit(as.character(name[1]), ";")[[1]][1]
    name <- sub("\\[EC:[^]]*\\]", "", name)
    name <- trimws(name)
    if (nchar(name) > MAX_LABEL) name <- paste0(substr(name, 1, MAX_LABEL - 3), "...")
    name
  }, error = function(e) {
    message("  Failed for ", kid, ": ", conditionMessage(e))
    kid   # fallback: return K-number
  })
}

fetch_ko_names <- function(ko_ids, cache_path) {
  if (file.exists(cache_path)) {
    cat("Loading cached KEGG names from", cache_path, "\n")
    cached <- readRDS(cache_path)
    # Fetch any IDs not yet in cache
    missing <- setdiff(ko_ids, names(cached))
    if (length(missing) == 0) return(cached)
    cat(sprintf("Fetching %d new KO names not in cache...\n", length(missing)))
    new_names <- setNames(vapply(missing, fetch_one_ko, character(1)), missing)
    result <- c(cached, new_names)
    saveRDS(result, cache_path)
    return(result)
  }
  
  cat(sprintf("Fetching KEGG names for %d KOs (first run only)...\n", length(ko_ids)))
  results <- setNames(vapply(ko_ids, function(k) {
    cat(sprintf("  %s\n", k))
    fetch_one_ko(k)
  }, character(1)), ko_ids)
  
  saveRDS(results, cache_path)
  cat("Saved cache to", cache_path, "\n")
  results
}

make_ko_label <- function(ko_id, name_lookup) {
  nm <- name_lookup[ko_id]
  if (is.na(nm) || nm == ko_id) return(ko_id)
  paste0(ko_id, ": ", nm)
}

# ---- Helpers --------------------------------------------------------

load_maaslin <- function(path, cohort) {
  if (!file.exists(path)) { warning("Missing: ", path); return(NULL) }
  read_tsv(path, show_col_types = FALSE) |>
    filter(metadata == "Disease") |>
    mutate(cohort = cohort, is_sig = qval < Q_THRESH)
}

make_coef_df <- function(all_dat, feature_list, name_lookup) {
  full_grid <- expand_grid(feature = feature_list, cohort = names(cohort_labels))
  dat <- all_dat |>
    filter(feature %in% feature_list) |>
    select(feature, cohort, coef, qval, is_sig) |>
    right_join(full_grid, by = c("feature", "cohort"))
  
  rank_df <- dat |>
    group_by(feature) |>
    summarise(n_sig         = sum(is_sig, na.rm = TRUE),
              min_q         = suppressWarnings(min(qval, na.rm = TRUE)),
              mean_abs_coef = mean(abs(coef), na.rm = TRUE),
              .groups = "drop") |>
    mutate(min_q = if_else(is.infinite(min_q), NA_real_, min_q)) |>
    arrange(desc(n_sig), min_q, desc(mean_abs_coef))
  
  feat_order <- rank_df$feature
  max_abs    <- max(abs(dat$coef), na.rm = TRUE)
  
  label_map <- setNames(
    vapply(feat_order, make_ko_label, character(1), name_lookup = name_lookup),
    feat_order
  )
  
  dat |>
    mutate(
      cohort_label  = factor(cohort_labels[cohort], levels = unname(cohort_labels)),
      display_label = label_map[feature],
      display_label = factor(display_label, levels = rev(unname(label_map))),
      is_sig_safe   = coalesce(is_sig, FALSE),
      sig_marker    = if_else(is_sig_safe & !is.na(coef), "*", ""),
      cell_label    = if_else(is.na(coef), "--",
                              paste0(sprintf("%.2f", coef), sig_marker)),
      label_col     = if_else(is.na(coef) | abs(coef) < 0.6 * max_abs,
                              "black", "white"),
      max_abs       = max_abs
    )
}

draw_heatmap <- function(df, title_text, subtitle_text,
                         out_file, fig_width = 13, fig_height = 8) {
  max_abs <- unique(df$max_abs)
  p <- ggplot(df, aes(x = cohort_label, y = display_label, fill = coef)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = cell_label, colour = label_col),
              size = 3.4, show.legend = FALSE) +
    scale_colour_identity() +
    scale_fill_gradient2(
      low = "#2166ac", mid = "white", high = "#d6604d",
      midpoint = 0, na.value = "grey88",
      limits = c(-max_abs, max_abs),
      name = "MaAsLin2\ncoefficient\n(disease vs ctrl)"
    ) +
    labs(title = title_text, subtitle = subtitle_text, x = NULL, y = NULL) +
    theme(
      axis.text.y     = element_text(size = 8.5),
      axis.text.x     = element_text(size = 11, lineheight = 0.9),
      panel.grid      = element_blank(),
      legend.position = "right",
      plot.title      = element_text(face = "bold"),
      plot.subtitle   = element_text(size = 10, colour = "grey35")
    )
  ggsave(out_file, p, width = fig_width, height = fig_height, dpi = 300)
  message("Saved: ", basename(out_file))
  invisible(p)
}

# ---- Load KO results ------------------------------------------------

cat("Loading KO MaAsLin2 results...\n")
ko_all <- imap_dfr(maaslin_paths, load_maaslin)

top_ko_features <- ko_all |>
  group_by(feature) |>
  summarise(n_sig = sum(is_sig, na.rm = TRUE),
            min_q = suppressWarnings(min(qval, na.rm = TRUE)),
            .groups = "drop") |>
  filter(n_sig >= MIN_COHORTS) |>
  arrange(desc(n_sig), min_q) |>
  slice_head(n = TOP_N_FEAT) |>
  pull(feature)

# Strip any description suffix ("K00194: something" -> "K00194")
ko_ids <- sub(":.*$", "", top_ko_features)

cat(sprintf("%d KOs in >=%d cohorts\n", length(ko_ids), MIN_COHORTS))

# ---- Fetch / load cached KEGG names ---------------------------------

ko_names <- fetch_ko_names(
  ko_ids     = ko_ids,
  cache_path = file.path(out_dir, "ko_name_cache.rds")
)

cat("\nKO annotations:\n")
for (k in names(ko_names)) cat(sprintf("  %s: %s\n", k, ko_names[k]))

# ---- Draw heatmap ---------------------------------------------------

df <- make_coef_df(ko_all, top_ko_features, name_lookup = ko_names)

draw_heatmap(
  df,
  title_text    = sprintf("Top KOs altered in >=2 cohorts"),
  subtitle_text = "Rows ranked by n cohorts significant then min q. * marks q < 0.25.",
  out_file      = file.path(out_dir, "04_ko_top_heatmap_annotated.png"),
  fig_height    = max(5, 0.38 * length(top_ko_features) + 2)
)

cat("Done.\n")
