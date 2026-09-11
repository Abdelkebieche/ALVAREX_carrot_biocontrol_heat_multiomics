############################################################
## 00_helpers_WGCNA_LCMS.R — Fonctions partagées (LC-MS)
## Sourcé par tous les scripts 01-05
# Repository configuration (portable; no HPC-specific paths)
if (!exists('CARROT_CONFIG_LOADED')) source('../../../config/config.R')


############################################################

# ── Packages ──────────────────────────────────────────────
load_pkgs <- function(pkgs, install_missing = FALSE, bioc = FALSE){
  for (p in pkgs){
    if (!requireNamespace(p, quietly = TRUE)){
      if (!install_missing) stop("Package manquant: ", p)
      if (bioc) BiocManager::install(p, ask = FALSE, update = FALSE)
      else install.packages(p, dependencies = TRUE)
    }
    suppressPackageStartupMessages(library(p, character.only = TRUE))
  }
}

# ── Metadata normalisation ────────────────────────────────
normalize_temperature <- function(x){
  x <- toupper(trimws(as.character(x)))
  x <- gsub("\\s+|°", "", x)
  x <- gsub("TEMP", "T", x)
  x <- gsub("^0$", "T0", x); x <- gsub("^1$", "T1", x); x <- gsub("^2$", "T2", x)
  x
}

normalize_treatment <- function(x){
  x <- toupper(trimws(as.character(x)))
  x <- gsub("\\s+", "", x)
  dplyr::case_when(
    x %in% c("EAU","WATER","H2O","CTRL","CONTROL") ~ "Control",
    x %in% c("SDP","TREATED","PRI") ~ "Treated",
    TRUE ~ x
  )
}

# ── LC-MS data reader ─────────────────────────────────────
#' Lit la matrice LC-MS et extrait les IDs d'échantillons (4 chiffres)
#' @param path Chemin vers le fichier CSV/TSV LC-MS
#' @param sep Séparateur (détecté automatiquement si NULL)
#' @return Liste: mat (metabolites × samples), feature_info (metadata des features)
read_lcms_matrix <- function(path, sep = NULL) {
  # Détection automatique du séparateur
  first_line <- readLines(path, n = 1)
  if (is.null(sep)) {
    # Compter les occurrences de chaque séparateur potentiel
    n_tab   <- nchar(gsub("[^\t]", "", first_line))
    n_semi  <- nchar(gsub("[^;]", "", first_line))
    n_comma <- nchar(gsub("[^,]", "", first_line))
    sep <- c("\t", ";", ",")[which.max(c(n_tab, n_semi, n_comma))]
    cat("  Séparateur détecté:", ifelse(sep == "\t", "TAB", sep), "\n")
  }
  
  raw <- utils::read.delim(path, header = TRUE, sep = sep,
                           check.names = FALSE, stringsAsFactors = FALSE,
                           quote = "\"")
  
  cat("  Colonnes brutes:", ncol(raw), "\n")
  cat("  Lignes brutes:", nrow(raw), "\n")
  
  all_cols <- colnames(raw)
  
  # ── Identifier les colonnes d'échantillons ──
  # Critère: contient ".raw" dans le nom (ex: "2025-F-0533_013.raw (F13)")
  # Exclure explicitement les colonnes Log2 Fold Change et Adj. P-value
  sample_cols <- grep("\\.raw\\b", all_cols, value = TRUE)
  
  # ── Colonnes à exclure (statistiques pré-calculées) ──
  exclude_pattern <- "^(Log2 Fold Change|Adj\\. P-value|P-value|Fold Change|T-test|ANOVA)"
  exclude_cols <- grep(exclude_pattern, all_cols, value = TRUE)
  
  # ── Colonnes info = tout le reste sauf samples et exclues ──
  info_cols <- setdiff(all_cols, c(sample_cols, exclude_cols))
  
  cat("  Colonnes échantillons (.raw):", length(sample_cols), "\n")
  cat("  Colonnes info:", length(info_cols), "\n")
  cat("  Colonnes statistiques exclues:", length(exclude_cols), "\n")
  
  stopifnot(length(sample_cols) > 0)
  
  # Extraire la matrice d'intensités
  mat <- as.matrix(raw[, sample_cols, drop = FALSE])
  mode(mat) <- "numeric"
  
  # Créer un ID de feature
  if ("ID_MT" %in% info_cols && "ID_MT" %in% colnames(raw)) {
    feat_ids <- as.character(raw[["ID_MT"]])
  } else if ("m/z" %in% colnames(raw) && "RT [min]" %in% colnames(raw)) {
    feat_ids <- paste0("M", round(as.numeric(raw[["m/z"]]), 4), "T",
                       round(as.numeric(raw[["RT [min]"]]), 3))
  } else if ("m/z" %in% colnames(raw) && "RT[min]" %in% colnames(raw)) {
    feat_ids <- paste0("M", round(as.numeric(raw[["m/z"]]), 4), "T",
                       round(as.numeric(raw[["RT[min]"]]), 3))
  } else {
    feat_ids <- paste0("Feature_", seq_len(nrow(raw)))
  }
  feat_ids <- make.unique(as.character(feat_ids))
  rownames(mat) <- feat_ids
  
  # ── Extraire les IDs d'échantillons (4 chiffres après le 2ème "F-") ──
  # Format: "2025-F-0533_013.raw (F13)"
  #   - "2025" = année (ignorer)
  #   - "0533" = ID échantillon (ce qu'on veut)
  #   - "013"  = numéro de run (ignorer)
  new_colnames <- sapply(sample_cols, function(x) {
    # Extraire tous les groupes de 4 chiffres après "F-"
    matches <- gregexpr("(?<=F-)[0-9]{4}", x, perl = TRUE)
    all_matches <- regmatches(x, matches)[[1]]
    # Le premier match après F- sera "0533" (pas "2025" car c'est avant F-)
    if (length(all_matches) >= 1) return(all_matches[1])
    # Fallback
    m <- regmatches(x, regexpr("[0-9]{4}", x))
    if (length(m) == 1) return(m)
    return(x)
  })
  colnames(mat) <- as.character(new_colnames)
  
  # Vérifier les doublons
  if (any(duplicated(new_colnames))) {
    cat("  ⚠ WARNING: IDs dupliqués détectés!\n")
    cat("  Doublons:", paste(new_colnames[duplicated(new_colnames)], collapse=", "), "\n")
    colnames(mat) <- make.unique(as.character(new_colnames))
  }
  
  # Diagnostic: afficher quelques mappings
  cat("  Mapping colonnes → IDs (premiers 5):\n")
  for (i in seq_len(min(5, length(sample_cols)))) {
    cat("    ", sample_cols[i], " → ", new_colnames[i], "\n")
  }
  
  # Feature info
  feature_info <- raw[, info_cols, drop = FALSE]
  feature_info$FeatureID <- feat_ids
  
  cat("  Matrice finale:", nrow(mat), "features x", ncol(mat), "échantillons\n")
  
  list(mat = mat, feature_info = feature_info)
}

# ── Theme publication (sized for journal figures) ─────────
theme_pub <- function(base_size = 14){
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "grey94", colour = NA),
      strip.text       = ggplot2::element_text(face = "bold", size = base_size),
      axis.title       = ggplot2::element_text(face = "bold", size = base_size),
      axis.text        = ggplot2::element_text(colour = "black", size = base_size - 2),
      plot.title       = ggplot2::element_text(face = "bold", size = base_size + 3),
      plot.subtitle    = ggplot2::element_text(colour = "grey35", size = base_size - 1),
      legend.title     = ggplot2::element_text(face = "bold", size = base_size - 1),
      legend.text      = ggplot2::element_text(size = base_size - 2),
      legend.key       = ggplot2::element_blank(),
      legend.position  = "top",
      panel.spacing    = ggplot2::unit(0.6, "lines")
    )
}

# ── Display label mapping (figures only) ──────────────────
temp_display <- c("T0" = "NH", "T1" = "HS-7", "T2" = "HS-2")
time_display <- c("P3" = "48h post-trt", "P4" = "10d post-trt")

# Relabel a vector of internal codes to display labels
relabel_temp <- function(x) { r <- temp_display[as.character(x)]; ifelse(is.na(r), as.character(x), r) }
relabel_time <- function(x) { r <- time_display[as.character(x)]; ifelse(is.na(r), as.character(x), r) }

# ── Save helpers (PDF + PNG + SVG) ────────────────────────
save_plot <- function(p, path_base, w = 7, h = 5, dpi = 350){
  ggplot2::ggsave(paste0(path_base, ".pdf"), p, width = w, height = h, bg = "white")
  ggplot2::ggsave(paste0(path_base, ".png"), p, width = w, height = h, dpi = dpi, bg = "white")
  ggplot2::ggsave(paste0(path_base, ".svg"), p, width = w, height = h, bg = "white",
                  device = svglite::svglite)
}

save_base_plot <- function(fn, path_base, w = 7, h = 5){
  pdf(paste0(path_base, ".pdf"), width = w, height = h, useDingbats = FALSE); fn(); dev.off()
  png(paste0(path_base, ".png"), width = w, height = h, units = "in", res = 300, bg = "white"); fn(); dev.off()
  svglite::svglite(paste0(path_base, ".svg"), width = w, height = h, bg = "white"); fn(); dev.off()
}

# ── Starify p-values ─────────────────────────────────────
starify <- function(p){
  ifelse(is.na(p), "",
         ifelse(p < 0.001, "***",
                ifelse(p < 0.01, "**",
                       ifelse(p < 0.05, "*",
                              ifelse(p < 0.10, ".", "")))))
}

# ── PCA centroïdes ────────────────────────────────────────
compute_centroids <- function(pca_df, pcx = "PC1", pcy = "PC2",
                              group_vars, trajectory_var = NULL,
                              trajectory_levels = NULL){
  cent <- pca_df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) |>
    dplyr::summarise(
      cx = mean(.data[[pcx]], na.rm = TRUE),
      cy = mean(.data[[pcy]], na.rm = TRUE),
      n = dplyr::n(), .groups = "drop"
    )
  
  seg <- NULL
  if (!is.null(trajectory_var) && trajectory_var %in% group_vars) {
    if (!is.null(trajectory_levels))
      cent[[trajectory_var]] <- factor(cent[[trajectory_var]], levels = trajectory_levels)
    
    seg_group <- setdiff(group_vars, trajectory_var)
    if (length(seg_group) > 0) {
      seg <- cent |>
        dplyr::group_by(dplyr::across(dplyr::all_of(seg_group))) |>
        dplyr::arrange(.data[[trajectory_var]]) |>
        dplyr::mutate(cxend = dplyr::lead(cx), cyend = dplyr::lead(cy)) |>
        dplyr::ungroup() |>
        dplyr::filter(!is.na(cxend))
    } else {
      seg <- cent |>
        dplyr::arrange(.data[[trajectory_var]]) |>
        dplyr::mutate(cxend = dplyr::lead(cx), cyend = dplyr::lead(cy)) |>
        dplyr::filter(!is.na(cxend))
    }
  }
  list(centroids = cent, segments = seg)
}

# ── Find latest run dir ──────────────────────────────────
find_latest_dir <- function(base, pattern){
  dirs <- list.files(base, full.names = TRUE, pattern = pattern)
  dirs <- dirs[file.info(dirs)$isdir]
  if (length(dirs) == 0) stop("No run dir found: ", base, " pattern=", pattern)
  dirs[order(file.info(dirs)$mtime, decreasing = TRUE)][1]
}

# ── Logger ────────────────────────────────────────────────
init_logger <- function(logfile){
  dir.create(dirname(logfile), recursive = TRUE, showWarnings = FALSE)
  function(...){
    msg <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", paste0(..., collapse = ""))
    cat(msg, "\n")
    cat(msg, "\n", file = logfile, append = TRUE)
  }
}

cat("✓ Helpers LC-MS loaded\n")