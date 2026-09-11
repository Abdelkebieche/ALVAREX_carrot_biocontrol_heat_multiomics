############################################################
## 00_helpers_WGCNA_v2.R — Fonctions partagées
## Sourcé par tous les scripts 01-07
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
clean_sample_names <- function(x){
  x <- gsub("^X", "", x)
  m <- regmatches(x, regexpr("\\d{4}", x))
  ifelse(m == "" | is.na(m), x, m)
}

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

# ── featureCounts reader ──────────────────────────────────
read_featurecounts_matrix <- function(path){
  fc <- utils::read.delim(path, header = TRUE, sep = "\t",
                          comment.char = "#", quote = "", check.names = FALSE,
                          stringsAsFactors = FALSE)
  stopifnot("Geneid" %in% colnames(fc))
  drop <- intersect(colnames(fc), c("Chr","Start","End","Strand","Length"))
  mat <- fc |>
    dplyr::select(-dplyr::all_of(drop)) |>
    tibble::column_to_rownames("Geneid") |>
    as.matrix()
  mode(mat) <- "numeric"
  colnames(mat) <- clean_sample_names(colnames(mat))
  colnames(mat) <- make.unique(colnames(mat))
  # extract 4-digit sample ID
  colnames(mat) <- sapply(colnames(mat), function(x){
    m <- regmatches(x, regexpr("[0-9]{4}", x))
    if (length(m) == 1 && nchar(m) == 4) m else x
  })
  mat
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
time_display <- c("P2" = "D2", "P3" = "D4")

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

cat("✓ Helpers loaded\n")