# NOTE: companion annotation-builder reference implementation.
# Paths have been sanitized. See docs/05_functional_dictionaries.md before running.

############################################################
## 09B_TF_panorama_complete_REVISED.R
## Figure 9B — Complete heat-responsive transcription-factor panorama
##
## SCIENTIFIC PRINCIPLE
##   1. TF identity/family is taken from the genome annotation / TF database,
##      not reconstructed from loose function keywords.
##   2. ALL annotated TF families are eligible. No "focus-family" filter is
##      used to define the panorama.
##   3. A TF is retained for the main heatmap when it is Up or Down in >=1 of
##      the FOUR intended intra-genotype heat-vs-NH contrasts only.
##   4. Family/subfamily labels are kept conservative. DREB is treated as an
##      ERF subgroup rather than as an independent top-level family when this
##      is explicit in the source annotation.
##   5. The main panorama contains ALL DE TFs (no top-N truncation).
##   6. Family interpretation is based on both raw counts and the proportion
##      of the expressed family that is DE; one-sided Fisher tests provide a
##      family-overrepresentation statistic with BH correction.
##
## FIXED DISPLAY ORDER REQUESTED
##   ROBILA NH -> ROBILA HS-7 -> ROBILA HS-2 ->
##   PRESTO NH -> PRESTO HS-7 -> PRESTO HS-2
##
## CONTRAST ORDER
##   ROBILA HS-7 vs NH
##   ROBILA HS-2 vs NH
##   PRESTO HS-7 vs NH
##   PRESTO HS-2 vs NH
##
## OUTPUTS
##   00_tables/
##     - complete TF census
##     - family enrichment
##     - per-contrast family counts and fractions
##     - complete selected DE-TF table
##     - temporal/condition-response patterns
##     - family completeness audit
##     - conservative curation table
##   01_panorama/
##     - COMPLETE DE-TF heatmap (all selected genes)
##     - optional top-variable all-expressed-TF context heatmap
##     - family count heatmap
##   02_family_heatmaps/
##     - one DE-only heatmap for EVERY family represented among DE TFs
##   03_QC/
##     - metadata, design-independent inputs, ID mapping/QC, RDS objects,
##       sessionInfo
##
## LITERATURE / DATABASE BASIS
##   - PlantRegMap / PlantTFDB v5 family assignment uses DNA-binding,
##     auxiliary and forbidden domains.
##     Tian et al. 2020, Nucleic Acids Research 48:D1104-D1113.
##     DOI: 10.1093/nar/gkz1020
##   - Current PlantTFDB Daucus carota reference repertoire:
##     1906 TF loci in 56 families (ARS-USDA v2.0 reference).
##   - ERF/DREB:
##     Nakano et al. 2006, Plant Physiology 140:411-432.
##     DOI: 10.1104/pp.105.073783
##     DREB proteins are sequence-defined AP2/ERF-domain factors and are
##     treated here as an ERF subgroup, not an independent top-level family.
##   - Heat-response transcriptomes routinely recover many TF families,
##     not HSF alone; therefore the panorama is annotation-driven rather than
##     limited to canonical heat-stress families.
############################################################

##############################
## 0) Configuration
##############################
source("00_builder_config.R")
source("00_builder_helpers.R")

# Same genome annotation used in the earlier TF panorama script.
TF_ANNOT_FILE <- DH13_ANNOT_FILE

# Four intended heat-vs-NH contrasts ONLY.
TF_CONTRAST_ORDER <- c(
  "ROBILA_T1_vs_T0",
  "ROBILA_T2_vs_T0",
  "PRESTO_T1_vs_T0",
  "PRESTO_T2_vs_T0"
)
TF_CONTRAST_LABELS <- c(
  "ROBILA HS-7",
  "ROBILA HS-2",
  "PRESTO HS-7",
  "PRESTO HS-2"
)

# Fixed expression display order requested by the user.
TF_GROUP_ORDER <- c(
  "ROBILA_T0", "ROBILA_T1", "ROBILA_T2",
  "PRESTO_T0", "PRESTO_T1", "PRESTO_T2"
)
TF_GROUP_LABELS <- c(
  "ROBILA NH", "ROBILA HS-7", "ROBILA HS-2",
  "PRESTO NH", "PRESTO HS-7", "PRESTO HS-2"
)

# DE definition. If the upstream DEG tables contain a Status column, it is
# deliberately reconstructed here from FDR/logFC so the selection criterion is
# explicit and reproducible in this script.
TF_FDR_CUTOFF <- PADJ_CUTOFF
TF_LFC_CUTOFF <- lfc_threshold

# Figure settings.
TF_Z_LIMIT              <- 2.5
TF_CLUSTER_ROWS          <- TRUE
TF_ROW_MM                <- 3.3
TF_FIG_WIDTH_MM          <- 235
TF_RASTER_DPI            <- 450
TF_WRITE_PDF             <- TRUE
TF_WRITE_SVG             <- TRUE
TF_WRITE_PNG             <- TRUE
TF_WRITE_TIFF            <- TRUE
TF_RASTER_MAX_HEIGHT_MM  <- 650  # skip giant raster; PDF/SVG remain complete

# ------------------------------------------------------------------
# GLOBAL COMPACT TF EXPRESSION PANORAMAS
# ------------------------------------------------------------------
# These figures are deliberately simple:
#   - expression only
#   - no TF-family split
#   - no Up/Down/NS panel
#   - no gene names
#   - hierarchical clustering of TF rows
#   - fixed biological column order:
#       ROBILA NH -> HS-7 -> HS-2 -> PRESTO NH -> HS-7 -> HS-2
#
# 1) ALL EXPRESSED TFs:
#    gives the global transcription-factor expression landscape.
# 2) ALL DE TFs:
#    same compact representation, restricted to TFs DE in >=1 of the
#    four intended heat-vs-NH contrasts.
#
# IMPORTANT:
# Row-wise z-scores are used only for pattern visualisation.
# They do not represent differential-expression significance.

TF_GLOBAL_Z_LIMIT                 <- 2
TF_GLOBAL_CLUSTER_ROWS            <- TRUE
TF_GLOBAL_SHOW_ROW_NAMES          <- FALSE
TF_GLOBAL_SHOW_ROW_DENDROGRAM     <- TRUE
TF_GLOBAL_WIDTH_MM                <- 185
TF_GLOBAL_ALL_HEIGHT_MM           <- 190
TF_GLOBAL_DE_HEIGHT_MM            <- 175
TF_GLOBAL_RASTER_QUALITY          <- 3

# Optional contextual heatmap: top-variable ALL EXPRESSED TFs.
# This does NOT affect selection for the complete DE-TF panorama.
TF_ALL_EXPRESSED_TOP_N <- 150L

# Family order in the complete panorama:
# "n_DE" = most heat-responsive families first; "alphabetical" = alphabetical.
TF_FAMILY_ORDER_MODE <- "n_DE"

# Keep source-annotated TFs by default. Suspicious/ambiguous cases are FLAGGED,
# not silently removed. Manual overrides can exclude a gene if required.
TF_EXCLUDE_AMBIGUOUS_BY_DEFAULT <- FALSE

# Current PlantTFDB v5 Daucus carota family repertoire (56 families).
# Used ONLY as an audit/reference list; it does NOT filter the DH13M14 analysis.
PLANTTFDB_DCA_FAMILIES <- c(
  "AP2", "ARF", "ARR-B", "B3", "BBR-BPC", "BES1",
  "C2H2", "C3H", "CAMTA", "CO-like", "CPP", "DBB", "Dof",
  "E2F/DP", "ERF", "FAR1", "G2-like", "GATA", "GRAS",
  "GRF", "GeBP", "HB-PHD", "HB-other", "HD-ZIP", "HRT-like",
  "HSF", "LBD", "LFY", "LSD", "M-type_MADS", "MIKC_MADS",
  "MYB", "MYB_related", "NAC", "NF-X1", "NF-YA", "NF-YB",
  "NF-YC", "Nin-like", "RAV", "S1Fa-like", "SAP", "SBP", "SRS",
  "STAT", "TALE", "TCP", "Trihelix", "VOZ", "WOX", "WRKY",
  "Whirly", "YABBY", "ZF-HD", "bHLH", "bZIP"
)

set.seed(SEED)
options(stringsAsFactors = FALSE)

##############################
## 1) Packages
##############################
load_pkgs(
  c("edgeR", "ComplexHeatmap"),
  install_missing = INSTALL_MISSING_PKGS,
  bioc = TRUE
)
load_pkgs(
  c("dplyr", "tidyr", "tibble", "readr", "stringr",
    "circlize", "svglite"),
  install_missing = INSTALL_MISSING_PKGS,
  bioc = FALSE
)

##############################
## 2) Utility functions
##############################
`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || all(is.na(x))) y else x
}

assert_columns <- function(x, required, object_name) {
  missing_cols <- setdiff(required, colnames(x))
  if (length(missing_cols) > 0L) {
    stop(
      object_name, " is missing required column(s): ",
      paste(missing_cols, collapse = ", ")
    )
  }
}

read_csv_flexible <- function(path) {
  if (!file.exists(path)) stop("Input file not found: ", path)
  x <- tryCatch(
    readr::read_csv(path, show_col_types = FALSE, progress = FALSE),
    error = function(e) NULL
  )
  if (is.null(x) || ncol(x) <= 1L) {
    x <- readr::read_delim(
      path,
      delim = ";",
      show_col_types = FALSE,
      progress = FALSE
    )
  }
  as.data.frame(x, stringsAsFactors = FALSE)
}

zscore_rows_safe <- function(mat) {
  mat <- as.matrix(mat)
  mu <- rowMeans(mat, na.rm = TRUE)
  sds <- apply(mat, 1L, stats::sd, na.rm = TRUE)
  out <- sweep(mat, 1L, mu, "-")
  valid <- is.finite(sds) & sds > 0
  out[valid, ] <- sweep(
    out[valid, , drop = FALSE],
    1L, sds[valid], "/"
  )
  out[!valid, ] <- 0
  out[!is.finite(out)] <- 0
  out
}

safe_file_stub <- function(x) {
  x <- iconv(x, to = "ASCII//TRANSLIT")
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x
}

clean_label <- function(x, max_chars = 68L) {
  x <- dplyr::coalesce(as.character(x), "")
  x <- trimws(gsub("\\s+", " ", x))
  x <- gsub("[.;]+$", "", x)
  too_long <- nchar(x) > max_chars
  x[too_long] <- paste0(substr(x[too_long], 1L, max_chars - 3L), "...")
  x
}

parse_logical_override <- function(x) {
  y <- trimws(tolower(as.character(x)))
  out <- rep(NA, length(y))
  out[y %in% c("true", "t", "1", "yes", "y", "oui")] <- TRUE
  out[y %in% c("false", "f", "0", "no", "n", "non")] <- FALSE
  out
}

status_from_de <- function(logFC, FDR) {
  dplyr::case_when(
    is.finite(FDR) &
      FDR <= TF_FDR_CUTOFF &
      logFC >= TF_LFC_CUTOFF ~ "Up",
    is.finite(FDR) &
      FDR <= TF_FDR_CUTOFF &
      logFC <= -TF_LFC_CUTOFF ~ "Down",
    TRUE ~ "NS"
  )
}

# Conservative harmonisation only.
# Source family is always retained separately in TF_family_source.
canonicalise_tf_family <- function(family, annotation = "") {
  f <- trimws(as.character(family))
  a <- stringr::str_to_lower(dplyr::coalesce(as.character(annotation), ""))
  
  f[is.na(f) | !nzchar(f)] <- "Unclassified_TF"
  
  # Safe spelling aliases.
  f[f %in% c("MYB-related", "MYB_related")] <- "MYB_related"
  f[f %in% c("M-type MADS", "M-type_MADS")] <- "M-type_MADS"
  f[f %in% c("MIKC MADS", "MIKC_MADS")] <- "MIKC_MADS"
  
  # DREB is an ERF/AP2-domain subgroup rather than an independent
  # PlantTFDB top-level family.
  f[f == "DREB"] <- "ERF"
  f[f %in% c("AP2/ERF-ERF", "AP2_ERF_ERF")] <- "ERF"
  f[f %in% c("AP2/ERF-AP2", "AP2_ERF_AP2")] <- "AP2"
  f[f %in% c("AP2/ERF-RAV", "AP2_ERF_RAV")] <- "RAV"
  
  # Generic AP2/ERF source labels are resolved only when the functional
  # annotation is explicit. Otherwise they remain visibly unresolved.
  generic <- f %in% c("AP2/ERF", "AP2_ERF")
  is_dreb <- stringr::str_detect(
    a,
    "\\bdreb[0-9a-z-]*\\b|dehydration[- ]responsive element[- ]binding"
  )
  is_erf <- stringr::str_detect(
    a,
    "ethylene[- ]responsive factor|ethylene response factor|\\berf[0-9a-z-]*\\b"
  )
  is_ap2 <- stringr::str_detect(
    a,
    "apetala2|\\bap2[- ]like\\b|ap2 transcription factor"
  )
  is_rav <- stringr::str_detect(a, "\\brav[0-9a-z-]*\\b")
  
  f[generic & is_dreb] <- "ERF"
  f[generic & !is_dreb & is_erf] <- "ERF"
  f[generic & !is_dreb & !is_erf & is_ap2] <- "AP2"
  f[generic & !is_dreb & !is_erf & !is_ap2 & is_rav] <- "RAV"
  f[generic & !(is_dreb | is_erf | is_ap2 | is_rav)] <- "AP2_ERF_unresolved"
  
  f
}

infer_tf_subgroup_conservative <- function(family_source, family_plot, annotation) {
  fsrc <- stringr::str_to_lower(dplyr::coalesce(as.character(family_source), ""))
  f <- stringr::str_to_lower(dplyr::coalesce(as.character(family_plot), ""))
  a <- stringr::str_to_lower(dplyr::coalesce(as.character(annotation), ""))
  
  dplyr::case_when(
    # ERF/AP2 superfamily — retain DREB only when explicit.
    f == "erf" &
      (fsrc == "dreb" |
         stringr::str_detect(
           a,
           "\\bdreb[0-9a-z-]*\\b|dehydration[- ]responsive element[- ]binding"
         )) ~ "DREB-like (ERF subgroup)",
    f == "erf" & stringr::str_detect(
      a,
      "ethylene[- ]responsive factor|ethylene response factor"
    ) ~ "ERF-like",
    f == "ap2" & stringr::str_detect(a, "apetala2|\\bap2\\b") ~ "AP2-like",
    f == "rav" & stringr::str_detect(a, "\\brav[0-9a-z-]*\\b") ~ "RAV-like",
    
    # HSF — A/B/C class only when explicitly named.
    f == "hsf" & stringr::str_detect(
      a,
      "hsfa|heat stress transcription factor a"
    ) ~ "HSFA-like",
    f == "hsf" & stringr::str_detect(
      a,
      "hsfb|heat stress transcription factor b"
    ) ~ "HSFB-like",
    f == "hsf" & stringr::str_detect(
      a,
      "hsfc|heat stress transcription factor c"
    ) ~ "HSFC-like",
    
    # MYB repeat architecture only when explicit.
    f %in% c("myb", "myb_related") &
      stringr::str_detect(a, "r2r3|r2-r3") ~ "R2R3-MYB",
    f %in% c("myb", "myb_related") &
      stringr::str_detect(a, "3r-myb|r1r2r3|myb3r") ~ "3R-MYB",
    f %in% c("myb", "myb_related") &
      stringr::str_detect(a, "1r-myb|shaqkyf") ~ "1R-MYB / SHAQKYF-like",
    
    # TCP class only if explicit.
    f == "tcp" & stringr::str_detect(a, "class i|\\bpcf\\b") ~ "TCP class I / PCF-like",
    f == "tcp" & stringr::str_detect(
      a,
      "class ii|\\bcin\\b|cycloidea|\\bcyc\\b|\\btb1\\b|branched1"
    ) ~ "TCP class II / CIN-CYC-TB1-like",
    
    # GRAS named clades only when explicit.
    f == "gras" & stringr::str_detect(
      a,
      "\\bdella\\b|protein gai|\\bgai\\b|\\brga\\b|\\brgl[0-9]*\\b"
    ) ~ "DELLA / GAI-like GRAS",
    f == "gras" & stringr::str_detect(a, "scarecrow|\\bscr\\b") ~ "SCARECROW-like GRAS",
    f == "gras" & stringr::str_detect(a, "short-root|\\bshr\\b") ~ "SHORT-ROOT-like GRAS",
    
    # Named orthologue-like annotations, without pretending to know
    # phylogenetic subgroup.
    f == "wrky" & stringr::str_detect(a, "wrky[0-9]+") ~ "Named WRKY-like",
    f == "nac" & stringr::str_detect(a, "nac[0-9]+") ~ "Named NAC-like",
    f == "bhlh" & stringr::str_detect(a, "bhlh[0-9]+|\\bmyc[0-9]*\\b") ~ "Named bHLH/MYC-like",
    f == "bzip" & stringr::str_detect(a, "bzip[0-9]+|\\bgbf[0-9]*\\b") ~ "Named bZIP-like",
    
    TRUE ~ NA_character_
  )
}

tf_annotation_confidence <- function(family_source, family_plot, annotation) {
  fs <- trimws(dplyr::coalesce(as.character(family_source), ""))
  fp <- trimws(dplyr::coalesce(as.character(family_plot), ""))
  a <- trimws(dplyr::coalesce(as.character(annotation), ""))
  
  dplyr::case_when(
    fp == "Unclassified_TF" ~ "Low",
    fp == "AP2_ERF_unresolved" ~ "Medium",
    nzchar(fs) ~ "High_family / subfamily_conservative",
    nzchar(a) ~ "Medium",
    TRUE ~ "Low"
  )
}

response_pattern_two <- function(status_hs7, status_hs2) {
  dplyr::case_when(
    status_hs7 == "NS" & status_hs2 == "NS" ~ "None",
    status_hs7 != "NS" & status_hs2 == "NS" ~ "HS-7 only",
    status_hs7 == "NS" & status_hs2 != "NS" ~ "HS-2 only",
    status_hs7 == status_hs2 & status_hs7 != "NS" ~ paste0("Both ", status_hs7),
    status_hs7 != "NS" & status_hs2 != "NS" &
      status_hs7 != status_hs2 ~ "Opposite directions",
    TRUE ~ "Other"
  )
}

##############################
## 3) Heatmap helpers
##############################
make_top_annotation <- function() {
  genotype <- c(rep("ROBILA", 3), rep("PRESTO", 3))
  condition <- c("NH", "HS-7", "HS-2", "NH", "HS-7", "HS-2")
  
  ComplexHeatmap::HeatmapAnnotation(
    Genotype = genotype,
    Condition = condition,
    col = list(
      Genotype = geno_colors,
      Condition = c(
        "NH" = unname(temp_colors["T0"]),
        "HS-7" = unname(temp_colors["T1"]),
        "HS-2" = unname(temp_colors["T2"])
      )
    ),
    annotation_name_gp = grid::gpar(fontsize = 8),
    simple_anno_size = grid::unit(3.5, "mm"),
    show_legend = TRUE
  )
}

draw_heatmap_object <- function(ht) {
  ComplexHeatmap::draw(
    ht,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    merge_legends = TRUE,
    padding = grid::unit(c(5, 5, 5, 5), "mm")
  )
}

save_complex_heatmap <- function(
    ht, file_base, width_mm, height_mm,
    allow_raster = TRUE) {
  
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  dir.create(dirname(file_base), recursive = TRUE, showWarnings = FALSE)
  
  render_device <- function(label, open_device) {
    device_before <- grDevices::dev.cur()
    tryCatch({
      open_device()
      draw_heatmap_object(ht)
      grDevices::dev.off()
      invisible(TRUE)
    }, error = function(e) {
      if (grDevices::dev.cur() != device_before) {
        try(grDevices::dev.off(), silent = TRUE)
      }
      warning(
        label, " export failed for ", file_base, ": ",
        conditionMessage(e)
      )
      invisible(FALSE)
    })
  }
  
  if (TF_WRITE_PDF) {
    render_device("PDF", function() {
      grDevices::pdf(
        paste0(file_base, ".pdf"),
        width = width_in,
        height = height_in,
        useDingbats = FALSE,
        onefile = TRUE
      )
    })
  }
  
  if (TF_WRITE_SVG) {
    render_device("SVG", function() {
      svglite::svglite(
        paste0(file_base, ".svg"),
        width = width_in,
        height = height_in,
        bg = "white"
      )
    })
  }
  
  raster_ok <- allow_raster && height_mm <= TF_RASTER_MAX_HEIGHT_MM
  
  if (TF_WRITE_PNG && raster_ok) {
    render_device("PNG", function() {
      grDevices::png(
        paste0(file_base, ".png"),
        width = width_in,
        height = height_in,
        units = "in",
        res = TF_RASTER_DPI,
        type = "cairo-png",
        bg = "white"
      )
    })
  }
  
  if (TF_WRITE_TIFF && raster_ok) {
    render_device("TIFF", function() {
      grDevices::tiff(
        paste0(file_base, ".tiff"),
        width = width_in,
        height = height_in,
        units = "in",
        res = TF_RASTER_DPI,
        compression = "lzw",
        type = "cairo",
        bg = "white"
      )
    })
  }
  
  if (!raster_ok && (TF_WRITE_PNG || TF_WRITE_TIFF)) {
    message(
      "Raster skipped for very tall heatmap (", round(height_mm, 1),
      " mm). Complete PDF/SVG were still written: ", file_base
    )
  }
}

get_row_order <- function(mat, cluster_rows = TRUE) {
  if (!cluster_rows || nrow(mat) <= 2L) return(seq_len(nrow(mat)))
  stats::hclust(stats::dist(mat), method = "complete")$order
}


# ------------------------------------------------------------------
# Compact expression-only heatmap used for the global TF panoramas.
#
# Unlike make_heatmap_pair(), this function intentionally does NOT:
#   - split rows by TF family,
#   - show DE status,
#   - show gene labels by default.
#
# This provides a condensed global view similar to a classical
# transcriptomic overview heatmap.
# ------------------------------------------------------------------
make_global_tf_expression_heatmap <- function(
    zmat,
    title,
    cluster_rows = TF_GLOBAL_CLUSTER_ROWS,
    show_row_names = TF_GLOBAL_SHOW_ROW_NAMES,
    show_row_dend = TF_GLOBAL_SHOW_ROW_DENDROGRAM) {
  
  if (nrow(zmat) < 2L) {
    stop("Global TF heatmap requires at least two TFs.")
  }
  
  zplot <- pmax(
    pmin(zmat, TF_GLOBAL_Z_LIMIT),
    -TF_GLOBAL_Z_LIMIT
  )
  
  # Remove rows that are entirely non-finite as a final safety check.
  keep_rows <- apply(
    zplot,
    1L,
    function(x) any(is.finite(x))
  )
  zplot <- zplot[keep_rows, , drop = FALSE]
  
  if (nrow(zplot) < 2L) {
    stop("Fewer than two finite TF-expression rows remain.")
  }
  
  z_col_fun <- circlize::colorRamp2(
    c(-TF_GLOBAL_Z_LIMIT, 0, TF_GLOBAL_Z_LIMIT),
    c("#2166AC", "#F7F7F7", "#B2182B")
  )
  
  ComplexHeatmap::Heatmap(
    zplot,
    name = "Row z-score",
    col = z_col_fun,
    na_col = "grey90",
    
    # Same genotype/temperature annotation as all other Figure 9B outputs.
    top_annotation = make_top_annotation(),
    
    column_title = title,
    column_title_gp = grid::gpar(
      fontsize = 11,
      fontface = "bold"
    ),
    
    # Biological conditions must never be clustered/reordered.
    cluster_columns = FALSE,
    
    # TFs themselves are clustered to reveal broad response patterns.
    cluster_rows = cluster_rows,
    clustering_distance_rows = "euclidean",
    clustering_method_rows = "complete",
    show_row_dend = show_row_dend,
    row_dend_width = grid::unit(12, "mm"),
    
    # Condensed overview: no gene labels and no family separation.
    show_row_names = show_row_names,
    row_names_gp = grid::gpar(fontsize = 5),
    show_column_names = TRUE,
    column_names_rot = 90,
    column_names_gp = grid::gpar(fontsize = 8),
    column_names_centered = TRUE,
    
    # No white cell borders: this keeps hundreds/thousands of TF rows
    # visually continuous, as in a global transcriptomic heatmap.
    rect_gp = grid::gpar(col = NA),
    border = FALSE,
    
    # Rasterise the heatmap body when large, while keeping annotations/text
    # handled by ComplexHeatmap.
    use_raster = nrow(zplot) > 200L,
    raster_quality = TF_GLOBAL_RASTER_QUALITY,
    
    heatmap_legend_param = list(
      title = "Expression",
      at = c(
        -TF_GLOBAL_Z_LIMIT,
        0,
        TF_GLOBAL_Z_LIMIT
      ),
      labels = c(
        paste0("<= -", TF_GLOBAL_Z_LIMIT),
        "0",
        paste0(">= ", TF_GLOBAL_Z_LIMIT)
      )
    )
  )
}

make_heatmap_pair <- function(
    zmat, status_mat, labels, title,
    row_split = NULL, cluster_rows = TRUE,
    show_row_names = TRUE) {
  
  stopifnot(identical(rownames(zmat), rownames(status_mat)))
  stopifnot(all(rownames(zmat) %in% names(labels)))
  
  zplot <- pmax(pmin(zmat, TF_Z_LIMIT), -TF_Z_LIMIT)
  
  if (is.null(row_split)) {
    ord <- get_row_order(zplot, cluster_rows)
    zplot <- zplot[ord, , drop = FALSE]
    status_mat <- status_mat[ord, , drop = FALSE]
    labels <- labels[rownames(zplot)]
  }
  
  z_col_fun <- circlize::colorRamp2(
    c(-TF_Z_LIMIT, 0, TF_Z_LIMIT),
    c("#2166AC", "#F7F7F7", "#B2182B")
  )
  
  hm_expr <- ComplexHeatmap::Heatmap(
    zplot,
    name = "Row z-score",
    col = z_col_fun,
    na_col = "grey90",
    top_annotation = make_top_annotation(),
    column_title = title,
    column_title_gp = grid::gpar(fontsize = 11, fontface = "bold"),
    cluster_columns = FALSE,
    cluster_rows = if (is.null(row_split)) FALSE else cluster_rows,
    cluster_row_slices = FALSE,
    row_split = row_split,
    row_title_gp = grid::gpar(fontsize = 8.5, fontface = "bold"),
    row_title_rot = 0,
    row_labels = labels,
    show_row_names = show_row_names,
    row_names_gp = grid::gpar(fontsize = 6.4),
    row_names_max_width = grid::unit(92, "mm"),
    column_names_rot = 45,
    column_names_gp = grid::gpar(fontsize = 8),
    column_names_centered = TRUE,
    rect_gp = grid::gpar(col = "white", lwd = 0.30),
    border = TRUE,
    heatmap_legend_param = list(
      title = "Expression",
      at = c(-TF_Z_LIMIT, 0, TF_Z_LIMIT),
      labels = c(
        paste0("<= -", TF_Z_LIMIT),
        "0",
        paste0(">= ", TF_Z_LIMIT)
      )
    )
  )
  
  hm_de <- ComplexHeatmap::Heatmap(
    status_mat,
    name = "DE status",
    col = c(
      "Down" = unname(deg_colors["Down"]),
      "NS" = unname(deg_colors["NS"]),
      "Up" = unname(deg_colors["Up"])
    ),
    na_col = unname(deg_colors["NS"]),
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    show_row_names = FALSE,
    column_names_rot = 45,
    column_names_gp = grid::gpar(fontsize = 7),
    rect_gp = grid::gpar(col = "white", lwd = 0.30),
    border = TRUE,
    width = grid::unit(38, "mm"),
    heatmap_legend_param = list(
      title = "Differential\nexpression",
      at = c("Up", "Down", "NS")
    )
  )
  
  hm_expr + hm_de
}

##############################
## 4) Load annotation and map to analysis gene IDs
##############################
if (!file.exists(TF_ANNOT_FILE)) {
  stop("TF annotation file not found: ", TF_ANNOT_FILE)
}
if (!file.exists(map_file)) {
  stop("Gene-ID map file not found: ", map_file)
}

annot <- utils::read.delim(
  TF_ANNOT_FILE,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  quote = "",
  check.names = FALSE
)
colnames(annot) <- trimws(colnames(annot))

assert_columns(
  annot,
  c("gene_ID", "tf_type", "tf_family", "function"),
  "DH13M14 annotation"
)

annot <- annot |>
  dplyr::mutate(
    tf_type = trimws(dplyr::coalesce(as.character(tf_type), "")),
    tf_family = trimws(dplyr::coalesce(as.character(tf_family), "")),
    func_clean = dplyr::coalesce(as.character(.data[["function"]]), ""),
    # A TF call is retained if the source provides either a TF type or family.
    is_TF = nzchar(tf_type) | nzchar(tf_family),
    TF_family_source = dplyr::if_else(
      nzchar(tf_family), tf_family, "Unclassified_TF"
    )
  )

map_df <- utils::read.csv(
  map_file,
  stringsAsFactors = FALSE
) |>
  dplyr::select(my_gene_id, plant2t_id) |>
  dplyr::mutate(
    my_gene_id = trimws(as.character(my_gene_id)),
    plant2t_id = trimws(as.character(plant2t_id))
  ) |>
  dplyr::filter(
    !is.na(my_gene_id), nzchar(my_gene_id),
    !is.na(plant2t_id), nzchar(plant2t_id)
  ) |>
  dplyr::distinct()

annot_map_raw <- annot |>
  dplyr::filter(is_TF) |>
  dplyr::select(
    gene_ID, tf_type, TF_family_source, func_clean
  ) |>
  dplyr::inner_join(
    map_df,
    by = c("gene_ID" = "plant2t_id")
  ) |>
  dplyr::mutate(
    TF_family = canonicalise_tf_family(
      TF_family_source,
      func_clean
    ),
    TF_subgroup = infer_tf_subgroup_conservative(
      TF_family_source,
      TF_family,
      func_clean
    ),
    Annotation_confidence = tf_annotation_confidence(
      TF_family_source,
      TF_family,
      func_clean
    )
  )

# Do not silently collapse a gene mapped to conflicting TF families.
family_conflicts <- annot_map_raw |>
  dplyr::group_by(my_gene_id) |>
  dplyr::summarise(
    n_families = dplyr::n_distinct(TF_family),
    families = paste(sort(unique(TF_family)), collapse = "; "),
    .groups = "drop"
  ) |>
  dplyr::filter(n_families > 1L)

if (nrow(family_conflicts) > 0L) {
  stop(
    "Conflicting TF-family assignments were detected for ",
    nrow(family_conflicts),
    " mapped gene(s). Resolve these mappings before plotting. Examples: ",
    paste(
      paste0(
        head(family_conflicts$my_gene_id, 10L),
        " [",
        head(family_conflicts$families, 10L),
        "]"
      ),
      collapse = ", "
    )
  )
}

annot_map <- annot_map_raw |>
  dplyr::distinct(my_gene_id, .keep_all = TRUE)

cat("Mapped annotated TFs:", nrow(annot_map), "\n")
cat(
  "Mapped TF families:",
  dplyr::n_distinct(annot_map$TF_family),
  "\n"
)

##############################
## 5) Load fitted expression object and the FOUR DE contrasts
##############################
pattern_run <- paste0(
  "^", SAMPLING_TIME_CHOSEN, "_Merged_BothGenotypes_"
)

main_dir <- find_latest_complete_run(
  base_out,
  pattern_run,
  required_relpaths = c(
    "01_edgeR_fit",
    "05_DEG_Contrasts/DEG_Tables"
  )
)

dge_file <- file.path(
  main_dir,
  "01_edgeR_fit",
  paste0("dge_MERGED_", SAMPLING_TIME_CHOSEN, ".rds")
)
if (!file.exists(dge_file)) {
  stop("Merged edgeR DGE object not found: ", dge_file)
}

y <- readRDS(dge_file)
logcpm <- edgeR::cpm(y, log = TRUE, prior.count = 2)

meta2 <- y$samples
meta2$Sample <- rownames(meta2)

if (!"GenoTemp" %in% colnames(meta2)) {
  assert_columns(meta2, c("Genotype", "Temperature"), "DGE sample metadata")
  meta2$GenoTemp <- paste(
    meta2$Genotype,
    meta2$Temperature,
    sep = "_"
  )
}

missing_groups <- setdiff(TF_GROUP_ORDER, unique(as.character(meta2$GenoTemp)))
if (length(missing_groups) > 0L) {
  stop(
    "Missing required GenoTemp group(s): ",
    paste(missing_groups, collapse = ", ")
  )
}

dir_deg <- file.path(
  main_dir,
  "05_DEG_Contrasts",
  "DEG_Tables"
)

read_one_contrast <- function(ct) {
  path <- file.path(dir_deg, paste0("DEG_", ct, ".csv"))
  if (!file.exists(path)) {
    stop("Required DE contrast file not found: ", path)
  }
  
  df <- read_csv_flexible(path)
  if (!"Gene" %in% colnames(df)) {
    if (!is.null(rownames(df)) && all(nzchar(rownames(df)))) {
      df$Gene <- rownames(df)
    } else {
      stop("No Gene column in: ", path)
    }
  }
  
  assert_columns(df, c("Gene", "logFC", "FDR"), basename(path))
  
  df |>
    dplyr::transmute(
      Gene = trimws(as.character(Gene)),
      Contrast = ct,
      logFC = as.numeric(logFC),
      FDR = dplyr::coalesce(as.numeric(FDR), 1),
      Status = status_from_de(
        as.numeric(logFC),
        dplyr::coalesce(as.numeric(FDR), 1)
      )
    )
}

all_deg <- lapply(
  TF_CONTRAST_ORDER,
  read_one_contrast
) |>
  dplyr::bind_rows() |>
  dplyr::left_join(
    annot_map |>
      dplyr::select(
        my_gene_id, gene_ID,
        TF_family_source, TF_family,
        TF_subgroup, Annotation_confidence,
        func_clean
      ),
    by = c("Gene" = "my_gene_id")
  ) |>
  dplyr::mutate(
    is_TF = !is.na(TF_family)
  )

all_sig <- all_deg |>
  dplyr::filter(Status %in% c("Up", "Down"))

all_tf_ids <- unique(annot_map$my_gene_id)
expressed_tf_ids <- intersect(all_tf_ids, rownames(logcpm))
de_tf_ids <- unique(all_sig$Gene[all_sig$is_TF])

cat("Expressed TFs:", length(expressed_tf_ids), "\n")
cat("DE TFs in >=1 intended heat contrast:", length(de_tf_ids), "\n")

##############################
## 6) Manual curation / overrides
##############################
tf_out <- file.path(
  main_dir,
  paste0(
    "09B_TF_panorama_complete_",
    SAMPLING_TIME_CHOSEN, "_", timestamp
  )
)
dir_tables <- file.path(tf_out, "00_tables")
dir_panorama <- file.path(tf_out, "01_panorama")
dir_family <- file.path(tf_out, "02_family_heatmaps")
dir_qc <- file.path(tf_out, "03_QC")
for (d in c(tf_out, dir_tables, dir_panorama, dir_family, dir_qc)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

TF_OVERRIDE_FILE <- file.path(
  dir_tables,
  "TF_manual_overrides.csv"
)

curation_base <- annot_map |>
  dplyr::filter(my_gene_id %in% de_tf_ids) |>
  dplyr::transmute(
    Gene = my_gene_id,
    gene_ID,
    TF_family_source,
    TF_family_auto = TF_family,
    TF_subgroup_auto = TF_subgroup,
    Annotation_confidence,
    func_clean,
    Ambiguous_family = TF_family %in% c(
      "Unclassified_TF",
      "AP2_ERF_unresolved"
    ),
    Include_auto = if (TF_EXCLUDE_AMBIGUOUS_BY_DEFAULT) {
      !Ambiguous_family
    } else {
      TRUE
    }
  )

if (!file.exists(TF_OVERRIDE_FILE)) {
  template <- curation_base |>
    dplyr::transmute(
      Gene,
      gene_ID,
      TF_family_source,
      TF_family_auto,
      TF_subgroup_auto,
      Annotation_confidence,
      Family_override = "",
      Subgroup_override = "",
      Short_name_override = "",
      Include_override = "",
      Curator_note = ""
    )
  readr::write_csv(template, TF_OVERRIDE_FILE)
  message("Created manual override template: ", TF_OVERRIDE_FILE)
}

ov <- read_csv_flexible(TF_OVERRIDE_FILE)
assert_columns(ov, "Gene", "TF manual override table")

for (nm in c(
  "Family_override", "Subgroup_override", "Short_name_override",
  "Include_override", "Curator_note"
)) {
  if (!nm %in% colnames(ov)) ov[[nm]] <- ""
}

ov <- ov |>
  dplyr::mutate(
    dplyr::across(
      c(
        Family_override, Subgroup_override,
        Short_name_override, Include_override, Curator_note
      ),
      ~ trimws(as.character(.x))
    ),
    Include_override_parsed = parse_logical_override(
      Include_override
    )
  ) |>
  dplyr::select(
    Gene,
    Family_override,
    Subgroup_override,
    Short_name_override,
    Include_override_parsed,
    Curator_note
  )

curation <- curation_base |>
  dplyr::left_join(ov, by = "Gene") |>
  dplyr::mutate(
    TF_family = dplyr::if_else(
      !is.na(Family_override) & nzchar(Family_override),
      Family_override,
      TF_family_auto
    ),
    TF_subgroup = dplyr::if_else(
      !is.na(Subgroup_override) & nzchar(Subgroup_override),
      Subgroup_override,
      dplyr::coalesce(TF_subgroup_auto, "")
    ),
    Short_name = dplyr::if_else(
      !is.na(Short_name_override) & nzchar(Short_name_override),
      Short_name_override,
      clean_label(func_clean)
    ),
    Include = dplyr::if_else(
      !is.na(Include_override_parsed),
      Include_override_parsed,
      Include_auto
    ),
    Curator_note = dplyr::coalesce(Curator_note, "")
  )

##############################
## 7) Restrict DE selection to the FOUR intended contrasts
##############################
selected_gene_stats <- all_deg |>
  dplyr::filter(
    is_TF,
    Contrast %in% TF_CONTRAST_ORDER
  ) |>
  dplyr::group_by(Gene) |>
  dplyr::summarise(
    n_DE_contrasts = sum(
      Status %in% c("Up", "Down"),
      na.rm = TRUE
    ),
    min_FDR = suppressWarnings(
      min(
        FDR[Status %in% c("Up", "Down")],
        na.rm = TRUE
      )
    ),
    max_abs_logFC = suppressWarnings(
      max(
        abs(logFC[Status %in% c("Up", "Down")]),
        na.rm = TRUE
      )
    ),
    response_direction = dplyr::case_when(
      any(Status == "Up", na.rm = TRUE) &
        any(Status == "Down", na.rm = TRUE) ~ "Mixed",
      any(Status == "Up", na.rm = TRUE) ~ "Up only",
      any(Status == "Down", na.rm = TRUE) ~ "Down only",
      TRUE ~ "NS"
    ),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    min_FDR = dplyr::if_else(
      is.infinite(min_FDR), NA_real_, min_FDR
    ),
    max_abs_logFC = dplyr::if_else(
      is.infinite(max_abs_logFC),
      NA_real_,
      max_abs_logFC
    )
  ) |>
  dplyr::filter(n_DE_contrasts >= 1L)

plot_tbl <- selected_gene_stats |>
  dplyr::inner_join(
    curation,
    by = "Gene"
  ) |>
  dplyr::filter(Include)

if (nrow(plot_tbl) == 0L) {
  stop("No DE TFs remain after curation.")
}

##############################
## 8) Complete six-group expression matrix
##############################
group_means_all <- vapply(
  TF_GROUP_ORDER,
  FUN.VALUE = numeric(nrow(logcpm)),
  FUN = function(g) {
    samps <- meta2$Sample[
      as.character(meta2$GenoTemp) == g
    ]
    if (length(samps) == 0L) {
      stop("No samples for group: ", g)
    }
    rowMeans(
      logcpm[, samps, drop = FALSE],
      na.rm = TRUE
    )
  }
)
rownames(group_means_all) <- rownames(logcpm)

plot_genes <- plot_tbl$Gene
missing_expr <- setdiff(plot_genes, rownames(group_means_all))
if (length(missing_expr) > 0L) {
  warning(
    length(missing_expr),
    " selected DE TF(s) are absent from the expression matrix and will be removed."
  )
}

plot_tbl <- plot_tbl |>
  dplyr::filter(Gene %in% rownames(group_means_all))

plot_genes <- plot_tbl$Gene
mean_logcpm <- group_means_all[
  plot_genes,
  TF_GROUP_ORDER,
  drop = FALSE
]
colnames(mean_logcpm) <- TF_GROUP_LABELS

if (any(!is.finite(mean_logcpm))) {
  bad <- rownames(mean_logcpm)[
    apply(mean_logcpm, 1L, function(x) any(!is.finite(x)))
  ]
  stop(
    "Non-finite group means for selected TF(s): ",
    paste(bad, collapse = ", ")
  )
}

zscore_mat <- zscore_rows_safe(mean_logcpm)

##############################
## 9) Four-contrast DE-status matrix
##############################
status_df <- all_deg |>
  dplyr::filter(
    Gene %in% plot_genes,
    Contrast %in% TF_CONTRAST_ORDER
  ) |>
  dplyr::select(Gene, Contrast, Status) |>
  dplyr::distinct() |>
  tidyr::complete(
    Gene = plot_genes,
    Contrast = TF_CONTRAST_ORDER,
    fill = list(Status = "NS")
  ) |>
  tidyr::pivot_wider(
    names_from = Contrast,
    values_from = Status,
    values_fill = "NS"
  ) |>
  dplyr::select(
    Gene,
    dplyr::all_of(TF_CONTRAST_ORDER)
  )

status_mat <- as.matrix(
  status_df[, TF_CONTRAST_ORDER, drop = FALSE]
)
rownames(status_mat) <- status_df$Gene
colnames(status_mat) <- TF_CONTRAST_LABELS
status_mat[
  !status_mat %in% c("Up", "Down", "NS")
] <- "NS"
status_mat <- status_mat[plot_genes, , drop = FALSE]

##############################
## 10) Condition-response patterns
##############################
status_named <- as.data.frame(status_mat, stringsAsFactors = FALSE) |>
  tibble::rownames_to_column("Gene")

response_patterns <- status_named |>
  dplyr::transmute(
    Gene,
    ROBILA_HS7_status = .data[["ROBILA HS-7"]],
    ROBILA_HS2_status = .data[["ROBILA HS-2"]],
    PRESTO_HS7_status = .data[["PRESTO HS-7"]],
    PRESTO_HS2_status = .data[["PRESTO HS-2"]],
    ROBILA_pattern = response_pattern_two(
      ROBILA_HS7_status,
      ROBILA_HS2_status
    ),
    PRESTO_pattern = response_pattern_two(
      PRESTO_HS7_status,
      PRESTO_HS2_status
    ),
    Cross_genotype_pattern = paste0(
      "ROBILA: ", ROBILA_pattern,
      " | PRESTO: ", PRESTO_pattern
    )
  )

plot_tbl <- plot_tbl |>
  dplyr::left_join(response_patterns, by = "Gene")

##############################
## 11) Family census — genome, expressed, DE, Up/Down
##############################
genome_tf <- annot_map |>
  dplyr::distinct(my_gene_id, TF_family)

expressed_tf <- annot_map |>
  dplyr::filter(my_gene_id %in% expressed_tf_ids) |>
  dplyr::distinct(my_gene_id, TF_family)

de_tf <- all_sig |>
  dplyr::filter(is_TF) |>
  dplyr::distinct(Gene, TF_family)

up_tf <- all_sig |>
  dplyr::filter(is_TF, Status == "Up") |>
  dplyr::distinct(Gene, TF_family)

down_tf <- all_sig |>
  dplyr::filter(is_TF, Status == "Down") |>
  dplyr::distinct(Gene, TF_family)

census <- genome_tf |>
  dplyr::count(TF_family, name = "n_genome") |>
  dplyr::full_join(
    expressed_tf |>
      dplyr::count(TF_family, name = "n_expressed"),
    by = "TF_family"
  ) |>
  dplyr::full_join(
    de_tf |>
      dplyr::count(TF_family, name = "n_DE"),
    by = "TF_family"
  ) |>
  dplyr::full_join(
    up_tf |>
      dplyr::count(TF_family, name = "n_Up"),
    by = "TF_family"
  ) |>
  dplyr::full_join(
    down_tf |>
      dplyr::count(TF_family, name = "n_Down"),
    by = "TF_family"
  ) |>
  dplyr::mutate(
    dplyr::across(
      dplyr::starts_with("n_"),
      ~ tidyr::replace_na(.x, 0L)
    ),
    pct_expressed = dplyr::if_else(
      n_genome > 0,
      100 * n_expressed / n_genome,
      NA_real_
    ),
    pct_DE_of_expressed = dplyr::if_else(
      n_expressed > 0,
      100 * n_DE / n_expressed,
      NA_real_
    )
  ) |>
  dplyr::arrange(
    dplyr::desc(n_DE),
    dplyr::desc(pct_DE_of_expressed),
    TF_family
  )

##############################
## 12) Family overrepresentation among heat-responsive TFs
##############################
n_expr_total <- length(unique(expressed_tf$my_gene_id))
n_de_total <- length(unique(de_tf$Gene))

family_enrichment <- census |>
  dplyr::filter(n_expressed > 0) |>
  dplyr::rowwise() |>
  dplyr::mutate(
    a_DE_in_family = n_DE,
    b_nonDE_in_family = n_expressed - n_DE,
    c_DE_other = n_de_total - n_DE,
    d_nonDE_other = (
      n_expr_total - n_expressed
    ) - c_DE_other,
    Fisher_p = {
      tab <- matrix(
        c(
          a_DE_in_family,
          b_nonDE_in_family,
          c_DE_other,
          d_nonDE_other
        ),
        nrow = 2L,
        byrow = TRUE
      )
      if (all(tab >= 0)) {
        stats::fisher.test(
          tab,
          alternative = "greater"
        )$p.value
      } else {
        NA_real_
      }
    },
    Odds_ratio = {
      tab <- matrix(
        c(
          a_DE_in_family,
          b_nonDE_in_family,
          c_DE_other,
          d_nonDE_other
        ),
        nrow = 2L,
        byrow = TRUE
      )
      if (all(tab >= 0)) {
        unname(
          stats::fisher.test(
            tab,
            alternative = "greater"
          )$estimate
        )
      } else {
        NA_real_
      }
    }
  ) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    Fisher_FDR = stats::p.adjust(
      Fisher_p,
      method = "BH"
    ),
    Enriched_FDR05 = !is.na(Fisher_FDR) &
      Fisher_FDR < 0.05
  ) |>
  dplyr::arrange(
    Fisher_FDR,
    dplyr::desc(Odds_ratio),
    dplyr::desc(n_DE)
  )

##############################
## 13) Per-contrast family counts and fractions
##############################
family_per_contrast <- all_deg |>
  dplyr::filter(
    is_TF,
    Contrast %in% TF_CONTRAST_ORDER,
    Status %in% c("Up", "Down")
  ) |>
  dplyr::distinct(
    Gene, Contrast, TF_family, Status
  ) |>
  dplyr::count(
    TF_family, Contrast, Status,
    name = "n"
  ) |>
  dplyr::left_join(
    census |>
      dplyr::select(
        TF_family,
        n_expressed
      ),
    by = "TF_family"
  ) |>
  dplyr::mutate(
    pct_of_expressed_family = dplyr::if_else(
      n_expressed > 0,
      100 * n / n_expressed,
      NA_real_
    )
  )

family_contrast_total <- all_deg |>
  dplyr::filter(
    is_TF,
    Contrast %in% TF_CONTRAST_ORDER,
    Status %in% c("Up", "Down")
  ) |>
  dplyr::distinct(
    Gene, Contrast, TF_family
  ) |>
  dplyr::count(
    TF_family, Contrast,
    name = "n_DE"
  ) |>
  tidyr::complete(
    TF_family = sort(unique(census$TF_family)),
    Contrast = TF_CONTRAST_ORDER,
    fill = list(n_DE = 0L)
  ) |>
  dplyr::left_join(
    census |>
      dplyr::select(
        TF_family, n_expressed
      ),
    by = "TF_family"
  ) |>
  dplyr::mutate(
    pct_DE_of_expressed_family = dplyr::if_else(
      n_expressed > 0,
      100 * n_DE / n_expressed,
      NA_real_
    )
  )

##############################
## 14) Family completeness audit
##############################
analysis_families <- sort(unique(annot_map$TF_family))

family_audit <- tibble::tibble(
  Family = sort(unique(c(
    PLANTTFDB_DCA_FAMILIES,
    analysis_families
  )))
) |>
  dplyr::mutate(
    In_current_PlantTFDB_Dca_reference = Family %in%
      PLANTTFDB_DCA_FAMILIES,
    In_DH13M14_source_annotation = Family %in%
      analysis_families,
    Status = dplyr::case_when(
      In_current_PlantTFDB_Dca_reference &
        In_DH13M14_source_annotation ~ "Shared",
      In_current_PlantTFDB_Dca_reference &
        !In_DH13M14_source_annotation ~ "Reference-only / naming-or-annotation difference",
      !In_current_PlantTFDB_Dca_reference &
        In_DH13M14_source_annotation ~ "DH13M14-only / naming-or-classification difference",
      TRUE ~ "Neither"
    )
  )

##############################
## 15) Family ordering for complete panorama
##############################
if (TF_FAMILY_ORDER_MODE == "n_DE") {
  family_order_base <- census |>
    dplyr::filter(n_DE > 0) |>
    dplyr::arrange(
      dplyr::desc(n_DE),
      dplyr::desc(pct_DE_of_expressed),
      TF_family
    ) |>
    dplyr::pull(TF_family)
} else if (TF_FAMILY_ORDER_MODE == "alphabetical") {
  family_order_base <- sort(
    unique(as.character(plot_tbl$TF_family))
  )
} else {
  stop(
    "TF_FAMILY_ORDER_MODE must be 'n_DE' or 'alphabetical'."
  )
}

# Manual family overrides must remain plottable even if their new family name
# was absent from the pre-override census.
families_after_curation <- unique(as.character(plot_tbl$TF_family))
family_order <- c(
  family_order_base[family_order_base %in% families_after_curation],
  sort(setdiff(families_after_curation, family_order_base))
)
family_order <- unique(family_order)

plot_tbl <- plot_tbl |>
  dplyr::mutate(
    TF_family = factor(
      as.character(TF_family),
      levels = family_order
    )
  ) |>
  dplyr::filter(!is.na(TF_family)) |>
  dplyr::arrange(
    TF_family,
    TF_subgroup,
    Gene
  )

plot_genes <- plot_tbl$Gene
mean_logcpm <- mean_logcpm[
  plot_genes, ,
  drop = FALSE
]
zscore_mat <- zscore_mat[
  plot_genes, ,
  drop = FALSE
]
status_mat <- status_mat[
  plot_genes, ,
  drop = FALSE
]

labels_all <- stats::setNames(
  ifelse(
    nzchar(plot_tbl$Short_name),
    paste0(
      plot_tbl$Gene,
      " | ",
      plot_tbl$Short_name
    ),
    plot_tbl$Gene
  ),
  plot_tbl$Gene
)

##############################
## 16) Export complete tables
##############################
readr::write_csv(
  census,
  file.path(
    dir_tables,
    "TF_census_all_families.csv"
  )
)

readr::write_csv(
  family_enrichment,
  file.path(
    dir_tables,
    "TF_family_overrepresentation_Fisher.csv"
  )
)

readr::write_csv(
  family_per_contrast,
  file.path(
    dir_tables,
    "TF_family_direction_per_contrast.csv"
  )
)

readr::write_csv(
  family_contrast_total,
  file.path(
    dir_tables,
    "TF_family_total_per_contrast.csv"
  )
)

readr::write_csv(
  family_audit,
  file.path(
    dir_tables,
    "TF_family_completeness_audit_vs_PlantTFDB_Dca.csv"
  )
)

readr::write_csv(
  curation,
  file.path(
    dir_tables,
    "TF_DE_curation_complete.csv"
  )
)

readr::write_csv(
  plot_tbl,
  file.path(
    dir_tables,
    "TF_DE_selected_complete.csv"
  )
)

readr::write_csv(
  response_patterns,
  file.path(
    dir_tables,
    "TF_DE_response_patterns.csv"
  )
)

readr::write_csv(
  tibble::rownames_to_column(
    as.data.frame(mean_logcpm),
    "Gene"
  ) |>
    dplyr::left_join(
      plot_tbl |>
        dplyr::mutate(
          TF_family = as.character(TF_family)
        ) |>
        dplyr::select(
          Gene, TF_family, TF_subgroup,
          Short_name
        ),
      by = "Gene"
    ),
  file.path(
    dir_tables,
    "TF_DE_meanLogCPM_six_groups.csv"
  )
)

readr::write_csv(
  tibble::rownames_to_column(
    as.data.frame(zscore_mat),
    "Gene"
  ) |>
    dplyr::left_join(
      plot_tbl |>
        dplyr::mutate(
          TF_family = as.character(TF_family)
        ) |>
        dplyr::select(
          Gene, TF_family, TF_subgroup,
          Short_name
        ),
      by = "Gene"
    ),
  file.path(
    dir_tables,
    "TF_DE_zscore_six_groups.csv"
  )
)

readr::write_csv(
  tibble::rownames_to_column(
    as.data.frame(status_mat),
    "Gene"
  ),
  file.path(
    dir_tables,
    "TF_DE_status_four_contrasts.csv"
  )
)

##############################
## 17) GLOBAL COMPACT TF EXPRESSION PANORAMAS
##
## 17A = ALL expressed TFs
## 17B = ALL heat-responsive / DE TFs
##
## These are expression-only overview figures:
## no family split, no gene names, no DE-status side panel.
##############################

# ================================================================
# 17A) ALL EXPRESSED TFs — complete global expression landscape
# ================================================================
global_all_tf_ids <- intersect(
  annot_map$my_gene_id,
  rownames(group_means_all)
)
global_all_tf_ids <- unique(global_all_tf_ids)

if (length(global_all_tf_ids) >= 2L) {
  
  global_all_mean_logcpm <- group_means_all[
    global_all_tf_ids,
    TF_GROUP_ORDER,
    drop = FALSE
  ]
  colnames(global_all_mean_logcpm) <- TF_GROUP_LABELS
  
  # Safety: remove any gene with a non-finite group mean.
  finite_global_all <- apply(
    global_all_mean_logcpm,
    1L,
    function(x) all(is.finite(x))
  )
  global_all_mean_logcpm <- global_all_mean_logcpm[
    finite_global_all,
    ,
    drop = FALSE
  ]
  
  global_all_zscore <- zscore_rows_safe(
    global_all_mean_logcpm
  )
  
  ht_global_all <- make_global_tf_expression_heatmap(
    zmat = global_all_zscore,
    title = paste0(
      "Global expression panorama of all expressed transcription factors (n=",
      nrow(global_all_zscore),
      ")"
    )
  )
  
  # Fixed compact dimensions on purpose:
  # rows are squeezed into a global landscape rather than displayed one by one.
  save_complex_heatmap(
    ht_global_all,
    file.path(
      dir_panorama,
      "Figure_9B_GLOBAL_ALL_EXPRESSED_TFs"
    ),
    width_mm = TF_GLOBAL_WIDTH_MM,
    height_mm = TF_GLOBAL_ALL_HEIGHT_MM,
    allow_raster = TRUE
  )
  
  # Export the exact matrices used in the figure.
  readr::write_csv(
    tibble::rownames_to_column(
      as.data.frame(global_all_mean_logcpm),
      "Gene"
    ) |>
      dplyr::left_join(
        annot_map |>
          dplyr::transmute(
            Gene = my_gene_id,
            TF_family,
            TF_family_source,
            TF_subgroup,
            func_clean
          ),
        by = "Gene"
      ),
    file.path(
      dir_tables,
      "TF_GLOBAL_ALL_EXPRESSED_meanLogCPM.csv"
    )
  )
  
  readr::write_csv(
    tibble::rownames_to_column(
      as.data.frame(global_all_zscore),
      "Gene"
    ) |>
      dplyr::left_join(
        annot_map |>
          dplyr::transmute(
            Gene = my_gene_id,
            TF_family,
            TF_family_source,
            TF_subgroup
          ),
        by = "Gene"
      ),
    file.path(
      dir_tables,
      "TF_GLOBAL_ALL_EXPRESSED_zscore.csv"
    )
  )
  
  cat(
    "Saved global ALL-expressed TF heatmap:",
    nrow(global_all_zscore),
    "TFs\n"
  )
  
} else {
  warning(
    "Fewer than two expressed TFs were available for the global ALL-TF heatmap."
  )
}


# ================================================================
# 17B) ALL DE TFs — compact expression-only heatmap
#
# This contains exactly the TF set retained in the detailed Figure 9B:
# DE (Up or Down) in >=1 of the four intended heat-vs-NH contrasts.
# It is useful when the family-split detailed panel is too dense for a slide
# or for an initial overview.
# ================================================================
global_de_tf_ids <- plot_tbl$Gene
global_de_tf_ids <- global_de_tf_ids[
  global_de_tf_ids %in% rownames(group_means_all)
]
global_de_tf_ids <- unique(global_de_tf_ids)

if (length(global_de_tf_ids) >= 2L) {
  
  global_de_mean_logcpm <- group_means_all[
    global_de_tf_ids,
    TF_GROUP_ORDER,
    drop = FALSE
  ]
  colnames(global_de_mean_logcpm) <- TF_GROUP_LABELS
  
  finite_global_de <- apply(
    global_de_mean_logcpm,
    1L,
    function(x) all(is.finite(x))
  )
  global_de_mean_logcpm <- global_de_mean_logcpm[
    finite_global_de,
    ,
    drop = FALSE
  ]
  
  global_de_zscore <- zscore_rows_safe(
    global_de_mean_logcpm
  )
  
  ht_global_de <- make_global_tf_expression_heatmap(
    zmat = global_de_zscore,
    title = paste0(
      "Global expression panorama of heat-responsive transcription factors (n=",
      nrow(global_de_zscore),
      ")"
    )
  )
  
  save_complex_heatmap(
    ht_global_de,
    file.path(
      dir_panorama,
      "Figure_9B_GLOBAL_DE_TFs_expression_only"
    ),
    width_mm = TF_GLOBAL_WIDTH_MM,
    height_mm = TF_GLOBAL_DE_HEIGHT_MM,
    allow_raster = TRUE
  )
  
  readr::write_csv(
    tibble::rownames_to_column(
      as.data.frame(global_de_mean_logcpm),
      "Gene"
    ) |>
      dplyr::left_join(
        plot_tbl |>
          dplyr::mutate(
            TF_family = as.character(TF_family)
          ) |>
          dplyr::select(
            Gene,
            TF_family,
            TF_subgroup,
            Short_name,
            n_DE_contrasts,
            min_FDR,
            max_abs_logFC,
            response_direction
          ),
        by = "Gene"
      ),
    file.path(
      dir_tables,
      "TF_GLOBAL_DE_meanLogCPM.csv"
    )
  )
  
  readr::write_csv(
    tibble::rownames_to_column(
      as.data.frame(global_de_zscore),
      "Gene"
    ) |>
      dplyr::left_join(
        plot_tbl |>
          dplyr::mutate(
            TF_family = as.character(TF_family)
          ) |>
          dplyr::select(
            Gene,
            TF_family,
            TF_subgroup,
            n_DE_contrasts,
            response_direction
          ),
        by = "Gene"
      ),
    file.path(
      dir_tables,
      "TF_GLOBAL_DE_zscore.csv"
    )
  )
  
  cat(
    "Saved global DE-TF expression-only heatmap:",
    nrow(global_de_zscore),
    "TFs\n"
  )
  
} else {
  warning(
    "Fewer than two DE TFs were available for the global DE-TF heatmap."
  )
}


##############################
## 18) COMPLETE DE-TF panorama heatmap — NO top-N truncation
##      Detailed version: TF-family split + Up/Down/NS panel
##############################
combined_split <- droplevels(
  plot_tbl$TF_family
)

ht_complete <- make_heatmap_pair(
  zmat = zscore_mat,
  status_mat = status_mat,
  labels = labels_all,
  title = paste0(
    "Figure 9B | Complete heat-responsive TF panorama — ",
    nrow(plot_tbl),
    " DE TFs in ",
    length(unique(as.character(plot_tbl$TF_family))),
    " families"
  ),
  row_split = combined_split,
  cluster_rows = TF_CLUSTER_ROWS,
  show_row_names = TRUE
)

complete_height_mm <- max(
  190,
  65 + TF_ROW_MM * nrow(plot_tbl)
)

save_complex_heatmap(
  ht_complete,
  file.path(
    dir_panorama,
    "Figure_9B_COMPLETE_DE_TF_panorama"
  ),
  width_mm = TF_FIG_WIDTH_MM + 25,
  height_mm = complete_height_mm,
  allow_raster = TRUE
)

##############################
## 19) Optional context: top-variable ALL EXPRESSED TFs
##     (does not affect main DE panorama)
##############################
all_expr_info <- annot_map |>
  dplyr::filter(
    my_gene_id %in% expressed_tf_ids
  ) |>
  dplyr::distinct(
    my_gene_id, .keep_all = TRUE
  )

if (
  TF_ALL_EXPRESSED_TOP_N > 0L &&
  nrow(all_expr_info) >= 2L
) {
  rv <- apply(
    logcpm[
      all_expr_info$my_gene_id,
      ,
      drop = FALSE
    ],
    1L,
    stats::var,
    na.rm = TRUE
  )
  
  top_ids <- names(
    sort(rv, decreasing = TRUE)
  )[
    seq_len(
      min(
        TF_ALL_EXPRESSED_TOP_N,
        length(rv)
      )
    )
  ]
  
  top_means <- group_means_all[
    top_ids,
    TF_GROUP_ORDER,
    drop = FALSE
  ]
  colnames(top_means) <- TF_GROUP_LABELS
  top_z <- zscore_rows_safe(top_means)
  
  top_status_df <- all_deg |>
    dplyr::filter(
      Gene %in% top_ids,
      Contrast %in% TF_CONTRAST_ORDER
    ) |>
    dplyr::select(
      Gene, Contrast, Status
    ) |>
    tidyr::complete(
      Gene = top_ids,
      Contrast = TF_CONTRAST_ORDER,
      fill = list(Status = "NS")
    ) |>
    tidyr::pivot_wider(
      names_from = Contrast,
      values_from = Status,
      values_fill = "NS"
    ) |>
    dplyr::select(
      Gene,
      dplyr::all_of(TF_CONTRAST_ORDER)
    )
  
  top_status <- as.matrix(
    top_status_df[
      , TF_CONTRAST_ORDER,
      drop = FALSE
    ]
  )
  rownames(top_status) <- top_status_df$Gene
  colnames(top_status) <- TF_CONTRAST_LABELS
  top_status <- top_status[top_ids, , drop = FALSE]
  
  top_info <- all_expr_info[
    match(top_ids, all_expr_info$my_gene_id),
    ,
    drop = FALSE
  ]
  
  top_family_order <- census |>
    dplyr::filter(
      TF_family %in% top_info$TF_family
    ) |>
    dplyr::arrange(
      dplyr::desc(n_DE),
      TF_family
    ) |>
    dplyr::pull(TF_family)
  
  top_info$TF_family_factor <- factor(
    top_info$TF_family,
    levels = unique(top_family_order)
  )
  
  ord <- order(
    top_info$TF_family_factor,
    top_info$my_gene_id
  )
  top_ids <- top_ids[ord]
  top_z <- top_z[top_ids, , drop = FALSE]
  top_status <- top_status[top_ids, , drop = FALSE]
  top_info <- top_info[ord, , drop = FALSE]
  
  top_labels <- stats::setNames(
    paste0(
      top_info$my_gene_id,
      " | ",
      clean_label(top_info$func_clean, 55L)
    ),
    top_info$my_gene_id
  )
  
  ht_all_context <- make_heatmap_pair(
    zmat = top_z,
    status_mat = top_status,
    labels = top_labels,
    title = paste0(
      "Top-variable expressed TFs — context only (n=",
      length(top_ids), ")"
    ),
    row_split = droplevels(
      top_info$TF_family_factor
    ),
    cluster_rows = TF_CLUSTER_ROWS,
    show_row_names = TRUE
  )
  
  context_height_mm <- max(
    180,
    60 + 2.7 * length(top_ids)
  )
  
  save_complex_heatmap(
    ht_all_context,
    file.path(
      dir_panorama,
      "Context_TOP_VARIABLE_ALL_EXPRESSED_TFs"
    ),
    width_mm = TF_FIG_WIDTH_MM + 20,
    height_mm = context_height_mm,
    allow_raster = TRUE
  )
}

##############################
## 20) Family x contrast count heatmap
##############################
count_wide <- family_contrast_total |>
  dplyr::select(
    TF_family, Contrast, n_DE
  ) |>
  tidyr::pivot_wider(
    names_from = Contrast,
    values_from = n_DE,
    values_fill = 0
  )

count_mat <- as.matrix(
  count_wide[
    , TF_CONTRAST_ORDER,
    drop = FALSE
  ]
)
rownames(count_mat) <- count_wide$TF_family
colnames(count_mat) <- TF_CONTRAST_LABELS

# Keep families with >=1 DE TF in any contrast.
count_mat <- count_mat[
  rowSums(count_mat) > 0,
  ,
  drop = FALSE
]

if (nrow(count_mat) >= 1L) {
  family_count_order <- family_order[
    family_order %in% rownames(count_mat)
  ]
  count_mat <- count_mat[
    family_count_order,
    ,
    drop = FALSE
  ]
  
  max_count <- max(count_mat, na.rm = TRUE)
  count_col_fun <- circlize::colorRamp2(
    c(0, max(1, max_count / 2), max(1, max_count)),
    c("#FFFFFF", "#FDBE85", "#D7301F")
  )
  
  ht_counts <- ComplexHeatmap::Heatmap(
    count_mat,
    name = "n DE TFs",
    col = count_col_fun,
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    rect_gp = grid::gpar(
      col = "white",
      lwd = 0.35
    ),
    border = TRUE,
    column_title = "DE transcription-factor counts per family and contrast",
    row_names_gp = grid::gpar(fontsize = 7.5),
    column_names_gp = grid::gpar(fontsize = 8),
    cell_fun = function(j, i, x, y, w, h, fill) {
      grid::grid.text(
        count_mat[i, j],
        x, y,
        gp = grid::gpar(fontsize = 7)
      )
    }
  )
  
  save_complex_heatmap(
    ht_counts,
    file.path(
      dir_panorama,
      "TF_family_counts_per_contrast"
    ),
    width_mm = 150,
    height_mm = max(
      120,
      45 + 5.0 * nrow(count_mat)
    ),
    allow_raster = TRUE
  )
}

##############################
## 21) One DE-only heatmap for EVERY represented TF family
##############################
families_present <- levels(
  droplevels(plot_tbl$TF_family)
)

for (fam in families_present) {
  fam_tbl <- plot_tbl |>
    dplyr::filter(
      as.character(TF_family) == fam
    ) |>
    dplyr::arrange(
      TF_subgroup,
      Gene
    )
  
  genes <- fam_tbl$Gene
  if (length(genes) == 0L) next
  
  z_fam <- zscore_mat[
    genes,
    ,
    drop = FALSE
  ]
  st_fam <- status_mat[
    genes,
    ,
    drop = FALSE
  ]
  lab_fam <- labels_all[genes]
  
  ht_fam <- make_heatmap_pair(
    zmat = z_fam,
    status_mat = st_fam,
    labels = lab_fam,
    title = paste0(
      fam,
      " — all DE TFs in >=1 heat-vs-NH contrast (n=",
      length(genes),
      ")"
    ),
    row_split = NULL,
    cluster_rows = TF_CLUSTER_ROWS,
    show_row_names = TRUE
  )
  
  height_mm <- max(
    95,
    42 + TF_ROW_MM * length(genes)
  )
  
  save_complex_heatmap(
    ht_fam,
    file.path(
      dir_family,
      paste0(
        "Figure_9B_",
        safe_file_stub(fam),
        "_DE_heatmap"
      )
    ),
    width_mm = TF_FIG_WIDTH_MM,
    height_mm = height_mm,
    allow_raster = TRUE
  )
  
  cat(
    "Saved TF family:",
    fam,
    "—",
    length(genes),
    "DE genes\n"
  )
}

##############################
## 22) Bibliography / methodological rules table
##############################
bibliography_rules <- tibble::tribble(
  ~Topic, ~Rule_used_here, ~Reason, ~Reference,
  "TF family definition",
  "Use source TF-family calls; do not reconstruct all families from free-text function keywords.",
  "Modern PlantTFDB family assignment is domain-based and can use DNA-binding, auxiliary and forbidden domains.",
  "Tian et al. 2020, Nucleic Acids Research, DOI:10.1093/nar/gkz1020; PlantTFDB v5 family-assignment rules",
  "Daucus carota completeness",
  "Do not prefilter to a small focus-family list.",
  "Current PlantTFDB D. carota reference contains 1906 TFs in 56 families.",
  "PlantTFDB v5 Daucus carota species page",
  "ERF / DREB",
  "Map explicit DREB source calls into top-level ERF and retain DREB as subgroup.",
  "DREB and ERF are AP2/ERF-domain factors; DREB is not treated as an independent PlantTFDB top-level family.",
  "Nakano et al. 2006, Plant Physiology, DOI:10.1104/pp.105.073783",
  "Subfamily inference",
  "Assign a subfamily only when it is explicit in the annotation; otherwise leave unresolved.",
  "Many TF subgroup definitions depend on protein sequence, motif composition or phylogeny.",
  "PlantTFDB v5 family-assignment rules; family-specific literature",
  "Heat-responsive TF panorama",
  "Select every annotated TF that is DE in >=1 of the four intended heat-vs-NH contrasts.",
  "Heat transcriptomes recover broad TF repertoires; restricting to canonical families can hide real regulatory responses.",
  "Genome-wide heat-stress transcriptomic studies in rice and other crops",
  "Family interpretation",
  "Report n_DE/n_expressed and test family overrepresentation, not n_DE alone.",
  "Large TF families can dominate raw counts simply because they contain more expressed genes.",
  "Statistical rationale: one-sided Fisher exact test with BH correction"
)
readr::write_csv(
  bibliography_rules,
  file.path(
    dir_tables,
    "TF_methodological_literature_rules.csv"
  )
)

##############################
## 23) QC / reproducibility
##############################
utils::write.csv(
  meta2,
  file.path(
    dir_qc,
    "metadata_from_DGE_used.csv"
  ),
  row.names = FALSE
)

readr::write_csv(
  annot_map,
  file.path(
    dir_qc,
    "TF_annotation_mapped_all.csv"
  )
)

readr::write_csv(
  all_deg,
  file.path(
    dir_qc,
    "TF_all_four_contrasts_with_annotation.csv"
  )
)

grand_summary <- tibble::tibble(
  Sampling_time = SAMPLING_TIME_CHOSEN,
  n_contrasts = length(TF_CONTRAST_ORDER),
  n_TF_mapped_genome = nrow(annot_map),
  n_TF_families_mapped = dplyr::n_distinct(
    annot_map$TF_family
  ),
  n_TF_expressed = length(expressed_tf_ids),
  n_TF_DE_any_of_four_contrasts = nrow(plot_tbl),
  n_TF_families_DE = length(
    unique(as.character(plot_tbl$TF_family))
  ),
  n_PlantTFDB_Dca_reference_families = length(
    PLANTTFDB_DCA_FAMILIES
  )
)
readr::write_csv(
  grand_summary,
  file.path(
    dir_qc,
    "TF_GRAND_SUMMARY.csv"
  )
)

saveRDS(
  list(
    main_dir = main_dir,
    annotation_file = TF_ANNOT_FILE,
    annotation_mapped = annot_map,
    metadata = meta2,
    DGE = y,
    logCPM = logcpm,
    all_DE_four_contrasts = all_deg,
    curation = curation,
    selected_TFs = plot_tbl,
    mean_logCPM = mean_logcpm,
    zscore = zscore_mat,
    status_matrix = status_mat,
    global_all_expressed_mean_logCPM = get0(
      "global_all_mean_logcpm",
      ifnotfound = NULL
    ),
    global_all_expressed_zscore = get0(
      "global_all_zscore",
      ifnotfound = NULL
    ),
    global_DE_mean_logCPM = get0(
      "global_de_mean_logcpm",
      ifnotfound = NULL
    ),
    global_DE_zscore = get0(
      "global_de_zscore",
      ifnotfound = NULL
    ),
    response_patterns = response_patterns,
    census = census,
    family_enrichment = family_enrichment,
    family_per_contrast = family_per_contrast,
    family_completeness_audit = family_audit,
    bibliography = bibliography_rules
  ),
  file.path(
    dir_qc,
    "Figure_9B_TF_complete_analysis_objects.rds"
  )
)

capture.output(
  utils::sessionInfo(),
  file = file.path(
    tf_out,
    "sessionInfo.txt"
  )
)

cat("\nDONE Figure 9B COMPLETE TF pipeline\n")
cat("Input run:", main_dir, "\n")
cat("Mapped TFs:", nrow(annot_map), "\n")
cat("Expressed TFs:", length(expressed_tf_ids), "\n")
cat("DE TFs plotted:", nrow(plot_tbl), "\n")
cat(
  "DE TF families plotted:",
  length(unique(as.character(plot_tbl$TF_family))),
  "\n"
)
cat(
  "Expression column order:",
  paste(TF_GROUP_LABELS, collapse = " | "),
  "\n"
)
cat(
  "Global ALL-expressed TF heatmap:",
  if (exists("global_all_zscore")) nrow(global_all_zscore) else 0L,
  "TFs\n"
)
cat(
  "Global DE-TF expression-only heatmap:",
  if (exists("global_de_zscore")) nrow(global_de_zscore) else 0L,
  "TFs\n"
)
cat(
  "Contrast order:",
  paste(TF_CONTRAST_LABELS, collapse = " | "),
  "\n"
)
cat("Output directory:", tf_out, "\n")
