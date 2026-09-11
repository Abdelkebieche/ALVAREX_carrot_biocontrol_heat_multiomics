############################################################
## 01_common_heatmap_helpers.R
## Common loader/helpers for final WGCNA biological heatmaps
##
## PRINCIPLE
##   - Gene selection comes ONLY from the WGCNA biological master.
##   - WGCNA networks were built on D2 + D4.
##   - For the SAME WGCNA-selected genes, expression is recovered at
##     D0 + D2 + D4 from the original count matrix for interpretation.
##   - D0 is PRE-TREATMENT: one baseline per temperature, NO Control/Treated split.
##   - Biological replicates are averaged BEFORE z-scoring.
##   - z-score is calculated gene-by-gene.
##   - For COMBINED PRESTO+ROBILA heatmaps, z-scores are calculated
##     independently within each genotype and then concatenated.
##   - NO family-by-family heatmaps. Detailed families/subfamilies stay in CSV.
##
## OUTPUT PER THEME
##   Temperature/
##   Defence/
##   Temperature_x_Defence/
##
## Each theme contains ONLY 3 heatmaps:
##   1. PRESTO
##   2. ROBILA
##   3. PRESTO + ROBILA combined
# Repository configuration (portable; no HPC-specific paths)
if (!exists('CARROT_CONFIG_LOADED')) source('../../../config/config.R')


############################################################

##############################
## 0. USER SETTINGS
##############################
PROJECT_ROOT <- OUTPUT_ROOT
BIO_ROOT <- BIOLOGICAL_HEATMAP_ROOT

COUNTS_FILE <- RNA_COUNTS_FILE
META_FILE <- RNA_METADATA_FILE

## Optional: set manually to a Biological_Master_... folder.
## Leave NA to use Biological_Heatmaps/LATEST_MASTER_PATH.txt.
MASTER_DIR <- NA_character_

## Same basic count filter as WGCNA script 01.
MIN_TOTAL_COUNTS <- 10
MIN_COUNT        <- 10
MIN_SAMPLES      <- 3

## Heatmap display only. Exact, unclipped z-scores are exported to CSV.
ZSCORE_DISPLAY_LIMIT <- 2.5
SHOW_GENE_NAMES_MAX  <- 120

## Themes to draw. Keep all three for the complete final analysis.
THEMES_TO_RUN <- c("Temperature", "Defence", "Temperature_x_Defence")

## Export resolution
PNG_DPI <- 400

options(stringsAsFactors = FALSE)

##############################
## 1. PACKAGES
##############################
cran_pkgs <- c("dplyr", "tibble", "readr", "stringr", "pheatmap", "svglite")
bioc_pkgs <- c("DESeq2", "SummarizedExperiment")
all_pkgs <- c(cran_pkgs, bioc_pkgs)
missing_pkgs <- all_pkgs[!vapply(all_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) {
  stop("Missing packages: ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
  library(stringr)
  library(pheatmap)
})

##############################
## 2. HELPERS
##############################
clean_sample_names <- function(x) {
  x <- gsub("^X", "", x)
  out <- vapply(x, function(xx) {
    m <- regmatches(xx, regexpr("[0-9]{4}", xx))
    if (length(m) == 1 && nchar(m) == 4) m else xx
  }, character(1))
  unname(out)
}

read_featurecounts_matrix <- function(path) {
  if (!file.exists(path)) stop("Counts file not found: ", path)
  fc <- utils::read.delim(
    path, header = TRUE, sep = "\t", comment.char = "#",
    quote = "", check.names = FALSE, stringsAsFactors = FALSE
  )
  if (!"Geneid" %in% colnames(fc)) stop("Geneid column not found in featureCounts file")
  drop_cols <- intersect(colnames(fc), c("Chr", "Start", "End", "Strand", "Length"))
  mat <- fc |>
    dplyr::select(-dplyr::all_of(drop_cols)) |>
    tibble::column_to_rownames("Geneid") |>
    as.matrix()
  mode(mat) <- "numeric"
  colnames(mat) <- clean_sample_names(colnames(mat))
  if (anyDuplicated(colnames(mat))) stop("Duplicated sample IDs after cleaning featureCounts columns")
  mat
}

normalize_temperature <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- gsub("\\s+|°", "", x)
  x <- gsub("TEMP", "T", x)
  x <- gsub("^0$", "T0", x)
  x <- gsub("^1$", "T1", x)
  x <- gsub("^2$", "T2", x)
  x
}

normalize_treatment <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- gsub("\\s+", "", x)
  dplyr::case_when(
    x %in% c("EAU", "WATER", "H2O", "CTRL", "CONTROL") ~ "Control",
    x %in% c("SDP", "TREATED", "PRI") ~ "Treated",
    TRUE ~ x
  )
}

normalize_time <- function(x) {
  x <- toupper(trimws(as.character(x)))
  dplyr::recode(x, "D0" = "P1", "D2" = "P2", "D4" = "P3", .default = x)
}

TEMP_DISPLAY <- c("T0" = "NH", "T1" = "HS-7", "T2" = "HS-2")
TIME_DISPLAY <- c("P1" = "D0", "P2" = "D2", "P3" = "D4")

as_flag <- function(x) {
  if (is.logical(x)) return(replace(x, is.na(x), FALSE))
  y <- tolower(trimws(as.character(x)))
  !is.na(y) & y %in% c("true", "t", "1", "yes", "y")
}

zscore_rows <- function(mat) {
  mat <- as.matrix(mat)
  storage.mode(mat) <- "numeric"
  mu <- rowMeans(mat, na.rm = TRUE)
  ss <- apply(mat, 1, stats::sd, na.rm = TRUE)
  ss[is.na(ss) | ss == 0] <- 1
  z <- sweep(mat, 1, mu, "-")
  z <- sweep(z, 1, ss, "/")
  z[!is.finite(z)] <- 0
  z
}

write_matrix_csv <- function(mat, path) {
  out <- as.data.frame(mat, check.names = FALSE) |>
    tibble::rownames_to_column("Gene")
  readr::write_csv(out, path)
}

safe_text <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x
}

module_colours <- function(modules) {
  mods <- unique(as.character(modules))
  mods <- mods[!is.na(mods) & nzchar(mods)]
  if (!length(mods)) return(character())
  vals <- vapply(mods, function(m) {
    ok <- tryCatch({ grDevices::col2rgb(m); TRUE }, error = function(e) FALSE)
    if (ok) m else "grey70"
  }, character(1))
  setNames(vals, mods)
}

save_pheatmap_all <- function(ph, base, width, height) {
  grDevices::pdf(paste0(base, ".pdf"), width = width, height = height, useDingbats = FALSE)
  grid::grid.newpage(); grid::grid.draw(ph$gtable); grDevices::dev.off()

  grDevices::png(paste0(base, ".png"), width = width, height = height,
                 units = "in", res = PNG_DPI, bg = "white")
  grid::grid.newpage(); grid::grid.draw(ph$gtable); grDevices::dev.off()

  svglite::svglite(paste0(base, ".svg"), width = width, height = height, bg = "white")
  grid::grid.newpage(); grid::grid.draw(ph$gtable); grDevices::dev.off()
}

##############################
## 3. LOCATE MASTER + OUTPUT
##############################
if (is.na(MASTER_DIR)) {
  latest_file <- file.path(BIO_ROOT, "LATEST_MASTER_PATH.txt")
  if (file.exists(latest_file)) {
    MASTER_DIR <- trimws(readLines(latest_file, warn = FALSE)[1])
  } else {
    candidates <- list.dirs(BIO_ROOT, recursive = FALSE, full.names = TRUE)
    candidates <- candidates[grepl("Biological_Master_", basename(candidates), fixed = TRUE)]
    if (!length(candidates)) stop("No Biological_Master_* directory found in: ", BIO_ROOT)
    MASTER_DIR <- candidates[order(file.info(candidates)$mtime, decreasing = TRUE)][1]
  }
}
if (!dir.exists(MASTER_DIR)) stop("MASTER_DIR does not exist: ", MASTER_DIR)

MASTER_FILE <- file.path(MASTER_DIR, "01_tables", "02_Master_gene_annotation_WGCNA_maSigPro.csv")
if (!file.exists(MASTER_FILE)) stop("Master annotation file not found: ", MASTER_FILE)

OUT_ROOT <- file.path(MASTER_DIR, "05_WGCNA_D0D2D4_Zscore_FINAL")
DIR_EXPR <- file.path(OUT_ROOT, "00_expression_mean_D0D2D4")
dir.create(DIR_EXPR, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(OUT_ROOT, "16_WGCNA_D0D2D4_heatmaps.log")
log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", paste0(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

log_msg("MASTER: ", MASTER_DIR)

##############################
## 4. LOAD WGCNA BIOLOGICAL MASTER
##############################
MASTER <- readr::read_csv(MASTER_FILE, show_col_types = FALSE, progress = FALSE) |>
  dplyr::mutate(
    Gene = as.character(Gene),
    Genotype = as.character(Genotype),
    In_WGCNA = as_flag(In_WGCNA),
    Candidate_Temperature = as_flag(Candidate_Temperature),
    Candidate_Defence = as_flag(Candidate_Defence),
    Candidate_Temp_x_Defence = as_flag(Candidate_Temp_x_Defence),
    Is_HSP = as_flag(Is_HSP),
    Is_TF = as_flag(Is_TF),
    Is_Immune = as_flag(Is_Immune),
    Is_FunctionalSet = as_flag(Is_FunctionalSet)
  ) |>
  dplyr::filter(In_WGCNA, Genotype %in% c("PRESTO", "ROBILA")) |>
  dplyr::distinct(Genotype, Gene, .keep_all = TRUE)

THEMES <- list(
  Temperature = list(
    flag = "Candidate_Temperature",
    title = "Temperature-associated WGCNA genes"
  ),
  Defence = list(
    flag = "Candidate_Defence",
    title = "Defence-associated WGCNA genes"
  ),
  Temperature_x_Defence = list(
    flag = "Candidate_Temp_x_Defence",
    title = "Temperature × defence WGCNA genes"
  )
)

all_selected <- unique(MASTER$Gene[
  MASTER$Candidate_Temperature |
    MASTER$Candidate_Defence |
    MASTER$Candidate_Temp_x_Defence
])
if (!length(all_selected)) stop("No selected WGCNA biological genes found in master table")
log_msg("Unique WGCNA-selected genes across themes: ", length(all_selected))

##############################
## 5. REBUILD D0+D2+D4 VST FOR VISUALISATION
##    FROM ORIGINAL COUNTS — NO maSigPro FILTERING
##############################
count_mat <- read_featurecounts_matrix(COUNTS_FILE)
meta_raw <- utils::read.csv(META_FILE, sep = ";", stringsAsFactors = FALSE)

required_meta <- c("Sample", "Genotype", "Temperature", "Treatment", "Sampling_time")
miss_meta <- setdiff(required_meta, colnames(meta_raw))
if (length(miss_meta)) stop("Missing metadata columns: ", paste(miss_meta, collapse = ", "))

meta_all <- meta_raw |>
  dplyr::mutate(
    Sample = sprintf("%04d", as.integer(as.character(Sample))),
    Genotype = trimws(as.character(Genotype)),
    Temperature = normalize_temperature(Temperature),
    Treatment = normalize_treatment(Treatment),
    Sampling_time = normalize_time(Sampling_time)
  ) |>
  dplyr::filter(
    Genotype %in% c("PRESTO", "ROBILA"),
    Sampling_time %in% c("P1", "P2", "P3"),
    Temperature %in% c("T0", "T1", "T2"),
    Sample %in% colnames(count_mat)
  ) |>
  ## D0 is pre-treatment. Any D0 row labelled Treated is not a real
  ## post-treatment condition and is excluded from visualisation.
  dplyr::filter(!(Sampling_time == "P1" & Treatment == "Treated"))

CONDITION_ORDER <- c(
  "NH_D0", "HS-7_D0", "HS-2_D0",
  "NH_Control_D2", "NH_Treated_D2",
  "HS-7_Control_D2", "HS-7_Treated_D2",
  "HS-2_Control_D2", "HS-2_Treated_D2",
  "NH_Control_D4", "NH_Treated_D4",
  "HS-7_Control_D4", "HS-7_Treated_D4",
  "HS-2_Control_D4", "HS-2_Treated_D4"
)

build_mean_vst <- function(genotype, selected_genes) {
  meta <- meta_all |>
    dplyr::filter(Genotype == genotype) |>
    dplyr::arrange(match(Sampling_time, c("P1", "P2", "P3")), Temperature, Treatment)

  if (!nrow(meta)) stop("No metadata rows for genotype: ", genotype)
  counts <- count_mat[, meta$Sample, drop = FALSE]
  if (!identical(colnames(counts), meta$Sample)) stop("Count/metadata alignment failed for ", genotype)

  keep_standard <-
    rowSums(counts, na.rm = TRUE) >= MIN_TOTAL_COUNTS &
    rowSums(counts >= MIN_COUNT, na.rm = TRUE) >= MIN_SAMPLES

  ## Force WGCNA-selected genes into the visualisation transform if they
  ## have any counts, even if the all-time-point filter changes slightly.
  keep_selected <- rownames(counts) %in% selected_genes & rowSums(counts, na.rm = TRUE) > 0
  counts_f <- counts[keep_standard | keep_selected, , drop = FALSE]

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = round(counts_f),
    colData = meta,
    design = ~ 1
  )
  dds <- DESeq2::estimateSizeFactors(dds)
  vsd <- DESeq2::vst(dds, blind = TRUE)
  expr <- SummarizedExperiment::assay(vsd)

  temp_disp <- unname(TEMP_DISPLAY[meta$Temperature])
  time_disp <- unname(TIME_DISPLAY[meta$Sampling_time])

  condition <- ifelse(
    meta$Sampling_time == "P1",
    paste(temp_disp, "D0", sep = "_"),
    paste(temp_disp, meta$Treatment, time_disp, sep = "_")
  )

  observed <- unique(condition)
  missing_conditions <- setdiff(CONDITION_ORDER, observed)
  if (length(missing_conditions)) {
    warning(genotype, " missing expected conditions: ", paste(missing_conditions, collapse = ", "))
  }

  means <- lapply(CONDITION_ORDER, function(cond) {
    samps <- meta$Sample[condition == cond]
    samps <- intersect(samps, colnames(expr))
    if (!length(samps)) return(rep(NA_real_, nrow(expr)))
    rowMeans(expr[, samps, drop = FALSE], na.rm = TRUE)
  })
  means <- do.call(cbind, means)
  rownames(means) <- rownames(expr)
  colnames(means) <- CONDITION_ORDER
  means
}

EXPR_PRESTO <- build_mean_vst("PRESTO", all_selected)
EXPR_ROBILA <- build_mean_vst("ROBILA", all_selected)

write_matrix_csv(EXPR_PRESTO, file.path(DIR_EXPR, "Mean_VST_D0_D2_D4_PRESTO.csv"))
write_matrix_csv(EXPR_ROBILA, file.path(DIR_EXPR, "Mean_VST_D0_D2_D4_ROBILA.csv"))

##############################
## 6. COLUMN ANNOTATIONS
##############################
make_condition_annotation <- function(genotype = NULL) {
  parse_one <- function(x) {
    parts <- strsplit(x, "_", fixed = TRUE)[[1]]
    if (length(parts) == 2 && parts[2] == "D0") {
      data.frame(
        Temperature = parts[1], Treatment = "Baseline", Time = "D0",
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(
        Temperature = parts[1], Treatment = parts[2], Time = parts[3],
        stringsAsFactors = FALSE
      )
    }
  }
  ann <- dplyr::bind_rows(lapply(CONDITION_ORDER, parse_one))
  rownames(ann) <- CONDITION_ORDER
  if (!is.null(genotype)) ann$Genotype <- genotype
  ann
}

COL_ANN_SINGLE <- make_condition_annotation()

TEMP_COLORS <- c("NH" = "#3B6FB6", "HS-7" = "#E69F00", "HS-2" = "#D55E00")
TRT_COLORS  <- c("Baseline" = "#9A9A9A", "Control" = "#56B4E9", "Treated" = "#009E73")
TIME_COLORS <- c("D0" = "#8C8C8C", "D2" = "#0072B2", "D4" = "#CC79A7")
GENO_COLORS <- c("PRESTO" = "#0072B2", "ROBILA" = "#D55E00")
SELECT_COLORS <- c("Both" = "#009E73", "PRESTO only" = "#0072B2", "ROBILA only" = "#D55E00")

HEAT_COLORS <- grDevices::colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(101)
HEAT_BREAKS <- seq(-ZSCORE_DISPLAY_LIMIT, ZSCORE_DISPLAY_LIMIT, length.out = 102)

##############################
## 7. HEATMAP WRITER
##############################
plot_single_genotype <- function(theme_name, theme_title, genotype, genes, master_sub, expr_all, out_dir) {
  genes <- unique(genes)
  genes <- intersect(genes, rownames(expr_all))
  if (length(genes) < 2) {
    log_msg(theme_name, " | ", genotype, ": fewer than 2 genes; heatmap skipped")
    return(invisible(NULL))
  }

  raw <- expr_all[genes, CONDITION_ORDER, drop = FALSE]
  keep <- rowSums(is.finite(raw)) >= 2
  raw <- raw[keep, , drop = FALSE]
  if (nrow(raw) < 2) return(invisible(NULL))

  z_exact <- zscore_rows(raw)
  z_plot <- pmax(pmin(z_exact, ZSCORE_DISPLAY_LIMIT), -ZSCORE_DISPLAY_LIMIT)

  ann_rows <- master_sub |>
    dplyr::filter(Genotype == genotype, Gene %in% rownames(raw)) |>
    dplyr::distinct(Gene, .keep_all = TRUE) |>
    dplyr::select(Gene, ModuleColor) |>
    as.data.frame()
  rownames(ann_rows) <- ann_rows$Gene
  ann_rows$Gene <- NULL
  ann_rows <- ann_rows[rownames(raw), , drop = FALSE]
  colnames(ann_rows) <- "WGCNA module"

  mod_pal <- module_colours(ann_rows[[1]])
  ann_colors <- list(
    Temperature = TEMP_COLORS,
    Treatment = TRT_COLORS,
    Time = TIME_COLORS,
    `WGCNA module` = mod_pal
  )

  show_names <- nrow(raw) <= SHOW_GENE_NAMES_MAX
  h <- max(6, min(18, 4 + nrow(raw) * 0.10))

  ph <- pheatmap::pheatmap(
    z_plot,
    color = HEAT_COLORS,
    breaks = HEAT_BREAKS,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    annotation_col = COL_ANN_SINGLE,
    annotation_row = ann_rows,
    annotation_colors = ann_colors,
    gaps_col = c(3, 9),
    show_rownames = show_names,
    show_colnames = TRUE,
    fontsize = 9,
    fontsize_row = if (show_names) 7 else 1,
    fontsize_col = 9,
    border_color = NA,
    main = paste0(theme_title, " — ", genotype, "\nWGCNA-selected genes; mean biological replicates; D0 added for interpretation"),
    silent = TRUE
  )

  base <- file.path(out_dir, paste0(theme_name, "_", genotype, "_Zscore"))
  save_pheatmap_all(ph, base, width = 12.5, height = h)
  write_matrix_csv(raw, paste0(base, "_MeanVST.csv"))
  write_matrix_csv(z_exact, paste0(base, "_Zscore_exact.csv"))
  log_msg(theme_name, " | ", genotype, ": ", nrow(raw), " genes")
  invisible(list(raw = raw, z = z_exact))
}

plot_combined <- function(theme_name, theme_title, genes_p, genes_r, master_sub, expr_p, expr_r, out_dir) {
  genes_union <- sort(unique(c(genes_p, genes_r)))
  genes_union <- intersect(genes_union, intersect(rownames(expr_p), rownames(expr_r)))
  if (length(genes_union) < 2) {
    log_msg(theme_name, " | COMBINED: fewer than 2 genes; heatmap skipped")
    return(invisible(NULL))
  }

  raw_p <- expr_p[genes_union, CONDITION_ORDER, drop = FALSE]
  raw_r <- expr_r[genes_union, CONDITION_ORDER, drop = FALSE]

  keep <- rowSums(is.finite(raw_p)) >= 2 & rowSums(is.finite(raw_r)) >= 2
  raw_p <- raw_p[keep, , drop = FALSE]
  raw_r <- raw_r[keep, , drop = FALSE]
  genes_union <- rownames(raw_p)
  if (length(genes_union) < 2) return(invisible(NULL))

  ## KEY RULE: z-score separately within genotype, then concatenate.
  z_p <- zscore_rows(raw_p)
  z_r <- zscore_rows(raw_r)
  colnames(z_p) <- paste0("PRESTO | ", CONDITION_ORDER)
  colnames(z_r) <- paste0("ROBILA | ", CONDITION_ORDER)
  z_exact <- cbind(z_p, z_r)
  z_plot <- pmax(pmin(z_exact, ZSCORE_DISPLAY_LIMIT), -ZSCORE_DISPLAY_LIMIT)

  raw_comb <- cbind(raw_p, raw_r)
  colnames(raw_comb) <- c(paste0("PRESTO | ", CONDITION_ORDER), paste0("ROBILA | ", CONDITION_ORDER))

  ann_p <- make_condition_annotation("PRESTO")
  ann_r <- make_condition_annotation("ROBILA")
  rownames(ann_p) <- paste0("PRESTO | ", CONDITION_ORDER)
  rownames(ann_r) <- paste0("ROBILA | ", CONDITION_ORDER)
  col_ann <- rbind(ann_p, ann_r)
  col_ann <- col_ann[colnames(z_exact), , drop = FALSE]

  mp <- master_sub |>
    dplyr::filter(Genotype == "PRESTO", Gene %in% genes_union) |>
    dplyr::distinct(Gene, .keep_all = TRUE) |>
    dplyr::select(Gene, PRESTO_module = ModuleColor)
  mr <- master_sub |>
    dplyr::filter(Genotype == "ROBILA", Gene %in% genes_union) |>
    dplyr::distinct(Gene, .keep_all = TRUE) |>
    dplyr::select(Gene, ROBILA_module = ModuleColor)

  row_ann <- data.frame(Gene = genes_union, stringsAsFactors = FALSE) |>
    dplyr::left_join(mp, by = "Gene") |>
    dplyr::left_join(mr, by = "Gene") |>
    dplyr::mutate(
      Selected_in = dplyr::case_when(
        Gene %in% genes_p & Gene %in% genes_r ~ "Both",
        Gene %in% genes_p ~ "PRESTO only",
        Gene %in% genes_r ~ "ROBILA only",
        TRUE ~ NA_character_
      )
    ) |>
    as.data.frame()
  rownames(row_ann) <- row_ann$Gene
  row_ann$Gene <- NULL
  row_ann <- row_ann[genes_union, , drop = FALSE]

  mods <- unique(c(safe_text(row_ann$PRESTO_module), safe_text(row_ann$ROBILA_module)))
  mods <- mods[nzchar(mods)]
  mod_pal <- module_colours(mods)

  ann_colors <- list(
    Temperature = TEMP_COLORS,
    Treatment = TRT_COLORS,
    Time = TIME_COLORS,
    Genotype = GENO_COLORS,
    PRESTO_module = mod_pal,
    ROBILA_module = mod_pal,
    Selected_in = SELECT_COLORS
  )

  show_names <- nrow(z_exact) <= SHOW_GENE_NAMES_MAX
  h <- max(6, min(18, 4 + nrow(z_exact) * 0.10))

  ph <- pheatmap::pheatmap(
    z_plot,
    color = HEAT_COLORS,
    breaks = HEAT_BREAKS,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    annotation_col = col_ann,
    annotation_row = row_ann,
    annotation_colors = ann_colors,
    gaps_col = c(3, 9, 15, 18, 24),
    show_rownames = show_names,
    show_colnames = TRUE,
    fontsize = 8,
    fontsize_row = if (show_names) 7 else 1,
    fontsize_col = 7.5,
    border_color = NA,
    main = paste0(theme_title, " — PRESTO + ROBILA\n",
                  "WGCNA selection; mean replicates; z-score calculated independently within each genotype"),
    silent = TRUE
  )

  base <- file.path(out_dir, paste0(theme_name, "_PRESTO_ROBILA_COMBINED_Zscore"))
  save_pheatmap_all(ph, base, width = 20, height = h)
  write_matrix_csv(raw_comb, paste0(base, "_MeanVST.csv"))
  write_matrix_csv(z_exact, paste0(base, "_Zscore_exact.csv"))
  log_msg(theme_name, " | COMBINED: ", nrow(z_exact), " union genes")
  invisible(list(raw = raw_comb, z = z_exact))
}


##############################
## 8. GENERIC WGCNA HEATMAP SET RUNNER
##############################
run_heatmap_set <- function(set_name, set_title, selected, out_root = NULL) {
  selected <- selected |>
    dplyr::filter(In_WGCNA, Genotype %in% c("PRESTO", "ROBILA")) |>
    dplyr::distinct(Genotype, Gene, .keep_all = TRUE)

  if (is.null(out_root)) out_root <- file.path(MASTER_DIR, "05_WGCNA_D0D2D4_Zscore_FINAL")
  out_dir <- file.path(out_root, set_name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  if (!nrow(selected)) {
    warning("No genes selected for set: ", set_name)
    readr::write_csv(data.frame(Set=set_name, N_PRESTO=0, N_ROBILA=0, N_union=0, N_shared=0),
                     file.path(out_dir, paste0(set_name, "_summary.csv")))
    return(invisible(NULL))
  }

  genes_p <- unique(selected$Gene[selected$Genotype == "PRESTO"])
  genes_r <- unique(selected$Gene[selected$Genotype == "ROBILA"])

  ## Keep detailed biological information in CSV only; no family-by-family heatmaps.
  keep_cols <- intersect(
    c(
      "Genotype", "Gene", "ModuleColor", "Statistical_module_class",
      "FDR_Temperature", "FDR_Treatment", "FDR_Time",
      "FDR_Temp_x_Treatment", "FDR_Temp_x_Time",
      "FDR_Treatment_x_Time", "FDR_Temp_x_Treatment_x_Time",
      "Temperature_recruited", "Treatment_recruited", "Time_recruited",
      "Temp_x_Treatment_recruited", "Temp_x_Time_recruited",
      "Treatment_x_Time_recruited", "Triple_interaction_recruited",
      "kME", "GS_Temperature", "GS_Treatment", "GS_Time",
      "IsHub", "IsHighGS_Temp", "IsHighGS_Trt",
      "Annotation_classes",
      "Is_HSP", "HSP_family", "HSP_subfamily", "func_clean",
      "Is_TF", "TF_family", "TF_subgroup", "TF_function",
      "Is_Immune", "Immune_module", "Immune_submodule", "Immune_evidence", "Immune_annotation",
      "Is_FunctionalSet", "Functional_families", "Functional_subfamilies",
      "Functional_KEGG_pathways", "Functional_annotations",
      "Candidate_Temperature", "Candidate_Defence", "Candidate_Temp_x_Defence",
      "Heatmap_primary_theme"
    ), colnames(selected)
  )

  readr::write_csv(
    selected |> dplyr::select(dplyr::all_of(keep_cols)) |> dplyr::arrange(Genotype, ModuleColor, Gene),
    file.path(out_dir, paste0(set_name, "_WGCNA_selected_genes_annotations.csv"))
  )

  module_summary <- selected |>
    dplyr::count(Genotype, ModuleColor, Statistical_module_class, name = "N_genes") |>
    dplyr::arrange(Genotype, dplyr::desc(N_genes))
  readr::write_csv(module_summary, file.path(out_dir, paste0(set_name, "_module_summary.csv")))

  ## Family/subfamily information remains as CSV summaries.
  if ("HSP_family" %in% names(selected)) {
    readr::write_csv(selected |> dplyr::filter(Is_HSP) |> dplyr::count(Genotype, HSP_family, HSP_subfamily, name="N_genes"),
                     file.path(out_dir, paste0(set_name, "_HSP_family_summary.csv")))
  }
  if ("TF_family" %in% names(selected)) {
    readr::write_csv(selected |> dplyr::filter(Is_TF) |> dplyr::count(Genotype, TF_family, TF_subgroup, name="N_genes"),
                     file.path(out_dir, paste0(set_name, "_TF_family_summary.csv")))
  }
  if ("Immune_module" %in% names(selected)) {
    readr::write_csv(selected |> dplyr::filter(Is_Immune) |> dplyr::count(Genotype, Immune_module, Immune_submodule, name="N_genes"),
                     file.path(out_dir, paste0(set_name, "_Immune_family_summary.csv")))
  }
  if ("Functional_families" %in% names(selected)) {
    readr::write_csv(selected |> dplyr::filter(Is_FunctionalSet) |>
                       dplyr::select(Genotype, Gene, Functional_families, Functional_subfamilies, Functional_KEGG_pathways),
                     file.path(out_dir, paste0(set_name, "_Functional_annotations.csv")))
  }

  ## Exactly 3 heatmaps: PRESTO, ROBILA, COMBINED.
  plot_single_genotype(set_name, set_title, "PRESTO", genes_p, selected, EXPR_PRESTO, out_dir)
  plot_single_genotype(set_name, set_title, "ROBILA", genes_r, selected, EXPR_ROBILA, out_dir)
  plot_combined(set_name, set_title, genes_p, genes_r, selected, EXPR_PRESTO, EXPR_ROBILA, out_dir)

  sm <- data.frame(
    Set = set_name,
    N_PRESTO = length(genes_p),
    N_ROBILA = length(genes_r),
    N_union = length(unique(c(genes_p, genes_r))),
    N_shared = length(intersect(genes_p, genes_r))
  )
  readr::write_csv(sm, file.path(out_dir, paste0(set_name, "_summary.csv")))
  invisible(sm)
}

cat("Common WGCNA D0-D2-D4 heatmap environment loaded.\n")
cat("MASTER_DIR: ", MASTER_DIR, "\n", sep="")
cat("D0 has 3 baseline columns only; D2/D4 have Control/Treated.\n")
