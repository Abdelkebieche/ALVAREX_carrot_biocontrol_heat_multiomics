############################################################
## 06_LCMS_FIG5_FIG6_ANNOTATED_HUBS_FINAL_WORKING.R
## POWERPOINT SVG REVISION — 2026-08-20
## SVG export now uses svglite::svglite(fix_text_size = FALSE)
## so text and graphical objects remain as simple SVG elements and
## are substantially easier to convert/ungroup/edit in PowerPoint.
##
## FIGURE 5
##   Heatmaps of manually structure-named historical hub features
##   from the technician's Excel workbooks:
##     historical Blue / Green / Purple / Magenta modules.
##
## FIGURE 6
##   Projection of the four HISTORICAL hub-feature programmes
##   across PRESTO + ROBILA and D0/D4/D10.
##
## IMPORTANT DESIGN DECISIONS
##   1) Historical module names are kept as the biological reference.
##   2) New WGCNA colours are used ONLY to verify concordance.
##   3) Figure 5 uses structure-named rows only (not generic "fgt",
##      "flavonoide" or "isotope" rows), avoiding over-representation
##      of multiple fragments/adducts from the same metabolite.
##   4) Figure 6 uses ALL features contained in each historical hub
##      workbook, because it represents the historical hub programme,
##      not only the chemically identified subset.
##   5) LC-MS abundance = log2(peak area + 1), then feature-wise Z-score.
##   6) Z-score is calculated ACROSS THE 120 BIOLOGICAL SAMPLES first;
##      biological replicates are retained for trajectory means/SE.
##
## This script follows the SAME metadata/sample-ID logic as the
## previously working LC-MS PLS-DA and WGCNA scripts:
##   - metadata separator = ";"
##   - Sample IDs are zero-padded to 4 digits (485 -> 0485)
##   - LC-MS columns provide the same 4-digit ID after "2025-F-"
# Repository configuration (portable; no HPC-specific paths)
if (!exists('CARROT_CONFIG_LOADED')) source('../../../config/config.R')


############################################################

##############################
## 0. PATHS
##############################

META_FILE <- LCMS_METADATA_FILE
LCMS_FILE <- LCMS_FILE

## Exact successful WGCNA run used for the technician-annotated hub files.
## Do NOT append the run name again: this path already IS the run directory.
WGCNA_DIR <- LCMS_REFERENCE_RUN

OUT_DIR <- file.path(LCMS_FIGURE_ROOT, "annotated_hubs_projection")

##############################
## 1. DISPLAY / BIOLOGICAL MAPPING
##############################

TIME_DISPLAY <- c("P1" = "D0", "P3" = "D4", "P4" = "D10")
TEMP_DISPLAY <- c("T0" = "NH", "T1" = "HS-7", "T2" = "HS-2")
TRT_DISPLAY  <- c("EAU" = "Control", "SDP" = "Treated")

GENO_LEVELS <- c("PRESTO", "ROBILA")
TIME_LEVELS <- c("D0", "D4", "D10")
TEMP_LEVELS <- c("NH", "HS-7", "HS-2")
TRT_LEVELS  <- c("Control", "Treated")
HIST_MODULE_LEVELS <- c("Blue", "Green", "Purple", "Magenta")

TEMP_COLORS <- c(
  "NH"   = "#313695",
  "HS-7" = "#74ADD1",
  "HS-2" = "#F46D43"
)

## Heatmap palette only; independent from experimental-factor colours.
HEAT_LOW  <- "#313695"
HEAT_MID  <- "#FFFFFF"
HEAT_HIGH <- "#A50026"

##############################
## 2. PACKAGES
##############################

required_pkgs <- c(
  "data.table", "dplyr", "tidyr", "readxl", "ggplot2", "scales",
  "svglite"
)
missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]
if (length(missing_pkgs)) {
  stop(
    "Missing R packages: ", paste(missing_pkgs, collapse = ", "),
    "\nInstall them, then rerun the script."
  )
}

##############################
## 3. HELPERS
##############################

msg <- function(...) message(sprintf(...))

zscore_rows <- function(mat) {
  z <- t(scale(t(mat)))
  z[!is.finite(z)] <- 0
  z
}

se_fun <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) <= 1L) return(NA_real_)
  stats::sd(x) / sqrt(length(x))
}

save_plot <- function(p, basename, width, height, out_dir) {
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  svg_file <- file.path(out_dir, paste0(basename, ".svg"))
  pdf_file <- file.path(out_dir, paste0(basename, ".pdf"))
  png_file <- file.path(out_dir, paste0(basename, ".png"))
  
  ############################################################
  ## POWERPOINT-EDITABLE SVG
  ##
  ## IMPORTANT:
  ##   Do NOT use grDevices::svg() here.
  ##   The Cairo-based R SVG device can encode text/glyphs in a way
  ##   that PowerPoint may display correctly but cannot reliably
  ##   convert/ungroup into editable Office shapes.
  ##
  ## svglite writes a cleaner SVG and keeps labels as SVG <text>
  ## elements.  fix_text_size = FALSE deliberately avoids
  ## textLength/lengthAdjust attributes, making subsequent editing
  ## in PowerPoint / Illustrator / Inkscape simpler.
  ############################################################
  
  svglite::svglite(
    filename = svg_file,
    width = width,
    height = height,
    bg = "white",
    standalone = TRUE,
    fix_text_size = FALSE
  )
  print(p)
  grDevices::dev.off()
  
  ## Basic SVG self-check: an editable ggplot SVG should contain
  ## ordinary text nodes rather than only Cairo glyph definitions.
  svg_txt <- paste(readLines(svg_file, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  if (!grepl("<text", svg_txt, fixed = TRUE)) {
    warning(
      "SVG export contains no <text> elements: ", svg_file,
      "\nText may not remain directly editable in PowerPoint."
    )
  }
  if (grepl("<symbol", svg_txt, fixed = TRUE) &&
      grepl("<use", svg_txt, fixed = TRUE)) {
    warning(
      "SVG contains <symbol>/<use> glyph references: ", svg_file,
      "\nThis is typical of Cairo SVG and can reduce PowerPoint editability."
    )
  }
  
  ############################################################
  ## PDF — unchanged publication/vector backup
  ############################################################
  grDevices::pdf(
    file = pdf_file,
    width = width,
    height = height,
    useDingbats = FALSE
  )
  print(p)
  grDevices::dev.off()
  
  ############################################################
  ## PNG — unchanged high-resolution raster backup
  ############################################################
  grDevices::png(
    filename = png_file,
    width = width,
    height = height,
    units = "in",
    res = 400,
    type = if (.Platform$OS.type == "windows") "windows" else "cairo"
  )
  print(p)
  grDevices::dev.off()
  
  message("Saved PowerPoint-editable SVG: ", svg_file)
  message("Saved PDF: ", pdf_file)
  message("Saved PNG: ", png_file)
}

theme_pub <- function(base_size = 11) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 2),
      plot.subtitle = ggplot2::element_text(size = base_size),
      strip.text = ggplot2::element_text(face = "bold"),
      strip.background = ggplot2::element_rect(fill = "white", colour = "black"),
      legend.position = "top",
      plot.margin = ggplot2::margin(8, 12, 8, 8)
    )
}

find_one_file <- function(directory, pattern) {
  hits <- list.files(
    directory,
    pattern = pattern,
    full.names = TRUE,
    ignore.case = TRUE
  )
  if (length(hits) != 1L) {
    stop(
      "Expected exactly one file matching pattern:\n", pattern,
      "\nin:\n", directory,
      "\nFound: ", length(hits),
      if (length(hits)) paste0("\n", paste(basename(hits), collapse = "\n")) else ""
    )
  }
  hits
}

clean_text <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "NaN", "<NA>")] <- NA_character_
  x
}

##############################
## 4. RESOLVE INPUT DIRECTORIES + OUTPUTS
##############################

if (!file.exists(META_FILE)) stop("META_FILE not found:\n", META_FILE)
if (!file.exists(LCMS_FILE)) stop("LCMS_FILE not found:\n", LCMS_FILE)

if (!dir.exists(WGCNA_DIR)) {
  stop(
    "WGCNA run directory not found:\n", WGCNA_DIR,
    "\n\nThis path must directly contain 04_hub_features, 03_module_trait, 02_network, etc."
  )
}

HUB_DIR <- file.path(WGCNA_DIR, "04_hub_features", "by_module")
if (!dir.exists(HUB_DIR)) {
  stop(
    "Hub directory not found:\n", HUB_DIR,
    "\nExpected the technician XLSX files and Module_*_central_features.csv here."
  )
}

DIR_TABLES <- file.path(OUT_DIR, "00_TABLES")
DIR_FIG5   <- file.path(OUT_DIR, "01_FIGURE5_HEATMAPS")
DIR_FIG6   <- file.path(OUT_DIR, "02_FIGURE6_TRAJECTORIES")
DIR_QC     <- file.path(OUT_DIR, "03_QC")
invisible(lapply(c(OUT_DIR, DIR_TABLES, DIR_FIG5, DIR_FIG6, DIR_QC), dir.create,
                 recursive = TRUE, showWarnings = FALSE))

msg("WGCNA run: %s", WGCNA_DIR)
msg("Hub directory: %s", HUB_DIR)

## Preflight: these 8 files are required before any LC-MS calculation starts.
required_hub_patterns <- c(
  "^Hub_metabolites_blue.*\\.xlsx$",
  "^Hub_metabolites_green.*\\.xlsx$",
  "^Hub_metabolites_purple.*\\.xlsx$",
  "^Hub_metabolites_magenta.*\\.xlsx$",
  "^Module_blue_central_features\\.csv$",
  "^Module_green_central_features\\.csv$",
  "^Module_greenyellow_central_features\\.csv$",
  "^Module_purple_central_features\\.csv$"
)
for (pat in required_hub_patterns) {
  hits <- list.files(HUB_DIR, pattern = pat, ignore.case = TRUE)
  if (length(hits) != 1L) {
    stop(
      "Preflight failed for HUB_DIR. Pattern: ", pat,
      "\nFound ", length(hits), " matching file(s) in:\n", HUB_DIR,
      "\nCurrent files:\n", paste(list.files(HUB_DIR), collapse = "\n")
    )
  }
}
msg("WGCNA/hub preflight OK.")

##############################
## 5. METADATA
##    EXACT LOGIC USED IN WORKING LC-MS/WGCNA SCRIPTS
##############################

meta <- data.table::fread(
  META_FILE,
  sep = ";",
  header = TRUE,
  data.table = FALSE,
  colClasses = "character",
  encoding = "UTF-8",
  showProgress = FALSE
)

names(meta) <- trimws(sub("^\\ufeff", "", names(meta)))

## Limited aliases, only as safeguards. Your actual file already uses
## the canonical names below.
rename_if_absent <- function(df, old, new) {
  if (!(new %in% names(df)) && old %in% names(df)) {
    names(df)[names(df) == old] <- new
  }
  df
}
meta <- rename_if_absent(meta, "Replicate", "Replicat")
meta <- rename_if_absent(meta, "Sampling", "Sampling_time")
meta <- rename_if_absent(meta, "Sampling time", "Sampling_time")

required_meta <- c(
  "Sample", "Genotype", "Temperature", "Treatment", "Replicat", "Sampling_time"
)
missing_meta <- setdiff(required_meta, names(meta))
if (length(missing_meta)) {
  stop(
    "Missing metadata columns: ", paste(missing_meta, collapse = ", "),
    "\nColumns detected: ", paste(names(meta), collapse = " | ")
  )
}

## CRITICAL: preserve the same 4-digit convention used by every working
## LC-MS/WGCNA script: 485 -> 0485, 996 -> 0996, 1021 -> 1021.
meta <- meta |>
  dplyr::mutate(
    Sample_num = suppressWarnings(as.integer(Sample)),
    Sample = sprintf("%04d", Sample_num),
    Genotype = toupper(trimws(Genotype)),
    Temperature = toupper(trimws(Temperature)),
    Treatment = toupper(trimws(Treatment)),
    Replicat = trimws(Replicat),
    Sampling_time = toupper(trimws(Sampling_time)),
    Sampling_D = unname(TIME_DISPLAY[Sampling_time]),
    Temperature_D = unname(TEMP_DISPLAY[Temperature]),
    Treatment_D = unname(TRT_DISPLAY[Treatment])
  )

if (anyNA(meta$Sample_num)) stop("Non-numeric Sample value found in metadata.")
if (nrow(meta) != 120L) stop("Expected 120 metadata rows; found ", nrow(meta), ".")
if (anyDuplicated(meta$Sample)) stop("Duplicated Sample IDs after 4-digit formatting.")

if (anyNA(meta$Sampling_D)) {
  stop("Unknown Sampling_time values: ", paste(unique(meta$Sampling_time[is.na(meta$Sampling_D)]), collapse = ", "))
}
if (anyNA(meta$Temperature_D)) {
  stop("Unknown Temperature values: ", paste(unique(meta$Temperature[is.na(meta$Temperature_D)]), collapse = ", "))
}
if (anyNA(meta$Treatment_D)) {
  stop("Unknown Treatment values: ", paste(unique(meta$Treatment[is.na(meta$Treatment_D)]), collapse = ", "))
}

if (any(meta$Sampling_time == "P1" & meta$Treatment != "EAU")) {
  stop("P1/D0 contains a non-EAU sample; design mismatch.")
}

meta <- meta |>
  dplyr::mutate(
    Genotype = factor(Genotype, levels = GENO_LEVELS),
    Sampling_D = factor(Sampling_D, levels = TIME_LEVELS),
    Temperature_D = factor(Temperature_D, levels = TEMP_LEVELS),
    Treatment_D = factor(Treatment_D, levels = TRT_LEVELS)
  )

msg("Metadata OK: %d rows", nrow(meta))
msg("Metadata sample IDs: %s ... %s", meta$Sample[1], meta$Sample[nrow(meta)])

utils::write.csv(meta, file.path(DIR_QC, "metadata_clean_120.csv"), row.names = FALSE)

##############################
## 6. FULL LC-MS INPUT + EXACT 120-SAMPLE MATCHING
##############################

msg("Reading LC-MS table...")
lcms <- data.table::fread(
  LCMS_FILE,
  sep = ",",
  header = TRUE,
  data.table = FALSE,
  check.names = FALSE,
  encoding = "Latin-1",
  showProgress = TRUE
)

if (nrow(lcms) != 17846L) {
  warning("Expected 17,846 LC-MS features; found ", nrow(lcms), ".")
}

## Exact format used in the canonical Compound Discoverer table.
sample_cols <- grep(
  "^2025-F-[0-9]{4}_[0-9]{3}\\.raw \\(F[0-9]+\\)$",
  names(lcms),
  value = TRUE
)

if (length(sample_cols) != 120L) {
  stop(
    "Expected exactly 120 LC-MS sample columns; detected ", length(sample_cols),
    "\nFirst detected columns: ", paste(utils::head(sample_cols, 5), collapse = " | ")
  )
}

sample_ids <- sub("^2025-F-([0-9]{4})_.*$", "\\1", sample_cols)
if (anyDuplicated(sample_ids)) stop("Duplicated sample IDs extracted from LC-MS headers.")

missing_in_lcms <- setdiff(meta$Sample, sample_ids)
extra_in_lcms   <- setdiff(sample_ids, meta$Sample)
if (length(missing_in_lcms) || length(extra_in_lcms)) {
  stop(
    "LC-MS / metadata sample matching failed after 4-digit normalization.",
    if (length(missing_in_lcms)) paste0("\nMissing in LC-MS: ", paste(missing_in_lcms, collapse = ", ")) else "",
    if (length(extra_in_lcms)) paste0("\nExtra in LC-MS: ", paste(extra_in_lcms, collapse = ", ")) else ""
  )
}

## Reorder LC-MS exactly to metadata.
ord <- match(meta$Sample, sample_ids)
sample_cols <- sample_cols[ord]
sample_ids  <- sample_ids[ord]
if (!identical(sample_ids, as.character(meta$Sample))) {
  stop("Internal sample-order mismatch after matching.")
}

feature_ids_original <- if ("ID_MT" %in% names(lcms)) {
  as.character(lcms[["ID_MT"]])
} else {
  stop("ID_MT column is required in LCMS_FILE.")
}

if (anyNA(feature_ids_original) || any(feature_ids_original == "")) {
  stop("Missing ID_MT values in LC-MS input.")
}

## Same safeguard used in the previous LC-MS/WGCNA scripts: Compound Discoverer
## contains one rounded duplicate ID_MT in this table, so make.unique() is used
## for matrix row names while the first occurrence keeps the original ID.
dup_mask <- duplicated(feature_ids_original) | duplicated(feature_ids_original, fromLast = TRUE)
if (any(dup_mask)) {
  dup_tbl <- data.frame(
    Row = which(dup_mask),
    ID_MT_original = feature_ids_original[dup_mask],
    stringsAsFactors = FALSE
  )
  utils::write.csv(
    dup_tbl,
    file.path(DIR_QC, "duplicate_ID_MT_original_input.csv"),
    row.names = FALSE
  )
  warning(
    "Duplicated rounded ID_MT detected in canonical LC-MS input; applying make.unique() exactly as in the working pipelines."
  )
}
feature_ids <- make.unique(feature_ids_original)

X_raw <- as.matrix(lcms[, sample_cols, drop = FALSE])
suppressWarnings(storage.mode(X_raw) <- "numeric")
rownames(X_raw) <- feature_ids
colnames(X_raw) <- sample_ids

if (any(!is.finite(X_raw))) stop("Non-finite intensity values detected in LC-MS matrix.")
if (any(X_raw < 0, na.rm = TRUE)) stop("Negative peak areas detected; log2(x+1) is not valid.")

X_log <- log2(X_raw + 1)

feature_info <- data.frame(
  ID_MT = feature_ids,
  mz = suppressWarnings(as.numeric(lcms[["m/z"]])),
  RT_min = suppressWarnings(as.numeric(lcms[["RT [min]"]])),
  stringsAsFactors = FALSE
)

msg("LC-MS matching OK: %d features x %d samples", nrow(X_log), ncol(X_log))

utils::write.csv(
  data.frame(
    MetadataSample = meta$Sample,
    LCMSHeader = sample_cols,
    ExtractedLCMS_ID = sample_ids,
    stringsAsFactors = FALSE
  ),
  file.path(DIR_QC, "sample_matching_120_exact.csv"),
  row.names = FALSE
)

##############################
## 7. TECHNICIAN'S HISTORICAL HUB WORKBOOKS
##############################

xlsx_blue <- find_one_file(HUB_DIR, "^Hub_metabolites_blue.*\\.xlsx$")
xlsx_green <- find_one_file(HUB_DIR, "^Hub_metabolites_green.*\\.xlsx$")
xlsx_purple <- find_one_file(HUB_DIR, "^Hub_metabolites_purple.*\\.xlsx$")
xlsx_magenta <- find_one_file(HUB_DIR, "^Hub_metabolites_magenta.*\\.xlsx$")

specs <- data.frame(
  HistoricalModule = HIST_MODULE_LEVELS,
  FilePath = c(xlsx_blue, xlsx_green, xlsx_purple, xlsx_magenta),
  NameColumn = c("Name", "Name", "Name2", "Name"),
  stringsAsFactors = FALSE
)

read_historical_hub <- function(module_name, path, name_col) {
  dat <- readxl::read_excel(path)
  names(dat) <- trimws(names(dat))
  
  needed <- c("ID", "kME", name_col)
  miss <- setdiff(needed, names(dat))
  if (length(miss)) {
    stop("Workbook ", basename(path), " misses columns: ", paste(miss, collapse = ", "))
  }
  
  data.frame(
    HistoricalModule = module_name,
    ID_MT = as.character(dat[["ID"]]),
    Historical_kME = suppressWarnings(as.numeric(dat[["kME"]])),
    ManualAnnotation = clean_text(dat[[name_col]]),
    Workbook = basename(path),
    stringsAsFactors = FALSE
  ) |>
    dplyr::filter(!is.na(ID_MT), ID_MT != "") |>
    dplyr::distinct(ID_MT, .keep_all = TRUE)
}

hist_list <- lapply(seq_len(nrow(specs)), function(i) {
  read_historical_hub(
    specs$HistoricalModule[i],
    specs$FilePath[i],
    specs$NameColumn[i]
  )
})
hist_hubs <- dplyr::bind_rows(hist_list) |>
  dplyr::mutate(
    HistoricalModule = factor(HistoricalModule, levels = HIST_MODULE_LEVELS)
  ) |>
  dplyr::left_join(feature_info, by = "ID_MT")

## Historical technician files are expected to contain exactly:
## Blue 300, Green 140, Purple 80, Magenta 98 = 618 unique LC-MS features.
hist_counts <- hist_hubs |>
  dplyr::count(HistoricalModule, name = "N_hub_features")

expected_counts <- c("Blue" = 300L, "Green" = 140L, "Purple" = 80L, "Magenta" = 98L)
observed_counts <- setNames(hist_counts$N_hub_features, as.character(hist_counts$HistoricalModule))
for (nm in names(expected_counts)) {
  if (is.na(observed_counts[nm]) || observed_counts[nm] != expected_counts[nm]) {
    warning(
      "Historical workbook count differs for ", nm,
      ": expected ", expected_counts[nm], ", observed ", observed_counts[nm]
    )
  }
}

if (anyDuplicated(hist_hubs$ID_MT)) {
  stop("An ID_MT occurs in more than one historical hub workbook.")
}

missing_hub_ids <- setdiff(hist_hubs$ID_MT, rownames(X_log))
if (length(missing_hub_ids)) {
  stop(
    length(missing_hub_ids),
    " historical hub IDs are absent from the canonical LC-MS matrix. First IDs: ",
    paste(utils::head(missing_hub_ids, 10), collapse = ", ")
  )
}

msg("Historical hub workbooks OK: %d unique features", nrow(hist_hubs))
print(hist_counts)

##############################
## 8. OLD <-> NEW WGCNA COLOUR CONCORDANCE
##############################

new_files <- list.files(
  HUB_DIR,
  pattern = "^Module_.*_central_features\\.csv$",
  full.names = TRUE
)
if (!length(new_files)) stop("No Module_*_central_features.csv files found in HUB_DIR.")

new_map <- dplyr::bind_rows(lapply(new_files, function(f) {
  x <- utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("ID_MT", "ModuleColor", "kME") %in% names(x))) {
    stop("Required columns absent from: ", basename(f))
  }
  data.frame(
    ID_MT = as.character(x$ID_MT),
    NewModuleColor = as.character(x$ModuleColor),
    New_kME = suppressWarnings(as.numeric(x$kME)),
    stringsAsFactors = FALSE
  )
})) |>
  dplyr::distinct(ID_MT, .keep_all = TRUE)

hist_hubs <- hist_hubs |>
  dplyr::left_join(new_map, by = "ID_MT")

concordance <- hist_hubs |>
  dplyr::count(HistoricalModule, NewModuleColor, name = "N") |>
  dplyr::group_by(HistoricalModule) |>
  dplyr::mutate(Percent = 100 * N / sum(N)) |>
  dplyr::arrange(dplyr::desc(N), .by_group = TRUE) |>
  dplyr::ungroup()

utils::write.csv(
  concordance,
  file.path(DIR_TABLES, "Old_to_new_WGCNA_module_concordance.csv"),
  row.names = FALSE
)
utils::write.csv(
  hist_hubs,
  file.path(DIR_TABLES, "Historical_618_hub_features_with_new_colour.csv"),
  row.names = FALSE
)

## Strong verification: for these technician files, every historical feature
## should map to one single dominant new colour.
con_check <- concordance |>
  dplyr::group_by(HistoricalModule) |>
  dplyr::summarise(
    N_new_colours = dplyr::n_distinct(NewModuleColor, na.rm = TRUE),
    BestNewColour = NewModuleColor[which.max(N)],
    BestN = max(N),
    TotalN = sum(N),
    BestPercent = 100 * BestN / TotalN,
    .groups = "drop"
  )

utils::write.csv(
  con_check,
  file.path(DIR_QC, "module_concordance_check.csv"),
  row.names = FALSE
)

msg("Historical -> new WGCNA concordance:")
print(con_check)

##############################
## 9. FIGURE-5 FEATURE SET = STRUCTURE-NAMED ROWS ONLY
##############################

## Generic technician codes are useful analytically but are not structure names.
## They are therefore excluded from the MAIN labelled heatmap to avoid dozens
## of indistinguishable rows labelled only "fgt" or "flavonoide".
generic_codes <- c(
  "fgt", "fragment", "isotope", "flavonoide", "flavonoïde", "flavonoid"
)

fig5_features <- hist_hubs |>
  dplyr::mutate(
    AnnotationLower = tolower(trimws(ManualAnnotation)),
    IsSpecificNamedStructure = !is.na(ManualAnnotation) &
      !(AnnotationLower %in% generic_codes)
  ) |>
  dplyr::filter(IsSpecificNamedStructure) |>
  dplyr::arrange(HistoricalModule, dplyr::desc(Historical_kME), ManualAnnotation) |>
  dplyr::mutate(
    DisplayLabel = paste0(
      ManualAnnotation,
      "  [m/z ", format(round(mz, 4), nsmall = 4, trim = TRUE),
      "; RT ", format(round(RT_min, 2), nsmall = 2, trim = TRUE), " min]"
    ),
    RowKey = paste(as.character(HistoricalModule), ID_MT, sep = "||")
  )

if (!nrow(fig5_features)) stop("No structure-named annotations detected in technician workbooks.")

fig5_counts <- fig5_features |>
  dplyr::count(HistoricalModule, name = "N_structure_named_features")

utils::write.csv(
  fig5_features,
  file.path(DIR_TABLES, "Figure5_structure_named_features.csv"),
  row.names = FALSE
)
utils::write.csv(
  fig5_counts,
  file.path(DIR_TABLES, "Figure5_structure_named_feature_counts.csv"),
  row.names = FALSE
)

msg("Figure 5 structure-named feature counts:")
print(fig5_counts)

##############################
## 10. FEATURE-WISE Z-SCORES ACROSS ALL 120 SAMPLES
##############################

all_needed_ids <- hist_hubs$ID_MT
X_hub_log <- X_log[all_needed_ids, , drop = FALSE]
X_hub_z   <- zscore_rows(X_hub_log)

## All 618 hub-feature sample-level Z-scores are preserved for reproducibility.
utils::write.csv(
  data.frame(ID_MT = rownames(X_hub_z), X_hub_z, check.names = FALSE),
  file.path(DIR_TABLES, "Historical_618_hub_features_sample_Zscores.csv"),
  row.names = FALSE
)

##############################
## 11. CONDITION ORDER
##############################

meta_for_plot <- meta |>
  dplyr::mutate(
    Genotype_chr = as.character(Genotype),
    Time_chr = as.character(Sampling_D),
    Temp_chr = as.character(Temperature_D),
    Trt_chr = as.character(Treatment_D),
    Condition = dplyr::if_else(
      Time_chr == "D0",
      paste(Genotype_chr, Time_chr, Temp_chr, sep = "_"),
      paste(Genotype_chr, Time_chr, Temp_chr, Trt_chr, sep = "_")
    )
  )

make_condition_order <- function(geno) {
  c(
    paste(geno, "D0", c("NH", "HS-7", "HS-2"), sep = "_"),
    as.vector(outer(
      c("NH", "HS-7", "HS-2"),
      c("Control", "Treated"),
      FUN = function(temp, trt) paste(geno, "D4", temp, trt, sep = "_")
    )),
    as.vector(outer(
      c("NH", "HS-7", "HS-2"),
      c("Control", "Treated"),
      FUN = function(temp, trt) paste(geno, "D10", temp, trt, sep = "_")
    ))
  )
}

## outer() is column-major; explicitly enforce NH-C, NH-T, HS7-C, HS7-T, HS2-C, HS2-T.
make_condition_order <- function(geno) {
  c(
    paste(geno, "D0", "NH", sep = "_"),
    paste(geno, "D0", "HS-7", sep = "_"),
    paste(geno, "D0", "HS-2", sep = "_"),
    paste(geno, "D4", "NH", "Control", sep = "_"),
    paste(geno, "D4", "NH", "Treated", sep = "_"),
    paste(geno, "D4", "HS-7", "Control", sep = "_"),
    paste(geno, "D4", "HS-7", "Treated", sep = "_"),
    paste(geno, "D4", "HS-2", "Control", sep = "_"),
    paste(geno, "D4", "HS-2", "Treated", sep = "_"),
    paste(geno, "D10", "NH", "Control", sep = "_"),
    paste(geno, "D10", "NH", "Treated", sep = "_"),
    paste(geno, "D10", "HS-7", "Control", sep = "_"),
    paste(geno, "D10", "HS-7", "Treated", sep = "_"),
    paste(geno, "D10", "HS-2", "Control", sep = "_"),
    paste(geno, "D10", "HS-2", "Treated", sep = "_")
  )
}

condition_order <- c(
  make_condition_order("PRESTO"),
  make_condition_order("ROBILA")
)

actual_conditions <- unique(meta_for_plot$Condition)
missing_conditions <- setdiff(condition_order, actual_conditions)
extra_conditions <- setdiff(actual_conditions, condition_order)
if (length(missing_conditions) || length(extra_conditions)) {
  stop(
    "Unexpected condition structure.",
    if (length(missing_conditions)) paste0("\nMissing: ", paste(missing_conditions, collapse = ", ")) else "",
    if (length(extra_conditions)) paste0("\nExtra: ", paste(extra_conditions, collapse = ", ")) else ""
  )
}

## Compact x labels. Genotype is indicated by a top line here and the strong separator.
condition_labels <- c(
  ## PRESTO
  "PRESTO\nD0 NH", "D0 HS-7", "D0 HS-2",
  "D4 NH C", "D4 NH T", "D4 HS-7 C", "D4 HS-7 T", "D4 HS-2 C", "D4 HS-2 T",
  "D10 NH C", "D10 NH T", "D10 HS-7 C", "D10 HS-7 T", "D10 HS-2 C", "D10 HS-2 T",
  ## ROBILA
  "ROBILA\nD0 NH", "D0 HS-7", "D0 HS-2",
  "D4 NH C", "D4 NH T", "D4 HS-7 C", "D4 HS-7 T", "D4 HS-2 C", "D4 HS-2 T",
  "D10 NH C", "D10 NH T", "D10 HS-7 C", "D10 HS-7 T", "D10 HS-2 C", "D10 HS-2 T"
)
names(condition_labels) <- condition_order

##############################
## 12. FIGURE 5: CONDITION-MEAN HEATMAPS
##     (Z-score was calculated before replicate averaging)
##############################

fig5_ids <- fig5_features$ID_MT
X_fig5_z <- X_hub_z[fig5_ids, , drop = FALSE]

condition_samples <- split(meta_for_plot$Sample, meta_for_plot$Condition)
condition_samples <- condition_samples[condition_order]
if (any(vapply(condition_samples, length, integer(1)) != 4L)) {
  stop("Every plotted condition should contain exactly 4 biological replicates.")
}

fig5_condition_mean <- sapply(condition_samples, function(samps) {
  rowMeans(X_fig5_z[, samps, drop = FALSE], na.rm = TRUE)
})
rownames(fig5_condition_mean) <- fig5_ids
colnames(fig5_condition_mean) <- condition_order

utils::write.csv(
  data.frame(ID_MT = rownames(fig5_condition_mean), fig5_condition_mean, check.names = FALSE),
  file.path(DIR_TABLES, "Figure5_condition_mean_Zscores.csv"),
  row.names = FALSE
)

fig5_long <- data.frame(
  ID_MT = rownames(fig5_condition_mean),
  fig5_condition_mean,
  check.names = FALSE
) |>
  tidyr::pivot_longer(
    cols = -ID_MT,
    names_to = "Condition",
    values_to = "MeanZ"
  ) |>
  dplyr::left_join(
    fig5_features |>
      dplyr::select(ID_MT, HistoricalModule, Historical_kME, ManualAnnotation, DisplayLabel, RowKey),
    by = "ID_MT"
  ) |>
  dplyr::mutate(
    HistoricalModule = factor(HistoricalModule, levels = HIST_MODULE_LEVELS),
    Condition = factor(Condition, levels = condition_order)
  )

## Explicit row order: historical module order, then decreasing historical kME.
row_order <- fig5_features |>
  dplyr::arrange(HistoricalModule, dplyr::desc(Historical_kME), ManualAnnotation)
row_levels <- rev(row_order$RowKey)
row_label_map <- stats::setNames(row_order$DisplayLabel, row_order$RowKey)
fig5_long$RowKey <- factor(fig5_long$RowKey, levels = row_levels)

p5 <- ggplot2::ggplot(
  fig5_long,
  ggplot2::aes(x = Condition, y = RowKey, fill = MeanZ)
) +
  ggplot2::geom_tile() +
  ## D0/D4/D10 separators within PRESTO and ROBILA.
  ggplot2::geom_vline(xintercept = c(3.5, 9.5, 18.5, 24.5), linewidth = 0.35, colour = "grey30") +
  ## Strong genotype separator.
  ggplot2::geom_vline(xintercept = 15.5, linewidth = 1.0, colour = "black") +
  ggplot2::facet_grid(
    HistoricalModule ~ .,
    scales = "free_y",
    space = "free_y"
  ) +
  ggplot2::scale_x_discrete(labels = function(x) unname(condition_labels[x])) +
  ggplot2::scale_y_discrete(labels = function(x) unname(row_label_map[x])) +
  ggplot2::scale_fill_gradient2(
    low = HEAT_LOW,
    mid = HEAT_MID,
    high = HEAT_HIGH,
    midpoint = 0,
    limits = c(-2.5, 2.5),
    oob = scales::squish,
    name = "Relative\nabundance\n(Z-score)"
  ) +
  ggplot2::labs(
    title = "Manually structure-annotated LC-MS hub features",
    subtitle = "Historical ROBILA-D4 WGCNA modules projected across PRESTO and ROBILA; C = Control, T = Treated",
    x = NULL,
    y = NULL
  ) +
  theme_pub(10) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 55, hjust = 1, vjust = 1, size = 7),
    axis.text.y = ggplot2::element_text(size = 6),
    strip.text.y = ggplot2::element_text(face = "bold", size = 10)
  )

save_plot(
  p5,
  "FIG05_Manually_annotated_hub_heatmaps_COMBINED",
  width = 16,
  height = 14,
  out_dir = DIR_FIG5
)

## One publication-quality heatmap per historical module.
for (mod in HIST_MODULE_LEVELS) {
  dd <- fig5_long |>
    dplyr::filter(as.character(HistoricalModule) == mod)
  
  nrows_mod <- dplyr::n_distinct(dd$ID_MT)
  if (!nrows_mod) next
  
  ## Rebuild module-specific labels/levels to eliminate unused rows.
  mod_info <- fig5_features |>
    dplyr::filter(as.character(HistoricalModule) == mod) |>
    dplyr::arrange(dplyr::desc(Historical_kME), ManualAnnotation)
  mod_levels <- rev(mod_info$RowKey)
  mod_labels <- stats::setNames(mod_info$DisplayLabel, mod_info$RowKey)
  dd$RowKey <- factor(as.character(dd$RowKey), levels = mod_levels)
  
  pp <- ggplot2::ggplot(dd, ggplot2::aes(x = Condition, y = RowKey, fill = MeanZ)) +
    ggplot2::geom_tile() +
    ggplot2::geom_vline(xintercept = c(3.5, 9.5, 18.5, 24.5), linewidth = 0.35, colour = "grey30") +
    ggplot2::geom_vline(xintercept = 15.5, linewidth = 1.0, colour = "black") +
    ggplot2::scale_x_discrete(labels = function(x) unname(condition_labels[x])) +
    ggplot2::scale_y_discrete(labels = function(x) unname(mod_labels[x])) +
    ggplot2::scale_fill_gradient2(
      low = HEAT_LOW, mid = HEAT_MID, high = HEAT_HIGH,
      midpoint = 0, limits = c(-2.5, 2.5), oob = scales::squish,
      name = "Relative\nabundance\n(Z-score)"
    ) +
    ggplot2::labs(
      title = paste0("Historical ", mod, " module"),
      subtitle = "Manually structure-named features; historical ROBILA-D4 hub ions projected across the full experiment",
      x = NULL, y = NULL
    ) +
    theme_pub(10) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 55, hjust = 1, vjust = 1, size = 7),
      axis.text.y = ggplot2::element_text(size = 7)
    )
  
  h <- max(5.5, min(13, 3.5 + 0.26 * nrows_mod))
  save_plot(
    pp,
    paste0("FIG05_Heatmap_", mod),
    width = 16,
    height = h,
    out_dir = DIR_FIG5
  )
}

##############################
## 13. FIGURE 6: SAMPLE-LEVEL PROJECTED HISTORICAL HUB SCORES
##     Uses ALL features delivered in each historical hub workbook.
##############################

module_sample_scores <- dplyr::bind_rows(lapply(HIST_MODULE_LEVELS, function(mod) {
  ids <- hist_hubs |>
    dplyr::filter(as.character(HistoricalModule) == mod) |>
    dplyr::pull(ID_MT)
  
  if (!length(ids)) return(NULL)
  
  data.frame(
    HistoricalModule = mod,
    Sample = colnames(X_hub_z),
    ProjectedHubScore = colMeans(X_hub_z[ids, , drop = FALSE], na.rm = TRUE),
    N_hub_features = length(ids),
    stringsAsFactors = FALSE
  )
})) |>
  dplyr::left_join(
    meta |>
      dplyr::select(Sample, Genotype, Sampling_D, Temperature_D, Treatment_D, Replicat),
    by = "Sample"
  ) |>
  dplyr::mutate(
    HistoricalModule = factor(HistoricalModule, levels = HIST_MODULE_LEVELS),
    Genotype = factor(Genotype, levels = GENO_LEVELS),
    Sampling_D = factor(Sampling_D, levels = TIME_LEVELS),
    Temperature_D = factor(Temperature_D, levels = TEMP_LEVELS),
    Treatment_D = factor(Treatment_D, levels = TRT_LEVELS)
  )

if (anyNA(module_sample_scores$Genotype)) stop("Failed to attach metadata to projected hub scores.")

utils::write.csv(
  module_sample_scores,
  file.path(DIR_TABLES, "Figure6_projected_hub_scores_SAMPLE_LEVEL.csv"),
  row.names = FALSE
)

score_summary <- module_sample_scores |>
  dplyr::group_by(
    HistoricalModule, Genotype, Sampling_D, Temperature_D, Treatment_D
  ) |>
  dplyr::summarise(
    MeanScore = mean(ProjectedHubScore, na.rm = TRUE),
    SE = se_fun(ProjectedHubScore),
    N = dplyr::n(),
    .groups = "drop"
  )

utils::write.csv(
  score_summary,
  file.path(DIR_TABLES, "Figure6_projected_hub_scores_CONDITION_SUMMARY.csv"),
  row.names = FALSE
)

##############################
## 14. FIGURE 6A: CONTROL TRAJECTORY D0 -> D4 -> D10
##############################

control_traj <- score_summary |>
  dplyr::filter(Treatment_D == "Control")

p6a <- ggplot2::ggplot(
  control_traj,
  ggplot2::aes(
    x = Sampling_D,
    y = MeanScore,
    group = Temperature_D,
    colour = Temperature_D
  )
) +
  ggplot2::geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.35) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::geom_point(size = 2.2) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = MeanScore - SE, ymax = MeanScore + SE),
    width = 0.08,
    linewidth = 0.45
  ) +
  ggplot2::facet_grid(HistoricalModule ~ Genotype) +
  ggplot2::scale_colour_manual(values = TEMP_COLORS, drop = FALSE) +
  ggplot2::labs(
    title = "Control trajectories of historical LC-MS hub programmes",
    subtitle = "Mean projected hub-feature score ± SE across four biological replicates",
    x = NULL,
    y = "Projected historical hub score",
    colour = "Temperature"
  ) +
  theme_pub(11)

save_plot(
  p6a,
  "FIG06A_Control_trajectories_D0_D4_D10",
  width = 11.5,
  height = 9.5,
  out_dir = DIR_FIG6
)

##############################
## 15. FIGURE 6B: TREATMENT EFFECT AT D4 AND D10
##     Delta = Treated - Control projected hub score
##############################

trt_wide <- score_summary |>
  dplyr::filter(Sampling_D %in% c("D4", "D10")) |>
  dplyr::select(
    HistoricalModule, Genotype, Sampling_D, Temperature_D,
    Treatment_D, MeanScore, SE
  ) |>
  tidyr::pivot_wider(
    names_from = Treatment_D,
    values_from = c(MeanScore, SE),
    names_sep = "_"
  )

needed_delta_cols <- c(
  "MeanScore_Control", "MeanScore_Treated",
  "SE_Control", "SE_Treated"
)
if (!all(needed_delta_cols %in% names(trt_wide))) {
  stop("Could not construct Treated-Control score table for D4/D10.")
}

trt_delta <- trt_wide |>
  dplyr::mutate(
    TreatmentDelta = MeanScore_Treated - MeanScore_Control,
    SE_Delta = sqrt(SE_Treated^2 + SE_Control^2)
  )

utils::write.csv(
  trt_delta,
  file.path(DIR_TABLES, "Figure6B_Treated_minus_Control_projected_hub_scores.csv"),
  row.names = FALSE
)

p6b <- ggplot2::ggplot(
  trt_delta,
  ggplot2::aes(
    x = Sampling_D,
    y = TreatmentDelta,
    group = Temperature_D,
    colour = Temperature_D
  )
) +
  ggplot2::geom_hline(yintercept = 0, colour = "black", linewidth = 0.45) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::geom_point(size = 2.2) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = TreatmentDelta - SE_Delta, ymax = TreatmentDelta + SE_Delta),
    width = 0.06,
    linewidth = 0.45
  ) +
  ggplot2::facet_grid(HistoricalModule ~ Genotype) +
  ggplot2::scale_colour_manual(values = TEMP_COLORS, drop = FALSE) +
  ggplot2::labs(
    title = "Treatment-dependent trajectories of historical LC-MS hub programmes",
    subtitle = "Difference in projected hub score (Treated - Control); positive values indicate higher abundance in Treated plants",
    x = NULL,
    y = "Treated - Control hub score",
    colour = "Temperature"
  ) +
  theme_pub(11)

save_plot(
  p6b,
  "FIG06B_Treatment_effect_trajectories_D4_D10",
  width = 11.5,
  height = 9.5,
  out_dir = DIR_FIG6
)

##############################
## 16. OPTIONAL SUPPORTING FIGURE: RAW CONTROL + TREATED SCORES D4/D10
##############################

raw_d4d10 <- score_summary |>
  dplyr::filter(Sampling_D %in% c("D4", "D10")) |>
  dplyr::mutate(
    TrajectoryGroup = interaction(Temperature_D, Treatment_D, drop = TRUE)
  )

p6supp <- ggplot2::ggplot(
  raw_d4d10,
  ggplot2::aes(
    x = Sampling_D,
    y = MeanScore,
    group = TrajectoryGroup,
    colour = Temperature_D,
    linetype = Treatment_D
  )
) +
  ggplot2::geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.35) +
  ggplot2::geom_line(linewidth = 0.75) +
  ggplot2::geom_point(size = 2.0) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = MeanScore - SE, ymax = MeanScore + SE),
    width = 0.05,
    linewidth = 0.4
  ) +
  ggplot2::facet_grid(HistoricalModule ~ Genotype) +
  ggplot2::scale_colour_manual(values = TEMP_COLORS, drop = FALSE) +
  ggplot2::scale_linetype_manual(
    values = c("Control" = "solid", "Treated" = "dashed"),
    drop = FALSE
  ) +
  ggplot2::labs(
    title = "Supporting view: projected hub scores at D4 and D10",
    x = NULL,
    y = "Projected historical hub score",
    colour = "Temperature",
    linetype = "Treatment"
  ) +
  theme_pub(10)

save_plot(
  p6supp,
  "SUPP_FIG_Projected_hub_scores_Control_vs_Treated_D4_D10",
  width = 11.5,
  height = 9.5,
  out_dir = DIR_FIG6
)

##############################
## 17. RUN SUMMARY / SELF-CHECK
##############################

summary_lines <- c(
  "FIGURE 5 + FIGURE 6 LC-MS annotated-hub analysis",
  "",
  paste0("META_FILE: ", META_FILE),
  paste0("LCMS_FILE: ", LCMS_FILE),
  paste0("WGCNA_DIR: ", WGCNA_DIR),
  paste0("OUT_DIR: ", OUT_DIR),
  "",
  paste0("Metadata rows: ", nrow(meta)),
  paste0("LC-MS features: ", nrow(X_log)),
  paste0("LC-MS samples: ", ncol(X_log)),
  paste0("Historical hub features: ", nrow(hist_hubs)),
  "",
  "Historical workbook counts:",
  paste0("  ", as.character(hist_counts$HistoricalModule), ": ", hist_counts$N_hub_features),
  "",
  "Figure-5 structure-named feature counts:",
  paste0("  ", as.character(fig5_counts$HistoricalModule), ": ", fig5_counts$N_structure_named_features),
  "",
  "Old -> new WGCNA concordance:",
  paste0(
    "  ", as.character(con_check$HistoricalModule), " -> ", con_check$BestNewColour,
    " : ", con_check$BestN, "/", con_check$TotalN,
    " (", round(con_check$BestPercent, 1), "%)"
  ),
  "",
  "Figure 5:",
  "  Z-score computed feature-wise across 120 samples BEFORE replicate averaging.",
  "  Main heatmap uses manually structure-named rows only.",
  "",
  "Figure 6:",
  "  Projected hub score uses ALL historical hub features in each technician workbook.",
  "  6A = Control trajectory D0-D4-D10.",
  "  6B = Treated-Control trajectory D4-D10.",
  "",
  "Analysis completed successfully."
)
writeLines(summary_lines, file.path(DIR_QC, "RUN_SUMMARY.txt"))

cat("\n============================================================\n")
cat("SUCCESS\n")
cat("============================================================\n")
cat(paste(summary_lines, collapse = "\n"), "\n")
cat("============================================================\n")

############################################################
## 16B. ADDITIONAL TRAJECTORY FIGURES
##
## FIGURE 6C:
##   Shared D0 baseline -> Control + Treated trajectories
##   separated by temperature.
##
## FIGURE 6D:
##   Shared D0 baseline -> Treated trajectory only
##   separated by temperature.
##
## GRAPHICAL CONVENTION
##   Treated = SOLID line + filled point
##   Control = DASHED line + open point
##
## IMPORTANT:
##   D0 is BEFORE treatment.
##   Therefore there is only ONE observed D0 condition.
##   D0 is duplicated below ONLY for graphical construction
##   of the Control and Treated trajectories.
##
##   It must NOT be interpreted as an experimentally observed
##   "Treated D0" group and must NOT be duplicated for statistics.
############################################################


############################################################
## 16B-1. FIGURE 6C
## SHARED D0 BASELINE ->
## CONTROL AND TREATED TRAJECTORIES
## D0 -> D4 -> D10
############################################################

##############################
## A. REAL D0 BASELINE
##############################

baseline_D0 <- score_summary |>
  dplyr::filter(
    Sampling_D == "D0",
    Treatment_D == "Control"
  )

## Check expected design:
## 4 modules x 2 genotypes x 3 temperatures = 24 rows
if (nrow(baseline_D0) !=
    length(HIST_MODULE_LEVELS) *
    length(GENO_LEVELS) *
    length(TEMP_LEVELS)) {
  
  warning(
    "Unexpected number of D0 baseline rows: ",
    nrow(baseline_D0)
  )
}


##############################
## B. DUPLICATE D0 ONLY FOR
##    GRAPHICAL TRAJECTORIES
##############################

## Control trajectory starts from the actual D0 baseline.
baseline_control <- baseline_D0 |>
  dplyr::mutate(
    PlotTreatment = "Control",
    BaselineType = "Shared pre-treatment D0"
  )

## Treated trajectory starts from the SAME actual D0 baseline.
## This is a graphical duplicate only.
baseline_treated <- baseline_D0 |>
  dplyr::mutate(
    PlotTreatment = "Treated",
    BaselineType = "Shared pre-treatment D0"
  )


##############################
## C. REAL D4 + D10 DATA
##############################

post_treatment <- score_summary |>
  dplyr::filter(
    Sampling_D %in% c("D4", "D10"),
    Treatment_D %in% c("Control", "Treated")
  ) |>
  dplyr::mutate(
    PlotTreatment = as.character(Treatment_D),
    BaselineType = "Observed post-treatment"
  )


##############################
## D. COMBINE THE TRAJECTORIES
##############################

trajectory_shared_D0 <- dplyr::bind_rows(
  baseline_control,
  baseline_treated,
  post_treatment
) |>
  dplyr::mutate(
    
    Sampling_D = factor(
      as.character(Sampling_D),
      levels = c("D0", "D4", "D10")
    ),
    
    Temperature_D = factor(
      as.character(Temperature_D),
      levels = c("NH", "HS-7", "HS-2")
    ),
    
    PlotTreatment = factor(
      PlotTreatment,
      levels = c("Control", "Treated")
    ),
    
    HistoricalModule = factor(
      as.character(HistoricalModule),
      levels = HIST_MODULE_LEVELS
    ),
    
    Genotype = factor(
      as.character(Genotype),
      levels = GENO_LEVELS
    ),
    
    ## One independent graphical trajectory for each:
    ## Temperature x Treatment
    TrajectoryGroup = interaction(
      Temperature_D,
      PlotTreatment,
      drop = TRUE
    )
  )


##############################
## E. EXPORT FIGURE 6C DATA
##############################

utils::write.csv(
  trajectory_shared_D0,
  file.path(
    DIR_TABLES,
    "Figure6C_SHARED_D0_Control_Treated_Temperature_trajectories.csv"
  ),
  row.names = FALSE
)


##############################
## F. FIGURE 6C
##
## Treated:
##   solid line
##   filled circle
##
## Control:
##   dashed line
##   open circle
##
## Temperature:
##   colour
##############################

p6c <- ggplot2::ggplot(
  trajectory_shared_D0,
  ggplot2::aes(
    x = Sampling_D,
    y = MeanScore,
    group = TrajectoryGroup,
    colour = Temperature_D,
    linetype = PlotTreatment,
    shape = PlotTreatment
  )
) +
  
  ## Zero reference
  ggplot2::geom_hline(
    yintercept = 0,
    colour = "grey70",
    linewidth = 0.35
  ) +
  
  ## Trajectories
  ggplot2::geom_line(
    linewidth = 0.95
  ) +
  
  ## Mean points
  ggplot2::geom_point(
    size = 2.8,
    stroke = 0.9
  ) +
  
  ## Biological replicate SE
  ggplot2::geom_errorbar(
    ggplot2::aes(
      ymin = MeanScore - SE,
      ymax = MeanScore + SE
    ),
    width = 0.07,
    linewidth = 0.45
  ) +
  
  ## Rows = historical metabolic programmes
  ## Columns = genotypes
  ggplot2::facet_grid(
    HistoricalModule ~ Genotype
  ) +
  
  ## Temperature colours
  ggplot2::scale_colour_manual(
    values = TEMP_COLORS,
    breaks = c("NH", "HS-7", "HS-2"),
    drop = FALSE
  ) +
  
  ## IMPORTANT GRAPHICAL CONVENTION
  ##
  ## Treated = SOLID
  ## Control = DASHED
  ggplot2::scale_linetype_manual(
    values = c(
      "Control" = "dashed",
      "Treated" = "solid"
    ),
    breaks = c("Treated", "Control"),
    drop = FALSE
  ) +
  
  ## Treated = filled point
  ## Control = open point
  ggplot2::scale_shape_manual(
    values = c(
      "Control" = 1,
      "Treated" = 16
    ),
    breaks = c("Treated", "Control"),
    drop = FALSE
  ) +
  
  ggplot2::labs(
    title =
      "Temperature- and treatment-dependent trajectories of historical LC-MS hub programmes",
    
    subtitle =
      "D0 is the shared pre-treatment baseline; Control and Treated trajectories diverge from D4 onwards",
    
    x = NULL,
    
    y =
      "Projected historical hub score",
    
    colour =
      "Temperature",
    
    linetype =
      "Treatment",
    
    shape =
      "Treatment"
  ) +
  
  theme_pub(11) +
  
  ggplot2::theme(
    
    ## Cleaner facets
    strip.text =
      ggplot2::element_text(
        face = "bold",
        size = 10
      ),
    
    ## Keep legend compact
    legend.position =
      "top",
    
    legend.box =
      "horizontal"
  )


save_plot(
  p6c,
  
  "FIG06C_SHARED_D0_Control_Treated_Temperature_trajectories",
  
  width = 12,
  height = 9.5,
  
  out_dir = DIR_FIG6
)



############################################################
## 16B-2. FIGURE 6D
## SHARED D0 BASELINE ->
## TREATED TRAJECTORY ONLY
##
## D0 -> D4 -> D10
## NH / HS-7 / HS-2
############################################################


##############################
## A. BUILD TREATED TRAJECTORY
##############################

treated_with_D0 <- dplyr::bind_rows(
  
  ################################
  ## Shared PRE-TREATMENT baseline
  ################################
  
  score_summary |>
    dplyr::filter(
      Sampling_D == "D0",
      Treatment_D == "Control"
    ) |>
    dplyr::mutate(
      PlotTreatment = "Treated",
      BaselineType = "Shared pre-treatment D0"
    ),
  
  ################################
  ## Actual Treated observations
  ################################
  
  score_summary |>
    dplyr::filter(
      Sampling_D %in% c("D4", "D10"),
      Treatment_D == "Treated"
    ) |>
    dplyr::mutate(
      PlotTreatment = "Treated",
      BaselineType = "Observed Treated"
    )
) |>
  
  dplyr::mutate(
    
    Sampling_D = factor(
      as.character(Sampling_D),
      levels = c("D0", "D4", "D10")
    ),
    
    Temperature_D = factor(
      as.character(Temperature_D),
      levels = c("NH", "HS-7", "HS-2")
    ),
    
    HistoricalModule = factor(
      as.character(HistoricalModule),
      levels = HIST_MODULE_LEVELS
    ),
    
    Genotype = factor(
      as.character(Genotype),
      levels = GENO_LEVELS
    )
  )


##############################
## B. EXPORT FIGURE 6D DATA
##############################

utils::write.csv(
  treated_with_D0,
  file.path(
    DIR_TABLES,
    "Figure6D_SHARED_D0_to_Treated_temperature_trajectories.csv"
  ),
  row.names = FALSE
)


##############################
## C. FIGURE 6D
##
## Treated only:
##   SOLID lines
##   filled points
##############################

p6d <- ggplot2::ggplot(
  treated_with_D0,
  ggplot2::aes(
    x = Sampling_D,
    y = MeanScore,
    group = Temperature_D,
    colour = Temperature_D
  )
) +
  
  ggplot2::geom_hline(
    yintercept = 0,
    colour = "grey70",
    linewidth = 0.35
  ) +
  
  ## Treated trajectories = SOLID
  ggplot2::geom_line(
    linewidth = 0.95,
    linetype = "solid"
  ) +
  
  ## Treated points = FILLED
  ggplot2::geom_point(
    size = 2.8,
    shape = 16
  ) +
  
  ggplot2::geom_errorbar(
    ggplot2::aes(
      ymin = MeanScore - SE,
      ymax = MeanScore + SE
    ),
    width = 0.07,
    linewidth = 0.45
  ) +
  
  ggplot2::facet_grid(
    HistoricalModule ~ Genotype
  ) +
  
  ggplot2::scale_colour_manual(
    values = TEMP_COLORS,
    breaks = c("NH", "HS-7", "HS-2"),
    drop = FALSE
  ) +
  
  ggplot2::labs(
    
    title =
      "Trajectories of treated plants from the shared pre-treatment baseline",
    
    subtitle =
      "D0 represents the shared pre-treatment state; D4 and D10 correspond to observed Treated plants",
    
    x = NULL,
    
    y =
      "Projected historical hub score",
    
    colour =
      "Temperature"
  ) +
  
  theme_pub(11) +
  
  ggplot2::theme(
    strip.text =
      ggplot2::element_text(
        face = "bold",
        size = 10
      ),
    legend.position =
      "top"
  )


save_plot(
  p6d,
  
  "FIG06D_SHARED_D0_to_Treated_Temperature_trajectories",
  
  width = 11.5,
  height = 9.5,
  
  out_dir = DIR_FIG6
)



############################################################
## 16B-3. OPTIONAL FIGURE 6E
## CONTROL TRAJECTORY ONLY
##
## Same graphical convention:
## Control = DASHED + open point
##
## This duplicates the biological information of Figure 6A,
## but uses the final graphical convention chosen for the article.
############################################################

control_with_D0 <- score_summary |>
  dplyr::filter(
    Treatment_D == "Control"
  ) |>
  dplyr::mutate(
    
    Sampling_D = factor(
      as.character(Sampling_D),
      levels = c("D0", "D4", "D10")
    ),
    
    Temperature_D = factor(
      as.character(Temperature_D),
      levels = c("NH", "HS-7", "HS-2")
    ),
    
    HistoricalModule = factor(
      as.character(HistoricalModule),
      levels = HIST_MODULE_LEVELS
    ),
    
    Genotype = factor(
      as.character(Genotype),
      levels = GENO_LEVELS
    )
  )


utils::write.csv(
  control_with_D0,
  file.path(
    DIR_TABLES,
    "Figure6E_Control_temperature_trajectories.csv"
  ),
  row.names = FALSE
)


p6e <- ggplot2::ggplot(
  control_with_D0,
  ggplot2::aes(
    x = Sampling_D,
    y = MeanScore,
    group = Temperature_D,
    colour = Temperature_D
  )
) +
  
  ggplot2::geom_hline(
    yintercept = 0,
    colour = "grey70",
    linewidth = 0.35
  ) +
  
  ## Control = DASHED
  ggplot2::geom_line(
    linewidth = 0.95,
    linetype = "dashed"
  ) +
  
  ## Control = OPEN POINT
  ggplot2::geom_point(
    size = 2.8,
    shape = 1,
    stroke = 0.9
  ) +
  
  ggplot2::geom_errorbar(
    ggplot2::aes(
      ymin = MeanScore - SE,
      ymax = MeanScore + SE
    ),
    width = 0.07,
    linewidth = 0.45
  ) +
  
  ggplot2::facet_grid(
    HistoricalModule ~ Genotype
  ) +
  
  ggplot2::scale_colour_manual(
    values = TEMP_COLORS,
    breaks = c("NH", "HS-7", "HS-2"),
    drop = FALSE
  ) +
  
  ggplot2::labs(
    
    title =
      "Control trajectories of historical LC-MS hub programmes",
    
    subtitle =
      "Observed Control plants from D0 to D10 across temperature regimes",
    
    x = NULL,
    
    y =
      "Projected historical hub score",
    
    colour =
      "Temperature"
  ) +
  
  theme_pub(11) +
  
  ggplot2::theme(
    strip.text =
      ggplot2::element_text(
        face = "bold",
        size = 10
      ),
    legend.position =
      "top"
  )


save_plot(
  p6e,
  
  "FIG06E_Control_Temperature_trajectories",
  
  width = 11.5,
  height = 9.5,
  
  out_dir = DIR_FIG6
)