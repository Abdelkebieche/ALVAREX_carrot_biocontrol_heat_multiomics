############################################################
## maSigPro_D0_D2_D4_SharedBaseline.R
##
## Design:
##   D0 (P1): Control only — shared baseline for both trajectories
##   D2 (P2): Control + Treated
##   D4 (P3): Control + Treated
##
## For each Temperature (NH, HS-7, HS-2):
##   Control trajectory: D0 → D2 → D4
##   Treated trajectory: D0 (=Control) → D2 → D4
##
## Pipeline:
##   1. Import + filter (genotype + P1/P2/P3)
##   2. Count filtering + normalization
##   3. MAD filtering
##   4. maSigPro design with shared baseline
##   5. maSigPro p.vector + T.fit
##   6. Extract significant genes
##   7. Mean profile matrix (display labels)
##   8. Hierarchical clustering
##   9. Heatmap (pheatmap)
##  10. Mean cluster profiles
##  10b. Individual gene profiles
##  11. Cluster size barplot
##  12. Top variable genes per cluster
##  13. GO ORA + Pathway ORA per cluster (PlanT2T)
##  14. Export gene lists compatible with intersection script
##  15. Summary
##
## Display labels:
##   T0 → NH, T1 → HS-7, T2 → HS-2
##   EAU → Control, SDP → Treated
##   P1 → D0, P2 → D2, P3 → D4
##
## Output format: PDF + PNG + SVG
# Repository configuration (portable; no HPC-specific paths)
if (!exists('CARROT_CONFIG_LOADED')) source('../../../config/config.R')


############################################################

##############################
## PARAMETERS
##############################
counts_file <- RNA_COUNTS_FILE
meta_file <- RNA_METADATA_FILE
base_out <- file.path(MASIGPRO_ROOT, paste0(GENOTYPE_CHOSEN, "_maSigPro"))

GENOTYPE_CHOSEN <- Sys.getenv("CARROT_GENOTYPE", unset = "PRESTO")

## Count filtering
MIN_TOTAL_COUNTS <- 30
MIN_COUNT        <- 5
MIN_SAMPLES      <- 3

## MAD filtering (keeps top 50% most variable genes)
MAD_QUANTILE <- 0.5

## maSigPro parameters
MASIGPRO_Q       <- 0.01
MASIGPRO_STEP    <- "backward"
MASIGPRO_DEGREE  <- 2       # degree 2 for 3 time points (quadratic possible)
MASIGPRO_NVAR    <- 2

## Clustering
K_CLUSTERS    <- 8
HCLUST_METHOD <- "ward.D2"
SEED          <- 123

## Enrichment
RUN_ENRICHMENT       <- TRUE
INSTALL_MISSING_PKGS <- FALSE

FDR_CUTOFF       <- 0.05
GO_ONT_LIST      <- c("BP", "MF", "CC")
GO_SHOWCATEGORY  <- 20
MIN_GENES_ORA    <- 10
MIN_GS_SIZE_PATH <- 5

ORGDB_URL <- "https://biobigdata.nju.edu.cn/plant2t/orgdb/org.Daucus.carota.DH13M14.eg.db_1.0.tar.gz"
ORGDB_PKG <- "org.Daucus.carota.DH13M14.eg.db"
map_file <- GENE_ID_MAP_FILE

##############################
## DISPLAY LABELS & PALETTES
##############################
temp_display <- c("T0" = "NH", "T1" = "HS-7", "T2" = "HS-2")
trt_display  <- c("EAU" = "Control", "SDP" = "Treated")
time_display <- c("P1" = "D0", "P2" = "D2", "P3" = "D4")

TEMP_COLORS <- c("T0" = "#313695", "T1" = "#74ADD1", "T2" = "#F46D43")
TRT_COLORS  <- c("EAU" = "#4575B4", "SDP" = "#D73027")
TIME_COLORS <- c("P1" = "#AAAAAA", "P2" = "#56B4E9", "P3" = "#009E73")

# Display-keyed palettes (for figure legends)
TEMP_COLORS_D <- setNames(TEMP_COLORS, temp_display[names(TEMP_COLORS)])
TRT_COLORS_D  <- setNames(TRT_COLORS, trt_display[names(TRT_COLORS)])
TIME_COLORS_D <- setNames(TIME_COLORS, time_display[names(TIME_COLORS)])

relabel_temp <- function(x) {
  r <- temp_display[as.character(x)]
  ifelse(is.na(r), as.character(x), r)
}
relabel_trt <- function(x) {
  r <- trt_display[as.character(x)]
  ifelse(is.na(r), as.character(x), r)
}
relabel_time <- function(x) {
  r <- time_display[as.character(x)]
  ifelse(is.na(r), as.character(x), r)
}

##############################
## INITIALIZATION
##############################
set.seed(SEED)
options(stringsAsFactors = FALSE)

dir.create(base_out, recursive = TRUE, showWarnings = FALSE)

dirs <- list(
  inputs     = file.path(base_out, "00_inputs"),
  masigpro   = file.path(base_out, "01_maSigPro"),
  clusters   = file.path(base_out, "02_clusters"),
  figures    = file.path(base_out, "03_figures"),
  tables     = file.path(base_out, "04_tables"),
  enrichment = file.path(base_out, "05_enrichment")
)
for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)

logfile <- file.path(base_out, "pipeline.log")
log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", paste0(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = logfile, append = TRUE)
}

##############################
## PACKAGE LOADING
##############################
load_pkgs <- function(pkgs, install_missing = FALSE, bioc = FALSE) {
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      if (!install_missing) stop("Missing package: ", p)
      if (bioc) {
        if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
        BiocManager::install(p, ask = FALSE, update = FALSE)
      } else {
        install.packages(p)
      }
    }
    suppressPackageStartupMessages(library(p, character.only = TRUE))
  }
}

load_pkgs(
  c("dplyr", "tidyr", "tibble", "ggplot2", "pheatmap", "RColorBrewer", "svglite"),
  install_missing = INSTALL_MISSING_PKGS
)

load_pkgs(
  c("DESeq2", "SummarizedExperiment", "edgeR", "matrixStats", "maSigPro",
    "AnnotationDbi", "clusterProfiler", "enrichplot"),
  install_missing = INSTALL_MISSING_PKGS,
  bioc = TRUE
)

##############################
## HELPER FUNCTIONS
##############################
read_featurecounts_matrix <- function(path) {
  fc <- read.delim(
    path, header = TRUE, sep = "\t",
    comment.char = "#", quote = "", check.names = FALSE,
    stringsAsFactors = FALSE
  )
  stopifnot("Geneid" %in% colnames(fc))
  drop <- intersect(colnames(fc), c("Chr", "Start", "End", "Strand", "Length"))
  mat <- fc |>
    dplyr::select(-dplyr::all_of(drop)) |>
    tibble::column_to_rownames("Geneid") |>
    as.matrix()
  mode(mat) <- "numeric"
  colnames(mat) <- sapply(colnames(mat), function(x) {
    m <- regmatches(x, regexpr("[0-9]{4}", x))
    if (length(m) == 1 && nchar(m) == 4) m else x
  })
  colnames(mat) <- make.unique(colnames(mat))
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
    x %in% c("EAU", "WATER", "H2O", "CTRL", "CONTROL") ~ "EAU",
    x %in% c("SDP") ~ "SDP",
    TRUE ~ x
  )
}

theme_pub <- function(base_size = 14) {
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

save_plot <- function(p, path_base, w = 8, h = 6, dpi = 600) {
  ggplot2::ggsave(paste0(path_base, ".pdf"), p, width = w, height = h, bg = "white")
  ggplot2::ggsave(paste0(path_base, ".png"), p, width = w, height = h, dpi = dpi, bg = "white")
  ggplot2::ggsave(paste0(path_base, ".svg"), p, width = w, height = h, bg = "white",
                  device = svglite::svglite)
}

save_pheatmap <- function(ph_obj, path_base, w = 11, h = 12) {
  grDevices::pdf(paste0(path_base, ".pdf"), width = w, height = h)
  grid::grid.draw(ph_obj$gtable)
  grDevices::dev.off()
  
  grDevices::png(paste0(path_base, ".png"), width = w, height = h, units = "in", res = 300, bg = "white")
  grid::grid.draw(ph_obj$gtable)
  grDevices::dev.off()
  
  svglite::svglite(paste0(path_base, ".svg"), width = w, height = h, bg = "white")
  grid::grid.draw(ph_obj$gtable)
  grDevices::dev.off()
}

extract_gene_ids <- function(x) {
  if (is.null(x)) return(character(0))
  x_df <- tryCatch(as.data.frame(x), error = function(e) NULL)
  if (is.null(x_df) || nrow(x_df) == 0) return(character(0))
  rn <- rownames(x_df)
  if (!is.null(rn) && length(rn) > 0) return(rn)
  for (col in c("Gene", "gene", "Geneid")) {
    if (col %in% colnames(x_df)) return(as.character(x_df[[col]]))
  }
  character(0)
}

##############################
## 1. IMPORT DATA
##############################
log_msg("========================================")
log_msg("maSigPro pipeline — Shared baseline D0")
log_msg("Genotype: ", GENOTYPE_CHOSEN)
log_msg("========================================")

log_msg("Importing counts matrix")
count_mat <- read_featurecounts_matrix(counts_file)
log_msg("Counts: ", nrow(count_mat), " genes x ", ncol(count_mat), " samples")

log_msg("Importing metadata")
meta_raw <- read.csv(meta_file, sep = ";", stringsAsFactors = FALSE)

required_cols <- c("Sample", "Genotype", "Temperature", "Treatment", "Sampling_time")
missing_cols <- setdiff(required_cols, colnames(meta_raw))
if (length(missing_cols) > 0) stop("Missing metadata columns: ", paste(missing_cols, collapse = ", "))

meta <- meta_raw |>
  dplyr::mutate(
    Sample        = sprintf("%04d", as.integer(as.character(Sample))),
    Genotype      = trimws(Genotype),
    Temperature   = normalize_temperature(Temperature),
    Treatment     = normalize_treatment(Treatment),
    Sampling_time = trimws(Sampling_time)
  ) |>
  dplyr::filter(
    Genotype == GENOTYPE_CHOSEN,
    Sampling_time %in% c("P1", "P2", "P3"),
    Sample %in% colnames(count_mat)
  ) |>
  dplyr::mutate(
    Temperature   = factor(Temperature, levels = c("T0", "T1", "T2")),
    Treatment     = factor(Treatment, levels = c("EAU", "SDP")),
    Sampling_time = factor(Sampling_time, levels = c("P1", "P2", "P3"))
  ) |>
  droplevels()

# ── KEY: P1 has only Control (EAU). Remove any P1 SDP if present ──
meta <- meta |>
  dplyr::filter(!(Sampling_time == "P1" & Treatment == "SDP"))

# Numeric time for maSigPro
meta$Time <- c(P1 = 0, P2 = 2, P3 = 4)[as.character(meta$Sampling_time)]
meta$Condition <- paste(meta$Temperature, meta$Treatment, sep = "_")

# Display columns
meta$Temp_disp <- relabel_temp(meta$Temperature)
meta$Trt_disp  <- relabel_trt(meta$Treatment)
meta$Time_disp <- relabel_time(meta$Sampling_time)

stopifnot(nrow(meta) >= 15)

write.csv(meta, file.path(dirs$inputs, "metadata_filtered.csv"), row.names = FALSE)

cat("\n", GENOTYPE_CHOSEN, " | D0+D2+D4 | Design table:\n")
print(table(
  Temperature = meta$Temp_disp,
  Treatment = meta$Trt_disp,
  Time = meta$Time_disp
))
cat("\n")

log_msg("Samples after filtering: ", nrow(meta))
log_msg("  D0 (P1, Control only): ", sum(meta$Sampling_time == "P1"))
log_msg("  D2 (P2): ", sum(meta$Sampling_time == "P2"))
log_msg("  D4 (P3): ", sum(meta$Sampling_time == "P3"))

##############################
## 2. SUBSET COUNT MATRIX + FILTERING
##############################
count_sub <- count_mat[, meta$Sample, drop = FALSE]
stopifnot(all(colnames(count_sub) == meta$Sample))

keep1 <- rowSums(count_sub, na.rm = TRUE) >= MIN_TOTAL_COUNTS
keep2 <- rowSums(count_sub >= MIN_COUNT, na.rm = TRUE) >= MIN_SAMPLES
count_f <- count_sub[keep1 & keep2, , drop = FALSE]

log_msg("Count filtering: ", nrow(count_mat), " -> ", nrow(count_f), " genes")

##############################
## 3. NORMALIZATION
##############################
log_msg("Normalization (VST + logCPM)")

dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = round(count_f),
  colData   = meta,
  design    = ~ Temperature + Sampling_time
)
dds <- DESeq2::estimateSizeFactors(dds)
vsd <- DESeq2::vst(dds, blind = FALSE)
expr_vst <- SummarizedExperiment::assay(vsd)

dge <- edgeR::DGEList(counts = count_f)
dge <- edgeR::calcNormFactors(dge, method = "TMM")
log_cpm <- edgeR::cpm(dge, log = TRUE, prior.count = 1)

saveRDS(expr_vst, file.path(dirs$inputs, "expr_vst_before_MAD.rds"))
saveRDS(log_cpm,  file.path(dirs$inputs, "log_cpm_before_MAD.rds"))

log_msg("Normalization OK: ", nrow(expr_vst), " genes")

##############################
## 3bis. MAD FILTERING
##############################
log_msg("MAD filtering")

n_before_mad <- nrow(expr_vst)
mad_vec  <- matrixStats::rowMads(expr_vst)
good_mad <- !is.na(mad_vec) & mad_vec > 0
mad_thr  <- quantile(mad_vec[good_mad], MAD_QUANTILE, na.rm = TRUE)
keep_mad <- good_mad & mad_vec >= mad_thr

expr_vst <- expr_vst[keep_mad, , drop = FALSE]
log_cpm  <- log_cpm[rownames(expr_vst), , drop = FALSE]

write.csv(
  data.frame(Step = c("Before_MAD", "After_MAD"), N_genes = c(n_before_mad, nrow(expr_vst))),
  file.path(dirs$tables, "MAD_filter_summary.csv"), row.names = FALSE
)

saveRDS(expr_vst, file.path(dirs$inputs, "expr_vst_after_MAD.rds"))
saveRDS(log_cpm,  file.path(dirs$inputs, "log_cpm_after_MAD.rds"))

log_msg("MAD filter: ", n_before_mad, " -> ", nrow(expr_vst), " genes")

##############################
## 4. maSigPro DESIGN — SHARED BASELINE
##############################
log_msg("Building maSigPro design with shared D0 baseline")

# ── Strategy: D0 Control samples are duplicated for each Treated trajectory ──
# Original samples at D0 belong to condition T0_EAU, T1_EAU, T2_EAU
# We duplicate them as T0_SDP, T1_SDP, T2_SDP at Time=0

meta_d0 <- meta |> dplyr::filter(Sampling_time == "P1")
meta_d2d4 <- meta |> dplyr::filter(Sampling_time %in% c("P2", "P3"))

# Create duplicated D0 rows for Treated trajectories
meta_d0_treated <- meta_d0 |>
  dplyr::mutate(
    Treatment = factor("SDP", levels = c("EAU", "SDP")),
    Condition = paste(Temperature, "SDP", sep = "_"),
    Trt_disp  = "Treated",
    # Mark as duplicated baseline
    Sample_orig = Sample,
    Sample = paste0(Sample, "_SDP_baseline")
  )

# Update original D0 to mark them
meta_d0 <- meta_d0 |>
  dplyr::mutate(Sample_orig = Sample)

meta_d2d4 <- meta_d2d4 |>
  dplyr::mutate(Sample_orig = Sample)

meta_d0_treated <- meta_d0_treated |>
  dplyr::mutate(Sample_orig = Sample_orig)

# Combine all metadata rows
meta_ext <- dplyr::bind_rows(meta_d0, meta_d0_treated, meta_d2d4)
meta_ext$Replicate <- seq_len(nrow(meta_ext))

log_msg("Extended design: ", nrow(meta_ext), " rows (", nrow(meta_d0), " D0 Control + ",
        nrow(meta_d0_treated), " D0 shared baseline + ", nrow(meta_d2d4), " D2+D4)")

# ── Expression matrix: duplicate D0 columns for the baseline ──
expr_for_masigpro <- log_cpm[, meta$Sample, drop = FALSE]

# Add duplicated columns
d0_dupl_cols <- meta_d0_treated$Sample
d0_orig_cols <- meta_d0$Sample

expr_dupl <- expr_for_masigpro[, d0_orig_cols, drop = FALSE]
colnames(expr_dupl) <- d0_dupl_cols

expr_ext <- cbind(expr_for_masigpro, expr_dupl)
expr_ext <- expr_ext[, meta_ext$Sample_orig[!grepl("_SDP_baseline$", meta_ext$Sample)]]

# Actually, simpler approach: build expression matrix indexed by row number
# since maSigPro uses rownames of edesign to match columns
rownames(meta_ext) <- paste0("S", seq_len(nrow(meta_ext)))

expr_masigpro <- matrix(NA_real_, nrow = nrow(log_cpm), ncol = nrow(meta_ext))
rownames(expr_masigpro) <- rownames(log_cpm)
colnames(expr_masigpro) <- rownames(meta_ext)

for (i in seq_len(nrow(meta_ext))) {
  orig_sample <- meta_ext$Sample_orig[i]
  if (orig_sample %in% colnames(log_cpm)) {
    expr_masigpro[, i] <- log_cpm[, orig_sample]
  }
}

# Remove genes with NA
keep_complete <- complete.cases(expr_masigpro)
expr_masigpro <- expr_masigpro[keep_complete, , drop = FALSE]

log_msg("Expression matrix for maSigPro: ", nrow(expr_masigpro), " genes x ", ncol(expr_masigpro), " columns")

# ── Build edesign ──
# Conditions (reference = T0_EAU)
# Dummy variables for each non-reference condition
conditions <- unique(meta_ext$Condition)
ref_cond <- "T0_EAU"
other_conds <- setdiff(conditions, ref_cond)

edesign_df <- data.frame(
  Time      = meta_ext$Time,
  Replicate = meta_ext$Replicate,
  row.names = rownames(meta_ext)
)

for (cond in other_conds) {
  edesign_df[[cond]] <- as.integer(meta_ext$Condition == cond)
}

edesign <- as.matrix(edesign_df)

write.csv(
  as.data.frame(edesign),
  file.path(dirs$inputs, "edesign_maSigPro_sharedBaseline.csv"),
  row.names = TRUE
)

# Save metadata extension for reference
write.csv(meta_ext, file.path(dirs$inputs, "metadata_extended_sharedBaseline.csv"), row.names = TRUE)

log_msg("edesign: ", nrow(edesign), " samples x ", ncol(edesign), " variables")
log_msg("Conditions: ", paste(colnames(edesign), collapse = ", "))

##############################
## 5. maSigPro
##############################
log_msg("Running maSigPro::make.design.matrix")
NBdesign <- maSigPro::make.design.matrix(
  edesign = edesign,
  degree  = MASIGPRO_DEGREE
)

log_msg("Running maSigPro::p.vector")
pvect <- maSigPro::p.vector(
  data      = expr_masigpro,
  design    = NBdesign,
  counts    = FALSE,
  Q         = MASIGPRO_Q,
  MT.adjust = "BH",
  min.obs   = ceiling(ncol(expr_masigpro) * 0.4)
)

saveRDS(pvect, file.path(dirs$masigpro, "pvect.rds"))
log_msg("p.vector: ", pvect$i, " initial significant genes")

log_msg("Running maSigPro::T.fit")
tstep <- maSigPro::T.fit(
  data        = pvect,
  step.method = MASIGPRO_STEP,
  alfa        = MASIGPRO_Q,
  nvar        = MASIGPRO_NVAR
)

saveRDS(tstep, file.path(dirs$masigpro, "tstep.rds"))
log_msg("T.fit completed")

##############################
## 6. EXTRACT SIGNIFICANT GENES
##############################
log_msg("Extracting significant genes")

sig_genes <- extract_gene_ids(pvect$SELEC)
if (length(sig_genes) == 0) sig_genes <- extract_gene_ids(pvect$Q)
if (length(sig_genes) == 0) sig_genes <- extract_gene_ids(tstep$sol)

sig_genes <- unique(sig_genes[!is.na(sig_genes) & sig_genes != ""])
sig_genes <- intersect(sig_genes, rownames(expr_vst))

log_msg("Temporal genes for clustering: ", length(sig_genes))

if (length(sig_genes) < 2) stop("Too few significant genes for clustering.")

write.csv(
  data.frame(Gene = sig_genes),
  file.path(dirs$tables, "maSigPro_significant_genes.csv"),
  row.names = FALSE
)

# Strict genes from T.fit
strict_genes <- character(0)
if (!is.null(tstep$sol)) {
  strict_genes <- extract_gene_ids(tstep$sol)
  strict_genes <- strict_genes[!is.na(strict_genes) & strict_genes != ""]
  strict_genes <- intersect(strict_genes, rownames(expr_vst))
  write.csv(
    data.frame(Gene = strict_genes),
    file.path(dirs$tables, "maSigPro_strict_genes_Tfit.csv"),
    row.names = FALSE
  )
  log_msg("Strict T.fit genes: ", length(strict_genes))
}

##############################
## 7. MEAN PROFILE MATRIX
##############################
log_msg("Computing mean expression profiles")

# 18 groups: 3 Temp x 2 Trt x 3 Time
# But D0 has no Treated → we use Control values as shared baseline
# Display labels for column names
group_levels <- c(
  "NH_Control_D0", "NH_Control_D2", "NH_Control_D4",
  "NH_Treated_D0", "NH_Treated_D2", "NH_Treated_D4",
  "HS-7_Control_D0", "HS-7_Control_D2", "HS-7_Control_D4",
  "HS-7_Treated_D0", "HS-7_Treated_D2", "HS-7_Treated_D4",
  "HS-2_Control_D0", "HS-2_Control_D2", "HS-2_Control_D4",
  "HS-2_Treated_D0", "HS-2_Treated_D2", "HS-2_Treated_D4"
)

# Build mapping: display group -> list of sample IDs
# Use original meta (not extended) for actual expression values
group_map_display <- paste(
  relabel_temp(meta$Temperature),
  relabel_trt(meta$Treatment),
  relabel_time(meta$Sampling_time),
  sep = "_"
)
names(group_map_display) <- meta$Sample

# For Treated_D0: use Control_D0 samples of same temperature
profile_mat <- sapply(group_levels, function(grp) {
  # Parse group
  parts <- strsplit(grp, "_")[[1]]
  temp_d <- parts[1]
  trt_d  <- parts[2]
  time_d <- parts[3]
  
  # Find matching samples
  samps <- names(group_map_display)[group_map_display == grp]
  
  # If Treated at D0: use Control at D0 for same temperature (shared baseline)
  if (trt_d == "Treated" && time_d == "D0") {
    baseline_grp <- paste(temp_d, "Control", "D0", sep = "_")
    samps <- names(group_map_display)[group_map_display == baseline_grp]
  }
  
  if (length(samps) == 0) {
    rep(NA_real_, length(sig_genes))
  } else if (length(samps) == 1) {
    expr_vst[sig_genes, samps]
  } else {
    rowMeans(expr_vst[sig_genes, samps, drop = FALSE], na.rm = TRUE)
  }
})

profile_mat <- as.matrix(profile_mat)
rownames(profile_mat) <- sig_genes

# Remove genes with any NA
keep_complete <- complete.cases(profile_mat)
profile_mat <- profile_mat[keep_complete, , drop = FALSE]

# Z-score
profile_z <- t(scale(t(profile_mat)))
profile_z[is.na(profile_z)] <- 0

# Clamp z-scores
profile_z[profile_z >  2] <-  2
profile_z[profile_z < -2] <- -2

saveRDS(profile_mat, file.path(dirs$clusters, "profile_mean_matrix.rds"))
saveRDS(profile_z,   file.path(dirs$clusters, "profile_mean_matrix_zscore.rds"))

log_msg("Profile matrix: ", nrow(profile_z), " genes x ", ncol(profile_z), " groups")

##############################
## 8. CLUSTERING
##############################
log_msg("Hierarchical clustering")

dist_mat <- dist(profile_z)
hc <- hclust(dist_mat, method = HCLUST_METHOD)
cluster_assign <- cutree(hc, k = K_CLUSTERS)

cluster_df <- data.frame(
  Gene    = names(cluster_assign),
  Cluster = as.integer(cluster_assign),
  row.names = NULL
) |> dplyr::arrange(Cluster, Gene)

write.csv(cluster_df, file.path(dirs$tables, "Gene_cluster_assignment.csv"), row.names = FALSE)

cluster_summary <- cluster_df |>
  dplyr::count(Cluster, name = "N_genes") |>
  dplyr::arrange(Cluster)

write.csv(cluster_summary, file.path(dirs$tables, "Cluster_summary.csv"), row.names = FALSE)

log_msg("Clustering: ", K_CLUSTERS, " clusters, ", nrow(profile_z), " genes")
for (i in seq_len(nrow(cluster_summary))) {
  log_msg("  Cluster ", cluster_summary$Cluster[i], ": ", cluster_summary$N_genes[i], " genes")
}

##############################
## 9. HEATMAP
##############################
log_msg("Generating heatmap")

ann_row <- data.frame(Cluster = factor(cluster_assign))
rownames(ann_row) <- names(cluster_assign)

cluster_colors <- setNames(
  colorRampPalette(brewer.pal(8, "Set2"))(K_CLUSTERS),
  as.character(seq_len(K_CLUSTERS))
)

# Column annotations
temp_vec <- sub("_.*", "", colnames(profile_z))
trt_vec  <- sub("^[^_]+_([^_]+)_.*$", "\\1", colnames(profile_z))
time_vec <- sub(".*_", "", colnames(profile_z))

ann_col <- data.frame(
  Temperature = factor(temp_vec, levels = c("NH", "HS-7", "HS-2")),
  Treatment   = factor(trt_vec, levels = c("Control", "Treated")),
  Time        = factor(time_vec, levels = c("D0", "D2", "D4"))
)
rownames(ann_col) <- colnames(profile_z)

hm_ann_colors <- list(
  Cluster     = cluster_colors,
  Temperature = TEMP_COLORS_D,
  Treatment   = TRT_COLORS_D,
  Time        = TIME_COLORS_D
)

hm_breaks <- seq(-2, 2, length.out = 101)
hm_colors <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100)

ph <- pheatmap::pheatmap(
  profile_z,
  cluster_rows    = hc,
  cluster_cols    = FALSE,
  annotation_row  = ann_row,
  annotation_col  = ann_col,
  annotation_colors = hm_ann_colors,
  show_rownames   = FALSE,
  fontsize_col    = 9,
  fontsize         = 10,
  breaks          = hm_breaks,
  color           = hm_colors,
  main = paste0(GENOTYPE_CHOSEN, " | D0→D2→D4 (shared baseline) | ",
                nrow(profile_z), " temporal genes"),
  silent = TRUE
)

save_pheatmap(ph, file.path(dirs$figures, "Heatmap_clusters"), w = 13, h = 13)

##############################
## 10. MEAN CLUSTER PROFILES
##############################
log_msg("Plotting cluster mean profiles")

cluster_labels <- setNames(
  paste0("Cluster ", cluster_summary$Cluster, " (n=", cluster_summary$N_genes, ")"),
  cluster_summary$Cluster
)

profile_long <- as.data.frame(profile_mat) |>
  tibble::rownames_to_column("Gene") |>
  tidyr::pivot_longer(-Gene, names_to = "Group", values_to = "Expression") |>
  dplyr::left_join(cluster_df, by = "Gene") |>
  dplyr::mutate(
    Temperature   = factor(sub("_.*", "", Group), levels = c("NH", "HS-7", "HS-2")),
    Treatment     = factor(sub("^[^_]+_([^_]+)_.*$", "\\1", Group), levels = c("Control", "Treated")),
    Time          = factor(sub(".*_", "", Group), levels = c("D0", "D2", "D4")),
    Cluster_label = factor(cluster_labels[as.character(Cluster)], levels = cluster_labels)
  )

cluster_mean <- profile_long |>
  dplyr::group_by(Cluster, Cluster_label, Temperature, Treatment, Time) |>
  dplyr::summarise(
    Mean = mean(Expression, na.rm = TRUE),
    SE   = sd(Expression, na.rm = TRUE) / sqrt(dplyr::n()),
    .groups = "drop"
  )

p_profiles <- ggplot2::ggplot(
  cluster_mean,
  ggplot2::aes(
    x = Time, y = Mean,
    group = interaction(Temperature, Treatment),
    colour = Temperature, linetype = Treatment
  )
) +
  ggplot2::geom_line(linewidth = 1) +
  ggplot2::geom_point(size = 2.5) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = Mean - SE, ymax = Mean + SE),
    width = 0.08, linewidth = 0.5
  ) +
  ggplot2::facet_wrap(~ Cluster_label, scales = "free_y") +
  ggplot2::scale_colour_manual(values = TEMP_COLORS_D) +
  ggplot2::scale_linetype_manual(values = c("Control" = "dashed", "Treated" = "solid")) +
  ggplot2::labs(
    title    = "Mean expression profiles by cluster",
    subtitle = paste0(GENOTYPE_CHOSEN, " | D0→D2→D4 (shared baseline) | ",
                      nrow(profile_z), " genes | Solid=Treated, Dashed=Control"),
    x = "Time point", y = "Mean VST expression"
  ) +
  theme_pub(14)

save_plot(p_profiles, file.path(dirs$figures, "Cluster_mean_profiles"), w = 15, h = 10)

##############################
## 10bis. INDIVIDUAL GENE PROFILES (z-score)
##############################
log_msg("Plotting individual gene profiles")

profile_gene_z <- as.data.frame(profile_z) |>
  tibble::rownames_to_column("Gene") |>
  tidyr::pivot_longer(-Gene, names_to = "Group", values_to = "Zscore") |>
  dplyr::left_join(cluster_df, by = "Gene") |>
  dplyr::mutate(
    Group = factor(Group, levels = group_levels),
    Cluster_label = factor(cluster_labels[as.character(Cluster)], levels = cluster_labels)
  )

p_gene_z <- ggplot2::ggplot(
  profile_gene_z,
  ggplot2::aes(x = Group, y = Zscore, group = Gene, colour = Gene)
) +
  ggplot2::geom_line(linewidth = 0.35, alpha = 0.7, show.legend = FALSE) +
  ggplot2::facet_wrap(~ Cluster_label, scales = "free_y", ncol = 3) +
  ggplot2::labs(
    title    = "Individual z-scored gene profiles by cluster",
    subtitle = paste0(GENOTYPE_CHOSEN, " | D0→D2→D4"),
    x = NULL, y = "Z-score expression"
  ) +
  theme_pub(10) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6)
  )

save_plot(p_gene_z, file.path(dirs$figures, "Cluster_individual_profiles_zscore"), w = 15, h = 11)

##############################
## 11. CLUSTER SIZE BARPLOT
##############################
p_sizes <- ggplot2::ggplot(cluster_summary, ggplot2::aes(x = factor(Cluster), y = N_genes)) +
  ggplot2::geom_col(width = 0.7, fill = "steelblue", alpha = 0.9, colour = "white") +
  ggplot2::geom_text(ggplot2::aes(label = N_genes), vjust = -0.3, size = 5, fontface = "bold") +
  ggplot2::labs(
    title    = "Number of genes per cluster",
    subtitle = paste0(GENOTYPE_CHOSEN, " | D0→D2→D4 | total = ", sum(cluster_summary$N_genes), " genes"),
    x = "Cluster", y = "Number of genes"
  ) +
  theme_pub(14)

save_plot(p_sizes, file.path(dirs$figures, "Cluster_sizes"), w = 9, h = 6)

##############################
## 12. TOP VARIABLE GENES PER CLUSTER
##############################
gene_var <- data.frame(
  Gene     = rownames(profile_mat),
  Variance = matrixStats::rowVars(profile_mat)
)

top_genes <- cluster_df |>
  dplyr::left_join(gene_var, by = "Gene") |>
  dplyr::group_by(Cluster) |>
  dplyr::arrange(dplyr::desc(Variance), .by_group = TRUE) |>
  dplyr::slice_head(n = 10) |>
  dplyr::ungroup()

write.csv(top_genes, file.path(dirs$tables, "Top10_variable_genes_per_cluster.csv"), row.names = FALSE)

##############################
## 13. GO ORA + PATHWAY ORA PER CLUSTER
##############################
if (RUN_ENRICHMENT) {
  log_msg("Starting GO / Pathway enrichment per cluster")
  
  # Install/load OrgDb
  if (!requireNamespace(ORGDB_PKG, quietly = TRUE)) {
    if (!INSTALL_MISSING_PKGS) stop("Missing OrgDb: ", ORGDB_PKG)
    install.packages(ORGDB_URL, repos = NULL, type = "source")
  }
  suppressPackageStartupMessages(library(ORGDB_PKG, character.only = TRUE))
  
  # Find OrgDb object
  OrgDb_obj <- NULL
  ns <- asNamespace(ORGDB_PKG)
  for (o in ls(ns, all.names = TRUE)) {
    x <- try(get(o, envir = ns), silent = TRUE)
    if (!inherits(x, "try-error") && inherits(x, "OrgDb")) { OrgDb_obj <- x; break }
  }
  if (is.null(OrgDb_obj)) stop("OrgDb object not found in: ", ORGDB_PKG)
  
  # Mapping
  stopifnot(file.exists(map_file))
  map_df <- read.csv(map_file) |>
    dplyr::select(my_gene_id, plant2t_id) |>
    dplyr::mutate(
      my_gene_id = trimws(as.character(my_gene_id)),
      plant2t_id = trimws(as.character(plant2t_id))
    ) |>
    dplyr::filter(!is.na(my_gene_id), my_gene_id != "",
                  !is.na(plant2t_id), plant2t_id != "") |>
    dplyr::distinct()
  
  write.csv(map_df, file.path(dirs$enrichment, "GeneID_mapping_used.csv"), row.names = FALSE)
  
  # Universe
  universe_ids <- data.frame(my_gene_id = rownames(expr_vst)) |>
    dplyr::left_join(map_df, by = "my_gene_id") |>
    dplyr::pull(plant2t_id) |>
    unique() |> na.omit()
  
  log_msg("Mapped universe: ", length(universe_ids), " PlanT2T IDs")
  
  # Pathway TERM2GENE
  p2g <- AnnotationDbi::select(
    OrgDb_obj, keys = universe_ids,
    columns = c("Pathway", "GID"), keytype = "GID"
  ) |>
    dplyr::filter(!is.na(Pathway), Pathway != "", !is.na(GID)) |>
    dplyr::distinct(Pathway, GID)
  
  # Pathway TERM2NAME
  ko_file <- system.file("extdata", "ko00001.PlanT2T.txt", package = ORGDB_PKG)
  if (file.exists(ko_file)) {
    ko_raw <- utils::read.delim(ko_file, header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)
    id_col <- intersect(c("pathway_id", "Pathway"), colnames(ko_raw))[1]
    nm_col <- intersect(c("level3", "Name", "name"), colnames(ko_raw))[1]
    p2n <- ko_raw |>
      dplyr::transmute(
        Pathway = trimws(as.character(.data[[id_col]])),
        Name    = trimws(as.character(.data[[nm_col]]))
      ) |>
      dplyr::filter(Pathway != "", Name != "") |>
      dplyr::distinct()
  } else {
    p2n <- data.frame(Pathway = unique(p2g$Pathway), Name = unique(p2g$Pathway))
  }
  
  write.csv(p2g, file.path(dirs$enrichment, "TERM2GENE_Pathway.csv"), row.names = FALSE)
  write.csv(p2n, file.path(dirs$enrichment, "TERM2NAME_Pathway.csv"), row.names = FALSE)
  
  all_go   <- data.frame()
  all_path <- data.frame()
  enrich_summary <- data.frame()
  
  for (cl in sort(unique(cluster_df$Cluster))) {
    cluster_label <- paste0("Cluster_", cl)
    dir_cl <- file.path(dirs$enrichment, cluster_label)
    dir.create(dir_cl, recursive = TRUE, showWarnings = FALSE)
    
    genes_cl <- cluster_df$Gene[cluster_df$Cluster == cl]
    
    mapped_ids <- data.frame(my_gene_id = genes_cl) |>
      dplyr::left_join(map_df, by = "my_gene_id") |>
      dplyr::pull(plant2t_id) |>
      unique() |> na.omit()
    
    enrich_summary <- dplyr::bind_rows(enrich_summary, data.frame(
      Cluster = cl, N_genes = length(genes_cl), N_mapped = length(mapped_ids)
    ))
    
    write.csv(data.frame(Gene = genes_cl),
              file.path(dir_cl, paste0(cluster_label, "_genes.csv")), row.names = FALSE)
    write.csv(data.frame(Plant2T_ID = mapped_ids),
              file.path(dir_cl, paste0(cluster_label, "_mapped_ids.csv")), row.names = FALSE)
    
    if (length(mapped_ids) < MIN_GENES_ORA) {
      log_msg("Cluster ", cl, " skipped for enrichment (", length(mapped_ids), " mapped)")
      next
    }
    
    ## GO ORA
    for (ont in GO_ONT_LIST) {
      ego <- tryCatch(
        suppressMessages(clusterProfiler::enrichGO(
          gene = mapped_ids, universe = universe_ids,
          OrgDb = OrgDb_obj, keyType = "GID", ont = ont,
          pAdjustMethod = "BH", pvalueCutoff = FDR_CUTOFF, qvalueCutoff = FDR_CUTOFF
        )),
        error = function(e) { log_msg("GO ", ont, " error cluster ", cl, ": ", conditionMessage(e)); NULL }
      )
      if (is.null(ego)) next
      
      ego_df <- as.data.frame(ego)
      if (nrow(ego_df) > 0) {
        ego_df$Cluster <- cl; ego_df$Ontology <- ont
        all_go <- dplyr::bind_rows(all_go, ego_df)
        write.csv(ego_df, file.path(dir_cl, paste0("GO_ORA_", ont, ".csv")), row.names = FALSE)
        
        p_go <- enrichplot::dotplot(ego, showCategory = GO_SHOWCATEGORY) +
          ggplot2::labs(title = paste0("Cluster ", cl, " — GO ", ont),
                        subtitle = paste0(GENOTYPE_CHOSEN, " | mapped = ", length(mapped_ids))) +
          theme_pub(10)
        save_plot(p_go, file.path(dir_cl, paste0("Dotplot_GO_", ont)),
                  w = 9, h = max(4, min(nrow(ego_df), GO_SHOWCATEGORY) * 0.25 + 2))
      }
    }
    
    ## Pathway ORA
    enr <- tryCatch(
      suppressMessages(clusterProfiler::enricher(
        gene = mapped_ids, universe = universe_ids,
        TERM2GENE = p2g, TERM2NAME = p2n,
        pAdjustMethod = "BH", pvalueCutoff = FDR_CUTOFF, qvalueCutoff = FDR_CUTOFF,
        minGSSize = MIN_GS_SIZE_PATH
      )),
      error = function(e) { log_msg("Pathway error cluster ", cl, ": ", conditionMessage(e)); NULL }
    )
    if (!is.null(enr)) {
      enr_df <- as.data.frame(enr)
      if (nrow(enr_df) > 0) {
        enr_df$Cluster <- cl
        all_path <- dplyr::bind_rows(all_path, enr_df)
        write.csv(enr_df, file.path(dir_cl, "Pathway_ORA.csv"), row.names = FALSE)
        
        p_pw <- enrichplot::dotplot(enr, showCategory = GO_SHOWCATEGORY) +
          ggplot2::labs(title = paste0("Cluster ", cl, " — Pathway ORA"),
                        subtitle = paste0(GENOTYPE_CHOSEN, " | mapped = ", length(mapped_ids))) +
          theme_pub(10)
        save_plot(p_pw, file.path(dir_cl, "Dotplot_Pathway"),
                  w = 9, h = max(4, min(nrow(enr_df), GO_SHOWCATEGORY) * 0.25 + 2))
      }
    }
    
    log_msg("Cluster ", cl, " | genes=", length(genes_cl), " | mapped=", length(mapped_ids),
            " | GO=", nrow(all_go |> dplyr::filter(Cluster == cl)),
            " | Pathways=", nrow(all_path |> dplyr::filter(Cluster == cl)))
  }
  
  write.csv(enrich_summary, file.path(dirs$enrichment, "Enrichment_summary.csv"), row.names = FALSE)
  if (nrow(all_go) > 0)   write.csv(all_go,   file.path(dirs$enrichment, "GO_ALL_combined.csv"), row.names = FALSE)
  if (nrow(all_path) > 0) write.csv(all_path, file.path(dirs$enrichment, "Pathway_ALL_combined.csv"), row.names = FALSE)
  
  log_msg("Enrichment completed")
}

##############################
## 14. EXPORT FOR INTERSECTION SCRIPT
##############################
# The Gene_cluster_assignment.csv is already compatible with the intersection script.
# It has columns: Gene, Cluster — same format expected by Intersection_maSigPro_WGCNA.R

log_msg("Gene_cluster_assignment.csv exported — compatible with intersection script")

##############################
## 15. FINAL SUMMARY
##############################
summary_df <- data.frame(
  Parameter = c(
    "Genotype", "Time_points", "Design",
    "Genes_after_MAD", "Genes_for_clustering", "Clustered_genes",
    "Strict_genes_Tfit", "K_clusters", "Degree", "Q",
    "Shared_baseline", "RUN_ENRICHMENT"
  ),
  Value = c(
    GENOTYPE_CHOSEN, "D0_D2_D4",
    "D0 Control shared as baseline for Treated trajectories",
    nrow(expr_vst), length(sig_genes), nrow(profile_z),
    length(strict_genes), K_CLUSTERS, MASIGPRO_DEGREE, MASIGPRO_Q,
    "YES — D0 Control duplicated for Treated conditions",
    RUN_ENRICHMENT
  )
)

write.csv(summary_df, file.path(dirs$tables, "Analysis_summary.csv"), row.names = FALSE)

log_msg("================================================")
log_msg("Pipeline completed")
log_msg("Genotype: ", GENOTYPE_CHOSEN)
log_msg("Time points: D0 (shared baseline) → D2 → D4")
log_msg("Genes clustered: ", nrow(profile_z))
log_msg("Clusters: ", K_CLUSTERS)
log_msg("Output: ", base_out)
log_msg("================================================")

cat("\n", strrep("=", 60), "\n")
cat("PIPELINE COMPLETED\n")
cat(strrep("=", 60), "\n")
cat("Genotype:      ", GENOTYPE_CHOSEN, "\n")
cat("Time points:   D0 (shared baseline) → D2 → D4\n")
cat("Genes:         ", nrow(profile_z), " clustered into ", K_CLUSTERS, " clusters\n")
cat("Output:        ", base_out, "\n")
cat("Intersection:  Gene_cluster_assignment.csv is ready\n")
cat(strrep("=", 60), "\n")