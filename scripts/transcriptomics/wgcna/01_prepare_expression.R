############################################################
## 01_import_filter_normalize.R
## Import counts + metadata → filtre génotype + P2+P3
## → DESeq2 normalisation + VST → filtrage MAD → datExpr
# Repository configuration (portable; no HPC-specific paths)
if (!exists('CARROT_CONFIG_LOADED')) source('../../../config/config.R')


############################################################

##############################
## PARAMÈTRES (modifier ici)
##############################
GENOTYPE_CHOSEN <- Sys.getenv("CARROT_GENOTYPE", unset = "PRESTO")
SEED <- 123
INSTALL_MISSING_PKGS <- FALSE

counts_file <- RNA_COUNTS_FILE
meta_file <- RNA_METADATA_FILE
base_out <- RNA_WGCNA_ROOT

# Filtrage
MIN_TOTAL_COUNTS <- 10
MIN_COUNT        <- 10
MIN_SAMPLES      <- 3
MAD_QUANTILE     <- 0.50  

set.seed(SEED)
options(stringsAsFactors = FALSE)

##############################
## SOURCE HELPERS
##############################
source("00_helpers.R")
load_pkgs(c("dplyr","tibble"), install_missing = INSTALL_MISSING_PKGS)
load_pkgs(c("DESeq2","SummarizedExperiment","matrixStats"), install_missing = INSTALL_MISSING_PKGS, bioc = TRUE)

GENO_TAG <- gsub("[^A-Za-z0-9]+", "", GENOTYPE_CHOSEN)
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_dir <- file.path(base_out, paste0("WGCNA_", GENO_TAG, "_P2P3_", timestamp))
dir_inputs <- file.path(out_dir, "00_inputs")
dir.create(dir_inputs, recursive = TRUE, showWarnings = FALSE)

log_msg <- init_logger(file.path(out_dir, "pipeline.log"))
log_msg("START Script 01 — ", GENOTYPE_CHOSEN)

##############################
## PARTIE 1: Import counts
##############################
cat("── IMPORT COUNTS ──\n")
count_mat <- read_featurecounts_matrix(counts_file)
log_msg("Counts: ", nrow(count_mat), " genes x ", ncol(count_mat), " samples")

##############################
## PARTIE 2: Import + clean metadata
##############################
cat("── IMPORT METADATA ──\n")
meta_raw <- utils::read.csv(meta_file, sep = ";", stringsAsFactors = FALSE)

meta <- meta_raw |>
  dplyr::mutate(
    Sample        = sprintf("%04d", as.integer(as.character(Sample))),
    Genotype      = trimws(Genotype),
    Sampling_time = trimws(Sampling_time),
    Temperature   = normalize_temperature(Temperature),
    Treatment     = normalize_treatment(Treatment)
  ) |>
  # ★ Filtrer: génotype choisi + P2 et P3 seulement
  dplyr::filter(Genotype == GENOTYPE_CHOSEN) |>
  dplyr::filter(Sampling_time %in% c("P2", "P3")) |>
  dplyr::filter(Sample %in% colnames(count_mat)) |>
  dplyr::mutate(
    Genotype      = factor(Genotype),
    Temperature   = factor(Temperature, levels = c("T0","T1","T2")),
    Treatment     = factor(Treatment, levels = c("Control","Treated")),
    Sampling_time = factor(Sampling_time, levels = c("P2","P3")),
    # Groupe complet pour les boxplots
    Group = paste(Temperature, Treatment, Sampling_time, sep = "_")
  ) |>
  droplevels()

stopifnot(nrow(meta) >= 10)
log_msg("Samples après filtrage (", GENOTYPE_CHOSEN, " D2+D4): ", nrow(meta))

cat("\nDesign table:\n")
print(table(meta$Temperature, meta$Treatment, meta$Sampling_time))
cat("\n")

##############################
## PARTIE 3: Aligner counts ↔ metadata
##############################
count_tp <- count_mat[, meta$Sample, drop = FALSE]
stopifnot(all(colnames(count_tp) == meta$Sample))
log_msg("Aligned: ", nrow(count_tp), " genes x ", ncol(count_tp), " samples")

##############################
## PARTIE 4: Pre-filtrage counts
##############################
cat("── PRE-FILTRAGE ──\n")
keep1 <- rowSums(count_tp, na.rm = TRUE) >= MIN_TOTAL_COUNTS
keep2 <- rowSums(count_tp >= MIN_COUNT) >= MIN_SAMPLES
count_f <- count_tp[keep1 & keep2, , drop = FALSE]
log_msg("Pre-filtre: ", nrow(count_tp), " → ", nrow(count_f), " genes")

##############################
## PARTIE 5: DESeq2 normalisation + VST
##############################
cat("── DESeq2 VST ──\n")
dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = round(count_f),
  colData   = meta,
  design    = ~ Temperature + Treatment + Sampling_time
)

# Filtre DESeq2-style
keep_deseq <- rowSums(DESeq2::counts(dds) >= MIN_COUNT) >= MIN_SAMPLES
dds <- dds[keep_deseq, ]
log_msg("DESeq2 filter: ", sum(keep_deseq), " genes kept")

dds <- DESeq2::estimateSizeFactors(dds)
vsd <- DESeq2::vst(dds, blind = FALSE)
expr_mat <- SummarizedExperiment::assay(vsd)

log_msg("VST done: ", nrow(expr_mat), " genes x ", ncol(expr_mat), " samples")

##############################
## PARTIE 6: Filtrage MAD (top 40%)
##############################
cat("── FILTRAGE MAD ──\n")
n_before <- nrow(expr_mat)
mad_vec <- matrixStats::rowMads(expr_mat)
names(mad_vec) <- rownames(expr_mat)

# Retirer gènes constants
var_vec <- matrixStats::rowVars(expr_mat)
good <- !is.na(mad_vec) & mad_vec > 0 & var_vec > 0
expr_mat <- expr_mat[good, , drop = FALSE]
mad_vec <- mad_vec[good]

# Garder top (1 - MAD_QUANTILE)
mad_thr <- quantile(mad_vec, MAD_QUANTILE)
expr_mat <- expr_mat[mad_vec >= mad_thr, , drop = FALSE]
n_after <- nrow(expr_mat)

log_msg("MAD filter (top ", round((1 - MAD_QUANTILE) * 100), "%): ",
        n_before, " → ", n_after, " genes")

# datExpr: samples × genes (format WGCNA)
datExpr <- as.data.frame(t(expr_mat))
stopifnot(rownames(datExpr) == meta$Sample)

log_msg("datExpr ready: ", nrow(datExpr), " samples x ", ncol(datExpr), " genes")

##############################
## PARTIE 7: Sauvegarder
##############################
cat("── EXPORT ──\n")
saveRDS(datExpr,  file.path(out_dir, "datExpr.rds"))
saveRDS(meta,     file.path(out_dir, "meta.rds"))
saveRDS(vsd,      file.path(out_dir, "vsd.rds"))
saveRDS(expr_mat, file.path(out_dir, "expr_mat_VST_MAD.rds"))

write.csv(meta, file.path(dir_inputs, "metadata_aligned.csv"), row.names = FALSE)

# Résumé du filtrage
filter_summary <- data.frame(
  Step = c("Raw genes", "Pre-filter (counts)", "DESeq2 filter", "MAD filter"),
  N_genes = c(nrow(count_mat), nrow(count_f), sum(keep_deseq), n_after)
)
write.csv(filter_summary, file.path(dir_inputs, "gene_filtering_summary.csv"), row.names = FALSE)
print(filter_summary)

# Sauver les paramètres
params <- data.frame(
  Parameter = c("GENOTYPE_CHOSEN","SEED","MIN_TOTAL_COUNTS","MIN_COUNT",
                "MIN_SAMPLES","MAD_QUANTILE","N_samples","N_genes_final"),
  Value = c(GENOTYPE_CHOSEN, SEED, MIN_TOTAL_COUNTS, MIN_COUNT,
            MIN_SAMPLES, MAD_QUANTILE, nrow(meta), ncol(datExpr))
)
write.csv(params, file.path(dir_inputs, "parameters.csv"), row.names = FALSE)

log_msg("Script 01 DONE — output: ", out_dir)
cat("\n✓ Script 01 terminé\n")
cat("  out_dir:", out_dir, "\n")