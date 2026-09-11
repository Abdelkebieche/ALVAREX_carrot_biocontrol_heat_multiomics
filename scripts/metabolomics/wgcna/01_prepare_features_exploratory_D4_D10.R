############################################################
## 01_import_filter_normalize.R  (LC-MS WGCNA)
## Import intensités LC-MS + metadata → filtre génotype + P3+P4
## → log2 + Pareto scaling → filtrage MAD → datExpr
# Repository configuration (portable; no HPC-specific paths)
if (!exists('CARROT_CONFIG_LOADED')) source('../../../config/config.R')


############################################################
## DIFFÉRENCES vs RNAseq:
##  - Pas de DESeq2 / VST (données déjà continues)
##  - Normalisation: log2(x+1) + Pareto scaling
##  - P3 + P4 au lieu de P2 + P3
##  - Features = métabolites (m/z_RT) au lieu de gènes
############################################################

##############################
## PARAMÈTRES (modifier ici)
##############################
GENOTYPE_CHOSEN <- Sys.getenv("CARROT_GENOTYPE", unset = "ROBILA")
SEED <- 123
INSTALL_MISSING_PKGS <- FALSE

# ★ CHEMINS À ADAPTER
lcms_file <- LCMS_FILE
meta_file <- LCMS_METADATA_FILE
base_out <- LCMS_WGCNA_ROOT

# Filtrage
MIN_INTENSITY    <- 1000     # Intensité minimale pour considérer un signal
MIN_DETECT_FRAC  <- 0.50     # Feature détectée dans au moins 50% des échantillons
MAD_QUANTILE     <- 0.2     # Garder top 50% les plus variables (MAD)

# Normalisation
LOG_TRANSFORM    <- TRUE     # log2(x + 1)
PARETO_SCALE     <- TRUE     # Pareto scaling (diviser par sqrt(SD))

set.seed(SEED)
options(stringsAsFactors = FALSE)

##############################
## SOURCE HELPERS
##############################
source("00_helpers.R")
load_pkgs(c("dplyr","tibble","matrixStats"), install_missing = INSTALL_MISSING_PKGS)

GENO_TAG <- gsub("[^A-Za-z0-9]+", "", GENOTYPE_CHOSEN)
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_dir <- file.path(base_out, paste0("WGCNA_LCMS_", GENO_TAG, "_P3P4_", timestamp))
dir_inputs <- file.path(out_dir, "00_inputs")
dir.create(dir_inputs, recursive = TRUE, showWarnings = FALSE)

log_msg <- init_logger(file.path(out_dir, "pipeline.log"))
log_msg("START Script 01 — ", GENOTYPE_CHOSEN, " (LC-MS P3+P4)")

##############################
## PARTIE 1: Import LC-MS
##############################
cat("── IMPORT LC-MS ──\n")
lcms_data <- read_lcms_matrix(lcms_file)
intensity_mat <- lcms_data$mat
feature_info  <- lcms_data$feature_info

log_msg("LC-MS brut: ", nrow(intensity_mat), " features x ", ncol(intensity_mat), " échantillons")

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
  # ★ Filtrer: génotype choisi + P3 et P4 seulement
  dplyr::filter(Genotype == GENOTYPE_CHOSEN) |>
  dplyr::filter(Sampling_time %in% c("P3", "P4")) |>
  dplyr::filter(Sample %in% colnames(intensity_mat)) |>
  dplyr::mutate(
    Genotype      = factor(Genotype),
    Temperature   = factor(Temperature, levels = c("T0","T1","T2")),
    Treatment     = factor(Treatment, levels = c("Control","Treated")),
    Sampling_time = factor(Sampling_time, levels = c("P3","P4")),
    # Groupe complet pour les boxplots
    Group = paste(Temperature, Treatment, Sampling_time, sep = "_")
  ) |>
  droplevels()

stopifnot(nrow(meta) >= 10)
log_msg("Samples après filtrage (", GENOTYPE_CHOSEN, " P3+P4): ", nrow(meta))

cat("\nDesign table:\n")
print(table(meta$Temperature, meta$Treatment, meta$Sampling_time))
cat("\n")

##############################
## PARTIE 3: Aligner intensités ↔ metadata
##############################
common_samples <- intersect(colnames(intensity_mat), meta$Sample)
stopifnot(length(common_samples) >= 10)

intensity_tp <- intensity_mat[, common_samples, drop = FALSE]
meta <- meta |> dplyr::filter(Sample %in% common_samples)
# Réordonner
intensity_tp <- intensity_tp[, meta$Sample, drop = FALSE]
stopifnot(all(colnames(intensity_tp) == meta$Sample))
log_msg("Aligned: ", nrow(intensity_tp), " features x ", ncol(intensity_tp), " samples")

##############################
## PARTIE 4: Pre-filtrage features
##############################
cat("── PRE-FILTRAGE ──\n")
n_raw <- nrow(intensity_tp)

# 4a) Remplacer 0 et NA par NA pour le filtrage
int_work <- intensity_tp
int_work[int_work == 0] <- NA

# 4b) Filtre de détection: feature présente dans au moins MIN_DETECT_FRAC des échantillons
detect_rate <- rowMeans(!is.na(int_work) & int_work >= MIN_INTENSITY)
keep_detect <- detect_rate >= MIN_DETECT_FRAC
cat("  Détection (>=", MIN_DETECT_FRAC*100, "% samples):", sum(keep_detect), "/", n_raw, "features\n")

int_work <- int_work[keep_detect, , drop = FALSE]

# 4c) Remplacer les NA restants par la moitié du minimum de chaque feature (imputation simple)
for (i in seq_len(nrow(int_work))) {
  row_vals <- int_work[i, ]
  na_idx <- is.na(row_vals) | row_vals == 0
  if (any(na_idx)) {
    min_val <- min(row_vals[!na_idx], na.rm = TRUE)
    int_work[i, na_idx] <- min_val / 2
  }
}

n_after_detect <- nrow(int_work)
log_msg("Pre-filtre détection: ", n_raw, " → ", n_after_detect, " features")

##############################
## PARTIE 5: Normalisation LC-MS
##############################
cat("── NORMALISATION ──\n")

# 5a) Log2 transformation
if (LOG_TRANSFORM) {
  cat("  Log2(x + 1) transformation...\n")
  expr_mat <- log2(int_work + 1)
  log_msg("Log2 transformation done")
} else {
  expr_mat <- int_work
}

# 5b) Pareto scaling (par feature: centrer + diviser par sqrt(SD))
if (PARETO_SCALE) {
  cat("  Pareto scaling...\n")
  row_means <- rowMeans(expr_mat, na.rm = TRUE)
  row_sds   <- matrixStats::rowSds(expr_mat, na.rm = TRUE)
  
  # Éviter division par zéro
  row_sds[row_sds == 0] <- 1
  
  expr_mat <- (expr_mat - row_means) / sqrt(row_sds)
  log_msg("Pareto scaling done")
}

log_msg("Normalisation done: ", nrow(expr_mat), " features x ", ncol(expr_mat), " samples")

##############################
## PARTIE 6: Filtrage MAD (top variable)
##############################
cat("── FILTRAGE MAD ──\n")
n_before <- nrow(expr_mat)
mad_vec <- matrixStats::rowMads(expr_mat)
names(mad_vec) <- rownames(expr_mat)

# Retirer features constantes
var_vec <- matrixStats::rowVars(expr_mat)
good <- !is.na(mad_vec) & mad_vec > 0 & var_vec > 0
expr_mat <- expr_mat[good, , drop = FALSE]
mad_vec <- mad_vec[good]

# Garder top (1 - MAD_QUANTILE)
mad_thr <- quantile(mad_vec, MAD_QUANTILE)
expr_mat <- expr_mat[mad_vec >= mad_thr, , drop = FALSE]
n_after <- nrow(expr_mat)

log_msg("MAD filter (top ", round((1 - MAD_QUANTILE) * 100), "%): ",
        n_before, " → ", n_after, " features")

# datExpr: samples × features (format WGCNA)
datExpr <- as.data.frame(t(expr_mat))
stopifnot(rownames(datExpr) == meta$Sample)

log_msg("datExpr ready: ", nrow(datExpr), " samples x ", ncol(datExpr), " features")

##############################
## PARTIE 7: Sauvegarder
##############################
cat("── EXPORT ──\n")
saveRDS(datExpr,      file.path(out_dir, "datExpr.rds"))
saveRDS(meta,         file.path(out_dir, "meta.rds"))
saveRDS(expr_mat,     file.path(out_dir, "expr_mat_log2_pareto_MAD.rds"))
saveRDS(feature_info, file.path(out_dir, "feature_info.rds"))

# Sauver aussi la matrice avant MAD (pour référence)
saveRDS(int_work, file.path(out_dir, "intensity_filtered_raw.rds"))

write.csv(meta, file.path(dir_inputs, "metadata_aligned.csv"), row.names = FALSE)

# Résumé du filtrage
filter_summary <- data.frame(
  Step = c("Raw features", "Detection filter", "After normalisation",
           "MAD filter (final)"),
  N_features = c(n_raw, n_after_detect, n_before, n_after)
)
write.csv(filter_summary, file.path(dir_inputs, "feature_filtering_summary.csv"), row.names = FALSE)
print(filter_summary)

# Sauver les paramètres
params <- data.frame(
  Parameter = c("GENOTYPE_CHOSEN","SEED","MIN_INTENSITY","MIN_DETECT_FRAC",
                "MAD_QUANTILE","LOG_TRANSFORM","PARETO_SCALE",
                "N_samples","N_features_final"),
  Value = c(GENOTYPE_CHOSEN, SEED, MIN_INTENSITY, MIN_DETECT_FRAC,
            MAD_QUANTILE, LOG_TRANSFORM, PARETO_SCALE,
            nrow(meta), ncol(datExpr))
)
write.csv(params, file.path(dir_inputs, "parameters.csv"), row.names = FALSE)

log_msg("Script 01 DONE — output: ", out_dir)
cat("\n✓ Script 01 terminé\n")
cat("  out_dir:", out_dir, "\n")