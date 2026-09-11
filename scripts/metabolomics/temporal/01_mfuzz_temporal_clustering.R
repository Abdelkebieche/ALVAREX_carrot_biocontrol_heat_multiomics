############################################################
## MFuzz_LCMS_GLOBAL_P1P3P4.R
## Clustering temporel global MFuzz des métabolites LC-MS
## Présentation finale : P1 -> P3 -> P4
## Couleur = Température | Ligne = Traitement
## P1 partagé par température (baseline commune)
# Repository configuration (portable; no HPC-specific paths)
if (!exists('CARROT_CONFIG_LOADED')) source('../../../config/config.R')


############################################################

##############################
## PARAMÈTRES
##############################
GENOTYPE_CHOSEN <- "PRESTO"
SEED <- 123
INSTALL_MISSING_PKGS <- FALSE

# CHEMINS
lcms_file <- LCMS_FILE
meta_lcms_file <- LCMS_METADATA_FILE
base_out <- file.path(OUTPUT_ROOT, "metabolomics", "mfuzz")
helper_lcms_file <- file.path("..", "wgcna", "00_helpers.R")

# Filtrage LC-MS
LCMS_MIN_INTENSITY <- 1000
LCMS_MIN_DETECT    <- 0.70

# Filtrage dynamique global
DYN_SD_THRESHOLD   <- 0.30
DYN_AMP_THRESHOLD  <- 0.80
MIN_FEATURES_FOR_CLUSTERING <- 30

# MFuzz
C_RANGE            <- 4:12
C_CHOSEN           <- NULL
FUZZIFIER_M        <- NULL
MEMBERSHIP_CUTOFF  <- 0.30

# Sorties par cluster
TOP_HEATMAP_FEATURES_PER_CLUSTER <- 150

# Sauvegarde graphique
PNG_DPI             <- 200
PNG_MAX_HEIGHT_IN   <- 30
PNG_MAX_WIDTH_IN    <- 16

set.seed(SEED)
options(stringsAsFactors = FALSE)

##############################
## OUTILS GÉNÉRAUX
##############################
load_pkgs <- function(pkgs, install_missing = FALSE, bioc = FALSE){
  for (p in pkgs){
    if (!requireNamespace(p, quietly = TRUE)){
      if (!install_missing) stop("Package manquant: ", p)
      if (bioc){
        if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
        BiocManager::install(p, ask = FALSE, update = FALSE)
      } else {
        install.packages(p, dependencies = TRUE)
      }
    }
    suppressPackageStartupMessages(library(p, character.only = TRUE))
  }
}

format_sample_id <- function(x){
  sx <- trimws(as.character(x))
  suppressWarnings(nx <- as.integer(sx))
  out <- ifelse(!is.na(nx), sprintf("%04d", nx), sx)
  out
}

safe_mean <- function(x){
  if (all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

safe_se <- function(x){
  x <- x[!is.na(x)]
  if (length(x) <= 1) return(0)
  stats::sd(x) / sqrt(length(x))
}

close_all_devices <- function(){
  while (grDevices::dev.cur() > 1) {
    try(grDevices::dev.off(), silent = TRUE)
  }
}

##############################
## PACKAGES
##############################
load_pkgs(
  c("dplyr","tidyr","tibble","ggplot2","patchwork","pheatmap",
    "RColorBrewer","grid","matrixStats","svglite"),
  install_missing = INSTALL_MISSING_PKGS
)

load_pkgs(
  c("Mfuzz","Biobase"),
  install_missing = INSTALL_MISSING_PKGS,
  bioc = TRUE
)

HAS_RAGG <- requireNamespace("ragg", quietly = TRUE)

##############################
## HELPERS PROJET
##############################
source(helper_lcms_file)

##############################
## OUTILS GRAPHIQUES ROBUSTES
##############################
safe_open_png <- function(filename, width, height, res = 200, bg = "white"){
  close_all_devices()
  if (HAS_RAGG) {
    ragg::agg_png(
      filename = filename,
      width = width,
      height = height,
      units = "in",
      res = res,
      background = bg
    )
  } else {
    grDevices::png(
      filename = filename,
      width = width,
      height = height,
      units = "in",
      res = res,
      bg = bg
    )
  }
}

safe_save_ggplot <- function(plot_obj, file_base,
                             width = 7, height = 5,
                             png_max_height = PNG_MAX_HEIGHT_IN,
                             png_max_width = PNG_MAX_WIDTH_IN,
                             dpi = PNG_DPI,
                             bg = "white"){
  width_png  <- min(width, png_max_width)
  height_png <- min(height, png_max_height)
  
  try({
    ggplot2::ggsave(
      filename = paste0(file_base, ".pdf"),
      plot = plot_obj,
      width = width,
      height = height,
      units = "in",
      device = grDevices::pdf,
      bg = bg
    )
  }, silent = TRUE)
  
  try({
    if (HAS_RAGG) {
      ggplot2::ggsave(
        filename = paste0(file_base, ".png"),
        plot = plot_obj,
        width = width_png,
        height = height_png,
        units = "in",
        dpi = dpi,
        device = ragg::agg_png,
        bg = bg
      )
    } else {
      ggplot2::ggsave(
        filename = paste0(file_base, ".png"),
        plot = plot_obj,
        width = width_png,
        height = height_png,
        units = "in",
        dpi = dpi,
        device = "png",
        bg = bg
      )
    }
  }, silent = TRUE)
  
  try({
    ggplot2::ggsave(
      filename = paste0(file_base, ".svg"),
      plot = plot_obj,
      width = width,
      height = height,
      units = "in",
      device = svglite::svglite,
      bg = bg
    )
  }, silent = TRUE)
}

safe_save_grob <- function(grob_obj, file_base,
                           width = 6,
                           height_pdf = 8,
                           height_png = 8,
                           png_max_height = PNG_MAX_HEIGHT_IN,
                           png_max_width = PNG_MAX_WIDTH_IN,
                           dpi = PNG_DPI,
                           bg = "white"){
  width_png  <- min(width, png_max_width)
  height_png <- min(height_png, png_max_height)
  
  tryCatch({
    close_all_devices()
    grDevices::pdf(paste0(file_base, ".pdf"), width = width, height = height_pdf, onefile = FALSE)
    grid::grid.newpage()
    grid::grid.draw(grob_obj)
    grDevices::dev.off()
  }, error = function(e){
    message("Erreur PDF pour ", file_base, " : ", conditionMessage(e))
    close_all_devices()
  })
  
  tryCatch({
    safe_open_png(paste0(file_base, ".png"),
                  width = width_png, height = height_png,
                  res = dpi, bg = bg)
    grid::grid.newpage()
    grid::grid.draw(grob_obj)
    grDevices::dev.off()
  }, error = function(e){
    message("PNG non généré pour ", file_base, " : ", conditionMessage(e))
    close_all_devices()
  })
  
  tryCatch({
    close_all_devices()
    svglite::svglite(paste0(file_base, ".svg"), width = width, height = height_pdf, bg = bg)
    grid::grid.newpage()
    grid::grid.draw(grob_obj)
    grDevices::dev.off()
  }, error = function(e){
    message("SVG non généré pour ", file_base, " : ", conditionMessage(e))
    close_all_devices()
  })
}

##############################
## DOSSIER DE SORTIE
##############################
GENO_TAG  <- gsub("[^A-Za-z0-9]+", "", GENOTYPE_CHOSEN)
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_dir   <- file.path(base_out, paste0("MFuzz_GLOBAL_P1P3P4_", GENO_TAG, "_", timestamp))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_msg <- init_logger(file.path(out_dir, "pipeline.log"))
log_msg("START MFuzz GLOBAL P1P3P4 — ", GENOTYPE_CHOSEN)

cat("═══════════════════════════════════════════════════════\n")
cat("  MFuzz GLOBAL LC-MS | Présentation P1 -> P3 -> P4\n")
cat("  ", GENOTYPE_CHOSEN, "\n")
cat("═══════════════════════════════════════════════════════\n\n")

############################################################
## 1. IMPORT
############################################################
cat("── 1. IMPORT ──\n")

sep_meta <- ifelse(grepl(";", readLines(meta_lcms_file, n = 1)), ";", ",")
meta_raw <- read.csv(meta_lcms_file, sep = sep_meta, stringsAsFactors = FALSE)

if (!"Sample" %in% colnames(meta_raw)) stop("Colonne 'Sample' absente dans metadata.csv")
if (!"Genotype" %in% colnames(meta_raw)) stop("Colonne 'Genotype' absente dans metadata.csv")
if (!"Sampling_time" %in% colnames(meta_raw)) stop("Colonne 'Sampling_time' absente dans metadata.csv")
if (!"Temperature" %in% colnames(meta_raw)) stop("Colonne 'Temperature' absente dans metadata.csv")
if (!"Treatment" %in% colnames(meta_raw)) stop("Colonne 'Treatment' absente dans metadata.csv")

if (!"Replicat" %in% colnames(meta_raw)) {
  meta_raw$Replicat <- NA
}

meta <- meta_raw |>
  dplyr::mutate(
    Sample = format_sample_id(Sample),
    Genotype = trimws(Genotype),
    Sampling_time = trimws(Sampling_time),
    Temperature = normalize_temperature(Temperature),
    Treatment = normalize_treatment(Treatment)
  ) |>
  dplyr::filter(Genotype == GENOTYPE_CHOSEN) |>
  dplyr::filter(Sampling_time %in% c("P1","P3","P4")) |>
  dplyr::mutate(
    Sampling_time = factor(Sampling_time, levels = c("P1","P3","P4")),
    Temperature   = factor(Temperature, levels = c("T0","T1","T2")),
    Treatment     = factor(Treatment, levels = c("Control","Treated"))
  )

meta <- meta |>
  dplyr::filter(!is.na(Temperature), !is.na(Treatment), !is.na(Sampling_time))

cat("  Échantillons retenus :", nrow(meta), "\n")
print(table(meta$Temperature, meta$Treatment, meta$Sampling_time))
cat("\n")

lcms_data <- read_lcms_matrix(lcms_file)
lcms_raw  <- lcms_data$mat
feat_info <- lcms_data$feature_info

common <- intersect(meta$Sample, colnames(lcms_raw))
if (length(common) == 0) stop("Aucun Sample commun entre metadata et matrice LC-MS.")

meta <- meta |>
  dplyr::filter(Sample %in% common) |>
  dplyr::arrange(Temperature, Treatment, Sampling_time, Replicat, Sample)

lcms_mat <- lcms_raw[, meta$Sample, drop = FALSE]

cat("  Matrice brute :", nrow(lcms_mat), "features x", ncol(lcms_mat), "samples\n")

############################################################
## 2. FILTRAGE + IMPUTATION + LOG2
############################################################
cat("── 2. FILTRAGE / IMPUTATION / LOG2 ──\n")

lcms_mat[lcms_mat == 0] <- NA

detect_rate <- rowMeans(!is.na(lcms_mat) & lcms_mat >= LCMS_MIN_INTENSITY)
lcms_mat <- lcms_mat[detect_rate >= LCMS_MIN_DETECT, , drop = FALSE]

cat("  Features après filtrage détection :", nrow(lcms_mat), "\n")

for (i in seq_len(nrow(lcms_mat))){
  rv <- lcms_mat[i, ]
  na_i <- is.na(rv)
  if (any(na_i)) {
    non_na <- rv[!na_i]
    if (length(non_na) > 0) {
      lcms_mat[i, na_i] <- min(non_na, na.rm = TRUE) / 2
    }
  }
}

lcms_log <- log2(lcms_mat + 1)

cat("  Matrice log2 :", nrow(lcms_log), "x", ncol(lcms_log), "\n\n")

############################################################
## 3. MATRICE GLOBALE POUR MFUZZ
## P1 partagé par température
############################################################
cat("── 3. MATRICE GLOBALE MFUZZ ──\n")

temps_levels <- c("T0","T1","T2")
trt_levels   <- c("Control","Treated")

global_mat_list <- list()

for (temp in temps_levels) {
  
  # P1 partagé par température
  samp_p1 <- meta$Sample[
    as.character(meta$Temperature) == temp &
      as.character(meta$Sampling_time) == "P1"
  ]
  
  if (length(samp_p1) > 0) {
    col_name <- paste(temp, "P1", sep = "_")
    global_mat_list[[col_name]] <- rowMeans(lcms_log[, samp_p1, drop = FALSE], na.rm = TRUE)
  }
  
  # P3/P4 par traitement
  for (trt in trt_levels) {
    for (tp in c("P3","P4")) {
      
      samp <- meta$Sample[
        as.character(meta$Temperature) == temp &
          as.character(meta$Treatment) == trt &
          as.character(meta$Sampling_time) == tp
      ]
      
      if (length(samp) > 0) {
        col_name <- paste(temp, trt, tp, sep = "_")
        global_mat_list[[col_name]] <- rowMeans(lcms_log[, samp, drop = FALSE], na.rm = TRUE)
      }
    }
  }
}

group_mat <- do.call(cbind, global_mat_list)
group_mat <- as.matrix(group_mat)
rownames(group_mat) <- rownames(lcms_log)

write.csv(group_mat, file.path(out_dir, "01_group_means_matrix_sharedP1.csv"))

cat("  Matrice globale :", nrow(group_mat), "features x", ncol(group_mat), "groupes\n")
print(colnames(group_mat))
cat("\n")

############################################################
## 4. FEATURES DYNAMIQUES
############################################################
cat("── 4. FEATURES DYNAMIQUES ──\n")

feature_sd  <- matrixStats::rowSds(group_mat, na.rm = TRUE)
feature_amp <- matrixStats::rowMaxs(group_mat, na.rm = TRUE) -
  matrixStats::rowMins(group_mat, na.rm = TRUE)

dynamic_df <- data.frame(
  Feature = rownames(group_mat),
  SD_global = feature_sd,
  AMP_global = feature_amp,
  Dynamic = feature_sd >= DYN_SD_THRESHOLD & feature_amp >= DYN_AMP_THRESHOLD,
  stringsAsFactors = FALSE
)

write.csv(dynamic_df, file.path(out_dir, "02_dynamic_filter_stats.csv"), row.names = FALSE)

group_mat_dyn <- group_mat[dynamic_df$Dynamic, , drop = FALSE]

cat("  Features dynamiques :", nrow(group_mat_dyn), "/", nrow(group_mat), "\n\n")

if (nrow(group_mat_dyn) < MIN_FEATURES_FOR_CLUSTERING) {
  stop("Trop peu de features dynamiques pour lancer MFuzz.")
}

############################################################
## 5. MFUZZ GLOBAL
############################################################
cat("── 5. MFUZZ GLOBAL ──\n")

eset <- Biobase::ExpressionSet(assayData = group_mat_dyn)
eset_std <- Mfuzz::standardise(eset)
expr_std <- Biobase::exprs(eset_std)

# m
if (is.null(FUZZIFIER_M)) {
  m_est <- tryCatch(
    Mfuzz::mestimate(eset_std),
    error = function(e) {
      cat("  ⚠ mestimate a échoué, m=1.25 utilisé\n")
      1.25
    }
  )
  cat("  m estimé :", round(m_est, 3), "\n")
} else {
  m_est <- FUZZIFIER_M
  cat("  m fixé :", round(m_est, 3), "\n")
}

# choix de c
if (is.null(C_CHOSEN)) {
  max_c_allowed <- max(4, min(max(C_RANGE), floor(nrow(expr_std) / 8)))
  c_range_use <- C_RANGE[C_RANGE <= max_c_allowed]
  c_range_use <- c_range_use[c_range_use >= 4]
  
  if (length(c_range_use) == 0) c_range_use <- 4
  
  if (length(c_range_use) == 1) {
    best_c <- c_range_use[1]
    dist_df <- data.frame(c = best_c, min_dist = NA_real_)
  } else {
    min_dist <- sapply(c_range_use, function(cc){
      set.seed(SEED)
      cl_tmp <- Mfuzz::mfuzz(eset_std, c = cc, m = m_est)
      d <- as.matrix(stats::dist(cl_tmp$centers))
      diag(d) <- NA
      min(d, na.rm = TRUE)
    })
    
    dist_df <- data.frame(c = c_range_use, min_dist = min_dist)
    dist_threshold <- max(min_dist, na.rm = TRUE) * 0.4
    cand <- c_range_use[min_dist >= dist_threshold]
    best_c <- if (length(cand) > 0) max(cand) else c_range_use[which.max(min_dist)]
    
    if (is.na(best_c) || best_c < 4) {
      best_c <- c_range_use[which.max(min_dist)]
    }
  }
  
  cat("  Best c =", best_c, "\n")
  
  p_elbow <- ggplot2::ggplot(dist_df, ggplot2::aes(c, min_dist)) +
    ggplot2::geom_line(colour = "grey50") +
    ggplot2::geom_point(size = 3, colour = "#2471A3") +
    ggplot2::geom_vline(xintercept = best_c, linetype = "dashed", colour = "#C0392B") +
    ggplot2::labs(
      title = paste0("MFuzz GLOBAL — ", GENOTYPE_CHOSEN),
      x = "c",
      y = "Min centroid distance"
    ) +
    theme_pub()
  
  safe_save_ggplot(
    p_elbow,
    file.path(out_dir, "03_elbow_c_global"),
    width = 6,
    height = 4
  )
  
} else {
  best_c <- C_CHOSEN
  cat("  Best c fixé =", best_c, "\n")
}

set.seed(SEED)
cl <- tryCatch(
  Mfuzz::mfuzz(eset_std, c = best_c, m = m_est),
  error = function(e) {
    stop("MFuzz a échoué : ", conditionMessage(e))
  }
)

membership <- cl$membership
cluster_assignment <- cl$cluster
centers <- cl$centers
cl_sizes <- table(cluster_assignment)

cat("  Tailles des clusters :\n")
print(cl_sizes)
cat("\n")

############################################################
## 6. TABLE RÉSULTATS GLOBALE
############################################################
cat("── 6. TABLES DE RÉSULTATS ──\n")

result_df <- data.frame(
  Feature = names(cluster_assignment),
  Cluster = cluster_assignment,
  Max_membership = apply(membership, 1, max),
  stringsAsFactors = FALSE
)

for (cc in seq_len(best_c)) {
  result_df[[paste0("Mem_C", cc)]] <- membership[, cc]
}

for (grp in colnames(group_mat_dyn)) {
  result_df[[paste0("Mean_", grp)]] <- group_mat_dyn[result_df$Feature, grp]
}

if (!is.null(feat_info) && "FeatureID" %in% colnames(feat_info)) {
  info_cols <- intersect(colnames(feat_info), c("FeatureID","Tags","Name","m/z","RT [min]","ID_MT"))
  result_df <- result_df |>
    dplyr::left_join(feat_info[, info_cols, drop = FALSE], by = c("Feature" = "FeatureID"))
}

result_df <- result_df |>
  dplyr::arrange(Cluster, dplyr::desc(Max_membership))

write.csv(result_df, file.path(out_dir, "04_mfuzz_global_results.csv"), row.names = FALSE)
write.csv(centers, file.path(out_dir, "05_centroids_global.csv"))

############################################################
## 7. PLOT MFUZZ NATIF
############################################################
cat("── 7. PLOT MFUZZ NATIF ──\n")

group_labels_short <- colnames(expr_std)
group_labels_short <- gsub("Control", "Ctl", group_labels_short)
group_labels_short <- gsub("Treated", "Trt", group_labels_short)

mfuzz_pdf <- file.path(out_dir, "06_mfuzz_global_clusters.pdf")
mfuzz_png <- file.path(out_dir, "06_mfuzz_global_clusters.png")

plot_w <- min(max(10, best_c * 2.5), PNG_MAX_WIDTH_IN)
plot_h <- min(max(6, ceiling(best_c / 3) * 3.5), PNG_MAX_HEIGHT_IN)

tryCatch({
  close_all_devices()
  grDevices::pdf(
    mfuzz_pdf,
    width = max(10, best_c * 2.5),
    height = max(6, ceiling(best_c / 3) * 3.5)
  )
  Mfuzz::mfuzz.plot2(
    eset_std,
    cl = cl,
    mfrow = c(ceiling(best_c / 3), min(3, best_c)),
    time.labels = group_labels_short,
    min.mem = MEMBERSHIP_CUTOFF,
    colo = "fancy",
    x11 = FALSE
  )
  grDevices::dev.off()
}, error = function(e){
  message("Erreur PDF mfuzz.plot2 : ", conditionMessage(e))
  close_all_devices()
})

tryCatch({
  safe_open_png(mfuzz_png, width = plot_w, height = plot_h, res = PNG_DPI, bg = "white")
  Mfuzz::mfuzz.plot2(
    eset_std,
    cl = cl,
    mfrow = c(ceiling(best_c / 3), min(3, best_c)),
    time.labels = group_labels_short,
    min.mem = MEMBERSHIP_CUTOFF,
    colo = "fancy",
    x11 = FALSE
  )
  grDevices::dev.off()
}, error = function(e){
  message("PNG mfuzz.plot2 non généré : ", conditionMessage(e))
  close_all_devices()
})

############################################################
## 8. PROFILS STANDARDISÉS GLOBAUX PAR CLUSTER
############################################################
cat("── 8. PROFILS STANDARDISÉS PAR CLUSTER ──\n")

centers_df <- data.frame()
indiv_df   <- data.frame()

for (cc in seq_len(best_c)) {
  
  centers_df <- dplyr::bind_rows(
    centers_df,
    data.frame(
      Cluster = paste0("Cluster ", cc, " (n=", cl_sizes[as.character(cc)], ")"),
      ClusterNum = cc,
      Group = factor(colnames(centers), levels = colnames(centers)),
      Value = as.numeric(centers[cc, ]),
      stringsAsFactors = FALSE
    )
  )
  
  cc_feat <- names(cluster_assignment)[cluster_assignment == cc]
  cc_mem  <- membership[cc_feat, cc]
  cc_ok   <- cc_feat[cc_mem >= MEMBERSHIP_CUTOFF]
  
  if (length(cc_ok) > 0) {
    std_mat <- expr_std[cc_ok, , drop = FALSE]
    
    for (f in cc_ok) {
      indiv_df <- dplyr::bind_rows(
        indiv_df,
        data.frame(
          Cluster = paste0("Cluster ", cc, " (n=", cl_sizes[as.character(cc)], ")"),
          ClusterNum = cc,
          Group = factor(colnames(std_mat), levels = colnames(std_mat)),
          Value = as.numeric(std_mat[f, ]),
          Membership = membership[f, cc],
          Feature = f,
          stringsAsFactors = FALSE
        )
      )
    }
  }
}

p_prof_std <- ggplot2::ggplot() +
  ggplot2::geom_line(
    data = indiv_df,
    ggplot2::aes(Group, Value, group = Feature, alpha = Membership),
    colour = "grey65",
    linewidth = 0.30
  ) +
  ggplot2::geom_line(
    data = centers_df,
    ggplot2::aes(Group, Value, group = 1),
    colour = "#C0392B",
    linewidth = 1.4
  ) +
  ggplot2::geom_point(
    data = centers_df,
    ggplot2::aes(Group, Value),
    colour = "#C0392B",
    size = 2.4
  ) +
  ggplot2::scale_alpha_continuous(range = c(0.05, 0.5), guide = "none") +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey45") +
  ggplot2::facet_wrap(~Cluster, scales = "free_y", ncol = 3) +
  ggplot2::labs(
    title = paste0("MFuzz GLOBAL — ", GENOTYPE_CHOSEN),
    subtitle = paste0("Profils standardisés | c=", best_c, " | m=", round(m_est, 2)),
    x = "Groupe global",
    y = "Expression standardisée"
  ) +
  theme_pub(base_size = 10) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

safe_save_ggplot(
  p_prof_std,
  file.path(out_dir, "07_mfuzz_global_profiles_standardised"),
  width = max(12, best_c * 2.5),
  height = max(6, ceiling(best_c / 3) * 3.8)
)

############################################################
## 9. HEATMAP GLOBALE DES MOYENNES DE GROUPES
############################################################
cat("── 9. HEATMAP GLOBALE ──\n")

ord <- order(result_df$Cluster, -result_df$Max_membership)
hm_mat <- expr_std[result_df$Feature[ord], , drop = FALSE]
hm_mat[hm_mat > 3]  <- 3
hm_mat[hm_mat < -3] <- -3

ann_row <- data.frame(
  Cluster = factor(paste0("C", result_df$Cluster[ord])),
  row.names = rownames(hm_mat)
)

if (best_c <= 12) {
  cl_pal <- RColorBrewer::brewer.pal(max(3, best_c), "Set3")[seq_len(best_c)]
} else {
  cl_pal <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(best_c)
}

ann_colors <- list(
  Cluster = setNames(cl_pal[seq_len(best_c)], paste0("C", seq_len(best_c)))
)

ph_global <- pheatmap::pheatmap(
  hm_mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_row = ann_row,
  annotation_colors = ann_colors,
  show_rownames = FALSE,
  color = grDevices::colorRampPalette(c("#2166AC","#F7F7F7","#B2182B"))(100),
  border_color = NA,
  fontsize = 9,
  main = paste0("MFuzz GLOBAL — ", GENOTYPE_CHOSEN),
  silent = TRUE
)

h_fig_pdf <- max(7, nrow(hm_mat) * 0.010 + 3)
h_fig_png <- min(h_fig_pdf, PNG_MAX_HEIGHT_IN)

safe_save_grob(
  grob_obj = ph_global$gtable,
  file_base = file.path(out_dir, "08_global_heatmap_group_means"),
  width = 9,
  height_pdf = h_fig_pdf,
  height_png = h_fig_png,
  dpi = PNG_DPI
)

############################################################
## 10. FIGURE FINALE : P1 -> P3 -> P4
## lignes par condition dans chaque cluster
############################################################
cat("── 10. FIGURE FINALE P1 -> P3 -> P4 ──\n")

plot_list <- list()

for (cc in seq_len(best_c)) {
  
  cc_feat <- result_df |>
    dplyr::filter(Cluster == cc, Max_membership >= MEMBERSHIP_CUTOFF) |>
    dplyr::pull(Feature)
  
  # fallback si cutoff trop strict
  if (length(cc_feat) < 2) {
    cc_feat <- result_df |>
      dplyr::filter(Cluster == cc) |>
      dplyr::pull(Feature)
  }
  
  if (length(cc_feat) < 2) next
  
  # score moyen du cluster par sample
  cluster_sample_score <- data.frame(
    Sample = colnames(lcms_log),
    ClusterScore = colMeans(lcms_log[cc_feat, colnames(lcms_log), drop = FALSE], na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  
  plot_df_raw <- meta |>
    dplyr::left_join(cluster_sample_score, by = "Sample")
  
  # P1 partagé par température, puis dupliqué visuellement pour Control/Treated
  p1_df <- plot_df_raw |>
    dplyr::filter(Sampling_time == "P1") |>
    dplyr::group_by(Temperature, Sampling_time) |>
    dplyr::summarise(
      MeanExpr = safe_mean(ClusterScore),
      SE = safe_se(ClusterScore),
      Nrep = sum(!is.na(ClusterScore)),
      .groups = "drop"
    ) |>
    tidyr::crossing(
      Treatment = factor(c("Control","Treated"), levels = c("Control","Treated"))
    )
  
  # P3/P4 par température x traitement
  p34_df <- plot_df_raw |>
    dplyr::filter(Sampling_time %in% c("P3","P4")) |>
    dplyr::group_by(Temperature, Treatment, Sampling_time) |>
    dplyr::summarise(
      MeanExpr = safe_mean(ClusterScore),
      SE = safe_se(ClusterScore),
      Nrep = sum(!is.na(ClusterScore)),
      .groups = "drop"
    )
  
  cc_plot_df <- dplyr::bind_rows(p1_df, p34_df) |>
    dplyr::mutate(
      Cluster = paste0("Cluster ", cc, " (n=", length(cc_feat), ")"),
      ClusterNum = cc,
      TimePoint = factor(as.character(Sampling_time), levels = c("P1","P3","P4"))
    ) |>
    dplyr::arrange(Temperature, Treatment, TimePoint)
  
  plot_list[[cc]] <- cc_plot_df
}

plot_clusters_df <- dplyr::bind_rows(plot_list)

write.csv(
  plot_clusters_df,
  file.path(out_dir, "09_cluster_profiles_for_plot.csv"),
  row.names = FALSE
)

# couleurs température
temp_colors <- c(
  T0 = "#3B3EAC",
  T1 = "#77B5D9",
  T2 = "#F46D43"
)

p_cluster_summary <- ggplot2::ggplot(
  plot_clusters_df,
  ggplot2::aes(
    x = TimePoint,
    y = MeanExpr,
    colour = Temperature,
    linetype = Treatment,
    group = interaction(Temperature, Treatment)
  )
) +
  ggplot2::geom_line(linewidth = 0.9) +
  ggplot2::geom_point(size = 2.5) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = MeanExpr - SE, ymax = MeanExpr + SE),
    width = 0.08,
    linewidth = 0.5
  ) +
  ggplot2::facet_wrap(~Cluster, scales = "free_y", ncol = 3) +
  ggplot2::scale_colour_manual(values = temp_colors) +
  ggplot2::scale_linetype_manual(values = c("Control" = "dashed", "Treated" = "solid")) +
  ggplot2::labs(
    title = "Mean expression profiles by cluster",
    subtitle = paste0(
      GENOTYPE_CHOSEN,
      " | P1->P3->P4 (shared baseline) | ",
      nrow(group_mat_dyn),
      " metabolites | Solid=Treated, Dashed=Control"
    ),
    x = "Time point",
    y = "Mean log2 expression",
    colour = "Temperature",
    linetype = "Treatment"
  ) +
  theme_pub(base_size = 12) +
  ggplot2::theme(
    strip.text = ggplot2::element_text(face = "bold", size = 12),
    axis.text.x = ggplot2::element_text(face = "bold"),
    legend.position = "top"
  )

safe_save_ggplot(
  p_cluster_summary,
  file.path(out_dir, "10_cluster_profiles_P1_P3_P4_by_condition"),
  width = 14,
  height = 10
)

############################################################
## 11. STRUCTURE INTERNE DES CLUSTERS : SAMPLES INDIVIDUELS
############################################################
cat("── 11. STRUCTURE INTERNE DES CLUSTERS ──\n")

dir_cl <- file.path(out_dir, "by_cluster")
dir.create(dir_cl, recursive = TRUE, showWarnings = FALSE)

sample_order <- meta |>
  dplyr::arrange(Temperature, Treatment, Sampling_time, Replicat, Sample) |>
  dplyr::pull(Sample)

indiv_mat <- lcms_log[, sample_order, drop = FALSE]

# standardisation par feature sur échantillons individuels
indiv_mat_std <- t(scale(t(indiv_mat)))
indiv_mat_std[is.na(indiv_mat_std)] <- 0

ann_col <- meta |>
  dplyr::select(Sample, Temperature, Treatment, Sampling_time, Replicat) |>
  dplyr::distinct() |>
  dplyr::slice(match(sample_order, Sample)) |>
  as.data.frame()

rownames(ann_col) <- ann_col$Sample
ann_col$Sample <- NULL

ann_col_colors <- list(
  Temperature = c(T0 = "#66C2A5", T1 = "#FC8D62", T2 = "#8DA0CB"),
  Treatment = c(Control = "#BDBDBD", Treated = "#E41A1C"),
  Sampling_time = c(P1 = "#4DAF4A", P3 = "#377EB8", P4 = "#984EA3")
)

for (cc in seq_len(best_c)) {
  
  cc_res <- result_df |>
    dplyr::filter(Cluster == cc) |>
    dplyr::arrange(dplyr::desc(Max_membership))
  
  write.csv(
    cc_res,
    file.path(dir_cl, paste0("cluster_", cc, "_features.csv")),
    row.names = FALSE
  )
  
  cc_feat <- cc_res |>
    dplyr::filter(Max_membership >= MEMBERSHIP_CUTOFF) |>
    dplyr::pull(Feature)
  
  if (length(cc_feat) == 0) {
    cc_feat <- cc_res$Feature
  }
  
  # table longue valeurs individuelles
  cc_long <- as.data.frame(indiv_mat[cc_feat, sample_order, drop = FALSE]) |>
    tibble::rownames_to_column("Feature") |>
    tidyr::pivot_longer(
      cols = -Feature,
      names_to = "Sample",
      values_to = "log2_value"
    ) |>
    dplyr::left_join(
      meta |>
        dplyr::select(Sample, Temperature, Treatment, Sampling_time, Replicat),
      by = "Sample"
    ) |>
    dplyr::left_join(
      result_df |>
        dplyr::select(Feature, Cluster, Max_membership),
      by = "Feature"
    )
  
  write.csv(
    cc_long,
    file.path(dir_cl, paste0("cluster_", cc, "_individual_values_long.csv")),
    row.names = FALSE
  )
  
  # heatmap cluster sur les samples individuels
  cc_feat_heat <- cc_res |>
    dplyr::filter(Max_membership >= MEMBERSHIP_CUTOFF) |>
    dplyr::slice_head(n = TOP_HEATMAP_FEATURES_PER_CLUSTER) |>
    dplyr::pull(Feature)
  
  if (length(cc_feat_heat) == 0) {
    cc_feat_heat <- cc_res |>
      dplyr::slice_head(n = min(TOP_HEATMAP_FEATURES_PER_CLUSTER, nrow(cc_res))) |>
      dplyr::pull(Feature)
  }
  
  hm_indiv <- indiv_mat_std[cc_feat_heat, sample_order, drop = FALSE]
  hm_indiv[hm_indiv > 3]  <- 3
  hm_indiv[hm_indiv < -3] <- -3
  
  ph_cc <- pheatmap::pheatmap(
    hm_indiv,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    annotation_col = ann_col,
    annotation_colors = ann_col_colors,
    show_rownames = FALSE,
    color = grDevices::colorRampPalette(c("#2166AC","#F7F7F7","#B2182B"))(100),
    border_color = NA,
    fontsize = 8,
    main = paste0("Cluster ", cc, " — individual samples"),
    silent = TRUE
  )
  
  h_pdf <- max(6, nrow(hm_indiv) * 0.012 + 3)
  h_png <- min(h_pdf, PNG_MAX_HEIGHT_IN)
  
  safe_save_grob(
    grob_obj = ph_cc$gtable,
    file_base = file.path(dir_cl, paste0("cluster_", cc, "_individual_heatmap")),
    width = 10,
    height_pdf = h_pdf,
    height_png = h_png,
    dpi = PNG_DPI
  )
}

############################################################
## 12. TABLES RÉSUMÉES PAR CLUSTER
############################################################
cat("── 12. TABLES RÉSUMÉES ──\n")

cluster_summary <- result_df |>
  dplyr::group_by(Cluster) |>
  dplyr::summarise(
    N_features = n(),
    N_membership_ge_cutoff = sum(Max_membership >= MEMBERSHIP_CUTOFF, na.rm = TRUE),
    Mean_membership = mean(Max_membership, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::arrange(Cluster)

write.csv(cluster_summary, file.path(out_dir, "11_cluster_summary.csv"), row.names = FALSE)

if ("Name" %in% colnames(result_df)) {
  named_summary <- result_df |>
    dplyr::filter(!is.na(Name), Name != "") |>
    dplyr::group_by(Cluster) |>
    dplyr::summarise(
      N_annotated = n(),
      Top_names = paste(head(Name[order(-Max_membership)], 10), collapse = " | "),
      .groups = "drop"
    )
  write.csv(named_summary, file.path(out_dir, "12_cluster_annotation_summary.csv"), row.names = FALSE)
}

############################################################
## 13. PARAMÈTRES + OBJET FINAL
############################################################
cat("── 13. SAUVEGARDE FINALE ──\n")

params <- data.frame(
  Parameter = c(
    "GENOTYPE",
    "N_samples",
    "N_features_total_after_filter",
    "N_features_dynamic",
    "LCMS_MIN_INTENSITY",
    "LCMS_MIN_DETECT",
    "DYN_SD_THRESHOLD",
    "DYN_AMP_THRESHOLD",
    "C_RANGE",
    "BEST_C",
    "FUZZIFIER_M",
    "MEMBERSHIP_CUTOFF",
    "N_groups_for_MFuzz"
  ),
  Value = c(
    GENOTYPE_CHOSEN,
    ncol(lcms_log),
    nrow(lcms_log),
    nrow(group_mat_dyn),
    LCMS_MIN_INTENSITY,
    LCMS_MIN_DETECT,
    DYN_SD_THRESHOLD,
    DYN_AMP_THRESHOLD,
    paste(range(C_RANGE), collapse = "-"),
    best_c,
    round(m_est, 4),
    MEMBERSHIP_CUTOFF,
    ncol(group_mat_dyn)
  )
)

write.csv(params, file.path(out_dir, "13_parameters.csv"), row.names = FALSE)

saveRDS(
  list(
    meta = meta,
    lcms_log = lcms_log,
    group_mat = group_mat,
    group_mat_dyn = group_mat_dyn,
    expr_std = expr_std,
    mfuzz = cl,
    membership = membership,
    result_df = result_df,
    plot_clusters_df = plot_clusters_df
  ),
  file.path(out_dir, "14_mfuzz_global_object.rds")
)

log_msg("MFuzz GLOBAL P1P3P4 DONE")
cat("\n✓ MFuzz GLOBAL terminé\n")
cat("  Résultats :", out_dir, "\n")