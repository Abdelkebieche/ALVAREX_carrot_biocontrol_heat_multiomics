############################################################
## 02_QC_PCA_dendro.R  (LC-MS WGCNA)
## PCA centroïdes (Temp, Trt, Time) + dendro + heatmap
# Repository configuration (portable; no HPC-specific paths)
if (!exists('CARROT_CONFIG_LOADED')) source('../../../config/config.R')


############################################################

##############################
## PARAMÈTRES
##############################
INSTALL_MISSING_PKGS <- FALSE
PCA_NTOP <- 500   # top features pour PCA (métabo = moins de features que RNAseq)

TEMP_COLORS <- c("T0" = "#313695", "T1" = "#74ADD1", "T2" = "#F46D43")
TRT_COLORS  <- c("Control" = "#4575B4", "Treated" = "#D73027")
TIME_COLORS <- c("P3" = "#56B4E9", "P4" = "#009E73")
TRT_SHAPES  <- c("Control" = 16, "Treated" = 17)
TIME_SHAPES <- c("P3" = 16, "P4" = 17)

base_out <- LCMS_WGCNA_ROOT

##############################
## SOURCE + LOAD
##############################
source("00_helpers.R")
load_pkgs(c("svglite","dplyr","ggplot2","ggrepel","patchwork"), install_missing = INSTALL_MISSING_PKGS)
load_pkgs(c("matrixStats","WGCNA"), install_missing = INSTALL_MISSING_PKGS, bioc = TRUE)

HAS_PHEATMAP <- requireNamespace("pheatmap", quietly = TRUE)
if (HAS_PHEATMAP) library(pheatmap)

# Trouver le dernier run
run_dir  <- find_latest_dir(base_out, "^WGCNA_LCMS_")
datExpr  <- readRDS(file.path(run_dir, "datExpr.rds"))
meta     <- readRDS(file.path(run_dir, "meta.rds"))
expr_mat <- readRDS(file.path(run_dir, "expr_mat_log2_pareto_MAD.rds"))

GENO <- as.character(meta$Genotype[1])
dir_qc <- file.path(run_dir, "01_QC")
dir.create(dir_qc, recursive = TRUE, showWarnings = FALSE)

cat("── Script 02: QC —", GENO, "──\n")
cat("  Samples:", nrow(meta), "| Features:", ncol(datExpr), "\n\n")

##############################
## PARTIE 1: PCA
##############################
cat("── PCA ──\n")

# Top variable features
rv <- matrixStats::rowVars(expr_mat)
n_top <- min(PCA_NTOP, length(rv))
idx <- order(rv, decreasing = TRUE)[seq_len(n_top)]
pca <- prcomp(t(expr_mat[idx, , drop = FALSE]), center = TRUE, scale. = FALSE)
pct <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

pca_df <- data.frame(pca$x[, 1:min(6, ncol(pca$x))]) |>
  dplyr::mutate(
    Sample        = meta$Sample,
    Temperature   = meta$Temperature,
    Treatment     = meta$Treatment,
    Sampling_time = meta$Sampling_time,
    Group         = meta$Group
  )

xlab <- paste0("PC1 (", pct[1], "%)")
ylab <- paste0("PC2 (", pct[2], "%)")

# ── PCA Vue 1: Température (centroïdes + trajectoires T0→T1→T2) ──
cs1 <- compute_centroids(pca_df, "PC1", "PC2",
                         group_vars = c("Temperature","Treatment","Sampling_time"),
                         trajectory_var = "Temperature",
                         trajectory_levels = c("T0","T1","T2"))

p_temp <- ggplot2::ggplot() +
  ggplot2::geom_segment(data = cs1$segments,
                        ggplot2::aes(x = cx, y = cy, xend = cxend, yend = cyend),
                        colour = "grey60", linewidth = 0.8,
                        arrow = ggplot2::arrow(length = ggplot2::unit(2, "mm"), type = "closed")) +
  ggplot2::geom_point(data = pca_df,
                      ggplot2::aes(PC1, PC2, colour = Temperature, shape = Treatment),
                      alpha = 0.35, size = 2) +
  ggplot2::geom_point(data = cs1$centroids,
                      ggplot2::aes(cx, cy, colour = Temperature, shape = Treatment),
                      size = 5, stroke = 1.2) +
  ggplot2::geom_text(data = cs1$centroids,
                     ggplot2::aes(cx, cy, label = paste(relabel_temp(Temperature), Treatment, relabel_time(Sampling_time), sep=".")),
                     size = 2.2, vjust = -1.2, colour = "grey30") +
  ggplot2::scale_colour_manual(values = TEMP_COLORS, labels = c("T0"="NH","T1"="HS-7","T2"="HS-2")) +
  ggplot2::scale_shape_manual(values = TRT_SHAPES) +
  ggplot2::labs(title = paste0("PCA — ", GENO, " (P3+P4) LC-MS"),
                subtitle = "Centroïdes par Température × Traitement × Time | trajectoire NH→HS-7→HS-2",
                x = xlab, y = ylab) +
  theme_pub()

save_plot(p_temp, file.path(dir_qc, "PCA_centroid_Temperature"), w = 9, h = 6.5)

# ── PCA Vue 2: Traitement (Control vs Treated) ──
cs2 <- compute_centroids(pca_df, "PC1", "PC2",
                         group_vars = c("Treatment","Temperature","Sampling_time"))

p_trt <- ggplot2::ggplot() +
  ggplot2::geom_point(data = pca_df,
                      ggplot2::aes(PC1, PC2, colour = Treatment, shape = Sampling_time),
                      alpha = 0.40, size = 2) +
  ggplot2::geom_point(data = cs2$centroids,
                      ggplot2::aes(cx, cy, colour = Treatment, shape = Sampling_time),
                      size = 5, stroke = 1.2) +
  ggplot2::facet_wrap(~ Temperature, labeller = ggplot2::labeller(Temperature = c("T0"="NH","T1"="HS-7","T2"="HS-2")), nrow = 1) +
  ggplot2::scale_colour_manual(values = TRT_COLORS) +
  ggplot2::scale_shape_manual(values = TIME_SHAPES) +
  ggplot2::labs(title = paste0("PCA — ", GENO, " | Control vs Treated par température (LC-MS)"),
                subtitle = "Facet = Température | Forme = Time point",
                x = xlab, y = ylab) +
  theme_pub()

save_plot(p_trt, file.path(dir_qc, "PCA_centroid_Treatment"), w = 12, h = 5)

# ── PCA Vue 3: Time point (P3 vs P4) ──
p_time <- ggplot2::ggplot() +
  ggplot2::geom_point(data = pca_df,
                      ggplot2::aes(PC1, PC2, colour = Sampling_time, shape = Treatment),
                      alpha = 0.40, size = 2) +
  ggplot2::facet_wrap(~ Temperature, labeller = ggplot2::labeller(Temperature = c("T0"="NH","T1"="HS-7","T2"="HS-2")), nrow = 1) +
  ggplot2::scale_colour_manual(values = TIME_COLORS, labels = c("P3"="48h post-trt","P4"="10d post-trt")) +
  ggplot2::scale_shape_manual(values = TRT_SHAPES) +
  ggplot2::labs(title = paste0("PCA — ", GENO, " | P3 vs P4 par température (LC-MS)"),
                x = xlab, y = ylab) +
  theme_pub()

save_plot(p_time, file.path(dir_qc, "PCA_centroid_TimePoint"), w = 12, h = 5)

# ── PCA Vue 4: Combinée (PC1-2 + PC3-4) ──
p_pc12 <- ggplot2::ggplot(pca_df, ggplot2::aes(PC1, PC2, colour = Temperature, shape = Treatment)) +
  ggplot2::geom_point(size = 2.5, alpha = 0.7) +
  ggplot2::scale_colour_manual(values = TEMP_COLORS, labels = c("T0"="NH","T1"="HS-7","T2"="HS-2")) +
  ggplot2::scale_shape_manual(values = TRT_SHAPES) +
  ggplot2::labs(x = xlab, y = ylab, title = "PC1 vs PC2") +
  theme_pub(base_size = 9)

p_pc34 <- ggplot2::ggplot(pca_df, ggplot2::aes(PC3, PC4, colour = Temperature, shape = Treatment)) +
  ggplot2::geom_point(size = 2.5, alpha = 0.7) +
  ggplot2::scale_colour_manual(values = TEMP_COLORS, labels = c("T0"="NH","T1"="HS-7","T2"="HS-2")) +
  ggplot2::scale_shape_manual(values = TRT_SHAPES) +
  ggplot2::labs(x = paste0("PC3 (", pct[3], "%)"), y = paste0("PC4 (", pct[4], "%)"),
                title = "PC3 vs PC4") +
  theme_pub(base_size = 9)

p_multi <- p_pc12 + p_pc34 + patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(title = paste0(GENO, " — Multi-PC view (P3+P4) LC-MS"))

save_plot(p_multi, file.path(dir_qc, "PCA_multiPC"), w = 13, h = 5.5)

##############################
## PARTIE 2: Dendrogramme échantillons
##############################
cat("── DENDROGRAMME ──\n")

sample_dist <- dist(datExpr)
sample_hc   <- hclust(sample_dist, method = "average")

sample_colors <- cbind(
  Temp = WGCNA::labels2colors(as.numeric(meta$Temperature)),
  Trt  = WGCNA::labels2colors(as.numeric(meta$Treatment)),
  Time = WGCNA::labels2colors(as.numeric(meta$Sampling_time))
)

save_base_plot(function(){
  WGCNA::plotDendroAndColors(
    dendro = sample_hc,
    colors = sample_colors,
    groupLabels = c("Temperature","Treatment","Time"),
    main = paste0("Sample clustering — ", GENO, " (P3+P4) LC-MS"),
    cex.dendroLabels = 0.6, cex.colorLabels = 0.8
  )
}, file.path(dir_qc, "Sample_dendrogram"), w = 14, h = 6)

##############################
## PARTIE 3: Heatmap top features
##############################
if (HAS_PHEATMAP) {
  cat("── HEATMAP TOP FEATURES ──\n")
  
  n_hm <- min(PCA_NTOP, nrow(expr_mat))
  hm_mat <- expr_mat[idx[1:n_hm], , drop = FALSE]
  
  # Z-score par feature
  hm_z <- t(scale(t(hm_mat)))
  hm_z[hm_z > 3] <- 3; hm_z[hm_z < -3] <- -3
  
  ann_col <- data.frame(
    Temperature = relabel_temp(meta$Temperature),
    Treatment   = meta$Treatment,
    Time        = relabel_time(meta$Sampling_time),
    row.names   = meta$Sample
  )
  
  ann_palette <- list(
    Temperature = setNames(TEMP_COLORS, relabel_temp(names(TEMP_COLORS))),
    Treatment   = TRT_COLORS,
    Time        = setNames(TIME_COLORS, relabel_time(names(TIME_COLORS)))
  )
  
  ph <- pheatmap::pheatmap(
    hm_z,
    annotation_col   = ann_col,
    annotation_colors = ann_palette,
    show_rownames    = FALSE,
    show_colnames    = FALSE,
    clustering_method = "ward.D2",
    color = colorRampPalette(c("#2166AC","#F7F7F7","#B2182B"))(100),
    border_color     = NA,
    fontsize         = 8,
    main = paste0("Top ", n_hm, " variable metabolites (z-score) — ", GENO, " LC-MS"),
    silent = TRUE
  )
  
  pdf(file.path(dir_qc, "Heatmap_top_features.pdf"), width = 10, height = 9)
  grid::grid.draw(ph$gtable)
  dev.off()
  
  png(file.path(dir_qc, "Heatmap_top_features.png"), width = 10, height = 9,
      units = "in", res = 300, bg = "white")
  grid::grid.draw(ph$gtable)
  dev.off()
  
  svglite::svglite(file.path(dir_qc, "Heatmap_top_features.svg"), width = 10, height = 9, bg = "white")
  grid::grid.draw(ph$gtable)
  dev.off()
}

cat("\n✓ Script 02 terminé\n")
cat("  Figures dans:", dir_qc, "\n")