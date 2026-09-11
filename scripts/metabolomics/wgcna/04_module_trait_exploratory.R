############################################################
## 04_module_trait.R — MODULE-TRAIT ANALYSIS (LC-MS)
## ─────────────────────────────────────────────────────────
## 1. Heatmap ME × ALL traits
## 2. ANOVA par module: Temp + Trt + Time + Temp:Trt
## 3. Heatmap ANOVA (p-values par facteur)
## 4. Contrastes spécifiques (SDP_in_T0, SDP_in_T1, etc.)
## 5. Boxplots ME par condition
## 6. Classification automatique des modules
# Repository configuration (portable; no HPC-specific paths)
if (!exists('CARROT_CONFIG_LOADED')) source('../../../config/config.R')


############################################################
## DIFFÉRENCES vs RNAseq:
##  - P3/P4 au lieu de P2/P3
##  - Labels: P3 = 48h post-trt, P4 = 10d post-trt
##  - Tout le reste est identique (analyse sur MEs)
############################################################

##############################
## PARAMÈTRES
##############################
INSTALL_MISSING_PKGS <- FALSE
FDR_CUTOFF <- 0.05

TEMP_COLORS <- c("T0" = "#313695", "T1" = "#74ADD1", "T2" = "#F46D43")
TRT_COLORS  <- c("Control" = "#4575B4", "Treated" = "#D73027")
TIME_COLORS <- c("P3" = "#56B4E9", "P4" = "#009E73")

base_out <- LCMS_WGCNA_ROOT

##############################
## SOURCE + LOAD
##############################
source("00_helpers.R")
load_pkgs(c("svglite","dplyr","tidyr","ggplot2","patchwork"), install_missing = INSTALL_MISSING_PKGS)
load_pkgs(c("WGCNA"), install_missing = INSTALL_MISSING_PKGS, bioc = TRUE)

run_dir <- find_latest_dir(base_out, "^WGCNA_LCMS_")
net     <- readRDS(file.path(run_dir, "02_network", "net.rds"))
MEs     <- readRDS(file.path(run_dir, "02_network", "MEs.rds"))
meta    <- readRDS(file.path(run_dir, "02_network", "meta_clean.rds"))
datExpr <- readRDS(file.path(run_dir, "02_network", "datExpr_clean.rds"))
me_map  <- read.csv(file.path(run_dir, "02_network", "ME_color_map.csv"))

GENO <- as.character(meta$Genotype[1])
stopifnot(rownames(MEs) == meta$Sample)

dir_mt <- file.path(run_dir, "03_module_trait")
dir_mt_box <- file.path(dir_mt, "boxplots")
dir_mt_bar <- file.path(dir_mt, "barplots_ME")
for (d in c(dir_mt, dir_mt_box, dir_mt_bar)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

log_msg <- init_logger(file.path(run_dir, "pipeline.log"))
log_msg("START Script 04 — Module-Trait (LC-MS)")
cat("── Script 04: Module-Trait —", GENO, "──\n\n")

# Exclude grey
MEs_clean <- MEs[, !grepl("^ME0$", colnames(MEs)), drop = FALSE]

# ME to color mapping
me_to_color <- setNames(me_map$ModuleColor, me_map$ME)


##############################################################
## PARTIE 1: MATRICE DE TRAITS BINAIRES
##############################################################
cat("── PARTIE 1: Matrice de traits ──\n")

# 1A) Indicatrices simples
traits <- data.frame(
  T0 = as.integer(meta$Temperature == "T0"),
  T1 = as.integer(meta$Temperature == "T1"),
  T2 = as.integer(meta$Temperature == "T2"),
  Control = as.integer(meta$Treatment == "Control"),
  Treated = as.integer(meta$Treatment == "Treated"),
  P3 = as.integer(meta$Sampling_time == "P3"),
  P4 = as.integer(meta$Sampling_time == "P4"),
  row.names = meta$Sample
)

# 1B) Contrastes codés (-1/+1)
x_SDP <- ifelse(meta$Treatment == "Treated", 1, ifelse(meta$Treatment == "Control", -1, NA))
x_T1  <- ifelse(meta$Temperature == "T1", 1, ifelse(meta$Temperature == "T0", -1, NA))
x_T2  <- ifelse(meta$Temperature == "T2", 1, ifelse(meta$Temperature == "T0", -1, NA))
x_P4  <- ifelse(meta$Sampling_time == "P4", 1, ifelse(meta$Sampling_time == "P3", -1, NA))

contrasts <- data.frame(
  # SDP effect global
  SDP_effect = x_SDP,
  # Treated dans chaque T°
  SDP_in_T0 = ifelse(meta$Temperature == "T0", x_SDP, NA),
  SDP_in_T1 = ifelse(meta$Temperature == "T1", x_SDP, NA),
  SDP_in_T2 = ifelse(meta$Temperature == "T2", x_SDP, NA),
  # Temp contrasts
  T1_vs_T0 = ifelse(meta$Temperature %in% c("T0","T1"), x_T1, NA),
  T2_vs_T0 = ifelse(meta$Temperature %in% c("T0","T2"), x_T2, NA),
  # Time contrast (P4 vs P3)
  P4_vs_P3 = x_P4,
  # Temp dans chaque traitement
  T1vsT0_in_EAU = ifelse(meta$Treatment == "Control" & meta$Temperature %in% c("T0","T1"), x_T1, NA),
  T1vsT0_in_SDP = ifelse(meta$Treatment == "Treated" & meta$Temperature %in% c("T0","T1"), x_T1, NA),
  T2vsT0_in_EAU = ifelse(meta$Treatment == "Control" & meta$Temperature %in% c("T0","T2"), x_T2, NA),
  T2vsT0_in_SDP = ifelse(meta$Treatment == "Treated" & meta$Temperature %in% c("T0","T2"), x_T2, NA),
  row.names = meta$Sample
)

traits_all <- cbind(traits, contrasts)
write.csv(traits_all, file.path(dir_mt, "Traits_matrix.csv"))


##############################################################
## PARTIE 2: CORRÉLATIONS bicor (ME × chaque trait)
##############################################################
cat("── PARTIE 2: Corrélations bicor ──\n")

bicor_safe <- function(MEs_mat, trait_vec) {
  keep <- !is.na(trait_vec)
  if (sum(keep) < 5) return(list(r = rep(NA, ncol(MEs_mat)), p = rep(NA, ncol(MEs_mat))))
  
  r_mat <- WGCNA::bicor(MEs_mat[keep, , drop = FALSE],
                        trait_vec[keep], use = "pairwise.complete.obs")
  r_vec <- as.numeric(r_mat[, 1])
  
  n <- sum(keep)
  p_vec <- sapply(r_vec, function(r) {
    if (is.na(r) || abs(r) >= 1 || n < 3) return(NA_real_)
    t_stat <- r * sqrt((n - 2) / (1 - r^2))
    2 * pt(-abs(t_stat), df = n - 2)
  })
  list(r = r_vec, p = as.numeric(p_vec))
}

cor_list <- list()
for (trait_name in colnames(traits_all)) {
  res <- bicor_safe(MEs_clean, as.numeric(traits_all[[trait_name]]))
  cor_list[[trait_name]] <- data.frame(
    ME = colnames(MEs_clean),
    bicor = res$r,
    pvalue = res$p
  )
}

cor_mat <- do.call(cbind, lapply(cor_list, function(x) x$bicor))
rownames(cor_mat) <- colnames(MEs_clean)
colnames(cor_mat) <- names(cor_list)

p_mat <- do.call(cbind, lapply(cor_list, function(x) x$pvalue))
rownames(p_mat) <- colnames(MEs_clean)
colnames(p_mat) <- names(cor_list)

fdr_mat <- matrix(p.adjust(as.vector(p_mat), "BH"),
                  nrow = nrow(p_mat), ncol = ncol(p_mat),
                  dimnames = dimnames(p_mat))

me_colors <- me_to_color[rownames(cor_mat)]
me_colors_u <- make.unique(as.character(me_colors))
rownames(cor_mat) <- me_colors_u
rownames(p_mat) <- me_colors_u
rownames(fdr_mat) <- me_colors_u

write.csv(cor_mat, file.path(dir_mt, "ME_bicor_AllTraits.csv"))
write.csv(p_mat, file.path(dir_mt, "ME_pvalue_AllTraits.csv"))
write.csv(fdr_mat, file.path(dir_mt, "ME_FDR_AllTraits.csv"))

cat("  Matrices:", nrow(cor_mat), "modules x", ncol(cor_mat), "traits\n")


##############################################################
## PARTIE 3: HEATMAPS
##############################################################
cat("── PARTIE 3: Heatmaps ──\n")

# 3A) Conditions simples
simple_cols <- c("T0","T1","T2","Control","Treated","P3","P4")
simple_cols_display <- c("NH","HS-7","HS-2","Control","Treated","P3 (48h)","P4 (10d)")
stars_simple <- matrix(starify(fdr_mat[, simple_cols]),
                       nrow = nrow(fdr_mat), ncol = length(simple_cols))
text_simple <- paste0(sprintf("%.2f", cor_mat[, simple_cols]), stars_simple)
dim(text_simple) <- dim(cor_mat[, simple_cols])

save_base_plot(function(){
  par(mar = c(10, 10, 4, 2))
  WGCNA::labeledHeatmap(
    Matrix = cor_mat[, simple_cols],
    xLabels = simple_cols_display,
    yLabels = rownames(cor_mat),
    ySymbols = rownames(cor_mat),
    colorLabels = FALSE,
    colors = WGCNA::blueWhiteRed(50),
    textMatrix = text_simple,
    setStdMargins = FALSE,
    cex.text = 0.55, zlim = c(-1, 1),
    xLabelsAngle = 45,
    main = paste0("ME vs conditions — ", GENO, " (P3+P4) LC-MS\nbicor | *FDR<0.05 **<0.01 ***<0.001")
  )
}, file.path(dir_mt, "Heatmap_ME_vs_Conditions"), w = 10,
h = max(8, nrow(cor_mat) * 0.30 + 3))

# 3B) Contrastes
contrast_cols <- c("SDP_in_T0","SDP_in_T1","SDP_in_T2","T1_vs_T0","T2_vs_T0","P4_vs_P3")
contrast_cols_display <- c("Treated vs Control (NH)","Treated vs Control (HS-7)","Treated vs Control (HS-2)","HS-7 vs NH","HS-2 vs NH","P4 vs P3")
contrast_cols <- intersect(contrast_cols, colnames(cor_mat))

if (length(contrast_cols) > 0) {
  stars_c <- matrix(starify(fdr_mat[, contrast_cols]),
                    nrow = nrow(fdr_mat), ncol = length(contrast_cols))
  text_c <- paste0(sprintf("%.2f", cor_mat[, contrast_cols]), stars_c)
  dim(text_c) <- dim(cor_mat[, contrast_cols])
  
  save_base_plot(function(){
    par(mar = c(12, 10, 4, 2))
    WGCNA::labeledHeatmap(
      Matrix = cor_mat[, contrast_cols],
      xLabels = contrast_cols_display[seq_along(contrast_cols)],
      yLabels = rownames(cor_mat),
      ySymbols = rownames(cor_mat),
      colorLabels = FALSE,
      colors = WGCNA::blueWhiteRed(50),
      textMatrix = text_c,
      setStdMargins = FALSE,
      cex.text = 0.55, zlim = c(-1, 1),
      xLabelsAngle = 45,
      main = paste0("ME vs effects — ", GENO, " (P3+P4) LC-MS\nbicor | *FDR<0.05")
    )
  }, file.path(dir_mt, "Heatmap_ME_vs_Effects"), w = 10,
  h = max(8, nrow(cor_mat) * 0.30 + 3))
}

# 3C) Détaillé
detailed_cols <- c("SDP_in_T0","SDP_in_T1","SDP_in_T2",
                   "T1vsT0_in_EAU","T1vsT0_in_SDP",
                   "T2vsT0_in_EAU","T2vsT0_in_SDP",
                   "P4_vs_P3")
detailed_cols <- intersect(detailed_cols, colnames(cor_mat))

detailed_display_map <- c(
  "SDP_in_T0" = "Treated vs Ctrl (NH)",
  "SDP_in_T1" = "Treated vs Ctrl (HS-7)",
  "SDP_in_T2" = "Treated vs Ctrl (HS-2)",
  "T1vsT0_in_EAU" = "HS-7 vs NH (Control)",
  "T1vsT0_in_SDP" = "HS-7 vs NH (Treated)",
  "T2vsT0_in_EAU" = "HS-2 vs NH (Control)",
  "T2vsT0_in_SDP" = "HS-2 vs NH (Treated)",
  "P4_vs_P3" = "P4 vs P3 (10d vs 48h)"
)
detailed_cols_display <- detailed_display_map[detailed_cols]

if (length(detailed_cols) > 2) {
  stars_d <- matrix(starify(fdr_mat[, detailed_cols]),
                    nrow = nrow(fdr_mat), ncol = length(detailed_cols))
  text_d <- paste0(sprintf("%.2f", cor_mat[, detailed_cols]), stars_d)
  dim(text_d) <- dim(cor_mat[, detailed_cols])
  
  save_base_plot(function(){
    par(mar = c(14, 10, 4, 2))
    WGCNA::labeledHeatmap(
      Matrix = cor_mat[, detailed_cols],
      xLabels = detailed_cols_display,
      yLabels = rownames(cor_mat),
      ySymbols = rownames(cor_mat),
      colorLabels = FALSE,
      colors = WGCNA::blueWhiteRed(50),
      textMatrix = text_d,
      setStdMargins = FALSE,
      cex.text = 0.50, zlim = c(-1, 1),
      xLabelsAngle = 60,
      main = paste0("ME vs detailed effects — ", GENO, " LC-MS\nTreated vs Control per T° | T° per Treatment | P4 vs P3")
    )
  }, file.path(dir_mt, "Heatmap_ME_vs_DetailedEffects"),
  w = max(10, length(detailed_cols) * 0.9 + 5),
  h = max(8, nrow(cor_mat) * 0.30 + 3))
}


##############################################################
## PARTIE 4: ANOVA
##############################################################
cat("── PARTIE 4: ANOVA ──\n")

anova_results <- data.frame()

for (me in colnames(MEs_clean)) {
  colr <- as.character(me_to_color[me])
  
  df <- data.frame(
    ME          = as.numeric(MEs_clean[, me]),
    Temperature = meta$Temperature,
    Treatment   = meta$Treatment,
    Time        = meta$Sampling_time
  )
  
  fit <- tryCatch(
    lm(ME ~ Temperature * Treatment + Time, data = df),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    anova_results <- dplyr::bind_rows(anova_results, data.frame(
      ME = me, ModuleColor = colr,
      p_Temperature = NA, p_Treatment = NA, p_Time = NA, p_TempxTrt = NA,
      R2 = NA, stringsAsFactors = FALSE
    ))
    next
  }
  
  a <- tryCatch({
    if (requireNamespace("car", quietly = TRUE)) {
      car::Anova(fit, type = "II")
    } else {
      anova(fit)
    }
  }, error = function(e) anova(fit))
  
  pT <- tryCatch(a["Temperature", "Pr(>F)"], error = function(e) NA)
  pTr <- tryCatch(a["Treatment", "Pr(>F)"], error = function(e) NA)
  pTi <- tryCatch(a["Time", "Pr(>F)"], error = function(e) NA)
  pTxTr <- tryCatch(a["Temperature:Treatment", "Pr(>F)"], error = function(e) NA)
  r2 <- summary(fit)$r.squared
  
  anova_results <- dplyr::bind_rows(anova_results, data.frame(
    ME = me, ModuleColor = colr,
    p_Temperature = as.numeric(pT),
    p_Treatment = as.numeric(pTr),
    p_Time = as.numeric(pTi),
    p_TempxTrt = as.numeric(pTxTr),
    R2 = as.numeric(r2),
    stringsAsFactors = FALSE
  ))
}

anova_results$FDR_Temperature <- p.adjust(anova_results$p_Temperature, "BH")
anova_results$FDR_Treatment   <- p.adjust(anova_results$p_Treatment, "BH")
anova_results$FDR_Time        <- p.adjust(anova_results$p_Time, "BH")
anova_results$FDR_TempxTrt    <- p.adjust(anova_results$p_TempxTrt, "BH")

anova_results <- anova_results |> dplyr::arrange(pmin(FDR_Temperature, FDR_Treatment, FDR_Time, na.rm = TRUE))
write.csv(anova_results, file.path(dir_mt, "Module_ANOVA_results.csv"), row.names = FALSE)

cat("ANOVA summary:\n")
cat("  Modules sig Temperature (FDR<0.05):", sum(anova_results$FDR_Temperature < FDR_CUTOFF, na.rm = TRUE), "\n")
cat("  Modules sig Treatment   (FDR<0.05):", sum(anova_results$FDR_Treatment < FDR_CUTOFF, na.rm = TRUE), "\n")
cat("  Modules sig Time        (FDR<0.05):", sum(anova_results$FDR_Time < FDR_CUTOFF, na.rm = TRUE), "\n")
cat("  Modules sig Temp×Trt    (FDR<0.05):", sum(anova_results$FDR_TempxTrt < FDR_CUTOFF, na.rm = TRUE), "\n")

# ── Figure: ANOVA heatmap ──
anova_pmat <- anova_results |>
  dplyr::select(ModuleColor, FDR_Temperature, FDR_Treatment, FDR_Time, FDR_TempxTrt) |>
  tibble::column_to_rownames("ModuleColor")
colnames(anova_pmat) <- c("Temperature", "Treatment", "Time (P3/P4)", "Temp × Trt")

log_pmat <- -log10(as.matrix(anova_pmat))
log_pmat[is.infinite(log_pmat)] <- 10
log_pmat[is.na(log_pmat)] <- 0

stars_anova <- apply(as.matrix(anova_pmat), c(1,2), starify)
text_anova <- paste0(sprintf("%.2f", as.matrix(anova_pmat)), "\n", stars_anova)
dim(text_anova) <- dim(log_pmat)

save_base_plot(function(){
  par(mar = c(10, 10, 4, 2))
  WGCNA::labeledHeatmap(
    Matrix = log_pmat,
    xLabels = colnames(log_pmat),
    yLabels = rownames(log_pmat),
    ySymbols = rownames(log_pmat),
    colorLabels = FALSE,
    colors = colorRampPalette(c("white","#FFF9C4","#FFAB00","#D50000"))(50),
    textMatrix = text_anova,
    setStdMargins = FALSE,
    cex.text = 0.50,
    zlim = c(0, max(log_pmat, na.rm = TRUE)),
    xLabelsAngle = 45,
    main = paste0("ANOVA — ", GENO, " (P3+P4) LC-MS\n-log10(FDR) | lm(ME ~ Temp * Trt + Time)")
  )
}, file.path(dir_mt, "ANOVA_summary_heatmap"),
w = 9, h = max(7, nrow(log_pmat) * 0.30 + 3))


##############################################################
## PARTIE 5: BOXPLOTS
##############################################################
cat("── PARTIE 5: Boxplots ──\n")

sig_any <- anova_results |>
  dplyr::filter(FDR_Temperature < FDR_CUTOFF | FDR_Treatment < FDR_CUTOFF |
                  FDR_Time < FDR_CUTOFF | FDR_TempxTrt < FDR_CUTOFF)

if (nrow(sig_any) == 0) {
  sig_any <- anova_results |> dplyr::slice_head(n = 8)
  cat("  Aucun module FDR<0.05 → plotting top 8 modules\n")
}

for (i in seq_len(nrow(sig_any))) {
  me   <- sig_any$ME[i]
  colr <- sig_any$ModuleColor[i]
  
  df <- data.frame(
    ME          = as.numeric(MEs_clean[, me]),
    Temperature = meta$Temperature,
    Treatment   = meta$Treatment,
    Time        = meta$Sampling_time,
    Group       = meta$Group
  )
  
  pT  <- sig_any$FDR_Temperature[i]
  pTr <- sig_any$FDR_Treatment[i]
  pTi <- sig_any$FDR_Time[i]
  pI  <- sig_any$FDR_TempxTrt[i]
  
  sub_txt <- paste0(
    "ANOVA FDR: Temp=", sprintf("%.2g", pT),
    starify(pT), " | Trt=", sprintf("%.2g", pTr),
    starify(pTr), " | Time=", sprintf("%.2g", pTi),
    starify(pTi), " | T×Trt=", sprintf("%.2g", pI), starify(pI)
  )
  
  # 5A) Boxplot par Température
  p_t <- ggplot2::ggplot(df, ggplot2::aes(Temperature, ME, fill = Temperature)) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6) +
    ggplot2::geom_jitter(width = 0.15, size = 1.5, alpha = 0.6) +
    ggplot2::scale_fill_manual(values = TEMP_COLORS) +
    ggplot2::labs(title = "par Température", y = "ME") +
    theme_pub(base_size = 9) + ggplot2::theme(legend.position = "none")
  
  # 5B) Boxplot par Traitement
  p_tr <- ggplot2::ggplot(df, ggplot2::aes(Treatment, ME, fill = Treatment)) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6) +
    ggplot2::geom_jitter(width = 0.15, size = 1.5, alpha = 0.6) +
    ggplot2::scale_fill_manual(values = TRT_COLORS) +
    ggplot2::labs(title = "par Traitement", y = NULL) +
    theme_pub(base_size = 9) + ggplot2::theme(legend.position = "none")
  
  # 5C) Boxplot par Time
  p_ti <- ggplot2::ggplot(df, ggplot2::aes(Time, ME, fill = Time)) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6) +
    ggplot2::geom_jitter(width = 0.15, size = 1.5, alpha = 0.6) +
    ggplot2::scale_fill_manual(values = TIME_COLORS) +
    ggplot2::labs(title = "par Time point", y = NULL) +
    theme_pub(base_size = 9) + ggplot2::theme(legend.position = "none")
  
  # 5D) Boxplot SDP vs Control dans chaque T°
  p_sdp <- ggplot2::ggplot(df, ggplot2::aes(Treatment, ME, fill = Treatment)) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    ggplot2::geom_jitter(width = 0.12, size = 1.5, alpha = 0.6) +
    ggplot2::facet_wrap(~ Temperature, labeller = ggplot2::labeller(Temperature = c("T0"="NH","T1"="HS-7","T2"="HS-2")), nrow = 1) +
    ggplot2::scale_fill_manual(values = TRT_COLORS) +
    ggplot2::labs(title = "Treated vs Control by Temperature", y = "ME") +
    theme_pub(base_size = 9) + ggplot2::theme(legend.position = "bottom")
  
  # 5E) Boxplot par groupe complet (Temp × Trt × Time)
  p_full <- ggplot2::ggplot(df, ggplot2::aes(Treatment, ME, fill = Treatment)) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    ggplot2::geom_jitter(width = 0.12, size = 1.2, alpha = 0.5) +
    ggplot2::facet_grid(Time ~ Temperature,
                        labeller = ggplot2::labeller(
                          Temperature = c("T0"="NH","T1"="HS-7","T2"="HS-2"),
                          Time = c("P3"="P3 (48h)","P4"="P4 (10d)")
                        )) +
    ggplot2::scale_fill_manual(values = TRT_COLORS) +
    ggplot2::labs(title = "Détail: Trt × Temp × Time", y = "ME") +
    theme_pub(base_size = 9) + ggplot2::theme(legend.position = "bottom")
  
  # Assembler
  p_top <- p_t + p_tr + p_ti + patchwork::plot_layout(ncol = 3)
  p_combo <- (p_top / p_sdp / p_full) +
    patchwork::plot_annotation(
      title = paste0("Module ", colr, " (", me, ") — ", GENO, " LC-MS"),
      subtitle = sub_txt,
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 14),
        plot.subtitle = ggplot2::element_text(size = 9, colour = "grey40")
      )
    )
  
  save_plot(p_combo, file.path(dir_mt_box, paste0(colr, "_", me, "_boxplots")),
            w = 10, h = 14)
  
  # ── Barplot moyen ME par groupe ──
  mean_df <- df |>
    dplyr::group_by(Temperature, Treatment, Time) |>
    dplyr::summarise(mean_ME = mean(ME), se_ME = sd(ME) / sqrt(dplyr::n()), .groups = "drop") |>
    dplyr::mutate(Group = paste(Temperature, Treatment, sep = "_"))
  
  p_bar <- ggplot2::ggplot(mean_df, ggplot2::aes(Group, mean_ME, fill = Treatment)) +
    ggplot2::geom_col(position = "dodge", width = 0.7, colour = "white", linewidth = 0.3) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = mean_ME - se_ME, ymax = mean_ME + se_ME),
                           width = 0.2, position = ggplot2::position_dodge(0.7)) +
    ggplot2::facet_wrap(~ Time, nrow = 1,
                        labeller = ggplot2::labeller(Time = c("P3"="P3 (48h)","P4"="P4 (10d)"))) +
    ggplot2::scale_fill_manual(values = TRT_COLORS) +
    ggplot2::labs(title = paste0(colr, " — Mean ME ± SE"),
                  x = NULL, y = "Module eigengene") +
    theme_pub(base_size = 10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  
  save_plot(p_bar, file.path(dir_mt_bar, paste0(colr, "_", me, "_barplot")), w = 9, h = 5)
}

cat("  Boxplots saved for", nrow(sig_any), "modules\n")


##############################################################
## PARTIE 6: TESTS PAIRÉS — Wilcoxon
##############################################################
cat("── PARTIE 6: Tests pairés ──\n")

pairwise_results <- data.frame()

for (me in colnames(MEs_clean)) {
  colr <- as.character(me_to_color[me])
  df <- data.frame(ME = as.numeric(MEs_clean[, me]),
                   Temperature = meta$Temperature,
                   Treatment = meta$Treatment,
                   Time = meta$Sampling_time)
  
  # SDP vs Control dans chaque T°
  for (temp in c("T0","T1","T2")) {
    sub <- df |> dplyr::filter(Temperature == temp)
    if (length(unique(sub$Treatment)) < 2 || nrow(sub) < 4) next
    pv <- tryCatch(
      wilcox.test(ME ~ Treatment, data = sub, exact = FALSE)$p.value,
      error = function(e) NA
    )
    pairwise_results <- dplyr::bind_rows(pairwise_results, data.frame(
      ME = me, ModuleColor = colr, Test = paste0("Treated_vs_Control_in_", temp),
      p_value = pv
    ))
  }
  
  # SDP vs Control dans chaque Time
  for (tp in c("P3","P4")) {
    sub <- df |> dplyr::filter(Time == tp)
    if (length(unique(sub$Treatment)) < 2 || nrow(sub) < 4) next
    pv <- tryCatch(
      wilcox.test(ME ~ Treatment, data = sub, exact = FALSE)$p.value,
      error = function(e) NA
    )
    pairwise_results <- dplyr::bind_rows(pairwise_results, data.frame(
      ME = me, ModuleColor = colr, Test = paste0("Treated_vs_Control_in_", tp),
      p_value = pv
    ))
  }
  
  # T1 vs T0, T2 vs T0 (global)
  for (ctr in list(c("T0","T1"), c("T0","T2"))) {
    sub <- df |> dplyr::filter(Temperature %in% ctr)
    if (nrow(sub) < 4) next
    pv <- tryCatch(
      wilcox.test(ME ~ Temperature, data = sub, exact = FALSE)$p.value,
      error = function(e) NA
    )
    pairwise_results <- dplyr::bind_rows(pairwise_results, data.frame(
      ME = me, ModuleColor = colr, Test = paste0(ctr[2], "_vs_", ctr[1]),
      p_value = pv
    ))
  }
}

pairwise_results <- pairwise_results |>
  dplyr::group_by(Test) |>
  dplyr::mutate(FDR = p.adjust(p_value, "BH")) |>
  dplyr::ungroup()

write.csv(pairwise_results, file.path(dir_mt, "Pairwise_tests_results.csv"), row.names = FALSE)


##############################################################
## PARTIE 7: CLASSIFICATION DES MODULES
##############################################################
cat("── PARTIE 7: Classification ──\n")

classify <- anova_results |>
  dplyr::mutate(
    Temp_responsive = FDR_Temperature < FDR_CUTOFF,
    Trt_responsive  = FDR_Treatment < FDR_CUTOFF,
    Time_responsive = FDR_Time < FDR_CUTOFF,
    Interaction     = FDR_TempxTrt < FDR_CUTOFF,
    
    Classification = dplyr::case_when(
      Interaction & Temp_responsive & Trt_responsive ~ "Interaction Temp×Treated",
      Interaction ~ "Interaction Temp×Treated",
      Temp_responsive & Trt_responsive & Time_responsive ~ "Multi-responsive",
      Temp_responsive & Trt_responsive ~ "Temp + Treated responsive",
      Temp_responsive & Time_responsive ~ "Temp + Time responsive",
      Trt_responsive & Time_responsive  ~ "Treated + Time responsive",
      Temp_responsive ~ "Temperature-responsive",
      Trt_responsive  ~ "Treated-responsive",
      Time_responsive ~ "Time-responsive (P3→P4)",
      TRUE ~ "Stable"
    )
  ) |>
  dplyr::select(ME, ModuleColor, Classification,
                FDR_Temperature, FDR_Treatment, FDR_Time, FDR_TempxTrt, R2)

write.csv(classify, file.path(dir_mt, "Module_Classification.csv"), row.names = FALSE)

cat("\nModule classification:\n")
print(table(classify$Classification))

# ── Figure: classification summary ──
class_count <- as.data.frame(table(classify$Classification)) |>
  setNames(c("Class", "N")) |>
  dplyr::arrange(dplyr::desc(N))

class_colors <- c(
  "Temperature-responsive" = "#D55E00",
  "Treated-responsive" = "#CC79A7",
  "Time-responsive (P3→P4)" = "#009E73",
  "Interaction Temp×Treated" = "#E69F00",
  "Temp + Treated responsive" = "#F0E442",
  "Multi-responsive" = "#56B4E9",
  "Temp + Time responsive" = "#0072B2",
  "Treated + Time responsive" = "#882255",
  "Stable" = "#999999"
)

p_class <- ggplot2::ggplot(class_count,
                           ggplot2::aes(x = reorder(Class, N), y = N, fill = Class)) +
  ggplot2::geom_col(width = 0.7) +
  ggplot2::geom_text(ggplot2::aes(label = N), hjust = -0.2, fontface = "bold", size = 3.5) +
  ggplot2::scale_fill_manual(values = class_colors, guide = "none") +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.15))) +
  ggplot2::coord_flip() +
  ggplot2::labs(title = paste0("Module classification — ", GENO, " LC-MS"),
                subtitle = paste0("ANOVA lm(ME ~ Temp * Trt + Time) | FDR < ", FDR_CUTOFF),
                x = NULL, y = "Number of modules") +
  theme_pub()

save_plot(p_class, file.path(dir_mt, "Module_Classification_barplot"), w = 9, h = 5)


log_msg("Script 04 DONE — ", nrow(anova_results), " modules analysés")
cat("\n✓ Script 04 terminé\n")
cat("  Heatmaps + ANOVA + boxplots + classification dans:", dir_mt, "\n")