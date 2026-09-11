############################################################
## 05_hub_metabolites.R  (LC-MS WGCNA)
## kME + Metabolite Significance → Hub metabolites → scatter
# Repository configuration (portable; no HPC-specific paths)
if (!exists('CARROT_CONFIG_LOADED')) source('../../../config/config.R')


############################################################
## DIFFÉRENCES vs RNAseq:
##  - "Hub genes" → "Hub metabolites"
##  - Feature IDs = m/z_RT au lieu de gene IDs
##  - GS = Metabolite Significance (MS)
##  - P3/P4 au lieu de P2/P3
############################################################

##############################
## PARAMÈTRES
##############################
INSTALL_MISSING_PKGS <- FALSE
KME_THRESHOLD  <- 0.70
MS_THRESHOLD   <- 0.20    # Metabolite Significance (equiv GS)
TOP_HUB_N      <- 30
FDR_CUTOFF     <- 0.05

base_out <- LCMS_WGCNA_ROOT

##############################
## SOURCE + LOAD
##############################
source("00_helpers.R")
load_pkgs(c("svglite","dplyr","ggplot2","ggrepel"), install_missing = INSTALL_MISSING_PKGS)
load_pkgs(c("WGCNA"), install_missing = INSTALL_MISSING_PKGS, bioc = TRUE)

run_dir <- find_latest_dir(base_out, "^WGCNA_LCMS_")
net      <- readRDS(file.path(run_dir, "02_network", "net.rds"))
MEs      <- readRDS(file.path(run_dir, "02_network", "MEs.rds"))
datExpr  <- readRDS(file.path(run_dir, "02_network", "datExpr_clean.rds"))
meta     <- readRDS(file.path(run_dir, "02_network", "meta_clean.rds"))
me_map   <- read.csv(file.path(run_dir, "02_network", "ME_color_map.csv"))
feature_mod <- read.csv(file.path(run_dir, "02_network", "Feature_to_Module.csv"))
anova_res   <- read.csv(file.path(run_dir, "03_module_trait", "Module_ANOVA_results.csv"))

# Charger feature_info si disponible (pour les annotations m/z, RT, Name)
feature_info <- tryCatch(
  readRDS(file.path(run_dir, "feature_info.rds")),
  error = function(e) NULL
)

GENO <- as.character(meta$Genotype[1])

dir_hubs <- file.path(run_dir, "04_hub_metabolites")
dir_hubs_plots <- file.path(dir_hubs, "scatter_plots")
dir_hubs_mod <- file.path(dir_hubs, "by_module")
for (d in c(dir_hubs, dir_hubs_plots, dir_hubs_mod)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

log_msg <- init_logger(file.path(run_dir, "pipeline.log"))
log_msg("START Script 05 — Hub metabolites (LC-MS)")
cat("── Script 05: Hub metabolites —", GENO, "──\n\n")

##############################
## PARTIE 1: Calcul kME
##############################
cat("── kME ──\n")

kME_mat <- WGCNA::signedKME(datExpr, MEs, corFnc = "bicor")
colnames(kME_mat) <- gsub("^kME", "ME", colnames(kME_mat))

# kME propre (vers son propre module)
kME_own <- rep(NA_real_, ncol(datExpr))
names(kME_own) <- colnames(datExpr)

for (i in seq_len(nrow(feature_mod))) {
  f <- feature_mod$Feature[i]
  me_col <- paste0("ME", feature_mod$ModuleLabel[i])
  if (f %in% rownames(kME_mat) && me_col %in% colnames(kME_mat)) {
    kME_own[f] <- kME_mat[f, me_col]
  }
}

##############################
## PARTIE 2: Metabolite Significance (MS) par facteur
##############################
cat("── MS (Metabolite Significance) ──\n")

# MS vs Température (numérique: T0=0, T1=1, T2=2)
temp_num <- as.numeric(meta$Temperature) - 1
MS_temp <- as.numeric(WGCNA::bicor(datExpr, as.matrix(temp_num), use = "pairwise.complete.obs"))
names(MS_temp) <- colnames(datExpr)

# MS vs Traitement (Treated=1, Control=0)
trt_num <- as.numeric(meta$Treatment == "Treated")
MS_trt <- as.numeric(WGCNA::bicor(datExpr, as.matrix(trt_num), use = "pairwise.complete.obs"))
names(MS_trt) <- colnames(datExpr)

# MS vs Time (P4=1, P3=0)
time_num <- as.numeric(meta$Sampling_time == "P4")
MS_time <- as.numeric(WGCNA::bicor(datExpr, as.matrix(time_num), use = "pairwise.complete.obs"))
names(MS_time) <- colnames(datExpr)

##############################
## PARTIE 3: Assembler table hub metabolites
##############################
cat("── Hub table ──\n")

hub_df <- feature_mod |>
  dplyr::mutate(
    kME = kME_own[Feature],
    MS_Temperature = MS_temp[Feature],
    MS_Treatment   = MS_trt[Feature],
    MS_Time        = MS_time[Feature],
    IsHub = !is.na(kME) & kME >= KME_THRESHOLD & ModuleColor != "grey",
    IsHighMS_Temp = abs(MS_Temperature) >= MS_THRESHOLD,
    IsHighMS_Trt  = abs(MS_Treatment) >= MS_THRESHOLD
  )

# Ajouter les annotations m/z, RT, Name si disponibles
if (!is.null(feature_info) && "FeatureID" %in% colnames(feature_info)) {
  # Mapper les infos
  anno_cols <- intersect(colnames(feature_info),
                         c("FeatureID","Name","m/z","mz","RT[min]","RT","RT.min.","ID_MT","Tags"))
  if (length(anno_cols) > 1) {
    hub_df <- hub_df |>
      dplyr::left_join(feature_info[, anno_cols, drop = FALSE],
                       by = c("Feature" = "FeatureID"))
  }
}

cat("  Hub metabolites (kME >=", KME_THRESHOLD, "):", sum(hub_df$IsHub, na.rm = TRUE), "\n")
cat("  High MS Temperature:", sum(hub_df$IsHighMS_Temp, na.rm = TRUE), "\n")
cat("  High MS Treatment:", sum(hub_df$IsHighMS_Trt, na.rm = TRUE), "\n")

write.csv(hub_df, file.path(dir_hubs, "AllFeatures_kME_MS.csv"), row.names = FALSE)

# Hub metabolites prioritaires
hub_priority <- hub_df |> dplyr::filter(IsHub & (IsHighMS_Temp | IsHighMS_Trt))
write.csv(hub_priority, file.path(dir_hubs, "HubMetabolites_PRIORITY.csv"), row.names = FALSE)
cat("  Priority hubs (Hub + high MS):", nrow(hub_priority), "\n\n")

# Top hub par module
top_per_mod <- hub_df |>
  dplyr::filter(ModuleColor != "grey") |>
  dplyr::group_by(ModuleColor) |>
  dplyr::arrange(dplyr::desc(kME), .by_group = TRUE) |>
  dplyr::slice_head(n = TOP_HUB_N) |>
  dplyr::ungroup()

write.csv(top_per_mod, file.path(dir_hubs, paste0("Top", TOP_HUB_N, "_HubMetabolites_perModule.csv")),
          row.names = FALSE)

# Exports par module
for (mc in unique(feature_mod$ModuleColor[feature_mod$ModuleColor != "grey"])) {
  mod_df <- hub_df |> dplyr::filter(ModuleColor == mc) |> dplyr::arrange(dplyr::desc(kME))
  write.csv(mod_df, file.path(dir_hubs_mod, paste0(mc, "_metabolites.csv")), row.names = FALSE)
}

##############################
## PARTIE 4: Scatter kME vs MS — modules significatifs
##############################
cat("── Scatter plots ──\n")

sig_modules <- anova_res |>
  dplyr::filter(FDR_Temperature < FDR_CUTOFF | FDR_Treatment < FDR_CUTOFF |
                  FDR_Time < FDR_CUTOFF | FDR_TempxTrt < FDR_CUTOFF)

if (nrow(sig_modules) == 0) sig_modules <- anova_res |> dplyr::slice_head(n = 6)

for (i in seq_len(nrow(sig_modules))) {
  me   <- sig_modules$ME[i]
  colr <- sig_modules$ModuleColor[i]
  
  mod_df <- hub_df |> dplyr::filter(ModuleColor == colr)
  if (nrow(mod_df) < 5) next
  
  pT <- sig_modules$FDR_Temperature[i]
  pTr <- sig_modules$FDR_Treatment[i]
  
  if (!is.na(pT) && (is.na(pTr) || pT < pTr)) {
    ms_col <- "MS_Temperature"
    ms_label <- "MS ~ Temperature (NH→HS-7→HS-2)"
  } else {
    ms_col <- "MS_Treatment"
    ms_label <- "MS ~ Treatment (Treated)"
  }
  
  mod_df$MS_plot <- mod_df[[ms_col]]
  
  mod_df$Category <- dplyr::case_when(
    mod_df$IsHub & abs(mod_df$MS_plot) >= MS_THRESHOLD ~ "Hub + high MS",
    mod_df$IsHub ~ "Hub",
    abs(mod_df$MS_plot) >= MS_THRESHOLD ~ "High MS",
    TRUE ~ "Other"
  )
  
  cat_colors <- c("Hub + high MS" = "#C0392B", "Hub" = "#E8820C",
                  "High MS" = "#2471A3", "Other" = "grey80")
  
  # Labels top metabolites
  top_lab <- mod_df |>
    dplyr::filter(Category == "Hub + high MS") |>
    dplyr::arrange(dplyr::desc(kME)) |>
    dplyr::slice_head(n = 15)
  
  p_sc <- ggplot2::ggplot(mod_df, ggplot2::aes(kME, MS_plot, colour = Category)) +
    ggplot2::geom_point(size = ifelse(mod_df$Category == "Other", 0.8, 2.5),
                        alpha = ifelse(mod_df$Category == "Other", 0.3, 0.85)) +
    ggrepel::geom_text_repel(data = top_lab, ggplot2::aes(label = Feature),
                             size = 2.0, fontface = "italic", colour = "#C0392B",
                             max.overlaps = 20, segment.colour = "grey55") +
    ggplot2::geom_vline(xintercept = KME_THRESHOLD, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_hline(yintercept = c(-MS_THRESHOLD, MS_THRESHOLD),
                        linetype = "dashed", colour = "grey50") +
    ggplot2::scale_colour_manual(values = cat_colors) +
    ggplot2::labs(title = paste0("Module ", colr, " — kME vs ", ms_label),
                  subtitle = paste0(GENO, " | n=", nrow(mod_df), " metabolites (LC-MS)"),
                  x = "Module Membership (kME)", y = ms_label) +
    theme_pub()
  
  save_plot(p_sc, file.path(dir_hubs_plots, paste0(colr, "_kME_vs_MS")), w = 8, h = 6)
}

log_msg("Script 05 DONE — ", nrow(hub_priority), " priority hub metabolites")
cat("\n✓ Script 05 terminé\n")