############################################################
## 05_hub_genes.R
## kME + Gene Significance → Hub genes → kME vs GS scatter
# Repository configuration (portable; no HPC-specific paths)
if (!exists('CARROT_CONFIG_LOADED')) source('../../../config/config.R')


############################################################

##############################
## PARAMÈTRES
##############################
INSTALL_MISSING_PKGS <- FALSE
KME_THRESHOLD  <- 0.70
GS_THRESHOLD   <- 0.20
TOP_HUB_N      <- 30
FDR_CUTOFF     <- 0.05

base_out <- RNA_WGCNA_ROOT

##############################
## SOURCE + LOAD
##############################
source("00_helpers.R")
load_pkgs(c("svglite","dplyr","ggplot2","ggrepel"), install_missing = INSTALL_MISSING_PKGS)
load_pkgs(c("WGCNA"), install_missing = INSTALL_MISSING_PKGS, bioc = TRUE)

GENOTYPE_CHOSEN <- Sys.getenv("CARROT_GENOTYPE", unset = "PRESTO")
GENO_TAG <- gsub("[^A-Za-z0-9]+", "", GENOTYPE_CHOSEN)
run_dir <- find_latest_dir(base_out, paste0("^WGCNA_", GENO_TAG, "_"))
net      <- readRDS(file.path(run_dir, "02_network", "net.rds"))
MEs      <- readRDS(file.path(run_dir, "02_network", "MEs.rds"))
datExpr  <- readRDS(file.path(run_dir, "02_network", "datExpr_clean.rds"))
meta     <- readRDS(file.path(run_dir, "02_network", "meta_clean.rds"))
me_map   <- read.csv(file.path(run_dir, "02_network", "ME_color_map.csv"))
gene_mod <- read.csv(file.path(run_dir, "02_network", "Gene_to_Module.csv"))
anova_res <- read.csv(file.path(run_dir, "03_module_trait", "Module_ANOVA_results.csv"))

GENO <- as.character(meta$Genotype[1])

dir_hubs <- file.path(run_dir, "04_hub_genes")
dir_hubs_plots <- file.path(dir_hubs, "scatter_plots")
dir_hubs_mod <- file.path(dir_hubs, "by_module")
for (d in c(dir_hubs, dir_hubs_plots, dir_hubs_mod)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

log_msg <- init_logger(file.path(run_dir, "pipeline.log"))
log_msg("START Script 05 — Hub genes")
cat("── Script 05: Hub genes —", GENO, "──\n\n")

##############################
## PARTIE 1: Calcul kME
##############################
cat("── kME ──\n")

kME_mat <- WGCNA::signedKME(datExpr, MEs, corFnc = "bicor")
# Harmoniser les noms de colonnes
colnames(kME_mat) <- gsub("^kME", "ME", colnames(kME_mat))

# kME propre (vers son propre module)
kME_own <- rep(NA_real_, ncol(datExpr))
names(kME_own) <- colnames(datExpr)

for (i in seq_len(nrow(gene_mod))) {
  g <- gene_mod$Gene[i]
  me_col <- paste0("ME", gene_mod$ModuleLabel[i])
  if (g %in% rownames(kME_mat) && me_col %in% colnames(kME_mat)) {
    kME_own[g] <- kME_mat[g, me_col]
  }
}

##############################
## PARTIE 2: Gene Significance (GS) par facteur
##############################
cat("── GS ──\n")

# GS vs Température (numérique: NH=0, HS-7=1, HS-2=2)
temp_num <- as.numeric(meta$Temperature) - 1  # 0,1,2
GS_temp <- as.numeric(WGCNA::bicor(datExpr, as.matrix(temp_num), use = "pairwise.complete.obs"))
names(GS_temp) <- colnames(datExpr)

# GS vs Traitement (SDP=1, EAU=0)
trt_num <- as.numeric(meta$Treatment == "Treated")
GS_trt <- as.numeric(WGCNA::bicor(datExpr, as.matrix(trt_num), use = "pairwise.complete.obs"))
names(GS_trt) <- colnames(datExpr)

# GS vs Time (P3=1, P2=0)
time_num <- as.numeric(meta$Sampling_time == "P3")
GS_time <- as.numeric(WGCNA::bicor(datExpr, as.matrix(time_num), use = "pairwise.complete.obs"))
names(GS_time) <- colnames(datExpr)

##############################
## PARTIE 3: Assembler table hub genes
##############################
cat("── Hub table ──\n")

hub_df <- gene_mod |>
  dplyr::mutate(
    kME = kME_own[Gene],
    GS_Temperature = GS_temp[Gene],
    GS_Treatment   = GS_trt[Gene],
    GS_Time        = GS_time[Gene],
    IsHub = !is.na(kME) & kME >= KME_THRESHOLD & ModuleColor != "grey",
    IsHighGS_Temp = abs(GS_Temperature) >= GS_THRESHOLD,
    IsHighGS_Trt  = abs(GS_Treatment) >= GS_THRESHOLD
  )

cat("  Hub genes (kME >=", KME_THRESHOLD, "):", sum(hub_df$IsHub, na.rm = TRUE), "\n")
cat("  High GS Temperature:", sum(hub_df$IsHighGS_Temp, na.rm = TRUE), "\n")
cat("  High GS Treatment:", sum(hub_df$IsHighGS_Trt, na.rm = TRUE), "\n")

write.csv(hub_df, file.path(dir_hubs, "AllGenes_kME_GS.csv"), row.names = FALSE)

# Hub genes prioritaires
hub_priority <- hub_df |> dplyr::filter(IsHub & (IsHighGS_Temp | IsHighGS_Trt))
write.csv(hub_priority, file.path(dir_hubs, "HubGenes_PRIORITY.csv"), row.names = FALSE)
cat("  Priority hubs (Hub + high GS):", nrow(hub_priority), "\n\n")

# Top hub par module
top_per_mod <- hub_df |>
  dplyr::filter(ModuleColor != "grey") |>
  dplyr::group_by(ModuleColor) |>
  dplyr::arrange(dplyr::desc(kME), .by_group = TRUE) |>
  dplyr::slice_head(n = TOP_HUB_N) |>
  dplyr::ungroup()

write.csv(top_per_mod, file.path(dir_hubs, paste0("Top", TOP_HUB_N, "_HubGenes_perModule.csv")),
          row.names = FALSE)

# Exports par module
for (mc in unique(gene_mod$ModuleColor[gene_mod$ModuleColor != "grey"])) {
  mod_df <- hub_df |> dplyr::filter(ModuleColor == mc) |> dplyr::arrange(dplyr::desc(kME))
  write.csv(mod_df, file.path(dir_hubs_mod, paste0(mc, "_genes.csv")), row.names = FALSE)
}

##############################
## PARTIE 4: Scatter kME vs GS — modules significatifs
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
  
  # Déterminer quel GS est le plus pertinent pour ce module
  pT <- sig_modules$FDR_Temperature[i]
  pTr <- sig_modules$FDR_Treatment[i]
  
  if (!is.na(pT) && (is.na(pTr) || pT < pTr)) {
    gs_col <- "GS_Temperature"
    gs_label <- "GS ~ Temperature (NH→HS-7→HS-2)"
  } else {
    gs_col <- "GS_Treatment"
    gs_label <- "GS ~ Treatment (Treated)"
  }
  
  mod_df$GS_plot <- mod_df[[gs_col]]
  
  # Points colorés par catégorie
  mod_df$Category <- dplyr::case_when(
    mod_df$IsHub & abs(mod_df$GS_plot) >= GS_THRESHOLD ~ "Hub + high GS",
    mod_df$IsHub ~ "Hub",
    abs(mod_df$GS_plot) >= GS_THRESHOLD ~ "High GS",
    TRUE ~ "Other"
  )
  
  cat_colors <- c("Hub + high GS" = "#C0392B", "Hub" = "#E8820C",
                  "High GS" = "#2471A3", "Other" = "grey80")
  
  # Labels top genes
  top_lab <- mod_df |>
    dplyr::filter(Category == "Hub + high GS") |>
    dplyr::arrange(dplyr::desc(kME)) |>
    dplyr::slice_head(n = 15)
  
  p_sc <- ggplot2::ggplot(mod_df, ggplot2::aes(kME, GS_plot, colour = Category)) +
    ggplot2::geom_point(size = ifelse(mod_df$Category == "Other", 0.8, 2.5),
                        alpha = ifelse(mod_df$Category == "Other", 0.3, 0.85)) +
    ggrepel::geom_text_repel(data = top_lab, ggplot2::aes(label = Gene),
                             size = 2.2, fontface = "italic", colour = "#C0392B",
                             max.overlaps = 20, segment.colour = "grey55") +
    ggplot2::geom_vline(xintercept = KME_THRESHOLD, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_hline(yintercept = c(-GS_THRESHOLD, GS_THRESHOLD),
                        linetype = "dashed", colour = "grey50") +
    ggplot2::scale_colour_manual(values = cat_colors) +
    ggplot2::labs(title = paste0("Module ", colr, " — kME vs ", gs_label),
                  subtitle = paste0(GENO, " | n=", nrow(mod_df), " genes"),
                  x = "Module Membership (kME)", y = gs_label) +
    theme_pub()
  
  save_plot(p_sc, file.path(dir_hubs_plots, paste0(colr, "_kME_vs_GS")), w = 8, h = 6)
}

log_msg("Script 05 DONE — ", nrow(hub_priority), " priority hubs")
cat("\n✓ Script 05 terminé\n")