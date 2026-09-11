############################################################
## 03_network_construction.R  (LC-MS WGCNA)
## Soft threshold → blockwiseModules → dendrogram + sizes
# Repository configuration (portable; no HPC-specific paths)
if (!exists('CARROT_CONFIG_LOADED')) source('../../../config/config.R')


############################################################
## DIFFÉRENCES vs RNAseq:
##  - Moins de features → maxBlockSize plus petit mais single block OK
##  - SFT_R2_TARGET peut être ajusté (métabo souvent plus bruité)
##  - MIN_MOD_SIZE réduit (moins de features total)
############################################################

##############################
## PARAMETERS
##############################
INSTALL_MISSING_PKGS <- FALSE

NETWORK_TYPE   <- "signed"
COR_TYPE       <- "bicor"
SFT_R2_TARGET  <- 0.7
DEEP_SPLIT     <- 3
MIN_MOD_SIZE   <- 100     # ★ réduit vs RNAseq (100) car moins de features
MERGE_CUT_H    <- 0.25

base_out <- LCMS_WGCNA_ROOT

##############################
## SOURCE + LOAD
##############################
source("00_helpers.R")
load_pkgs(c("dplyr","ggplot2"), install_missing = INSTALL_MISSING_PKGS)
load_pkgs(c("WGCNA"), install_missing = INSTALL_MISSING_PKGS, bioc = TRUE)
WGCNA::allowWGCNAThreads()

run_dir <- find_latest_dir(base_out, "^WGCNA_LCMS_")
datExpr <- readRDS(file.path(run_dir, "datExpr.rds"))
meta    <- readRDS(file.path(run_dir, "meta.rds"))

GENO <- as.character(meta$Genotype[1])
dir_net <- file.path(run_dir, "02_network")
dir.create(dir_net, recursive = TRUE, showWarnings = FALSE)

log_msg <- init_logger(file.path(run_dir, "pipeline.log"))
log_msg("START Script 03 — Network construction (LC-MS)")
cat("── Script 03: Network —", GENO, "──\n")
cat("  datExpr:", nrow(datExpr), "samples x", ncol(datExpr), "features\n\n")

##############################
## PART 1: QC samples/features
##############################
cat("── QC WGCNA ──\n")
gsg <- WGCNA::goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  if (sum(!gsg$goodSamples) > 0) {
    cat("  Removing", sum(!gsg$goodSamples), "bad samples\n")
    datExpr <- datExpr[gsg$goodSamples, , drop = FALSE]
    meta <- meta[match(rownames(datExpr), meta$Sample), ]
  }
  if (sum(!gsg$goodGenes) > 0) {
    cat("  Removing", sum(!gsg$goodGenes), "bad features\n")
    datExpr <- datExpr[, gsg$goodGenes, drop = FALSE]
  }
}

# Remove remaining constant features
mad_v <- apply(datExpr, 2, mad, na.rm = TRUE)
var_v <- apply(datExpr, 2, var, na.rm = TRUE)
bad <- is.na(mad_v) | mad_v == 0 | is.na(var_v) | var_v == 0
if (sum(bad) > 0) {
  cat("  Removing", sum(bad), "constant features\n")
  datExpr <- datExpr[, !bad, drop = FALSE]
}

log_msg("QC done: ", nrow(datExpr), " samples x ", ncol(datExpr), " features")

##############################
## PART 2: Soft threshold
##############################
cat("── SOFT THRESHOLD ──\n")
Powers <- 1:20
sft <- WGCNA::pickSoftThreshold(
  datExpr, powerVector = Powers,
  networkType = NETWORK_TYPE,
  corFnc = COR_TYPE, verbose = 3
)

fit_r2 <- -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2]
ok <- which(fit_r2 >= SFT_R2_TARGET)
softPower_auto <- if (length(ok) > 0) Powers[min(ok)] else Powers[which.max(fit_r2)]
softPower <- if (tolower(Sys.getenv("CARROT_PUBLICATION_MODE", unset="false")) == "true") LCMS_SOFT_POWER else softPower_auto
log_msg("Soft power used: ", softPower, " | diagnostic auto-choice: ", softPower_auto)

log_msg("Soft power: ", softPower, " | max R²=", round(max(fit_r2), 3))

sft_df <- data.frame(Power = Powers, SFT_R2 = fit_r2, MeanK = sft$fitIndices[, 5])
write.csv(sft_df, file.path(dir_net, "SoftThreshold_fitIndices.csv"), row.names = FALSE)

# Figure
p_sft <- ggplot2::ggplot(sft_df, ggplot2::aes(Power, SFT_R2)) +
  ggplot2::geom_line(colour = "grey70") +
  ggplot2::geom_point(size = 2.5, colour = "#2471A3") +
  ggplot2::geom_text(ggplot2::aes(label = Power), vjust = -0.8, size = 2.8) +
  ggplot2::geom_hline(yintercept = SFT_R2_TARGET, linetype = "dashed", colour = "grey45") +
  ggplot2::geom_vline(xintercept = softPower, linetype = "dashed", colour = "#C0392B") +
  ggplot2::annotate("text", x = softPower + 0.5, y = min(sft_df$SFT_R2),
                    label = paste0("β = ", softPower), colour = "#C0392B",
                    hjust = 0, fontface = "bold", size = 3.5) +
  ggplot2::labs(title = "Scale-free topology fit (LC-MS)",
                subtitle = paste0("Target R² ≥ ", SFT_R2_TARGET, " | ", NETWORK_TYPE, " ", COR_TYPE),
                x = "Soft-threshold power (β)", y = "Scale-free fit (R²)") +
  theme_pub()

save_plot(p_sft, file.path(dir_net, "SoftThreshold_choice"), w = 7, h = 5)

##############################
## PART 3: Network construction
##############################
cat("── NETWORK CONSTRUCTION ──\n")
cat("  Single block: nFeatures =", ncol(datExpr), "\n")

net <- WGCNA::blockwiseModules(
  datExpr,
  power          = softPower,
  networkType    = "signed",
  TOMType        = "signed",
  corType        = COR_TYPE,
  maxBlockSize   = ncol(datExpr) + 100,
  deepSplit      = DEEP_SPLIT,
  minModuleSize  = MIN_MOD_SIZE,
  mergeCutHeight = MERGE_CUT_H,
  numericLabels  = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs       = FALSE,
  verbose        = 5
)

colors <- WGCNA::labels2colors(net$colors)
MEs <- WGCNA::orderMEs(net$MEs)

n_modules <- length(unique(colors[colors != "grey"]))
log_msg("Modules: ", n_modules, " (+ grey)")
cat("Top module sizes:\n")
print(head(sort(table(colors), decreasing = TRUE), 15))

##############################
## PART 4: Save outputs
##############################
saveRDS(net,       file.path(dir_net, "net.rds"))
saveRDS(MEs,       file.path(dir_net, "MEs.rds"))
saveRDS(softPower, file.path(dir_net, "softPower.rds"))
saveRDS(datExpr,   file.path(dir_net, "datExpr_clean.rds"))
saveRDS(meta,      file.path(dir_net, "meta_clean.rds"))

# Feature → Module table
feature_mod <- data.frame(
  Feature     = colnames(datExpr),
  ModuleLabel = as.integer(net$colors),
  ModuleColor = colors,
  stringsAsFactors = FALSE
)
write.csv(feature_mod, file.path(dir_net, "Feature_to_Module.csv"), row.names = FALSE)

# ME → Color map
me_map <- data.frame(
  ME          = colnames(MEs),
  ModuleLabel = suppressWarnings(as.integer(gsub("^ME", "", colnames(MEs)))),
  ModuleColor = WGCNA::labels2colors(suppressWarnings(as.integer(gsub("^ME", "", colnames(MEs)))))
)
write.csv(me_map, file.path(dir_net, "ME_color_map.csv"), row.names = FALSE)

##############################
## PART 5: Figures
##############################
cat("── FIGURES ──\n")

# Feature dendrogram + modules
if (!is.null(net$dendrograms) && length(net$dendrograms) > 0) {
  save_base_plot(function(){
    WGCNA::plotDendroAndColors(
      net$dendrograms[[1]],
      WGCNA::labels2colors(net$colors[net$blockGenes[[1]]]),
      "Modules", dendroLabels = FALSE, hang = 0.03,
      main = paste0("Metabolite dendrogram — ", GENO, " | β=", softPower, " (LC-MS)")
    )
  }, file.path(dir_net, "Feature_dendrogram_modules"), w = 14, h = 6)
}

# Module sizes barplot
mod_sizes <- as.data.frame(table(colors)) |>
  setNames(c("Module", "N")) |>
  dplyr::filter(Module != "grey") |>
  dplyr::arrange(dplyr::desc(N)) |>
  dplyr::mutate(Module = factor(Module, levels = Module[order(N)]))

p_sizes <- ggplot2::ggplot(mod_sizes, ggplot2::aes(Module, N, fill = as.character(Module))) +
  ggplot2::geom_col(colour = "white", linewidth = 0.3) +
  ggplot2::geom_text(ggplot2::aes(label = N), hjust = -0.15, size = 3, fontface = "bold") +
  ggplot2::scale_fill_identity() +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.15))) +
  ggplot2::coord_flip() +
  ggplot2::labs(title = paste0("Module sizes — ", GENO, " (LC-MS)"),
                subtitle = paste0(n_modules, " modules | β=", softPower, " | ", NETWORK_TYPE),
                x = NULL, y = "Number of metabolites") +
  theme_pub() + ggplot2::theme(legend.position = "none")

save_plot(p_sizes, file.path(dir_net, "ModuleSizes_barplot"), w = 7,
          h = max(4, nrow(mod_sizes) * 0.35 + 2))

log_msg("Script 03 DONE")
cat("\n✓ Script 03 completed\n")