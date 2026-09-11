############################################################
## 07_module_profiles_classification.R
## Profils temporels ME par condition (P2→P3)
## Classification finale + résumé + Cytoscape export
# Repository configuration (portable; no HPC-specific paths)
if (!exists('CARROT_CONFIG_LOADED')) source('../../../config/config.R')


############################################################

##############################
## PARAMÈTRES
##############################
INSTALL_MISSING_PKGS <- FALSE
FDR_CUTOFF <- 0.05
TOP_CYTO_MODULES <- 4
CYTO_TOP_GENES <- 80

TEMP_COLORS <- c("T0" = "#313695", "T1" = "#74ADD1", "T2" = "#F46D43")
TRT_COLORS  <- c("EAU" = "#4D4D4D", "SDP" = "#CC79A7")
TRT_LINETYPES <- c("Control" = "dashed", "Treated" = "solid")

base_out <- RNA_WGCNA_ROOT

##############################
## SOURCE + LOAD
##############################
source("00_helpers.R")
load_pkgs(c("svglite","dplyr","tidyr","ggplot2","patchwork"), install_missing = INSTALL_MISSING_PKGS)
load_pkgs(c("WGCNA"), install_missing = INSTALL_MISSING_PKGS, bioc = TRUE)

GENOTYPE_CHOSEN <- Sys.getenv("CARROT_GENOTYPE", unset = "PRESTO")
GENO_TAG <- gsub("[^A-Za-z0-9]+", "", GENOTYPE_CHOSEN)
run_dir  <- find_latest_dir(base_out, paste0("^WGCNA_", GENO_TAG, "_"))
MEs      <- readRDS(file.path(run_dir, "02_network", "MEs.rds"))
meta     <- readRDS(file.path(run_dir, "02_network", "meta_clean.rds"))
datExpr  <- readRDS(file.path(run_dir, "02_network", "datExpr_clean.rds"))
net      <- readRDS(file.path(run_dir, "02_network", "net.rds"))
me_map   <- read.csv(file.path(run_dir, "02_network", "ME_color_map.csv"))
gene_mod <- read.csv(file.path(run_dir, "02_network", "Gene_to_Module.csv"))
classify <- read.csv(file.path(run_dir, "03_module_trait", "Module_Classification.csv"))
anova_res <- read.csv(file.path(run_dir, "03_module_trait", "Module_ANOVA_results.csv"))

# Hub genes si disponibles
hub_file <- file.path(run_dir, "04_hub_genes", "AllGenes_kME_GS.csv")
hub_df <- if (file.exists(hub_file)) read.csv(hub_file) else NULL

# GO enrichment si disponible
go_file <- file.path(run_dir, "05_enrichment", "GO_ALL_modules_combined.csv")
go_all <- if (file.exists(go_file)) read.csv(go_file) else NULL

GENO <- as.character(meta$Genotype[1])
me_to_color <- setNames(me_map$ModuleColor, me_map$ME)

dir_prof <- file.path(run_dir, "06_profiles")
dir_cyto <- file.path(run_dir, "06_profiles", "cytoscape")
for (d in c(dir_prof, dir_cyto)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

log_msg <- init_logger(file.path(run_dir, "pipeline.log"))
log_msg("START Script 07 — Profiles + Final summary")
cat("── Script 07: Profiles + Summary —", GENO, "──\n\n")


##############################################################
## PARTIE 1: Profils ME par condition (P2 → P3)
##############################################################
cat("── PARTIE 1: Profils temporels ──\n")

MEs_clean <- MEs[, !grepl("^ME0$", colnames(MEs)), drop = FALSE]

# Long format: ME × condition
me_long <- as.data.frame(MEs_clean) |>
  dplyr::mutate(
    Sample        = meta$Sample,
    Temperature   = meta$Temperature,
    Treatment     = meta$Treatment,
    Sampling_time = meta$Sampling_time
  ) |>
  tidyr::pivot_longer(
    cols = starts_with("ME"),
    names_to = "ME",
    values_to = "Eigengene"
  ) |>
  dplyr::mutate(
    ModuleColor = me_to_color[ME],
    Group = paste(Temperature, Treatment, sep = "_")
  )

# Moyennes par groupe
me_means <- me_long |>
  dplyr::group_by(ME, ModuleColor, Temperature, Treatment, Sampling_time, Group) |>
  dplyr::summarise(
    mean_ME = mean(Eigengene, na.rm = TRUE),
    se_ME = sd(Eigengene, na.rm = TRUE) / sqrt(dplyr::n()),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    Time_num = ifelse(Sampling_time == "P2", 1, 2)
  )

# Plot profils pour chaque module significatif
sig_mods <- classify |> dplyr::filter(Classification != "Stable")
if (nrow(sig_mods) == 0) sig_mods <- classify |> dplyr::slice_head(n = 8)

profile_plots <- list()

for (i in seq_len(nrow(sig_mods))) {
  me <- sig_mods$ME[i]
  colr <- sig_mods$ModuleColor[i]
  class_label <- sig_mods$Classification[i]
  
  sub <- me_means |> dplyr::filter(ME == me)
  if (nrow(sub) == 0) next
  
  p <- ggplot2::ggplot(sub, ggplot2::aes(x = Time_num, y = mean_ME,
                                         colour = Temperature, linetype = Treatment,
                                         group = Group)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = mean_ME - se_ME, ymax = mean_ME + se_ME),
                           width = 0.08, linewidth = 0.5) +
    ggplot2::scale_colour_manual(values = TEMP_COLORS, labels = c("T0"="NH","T1"="HS-7","T2"="HS-2")) +
    ggplot2::scale_linetype_manual(values = TRT_LINETYPES) +
    ggplot2::scale_x_continuous(breaks = c(1, 2), labels = c("D2", "D4")) +
    ggplot2::labs(title = paste0(colr, " — ", class_label),
                  x = NULL, y = "ME") +
    theme_pub(base_size = 9) +
    ggplot2::theme(legend.position = "bottom",
                   legend.key.width = ggplot2::unit(15, "mm"))
  
  profile_plots[[colr]] <- p
}

# Combiner en une grande figure
if (length(profile_plots) > 0) {
  n_cols <- min(3, length(profile_plots))
  p_all <- patchwork::wrap_plots(profile_plots, ncol = n_cols) +
    patchwork::plot_annotation(
      title = paste0("Module eigengene profiles — ", GENO, " (P2→P3)"),
      subtitle = "Solid = Treated | Dashed = Control | Color = Temperature",
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 14),
        plot.subtitle = ggplot2::element_text(colour = "grey40", size = 11)
      )
    )
  
  n_rows <- ceiling(length(profile_plots) / n_cols)
  save_plot(p_all, file.path(dir_prof, "ModuleProfiles_allModules"),
            w = n_cols * 4.5, h = n_rows * 3.5 + 1)
  
  # Aussi individuellement
  for (colr in names(profile_plots)) {
    save_plot(profile_plots[[colr]],
              file.path(dir_prof, paste0("Profile_", colr)),
              w = 5, h = 4)
  }
}


##############################################################
## PARTIE 2: Tableau final résumé (1 ligne par module)
##############################################################
cat("── PARTIE 2: Final summary table ──\n")

final_summary <- classify |>
  dplyr::left_join(
    gene_mod |> dplyr::filter(ModuleColor != "grey") |>
      dplyr::count(ModuleColor, name = "ModuleSize"),
    by = "ModuleColor"
  )

# Ajouter top hub gene
if (!is.null(hub_df)) {
  top_hub <- hub_df |>
    dplyr::filter(ModuleColor != "grey") |>
    dplyr::group_by(ModuleColor) |>
    dplyr::arrange(dplyr::desc(kME), .by_group = TRUE) |>
    dplyr::slice_head(n = 1) |>
    dplyr::ungroup() |>
    dplyr::select(ModuleColor, HubGene = Gene, Hub_kME = kME)
  
  final_summary <- final_summary |> dplyr::left_join(top_hub, by = "ModuleColor")
}

# Ajouter top GO BP term
if (!is.null(go_all)) {
  desc_col <- intersect(c("Description","Term"), colnames(go_all))[1]
  padj_col <- intersect(c("p.adjust","qvalue"), colnames(go_all))[1]
  ont_col <- intersect(c("Ontology","ONTOLOGY","ont"), colnames(go_all))[1]
  
  if (!is.na(desc_col) && !is.na(padj_col)) {
    top_go <- go_all
    if (!is.na(ont_col)) top_go <- top_go |> dplyr::filter(.data[[ont_col]] == "BP")
    
    top_go <- top_go |>
      dplyr::group_by(Module) |>
      dplyr::arrange(.data[[padj_col]], .by_group = TRUE) |>
      dplyr::slice_head(n = 3) |>
      dplyr::summarise(
        Top_GO_BP = paste(.data[[desc_col]], collapse = " | "),
        .groups = "drop"
      ) |>
      dplyr::rename(ModuleColor = Module)
    
    final_summary <- final_summary |> dplyr::left_join(top_go, by = "ModuleColor")
  }
}

final_summary <- final_summary |> dplyr::arrange(Classification, dplyr::desc(ModuleSize))
write.csv(final_summary, file.path(dir_prof, "Module_Final_Summary.csv"), row.names = FALSE)

cat("  Final summary:", nrow(final_summary), "modules\n")
cat("  Classification breakdown:\n")
print(table(final_summary$Classification))


##############################################################
## PARTIE 3: Cytoscape export (top modules)
##############################################################
cat("\n── PARTIE 3: Cytoscape export ──\n")

# Sélectionner les modules les plus intéressants
# ★ FIX: vérifier que R2 existe dans anova_res, sinon utiliser p-values
top_mods <- classify |>
  dplyr::filter(Classification != "Stable")

if (nrow(top_mods) == 0) {
  # Fallback: prendre les premiers modules par taille
  top_mods <- classify |> dplyr::slice_head(n = TOP_CYTO_MODULES)
} else {
  # Joindre anova_res pour trier par importance
  cols_to_join <- intersect(c("ME", "R2", "FDR_Temperature", "FDR_Treatment"), colnames(anova_res))
  if (length(cols_to_join) > 1) {
    top_mods <- top_mods |>
      dplyr::left_join(anova_res |> dplyr::select(dplyr::all_of(cols_to_join)), by = "ME")
  }
  
  # Trier: R2 si disponible, sinon par min(FDR)
  if ("R2" %in% colnames(top_mods)) {
    top_mods <- top_mods |> dplyr::arrange(dplyr::desc(R2))
  } else if ("FDR_Temperature" %in% colnames(top_mods)) {
    top_mods <- top_mods |>
      dplyr::mutate(.min_fdr = pmin(FDR_Temperature, FDR_Treatment, na.rm = TRUE)) |>
      dplyr::arrange(.min_fdr) |>
      dplyr::select(-.min_fdr)
  }
  
  top_mods <- top_mods |> dplyr::slice_head(n = TOP_CYTO_MODULES)
}

cat("  Modules for Cytoscape:", nrow(top_mods), "\n")

if (nrow(top_mods) == 0) {
  cat("  No modules for Cytoscape export\n")
} else {
  for (i in seq_len(nrow(top_mods))) {
    colr <- top_mods$ModuleColor[i]
    lab  <- as.integer(gsub("^ME", "", top_mods$ME[i]))
    
    mod_genes <- gene_mod$Gene[gene_mod$ModuleLabel == lab]
    if (length(mod_genes) < 3) { cat("  SKIP", colr, "(too few genes)\n"); next }
    
    # Prendre les top genes par kME si hub_df est disponible
    if (!is.null(hub_df) && "kME" %in% colnames(hub_df)) {
      top_genes <- hub_df |>
        dplyr::filter(ModuleColor == colr, !is.na(kME)) |>
        dplyr::arrange(dplyr::desc(kME)) |>
        dplyr::slice_head(n = CYTO_TOP_GENES) |>
        dplyr::pull(Gene)
    } else {
      top_genes <- mod_genes[seq_len(min(CYTO_TOP_GENES, length(mod_genes)))]
    }
    
    if (length(top_genes) < 3) { cat("  SKIP", colr, "(too few top genes)\n"); next }
    
    # Calculer les corrélations entre les top genes
    genes_ok <- intersect(top_genes, colnames(datExpr))
    if (length(genes_ok) < 3) next
    
    expr_sub <- datExpr[, genes_ok, drop = FALSE]
    cor_sub <- WGCNA::bicor(expr_sub, use = "pairwise.complete.obs")
    
    # Edge list (top correlations)
    edges <- which(upper.tri(cor_sub) & abs(cor_sub) >= 0.5, arr.ind = TRUE)
    if (nrow(edges) == 0) {
      # Relax threshold if no edges
      edges <- which(upper.tri(cor_sub) & abs(cor_sub) >= 0.3, arr.ind = TRUE)
      if (nrow(edges) == 0) { cat("  SKIP", colr, "(no edges)\n"); next }
    }
    
    edge_df <- data.frame(
      Source = rownames(cor_sub)[edges[, 1]],
      Target = colnames(cor_sub)[edges[, 2]],
      Weight = cor_sub[edges],
      stringsAsFactors = FALSE
    ) |> dplyr::arrange(dplyr::desc(abs(Weight)))
    
    # Node table
    nodes_in <- unique(c(edge_df$Source, edge_df$Target))
    node_df <- data.frame(Gene = nodes_in, stringsAsFactors = FALSE) |>
      dplyr::left_join(gene_mod, by = "Gene")
    
    if (!is.null(hub_df)) {
      hub_cols <- intersect(c("Gene","kME","GS_Temperature","GS_Treatment"), colnames(hub_df))
      if (length(hub_cols) > 1) {
        node_df <- node_df |>
          dplyr::left_join(hub_df |> dplyr::select(dplyr::all_of(hub_cols)), by = "Gene")
      }
    }
    
    write.csv(edge_df, file.path(dir_cyto, paste0(colr, "_edges.csv")), row.names = FALSE)
    write.csv(node_df, file.path(dir_cyto, paste0(colr, "_nodes.csv")), row.names = FALSE)
    
    cat("  Cytoscape:", colr, "—", nrow(node_df), "nodes,", nrow(edge_df), "edges\n")
  }
}


##############################################################
## PARTIE 4: Session info
##############################################################
sink(file.path(run_dir, "sessionInfo.txt"))
cat("WGCNA v2 Pipeline —", GENO, "(P2+P3)\n")
cat("Date:", format(Sys.time()), "\n\n")
print(sessionInfo())
sink()

log_msg("Script 07 DONE — Pipeline complete")
cat("\n══════════════════════════════════════\n")
cat("  WGCNA v2 PIPELINE COMPLETE\n")
cat("  Génotype:", GENO, "\n")
cat("  Output:", run_dir, "\n")
cat("══════════════════════════════════════\n\n")

cat("Dossiers:\n")
cat("  00_inputs/        — metadata, paramètres\n")
cat("  01_QC/            — PCA, dendro, heatmap\n")
cat("  02_network/       — net, MEs, soft threshold\n")
cat("  03_module_trait/  — heatmaps, ANOVA, boxplots, classification\n")
cat("  04_hub_genes/     — kME, scatter, by_module\n")
cat("  05_enrichment/    — GO, Pathway par module\n")
cat("  06_profiles/      — profils, résumé final, Cytoscape\n")