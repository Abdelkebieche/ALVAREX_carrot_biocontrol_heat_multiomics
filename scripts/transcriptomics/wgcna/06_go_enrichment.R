############################################################
## 06_enrichment_GO.R
## GO ORA + Pathway ORA (PlanT2T) par module significatif
## clusterProfiler + PlanT2T OrgDb
# Repository configuration (portable; no HPC-specific paths)
if (!exists('CARROT_CONFIG_LOADED')) source('../../../config/config.R')


############################################################

##############################
## PARAMÈTRES
##############################
INSTALL_MISSING_PKGS <- FALSE
FDR_CUTOFF     <- 0.05
KME_CUTOFF_ORA <- 0.40
GO_ONT_LIST    <- c("BP","MF","CC")
GO_SHOWCATEGORY <- 20
MIN_GENES_ORA  <- 10

ORGDB_URL <- "https://biobigdata.nju.edu.cn/plant2t/orgdb/org.Daucus.carota.DH13M14.eg.db_1.0.tar.gz"
ORGDB_PKG <- "org.Daucus.carota.DH13M14.eg.db"

map_file <- GENE_ID_MAP_FILE
base_out <- RNA_WGCNA_ROOT

##############################
## SOURCE + LOAD
##############################
source("00_helpers.R")
load_pkgs(c("dplyr","ggplot2"), install_missing = INSTALL_MISSING_PKGS)
load_pkgs(c("WGCNA","AnnotationDbi","clusterProfiler","enrichplot"),
          install_missing = INSTALL_MISSING_PKGS, bioc = TRUE)

GENOTYPE_CHOSEN <- Sys.getenv("CARROT_GENOTYPE", unset = "PRESTO")
GENO_TAG <- gsub("[^A-Za-z0-9]+", "", GENOTYPE_CHOSEN)
run_dir  <- find_latest_dir(base_out, paste0("^WGCNA_", GENO_TAG, "_"))
datExpr  <- readRDS(file.path(run_dir, "02_network", "datExpr_clean.rds"))
net      <- readRDS(file.path(run_dir, "02_network", "net.rds"))
MEs      <- readRDS(file.path(run_dir, "02_network", "MEs.rds"))
gene_mod <- read.csv(file.path(run_dir, "02_network", "Gene_to_Module.csv"))
anova_res <- read.csv(file.path(run_dir, "03_module_trait", "Module_ANOVA_results.csv"))

GENO <- as.character(readRDS(file.path(run_dir, "02_network", "meta_clean.rds"))$Genotype[1])

dir_enrich <- file.path(run_dir, "05_enrichment")
dir.create(dir_enrich, recursive = TRUE, showWarnings = FALSE)

log_msg <- init_logger(file.path(run_dir, "pipeline.log"))
log_msg("START Script 06 — Enrichment")
cat("── Script 06: Enrichment —", GENO, "──\n\n")

##############################
## PARTIE 1: Load OrgDb + mapping
##############################
cat("── Load OrgDb + mapping ──\n")

if (!requireNamespace(ORGDB_PKG, quietly = TRUE)) {
  if (!INSTALL_MISSING_PKGS) stop("OrgDb manquant: ", ORGDB_PKG)
  install.packages(ORGDB_URL, repos = NULL, type = "source")
}
suppressPackageStartupMessages(library(ORGDB_PKG, character.only = TRUE))

# Trouver l'objet OrgDb
OrgDb_obj <- NULL
ns <- asNamespace(ORGDB_PKG)
for (o in ls(ns, all.names = TRUE)) {
  x <- try(get(o, envir = ns), silent = TRUE)
  if (!inherits(x, "try-error") && inherits(x, "OrgDb")) { OrgDb_obj <- x; break }
}
if (is.null(OrgDb_obj)) stop("OrgDb object not found in: ", ORGDB_PKG)

stopifnot(file.exists(map_file))
map_df <- read.csv(map_file) |>
  dplyr::select(my_gene_id, plant2t_id) |>
  dplyr::mutate(my_gene_id = trimws(as.character(my_gene_id)),
                plant2t_id = trimws(as.character(plant2t_id))) |>
  dplyr::filter(!is.na(my_gene_id), my_gene_id != "") |>
  dplyr::distinct()

# Universe
universe_ids <- data.frame(my_gene_id = colnames(datExpr)) |>
  dplyr::left_join(map_df, by = "my_gene_id") |>
  dplyr::pull(plant2t_id) |> unique() |> na.omit()

cat("  Universe:", length(universe_ids), "PlantT2T IDs\n")

##############################
## PARTIE 2: Pathway mapping
##############################
cat("── Pathway mapping ──\n")

p2g <- AnnotationDbi::select(OrgDb_obj, keys = universe_ids,
                             columns = c("Pathway","GID"), keytype = "GID") |>
  dplyr::filter(!is.na(Pathway), Pathway != "", !is.na(GID)) |>
  dplyr::distinct(Pathway, GID)

ko_file <- system.file("extdata", "ko00001.PlanT2T.txt", package = ORGDB_PKG)
if (file.exists(ko_file)) {
  ko_raw <- utils::read.delim(ko_file, header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)
  id_col <- intersect(c("pathway_id","Pathway"), colnames(ko_raw))[1]
  nm_col <- intersect(c("level3","Name","name"), colnames(ko_raw))[1]
  p2n <- ko_raw |>
    dplyr::transmute(Pathway = trimws(as.character(.data[[id_col]])),
                     Name = trimws(as.character(.data[[nm_col]]))) |>
    dplyr::filter(Pathway != "", Name != "") |> dplyr::distinct()
} else {
  p2n <- data.frame(Pathway = unique(p2g$Pathway), Name = unique(p2g$Pathway))
}

cat("  Pathway pairs:", nrow(p2g), "\n")

##############################
## PARTIE 3: Sélection modules à enrichir
##############################
# Modules significatifs (au moins 1 facteur ANOVA)
sig_mods <- anova_res |>
  dplyr::filter(FDR_Temperature < 0.20 | FDR_Treatment < 0.20 |
                  FDR_Time < 0.20 | FDR_TempxTrt < 0.20)

if (nrow(sig_mods) < 3) {
  sig_mods <- anova_res |> dplyr::slice_head(n = min(10, nrow(anova_res)))
}

cat("  Modules to enrich:", nrow(sig_mods), "\n\n")

##############################
## PARTIE 4: GO ORA + Pathway ORA par module
##############################
all_go <- data.frame()
all_path <- data.frame()

for (idx in seq_len(nrow(sig_mods))) {
  me   <- sig_mods$ME[idx]
  colr <- sig_mods$ModuleColor[idx]
  lab  <- as.integer(gsub("^ME", "", me))
  
  cat("── Module:", colr, "(", me, ") ──\n")
  
  # Gènes du module avec kME >= cutoff
  mod_genes <- gene_mod$Gene[gene_mod$ModuleLabel == lab]
  
  # kME pour filtrer les gènes core
  kME_mat_mod <- WGCNA::bicor(datExpr[, mod_genes, drop = FALSE],
                              as.numeric(MEs[, me]),
                              use = "pairwise.complete.obs")
  kme_vec <- as.numeric(kME_mat_mod[, 1])
  names(kme_vec) <- mod_genes
  core_genes <- names(kme_vec[kme_vec >= KME_CUTOFF_ORA])
  
  cat("  All genes:", length(mod_genes), "| Core (kME>=", KME_CUTOFF_ORA, "):", length(core_genes), "\n")
  
  # Map to PlanT2T
  core_ids <- data.frame(my_gene_id = core_genes) |>
    dplyr::left_join(map_df, by = "my_gene_id") |>
    dplyr::pull(plant2t_id) |> unique() |> na.omit()
  
  if (length(core_ids) < MIN_GENES_ORA) {
    cat("  SKIP: too few mapped genes (", length(core_ids), ")\n")
    next
  }
  
  dir_mod <- file.path(dir_enrich, paste0(colr, "_", me))
  dir.create(dir_mod, recursive = TRUE, showWarnings = FALSE)
  
  # ── GO ORA ──
  for (ont in GO_ONT_LIST) {
    ego <- tryCatch(
      suppressMessages(clusterProfiler::enrichGO(
        gene = core_ids, universe = universe_ids,
        OrgDb = OrgDb_obj, keyType = "GID",
        ont = ont, pAdjustMethod = "BH",
        pvalueCutoff = 0.05, qvalueCutoff = 0.05
      )),
      error = function(e) { cat("  GO", ont, "error:", conditionMessage(e), "\n"); NULL }
    )
    
    if (is.null(ego)) next
    ego_df <- as.data.frame(ego)
    
    if (nrow(ego_df) > 0) {
      ego_df$Module <- colr; ego_df$Ontology <- ont
      write.csv(ego_df, file.path(dir_mod, paste0("GO_ORA_", ont, ".csv")), row.names = FALSE)
      all_go <- dplyr::bind_rows(all_go, ego_df)
      
      p <- enrichplot::dotplot(ego, showCategory = GO_SHOWCATEGORY) +
        ggplot2::labs(title = paste0(colr, " — GO ", ont, " (", GENO, ")")) +
        theme_pub(base_size = 9)
      save_plot(p, file.path(dir_mod, paste0("Dotplot_GO_", ont)), w = 8,
                h = max(4, min(nrow(ego_df), GO_SHOWCATEGORY) * 0.25 + 2))
    }
  }
  
  # ── Pathway ORA ──
  enr <- tryCatch(
    suppressMessages(clusterProfiler::enricher(
      gene = core_ids, universe = universe_ids,
      TERM2GENE = p2g, TERM2NAME = p2n,
      pAdjustMethod = "BH", pvalueCutoff = 0.05,
      qvalueCutoff = 0.05, minGSSize = 5
    )),
    error = function(e) { cat("  Pathway error:", conditionMessage(e), "\n"); NULL }
  )
  
  if (!is.null(enr)) {
    enr_df <- as.data.frame(enr)
    if (nrow(enr_df) > 0) {
      enr_df$Module <- colr
      write.csv(enr_df, file.path(dir_mod, "Pathway_ORA.csv"), row.names = FALSE)
      all_path <- dplyr::bind_rows(all_path, enr_df)
      
      p_pw <- enrichplot::dotplot(enr, showCategory = GO_SHOWCATEGORY) +
        ggplot2::labs(title = paste0(colr, " — Pathway (", GENO, ")")) +
        theme_pub(base_size = 9)
      save_plot(p_pw, file.path(dir_mod, "Dotplot_Pathway"), w = 8,
                h = max(4, min(nrow(enr_df), GO_SHOWCATEGORY) * 0.25 + 2))
    }
  }
  
  cat("  GO terms:", nrow(all_go |> dplyr::filter(Module == colr)),
      "| Pathways:", nrow(all_path |> dplyr::filter(Module == colr)), "\n")
}

# Save combined tables
if (nrow(all_go) > 0) write.csv(all_go, file.path(dir_enrich, "GO_ALL_modules_combined.csv"), row.names = FALSE)
if (nrow(all_path) > 0) write.csv(all_path, file.path(dir_enrich, "Pathway_ALL_modules_combined.csv"), row.names = FALSE)

log_msg("Script 06 DONE")
cat("\n✓ Script 06 terminé\n")