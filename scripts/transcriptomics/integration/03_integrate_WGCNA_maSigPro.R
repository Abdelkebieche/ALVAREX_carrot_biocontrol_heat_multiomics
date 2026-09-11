
# Repository configuration (portable; no HPC-specific paths)
if (!exists('CARROT_CONFIG_LOADED')) source('../../../config/config.R')

# ======================================================================
# WGCNA × maSigPro INTEGRATION — PRESTO + ROBILA
# Publication-oriented integration of co-expression modules and temporal
# expression clusters.
#
# Main questions
#   1) Which maSigPro temporal clusters are preferentially embedded in
#      which WGCNA modules?
#   2) What fraction of each WGCNA module is dynamically regulated?
#   3) Which biological functions are supported independently by both
#      WGCNA module enrichment and maSigPro cluster enrichment?
#   4) Which WGCNA hub genes also belong to significant maSigPro profiles?
#   5) Which dynamic genes / cluster relationships are conserved between
#      PRESTO and ROBILA?
#
# IMPORTANT BIOLOGICAL INTERPRETATION
#   - WGCNA was built from D2/D4 samples and describes co-expression structure.
#   - maSigPro uses D0/D2/D4 and describes expression trajectories.
#   - Therefore this script DOES NOT call D2→D4 changes "time effects" by itself.
#     It integrates network membership with full maSigPro trajectories.
#
# DEFAULT
#   STRICT_ONLY = TRUE
#   maSigPro clusters are retained, but only genes passing the strict T.fit
#   selection are used for inferential overlap analyses.
#   Set STRICT_ONLY = FALSE to reproduce the broader p.vector-selected set.
# ======================================================================

# ----------------------------------------------------------------------
# 0. PACKAGES
# ----------------------------------------------------------------------

required_pkgs <- c(
  "dplyr", "tidyr", "readr", "stringr", "ggplot2", "forcats", "scales"
)

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0) {
  stop(
    "Missing R packages: ", paste(missing_pkgs, collapse = ", "),
    "\nInstall them before running this script."
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(forcats)
  library(scales)
})

# Explicit namespace protection against count()/select() conflicts
count     <- dplyr::count
select    <- dplyr::select
filter    <- dplyr::filter
mutate    <- dplyr::mutate
summarise <- dplyr::summarise
group_by  <- dplyr::group_by
ungroup   <- dplyr::ungroup

# ----------------------------------------------------------------------
# 1. USER SETTINGS
# ----------------------------------------------------------------------

FDR_CUTOFF <- 0.05
STRICT_ONLY <- TRUE
EXCLUDE_GREY_FROM_MAIN_PLOTS <- TRUE

# Minimum overlap for a cluster × module pair to be highlighted as a
# biologically interpretable association.
MIN_OVERLAP_GENES <- 10

# Hub definitions
KME_HUB_CUTOFF <- 0.85
GS_TREATMENT_CUTOFF <- 0.30

# Functional enrichment cut-off
FUNCTION_FDR <- 0.05

# Project root
PROJECT_ROOT <- OUTPUT_ROOT

# Root directory containing the WGCNA runs.
WGCNA_ROOT <- RNA_WGCNA_ROOT

# Preferred root directory containing PRESTO_maSigPro and ROBILA_maSigPro.
# If this exact folder does not exist, the script also searches PROJECT_ROOT
# recursively for directories named PRESTO_maSigPro / ROBILA_maSigPro.
MASIGPRO_ROOT <- MASIGPRO_ROOT

# Optional manual run paths.
# Leave NA to auto-detect.
PRESTO_WGCNA <- NA_character_
ROBILA_WGCNA <- NA_character_
PRESTO_MASIG <- NA_character_
ROBILA_MASIG <- NA_character_

# Output directory
OUTPUT_ROOT <- file.path(
  WGCNA_ROOT,
  paste0(
    "WGCNA_maSigPro_Integration_",
    format(Sys.time(), "%Y%m%d_%H%M%S")
  )
)

DIR_TAB <- file.path(OUTPUT_ROOT, "01_tables")
DIR_FIG <- file.path(OUTPUT_ROOT, "02_figures")
DIR_RDS <- file.path(OUTPUT_ROOT, "03_rds")
DIR_REPORT <- file.path(OUTPUT_ROOT, "04_report")

for (d in c(DIR_TAB, DIR_FIG, DIR_RDS, DIR_REPORT)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ----------------------------------------------------------------------
# 2. HELPERS — PATH DISCOVERY
# ----------------------------------------------------------------------

latest_matching_dir <- function(root, pattern, recursive = TRUE) {
  if (!dir.exists(root)) return(NA_character_)
  dd <- list.dirs(root, recursive = recursive, full.names = TRUE)
  hit <- dd[grepl(pattern, basename(dd))]
  if (length(hit) == 0) return(NA_character_)
  hit[order(file.info(hit)$mtime, decreasing = TRUE)][1]
}

find_masig_dir <- function(root, genotype) {
  target <- paste0(genotype, "_maSigPro")
  
  # First try exact immediate directory
  exact <- file.path(root, target)
  if (dir.exists(exact)) return(exact)
  
  # Then search preferred root recursively
  z <- latest_matching_dir(root, paste0("^", target, "$"))
  if (!is.na(z)) return(z)
  
  # Final fallback: search whole project root recursively
  latest_matching_dir(PROJECT_ROOT, paste0("^", target, "$"))
}

if (is.na(PRESTO_WGCNA)) {
  PRESTO_WGCNA <- latest_matching_dir(WGCNA_ROOT, "^WGCNA_PRESTO_")
}
if (is.na(ROBILA_WGCNA)) {
  ROBILA_WGCNA <- latest_matching_dir(WGCNA_ROOT, "^WGCNA_ROBILA_")
}
if (is.na(PRESTO_MASIG)) {
  PRESTO_MASIG <- find_masig_dir(MASIGPRO_ROOT, "PRESTO")
}
if (is.na(ROBILA_MASIG)) {
  ROBILA_MASIG <- find_masig_dir(MASIGPRO_ROOT, "ROBILA")
}

paths_required <- c(
  PRESTO_WGCNA = PRESTO_WGCNA,
  ROBILA_WGCNA = ROBILA_WGCNA,
  PRESTO_MASIG = PRESTO_MASIG,
  ROBILA_MASIG = ROBILA_MASIG
)

for (nm in names(paths_required)) {
  if (is.na(paths_required[[nm]]) || !dir.exists(paths_required[[nm]])) {
    stop(
      nm, " was not found.\n",
      "Set the corresponding path manually in section 1."
    )
  }
}

message("PRESTO WGCNA : ", PRESTO_WGCNA)
message("ROBILA WGCNA : ", ROBILA_WGCNA)
message("PRESTO maSig : ", PRESTO_MASIG)
message("ROBILA maSig : ", ROBILA_MASIG)
message("OUTPUT        : ", OUTPUT_ROOT)

# ----------------------------------------------------------------------
# 3. GENERIC HELPERS
# ----------------------------------------------------------------------

safe_read_csv <- function(path) {
  if (!file.exists(path)) {
    stop("Required file not found:\n", path)
  }
  readr::read_csv(path, show_col_types = FALSE)
}

safe_read_optional <- function(path) {
  if (!file.exists(path)) return(tibble::tibble())
  readr::read_csv(path, show_col_types = FALSE)
}

bh <- function(p) p.adjust(p, method = "BH")

save_plot <- function(p, basename, w = 9, h = 6) {
  ggplot2::ggsave(
    file.path(DIR_FIG, paste0(basename, ".png")),
    p, width = w, height = h, dpi = 400
  )
  ggplot2::ggsave(
    file.path(DIR_FIG, paste0(basename, ".pdf")),
    p, width = w, height = h
  )
  if (requireNamespace("svglite", quietly = TRUE)) {
    ggplot2::ggsave(
      file.path(DIR_FIG, paste0(basename, ".svg")),
      p, width = w, height = h,
      device = svglite::svglite
    )
  }
}

theme_pub <- function(base_size = 10) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey94", colour = NA),
      strip.text = element_text(face = "bold"),
      legend.title = element_text(face = "bold")
    )
}

# ----------------------------------------------------------------------
# 4. BIOLOGICAL THEME CLASSIFIER
# ----------------------------------------------------------------------

classify_theme <- function(description) {
  x <- stringr::str_to_lower(description)
  
  dplyr::case_when(
    str_detect(x, "defen[cs]e|immune|biotic stimulus|pathogen|systemic acquired") ~
      "Defence & immunity",
    str_detect(x, "hydrogen peroxide|reactive oxygen|oxidative stress|redox") ~
      "ROS & redox",
    str_detect(x, "mapk|signal transduction|signaling|signalling|receptor|cell recognition|cell communication") ~
      "Signalling & recognition",
    str_detect(x, "abscisic|auxin|jasmon|salicyl|ethylene|hormone") ~
      "Hormone signalling",
    str_detect(x, "phenylpropanoid|flavonoid|lignin|secondary metabol") ~
      "Phenylpropanoid & secondary metabolism",
    str_detect(x, "lipid|fatty acid|glycerolipid|phospholipid|wax|steroid") ~
      "Lipid metabolism",
    str_detect(x, "photosynth|light harvesting|light reaction") ~
      "Photosynthesis",
    str_detect(x, "translation|ribosom|peptide biosynt|amide biosynt") ~
      "Translation & ribosome",
    str_detect(x, "chaperone|protein folding|endoplasmic reticulum|protein processing") ~
      "Proteostasis & folding",
    str_detect(x, "proteasom|ubiquitin|deubiquit") ~
      "Protein turnover",
    str_detect(x, "rna processing|mrna|trna|spliceos|rna degrad|rna biosynt") ~
      "RNA processing",
    str_detect(x, "chromatin|chromosome|dna replication|dna repair|recombination") ~
      "Chromatin & DNA",
    str_detect(x, "transcription|gene expression") ~
      "Transcriptional regulation",
    str_detect(x, "glycol|gluconeogenesis|pentose phosphate|carbohydrate|starch|sucrose|calvin|carbon fixation") ~
      "Carbon & carbohydrate metabolism",
    str_detect(x, "organic acid|oxoacid|carboxylic acid|amino acid|pyruvate") ~
      "Organic-acid & amino-acid metabolism",
    str_detect(x, "oxidative phosphorylation|electron transport|atp synth") ~
      "Energy metabolism",
    str_detect(x, "mitochond") ~
      "Mitochondrial functions",
    str_detect(x, "transport|trafficking|localization|localisation|snare|transmembrane|ion channel") ~
      "Transport & trafficking",
    str_detect(x, "cytochrome p450|detox|xenobiotic") ~
      "Detoxification",
    TRUE ~ "Other"
  )
}

# ----------------------------------------------------------------------
# 5. LOAD ONE GENOTYPE
# ----------------------------------------------------------------------

load_genotype <- function(genotype, wgcna_dir, masig_dir) {
  
  # ---- WGCNA membership
  module_df <- safe_read_csv(
    file.path(wgcna_dir, "02_network", "Gene_to_Module.csv")
  ) |>
    dplyr::transmute(
      Gene = as.character(Gene),
      ModuleLabel = ModuleLabel,
      ModuleColor = as.character(ModuleColor)
    )
  
  # ---- WGCNA connectivity / GS
  hub_df <- safe_read_csv(
    file.path(wgcna_dir, "04_hub_genes", "AllGenes_kME_GS.csv")
  ) |>
    dplyr::mutate(
      Gene = as.character(Gene),
      ModuleColor = as.character(ModuleColor),
      kME = as.numeric(kME),
      GS_Temperature = as.numeric(GS_Temperature),
      GS_Treatment = as.numeric(GS_Treatment),
      GS_Time = as.numeric(GS_Time)
    )
  
  # ---- maSigPro cluster assignment
  cluster_df <- safe_read_csv(
    file.path(masig_dir, "04_tables", "Gene_cluster_assignment.csv")
  ) |>
    dplyr::transmute(
      Gene = as.character(Gene),
      Cluster = as.integer(Cluster)
    )
  
  strict_df <- safe_read_csv(
    file.path(masig_dir, "04_tables", "maSigPro_strict_genes_Tfit.csv")
  ) |>
    dplyr::transmute(Gene = as.character(Gene))
  
  # ---- maSigPro analysis universe after MAD
  expr_after_mad <- readRDS(
    file.path(masig_dir, "00_inputs", "expr_vst_after_MAD.rds")
  )
  masig_universe <- rownames(expr_after_mad)
  
  if (is.null(masig_universe)) {
    stop(
      "rownames are missing from expr_vst_after_MAD.rds for ", genotype
    )
  }
  
  # ---- Strict set for integration
  if (STRICT_ONLY) {
    cluster_used <- cluster_df |>
      dplyr::semi_join(strict_df, by = "Gene")
  } else {
    cluster_used <- cluster_df
  }
  
  # ---- Common eligible universe
  common_universe <- intersect(module_df$Gene, masig_universe)
  
  # Restrict all objects to common analytical universe
  module_u <- module_df |>
    dplyr::filter(Gene %in% common_universe)
  
  hub_u <- hub_df |>
    dplyr::filter(Gene %in% common_universe)
  
  cluster_u <- cluster_used |>
    dplyr::filter(Gene %in% common_universe)
  
  # Join dynamic genes with WGCNA membership and hub statistics
  integrated_genes <- cluster_u |>
    dplyr::left_join(module_u, by = "Gene") |>
    dplyr::left_join(
      hub_u |>
        dplyr::select(
          Gene, kME, GS_Temperature, GS_Treatment, GS_Time,
          IsHub, IsHighGS_Temp, IsHighGS_Trt
        ),
      by = "Gene"
    ) |>
    dplyr::mutate(
      Genotype = genotype,
      DynamicHub = !is.na(kME) & kME >= KME_HUB_CUTOFF,
      DynamicTreatmentHub =
        !is.na(kME) &
        kME >= KME_HUB_CUTOFF &
        !is.na(GS_Treatment) &
        abs(GS_Treatment) >= GS_TREATMENT_CUTOFF
    )
  
  # ---- Functional enrichment
  w_go <- safe_read_optional(
    file.path(wgcna_dir, "05_enrichment", "GO_ALL_modules_combined.csv")
  )
  w_path <- safe_read_optional(
    file.path(wgcna_dir, "05_enrichment", "Pathway_ALL_modules_combined.csv")
  )
  m_go <- safe_read_optional(
    file.path(masig_dir, "05_enrichment", "GO_ALL_combined.csv")
  )
  m_path <- safe_read_optional(
    file.path(masig_dir, "05_enrichment", "Pathway_ALL_combined.csv")
  )
  
  if (nrow(w_go) > 0) {
    w_go <- w_go |>
      dplyr::mutate(
        Genotype = genotype,
        Analysis = "WGCNA",
        Source = "GO_BP",
        ModuleColor = as.character(Module),
        p.adjust = as.numeric(p.adjust)
      ) |>
      dplyr::filter(
        Ontology == "BP",
        !is.na(p.adjust),
        p.adjust < FUNCTION_FDR
      )
  }
  
  if (nrow(m_go) > 0) {
    m_go <- m_go |>
      dplyr::mutate(
        Genotype = genotype,
        Analysis = "maSigPro",
        Source = "GO_BP",
        Cluster = as.integer(Cluster),
        p.adjust = as.numeric(p.adjust)
      ) |>
      dplyr::filter(
        Ontology == "BP",
        !is.na(p.adjust),
        p.adjust < FUNCTION_FDR
      )
  }
  
  if (nrow(w_path) > 0) {
    w_path <- w_path |>
      dplyr::mutate(
        Genotype = genotype,
        Analysis = "WGCNA",
        Source = "Pathway",
        ModuleColor = as.character(Module),
        p.adjust = as.numeric(p.adjust)
      ) |>
      dplyr::filter(
        !is.na(p.adjust),
        p.adjust < FUNCTION_FDR
      )
  }
  
  if (nrow(m_path) > 0) {
    m_path <- m_path |>
      dplyr::mutate(
        Genotype = genotype,
        Analysis = "maSigPro",
        Source = "Pathway",
        Cluster = as.integer(Cluster),
        p.adjust = as.numeric(p.adjust)
      ) |>
      dplyr::filter(
        !is.na(p.adjust),
        p.adjust < FUNCTION_FDR
      )
  }
  
  list(
    genotype = genotype,
    wgcna_dir = wgcna_dir,
    masig_dir = masig_dir,
    module = module_u,
    hub = hub_u,
    clusters = cluster_u,
    integrated_genes = integrated_genes,
    universe = common_universe,
    w_go = w_go,
    w_path = w_path,
    m_go = m_go,
    m_path = m_path
  )
}

P <- load_genotype("PRESTO", PRESTO_WGCNA, PRESTO_MASIG)
R <- load_genotype("ROBILA", ROBILA_WGCNA, ROBILA_MASIG)

# ----------------------------------------------------------------------
# 6. CLUSTER × MODULE OVERLAP AND HYPERGEOMETRIC ENRICHMENT
# ----------------------------------------------------------------------

cluster_module_overlap <- function(x) {
  
  N <- length(x$universe)
  
  cluster_sizes <- x$clusters |>
    dplyr::count(Cluster, name = "ClusterSize")
  
  module_sizes <- x$module |>
    dplyr::count(ModuleColor, name = "ModuleSize")
  
  observed <- x$integrated_genes |>
    dplyr::count(Cluster, ModuleColor, name = "Overlap")
  
  grid <- tidyr::crossing(
    Cluster = sort(unique(x$clusters$Cluster)),
    ModuleColor = sort(unique(x$module$ModuleColor))
  ) |>
    dplyr::left_join(cluster_sizes, by = "Cluster") |>
    dplyr::left_join(module_sizes, by = "ModuleColor") |>
    dplyr::left_join(observed, by = c("Cluster", "ModuleColor")) |>
    dplyr::mutate(
      Overlap = tidyr::replace_na(Overlap, 0L),
      Expected = ClusterSize * ModuleSize / N,
      FoldEnrichment = dplyr::if_else(
        Expected > 0,
        Overlap / Expected,
        NA_real_
      ),
      Fraction_of_cluster = Overlap / ClusterSize,
      Fraction_of_module_dynamic = Overlap / ModuleSize,
      Jaccard = dplyr::if_else(
        (ClusterSize + ModuleSize - Overlap) > 0,
        Overlap / (ClusterSize + ModuleSize - Overlap),
        NA_real_
      ),
      p_value = dplyr::if_else(
        Overlap == 0,
        1,
        stats::phyper(
          q = Overlap - 1,
          m = ModuleSize,
          n = N - ModuleSize,
          k = ClusterSize,
          lower.tail = FALSE
        )
      ),
      Genotype = x$genotype
    ) |>
    dplyr::mutate(
      FDR = bh(p_value),
      Significant =
        FDR < FDR_CUTOFF &
        FoldEnrichment > 1 &
        Overlap >= MIN_OVERLAP_GENES,
      neglogFDR = -log10(pmax(FDR, 1e-300)),
      log2FE = log2(pmax(FoldEnrichment, 1e-12))
    )
  
  grid
}

OV_P <- cluster_module_overlap(P)
OV_R <- cluster_module_overlap(R)
OV_ALL <- dplyr::bind_rows(OV_P, OV_R)

readr::write_csv(
  OV_ALL,
  file.path(DIR_TAB, "01_Cluster_x_WGCNA_Module_overlap_enrichment.csv")
)

# Best module for each maSigPro cluster
BEST_CLUSTER_MODULE <- OV_ALL |>
  dplyr::filter(ModuleColor != "grey") |>
  dplyr::group_by(Genotype, Cluster) |>
  dplyr::arrange(FDR, dplyr::desc(FoldEnrichment), dplyr::desc(Overlap)) |>
  dplyr::slice_head(n = 1) |>
  dplyr::ungroup()

readr::write_csv(
  BEST_CLUSTER_MODULE,
  file.path(DIR_TAB, "01_Best_WGCNA_module_per_maSigPro_cluster.csv")
)

# Best maSigPro cluster for each WGCNA module
BEST_MODULE_CLUSTER <- OV_ALL |>
  dplyr::filter(ModuleColor != "grey") |>
  dplyr::group_by(Genotype, ModuleColor) |>
  dplyr::arrange(FDR, dplyr::desc(FoldEnrichment), dplyr::desc(Overlap)) |>
  dplyr::slice_head(n = 1) |>
  dplyr::ungroup()

readr::write_csv(
  BEST_MODULE_CLUSTER,
  file.path(DIR_TAB, "01_Best_maSigPro_cluster_per_WGCNA_module.csv")
)

# ----------------------------------------------------------------------
# 7. DYNAMIC FRACTION OF EACH WGCNA MODULE
# ----------------------------------------------------------------------

module_dynamic_summary <- function(x) {
  
  x$module |>
    dplyr::count(ModuleColor, name = "WGCNA_module_size") |>
    dplyr::left_join(
      x$integrated_genes |>
        dplyr::count(ModuleColor, name = "Dynamic_genes"),
      by = "ModuleColor"
    ) |>
    dplyr::left_join(
      x$integrated_genes |>
        dplyr::filter(DynamicHub) |>
        dplyr::count(ModuleColor, name = "Dynamic_hubs"),
      by = "ModuleColor"
    ) |>
    dplyr::left_join(
      x$integrated_genes |>
        dplyr::filter(DynamicTreatmentHub) |>
        dplyr::count(ModuleColor, name = "Dynamic_treatment_hubs"),
      by = "ModuleColor"
    ) |>
    dplyr::mutate(
      Dynamic_genes = tidyr::replace_na(Dynamic_genes, 0L),
      Dynamic_hubs = tidyr::replace_na(Dynamic_hubs, 0L),
      Dynamic_treatment_hubs =
        tidyr::replace_na(Dynamic_treatment_hubs, 0L),
      Dynamic_fraction = Dynamic_genes / WGCNA_module_size,
      Dynamic_percent = 100 * Dynamic_fraction,
      Genotype = x$genotype
    )
}

MOD_DYN_P <- module_dynamic_summary(P)
MOD_DYN_R <- module_dynamic_summary(R)
MOD_DYN <- dplyr::bind_rows(MOD_DYN_P, MOD_DYN_R)

readr::write_csv(
  MOD_DYN,
  file.path(DIR_TAB, "02_WGCNA_module_dynamic_fraction_summary.csv")
)

# Cluster composition by module
CLUSTER_COMPOSITION <- dplyr::bind_rows(
  P$integrated_genes,
  R$integrated_genes
) |>
  dplyr::count(Genotype, Cluster, ModuleColor, name = "N") |>
  dplyr::group_by(Genotype, Cluster) |>
  dplyr::mutate(Fraction = N / sum(N)) |>
  dplyr::ungroup()

readr::write_csv(
  CLUSTER_COMPOSITION,
  file.path(DIR_TAB, "02_maSigPro_cluster_composition_by_WGCNA_module.csv")
)

# ----------------------------------------------------------------------
# 8. DYNAMIC HUB GENES
# ----------------------------------------------------------------------

DYNAMIC_GENES_ALL <- dplyr::bind_rows(
  P$integrated_genes,
  R$integrated_genes
)

readr::write_csv(
  DYNAMIC_GENES_ALL,
  file.path(DIR_TAB, "03_All_maSigPro_dynamic_genes_with_WGCNA_membership.csv")
)

DYNAMIC_HUBS <- DYNAMIC_GENES_ALL |>
  dplyr::filter(DynamicHub) |>
  dplyr::arrange(
    Genotype,
    ModuleColor,
    Cluster,
    dplyr::desc(kME),
    dplyr::desc(abs(GS_Treatment))
  )

readr::write_csv(
  DYNAMIC_HUBS,
  file.path(DIR_TAB, "03_Dynamic_WGCNA_hubs_kME085.csv")
)

DYNAMIC_TREATMENT_HUBS <- DYNAMIC_GENES_ALL |>
  dplyr::filter(DynamicTreatmentHub) |>
  dplyr::arrange(
    Genotype,
    ModuleColor,
    Cluster,
    dplyr::desc(kME),
    dplyr::desc(abs(GS_Treatment))
  )

readr::write_csv(
  DYNAMIC_TREATMENT_HUBS,
  file.path(DIR_TAB, "03_Dynamic_treatment_priority_hubs.csv")
)

HUB_SUMMARY <- DYNAMIC_GENES_ALL |>
  dplyr::group_by(Genotype, ModuleColor, Cluster) |>
  dplyr::summarise(
    DynamicGenes = dplyr::n(),
    DynamicHubs = sum(DynamicHub, na.rm = TRUE),
    DynamicTreatmentHubs = sum(DynamicTreatmentHub, na.rm = TRUE),
    Median_kME = median(kME, na.rm = TRUE),
    Median_abs_GS_Treatment = median(abs(GS_Treatment), na.rm = TRUE),
    .groups = "drop"
  )

readr::write_csv(
  HUB_SUMMARY,
  file.path(DIR_TAB, "03_Dynamic_hub_summary_by_module_cluster.csv")
)

# ----------------------------------------------------------------------
# 9. FUNCTIONAL CONVERGENCE — EXACT GO / PATHWAY TERMS
# ----------------------------------------------------------------------

exact_functional_overlap_one <- function(x, overlap_df) {
  
  sig_pairs <- overlap_df |>
    dplyr::filter(
      Significant,
      ModuleColor != "grey"
    ) |>
    dplyr::select(Genotype, Cluster, ModuleColor, Overlap, FoldEnrichment, FDR)
  
  if (nrow(sig_pairs) == 0) return(tibble::tibble())
  
  out <- list()
  ii <- 1
  
  for (i in seq_len(nrow(sig_pairs))) {
    
    cl <- sig_pairs$Cluster[i]
    mo <- sig_pairs$ModuleColor[i]
    
    # GO BP
    wgo <- x$w_go |>
      dplyr::filter(ModuleColor == mo) |>
      dplyr::select(ID, Description, WGCNA_FDR = p.adjust)
    
    mgo <- x$m_go |>
      dplyr::filter(Cluster == cl) |>
      dplyr::select(ID, Description, maSigPro_FDR = p.adjust)
    
    shared_go <- dplyr::inner_join(
      wgo, mgo,
      by = c("ID", "Description")
    )
    
    # Pathway
    wp <- x$w_path |>
      dplyr::filter(ModuleColor == mo) |>
      dplyr::select(ID, Description, WGCNA_FDR = p.adjust)
    
    mp <- x$m_path |>
      dplyr::filter(Cluster == cl) |>
      dplyr::select(ID, Description, maSigPro_FDR = p.adjust)
    
    shared_path <- dplyr::inner_join(
      wp, mp,
      by = c("ID", "Description")
    )
    
    go_union <- union(wgo$ID, mgo$ID)
    path_union <- union(wp$ID, mp$ID)
    
    out[[ii]] <- tibble::tibble(
      Genotype = x$genotype,
      Cluster = cl,
      ModuleColor = mo,
      GeneOverlap = sig_pairs$Overlap[i],
      GeneFoldEnrichment = sig_pairs$FoldEnrichment[i],
      GeneOverlap_FDR = sig_pairs$FDR[i],
      N_WGCNA_GO = dplyr::n_distinct(wgo$ID),
      N_maSigPro_GO = dplyr::n_distinct(mgo$ID),
      N_shared_GO = dplyr::n_distinct(shared_go$ID),
      GO_Jaccard = ifelse(
        length(go_union) > 0,
        length(intersect(wgo$ID, mgo$ID)) / length(go_union),
        NA_real_
      ),
      Shared_GO_terms = paste(
        unique(shared_go$Description),
        collapse = "; "
      ),
      N_WGCNA_Pathways = dplyr::n_distinct(wp$ID),
      N_maSigPro_Pathways = dplyr::n_distinct(mp$ID),
      N_shared_Pathways = dplyr::n_distinct(shared_path$ID),
      Pathway_Jaccard = ifelse(
        length(path_union) > 0,
        length(intersect(wp$ID, mp$ID)) / length(path_union),
        NA_real_
      ),
      Shared_Pathways = paste(
        unique(shared_path$Description),
        collapse = "; "
      )
    )
    
    ii <- ii + 1
  }
  
  dplyr::bind_rows(out)
}

FUNC_EXACT <- dplyr::bind_rows(
  exact_functional_overlap_one(P, OV_P),
  exact_functional_overlap_one(R, OV_R)
)

readr::write_csv(
  FUNC_EXACT,
  file.path(DIR_TAB, "04_Functional_convergence_exact_GO_pathway_terms.csv")
)

# ----------------------------------------------------------------------
# 10. FUNCTIONAL CONVERGENCE — BIOLOGICAL THEMES
# ----------------------------------------------------------------------

theme_table_one <- function(x) {
  
  # WGCNA themes
  w_terms <- dplyr::bind_rows(
    x$w_go |>
      dplyr::transmute(
        Unit = ModuleColor,
        UnitType = "WGCNA module",
        Description,
        p.adjust,
        Source = "GO"
      ),
    x$w_path |>
      dplyr::transmute(
        Unit = ModuleColor,
        UnitType = "WGCNA module",
        Description,
        p.adjust,
        Source = "Pathway"
      )
  ) |>
    dplyr::mutate(Theme = classify_theme(Description)) |>
    dplyr::filter(Theme != "Other") |>
    dplyr::group_by(Unit, UnitType, Theme) |>
    dplyr::summarise(
      N_terms = dplyr::n(),
      Best_FDR = min(p.adjust, na.rm = TRUE),
      Representative_terms =
        paste(head(Description[order(p.adjust)], 5), collapse = "; "),
      .groups = "drop"
    )
  
  # maSigPro themes
  m_terms <- dplyr::bind_rows(
    x$m_go |>
      dplyr::transmute(
        Unit = as.character(Cluster),
        UnitType = "maSigPro cluster",
        Description,
        p.adjust,
        Source = "GO"
      ),
    x$m_path |>
      dplyr::transmute(
        Unit = as.character(Cluster),
        UnitType = "maSigPro cluster",
        Description,
        p.adjust,
        Source = "Pathway"
      )
  ) |>
    dplyr::mutate(Theme = classify_theme(Description)) |>
    dplyr::filter(Theme != "Other") |>
    dplyr::group_by(Unit, UnitType, Theme) |>
    dplyr::summarise(
      N_terms = dplyr::n(),
      Best_FDR = min(p.adjust, na.rm = TRUE),
      Representative_terms =
        paste(head(Description[order(p.adjust)], 5), collapse = "; "),
      .groups = "drop"
    )
  
  # Cross only significant enriched cluster-module relationships
  pairs <- if (x$genotype == "PRESTO") OV_P else OV_R
  
  pairs <- pairs |>
    dplyr::filter(
      Significant,
      ModuleColor != "grey"
    )
  
  out <- list()
  ii <- 1
  
  for (i in seq_len(nrow(pairs))) {
    cl <- as.character(pairs$Cluster[i])
    mo <- pairs$ModuleColor[i]
    
    wt <- w_terms |>
      dplyr::filter(Unit == mo) |>
      dplyr::rename(
        WGCNA_N_terms = N_terms,
        WGCNA_Best_FDR = Best_FDR,
        WGCNA_terms = Representative_terms
      ) |>
      dplyr::select(
        Theme, WGCNA_N_terms,
        WGCNA_Best_FDR, WGCNA_terms
      )
    
    mt <- m_terms |>
      dplyr::filter(Unit == cl) |>
      dplyr::rename(
        maSigPro_N_terms = N_terms,
        maSigPro_Best_FDR = Best_FDR,
        maSigPro_terms = Representative_terms
      ) |>
      dplyr::select(
        Theme, maSigPro_N_terms,
        maSigPro_Best_FDR, maSigPro_terms
      )
    
    both <- dplyr::inner_join(wt, mt, by = "Theme")
    
    if (nrow(both) > 0) {
      out[[ii]] <- both |>
        dplyr::mutate(
          Genotype = x$genotype,
          Cluster = as.integer(cl),
          ModuleColor = mo,
          GeneOverlap = pairs$Overlap[i],
          GeneFoldEnrichment = pairs$FoldEnrichment[i],
          GeneOverlap_FDR = pairs$FDR[i],
          CombinedFunctionalScore =
            -log10(
              pmax(
                sqrt(WGCNA_Best_FDR * maSigPro_Best_FDR),
                1e-300
              )
            )
        )
      ii <- ii + 1
    }
  }
  
  if (length(out) == 0) return(tibble::tibble())
  dplyr::bind_rows(out)
}

THEME_CONV <- dplyr::bind_rows(
  theme_table_one(P),
  theme_table_one(R)
)

readr::write_csv(
  THEME_CONV,
  file.path(DIR_TAB, "04_Biological_theme_convergence_WGCNA_maSigPro.csv")
)

# ----------------------------------------------------------------------
# 11. CROSS-GENOTYPE SHARED DYNAMIC GENES
# ----------------------------------------------------------------------

P_dyn <- P$integrated_genes |>
  dplyr::select(
    Gene,
    PRESTO_Cluster = Cluster,
    PRESTO_Module = ModuleColor,
    PRESTO_kME = kME,
    PRESTO_GS_Treatment = GS_Treatment
  )

R_dyn <- R$integrated_genes |>
  dplyr::select(
    Gene,
    ROBILA_Cluster = Cluster,
    ROBILA_Module = ModuleColor,
    ROBILA_kME = kME,
    ROBILA_GS_Treatment = GS_Treatment
  )

SHARED_DYNAMIC <- dplyr::inner_join(P_dyn, R_dyn, by = "Gene")

readr::write_csv(
  SHARED_DYNAMIC,
  file.path(DIR_TAB, "05_Shared_dynamic_genes_PRESTO_ROBILA.csv")
)

CLUSTER_FLOW_PR <- SHARED_DYNAMIC |>
  dplyr::count(
    PRESTO_Cluster, ROBILA_Cluster,
    name = "SharedGenes"
  ) |>
  dplyr::group_by(PRESTO_Cluster) |>
  dplyr::mutate(
    Fraction_of_PRESTO_cluster_shared_set =
      SharedGenes / sum(SharedGenes)
  ) |>
  dplyr::ungroup()

readr::write_csv(
  CLUSTER_FLOW_PR,
  file.path(DIR_TAB, "05_PRESTO_vs_ROBILA_maSigPro_cluster_mapping.csv")
)

MODULE_FLOW_PR <- SHARED_DYNAMIC |>
  dplyr::count(
    PRESTO_Module, ROBILA_Module,
    name = "SharedGenes"
  ) |>
  dplyr::group_by(PRESTO_Module) |>
  dplyr::mutate(
    Fraction_of_PRESTO_module_shared_set =
      SharedGenes / sum(SharedGenes)
  ) |>
  dplyr::ungroup()

readr::write_csv(
  MODULE_FLOW_PR,
  file.path(DIR_TAB, "05_PRESTO_vs_ROBILA_WGCNA_module_mapping_dynamicGenes.csv")
)

# ----------------------------------------------------------------------
# 12. PUBLICATION SUMMARY TABLE
# ----------------------------------------------------------------------

THEME_PAIR_SUMMARY <- if (nrow(THEME_CONV) > 0) {
  THEME_CONV |>
    dplyr::group_by(Genotype, Cluster, ModuleColor) |>
    dplyr::summarise(
      Shared_biological_themes = paste(
        sort(unique(Theme)),
        collapse = "; "
      ),
      N_shared_biological_themes = dplyr::n_distinct(Theme),
      .groups = "drop"
    )
} else {
  tibble::tibble(
    Genotype = character(),
    Cluster = integer(),
    ModuleColor = character(),
    Shared_biological_themes = character(),
    N_shared_biological_themes = integer()
  )
}

PUBLICATION_SUMMARY <- OV_ALL |>
  dplyr::filter(
    ModuleColor != "grey",
    Significant
  ) |>
  dplyr::left_join(
    HUB_SUMMARY,
    by = c("Genotype", "ModuleColor", "Cluster")
  ) |>
  dplyr::left_join(
    FUNC_EXACT |>
      dplyr::select(
        Genotype, Cluster, ModuleColor,
        N_shared_GO, GO_Jaccard,
        N_shared_Pathways, Pathway_Jaccard,
        Shared_GO_terms, Shared_Pathways
      ),
    by = c("Genotype", "Cluster", "ModuleColor")
  ) |>
  dplyr::left_join(
    THEME_PAIR_SUMMARY,
    by = c("Genotype", "Cluster", "ModuleColor")
  ) |>
  dplyr::arrange(
    Genotype,
    Cluster,
    FDR,
    dplyr::desc(FoldEnrichment)
  )

readr::write_csv(
  PUBLICATION_SUMMARY,
  file.path(DIR_TAB, "06_Publication_WGCNA_maSigPro_integration_summary.csv")
)

# ----------------------------------------------------------------------
# 13. FIGURES — CLUSTER × MODULE ENRICHMENT
# ----------------------------------------------------------------------

plot_overlap_heatmap <- function(df, genotype) {
  
  dd <- df |>
    dplyr::filter(Genotype == genotype)
  
  if (EXCLUDE_GREY_FROM_MAIN_PLOTS) {
    dd <- dd |>
      dplyr::filter(ModuleColor != "grey")
  }
  
  # Stable module order based on WGCNA colour conventions and observed modules
  preferred <- c(
    "brown", "blue", "turquoise",
    "yellow", "green", "red", "black", "pink", "grey"
  )
  present <- unique(dd$ModuleColor)
  ord <- c(preferred[preferred %in% present], setdiff(present, preferred))
  
  dd <- dd |>
    dplyr::mutate(
      ModuleColor = factor(ModuleColor, levels = ord),
      Cluster = factor(Cluster, levels = sort(unique(Cluster)))
    )
  
  ggplot(
    dd,
    aes(ModuleColor, Cluster)
  ) +
    geom_tile(
      aes(fill = log2FE),
      colour = "white",
      linewidth = 0.4
    ) +
    geom_text(
      aes(
        label = ifelse(Significant, paste0(Overlap, "*"), Overlap)
      ),
      size = 3.2
    ) +
    scale_fill_gradient2(
      low = "steelblue4",
      mid = "white",
      high = "firebrick3",
      midpoint = 0,
      name = "log2 fold\nenrichment"
    ) +
    labs(
      title = paste0(genotype, ": maSigPro clusters within WGCNA modules"),
      subtitle = paste0(
        "Numbers = overlapping genes; * = enrichment FDR < ",
        FDR_CUTOFF,
        " and overlap ≥ ", MIN_OVERLAP_GENES
      ),
      x = "WGCNA module",
      y = "maSigPro cluster"
    ) +
    theme_pub(10) +
    theme(
      axis.text.x = element_text(
        angle = 45, hjust = 1, face = "bold"
      )
    )
}

p_ov_P <- plot_overlap_heatmap(OV_ALL, "PRESTO")
p_ov_R <- plot_overlap_heatmap(OV_ALL, "ROBILA")

save_plot(
  p_ov_P,
  "01_ClusterModule_enrichment_PRESTO",
  8.5, 6.5
)
save_plot(
  p_ov_R,
  "01_ClusterModule_enrichment_ROBILA",
  9.5, 6.5
)

# ----------------------------------------------------------------------
# 14. FIGURE — DYNAMIC FRACTION OF WGCNA MODULES
# ----------------------------------------------------------------------

plot_mod_dyn <- MOD_DYN

if (EXCLUDE_GREY_FROM_MAIN_PLOTS) {
  plot_mod_dyn <- plot_mod_dyn |>
    dplyr::filter(ModuleColor != "grey")
}

p_dyn <- ggplot(
  plot_mod_dyn,
  aes(
    x = reorder(ModuleColor, Dynamic_percent),
    y = Dynamic_percent
  )
) +
  geom_col(width = 0.7) +
  geom_text(
    aes(
      label = paste0(
        Dynamic_genes, "/",
        WGCNA_module_size,
        "\n(",
        round(Dynamic_percent, 1),
        "%)"
      )
    ),
    vjust = -0.25,
    size = 3
  ) +
  facet_wrap(~ Genotype, scales = "free_x") +
  coord_cartesian(
    ylim = c(
      0,
      max(plot_mod_dyn$Dynamic_percent, na.rm = TRUE) * 1.18
    )
  ) +
  labs(
    title = "Fraction of each WGCNA module represented by maSigPro-dynamic genes",
    subtitle = ifelse(
      STRICT_ONLY,
      "Strict T.fit maSigPro gene set",
      "p.vector-selected maSigPro gene set"
    ),
    x = "WGCNA module",
    y = "Dynamic genes (%)"
  ) +
  theme_pub(10) +
  theme(
    axis.text.x = element_text(
      angle = 45, hjust = 1, face = "bold"
    )
  )

save_plot(
  p_dyn,
  "02_WGCNA_module_dynamic_fraction_bothGenotypes",
  11, 5.5
)

# ----------------------------------------------------------------------
# 15. FIGURE — CLUSTER COMPOSITION BY WGCNA MODULE
# ----------------------------------------------------------------------

cc_plot <- CLUSTER_COMPOSITION

if (EXCLUDE_GREY_FROM_MAIN_PLOTS) {
  cc_plot <- cc_plot |>
    dplyr::filter(ModuleColor != "grey") |>
    dplyr::group_by(Genotype, Cluster) |>
    dplyr::mutate(Fraction_nonGrey = N / sum(N)) |>
    dplyr::ungroup()
} else {
  cc_plot <- cc_plot |>
    dplyr::mutate(Fraction_nonGrey = Fraction)
}

# Use module colour names as actual fill colours when valid
module_cols <- unique(cc_plot$ModuleColor)
fill_values <- setNames(module_cols, module_cols)

p_comp <- ggplot(
  cc_plot,
  aes(
    x = factor(Cluster),
    y = Fraction_nonGrey,
    fill = ModuleColor
  )
) +
  geom_col(width = 0.78) +
  facet_wrap(~ Genotype, scales = "free_x") +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1)
  ) +
  scale_fill_manual(
    values = fill_values,
    name = "WGCNA module"
  ) +
  labs(
    title = "WGCNA module composition of maSigPro temporal clusters",
    x = "maSigPro cluster",
    y = "Fraction of cluster"
  ) +
  theme_pub(10)

save_plot(
  p_comp,
  "02_maSigPro_cluster_composition_by_WGCNA_module",
  11.5, 5.5
)

# ----------------------------------------------------------------------
# 16. FIGURE — BIOLOGICAL THEME CONVERGENCE
# ----------------------------------------------------------------------

if (nrow(THEME_CONV) > 0) {
  
  # Keep only pairs with strongest supported gene overlap and all shared themes
  tc <- THEME_CONV |>
    dplyr::mutate(
      Pair = paste0(
        "C", Cluster, " ↔ ", ModuleColor
      ),
      Pair = forcats::fct_reorder(
        Pair,
        GeneOverlap,
        .fun = max
      )
    )
  
  p_theme <- ggplot(
    tc,
    aes(
      x = Theme,
      y = Pair
    )
  ) +
    geom_point(
      aes(
        size = pmin(
          WGCNA_N_terms,
          maSigPro_N_terms
        ),
        fill = CombinedFunctionalScore
      ),
      shape = 21,
      colour = "grey20",
      stroke = 0.35
    ) +
    facet_wrap(
      ~ Genotype,
      scales = "free_y"
    ) +
    scale_size_area(
      max_size = 11,
      name = "Min. no. enriched\nterms"
    ) +
    scale_fill_gradient(
      low = "white",
      high = "firebrick",
      name = "Functional\nsupport"
    ) +
    labs(
      title = "Biological convergence between WGCNA modules and maSigPro clusters",
      subtitle = "A theme is shown only when independently enriched in both analyses",
      x = NULL,
      y = "maSigPro cluster ↔ WGCNA module"
    ) +
    theme_pub(9) +
    theme(
      axis.text.x = element_text(
        angle = 45, hjust = 1
      )
    )
  
  save_plot(
    p_theme,
    "04_Biological_theme_convergence_WGCNA_maSigPro",
    14, 8
  )
}

# ----------------------------------------------------------------------
# 17. FIGURE — DYNAMIC HUBS
# ----------------------------------------------------------------------

hub_plot <- HUB_SUMMARY |>
  dplyr::filter(ModuleColor != "grey")

p_hubs <- ggplot(
  hub_plot,
  aes(
    x = factor(Cluster),
    y = ModuleColor,
    size = DynamicHubs,
    fill = DynamicTreatmentHubs
  )
) +
  geom_point(
    shape = 21,
    colour = "grey20",
    stroke = 0.4
  ) +
  facet_wrap(
    ~ Genotype,
    scales = "free"
  ) +
  scale_size_area(
    max_size = 12,
    name = paste0("Dynamic hubs\nkME ≥ ", KME_HUB_CUTOFF)
  ) +
  scale_fill_gradient(
    low = "white",
    high = "firebrick",
    name = paste0(
      "Treatment-priority\nhubs |GS| ≥ ",
      GS_TREATMENT_CUTOFF
    )
  ) +
  labs(
    title = "Dynamic WGCNA hubs captured by maSigPro trajectories",
    x = "maSigPro cluster",
    y = "WGCNA module"
  ) +
  theme_pub(10)

save_plot(
  p_hubs,
  "03_Dynamic_hubs_by_WGCNA_module_and_maSigPro_cluster",
  11, 6.5
)

# ----------------------------------------------------------------------
# 18. FIGURE — CROSS-GENOTYPE maSigPro CLUSTER MAPPING
# ----------------------------------------------------------------------

if (nrow(SHARED_DYNAMIC) > 0) {
  
  p_cross_cluster <- ggplot(
    CLUSTER_FLOW_PR,
    aes(
      x = factor(ROBILA_Cluster),
      y = factor(PRESTO_Cluster)
    )
  ) +
    geom_tile(
      aes(fill = SharedGenes),
      colour = "white",
      linewidth = 0.4
    ) +
    geom_text(
      aes(label = SharedGenes),
      size = 3
    ) +
    scale_fill_gradient(
      low = "white",
      high = "firebrick",
      name = "Shared genes"
    ) +
    labs(
      title = "Cross-genotype mapping of shared maSigPro-dynamic genes",
      subtitle = paste0(
        nrow(SHARED_DYNAMIC),
        " genes dynamic in both PRESTO and ROBILA"
      ),
      x = "ROBILA maSigPro cluster",
      y = "PRESTO maSigPro cluster"
    ) +
    theme_pub(10)
  
  save_plot(
    p_cross_cluster,
    "05_SharedDynamicGenes_PRESTO_vs_ROBILA_cluster_mapping",
    7.5, 6.5
  )
}

# ----------------------------------------------------------------------
# 19. OPTIONAL COMPOSITE PUBLICATION FIGURE
# ----------------------------------------------------------------------

if (requireNamespace("patchwork", quietly = TRUE)) {
  
  suppressPackageStartupMessages(library(patchwork))
  
  panels <- list(
    p_ov_P + labs(title = "A  PRESTO cluster–module integration"),
    p_ov_R + labs(title = "B  ROBILA cluster–module integration"),
    p_dyn + labs(title = "C  Dynamic fraction of WGCNA modules"),
    p_hubs + labs(title = "D  Dynamic hubs")
  )
  
  if (exists("p_theme")) {
    composite <- (
      panels[[1]] | panels[[2]]
    ) / (
      panels[[3]] | panels[[4]]
    ) / (
      p_theme + labs(title = "E  Functional convergence")
    ) +
      patchwork::plot_layout(
        heights = c(1, 0.85, 1.05)
      ) +
      patchwork::plot_annotation(
        title = paste0(
          "WGCNA × maSigPro integration",
          ifelse(
            STRICT_ONLY,
            " (strict maSigPro T.fit genes)",
            " (p.vector-selected genes)"
          )
        )
      )
    
    ggplot2::ggsave(
      file.path(
        DIR_FIG,
        "06_Publication_WGCNA_maSigPro_integration.png"
      ),
      composite,
      width = 17,
      height = 17,
      dpi = 400
    )
    ggplot2::ggsave(
      file.path(
        DIR_FIG,
        "06_Publication_WGCNA_maSigPro_integration.pdf"
      ),
      composite,
      width = 17,
      height = 17
    )
  }
}

# ----------------------------------------------------------------------
# 20. ANALYSIS SUMMARY
# ----------------------------------------------------------------------

summary_tbl <- dplyr::bind_rows(
  tibble::tibble(
    Genotype = "PRESTO",
    WGCNA_common_universe = length(P$universe),
    maSigPro_dynamic_used = nrow(P$clusters),
    Significant_cluster_module_pairs =
      sum(
        OV_P$Significant &
          OV_P$ModuleColor != "grey",
        na.rm = TRUE
      ),
    Dynamic_hubs = sum(
      P$integrated_genes$DynamicHub,
      na.rm = TRUE
    ),
    Dynamic_treatment_hubs = sum(
      P$integrated_genes$DynamicTreatmentHub,
      na.rm = TRUE
    )
  ),
  tibble::tibble(
    Genotype = "ROBILA",
    WGCNA_common_universe = length(R$universe),
    maSigPro_dynamic_used = nrow(R$clusters),
    Significant_cluster_module_pairs =
      sum(
        OV_R$Significant &
          OV_R$ModuleColor != "grey",
        na.rm = TRUE
      ),
    Dynamic_hubs = sum(
      R$integrated_genes$DynamicHub,
      na.rm = TRUE
    ),
    Dynamic_treatment_hubs = sum(
      R$integrated_genes$DynamicTreatmentHub,
      na.rm = TRUE
    )
  )
)

readr::write_csv(
  summary_tbl,
  file.path(DIR_TAB, "00_Integration_analysis_summary.csv")
)

# ----------------------------------------------------------------------
# 21. SAVE R OBJECTS
# ----------------------------------------------------------------------

saveRDS(
  list(
    PRESTO = P,
    ROBILA = R,
    Overlap = OV_ALL,
    ModuleDynamic = MOD_DYN,
    FunctionalExact = FUNC_EXACT,
    FunctionalThemes = THEME_CONV,
    SharedDynamicGenes = SHARED_DYNAMIC,
    PublicationSummary = PUBLICATION_SUMMARY
  ),
  file.path(
    DIR_RDS,
    "WGCNA_maSigPro_integration_objects.rds"
  )
)

# ----------------------------------------------------------------------
# 22. HUMAN-READABLE REPORT
# ----------------------------------------------------------------------

report <- c(
  "WGCNA × maSigPro integration",
  "========================================",
  "",
  paste0("STRICT_ONLY = ", STRICT_ONLY),
  paste0("FDR_CUTOFF = ", FDR_CUTOFF),
  paste0("MIN_OVERLAP_GENES = ", MIN_OVERLAP_GENES),
  paste0("KME_HUB_CUTOFF = ", KME_HUB_CUTOFF),
  paste0("GS_TREATMENT_CUTOFF = ", GS_TREATMENT_CUTOFF),
  "",
  "PRESTO",
  paste0("  Common analytical universe: ", length(P$universe)),
  paste0("  maSigPro dynamic genes used: ", nrow(P$clusters)),
  paste0(
    "  Significant cluster-module enrichments: ",
    sum(
      OV_P$Significant & OV_P$ModuleColor != "grey",
      na.rm = TRUE
    )
  ),
  paste0(
    "  Dynamic hubs: ",
    sum(P$integrated_genes$DynamicHub, na.rm = TRUE)
  ),
  paste0(
    "  Dynamic treatment-priority hubs: ",
    sum(P$integrated_genes$DynamicTreatmentHub, na.rm = TRUE)
  ),
  "",
  "ROBILA",
  paste0("  Common analytical universe: ", length(R$universe)),
  paste0("  maSigPro dynamic genes used: ", nrow(R$clusters)),
  paste0(
    "  Significant cluster-module enrichments: ",
    sum(
      OV_R$Significant & OV_R$ModuleColor != "grey",
      na.rm = TRUE
    )
  ),
  paste0(
    "  Dynamic hubs: ",
    sum(R$integrated_genes$DynamicHub, na.rm = TRUE)
  ),
  paste0(
    "  Dynamic treatment-priority hubs: ",
    sum(R$integrated_genes$DynamicTreatmentHub, na.rm = TRUE)
  ),
  "",
  paste0(
    "Shared dynamic genes between genotypes: ",
    nrow(SHARED_DYNAMIC)
  ),
  "",
  "Key output tables:",
  "  01_Cluster_x_WGCNA_Module_overlap_enrichment.csv",
  "  01_Best_WGCNA_module_per_maSigPro_cluster.csv",
  "  02_WGCNA_module_dynamic_fraction_summary.csv",
  "  03_Dynamic_treatment_priority_hubs.csv",
  "  04_Functional_convergence_exact_GO_pathway_terms.csv",
  "  04_Biological_theme_convergence_WGCNA_maSigPro.csv",
  "  05_Shared_dynamic_genes_PRESTO_ROBILA.csv",
  "  06_Publication_WGCNA_maSigPro_integration_summary.csv",
  "",
  "Key figures:",
  "  01_ClusterModule_enrichment_PRESTO.*",
  "  01_ClusterModule_enrichment_ROBILA.*",
  "  02_WGCNA_module_dynamic_fraction_bothGenotypes.*",
  "  02_maSigPro_cluster_composition_by_WGCNA_module.*",
  "  03_Dynamic_hubs_by_WGCNA_module_and_maSigPro_cluster.*",
  "  04_Biological_theme_convergence_WGCNA_maSigPro.*",
  "  05_SharedDynamicGenes_PRESTO_vs_ROBILA_cluster_mapping.*",
  "  06_Publication_WGCNA_maSigPro_integration.*"
)

writeLines(
  report,
  file.path(
    DIR_REPORT,
    "WGCNA_maSigPro_integration_summary.txt"
  )
)

# ----------------------------------------------------------------------
# 23. PARAMETERS / SESSION
# ----------------------------------------------------------------------

param_tbl <- tibble::tibble(
  Parameter = c(
    "STRICT_ONLY",
    "FDR_CUTOFF",
    "MIN_OVERLAP_GENES",
    "KME_HUB_CUTOFF",
    "GS_TREATMENT_CUTOFF",
    "FUNCTION_FDR",
    "PRESTO_WGCNA",
    "ROBILA_WGCNA",
    "PRESTO_MASIG",
    "ROBILA_MASIG"
  ),
  Value = c(
    as.character(STRICT_ONLY),
    as.character(FDR_CUTOFF),
    as.character(MIN_OVERLAP_GENES),
    as.character(KME_HUB_CUTOFF),
    as.character(GS_TREATMENT_CUTOFF),
    as.character(FUNCTION_FDR),
    PRESTO_WGCNA,
    ROBILA_WGCNA,
    PRESTO_MASIG,
    ROBILA_MASIG
  )
)

readr::write_csv(
  param_tbl,
  file.path(DIR_TAB, "99_Integration_parameters.csv")
)

capture.output(
  sessionInfo(),
  file = file.path(OUTPUT_ROOT, "sessionInfo.txt")
)

message("\n============================================================")
message("WGCNA × maSigPro integration completed.")
message("Output: ", OUTPUT_ROOT)
message("============================================================\n")
