############################################################
## 10_integration_PRESTO_ROBILA_v2_CORRECTED.R
## Publication-oriented integration of independent WGCNA networks
## PRESTO vs ROBILA
## Cleaned version: legacy duplicate removed; figure colours harmonised.
# Repository configuration (portable; no HPC-specific paths)
if (!exists('CARROT_CONFIG_LOADED')) source('../../../config/config.R')


############################################################

## IMPORTANT
## - Run AFTER scripts 01-08 have been completed for BOTH genotypes.
## - This script DOES NOT rebuild the WGCNA networks.
## - It corrects the main limitations of the former integration:
##   (1) module overlap is calculated on the SAME common gene universe;
##   (2) overlap significance is tested by hypergeometric test + BH-FDR;
##   (3) nested/refined submodules are distinguished from truly absent modules;
##   (4) WGCNA modulePreservation is run in BOTH directions;
##   (5) module eigengenes are re-analysed with Temp × Treatment × Time,
##       with Block included as a fixed blocking factor;
##   (6) temporal Treated responses (P2/P3) are classified formally;
##   (7) ALL significant GO BP/pathway terms are compared, not only top-30;
##   (8) shared hubs are defined within homologous/nested module pairs,
##       rather than from an arbitrary global top-50 list.
############################################################

##############################
## 0. USER PARAMETERS
##############################
INSTALL_MISSING_PKGS <- FALSE
AUTO_DETECT_RUNS     <- TRUE

BASE_WGCNA <- RNA_WGCNA_ROOT
RUN_DIR_PRESTO <- file.path(BASE_WGCNA, "WGCNA_PRESTO_P2P3_YYYYMMDD_HHMMSS")
RUN_DIR_ROBILA <- file.path(BASE_WGCNA, "WGCNA_ROBILA_P2P3_YYYYMMDD_HHMMSS")
BASE_OUT      <- file.path(BASE_WGCNA, "Integration_PRESTO_vs_ROBILA_v2")

FDR_CUTOFF       <- 0.05
JACCARD_CORE_THR <- 0.30   # reciprocal homologous pair
NESTED_SHARE_THR <- 0.50   # >=50% of one module nested in the other
KME_SHARED_THR   <- 0.85   # structural shared hub threshold
GS_PRIORITY_THR  <- 0.30   # optional trait-priority threshold

RUN_MODULE_PRESERVATION <- TRUE
N_PRESERVATION_PERM     <- 10   # 200 recommended for final run; use 50-100 for testing
PRESERVATION_SEED       <- 123
PRESERVATION_MAX_SIZE   <- 5000

BLOCK_COL <- "Block"     # B1/B2 in current metadata; fitted as FIXED blocking factor
TIME_COL  <- "Sampling_time"

options(stringsAsFactors = FALSE)
options(contrasts = c("contr.sum", "contr.poly"))
set.seed(PRESERVATION_SEED)

## Publication-ready colour system (Okabe-Ito inspired + sequential palettes)
## Genotype colours are kept identical across every categorical figure.
COL_PRESTO <- "#0072B2"   # blue
COL_ROBILA <- "#D55E00"   # vermillion
COL_SHARED <- "#009E73"   # bluish green
COL_PURPLE <- "#CC79A7"   # reddish purple
COL_GREY   <- "#6B7280"

PAL_GENOTYPE <- c(PRESTO = COL_PRESTO, ROBILA = COL_ROBILA)
PAL_FUNCTION <- c(GO_BP = COL_SHARED, Pathway = COL_PURPLE)
PAL_JACCARD  <- grDevices::colorRampPalette(c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B"))(100)
PAL_FRACTION <- grDevices::colorRampPalette(c("#F7FBFF", "#DEEBF7", "#9ECAE1", "#4292C6", "#08519C"))(100)
PAL_FDR      <- grDevices::colorRampPalette(c("#FCFBFD", "#EFEDF5", "#BCBDDC", "#756BB1", "#54278F"))(100)

##############################
## 1. PACKAGES + HELPERS
##############################
load_pkg <- function(p, bioc = FALSE) {
  if (!requireNamespace(p, quietly = TRUE)) {
    if (!INSTALL_MISSING_PKGS) stop("Missing package: ", p)
    if (bioc) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install(p, ask = FALSE, update = FALSE)
    } else install.packages(p, dependencies = TRUE)
  }
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}

for (p in c("WGCNA","dplyr","tidyr","purrr","tibble","stringr",
            "ggplot2","pheatmap","car","emmeans")) load_pkg(p)

WGCNA::enableWGCNAThreads()

safe_read_csv <- function(path) {
  if (!file.exists(path)) return(NULL)
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

find_latest_dir <- function(base, pattern) {
  x <- list.files(base, full.names = TRUE, pattern = pattern)
  x <- x[file.info(x)$isdir]
  if (!length(x)) stop("No run directory found in ", base, " matching ", pattern)
  x[order(file.info(x)$mtime, decreasing = TRUE)][1]
}

save_plot <- function(p, base, w = 8, h = 6, dpi = 320) {
  ggplot2::ggsave(paste0(base, ".pdf"), p, width = w, height = h, bg = "white")
  ggplot2::ggsave(paste0(base, ".svg"), p, width = w, height = h, dpi = dpi, bg = "white")
}

theme_pub <- function(base_size = 10) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.4),
      axis.line = element_blank(),
      axis.text = element_text(colour = "black"),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold"),
      strip.background = element_rect(fill = "#F3F4F6", colour = "#4B5563", linewidth = 0.3),
      strip.text = element_text(face = "bold", colour = "#111827"),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(colour = "#111827"),
      plot.title.position = "plot"
    )
}

bh <- function(x) p.adjust(x, method = "BH")

##############################
## 2. DETECT RUNS + OUTPUT
##############################
if (AUTO_DETECT_RUNS) {
  RUN_DIR_PRESTO <- find_latest_dir(BASE_WGCNA, "^WGCNA_PRESTO_P2P3_")
  RUN_DIR_ROBILA <- find_latest_dir(BASE_WGCNA, "^WGCNA_ROBILA_P2P3_")
}
stopifnot(dir.exists(RUN_DIR_PRESTO), dir.exists(RUN_DIR_ROBILA))

ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUT_DIR <- file.path(BASE_OUT, paste0("Integration_v2_", ts))
DIR_TAB <- file.path(OUT_DIR, "01_tables")
DIR_FIG <- file.path(OUT_DIR, "02_figures")
DIR_RDS <- file.path(OUT_DIR, "03_rds")
DIR_REP <- file.path(OUT_DIR, "04_report")
for (d in c(OUT_DIR, DIR_TAB, DIR_FIG, DIR_RDS, DIR_REP)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

cat("\n============================================================\n")
cat(" WGCNA INTEGRATION v2 — PRESTO vs ROBILA\n")
cat("============================================================\n")
cat("PRESTO:", RUN_DIR_PRESTO, "\n")
cat("ROBILA:", RUN_DIR_ROBILA, "\n")
cat("OUTPUT :", OUT_DIR, "\n\n")

##############################
## 3. LOAD A COMPLETE RUN
##############################
load_run <- function(run_dir, genotype) {
  x <- list(genotype = genotype, run_dir = run_dir)
  
  x$gene_mod <- safe_read_csv(file.path(run_dir, "02_network", "Gene_to_Module.csv"))
  x$me_map   <- safe_read_csv(file.path(run_dir, "02_network", "ME_color_map.csv"))
  x$MEs      <- readRDS(file.path(run_dir, "02_network", "MEs.rds"))
  x$datExpr  <- readRDS(file.path(run_dir, "02_network", "datExpr_clean.rds"))
  x$meta     <- readRDS(file.path(run_dir, "02_network", "meta_clean.rds"))
  
  x$anova_old <- safe_read_csv(file.path(run_dir, "03_module_trait", "Module_ANOVA_results.csv"))
  x$hubs      <- safe_read_csv(file.path(run_dir, "04_hub_genes", "AllGenes_kME_GS.csv"))
  x$go        <- safe_read_csv(file.path(run_dir, "05_enrichment", "GO_ALL_modules_combined.csv"))
  x$pathway   <- safe_read_csv(file.path(run_dir, "05_enrichment", "Pathway_ALL_modules_combined.csv"))
  x$time_old  <- safe_read_csv(file.path(run_dir, "07_sampling_effect", "P2vsP3_all_tests.csv"))
  
  stopifnot(!is.null(x$gene_mod), !is.null(x$me_map), !is.null(x$MEs),
            !is.null(x$datExpr), !is.null(x$meta), !is.null(x$hubs))
  
  ## Align metadata to eigengenes exactly
  x$meta$Sample <- as.character(x$meta$Sample)
  if (is.null(rownames(x$MEs))) stop("MEs have no sample rownames in ", genotype)
  idx <- match(rownames(x$MEs), x$meta$Sample)
  if (anyNA(idx)) stop("Could not align all MEs to metadata in ", genotype)
  x$meta <- x$meta[idx, , drop = FALSE]
  stopifnot(as.character(x$meta$Sample) == rownames(x$MEs))
  
  ## Factor order
  x$meta$Temperature   <- factor(x$meta$Temperature, levels = c("T0","T1","T2"))
  x$meta$Treatment     <- factor(x$meta$Treatment, levels = c("Control","Treated"))
  x$meta[[TIME_COL]]   <- factor(x$meta[[TIME_COL]], levels = c("P2","P3"))
  if (BLOCK_COL %in% colnames(x$meta)) x$meta[[BLOCK_COL]] <- factor(x$meta[[BLOCK_COL]])
  
  ## ME column -> biological module colour
  x$me_to_color <- setNames(x$me_map$ModuleColor, x$me_map$ME)
  missing_map <- setdiff(colnames(x$MEs), names(x$me_to_color))
  if (length(missing_map)) stop("ME map missing columns: ", paste(missing_map, collapse=", "))
  
  ## Named gene -> module vector
  x$gene_color <- setNames(as.character(x$gene_mod$ModuleColor), x$gene_mod$Gene)
  
  x
}

P <- load_run(RUN_DIR_PRESTO, "PRESTO")
R <- load_run(RUN_DIR_ROBILA, "ROBILA")

##############################
## 4. COMMON GENE UNIVERSE
##############################
common_genes <- Reduce(intersect, list(
  colnames(P$datExpr), colnames(R$datExpr),
  names(P$gene_color), names(R$gene_color)
))
common_genes <- sort(unique(common_genes))
if (length(common_genes) < 1000) stop("Too few common genes: ", length(common_genes))

common_summary <- data.frame(
  Metric = c("PRESTO genes in WGCNA", "ROBILA genes in WGCNA", "Common gene universe"),
  N = c(ncol(P$datExpr), ncol(R$datExpr), length(common_genes))
)
write.csv(common_summary, file.path(DIR_TAB, "00_Common_gene_universe_summary.csv"), row.names = FALSE)
write.csv(data.frame(Gene = common_genes), file.path(DIR_TAB, "00_Common_gene_universe.csv"), row.names = FALSE)
cat("Common gene universe:", length(common_genes), "genes\n")

##############################
## 5. MODULE OVERLAP ON COMMON UNIVERSE
##    + hypergeometric significance + reciprocal best match
##############################
assign_common <- data.frame(
  Gene = common_genes,
  PRESTO = unname(P$gene_color[common_genes]),
  ROBILA = unname(R$gene_color[common_genes]),
  stringsAsFactors = FALSE
)

modP <- sort(setdiff(unique(assign_common$PRESTO), "grey"))
modR <- sort(setdiff(unique(assign_common$ROBILA), "grey"))
M <- nrow(assign_common)

ov <- purrr::map_dfr(modP, function(mp) {
  purrr::map_dfr(modR, function(mr) {
    gp <- assign_common$Gene[assign_common$PRESTO == mp]
    gr <- assign_common$Gene[assign_common$ROBILA == mr]
    x  <- length(intersect(gp, gr))
    np <- length(gp); nr <- length(gr)
    data.frame(
      PRESTO_mod = mp, ROBILA_mod = mr,
      N_PRESTO_common = np, N_ROBILA_common = nr,
      N_overlap = x,
      Jaccard_common = ifelse(np + nr - x > 0, x/(np + nr - x), 0),
      Fraction_PRESTO_captured_by_ROBILA = ifelse(np > 0, x/np, NA_real_),
      Fraction_ROBILA_nested_in_PRESTO   = ifelse(nr > 0, x/nr, NA_real_),
      p_hyper = stats::phyper(x - 1, np, M - np, nr, lower.tail = FALSE)
    )
  })
}) |>
  dplyr::mutate(FDR_hyper = bh(p_hyper))

bestP <- ov |>
  group_by(PRESTO_mod) |>
  slice_max(Jaccard_common, n = 1, with_ties = FALSE) |>
  ungroup() |>
  transmute(PRESTO_mod, Best_R_by_Jaccard = ROBILA_mod)

bestR <- ov |>
  group_by(ROBILA_mod) |>
  slice_max(Jaccard_common, n = 1, with_ties = FALSE) |>
  ungroup() |>
  transmute(ROBILA_mod, Best_P_by_Jaccard = PRESTO_mod)

ov <- ov |>
  left_join(bestP, by = "PRESTO_mod") |>
  left_join(bestR, by = "ROBILA_mod") |>
  mutate(
    Reciprocal_best = (ROBILA_mod == Best_R_by_Jaccard & PRESTO_mod == Best_P_by_Jaccard),
    Pair_class = case_when(
      FDR_hyper < FDR_CUTOFF & Reciprocal_best & Jaccard_common >= JACCARD_CORE_THR ~ "Conserved homologous",
      FDR_hyper < FDR_CUTOFF &
        pmax(Fraction_PRESTO_captured_by_ROBILA, Fraction_ROBILA_nested_in_PRESTO, na.rm = TRUE) >= NESTED_SHARE_THR ~ "Nested/refined submodule",
      FDR_hyper < FDR_CUTOFF & Jaccard_common >= 0.15 ~ "Partial conserved overlap",
      TRUE ~ "Weak/no overlap"
    )
  ) |>
  arrange(desc(Jaccard_common))

write.csv(ov, file.path(DIR_TAB, "01_Module_overlap_common_universe.csv"), row.names = FALSE)

## ROBILA modules: which PRESTO parent contains most of their genes?
robila_parent <- ov |>
  group_by(ROBILA_mod) |>
  slice_max(Fraction_ROBILA_nested_in_PRESTO, n = 1, with_ties = FALSE) |>
  ungroup() |>
  transmute(
    ROBILA_mod,
    Dominant_PRESTO_parent = PRESTO_mod,
    Parent_fraction = Fraction_ROBILA_nested_in_PRESTO,
    Jaccard_common,
    FDR_hyper,
    Interpretation = case_when(
      Parent_fraction >= 0.75 ~ "Strong subdivision of PRESTO parent",
      Parent_fraction >= 0.50 ~ "Dominant PRESTO parent / refined submodule",
      Parent_fraction >= 0.30 ~ "Mixed ancestry with partial PRESTO parent",
      TRUE ~ "Composite / weakly nested"
    )
  )
write.csv(robila_parent, file.path(DIR_TAB, "01_ROBILA_modules_as_refinements_of_PRESTO.csv"), row.names = FALSE)

presto_parent <- ov |>
  group_by(PRESTO_mod) |>
  slice_max(Fraction_PRESTO_captured_by_ROBILA, n = 1, with_ties = FALSE) |>
  ungroup() |>
  transmute(
    PRESTO_mod,
    Dominant_ROBILA_match = ROBILA_mod,
    Captured_fraction = Fraction_PRESTO_captured_by_ROBILA,
    Jaccard_common,
    FDR_hyper
  )
write.csv(presto_parent, file.path(DIR_TAB, "01_PRESTO_modules_best_ROBILA_match.csv"), row.names = FALSE)

## Heatmap: common-universe Jaccard
jac_mat <- ov |>
  dplyr::select(PRESTO_mod, ROBILA_mod, Jaccard_common) |>
  pivot_wider(names_from = ROBILA_mod, values_from = Jaccard_common, values_fill = 0) |>
  tibble::column_to_rownames("PRESTO_mod") |>
  as.matrix()

pdf(file.path(DIR_FIG, "01_ModuleOverlap_Jaccard_commonUniverse.pdf"), width = 10, height = 5)
pheatmap::pheatmap(jac_mat,
                   color = PAL_JACCARD,
                   cluster_rows = TRUE, cluster_cols = TRUE,
                   display_numbers = matrix(sprintf("%.2f", jac_mat), nrow=nrow(jac_mat), dimnames=dimnames(jac_mat)),
                   number_color = "#111827",
                   border_color = NA,
                   main = paste0("Module overlap on common gene universe (n=", M, ")\nJaccard index"))
dev.off()
grDevices::svg(file.path(DIR_FIG, "01_ModuleOverlap_Jaccard_commonUniverse.svg"), width = 10, height = 5)
pheatmap::pheatmap(jac_mat,
                   color = PAL_JACCARD,
                   cluster_rows = TRUE, cluster_cols = TRUE,
                   display_numbers = matrix(sprintf("%.2f", jac_mat), nrow=nrow(jac_mat), dimnames=dimnames(jac_mat)),
                   number_color = "#111827",
                   border_color = NA,
                   main = paste0("Module overlap on common gene universe (n=", M, ")\nJaccard index"))
dev.off()

## Refinement heatmap: each ROBILA module sums to 1 across PRESTO modules + grey
flow <- assign_common |>
  dplyr::count(ROBILA, PRESTO, name = "N") |>
  dplyr::group_by(ROBILA) |>
  dplyr::mutate(Fraction_of_ROBILA = N / sum(N)) |>
  dplyr::ungroup()
write.csv(flow, file.path(DIR_TAB, "01_Module_assignment_flow_commonUniverse.csv"), row.names = FALSE)

flow_plot <- flow |>
  filter(ROBILA != "grey") |>
  ggplot(aes(PRESTO, ROBILA, fill = Fraction_of_ROBILA)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(Fraction_of_ROBILA >= 0.08,
                               paste0(round(100*Fraction_of_ROBILA), "%"), ""),
                colour = Fraction_of_ROBILA >= 0.55), size = 3) +
  scale_colour_manual(values = c(`FALSE` = "#111827", `TRUE` = "white"), guide = "none") +
  scale_fill_gradientn(colours = PAL_FRACTION, limits = c(0,1), name = "Fraction\nof ROBILA") +
  labs(title = "ROBILA modular refinement of PRESTO-like programmes",
       subtitle = "Fraction of each ROBILA module assigned to PRESTO modules on the common gene universe",
       x = "PRESTO module", y = "ROBILA module") +
  theme_pub(10)
save_plot(flow_plot, file.path(DIR_FIG, "01_ROBILA_refinement_heatmap"), w = 7.5, h = 5.5)

##############################
## 6. FORMAL MODULE PRESERVATION — BOTH DIRECTIONS
##############################
preservation_tbl <- NULL
if (RUN_MODULE_PRESERVATION) {
  cat("Running modulePreservation with", N_PRESERVATION_PERM, "permutations...\n")
  
  datP <- as.data.frame(P$datExpr[, common_genes, drop = FALSE])
  datR <- as.data.frame(R$datExpr[, common_genes, drop = FALSE])
  stopifnot(identical(colnames(datP), colnames(datR)))
  
  colP <- unname(P$gene_color[common_genes])
  colR <- unname(R$gene_color[common_genes])
  
  multiData  <- list(PRESTO = list(data = datP), ROBILA = list(data = datR))
  multiColor <- list(PRESTO = colP, ROBILA = colR)
  
  mp <- tryCatch(
    WGCNA::modulePreservation(
      multiData = multiData,
      multiColor = multiColor,
      referenceNetworks = c(1,2),
      nPermutations = N_PRESERVATION_PERM,
      randomSeed = PRESERVATION_SEED,
      quickCor = 0,
      networkType = "signed",
      corFnc = "bicor",
      corOptions = "use = 'p'",
      maxModuleSize = PRESERVATION_MAX_SIZE,
      verbose = 3
    ),
    error = function(e) e
  )
  
  if (inherits(mp, "error")) {
    writeLines(paste("modulePreservation ERROR:", conditionMessage(mp)),
               file.path(DIR_REP, "modulePreservation_ERROR.txt"))
    warning("modulePreservation failed: ", conditionMessage(mp))
  } else {
    saveRDS(mp, file.path(DIR_RDS, "modulePreservation_bothDirections.rds"))
    
    extract_preservation <- function(mp) {
      Z <- mp$preservation$Z
      O <- mp$preservation$observed
      purrr::map_dfr(names(Z), function(refnm) {
        purrr::map_dfr(names(Z[[refnm]]), function(testnm) {
          z <- as.data.frame(Z[[refnm]][[testnm]])
          z$Module <- rownames(z)
          z$Reference <- sub("^ref\\.", "", refnm)
          z$Test <- sub("^inColumnsAlsoPresentIn\\.", "", testnm)
          rownames(z) <- NULL
          
          ## medianRank is usually stored in preservation$observed
          o <- NULL
          if (!is.null(O[[refnm]]) && !is.null(O[[refnm]][[testnm]])) {
            o <- as.data.frame(O[[refnm]][[testnm]])
            o$Module <- rownames(o)
            rownames(o) <- NULL
            keep_obs <- intersect(c("Module","medianRank.pres","medianRank.qual"), colnames(o))
            o <- o[, keep_obs, drop=FALSE]
            z <- dplyr::left_join(z, o, by="Module")
          }
          z
        })
      })
    }
    
    preservation_tbl <- extract_preservation(mp)
    
    ## WGCNA may label networks as 1/2 or by list names depending on version.
    normalise_network_name <- function(z) {
      z <- as.character(z)
      dplyr::recode(z, `1` = "PRESTO", `2` = "ROBILA",
                    PRESTO = "PRESTO", ROBILA = "ROBILA", .default = z)
    }
    preservation_tbl <- preservation_tbl |>
      mutate(Reference = normalise_network_name(Reference),
             Test = normalise_network_name(Test))
    
    zcol <- intersect(c("Zsummary.pres", "Zsummary"), colnames(preservation_tbl))[1]
    mcol <- intersect(c("medianRank.pres", "medianRank"), colnames(preservation_tbl))[1]
    
    if (!is.na(zcol)) {
      preservation_tbl <- preservation_tbl |>
        mutate(
          Zsummary = .data[[zcol]],
          medianRank = if (!is.na(mcol)) .data[[mcol]] else NA_real_,
          Preservation_class = case_when(
            Zsummary >= 10 ~ "Strong preservation",
            Zsummary >= 2  ~ "Moderate preservation",
            TRUE ~ "Weak/no evidence"
          )
        )
    }
    
    write.csv(preservation_tbl, file.path(DIR_TAB, "02_ModulePreservation_bothDirections.csv"), row.names = FALSE)
    
    if ("Zsummary" %in% colnames(preservation_tbl)) {
      pp <- preservation_tbl |>
        filter(!tolower(Module) %in% c("grey","gold")) |>
        ggplot(aes(reorder(Module, Zsummary), Zsummary, shape = Reference, colour = Reference)) +
        geom_hline(yintercept = 2, linetype = 2) +
        geom_hline(yintercept = 10, linetype = 2) +
        geom_point(size = 3) +
        scale_colour_manual(values = PAL_GENOTYPE) +
        coord_flip() +
        facet_wrap(~ paste0(Reference, " → ", Test), scales = "free_y") +
        labs(title = "Formal WGCNA module preservation",
             subtitle = paste0("Common gene universe; ", N_PRESERVATION_PERM, " permutations"),
             x = "Reference module", y = "Zsummary preservation") +
        theme_pub(10) + theme(legend.position = "none")
      save_plot(pp, file.path(DIR_FIG, "02_ModulePreservation_Zsummary"), w = 9, h = 6)
    }
  }
}

##############################
## 7. RE-ANALYSE MODULE EIGENGENES
##    FULL: Temperature × Treatment × Sampling_time + Block
##############################
extract_type3 <- function(fit) {
  a <- as.data.frame(car::Anova(fit, type = 3))
  a$Effect <- rownames(a)
  rownames(a) <- NULL
  pcol <- grep("Pr\\(>F\\)", colnames(a), value = TRUE)[1]
  if (is.na(pcol)) stop("Cannot find Type-III p-value column")
  a |>
    transmute(Effect, p_value = .data[[pcol]])
}

fit_all_ME_models <- function(x) {
  out_anova <- list()
  out_trt_TT <- list()
  out_time_TT <- list()
  out_trt_by_time <- list()
  out_temp <- list()
  
  for (me in colnames(x$MEs)) {
    mod <- unname(x$me_to_color[me])
    if (is.na(mod) || mod == "grey") next
    
    df <- x$meta
    df$ME <- as.numeric(x$MEs[, me])
    
    use_block <- BLOCK_COL %in% colnames(df) && nlevels(factor(df[[BLOCK_COL]])) > 1
    fml <- if (use_block) {
      stats::as.formula(paste0("ME ~ ", BLOCK_COL, " + Temperature * Treatment * ", TIME_COL))
    } else {
      stats::as.formula(paste0("ME ~ Temperature * Treatment * ", TIME_COL))
    }
    
    fit <- lm(fml, data = df)
    aa <- extract_type3(fit) |>
      filter(!Effect %in% c("(Intercept)", BLOCK_COL, "Residuals")) |>
      mutate(Genotype = x$genotype, ME_name = me, ModuleColor = mod,
             R2 = summary(fit)$r.squared, Adj_R2 = summary(fit)$adj.r.squared)
    out_anova[[me]] <- aa
    
    ## Treated - Control within each Temperature × Time
    em1_fml <- stats::as.formula(paste0("~ Treatment | Temperature * ", TIME_COL))
    em1 <- emmeans::emmeans(fit, em1_fml)
    ct1 <- as.data.frame(emmeans::contrast(em1, method = list("Treated-Control" = c(-1,1)), adjust = "none")) |>
      mutate(Genotype=x$genotype, ME_name=me, ModuleColor=mod,
             Family="Treatment_within_Temperature_x_Time")
    out_trt_TT[[me]] <- ct1
    
    ## P3 - P2 within each Temperature × Treatment
    em2_fml <- stats::as.formula(paste0("~ ", TIME_COL, " | Temperature * Treatment"))
    em2 <- emmeans::emmeans(fit, em2_fml)
    ct2 <- as.data.frame(emmeans::contrast(em2, method = list("P3-P2" = c(-1,1)), adjust = "none")) |>
      mutate(Genotype=x$genotype, ME_name=me, ModuleColor=mod,
             Family="Time_within_Temperature_x_Treatment")
    out_time_TT[[me]] <- ct2
    
    ## Treated - Control at each Time, averaged over Temperature
    em3_fml <- stats::as.formula(paste0("~ Treatment | ", TIME_COL))
    em3 <- emmeans::emmeans(fit, em3_fml)
    ct3 <- as.data.frame(emmeans::contrast(em3, method = list("Treated-Control" = c(-1,1)), adjust = "none")) |>
      mutate(Genotype=x$genotype, ME_name=me, ModuleColor=mod,
             Family="Treatment_by_Time_global")
    out_trt_by_time[[me]] <- ct3
    
    ## Temperature EMMs (for direction/profile; Tukey within module)
    em4 <- as.data.frame(emmeans::emmeans(fit, ~ Temperature)) |>
      mutate(Genotype=x$genotype, ME_name=me, ModuleColor=mod)
    out_temp[[me]] <- em4
  }
  
  an <- bind_rows(out_anova) |>
    group_by(Effect) |>
    mutate(FDR = bh(p_value)) |>
    ungroup()
  
  add_global_FDR <- function(z) {
    bind_rows(z) |>
      group_by(Family) |>
      mutate(FDR = bh(p.value)) |>
      ungroup()
  }
  
  list(
    anova = an,
    trt_temp_time = add_global_FDR(out_trt_TT),
    time_temp_trt = add_global_FDR(out_time_TT),
    trt_by_time = add_global_FDR(out_trt_by_time),
    temp_emm = bind_rows(out_temp)
  )
}

MEP <- fit_all_ME_models(P)
MER <- fit_all_ME_models(R)
ME_all <- bind_rows(MEP$anova, MER$anova)
TRT_TT <- bind_rows(MEP$trt_temp_time, MER$trt_temp_time)
TIME_TT <- bind_rows(MEP$time_temp_trt, MER$time_temp_trt)
TRT_TIME <- bind_rows(MEP$trt_by_time, MER$trt_by_time)
TEMP_EMM <- bind_rows(MEP$temp_emm, MER$temp_emm)

write.csv(ME_all, file.path(DIR_TAB, "03_ME_TypeIII_full_Temp_Treatment_Time.csv"), row.names = FALSE)
write.csv(TRT_TT, file.path(DIR_TAB, "03_ME_Treated_vs_Control_by_Temp_and_Time.csv"), row.names = FALSE)
write.csv(TIME_TT, file.path(DIR_TAB, "03_ME_P3_vs_P2_by_Temp_and_Treatment.csv"), row.names = FALSE)
write.csv(TRT_TIME, file.path(DIR_TAB, "03_ME_Treatment_by_Time_global.csv"), row.names = FALSE)
write.csv(TEMP_EMM, file.path(DIR_TAB, "03_ME_Temperature_EMMs.csv"), row.names = FALSE)

## Temporal classification of treatment effect using GLOBAL Treated-Control at P2/P3
make_temporal_class <- function(trt_time_df) {
  w <- trt_time_df |>
    dplyr::select(Genotype, ModuleColor, TimeValue = all_of(TIME_COL), estimate, FDR) |>
    mutate(TimeValue = as.character(TimeValue)) |>
    pivot_wider(names_from = TimeValue, values_from = c(estimate,FDR), names_sep = "_")
  
  for (nm in c("estimate_P2","estimate_P3","FDR_P2","FDR_P3")) if (!nm %in% names(w)) w[[nm]] <- NA_real_
  
  w |>
    mutate(
      sig_P2 = FDR_P2 < FDR_CUTOFF,
      sig_P3 = FDR_P3 < FDR_CUTOFF,
      same_direction = sign(estimate_P2) == sign(estimate_P3),
      Treatment_temporal_class = case_when(
        sig_P2 & sig_P3 & same_direction ~ "Persistent P2+P3",
        sig_P2 & sig_P3 & !same_direction ~ "Direction reversal",
        sig_P2 & !sig_P3 ~ "Transient P2 only",
        !sig_P2 & sig_P3 ~ "Late P3 only",
        TRUE ~ "No global treatment effect"
      )
    )
}

temporal_class <- make_temporal_class(TRT_TIME)
write.csv(temporal_class, file.path(DIR_TAB, "03_Treatment_temporal_classification.csv"), row.names = FALSE)

## Significance heatmap for FULL model
wanted_effects <- c("Temperature", "Treatment", TIME_COL,
                    "Temperature:Treatment", paste0("Temperature:", TIME_COL),
                    paste0("Treatment:", TIME_COL), paste0("Temperature:Treatment:", TIME_COL))

hm_effect <- ME_all |>
  filter(Effect %in% wanted_effects) |>
  mutate(
    Effect = factor(Effect, levels = wanted_effects,
                    labels = c("Temperature","Treatment","Time",
                               "Temp×Trt","Temp×Time","Trt×Time","Temp×Trt×Time")),
    Row = paste(Genotype, ModuleColor, sep = " | "),
    Score = pmin(-log10(pmax(FDR, 1e-300)), 12),
    Mark = case_when(FDR < 0.001 ~ "***", FDR < 0.01 ~ "**", FDR < 0.05 ~ "*", TRUE ~ "")
  )

p_eff <- ggplot(hm_effect, aes(Effect, Row, fill = Score)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = Mark, colour = Score >= 7), size = 3.5) +
  scale_colour_manual(values = c(`FALSE` = "#111827", `TRUE` = "white"), guide = "none") +
  scale_fill_gradientn(colours = PAL_FDR, limits = c(0, 12), name = "-log10(FDR)") +
  labs(title = "Module eigengene responses: full factorial model",
       subtitle = paste0("ME ~ ", BLOCK_COL, " + Temperature × Treatment × Time; BH-FDR across modules per effect"),
       x = NULL, y = NULL) +
  theme_pub(9) + theme(axis.text.x = element_text(angle = 35, hjust = 1))
save_plot(p_eff, file.path(DIR_FIG, "03_ME_fullFactorial_significance_heatmap"), w = 9.5, h = 6.5)

## Temporal classification barplot
p_tc <- temporal_class |>
  dplyr::count(Genotype, Treatment_temporal_class, name = "N") |>
  ggplot(aes(Treatment_temporal_class, N, fill = Genotype)) +
  geom_col(position = position_dodge(width=0.75), width=0.65) +
  geom_text(aes(label=N), position=position_dodge(width=0.75), vjust=-0.25, size=3.5) +
  scale_fill_manual(values = PAL_GENOTYPE) +
  coord_flip() +
  labs(title="Temporal classification of the treatment-responsive modules",
       subtitle="Global Treated-Control contrasts at P2 and P3; BH-FDR across modules",
       x=NULL, y="Number of modules") +
  theme_pub(10)
save_plot(p_tc, file.path(DIR_FIG, "03_Treatment_temporal_classification"), w = 8, h = 5)

##############################
## 8. EFFECT CONSERVATION IN MATCHED MODULE PAIRS
##############################
matched_pairs <- ov |>
  filter(Pair_class %in% c("Conserved homologous","Nested/refined submodule","Partial conserved overlap")) |>
  filter(N_overlap >= 10)

compare_effects <- function(pair_df, effects_df) {
  purrr::pmap_dfr(pair_df[,c("PRESTO_mod","ROBILA_mod","Pair_class","Jaccard_common","N_overlap")],
                  function(PRESTO_mod, ROBILA_mod, Pair_class, Jaccard_common, N_overlap) {
                    p1 <- effects_df |> filter(Genotype=="PRESTO", ModuleColor==PRESTO_mod) |> dplyr::select(Effect, FDR_PRESTO=FDR)
                    r1 <- effects_df |> filter(Genotype=="ROBILA", ModuleColor==ROBILA_mod) |> dplyr::select(Effect, FDR_ROBILA=FDR)
                    full_join(p1,r1,by="Effect") |>
                      mutate(PRESTO_mod=PRESTO_mod, ROBILA_mod=ROBILA_mod, Pair_class=Pair_class,
                             Jaccard_common=Jaccard_common, N_overlap=N_overlap,
                             Effect_status = case_when(
                               FDR_PRESTO < FDR_CUTOFF & FDR_ROBILA < FDR_CUTOFF ~ "Significant in both",
                               FDR_PRESTO < FDR_CUTOFF & !(FDR_ROBILA < FDR_CUTOFF) ~ "PRESTO only",
                               !(FDR_PRESTO < FDR_CUTOFF) & FDR_ROBILA < FDR_CUTOFF ~ "ROBILA only",
                               TRUE ~ "Neither"
                             ))
                  })
}

effect_conservation <- compare_effects(matched_pairs, ME_all)
write.csv(effect_conservation, file.path(DIR_TAB, "04_Effect_conservation_matched_modules.csv"), row.names = FALSE)

##############################
## 9. FUNCTIONAL CONSERVATION — ALL SIGNIFICANT GO BP + PATHWAYS
##############################
prep_go <- function(df) {
  if (is.null(df) || !nrow(df)) return(data.frame())
  df |>
    mutate(Module = as.character(Module), ID = as.character(ID),
           Description = as.character(Description), p.adjust = as.numeric(p.adjust),
           Ontology = as.character(Ontology)) |>
    filter(Ontology == "BP", !is.na(p.adjust), p.adjust < FDR_CUTOFF, !is.na(ID))
}
prep_pw <- function(df) {
  if (is.null(df) || !nrow(df)) return(data.frame())
  df |>
    mutate(Module = as.character(Module), ID = as.character(ID),
           Description = as.character(Description), p.adjust = as.numeric(p.adjust)) |>
    filter(!is.na(p.adjust), p.adjust < FDR_CUTOFF, !is.na(ID))
}

Pgo <- prep_go(P$go); Rgo <- prep_go(R$go)
Ppw <- prep_pw(P$pathway); Rpw <- prep_pw(R$pathway)

compare_terms_pair <- function(pair_df, A, B, label) {
  summary <- list(); details <- list()
  for (i in seq_len(nrow(pair_df))) {
    mp <- pair_df$PRESTO_mod[i]; mr <- pair_df$ROBILA_mod[i]
    a <- A |> filter(Module == mp) |> arrange(p.adjust) |> distinct(ID, .keep_all = TRUE)
    b <- B |> filter(Module == mr) |> arrange(p.adjust) |> distinct(ID, .keep_all = TRUE)
    idsA <- a$ID; idsB <- b$ID; shared <- intersect(idsA, idsB); uni <- union(idsA, idsB)
    summary[[i]] <- data.frame(
      PRESTO_mod=mp, ROBILA_mod=mr, Pair_class=pair_df$Pair_class[i],
      N_PRESTO=length(idsA), N_ROBILA=length(idsB), N_shared=length(shared),
      Functional_Jaccard=ifelse(length(uni)>0,length(shared)/length(uni),NA_real_),
      Type=label
    )
    if (length(uni)) {
      desc_map <- bind_rows(
        a |> dplyr::select(ID, Description, p.adjust) |> rename(PRESTO_padj=p.adjust),
        b |> dplyr::select(ID, Description, p.adjust) |> rename(ROBILA_padj=p.adjust)
      ) |>
        group_by(ID) |>
        summarise(Description=first(na.omit(Description)),
                  PRESTO_padj=min(PRESTO_padj,na.rm=TRUE),
                  ROBILA_padj=min(ROBILA_padj,na.rm=TRUE), .groups="drop") |>
        mutate(across(c(PRESTO_padj,ROBILA_padj), ~ifelse(is.infinite(.x),NA_real_,.x)),
               PRESTO_mod=mp, ROBILA_mod=mr, Type=label,
               Status=case_when(
                 !is.na(PRESTO_padj) & !is.na(ROBILA_padj) ~ "Shared",
                 !is.na(PRESTO_padj) ~ "PRESTO only enriched",
                 TRUE ~ "ROBILA only enriched"))
      details[[i]] <- desc_map
    }
  }
  list(summary=bind_rows(summary), details=bind_rows(details))
}

func_GO <- compare_terms_pair(matched_pairs, Pgo, Rgo, "GO_BP")
func_PW <- compare_terms_pair(matched_pairs, Ppw, Rpw, "Pathway")
func_summary <- bind_rows(func_GO$summary, func_PW$summary)
func_details <- bind_rows(func_GO$details, func_PW$details)

write.csv(func_summary, file.path(DIR_TAB, "05_Functional_conservation_matched_modules_summary.csv"), row.names = FALSE)
write.csv(func_details, file.path(DIR_TAB, "05_Functional_conservation_matched_modules_allTerms.csv"), row.names = FALSE)

## Global ALL-significant terms — descriptive overlap only (NOT 'exclusive biology')
global_overlap <- bind_rows(
  data.frame(Type="GO_BP",
             N_PRESTO=n_distinct(Pgo$ID), N_ROBILA=n_distinct(Rgo$ID),
             N_shared=length(intersect(unique(Pgo$ID), unique(Rgo$ID)))),
  data.frame(Type="Pathway",
             N_PRESTO=n_distinct(Ppw$ID), N_ROBILA=n_distinct(Rpw$ID),
             N_shared=length(intersect(unique(Ppw$ID), unique(Rpw$ID))))
) |>
  mutate(Jaccard = N_shared/(N_PRESTO + N_ROBILA - N_shared))
write.csv(global_overlap, file.path(DIR_TAB, "05_Global_allSignificant_functional_overlap.csv"), row.names = FALSE)

p_func <- func_summary |>
  filter(!is.na(Functional_Jaccard)) |>
  mutate(Pair = paste(PRESTO_mod, "↔", ROBILA_mod)) |>
  ggplot(aes(Functional_Jaccard, reorder(Pair, Functional_Jaccard), shape=Type, colour=Type)) +
  geom_point(size=3) +
  scale_colour_manual(values = PAL_FUNCTION) +
  facet_wrap(~Type, scales="free_y") +
  labs(title="Functional conservation of matched module pairs",
       subtitle="All significant terms (FDR < 0.05), not top-N lists",
       x="Functional Jaccard", y=NULL) +
  theme_pub(10)
save_plot(p_func, file.path(DIR_FIG, "05_Functional_conservation_matchedModules"), w=8, h=6)

##############################
## 10. SHARED HUBS WITHIN MATCHED MODULE PAIRS
##############################
hubP <- P$hubs |>
  transmute(Gene, PRESTO_mod=ModuleColor,
            kME_PRESTO=as.numeric(kME),
            GS_Treatment_PRESTO=as.numeric(GS_Treatment),
            GS_Temperature_PRESTO=as.numeric(GS_Temperature),
            GS_Time_PRESTO=as.numeric(GS_Time))
hubR <- R$hubs |>
  transmute(Gene, ROBILA_mod=ModuleColor,
            kME_ROBILA=as.numeric(kME),
            GS_Treatment_ROBILA=as.numeric(GS_Treatment),
            GS_Temperature_ROBILA=as.numeric(GS_Temperature),
            GS_Time_ROBILA=as.numeric(GS_Time))

hub_pairs <- purrr::pmap_dfr(matched_pairs[,c("PRESTO_mod","ROBILA_mod","Pair_class")],
                             function(PRESTO_mod, ROBILA_mod, Pair_class) {
                               mp <- PRESTO_mod; mr <- ROBILA_mod; pc <- Pair_class
                               a <- hubP[hubP$PRESTO_mod == mp, , drop=FALSE]
                               b <- hubR[hubR$ROBILA_mod == mr, , drop=FALSE]
                               inner_join(a,b,by="Gene") |>
                                 mutate(Pair_class=pc,
                                        min_kME=pmin(kME_PRESTO,kME_ROBILA,na.rm=TRUE),
                                        mean_kME=rowMeans(cbind(kME_PRESTO,kME_ROBILA),na.rm=TRUE),
                                        Shared_structural_hub=(kME_PRESTO>=KME_SHARED_THR & kME_ROBILA>=KME_SHARED_THR),
                                        Shared_treatment_priority=(Shared_structural_hub &
                                                                     abs(GS_Treatment_PRESTO)>=GS_PRIORITY_THR & abs(GS_Treatment_ROBILA)>=GS_PRIORITY_THR &
                                                                     sign(GS_Treatment_PRESTO)==sign(GS_Treatment_ROBILA)))
                             })

write.csv(hub_pairs, file.path(DIR_TAB, "06_AllGenes_in_matched_module_pairs_kME_GS.csv"), row.names = FALSE)
shared_hubs <- hub_pairs |>
  filter(Shared_structural_hub) |>
  arrange(desc(min_kME), desc(mean_kME))
write.csv(shared_hubs, file.path(DIR_TAB, "06_Shared_structural_hubs_kME085.csv"), row.names = FALSE)
write.csv(shared_hubs |> filter(Shared_treatment_priority),
          file.path(DIR_TAB, "06_Shared_treatment_priority_hubs.csv"), row.names = FALSE)

hub_summary <- shared_hubs |>
  dplyr::count(PRESTO_mod, ROBILA_mod, Pair_class, name="N_shared_hubs") |>
  arrange(desc(N_shared_hubs))
write.csv(hub_summary, file.path(DIR_TAB, "06_Shared_hubs_summary_by_pair.csv"), row.names = FALSE)

##############################
## 11. PUBLICATION-ORIENTED INTEGRATION TABLE
##############################
## Keep one dominant/biologically relevant match per ROBILA module,
## but do NOT call low-Jaccard modules 'genotype-specific' if they are nested.
pub_pairs <- robila_parent |>
  dplyr::select(ROBILA_mod, Dominant_PRESTO_parent, Parent_fraction, Interpretation) |>
  left_join(
    ov |> dplyr::select(PRESTO_mod, ROBILA_mod, Pair_class, Reciprocal_best,
                        N_PRESTO_common,N_ROBILA_common,N_overlap,Jaccard_common,
                        Fraction_PRESTO_captured_by_ROBILA,Fraction_ROBILA_nested_in_PRESTO,FDR_hyper),
    by=c("ROBILA_mod"="ROBILA_mod","Dominant_PRESTO_parent"="PRESTO_mod")
  ) |>
  left_join(
    func_summary |> filter(Type=="GO_BP") |>
      dplyr::select(PRESTO_mod,ROBILA_mod,GO_BP_Jaccard=Functional_Jaccard,N_shared_GO=N_shared),
    by=c("Dominant_PRESTO_parent"="PRESTO_mod","ROBILA_mod"="ROBILA_mod")
  ) |>
  left_join(
    func_summary |> filter(Type=="Pathway") |>
      dplyr::select(PRESTO_mod,ROBILA_mod,Pathway_Jaccard=Functional_Jaccard,N_shared_pathways=N_shared),
    by=c("Dominant_PRESTO_parent"="PRESTO_mod","ROBILA_mod"="ROBILA_mod")
  )

write.csv(pub_pairs, file.path(DIR_TAB, "07_Publication_module_integration_summary.csv"), row.names = FALSE)

##############################
## 12. AUTOMATIC TEXT REPORT
##############################
core_pairs <- ov |> filter(Pair_class=="Conserved homologous")
nested_R <- robila_parent |> filter(Parent_fraction >= NESTED_SHARE_THR)

count_effect <- function(geno, eff) {
  sum(ME_all$Genotype==geno & ME_all$Effect==eff & ME_all$FDR<FDR_CUTOFF, na.rm=TRUE)
}

report <- c(
  "WGCNA comparative integration v2 — PRESTO vs ROBILA",
  "====================================================",
  paste0("PRESTO run: ", RUN_DIR_PRESTO),
  paste0("ROBILA run: ", RUN_DIR_ROBILA),
  "",
  paste0("Common WGCNA gene universe: ", length(common_genes), " genes."),
  paste0("Conserved reciprocal homologous module pairs: ", nrow(core_pairs), "."),
  paste0("ROBILA modules with >=", round(100*NESTED_SHARE_THR), "% of genes nested in one PRESTO parent: ", nrow(nested_R), "."),
  "",
  "Full eigengene model (BH-FDR < 0.05):",
  paste0("  PRESTO Time main effect: ", count_effect("PRESTO", TIME_COL), " modules"),
  paste0("  ROBILA Time main effect: ", count_effect("ROBILA", TIME_COL), " modules"),
  paste0("  PRESTO Treatment×Time: ", count_effect("PRESTO", paste0("Treatment:", TIME_COL)), " modules"),
  paste0("  ROBILA Treatment×Time: ", count_effect("ROBILA", paste0("Treatment:", TIME_COL)), " modules"),
  paste0("  PRESTO Temp×Treatment×Time: ", count_effect("PRESTO", paste0("Temperature:Treatment:", TIME_COL)), " modules"),
  paste0("  ROBILA Temp×Treatment×Time: ", count_effect("ROBILA", paste0("Temperature:Treatment:", TIME_COL)), " modules"),
  "",
  "IMPORTANT INTERPRETATION RULE:",
  "A ROBILA module is not labelled genotype-specific solely because its Jaccard is low.",
  "If a large fraction of that module is nested in a PRESTO parent, it is classified as modular refinement/subdivision.",
  "GO/pathway overlap uses ALL FDR-significant terms; 'PRESTO only enriched' and 'ROBILA only enriched' mean",
  "not significantly enriched in the other module under the current ORA, not biological absence.",
  ""
)
writeLines(report, file.path(DIR_REP, "Integration_v2_summary.txt"))

##############################
## 13. SESSION + PARAMETERS
##############################
params <- data.frame(
  Parameter=c("FDR_CUTOFF","JACCARD_CORE_THR","NESTED_SHARE_THR","KME_SHARED_THR",
              "GS_PRIORITY_THR","N_PRESERVATION_PERM","BLOCK_COL","TIME_COL",
              "N_common_genes"),
  Value=c(FDR_CUTOFF,JACCARD_CORE_THR,NESTED_SHARE_THR,KME_SHARED_THR,
          GS_PRIORITY_THR,N_PRESERVATION_PERM,BLOCK_COL,TIME_COL,length(common_genes))
)
write.csv(params, file.path(DIR_TAB, "99_Integration_parameters.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))

cat("\n============================================================\n")
cat(" INTEGRATION v2 COMPLETE\n")
cat("============================================================\n")
cat("Output:", OUT_DIR, "\n")
cat("Key files:\n")
cat("  01_Module_overlap_common_universe.csv\n")
cat("  01_ROBILA_modules_as_refinements_of_PRESTO.csv\n")
cat("  02_ModulePreservation_bothDirections.csv\n")
cat("  03_ME_TypeIII_full_Temp_Treatment_Time.csv\n")
cat("  03_Treatment_temporal_classification.csv\n")
cat("  04_Effect_conservation_matched_modules.csv\n")
cat("  05_Functional_conservation_matched_modules_summary.csv\n")
cat("  06_Shared_structural_hubs_kME085.csv\n")
cat("  07_Publication_module_integration_summary.csv\n")
cat("\n")


############################################################
## NOTE
## The legacy script that had been concatenated after this v2 script was
## intentionally removed. It redefined parameters/functions and reverted to
## older comparison rules (e.g. top-N GO terms and non-common gene universes).
############################################################

##############################
## 11B. EXTRA PUBLICATION FIGURES
## Add this block AFTER section 11 and BEFORE the automatic text report.
##############################

## Optional PNG export helper (keeps the existing PDF/SVG export intact)
save_plot_all <- function(p, base, w = 8, h = 6, dpi = 320) {
  save_plot(p, base, w = w, h = h, dpi = dpi)
  ggplot2::ggsave(paste0(base, ".png"), p, width = w, height = h, dpi = dpi, bg = "white")
}

shorten_label <- function(x, width = 40) {
  x <- stringr::str_replace_all(x, "_", " ")
  stringr::str_wrap(x, width = width)
}

## ----- Figure 04. Conservation/divergence of factor recruitment -----
recruitment_effects <- c(
  "Temperature",
  "Treatment",
  TIME_COL,
  "Temperature:Treatment",
  paste0("Temperature:", TIME_COL),
  paste0("Treatment:", TIME_COL),
  paste0("Temperature:Treatment:", TIME_COL)
)

recruitment_labels <- c(
  "Temperature",
  "Treatment",
  "Time",
  "Temp×Trt",
  "Temp×Time",
  "Trt×Time",
  "Temp×Trt×Time"
)

pal_status <- c(
  "Significant in both" = COL_SHARED,
  "PRESTO only"        = COL_PRESTO,
  "ROBILA only"        = COL_ROBILA,
  "Neither"            = "grey90"
)

plot_effect_conservation <- effect_conservation |>
  dplyr::filter(Effect %in% recruitment_effects) |>
  dplyr::mutate(
    Pair = paste(PRESTO_mod, "↔", ROBILA_mod),
    Pair = reorder(Pair, Jaccard_common),
    Effect = factor(Effect, levels = recruitment_effects, labels = recruitment_labels),
    Effect_status = factor(Effect_status,
                           levels = c("Significant in both", "PRESTO only", "ROBILA only", "Neither"))
  )

p_effect_conservation <- ggplot(plot_effect_conservation,
                                aes(x = Effect, y = Pair, fill = Effect_status)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  scale_fill_manual(values = pal_status, drop = FALSE, name = "Recruitment status") +
  labs(
    title = "Conservation and divergence of factor recruitment in matched module pairs",
    subtitle = "Structural overlap is separated from differential Temperature/Treatment/Time recruitment",
    x = NULL, y = NULL
  ) +
  theme_pub(9) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

save_plot_all(
  p_effect_conservation,
  file.path(DIR_FIG, "04_EffectRecruitment_conservation_divergence"),
  w = 9.5, h = 6.5
)

## ----- Figure 06. Biological schematic: PRESTO-blue -> ROBILA children -----
get_top_terms <- function(go_df, pw_df, module, n = 2) {
  tb <- dplyr::bind_rows(
    go_df |>
      dplyr::filter(Module == module) |>
      dplyr::transmute(Source = "GO", Description, padj = p.adjust),
    pw_df |>
      dplyr::filter(Module == module) |>
      dplyr::transmute(Source = "Pathway", Description, padj = p.adjust)
  ) |>
    dplyr::filter(!is.na(Description), !is.na(padj)) |>
    dplyr::arrange(padj) |>
    dplyr::distinct(Description, .keep_all = TRUE) |>
    dplyr::slice_head(n = n)
  
  if (!nrow(tb)) return("No enriched term")
  paste(shorten_label(tb$Description, width = 26), collapse = "\n")
}

blue_children <- pub_pairs |>
  dplyr::filter(Dominant_PRESTO_parent == "blue", Parent_fraction >= NESTED_SHARE_THR) |>
  dplyr::arrange(desc(Parent_fraction), desc(Jaccard_common)) |>
  dplyr::select(ROBILA_mod, Parent_fraction, Pair_class, Jaccard_common)

if (nrow(blue_children) > 0) {
  presto_blue_terms <- get_top_terms(Pgo, Ppw, module = "blue", n = 2)
  
  blue_children <- blue_children |>
    dplyr::mutate(
      y = rev(seq_len(dplyr::n())),
      ROBILA_label = purrr::map_chr(ROBILA_mod, ~ get_top_terms(Rgo, Rpw, module = .x, n = 2))
    )
  
  schematic_nodes <- dplyr::bind_rows(
    data.frame(
      side = "PRESTO",
      x = 1,
      y = mean(blue_children$y),
      module = "blue",
      label = paste0("PRESTO blue\n", presto_blue_terms),
      stringsAsFactors = FALSE
    ),
    blue_children |>
      dplyr::transmute(
        side = "ROBILA",
        x = 2,
        y = y,
        module = ROBILA_mod,
        label = paste0("ROBILA ", ROBILA_mod, "\n", ROBILA_label)
      )
  )
  
  schematic_links <- blue_children |>
    dplyr::transmute(
      x = 1.08,
      y = mean(blue_children$y),
      xend = 1.92,
      yend = y,
      Parent_fraction = Parent_fraction,
      Pair_class = Pair_class
    )
  
  p_schematic <- ggplot() +
    geom_segment(
      data = schematic_links,
      aes(x = x, y = y, xend = xend, yend = yend, linewidth = Parent_fraction),
      colour = COL_GREY,
      alpha = 0.8,
      lineend = "round"
    ) +
    geom_label(
      data = schematic_nodes,
      aes(x = x, y = y, label = label, fill = module),
      colour = "black",
      label.size = 0.25,
      label.r = grid::unit(0.15, "lines"),
      size = 3.2,
      lineheight = 0.95
    ) +
    scale_fill_identity() +
    scale_linewidth(range = c(0.8, 2.5), guide = "none") +
    scale_x_continuous(breaks = c(1, 2), labels = c("PRESTO", "ROBILA"), limits = c(0.7, 2.3)) +
    labs(
      title = "Biological refinement of the PRESTO blue backbone in ROBILA",
      subtitle = "ROBILA child modules are shown as subdivisions of the broad PRESTO-blue programme",
      x = NULL, y = NULL
    ) +
    theme_pub(10) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.border = element_blank()
    )
  
  save_plot_all(
    p_schematic,
    file.path(DIR_FIG, "06_Biological_schematic_PRESTOblue_to_ROBILAchildren"),
    w = 10, h = 5.5
  )
}

## ----- Figure 07. Conserved brown defence core + selected shared hubs -----
brown_pair <- matched_pairs |>
  dplyr::filter(PRESTO_mod == "brown" | ROBILA_mod == "brown") |>
  dplyr::arrange(desc(Jaccard_common), desc(N_overlap)) |>
  dplyr::slice_head(n = 1)

if (nrow(brown_pair) == 1) {
  bp <- brown_pair$PRESTO_mod[1]
  br <- brown_pair$ROBILA_mod[1]
  
  brown_terms <- func_details |>
    dplyr::filter(PRESTO_mod == bp, ROBILA_mod == br, Status == "Shared") |>
    dplyr::mutate(best_padj = pmin(PRESTO_padj, ROBILA_padj, na.rm = TRUE)) |>
    dplyr::arrange(best_padj) |>
    dplyr::distinct(ID, .keep_all = TRUE) |>
    dplyr::slice_head(n = 6)
  
  brown_term_line <- if (nrow(brown_terms) > 0) {
    paste(shorten_label(brown_terms$Description, width = 18), collapse = " | ")
  } else {
    "Shared defence/signalling terms"
  }
  
  brown_hubs <- shared_hubs |>
    dplyr::filter(PRESTO_mod == bp, ROBILA_mod == br) |>
    dplyr::mutate(
      Hub_class = ifelse(Shared_treatment_priority,
                         "Shared treatment-priority hub",
                         "Shared structural hub")
    ) |>
    dplyr::arrange(desc(Shared_treatment_priority), desc(mean_kME), desc(min_kME)) |>
    dplyr::slice_head(n = 15)
  
  if (nrow(brown_hubs) > 0) {
    p_brown_core <- ggplot(brown_hubs,
                           aes(x = mean_kME, y = reorder(Gene, mean_kME), colour = Hub_class)) +
      geom_segment(aes(x = 0, xend = mean_kME,
                       y = reorder(Gene, mean_kME), yend = reorder(Gene, mean_kME)),
                   linewidth = 0.6, alpha = 0.7) +
      geom_point(size = 3) +
      scale_colour_manual(values = c(
        "Shared structural hub" = COL_GREY,
        "Shared treatment-priority hub" = COL_SHARED
      )) +
      labs(
        title = paste0("Conserved brown defence/signalling core: ", bp, " ↔ ", br),
        subtitle = paste0("Shared enriched terms: ", brown_term_line),
        x = "Mean intramodular connectivity (mean kME)",
        y = NULL,
        colour = NULL
      ) +
      theme_pub(10)
    
    save_plot_all(
      p_brown_core,
      file.path(DIR_FIG, "07_Conserved_brown_defence_core_sharedHubs"),
      w = 9.5, h = 5.5
    )
  }
}




