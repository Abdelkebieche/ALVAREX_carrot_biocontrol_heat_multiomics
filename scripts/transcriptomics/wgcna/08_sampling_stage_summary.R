############################################################
## 08_sampling_effect.R — SAMPLING EFFECT (D2 vs D4)
## ────────────────────────────────────────────────────────
## Dedicated script for:
##   "Which modules change between P2 and P3?"
##   "Does the SDP effect persist at P3?"
##   "Does heat stress effect increase or decrease?"
##
## Analyses:
##   1. Paired-like analysis: ΔME = ME(P3) - ME(P2) per condition
##   2. Heatmap ME × (D2 vs D4) within each Temp×Trt
##   3. Temporal profiles D2→D4 with statistical tests
##   4. Temporal classification: transient, persistent, late
##   5. Time × Temperature interaction
##   6. Time × Treatment (SDP) interaction
##   7. Genes changing between P2 and P3 (module-level)
# Repository configuration (portable; no HPC-specific paths)
if (!exists('CARROT_CONFIG_LOADED')) source('../../../config/config.R')


############################################################

##############################
## PARAMETERS
##############################
INSTALL_MISSING_PKGS <- FALSE
FDR_CUTOFF <- 0.05

TEMP_COLORS <- c("T0" = "#313695", "T1" = "#74ADD1", "T2" = "#F46D43")
TRT_COLORS  <- c("EAU" = "#4575B4", "SDP" = "#D73027")
TIME_COLORS <- c("P2" = "#56B4E9", "P3" = "#009E73")

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
me_map   <- read.csv(file.path(run_dir, "02_network", "ME_color_map.csv"))
gene_mod <- read.csv(file.path(run_dir, "02_network", "Gene_to_Module.csv"))

GENO <- as.character(meta$Genotype[1])
stopifnot(rownames(MEs) == meta$Sample)

MEs_clean <- MEs[, !grepl("^ME0$", colnames(MEs)), drop = FALSE]
me_to_color <- setNames(me_map$ModuleColor, me_map$ME)

dir_time <- file.path(run_dir, "07_sampling_effect")
dir_time_box   <- file.path(dir_time, "boxplots")
dir_time_prof  <- file.path(dir_time, "profiles")
dir_time_delta <- file.path(dir_time, "delta_ME")
for (d in c(dir_time, dir_time_box, dir_time_prof, dir_time_delta))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

log_msg <- init_logger(file.path(run_dir, "pipeline.log"))
log_msg("START Script 08 — Sampling effect (D2 vs D4)")

cat("════════════════════════════════════════════════════\n")
cat("  Script 08: Sampling effect (D2 vs D4)\n")
cat("  Genotype:", GENO, "\n")
cat("  P2 samples:", sum(meta$Sampling_time == "P2"), "\n")
cat("  P3 samples:", sum(meta$Sampling_time == "P3"), "\n")
cat("════════════════════════════════════════════════════\n\n")

##############################################################
## PART 1: GLOBAL TEST — D2 vs D4 per module
##            within each sub-condition (Temp × Trt)
##############################################################
cat("── PART 1: D2 vs D4 tests per condition ──\n")

# Build all combinations to test
conditions <- expand.grid(
  Temperature = levels(meta$Temperature),
  Treatment   = levels(meta$Treatment),
  stringsAsFactors = FALSE
) |>
  dplyr::mutate(Condition = paste(Temperature, Treatment, sep="_"))

# Add global tests (all Temp, all Trt)
conditions_extended <- dplyr::bind_rows(
  conditions |> dplyr::mutate(Scope = "within_TempTrt"),
  data.frame(Temperature = "ALL", Treatment = levels(meta$Treatment),
             Condition = paste0("ALL_", levels(meta$Treatment)), Scope = "within_Trt"),
  data.frame(Temperature = levels(meta$Temperature), Treatment = "ALL",
             Condition = paste0(levels(meta$Temperature), "_ALL"), Scope = "within_Temp"),
  data.frame(Temperature = "ALL", Treatment = "ALL",
             Condition = "GLOBAL", Scope = "global")
)

# Wilcoxon tests D2 vs D4 for each ME × condition
time_tests <- data.frame()

for (me in colnames(MEs_clean)) {
  colr <- as.character(me_to_color[me])
  
  for (ci in seq_len(nrow(conditions_extended))) {
    temp_f <- conditions_extended$Temperature[ci]
    trt_f  <- conditions_extended$Treatment[ci]
    cond   <- conditions_extended$Condition[ci]
    scope  <- conditions_extended$Scope[ci]
    
    # Filter samples
    idx <- rep(TRUE, nrow(meta))
    if (temp_f != "ALL") idx <- idx & meta$Temperature == temp_f
    if (trt_f != "ALL")  idx <- idx & meta$Treatment == trt_f
    
    df_sub <- data.frame(
      ME   = as.numeric(MEs_clean[idx, me]),
      Time = meta$Sampling_time[idx]
    )
    
    n_p2 <- sum(df_sub$Time == "P2")
    n_p3 <- sum(df_sub$Time == "P3")
    
    if (n_p2 < 2 || n_p3 < 2) next
    
    # Wilcoxon
    pv <- tryCatch(
      wilcox.test(ME ~ Time, data=df_sub, exact=FALSE)$p.value,
      error=function(e) NA_real_
    )
    
    # Means and delta
    mean_p2 <- mean(df_sub$ME[df_sub$Time == "P2"], na.rm = TRUE)
    mean_p3 <- mean(df_sub$ME[df_sub$Time == "P3"], na.rm = TRUE)
    delta <- mean_p3 - mean_p2
    
    # Effect size (rank-biserial r)
    z_stat <- tryCatch({
      wt <- wilcox.test(ME ~ Time, data=df_sub, exact=FALSE)
      qnorm(wt$p.value / 2, lower.tail=FALSE) * sign(delta)
    }, error=function(e) NA_real_)
    
    r_effect <- z_stat / sqrt(n_p2 + n_p3)
    
    time_tests <- dplyr::bind_rows(time_tests, data.frame(
      ME = me, ModuleColor = colr, Condition = cond, Scope = scope,
      Temperature = temp_f, Treatment = trt_f,
      n_P2 = n_p2, n_P3 = n_p3,
      mean_P2 = mean_p2, mean_P3 = mean_p3,
      delta_ME = delta, p_value = pv, r_effect = r_effect,
      stringsAsFactors = FALSE
    ))
  }
}

# FDR per scope
time_tests <- time_tests |>
  dplyr::group_by(Scope) |>
  dplyr::mutate(FDR = p.adjust(p_value, "BH")) |>
  dplyr::ungroup() |>
  dplyr::arrange(FDR)

write.csv(time_tests, file.path(dir_time, "P2vsP3_all_tests.csv"), row.names=FALSE)

cat("  Tests computed:", nrow(time_tests), "\n")
cat("  Significant (FDR<0.05):", sum(time_tests$FDR < FDR_CUTOFF, na.rm=TRUE), "\n\n")

for (sc in unique(time_tests$Scope)) {
  sub <- time_tests |> dplyr::filter(Scope==sc, FDR<FDR_CUTOFF)
  cat("  ", sc, ": ", nrow(sub), " significant module×condition\n")
}
cat("\n")