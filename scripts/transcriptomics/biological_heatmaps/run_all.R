############################################################
## 17_run_ALL_WGCNA_biological_heatmaps_FINAL.R
## Complete runner for the FINAL biological WGCNA heatmap part.
##
## 1. Build/update biological master (13_00)
## 2. Load D0+D2+D4 expression once (helper)
## 3. HSP (13A)
## 4. TF (13B)
## 5. Immunity/Defence (13C)
## 6. Functional pathways (14)
## 7. Global Temperature / Defence / Temperature×Defence panorama (16)
# Repository configuration (portable; no HPC-specific paths)
if (!exists('CARROT_CONFIG_LOADED')) source('../../../config/config.R')


############################################################

required <- c(
  "00_prepare_master_annotation.R",
  "01_common_heatmap_helpers.R",
  "02_hsp_heatmaps.R",
  "03_tf_heatmaps.R",
  "04_immune_defence_heatmaps.R",
  "05_functional_heatmaps.R",
  "06_global_themes.R"
)
missing <- required[!file.exists(required)]
if (length(missing)) stop("Missing scripts: ", paste(missing, collapse=", "))



cat("\nSTEP 1/7 — Build/update biological annotation master\n")
source("00_prepare_master_annotation.R")

cat("\nSTEP 2/7 — D0+D2+D4 expression helper\n")
source("01_common_heatmap_helpers.R")

cat("\nSTEP 3/7 — HSP\n")
source("02_hsp_heatmaps.R")

cat("\nSTEP 4/7 — TF\n")
source("03_tf_heatmaps.R")

cat("\nSTEP 5/7 — Immunity/Defence\n")
source("04_immune_defence_heatmaps.R")

cat("\nSTEP 6/7 — Functional\n")
source("05_functional_heatmaps.R")

cat("\nSTEP 7/7 — Global themes\n")
source("06_global_themes.R")

cat("\n============================================================\n")
cat("ALL FINAL WGCNA BIOLOGICAL HEATMAP SCRIPTS COMPLETED\n")
cat("============================================================\n")
