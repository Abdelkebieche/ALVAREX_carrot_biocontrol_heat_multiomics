############################################################
## 16_Global_Themes_WGCNA_D0D2D4_FINAL.R
## Global WGCNA biological panorama.
##
## Produces 3 broad themes:
##   Temperature
##   Defence
##   Temperature_x_Defence
##
## For EACH theme: exactly 3 heatmaps
##   PRESTO / ROBILA / COMBINED
## No family-by-family figures.
# Repository configuration (portable; no HPC-specific paths)
if (!exists('CARROT_CONFIG_LOADED')) source('../../../config/config.R')


############################################################

helper <- "01_common_heatmap_helpers.R"
if (!exists("run_heatmap_set", mode = "function")) {
  if (!file.exists(helper)) stop("Put ", helper, " in the same working directory as this script.")
  source(helper)
}

GLOBAL_SETS <- list(
  Temperature = MASTER |> dplyr::filter(In_WGCNA, Candidate_Temperature),
  Defence = MASTER |> dplyr::filter(In_WGCNA, Candidate_Defence),
  Temperature_x_Defence = MASTER |> dplyr::filter(In_WGCNA, Candidate_Temp_x_Defence)
)
GLOBAL_TITLES <- c(
  Temperature = "Temperature-associated WGCNA genes",
  Defence = "Defence-associated WGCNA genes",
  Temperature_x_Defence = "Temperature × defence WGCNA genes"
)

summary_list <- list()
for (nm in names(GLOBAL_SETS)) {
  summary_list[[nm]] <- run_heatmap_set(
    set_name = paste0("16_", nm),
    set_title = GLOBAL_TITLES[[nm]],
    selected = GLOBAL_SETS[[nm]]
  )
}

sm <- dplyr::bind_rows(summary_list)
if (nrow(sm)) readr::write_csv(sm, file.path(MASTER_DIR, "05_WGCNA_D0D2D4_Zscore_FINAL", "16_Global_theme_summary.csv"))
cat("16 Global themes complete.\n")
