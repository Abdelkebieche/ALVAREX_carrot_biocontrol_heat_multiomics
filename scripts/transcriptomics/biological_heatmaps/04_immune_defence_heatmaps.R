############################################################
## 13C_Immune_Defence_WGCNA_D0D2D4_FINAL.R
## Immunity/defence heatmap based ONLY on WGCNA-selected genes.
##
## - Gene selection: immune-annotated genes from Defence and
##   Temperature×Defence WGCNA candidate modules.
## - Temperature-only modules are NOT automatically called defence.
## - Expression shown: D0 + D2 + D4 for the SAME genes.
## - Exactly 3 heatmaps: PRESTO, ROBILA, COMBINED.
## - Receptor/NLR/PTI/ROS/hormone/etc. subfamilies remain in CSV.
# Repository configuration (portable; no HPC-specific paths)
if (!exists('CARROT_CONFIG_LOADED')) source('../../../config/config.R')


############################################################

helper <- "01_common_heatmap_helpers.R"
if (!exists("run_heatmap_set", mode = "function")) {
  if (!file.exists(helper)) stop("Put ", helper, " in the same working directory as this script.")
  source(helper)
}

IMMUNE_SELECTED <- MASTER |>
  dplyr::filter(
    In_WGCNA,
    Is_Immune,
    Candidate_Defence | Candidate_Temp_x_Defence
  ) |>
  dplyr::mutate(
    Defence_context = dplyr::case_when(
      Candidate_Defence & Candidate_Temp_x_Defence ~ "Defence + Temperature×Defence",
      Candidate_Temp_x_Defence ~ "Temperature×Defence",
      Candidate_Defence ~ "Defence",
      TRUE ~ "Other"
    )
  )

run_heatmap_set(
  set_name = "13C_Immune_Defence",
  set_title = "Immune/defence genes recruited in defence-associated WGCNA modules",
  selected = IMMUNE_SELECTED
)

cat("13C Immune/Defence complete. Immune submodules are retained in CSV only.\n")
