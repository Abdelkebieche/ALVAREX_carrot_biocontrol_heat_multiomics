############################################################
## 13A_HSP_WGCNA_D0D2D4_FINAL.R
## HSP/chaperone heatmap based ONLY on WGCNA-selected genes.
##
## - Gene selection: WGCNA Temperature and Temperature×Defence candidates.
## - Expression shown: D0 + D2 + D4 for the SAME genes.
## - D0 = pre-treatment baseline only (NH, HS-7, HS-2).
## - Biological replicates are averaged before row z-score.
## - NO HSP-family-by-family heatmaps.
## - HSP family/subfamily details are exported to CSV.
## - Exactly 3 heatmaps: PRESTO, ROBILA, COMBINED.
# Repository configuration (portable; no HPC-specific paths)
if (!exists('CARROT_CONFIG_LOADED')) source('../../../config/config.R')


############################################################

helper <- "01_common_heatmap_helpers.R"
if (!exists("run_heatmap_set", mode = "function")) {
  if (!file.exists(helper)) stop("Put ", helper, " in the same working directory as this script.")
  source(helper)
}

HSP_SELECTED <- MASTER |>
  dplyr::filter(
    In_WGCNA,
    Is_HSP,
    Candidate_Temperature | Candidate_Temp_x_Defence
  ) |>
  dplyr::mutate(
    HSP_context = dplyr::case_when(
      Candidate_Temperature & Candidate_Temp_x_Defence ~ "Temperature + Temperature×Defence",
      Candidate_Temp_x_Defence ~ "Temperature×Defence",
      Candidate_Temperature ~ "Temperature",
      TRUE ~ "Other"
    )
  )

run_heatmap_set(
  set_name = "13A_HSP",
  set_title = "HSP/chaperone genes recruited in temperature-associated WGCNA modules",
  selected = HSP_SELECTED
)

cat("13A HSP complete. No family-by-family figures were generated.\n")
