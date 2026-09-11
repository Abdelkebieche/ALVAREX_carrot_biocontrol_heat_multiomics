############################################################
## 14_Functional_WGCNA_D0D2D4_FINAL.R
## Functional-gene heatmap based ONLY on WGCNA-selected genes.
##
## - Uses the functional categories inherited from the former script 10.
## - Includes functional genes recruited by Temperature, Defence, or
##   Temperature×Defence WGCNA candidate modules.
## - No pathway/family-by-family heatmaps: details remain in CSV.
## - Expression shown: D0 + D2 + D4 for the SAME genes.
## - Exactly 3 heatmaps: PRESTO, ROBILA, COMBINED.
# Repository configuration (portable; no HPC-specific paths)
if (!exists('CARROT_CONFIG_LOADED')) source('../../../config/config.R')


############################################################

helper <- "01_common_heatmap_helpers.R"
if (!exists("run_heatmap_set", mode = "function")) {
  if (!file.exists(helper)) stop("Put ", helper, " in the same working directory as this script.")
  source(helper)
}

FUNCTIONAL_SELECTED <- MASTER |>
  dplyr::filter(
    In_WGCNA,
    Is_FunctionalSet,
    Candidate_Temperature | Candidate_Defence | Candidate_Temp_x_Defence
  ) |>
  dplyr::mutate(
    Functional_context = dplyr::case_when(
      Candidate_Temp_x_Defence ~ "Temperature×Defence",
      Candidate_Defence & Candidate_Temperature ~ "Temperature + Defence",
      Candidate_Defence ~ "Defence",
      Candidate_Temperature ~ "Temperature",
      TRUE ~ "Other"
    )
  )

run_heatmap_set(
  set_name = "14_Functional",
  set_title = "Functional pathways represented in biologically relevant WGCNA modules",
  selected = FUNCTIONAL_SELECTED
)

cat("14 Functional complete. Functional families/pathways are retained in CSV only.\n")
