############################################################
## 13B_TF_WGCNA_D0D2D4_FINAL.R
## Transcription-factor heatmap based ONLY on WGCNA-selected genes.
##
## - Includes TFs recruited by Temperature, Defence, or Temperature×Defence modules.
## - Biological context remains in CSV; no TF-family-by-family heatmaps.
## - Expression shown: D0 + D2 + D4 for the SAME WGCNA-selected TFs.
## - Exactly 3 heatmaps: PRESTO, ROBILA, COMBINED.
# Repository configuration (portable; no HPC-specific paths)
if (!exists('CARROT_CONFIG_LOADED')) source('../../../config/config.R')


############################################################

helper <- "01_common_heatmap_helpers.R"
if (!exists("run_heatmap_set", mode = "function")) {
  if (!file.exists(helper)) stop("Put ", helper, " in the same working directory as this script.")
  source(helper)
}

TF_SELECTED <- MASTER |>
  dplyr::filter(
    In_WGCNA,
    Is_TF,
    Candidate_Temperature | Candidate_Defence | Candidate_Temp_x_Defence
  ) |>
  dplyr::mutate(
    TF_context = dplyr::case_when(
      Candidate_Temp_x_Defence ~ "Temperature×Defence",
      Candidate_Defence & Candidate_Temperature ~ "Temperature + Defence",
      Candidate_Defence ~ "Defence",
      Candidate_Temperature ~ "Temperature",
      TRUE ~ "Other"
    )
  )

run_heatmap_set(
  set_name = "13B_TF",
  set_title = "Transcription factors recruited in biologically relevant WGCNA modules",
  selected = TF_SELECTED
)

cat("13B TF complete. TF families/subgroups are retained in CSV only.\n")
