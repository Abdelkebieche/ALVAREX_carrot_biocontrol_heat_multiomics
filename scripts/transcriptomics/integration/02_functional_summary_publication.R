
# Repository configuration (portable; no HPC-specific paths)
if (!exists('CARROT_CONFIG_LOADED')) source('../../../config/config.R')

# ======================================================================
# WGCNA FUNCTIONAL SUMMARY — PRESTO + ROBILA
# GO-BP / pathway summary tables + publication-style bubble matrices
#
# PURPOSE
#   1) Combine ALL WGCNA modules from PRESTO and ROBILA.
#   2) Summarise significant GO Biological Process enrichments by
#      biologically interpretable categories/subcategories.
#   3) Produce a compact main-figure matrix and a detailed term-level matrix.
#   4) Optionally combine GO-BP + enriched pathways to retain biological
#      messages that GO alone may miss (e.g. chaperones, MAPK, pathways).
#   5) Add module-level significant experimental effects when the integration
#      model output is available.
#
# IMPORTANT
#   WGCNA module enrichment has NO intrinsic "up/down" direction.
#   Therefore bubble colour = enrichment strength (-log10 FDR), NOT up/down.
#   Bubble labels = gene count for a term, or number of significant terms in
#   the summary panel. Direction should only be added later for a defined
#   contrast (e.g. Treated-Control within NH) using module eigengene contrasts.
# ======================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(forcats)
  library(scales)
})

# ----------------------------------------------------------------------
# 0. USER SETTINGS
# ----------------------------------------------------------------------

FDR_CUTOFF <- 0.05

# Number of representative terms retained per biological subcategory
# for the detailed publication-style figure.
TOP_TERMS_PER_SUBCATEGORY <- 3

# Minimum gene count for a term to be eligible for the publication matrix.
MIN_TERM_COUNT <- 3

# Set the root containing WGCNA_PRESTO_* and WGCNA_ROBILA_* directories.
# Change this path if needed.
WGCNA_ROOT <- RNA_WGCNA_ROOT

# OPTIONAL manual paths.
# Leave NA to auto-detect the latest matching directory.
PRESTO_RUN <- NA_character_
ROBILA_RUN <- NA_character_
INTEGRATION_RUN <- NA_character_

# Output folder.
OUT_DIR <- file.path(
  WGCNA_ROOT,
  paste0("Functional_GO_summary_PRESTO_ROBILA_", format(Sys.time(), "%Y%m%d_%H%M%S"))
)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
DIR_TAB <- file.path(OUT_DIR, "01_tables")
DIR_FIG <- file.path(OUT_DIR, "02_figures")
dir.create(DIR_TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_FIG, recursive = TRUE, showWarnings = FALSE)

# ----------------------------------------------------------------------
# 1. HELPERS: locate latest runs
# ----------------------------------------------------------------------

latest_matching_dir <- function(root, pattern) {
  dd <- list.dirs(root, recursive = TRUE, full.names = TRUE)
  hit <- dd[grepl(pattern, basename(dd))]
  if (length(hit) == 0) return(NA_character_)
  info <- file.info(hit)
  hit[order(info$mtime, decreasing = TRUE)][1]
}

if (is.na(PRESTO_RUN)) {
  PRESTO_RUN <- latest_matching_dir(WGCNA_ROOT, "^WGCNA_PRESTO_")
}
if (is.na(ROBILA_RUN)) {
  ROBILA_RUN <- latest_matching_dir(WGCNA_ROOT, "^WGCNA_ROBILA_")
}
if (is.na(INTEGRATION_RUN)) {
  INTEGRATION_RUN <- latest_matching_dir(WGCNA_ROOT, "^Integration_v2_")
}

if (is.na(PRESTO_RUN) || !dir.exists(PRESTO_RUN)) {
  stop("PRESTO WGCNA run not found. Set PRESTO_RUN manually.")
}
if (is.na(ROBILA_RUN) || !dir.exists(ROBILA_RUN)) {
  stop("ROBILA WGCNA run not found. Set ROBILA_RUN manually.")
}

message("PRESTO_RUN: ", PRESTO_RUN)
message("ROBILA_RUN: ", ROBILA_RUN)
message("INTEGRATION_RUN: ", INTEGRATION_RUN)
message("OUT_DIR: ", OUT_DIR)

# ----------------------------------------------------------------------
# 2. READ ENRICHMENT OUTPUTS
# ----------------------------------------------------------------------

read_go <- function(run_dir, genotype) {
  f <- file.path(run_dir, "05_enrichment", "GO_ALL_modules_combined.csv")
  if (!file.exists(f)) stop("Missing GO file: ", f)
  
  readr::read_csv(f, show_col_types = FALSE) |>
    dplyr::mutate(
      Genotype = genotype,
      Source = "GO_BP",
      p.adjust = as.numeric(p.adjust),
      Count = as.numeric(Count)
    ) |>
    dplyr::filter(
      Ontology == "BP",
      !is.na(p.adjust),
      p.adjust < FDR_CUTOFF,
      !is.na(Description),
      Count >= MIN_TERM_COUNT
    )
}

read_pathway <- function(run_dir, genotype) {
  f <- file.path(run_dir, "05_enrichment", "Pathway_ALL_modules_combined.csv")
  if (!file.exists(f)) {
    warning("Missing pathway file: ", f)
    return(tibble())
  }
  
  readr::read_csv(f, show_col_types = FALSE) |>
    dplyr::mutate(
      Genotype = genotype,
      Source = "Pathway",
      Ontology = "Pathway",
      p.adjust = as.numeric(p.adjust),
      Count = as.numeric(Count)
    ) |>
    dplyr::filter(
      !is.na(p.adjust),
      p.adjust < FDR_CUTOFF,
      !is.na(Description),
      Count >= MIN_TERM_COUNT
    )
}

GO_ALL <- dplyr::bind_rows(
  read_go(PRESTO_RUN, "PRESTO"),
  read_go(ROBILA_RUN, "ROBILA")
)

PATH_ALL <- dplyr::bind_rows(
  read_pathway(PRESTO_RUN, "PRESTO"),
  read_pathway(ROBILA_RUN, "ROBILA")
)

FUNC_ALL <- dplyr::bind_rows(GO_ALL, PATH_ALL)

# ----------------------------------------------------------------------
# 3. BIOLOGICAL CLASSIFICATION
# ----------------------------------------------------------------------
# Keyword-based classification is used ONLY to organise the figure.
# All original GO/pathway descriptions and FDR values are retained in tables.
# Terms not captured by the dictionary remain in "Other".

classify_subcategory <- function(description) {
  x <- stringr::str_to_lower(description)
  
  dplyr::case_when(
    # Defence / stress / signalling
    str_detect(x, "defen[cs]e|immune|biotic stimulus|pathogen|systemic acquired") ~
      "Defence & immunity",
    str_detect(x, "hydrogen peroxide|reactive oxygen|oxidative stress|redox") ~
      "ROS & oxidative stress",
    str_detect(x, "mapk|signal transduction|signaling|signalling|receptor|cell communication|cell recognition") ~
      "Signalling & recognition",
    str_detect(x, "auxin|abscisic|jasmon|salicyl|ethylene|hormone") ~
      "Hormone response",
    
    # Photosynthesis / energy
    str_detect(x, "photosynth|light harvesting|light reaction|pigment") ~
      "Photosynthesis & light",
    str_detect(x, "electron transport|oxidative phosphorylation|precursor metabolites and energy|carbon fixation|atp synth") ~
      "Energy & electron transport",
    str_detect(x, "mitochond") ~
      "Mitochondrial functions",
    
    # Protein synthesis / proteostasis
    str_detect(x, "translation|ribosom|peptide biosynt|amide biosynt") ~
      "Translation & ribosome",
    str_detect(x, "chaperone|protein folding|endoplasmic reticulum|protein processing|protein maturation") ~
      "Proteostasis & folding",
    str_detect(x, "ubiquitin|proteasom|deubiquit|small protein conjugation|small protein removal") ~
      "Protein turnover",
    
    # Genetic information
    str_detect(x, "rna processing|mrna|trna|spliceos|rna degrad|rna biosynt|rna metabolic") ~
      "RNA processing",
    str_detect(x, "chromatin|chromosome|dna replication|dna repair|recombination") ~
      "Chromatin & DNA",
    str_detect(x, "transcription|gene expression|nucleobase-containing compound") ~
      "Transcriptional regulation",
    
    # Metabolism
    str_detect(x, "phenylpropanoid|flavonoid|terpenoid|secondary metabol") ~
      "Secondary metabolism",
    str_detect(x, "lipid|fatty acid|glycerolipid|phospholipid|steroid|wax biosynth") ~
      "Lipid metabolism",
    str_detect(x, "carbohydrate|glycol|gluconeogenesis|pentose phosphate|starch|sucrose|sugar|calvin") ~
      "Carbon & carbohydrate metabolism",
    str_detect(x, "organic acid|oxoacid|carboxylic acid|amino acid|nitrogen") ~
      "Amino/organic-acid metabolism",
    str_detect(x, "purine|pyrimidine|nucleotide|nucleoside") ~
      "Nucleotide metabolism",
    str_detect(x, "cytochrome p450|xenobiotic|drug metabolism") ~
      "Detoxification",
    
    # Transport / cellular organisation
    str_detect(x, "transport|localization|localisation|trafficking|snare|transmembrane|ion channel") ~
      "Transport & membrane trafficking",
    str_detect(x, "organelle organization|organelle organisation|cellular component organization|cellular component organisation|microtubule|cytoskeleton") ~
      "Cellular organisation",
    
    # Development / reproduction
    str_detect(x, "pollen|pollination|reproductive|development") ~
      "Development & reproduction",
    
    TRUE ~ "Other"
  )
}

subcategory_to_category <- function(subcat) {
  dplyr::case_when(
    subcat %in% c("Defence & immunity", "ROS & oxidative stress",
                  "Signalling & recognition", "Hormone response") ~
      "Stress, defence & signalling",
    
    subcat %in% c("Photosynthesis & light", "Energy & electron transport",
                  "Mitochondrial functions") ~
      "Photosynthesis & energy",
    
    subcat %in% c("Translation & ribosome", "Proteostasis & folding",
                  "Protein turnover") ~
      "Protein synthesis & proteostasis",
    
    subcat %in% c("RNA processing", "Chromatin & DNA",
                  "Transcriptional regulation") ~
      "Genetic information processing",
    
    subcat %in% c("Secondary metabolism", "Lipid metabolism",
                  "Carbon & carbohydrate metabolism",
                  "Amino/organic-acid metabolism",
                  "Nucleotide metabolism", "Detoxification") ~
      "Metabolism",
    
    subcat %in% c("Transport & membrane trafficking", "Cellular organisation") ~
      "Transport & cellular organisation",
    
    subcat == "Development & reproduction" ~
      "Development & reproduction",
    
    TRUE ~ "Other"
  )
}

GO_ALL <- GO_ALL |>
  dplyr::mutate(
    Subcategory = classify_subcategory(Description),
    Category = subcategory_to_category(Subcategory),
    Cell = paste(Genotype, Module, sep = " | "),
    neglogFDR = -log10(pmax(p.adjust, 1e-300))
  )

FUNC_ALL <- FUNC_ALL |>
  dplyr::mutate(
    Subcategory = classify_subcategory(Description),
    Category = subcategory_to_category(Subcategory),
    Cell = paste(Genotype, Module, sep = " | "),
    neglogFDR = -log10(pmax(p.adjust, 1e-300))
  )

# ----------------------------------------------------------------------
# 4. MODULE ORDER — keep ALL modules, even if a module has no significant GO
# ----------------------------------------------------------------------

PRESTO_MODULE_ORDER <- c("brown", "blue", "turquoise")
ROBILA_MODULE_ORDER <- c("brown", "blue", "turquoise",
                         "yellow", "green", "red", "black", "pink")

all_modules <- tibble::tibble(
  Genotype = c(rep("PRESTO", length(PRESTO_MODULE_ORDER)),
               rep("ROBILA", length(ROBILA_MODULE_ORDER))),
  Module = c(PRESTO_MODULE_ORDER, ROBILA_MODULE_ORDER)
)

# ----------------------------------------------------------------------
# 5. UNIQUE-GENE COUNT HELPER
# ----------------------------------------------------------------------

unique_gene_count <- function(x) {
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0) return(0L)
  
  genes <- unique(unlist(strsplit(paste(x, collapse = "/"), "/", fixed = TRUE)))
  genes <- genes[!is.na(genes) & nzchar(genes)]
  length(genes)
}

# ----------------------------------------------------------------------
# 6. COMPLETE GO TABLE — ALL SIGNIFICANT GO-BP TERMS
# ----------------------------------------------------------------------

GO_TABLE_ALL <- GO_ALL |>
  dplyr::arrange(Category, Subcategory, Genotype, Module, p.adjust)

readr::write_csv(
  GO_TABLE_ALL,
  file.path(DIR_TAB, "00_GO_BP_allSignificant_allModules_bothGenotypes.csv")
)

# ----------------------------------------------------------------------
# 7. BIOLOGICAL SUBCATEGORY SUMMARY
#
# One row = Genotype x Module x Subcategory.
# N_GO_terms        = number of significant GO-BP terms
# Unique_gene_count = union of genes across all significant terms
# Best_FDR          = strongest enrichment in that subcategory
# Representative_GO = best-FDR term
# ----------------------------------------------------------------------

GO_SUBCAT <- GO_ALL |>
  dplyr::group_by(Genotype, Module, Category, Subcategory) |>
  dplyr::summarise(
    N_GO_terms = dplyr::n_distinct(ID),
    Unique_gene_count = unique_gene_count(geneID),
    Best_FDR = min(p.adjust, na.rm = TRUE),
    Best_neglogFDR = max(neglogFDR, na.rm = TRUE),
    Representative_GO = Description[which.min(p.adjust)][1],
    Representative_GO_ID = ID[which.min(p.adjust)][1],
    .groups = "drop"
  )

readr::write_csv(
  GO_SUBCAT,
  file.path(DIR_TAB, "01_GO_BP_biologicalSubcategory_summary_long.csv")
)

# Wide numeric matrix: unique genes per subcategory x module
GO_SUBCAT_COUNTS_WIDE <- GO_SUBCAT |>
  dplyr::mutate(Cell = paste(Genotype, Module, sep = "__")) |>
  dplyr::select(Category, Subcategory, Cell, Unique_gene_count) |>
  tidyr::pivot_wider(
    names_from = Cell,
    values_from = Unique_gene_count,
    values_fill = 0
  )

readr::write_csv(
  GO_SUBCAT_COUNTS_WIDE,
  file.path(DIR_TAB, "01_GO_BP_biologicalSubcategory_uniqueGeneCounts_wide.csv")
)

# Wide matrix: number of significant GO terms
GO_SUBCAT_NTERMS_WIDE <- GO_SUBCAT |>
  dplyr::mutate(Cell = paste(Genotype, Module, sep = "__")) |>
  dplyr::select(Category, Subcategory, Cell, N_GO_terms) |>
  tidyr::pivot_wider(
    names_from = Cell,
    values_from = N_GO_terms,
    values_fill = 0
  )

readr::write_csv(
  GO_SUBCAT_NTERMS_WIDE,
  file.path(DIR_TAB, "01_GO_BP_biologicalSubcategory_Nterms_wide.csv")
)

# ----------------------------------------------------------------------
# 8. MAIN-FIGURE FRIENDLY SUMMARY MATRIX
#
# Bubble size  = unique genes represented in significant GO terms
# Bubble fill  = best -log10(FDR)
# Number inside bubble = number of significant GO terms
# ----------------------------------------------------------------------

category_order <- c(
  "Stress, defence & signalling",
  "Photosynthesis & energy",
  "Protein synthesis & proteostasis",
  "Genetic information processing",
  "Metabolism",
  "Transport & cellular organisation",
  "Development & reproduction",
  "Other"
)

subcategory_order <- c(
  "Defence & immunity",
  "ROS & oxidative stress",
  "Signalling & recognition",
  "Hormone response",
  "Photosynthesis & light",
  "Energy & electron transport",
  "Mitochondrial functions",
  "Translation & ribosome",
  "Proteostasis & folding",
  "Protein turnover",
  "RNA processing",
  "Chromatin & DNA",
  "Transcriptional regulation",
  "Secondary metabolism",
  "Lipid metabolism",
  "Carbon & carbohydrate metabolism",
  "Amino/organic-acid metabolism",
  "Nucleotide metabolism",
  "Detoxification",
  "Transport & membrane trafficking",
  "Cellular organisation",
  "Development & reproduction",
  "Other"
)

plot_summary_data <- GO_SUBCAT |>
  dplyr::mutate(
    Category = factor(Category, levels = category_order),
    Subcategory = factor(Subcategory, levels = subcategory_order),
    Module = dplyr::case_when(
      Genotype == "PRESTO" ~ factor(Module, levels = PRESTO_MODULE_ORDER),
      TRUE ~ factor(Module, levels = ROBILA_MODULE_ORDER)
    )
  )

p_summary <- ggplot(
  plot_summary_data,
  aes(x = Module, y = forcats::fct_rev(Subcategory))
) +
  geom_point(
    aes(size = Unique_gene_count, fill = Best_neglogFDR),
    shape = 21, colour = "grey25", stroke = 0.35
  ) +
  geom_text(
    aes(label = N_GO_terms),
    size = 2.6, fontface = "bold"
  ) +
  facet_grid(
    Category ~ Genotype,
    scales = "free",
    space = "free",
    switch = "y"
  ) +
  scale_size_area(max_size = 13, name = "Unique genes") +
  scale_fill_gradient(
    low = "white", high = "firebrick",
    name = expression(-log[10]("FDR"))
  ) +
  labs(
    title = "Functional summary of WGCNA modules",
    subtitle = "GO Biological Process enrichment; number = significant GO terms",
    x = "WGCNA module",
    y = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(
    strip.background = element_rect(fill = "grey94", colour = "grey75"),
    strip.text = element_text(face = "bold"),
    strip.placement = "outside",
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(size = 8.5),
    legend.position = "right",
    panel.spacing = grid::unit(0.55, "lines")
  )

ggsave(
  file.path(DIR_FIG, "Fig_GO_BP_SubcategorySummary_bothGenotypes.png"),
  p_summary, width = 12.5, height = 10.5, dpi = 400
)
ggsave(
  file.path(DIR_FIG, "Fig_GO_BP_SubcategorySummary_bothGenotypes.pdf"),
  p_summary, width = 12.5, height = 10.5
)

if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(
    file.path(DIR_FIG, "Fig_GO_BP_SubcategorySummary_bothGenotypes.svg"),
    p_summary, width = 12.5, height = 10.5,
    device = svglite::svglite
  )
}

# ----------------------------------------------------------------------
# 9. SELECT REPRESENTATIVE GO TERMS FOR A REFERENCE-STYLE MATRIX
#
# Selection is performed across BOTH genotypes/modules:
#   1) terms enriched in more cells are prioritised
#   2) then strongest minimum FDR
#   3) then largest gene count
#
# This avoids choosing terms separately for each module and makes the
# cross-genotype visual comparison fairer.
# ----------------------------------------------------------------------

GO_TERM_RANK <- GO_ALL |>
  dplyr::group_by(Category, Subcategory, ID, Description) |>
  dplyr::summarise(
    N_cells = dplyr::n_distinct(Cell),
    Min_FDR = min(p.adjust, na.rm = TRUE),
    Max_Count = max(Count, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::arrange(Category, Subcategory, dplyr::desc(N_cells), Min_FDR, dplyr::desc(Max_Count)) |>
  dplyr::group_by(Category, Subcategory) |>
  dplyr::slice_head(n = TOP_TERMS_PER_SUBCATEGORY) |>
  dplyr::ungroup()

GO_REP <- GO_ALL |>
  dplyr::semi_join(
    GO_TERM_RANK,
    by = c("Category", "Subcategory", "ID", "Description")
  )

readr::write_csv(
  GO_TERM_RANK,
  file.path(DIR_TAB, "02_GO_BP_representativeTerms_selection.csv")
)
readr::write_csv(
  GO_REP,
  file.path(DIR_TAB, "02_GO_BP_representativeTerms_allCells.csv")
)

# Create stable y-axis order
term_order <- GO_TERM_RANK |>
  dplyr::mutate(
    Category = factor(Category, levels = category_order),
    Subcategory = factor(Subcategory, levels = subcategory_order)
  ) |>
  dplyr::arrange(Category, Subcategory, Min_FDR) |>
  dplyr::pull(Description) |>
  unique()

GO_REP <- GO_REP |>
  dplyr::mutate(
    Category = factor(Category, levels = category_order),
    Subcategory = factor(Subcategory, levels = subcategory_order),
    TermLabel = paste0(Subcategory, " — ", Description),
    TermLabel = factor(
      TermLabel,
      levels = unique(paste0(
        GO_TERM_RANK$Subcategory[match(term_order, GO_TERM_RANK$Description)],
        " — ", term_order
      ))
    )
  )

# Safer order if duplicate Description strings exist
GO_REP <- GO_REP |>
  dplyr::mutate(
    TermKey = paste(Category, Subcategory, Description, sep = "|||")
  )

term_key_order <- GO_TERM_RANK |>
  dplyr::mutate(
    Category = factor(Category, levels = category_order),
    Subcategory = factor(Subcategory, levels = subcategory_order),
    TermKey = paste(Category, Subcategory, Description, sep = "|||")
  ) |>
  dplyr::arrange(Category, Subcategory, Min_FDR) |>
  dplyr::pull(TermKey)

GO_REP <- GO_REP |>
  dplyr::mutate(
    TermKey = factor(TermKey, levels = unique(term_key_order)),
    PrettyTerm = paste0(Subcategory, " — ", Description)
  )

p_terms <- ggplot(
  GO_REP,
  aes(x = Module, y = forcats::fct_rev(TermKey))
) +
  geom_point(
    aes(size = Count, fill = neglogFDR),
    shape = 21, colour = "grey20", stroke = 0.35
  ) +
  geom_text(
    aes(label = Count),
    size = 2.2, fontface = "bold"
  ) +
  facet_grid(
    Category ~ Genotype,
    scales = "free",
    space = "free",
    switch = "y"
  ) +
  scale_y_discrete(labels = function(x) {
    # x is Category|||Subcategory|||Description
    vapply(strsplit(x, "\\|\\|\\|"), function(z) {
      paste0(z[2], " — ", z[3])
    }, character(1))
  }) +
  scale_size_area(max_size = 12, name = "Gene count") +
  scale_fill_gradient(
    low = "white", high = "firebrick",
    name = expression(-log[10]("FDR"))
  ) +
  labs(
    title = "Representative biological processes across WGCNA modules",
    subtitle = paste0(
      "Top ", TOP_TERMS_PER_SUBCATEGORY,
      " non-redundant representative GO-BP terms per biological subcategory; ",
      "number inside bubble = genes"
    ),
    x = "WGCNA module",
    y = NULL
  ) +
  theme_classic(base_size = 9) +
  theme(
    strip.background = element_rect(fill = "grey94", colour = "grey75"),
    strip.text = element_text(face = "bold"),
    strip.placement = "outside",
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(size = 7.5),
    legend.position = "right",
    panel.spacing = grid::unit(0.45, "lines")
  )

ggsave(
  file.path(DIR_FIG, "Fig_GO_BP_RepresentativeTerms_bothGenotypes.png"),
  p_terms, width = 14.5, height = 14, dpi = 400
)
ggsave(
  file.path(DIR_FIG, "Fig_GO_BP_RepresentativeTerms_bothGenotypes.pdf"),
  p_terms, width = 14.5, height = 14
)

if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(
    file.path(DIR_FIG, "Fig_GO_BP_RepresentativeTerms_bothGenotypes.svg"),
    p_terms, width = 14.5, height = 14,
    device = svglite::svglite
  )
}

# ----------------------------------------------------------------------
# 10. OPTIONAL: COMBINED GO-BP + PATHWAY BIOLOGICAL MATRIX
#
# This is often more useful for the MAIN article because some modules
# have a strong pathway signal but few/no significant GO-BP terms.
# Example: chaperones/protein folding in ROBILA green.
# ----------------------------------------------------------------------

FUNC_TERM_RANK <- FUNC_ALL |>
  dplyr::group_by(Source, Category, Subcategory, ID, Description) |>
  dplyr::summarise(
    N_cells = dplyr::n_distinct(Cell),
    Min_FDR = min(p.adjust, na.rm = TRUE),
    Max_Count = max(Count, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::arrange(
    Category, Subcategory,
    dplyr::desc(N_cells), Min_FDR, dplyr::desc(Max_Count)
  ) |>
  dplyr::group_by(Category, Subcategory, Source) |>
  dplyr::slice_head(n = 2) |>
  dplyr::ungroup()

FUNC_REP <- FUNC_ALL |>
  dplyr::semi_join(
    FUNC_TERM_RANK,
    by = c("Source", "Category", "Subcategory", "ID", "Description")
  ) |>
  dplyr::mutate(
    Category = factor(Category, levels = category_order),
    Subcategory = factor(Subcategory, levels = subcategory_order),
    TermKey = paste(Category, Subcategory, Source, Description, sep = "|||")
  )

func_key_order <- FUNC_TERM_RANK |>
  dplyr::mutate(
    Category = factor(Category, levels = category_order),
    Subcategory = factor(Subcategory, levels = subcategory_order),
    TermKey = paste(Category, Subcategory, Source, Description, sep = "|||")
  ) |>
  dplyr::arrange(Category, Subcategory, Source, Min_FDR) |>
  dplyr::pull(TermKey)

FUNC_REP <- FUNC_REP |>
  dplyr::mutate(TermKey = factor(TermKey, levels = unique(func_key_order)))

readr::write_csv(
  FUNC_TERM_RANK,
  file.path(DIR_TAB, "03_Functional_GOplusPathway_representativeTerms_selection.csv")
)
readr::write_csv(
  FUNC_REP,
  file.path(DIR_TAB, "03_Functional_GOplusPathway_representativeTerms_allCells.csv")
)

p_func <- ggplot(
  FUNC_REP,
  aes(x = Module, y = forcats::fct_rev(TermKey))
) +
  geom_point(
    aes(size = Count, fill = neglogFDR),
    shape = 21, colour = "grey20", stroke = 0.35
  ) +
  geom_text(
    aes(label = Count),
    size = 2.05, fontface = "bold"
  ) +
  facet_grid(
    Category ~ Genotype,
    scales = "free",
    space = "free",
    switch = "y"
  ) +
  scale_y_discrete(labels = function(x) {
    vapply(strsplit(x, "\\|\\|\\|"), function(z) {
      src <- ifelse(z[3] == "Pathway", "[Path]", "[GO]")
      paste0(src, " ", z[2], " — ", z[4])
    }, character(1))
  }) +
  scale_size_area(max_size = 12, name = "Gene count") +
  scale_fill_gradient(
    low = "white", high = "firebrick",
    name = expression(-log[10]("FDR"))
  ) +
  labs(
    title = "Biological functions associated with PRESTO and ROBILA WGCNA modules",
    subtitle = "Representative GO Biological Process and pathway enrichments; number inside bubble = genes",
    x = "WGCNA module",
    y = NULL
  ) +
  theme_classic(base_size = 9) +
  theme(
    strip.background = element_rect(fill = "grey94", colour = "grey75"),
    strip.text = element_text(face = "bold"),
    strip.placement = "outside",
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(size = 7.2),
    legend.position = "right",
    panel.spacing = grid::unit(0.45, "lines")
  )

ggsave(
  file.path(DIR_FIG, "Fig_Functional_GOplusPathway_RepresentativeTerms.png"),
  p_func, width = 15.5, height = 15, dpi = 400
)
ggsave(
  file.path(DIR_FIG, "Fig_Functional_GOplusPathway_RepresentativeTerms.pdf"),
  p_func, width = 15.5, height = 15
)

if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(
    file.path(DIR_FIG, "Fig_Functional_GOplusPathway_RepresentativeTerms.svg"),
    p_func, width = 15.5, height = 15,
    device = svglite::svglite
  )
}

# ----------------------------------------------------------------------
# 11. MODULE-LEVEL BIOLOGICAL SUMMARY TABLE
#
# One row per genotype/module:
#   - number of significant GO terms
#   - unique GO genes
#   - dominant biological subcategories
#   - top GO terms
#   - top pathways
#   - significant experimental effects (if integration output is found)
# ----------------------------------------------------------------------

top_concat <- function(x, score, n = 5) {
  oo <- order(score, decreasing = FALSE, na.last = NA)
  x <- unique(x[oo])
  x <- x[!is.na(x) & nzchar(x)]
  paste(head(x, n), collapse = "; ")
}

MOD_GO <- GO_ALL |>
  dplyr::group_by(Genotype, Module) |>
  dplyr::summarise(
    N_sig_GO_BP = dplyr::n_distinct(ID),
    Unique_GO_genes = unique_gene_count(geneID),
    Dominant_subcategories = paste(
      names(sort(table(Subcategory), decreasing = TRUE))[1:min(4, length(unique(Subcategory)))],
      collapse = "; "
    ),
    Top_GO_terms = top_concat(Description, p.adjust, n = 5),
    Best_GO_FDR = min(p.adjust, na.rm = TRUE),
    .groups = "drop"
  )

if (nrow(PATH_ALL) > 0) {
  PATH_CLASS <- PATH_ALL |>
    dplyr::mutate(
      Subcategory = classify_subcategory(Description),
      Category = subcategory_to_category(Subcategory)
    )
  
  MOD_PATH <- PATH_CLASS |>
    dplyr::group_by(Genotype, Module) |>
    dplyr::summarise(
      N_sig_Pathways = dplyr::n_distinct(ID),
      Top_pathways = top_concat(Description, p.adjust, n = 5),
      Best_Pathway_FDR = min(p.adjust, na.rm = TRUE),
      .groups = "drop"
    )
} else {
  MOD_PATH <- tibble()
}

MODULE_SUMMARY <- all_modules |>
  dplyr::left_join(MOD_GO, by = c("Genotype", "Module")) |>
  dplyr::left_join(MOD_PATH, by = c("Genotype", "Module"))

# Optional: add significant factors from integration model
if (!is.na(INTEGRATION_RUN)) {
  effect_file <- file.path(
    INTEGRATION_RUN,
    "01_tables",
    "03_ME_TypeIII_full_Temp_Treatment_Time.csv"
  )
  
  if (file.exists(effect_file)) {
    EFF <- readr::read_csv(effect_file, show_col_types = FALSE) |>
      dplyr::mutate(
        Effect_clean = dplyr::recode(
          Effect,
          "Sampling_time" = "Sampling stage",
          "Temperature" = "Thermal history",
          "Temperature:Treatment" = "Thermal history × Treatment",
          "Temperature:Sampling_time" = "Thermal history × Sampling stage",
          "Treatment:Sampling_time" = "Treatment × Sampling stage",
          "Temperature:Treatment:Sampling_time" =
            "Thermal history × Treatment × Sampling stage"
        )
      ) |>
      dplyr::group_by(Genotype, ModuleColor) |>
      dplyr::summarise(
        Significant_effects = paste(
          Effect_clean[!is.na(FDR) & FDR < FDR_CUTOFF],
          collapse = "; "
        ),
        .groups = "drop"
      ) |>
      dplyr::rename(Module = ModuleColor)
    
    MODULE_SUMMARY <- MODULE_SUMMARY |>
      dplyr::left_join(EFF, by = c("Genotype", "Module"))
  }
}

MODULE_SUMMARY <- MODULE_SUMMARY |>
  dplyr::mutate(
    dplyr::across(
      c(N_sig_GO_BP, Unique_GO_genes, N_sig_Pathways),
      ~ tidyr::replace_na(.x, 0)
    ),
    Significant_effects = tidyr::replace_na(Significant_effects, "None at FDR < 0.05"),
    Top_GO_terms = tidyr::replace_na(Top_GO_terms, "No significant GO-BP enrichment"),
    Top_pathways = tidyr::replace_na(Top_pathways, "No significant pathway enrichment"),
    Dominant_subcategories = tidyr::replace_na(
      Dominant_subcategories, "No significant GO-BP enrichment"
    )
  )

readr::write_csv(
  MODULE_SUMMARY,
  file.path(DIR_TAB, "04_WGCNA_Module_BiologicalSummary_PRESTO_ROBILA.csv")
)

# ----------------------------------------------------------------------
# 12. README / METHOD NOTE
# ----------------------------------------------------------------------

readme <- c(
  "WGCNA functional summary — PRESTO + ROBILA",
  "",
  paste0("FDR threshold: ", FDR_CUTOFF),
  paste0("Minimum genes per displayed term: ", MIN_TERM_COUNT),
  paste0("Representative GO terms per subcategory: ", TOP_TERMS_PER_SUBCATEGORY),
  "",
  "Main outputs:",
  "01_tables/00_GO_BP_allSignificant_allModules_bothGenotypes.csv",
  "01_tables/01_GO_BP_biologicalSubcategory_summary_long.csv",
  "01_tables/01_GO_BP_biologicalSubcategory_uniqueGeneCounts_wide.csv",
  "01_tables/01_GO_BP_biologicalSubcategory_Nterms_wide.csv",
  "01_tables/02_GO_BP_representativeTerms_selection.csv",
  "01_tables/03_Functional_GOplusPathway_representativeTerms_selection.csv",
  "01_tables/04_WGCNA_Module_BiologicalSummary_PRESTO_ROBILA.csv",
  "",
  "Figures:",
  "02_figures/Fig_GO_BP_SubcategorySummary_bothGenotypes.*",
  "02_figures/Fig_GO_BP_RepresentativeTerms_bothGenotypes.*",
  "02_figures/Fig_Functional_GOplusPathway_RepresentativeTerms.*",
  "",
  "Interpretation:",
  "Bubble size = number of genes represented by the enrichment.",
  "Bubble colour = enrichment strength (-log10 FDR).",
  "Number inside the bubble = gene Count for term-level figures;",
  "for the subcategory summary it is the number of significant GO terms.",
  "",
  "IMPORTANT: WGCNA module enrichment is not intrinsically up/down.",
  "Red/blue direction must be defined from an explicit module-eigengene contrast",
  "(e.g. Treated-Control within a given thermal history), not from GO enrichment itself."
)

writeLines(readme, file.path(OUT_DIR, "README_GO_functional_summary.txt"))

message("\nDONE.")
message("Results written to: ", OUT_DIR)
