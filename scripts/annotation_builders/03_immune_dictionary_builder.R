# NOTE: companion annotation-builder reference implementation.
# Paths have been sanitized. See docs/05_functional_dictionaries.md before running.

############################################################
## 09C_Immune_receptors_microbe_heatstress_panorama.R
##
## Figure 9C — Heat-stress profile of immune receptors and
## plant–microbe / plant–pathogen interaction genes in carrot
##
## BIOLOGICAL QUESTION
## ------------------------------------------------------------
## How does heat stress reprogramme genes involved in:
##   1) cell-surface immune perception,
##   2) receptor-proximal PTI signalling,
##   3) Ca2+ / ROS / MAPK early immune outputs,
##   4) intracellular NLR / ETI machinery,
##   5) SA/NHP and defence-hormone signalling,
##   6) downstream defence outputs,
## without pathogen inoculation?
##
## IMPORTANT INTERPRETIVE RULE
## ------------------------------------------------------------
## This analysis asks whether heat ALONE changes the basal
## transcriptional state / "immune competence" of genes that are
## known or candidate components of plant–microbe interactions.
##
## It does NOT demonstrate:
##   - receptor activation,
##   - ligand perception,
##   - PTI/ETI activation,
##   - pathogen resistance,
## unless supported by an independent pathogen/elicitor experiment.
##
## FIXED DISPLAY ORDER
## ------------------------------------------------------------
## ROBILA NH -> ROBILA HS-7 -> ROBILA HS-2 ->
## PRESTO NH -> PRESTO HS-7 -> PRESTO HS-2
##
## FOUR INTENDED DE CONTRASTS ONLY
## ------------------------------------------------------------
## ROBILA HS-7 vs NH
## ROBILA HS-2 vs NH
## PRESTO HS-7 vs NH
## PRESTO HS-2 vs NH
##
## SCIENTIFIC BASIS
## ------------------------------------------------------------
## CARROT / Alternaria dauci:
## Liu et al. 2024, The Plant Journal 120:1643-1661
## DOI: 10.1111/tpj.17049
##
## The DH13M14 T2T study identified:
##   - 283 LRR-RLKs
##   - 10 LysM-RLKs
##   - 202 NLRs: 79 CNL, 34 TNL, 6 RNL, 83 other
## and reported infection-responsive RLK/NLR/PTI/ETI components.
##
## Their PTI/ETI transcriptomic analysis highlighted, among others:
##   CNGCs, CDPK/CPKs, CaM/CML,
##   RBOH, OXI1, NDPK2,
##   WRKY22, ACS6, PR1,
##   EDS1, PBS1, RIN4, RPS2, SGT1, HSP90.
##
## CELL-SURFACE IMMUNITY:
##   FLS2–BAK1:
##     Sun et al. 2013, Science 342:624-628
##     DOI: 10.1126/science.1243825
##
##   LYK5–CERK1 chitin perception:
##     Cao et al. 2014, eLife 3:e03766
##     DOI: 10.7554/eLife.03766
##
##   BAK1 as common PTI co-receptor:
##     Chinchilla et al./Heese et al. era; key BAK1 immune study:
##     DOI / PMID represented in bibliography table below.
##
## HEAT x IMMUNITY INTERSECTION:
##   Temporary heat stress suppresses PAMP-triggered immunity:
##     Mühlenbeck et al. 2019, Molecular Plant Pathology
##     PMID: 30924595
##
##   Elevated temperature suppresses SA immunity via
##   GBPL3 -> CBP60g/SARD1:
##     Kim et al. 2022, Nature 607:339-346
##     DOI: 10.1038/s41586-022-04902-y
##
##   Warm temperature suppresses SAR/NHP:
##     2025 study, PMID: 40758777
##
##   BAM1–PBS1 membrane complex transduces heat to RBOHD:
##     2025, Molecular Plant
##     PMID: 41174877
##     DOI: 10.1016/j.molp.2025.10.021
##
##   FERONIA plasma-membrane nanoclusters control thermotolerance:
##     2026, Science
##     PMID: 42166587
##     DOI: 10.1126/science.aeb1752
##
##   FER also scaffolds FLS2/EFR–BAK1 immune receptor complexes:
##     Stegmann et al. 2017, Science
##     DOI: 10.1126/science.aal2541
##
## CROP EXAMPLES OF DIRECT HEAT–BIOTIC CROSSTALK:
##   TaSERK1 — high-temperature resistance to stripe rust:
##     Shi et al. 2023, Phytopathology
##     DOI: 10.1094/PHYTO-11-22-0429-R
##
##   TaCRK10 — high-temperature resistance to stripe rust:
##     DOI: 10.1111/tpj.15513
##
##   Rice PWL1 G-type LecRLK — heat tolerance / bacterial resistance tradeoff:
##     DOI: 10.1111/pbi.14150
##
##   Rice CRK cluster — heat and pathogen stress:
##     DOI: 10.1111/pbi.14381
##
##   CaMAPK1 — heat tolerance + Ralstonia resistance:
##     DOI: 10.3390/plants13131775
##
## METHODOLOGICAL PRINCIPLE
## ------------------------------------------------------------
## We distinguish:
##
## HIGH confidence:
##   explicit known immune/microbe-interaction gene names or
##   explicit immune-function annotation.
##
## MEDIUM confidence:
##   structurally/functionally relevant families in which many
##   members participate in immunity (e.g. LRR-RLK, LysM-RLK,
##   CRK, LecRLK, WAK, NLR, RLCK), but where the individual
##   carrot gene is not functionally validated.
##
## LOW confidence:
##   generic "receptor-like kinase" annotations without sufficient
##   structural or immune-specific information.
##
## LOW-confidence generic RLKs are exported for curation but are
## excluded from the main immune-network figure by default.
##
## OUTPUTS
## ------------------------------------------------------------
## 00_tables/
##   Immune_gene_catalogue_complete.csv
##   Immune_gene_catalogue_main.csv
##   Immune_DE_selected.csv
##   Immune_module_census.csv
##   Immune_module_overrepresentation_Fisher.csv
##   Immune_response_patterns.csv
##   Immune_carrot_Liu2024_reference_audit.csv
##   Immune_literature_rules.csv
##
## 01_global_heatmaps/
##   Figure_9C_GLOBAL_ALL_IMMUNE_GENES
##   Figure_9C_GLOBAL_CELL_SURFACE_RECEPTORS
##   Figure_9C_GLOBAL_DIRECT_HEAT_IMMUNITY_INTERSECTION
##
## 02_DE_heatmaps/
##   Figure_9C_COMPLETE_DE_IMMUNE_NETWORK
##   Figure_9C_DE_CELL_SURFACE_RECEPTORS
##   Figure_9C_DE_NLR_ETI
##   Figure_9C_DE_EARLY_PTI
##
## 03_module_heatmaps/
##   one DE heatmap per immune module
##
## 04_QC/
##   mapping, annotation, rejected/ambiguous candidates, RDS,
##   sessionInfo
############################################################


##############################
## 0) Configuration
##############################
source("00_builder_config.R")
source("00_builder_helpers.R")

IMMUNE_ANNOT_FILE <- Sys.getenv(
  "IMMUNE_ANNOT_FILE",
  unset = DH13_ANNOT_FILE
)

# Exact four heat-vs-NH contrasts.
IMMUNE_CONTRAST_ORDER <- c(
  "ROBILA_T1_vs_T0",
  "ROBILA_T2_vs_T0",
  "PRESTO_T1_vs_T0",
  "PRESTO_T2_vs_T0"
)

IMMUNE_CONTRAST_LABELS <- c(
  "ROBILA HS-7",
  "ROBILA HS-2",
  "PRESTO HS-7",
  "PRESTO HS-2"
)

# Fixed expression-column order.
IMMUNE_GROUP_ORDER <- c(
  "ROBILA_T0", "ROBILA_T1", "ROBILA_T2",
  "PRESTO_T0", "PRESTO_T1", "PRESTO_T2"
)

IMMUNE_GROUP_LABELS <- c(
  "ROBILA NH", "ROBILA HS-7", "ROBILA HS-2",
  "PRESTO NH", "PRESTO HS-7", "PRESTO HS-2"
)

# DE definition.
IMMUNE_FDR_CUTOFF <- PADJ_CUTOFF
IMMUNE_LFC_CUTOFF <- lfc_threshold

# Inclusion policy.
INCLUDE_LOW_CONFIDENCE_GENERIC_RLK_IN_MAIN <- FALSE

# Keep TF/HSP immune nodes in the broad immune panorama.
# They can be excluded if one wants strictly non-TF/non-HSP modules.
INCLUDE_IMMUNE_TF_HSP_NODES <- TRUE

# Heatmap settings.
IMMUNE_Z_LIMIT       <- 2.5
IMMUNE_GLOBAL_Z_LIMIT <- 2.0
IMMUNE_CLUSTER_ROWS  <- TRUE
IMMUNE_ROW_MM        <- 3.5
IMMUNE_FIG_WIDTH_MM  <- 235
IMMUNE_GLOBAL_WIDTH_MM <- 190
IMMUNE_RASTER_DPI    <- 450
IMMUNE_WRITE_PDF     <- TRUE
IMMUNE_WRITE_SVG     <- TRUE
IMMUNE_WRITE_PNG     <- TRUE
IMMUNE_WRITE_TIFF    <- TRUE
IMMUNE_RASTER_MAX_HEIGHT_MM <- 650

# Main module order reflects the signalling architecture.
IMMUNE_MODULE_ORDER <- c(
  "01_Cell-surface immune receptors",
  "02_Co-receptors / receptor scaffolds",
  "03_Receptor-proximal RLCKs",
  "04_Ca2+ signalling",
  "05_ROS / redox signalling",
  "06_MAPK signalling",
  "07_Intracellular NLR receptors",
  "08_ETI support / guardee machinery",
  "09_SA / NHP / systemic immunity",
  "10_JA / ethylene defence crosstalk",
  "11_Immune transcriptional regulators",
  "12_Defence outputs / PR / cell wall",
  "13_Beneficial-microbe / symbiosis perception",
  "14_Other immune-associated candidate"
)

# Published carrot reference counts from Liu et al. 2024.
# These are an AUDIT reference, not hard expected counts for the local
# annotation-text classifier.
LIU2024_REFERENCE_COUNTS <- tibble::tribble(
  ~Reference_class, ~Published_n,
  "LRR-RLK", 283L,
  "LysM-RLK", 10L,
  "NLR_total", 202L,
  "CNL", 79L,
  "TNL", 34L,
  "RNL", 6L,
  "Other_NLR", 83L
)

set.seed(SEED)
options(stringsAsFactors = FALSE)


##############################
## 1) Packages
##############################
load_pkgs(
  c("edgeR", "ComplexHeatmap"),
  install_missing = INSTALL_MISSING_PKGS,
  bioc = TRUE
)

load_pkgs(
  c(
    "dplyr", "tidyr", "tibble", "readr", "stringr",
    "circlize", "svglite"
  ),
  install_missing = INSTALL_MISSING_PKGS,
  bioc = FALSE
)


##############################
## 2) Utility functions
##############################
`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || all(is.na(x))) y else x
}

assert_columns <- function(x, required, object_name) {
  missing_cols <- setdiff(required, colnames(x))
  if (length(missing_cols) > 0L) {
    stop(
      object_name,
      " is missing required column(s): ",
      paste(missing_cols, collapse = ", ")
    )
  }
}

read_csv_flexible <- function(path) {
  if (!file.exists(path)) {
    stop("Input file not found: ", path)
  }
  
  x <- tryCatch(
    readr::read_csv(
      path,
      show_col_types = FALSE,
      progress = FALSE
    ),
    error = function(e) NULL
  )
  
  if (is.null(x) || ncol(x) <= 1L) {
    x <- readr::read_delim(
      path,
      delim = ";",
      show_col_types = FALSE,
      progress = FALSE
    )
  }
  
  as.data.frame(
    x,
    stringsAsFactors = FALSE
  )
}

zscore_rows_safe <- function(mat) {
  mat <- as.matrix(mat)
  
  mu <- rowMeans(
    mat,
    na.rm = TRUE
  )
  
  sds <- apply(
    mat,
    1L,
    stats::sd,
    na.rm = TRUE
  )
  
  out <- sweep(
    mat,
    1L,
    mu,
    "-"
  )
  
  valid <- is.finite(sds) & sds > 0
  
  out[valid, ] <- sweep(
    out[valid, , drop = FALSE],
    1L,
    sds[valid],
    "/"
  )
  
  out[!valid, ] <- 0
  out[!is.finite(out)] <- 0
  
  out
}

safe_file_stub <- function(x) {
  x <- iconv(
    x,
    to = "ASCII//TRANSLIT"
  )
  
  x <- gsub(
    "[^A-Za-z0-9]+",
    "_",
    x
  )
  
  x <- gsub(
    "^_+|_+$",
    "",
    x
  )
  
  x
}

clean_label <- function(
    x,
    max_chars = 68L) {
  
  x <- dplyr::coalesce(
    as.character(x),
    ""
  )
  
  x <- trimws(
    gsub(
      "\\s+",
      " ",
      x
    )
  )
  
  x <- gsub(
    "[.;]+$",
    "",
    x
  )
  
  too_long <- nchar(x) > max_chars
  
  x[too_long] <- paste0(
    substr(
      x[too_long],
      1L,
      max_chars - 3L
    ),
    "..."
  )
  
  x
}

parse_logical_override <- function(x) {
  y <- trimws(
    tolower(
      as.character(x)
    )
  )
  
  out <- rep(
    NA,
    length(y)
  )
  
  out[
    y %in% c(
      "true", "t", "1",
      "yes", "y", "oui"
    )
  ] <- TRUE
  
  out[
    y %in% c(
      "false", "f", "0",
      "no", "n", "non"
    )
  ] <- FALSE
  
  out
}

detect_any <- function(
    text,
    pattern) {
  
  stringr::str_detect(
    dplyr::coalesce(
      as.character(text),
      ""
    ),
    stringr::regex(
      pattern,
      ignore_case = TRUE
    )
  )
}

status_from_de <- function(
    logFC,
    FDR) {
  
  dplyr::case_when(
    is.finite(FDR) &
      FDR <= IMMUNE_FDR_CUTOFF &
      logFC >= IMMUNE_LFC_CUTOFF ~ "Up",
    
    is.finite(FDR) &
      FDR <= IMMUNE_FDR_CUTOFF &
      logFC <= -IMMUNE_LFC_CUTOFF ~ "Down",
    
    TRUE ~ "NS"
  )
}

response_pattern_two <- function(
    status_hs7,
    status_hs2) {
  
  dplyr::case_when(
    status_hs7 == "NS" &
      status_hs2 == "NS" ~ "None",
    
    status_hs7 != "NS" &
      status_hs2 == "NS" ~ "HS-7 only",
    
    status_hs7 == "NS" &
      status_hs2 != "NS" ~ "HS-2 only",
    
    status_hs7 == status_hs2 &
      status_hs7 != "NS" ~ paste0(
        "Both ",
        status_hs7
      ),
    
    status_hs7 != "NS" &
      status_hs2 != "NS" &
      status_hs7 != status_hs2 ~
      "Opposite directions",
    
    TRUE ~ "Other"
  )
}


##############################
## 3) Annotation-text preparation
##
## The local DH13M14 annotation can contain useful information in
## columns other than "function". We concatenate plausible
## annotation/domain columns when present.
##############################
build_annotation_text <- function(annot) {
  
  candidate_cols <- unique(
    c(
      "function",
      "Function",
      "description",
      "Description",
      "protein_name",
      "product",
      "InterPro",
      "interpro",
      "InterProScan",
      "Pfam",
      "pfam",
      "domain",
      "domains",
      "GO",
      "go",
      "KEGG",
      "kegg",
      "KOG",
      "SwissProt",
      "NR",
      "nr"
    )
  )
  
  existing <- intersect(
    candidate_cols,
    colnames(annot)
  )
  
  if (length(existing) == 0L) {
    
    message(
      "No additional annotation/domain columns detected; using function only."
    )
    
    if (!"function" %in% colnames(annot)) {
      stop(
        "No usable annotation column was found."
      )
    }
    
    return(
      as.character(
        annot[["function"]]
      )
    )
  }
  
  x <- annot[
    ,
    existing,
    drop = FALSE
  ]
  
  x[] <- lapply(
    x,
    function(v) {
      v <- as.character(v)
      v[is.na(v)] <- ""
      v
    }
  )
  
  apply(
    x,
    1L,
    function(v) {
      paste(
        unique(
          v[nzchar(v)]
        ),
        collapse = " | "
      )
    }
  )
}


##############################
## 4) Literature-guided immune classification
##############################
classify_immune_annotations <- function(tbl) {
  
  assert_columns(
    tbl,
    c(
      "Gene",
      "gene_ID",
      "annotation_text"
    ),
    "Mapped carrot annotation"
  )
  
  tbl |>
    dplyr::mutate(
      annotation_raw = dplyr::coalesce(
        as.character(annotation_text),
        ""
      ),
      
      annotation_lc = stringr::str_to_lower(
        annotation_raw
      ),
      
      ########################################################
      ## A) Exact / high-confidence cell-surface immune nodes
      ########################################################
      hit_FLS2 = detect_any(
        annotation_lc,
        "\\bfls2\\b|flagellin[- ]sensing 2|flagellin receptor"
      ),
      
      hit_EFR = detect_any(
        annotation_lc,
        "\\befr\\b|ef-tu receptor|elongation factor[- ]tu receptor"
      ),
      
      hit_CERK1 = detect_any(
        annotation_lc,
        "\\bcerk1\\b|chitin elicitor receptor kinase"
      ),
      
      hit_LYK = detect_any(
        annotation_lc,
        "\\blyk[0-9]+\\b|lysin motif.*receptor[- ]like kinase|lysm.*receptor[- ]like kinase"
      ),
      
      hit_LYM_CEBiP = detect_any(
        annotation_lc,
        "\\blym[0-9]+\\b|cebip|chitin elicitor[- ]binding protein|lysm.*receptor[- ]like protein"
      ),
      
      hit_PEPR = detect_any(
        annotation_lc,
        "\\bpepr[12]\\b|pep receptor [12]|plant elicitor peptide receptor"
      ),
      
      hit_MIK2 = detect_any(
        annotation_lc,
        "\\bmik2\\b|male discoverer 1[- ]interacting receptor[- ]like kinase 2"
      ),
      
      hit_RLK7 = detect_any(
        annotation_lc,
        "\\brlk7\\b|receptor[- ]like kinase 7"
      ),
      
      hit_RLP_immune = detect_any(
        annotation_lc,
        "\\brlp23\\b|\\brlp42\\b|\\brlp30\\b|\\bcf[- ]?[249]\\b|cladosporium resistance|\\bve1\\b|elicitin response|\\brxeg1\\b|nlp20 receptor"
      ),
      
      hit_LORE = detect_any(
        annotation_lc,
        "\\blore\\b|lipooligosaccharide[- ]specific reduced elicitation|\\bcard1\\b|cannot respond to dmbq"
      ),
      
      hit_P2K_DORN = detect_any(
        annotation_lc,
        "\\bdorn1\\b|\\bp2k1\\b|does not respond to nucleotides|extracellular atp receptor"
      ),
      
      hit_HSL3 = detect_any(
        annotation_lc,
        "\\bhsl3\\b|haesa[- ]like 3"
      ),
      
      hit_NILR1 = detect_any(
        annotation_lc,
        "\\bnilr1\\b|nematode[- ]induced lrr[- ]rlk"
      ),
      
      hit_BAM1 = detect_any(
        annotation_lc,
        "\\bbam1\\b|barely any meristem 1"
      ),
      
      hit_FER = detect_any(
        annotation_lc,
        "\\bferonia\\b|\\bfer\\b.*receptor|catharanthus roseus receptor[- ]like kinase 1[- ]like|crrlk1l"
      ),
      
      ########################################################
      ## B) Broader receptor families relevant to immunity
      ##    Individual function is NOT assumed.
      ########################################################
      hit_LRR_RLK = detect_any(
        annotation_lc,
        "leucine[- ]rich repeat receptor[- ]like kinase|lrr[- ]rlk|lrr receptor kinase|leucine[- ]rich repeat.*protein kinase"
      ),
      
      hit_LysM_RLK = detect_any(
        annotation_lc,
        "lysm.*receptor[- ]like kinase|lysin motif.*receptor[- ]like kinase"
      ),
      
      hit_LRR_RLP = detect_any(
        annotation_lc,
        "leucine[- ]rich repeat receptor[- ]like protein|lrr[- ]rlp|receptor[- ]like protein.*leucine[- ]rich"
      ),
      
      hit_LecRLK = detect_any(
        annotation_lc,
        "lectin receptor[- ]like kinase|lecrk|g[- ]type lectin.*receptor[- ]like kinase|s[- ]type lectin.*receptor[- ]like kinase|c[- ]type lectin.*receptor[- ]like kinase"
      ),
      
      hit_CRK = detect_any(
        annotation_lc,
        "cysteine[- ]rich receptor[- ]like kinase|\\bcrk[0-9]+\\b|duf26.*receptor[- ]like kinase"
      ),
      
      hit_WAK = detect_any(
        annotation_lc,
        "wall[- ]associated kinase|\\bwakl?[0-9]*\\b"
      ),
      
      hit_CrRLK1L = detect_any(
        annotation_lc,
        "crrlk1l|catharanthus roseus receptor[- ]like kinase"
      ),
      
      hit_generic_RLK = detect_any(
        annotation_lc,
        "receptor[- ]like kinase|receptor kinase"
      ),
      
      hit_generic_RLP = detect_any(
        annotation_lc,
        "receptor[- ]like protein"
      ),
      
      ########################################################
      ## C) Co-receptors / scaffolds / receptor regulators
      ########################################################
      hit_BAK1_SERK = detect_any(
        annotation_lc,
        "\\bbak1\\b|serk3|bri1[- ]associated receptor kinase 1|somatic embryogenesis receptor kinase"
      ),
      
      hit_SOBIR1 = detect_any(
        annotation_lc,
        "\\bsobir1\\b|suppressor of bir1"
      ),
      
      hit_BIR = detect_any(
        annotation_lc,
        "\\bbir[1234]\\b|bak1[- ]interacting receptor[- ]like kinase"
      ),
      
      hit_IOS1 = detect_any(
        annotation_lc,
        "\\bios1\\b|impaired oomycete susceptibility 1"
      ),
      
      hit_QSK1 = detect_any(
        annotation_lc,
        "\\bqsk1\\b|qian shou kinase 1"
      ),
      
      ########################################################
      ## D) Receptor-like cytoplasmic kinases / proximal PTI
      ########################################################
      hit_RLCK = detect_any(
        annotation_lc,
        "receptor[- ]like cytoplasmic kinase|\\brlck\\b"
      ),
      
      hit_BIK1 = detect_any(
        annotation_lc,
        "\\bbik1\\b|botrytis[- ]induced kinase 1"
      ),
      
      hit_PBL = detect_any(
        annotation_lc,
        "\\bpbl[0-9]+\\b|pbs1[- ]like"
      ),
      
      hit_PBS1 = detect_any(
        annotation_lc,
        "\\bpbs1\\b|avrpphb susceptible 1"
      ),
      
      hit_BSK1 = detect_any(
        annotation_lc,
        "\\bbsk1\\b|brassinosteroid[- ]signaling kinase 1"
      ),
      
      hit_RIPK = detect_any(
        annotation_lc,
        "\\bripk\\b|rpm1[- ]induced protein kinase"
      ),
      
      hit_PCRK = detect_any(
        annotation_lc,
        "\\bpcrk[12]\\b|pti compromised receptor[- ]like cytoplasmic kinase"
      ),
      
      hit_OXI1 = detect_any(
        annotation_lc,
        "\\boxi1\\b|oxidative signal[- ]inducible 1"
      ),
      
      ########################################################
      ## E) Ca2+ signalling / channels
      ########################################################
      hit_CNGC = detect_any(
        annotation_lc,
        "cyclic nucleotide[- ]gated channel|\\bcngc[0-9]+\\b"
      ),
      
      hit_CPK_CDPK = detect_any(
        annotation_lc,
        "calcium[- ]dependent protein kinase|\\bcdpk\\b|\\bcpk[0-9]+\\b"
      ),
      
      hit_CaM_CML = detect_any(
        annotation_lc,
        "calmodulin[- ]like|calmodulin protein|\\bcml[0-9]+\\b"
      ),
      
      hit_ANN = detect_any(
        annotation_lc,
        "\\bannexin\\b|\\bann[14]\\b"
      ),
      
      hit_OSCA = detect_any(
        annotation_lc,
        "\\bosca1\\.3\\b|hyperosmolality[- ]gated calcium[- ]permeable channel"
      ),
      
      ########################################################
      ## F) ROS / redox / lipid signalling
      ########################################################
      hit_RBOH = detect_any(
        annotation_lc,
        "respiratory burst oxidase homolog|\\brboh[a-j]\\b|nadph oxidase.*respiratory burst"
      ),
      
      hit_NDPK2 = detect_any(
        annotation_lc,
        "\\bndpk2\\b|nucleoside diphosphate kinase 2"
      ),
      
      hit_DGK5 = detect_any(
        annotation_lc,
        "\\bdgk5\\b|diacylglycerol kinase 5"
      ),
      
      hit_HPCA1 = detect_any(
        annotation_lc,
        "\\bhpca1\\b|hydrogen[- ]peroxide[- ]induced ca2\\+ increases 1"
      ),
      
      ########################################################
      ## G) MAPK cascade
      ########################################################
      hit_MAPK = detect_any(
        annotation_lc,
        "mitogen[- ]activated protein kinase|map kinase|\\bmpk[0-9]+\\b|\\bmapk[0-9]+\\b"
      ),
      
      hit_MKK = detect_any(
        annotation_lc,
        "mitogen[- ]activated protein kinase kinase|\\bmkk[0-9]+\\b"
      ),
      
      hit_MAPKKK = detect_any(
        annotation_lc,
        "mitogen[- ]activated protein kinase kinase kinase|\\bmapkkk[0-9]+\\b"
      ),
      
      ########################################################
      ## H) NLR / intracellular immune receptors
      ########################################################
      hit_NLR_domain = detect_any(
        annotation_lc,
        "\\bnb[- ]arc\\b|nucleotide[- ]binding.*leucine[- ]rich repeat|nbs[- ]lrr|nb[- ]lrr|tir[- ]nbs|cc[- ]nbs|rpw8.*nbs|nucleotide[- ]binding site.*leucine[- ]rich"
      ),
      
      hit_NLR_named = detect_any(
        annotation_lc,
        "\\brps[2456]\\b|\\brpm1\\b|\\bsnc1\\b|\\brpp[0-9]+\\b|resistance to pseudomonas|disease resistance protein.*nbs|disease resistance protein.*nb[- ]arc"
      ),
      
      hit_TNL = detect_any(
        annotation_lc,
        "tir[- ]nbs[- ]lrr|tir.*nb[- ]arc.*lrr|toll.*interleukin.*nucleotide[- ]binding"
      ),
      
      hit_CNL = detect_any(
        annotation_lc,
        "cc[- ]nbs[- ]lrr|coiled[- ]coil.*nb[- ]arc.*lrr"
      ),
      
      hit_RNL = detect_any(
        annotation_lc,
        "rpw8.*nb[- ]arc|helper nlr|\\badr1\\b.*nlr|\\bnrg1\\b.*nlr"
      ),
      
      ########################################################
      ## I) ETI support / guardee machinery
      ########################################################
      hit_EDS1 = detect_any(
        annotation_lc,
        "\\beds1\\b|enhanced disease susceptibility 1"
      ),
      
      hit_PAD4 = detect_any(
        annotation_lc,
        "\\bpad4\\b|phytoalexin deficient 4"
      ),
      
      hit_SAG101 = detect_any(
        annotation_lc,
        "\\bsag101\\b|senescence[- ]associated gene 101"
      ),
      
      hit_ADR1 = detect_any(
        annotation_lc,
        "\\badr1\\b|activated disease resistance 1"
      ),
      
      hit_NRG1 = detect_any(
        annotation_lc,
        "\\bnrg1\\b|n requirement gene 1"
      ),
      
      hit_NDR1 = detect_any(
        annotation_lc,
        "\\bndr1\\b|non[- ]race specific disease resistance 1"
      ),
      
      hit_RIN4 = detect_any(
        annotation_lc,
        "\\brin4\\b|rpm1[- ]interacting protein 4"
      ),
      
      hit_SGT1 = detect_any(
        annotation_lc,
        "\\bsgt1\\b|suppressor of g2 allele of skp1"
      ),
      
      hit_RAR1 = detect_any(
        annotation_lc,
        "\\brar1\\b|required for mla12 resistance"
      ),
      
      hit_HSP90 = detect_any(
        annotation_lc,
        "heat shock protein 90|\\bhsp90\\b"
      ),
      
      ########################################################
      ## J) SA / NHP / systemic immunity
      ########################################################
      hit_CBP60g = detect_any(
        annotation_lc,
        "\\bcbp60g\\b|calmodulin[- ]binding protein 60[- ]like g"
      ),
      
      hit_SARD1 = detect_any(
        annotation_lc,
        "\\bsard1\\b|systemic acquired resistance deficient 1"
      ),
      
      hit_ICS1 = detect_any(
        annotation_lc,
        "\\bics1\\b|\\bsid2\\b|isochorismate synthase 1"
      ),
      
      hit_EDS5 = detect_any(
        annotation_lc,
        "\\beds5\\b|enhanced disease susceptibility 5"
      ),
      
      hit_PBS3 = detect_any(
        annotation_lc,
        "\\bpbs3\\b|avrpphb susceptible 3"
      ),
      
      hit_NPR = detect_any(
        annotation_lc,
        "\\bnpr[134]\\b|nonexpressor of pathogenesis[- ]related genes"
      ),
      
      hit_ALD1 = detect_any(
        annotation_lc,
        "\\bald1\\b|ag2[- ]like defense response protein 1"
      ),
      
      hit_FMO1 = detect_any(
        annotation_lc,
        "\\bfmo1\\b|flavin[- ]dependent monooxygenase 1"
      ),
      
      ########################################################
      ## K) JA / ethylene defence crosstalk
      ########################################################
      hit_COI1 = detect_any(
        annotation_lc,
        "\\bcoi1\\b|coronatine insensitive 1"
      ),
      
      hit_JAZ = detect_any(
        annotation_lc,
        "\\bjaz[0-9]+\\b|jasmonate zim[- ]domain"
      ),
      
      hit_MYC2 = detect_any(
        annotation_lc,
        "\\bmyc2\\b"
      ),
      
      hit_EIN2 = detect_any(
        annotation_lc,
        "\\bein2\\b|ethylene insensitive 2"
      ),
      
      hit_EIN3_EIL = detect_any(
        annotation_lc,
        "\\bein3\\b|ein3[- ]like|ethylene insensitive 3"
      ),
      
      hit_ACS6 = detect_any(
        annotation_lc,
        "\\bacs6\\b|1[- ]aminocyclopropane[- ]1[- ]carboxylate synthase 6"
      ),
      
      ########################################################
      ## L) Transcriptional immune nodes
      ########################################################
      hit_WRKY22_29 = detect_any(
        annotation_lc,
        "\\bwrky22\\b|\\bwrky29\\b"
      ),
      
      hit_TGA = detect_any(
        annotation_lc,
        "\\btga[0-9]+\\b|tga transcription factor"
      ),
      
      ########################################################
      ## M) Defence outputs
      ########################################################
      hit_PR1 = detect_any(
        annotation_lc,
        "pathogenesis[- ]related protein 1|\\bpr1\\b"
      ),
      
      hit_chitinase = detect_any(
        annotation_lc,
        "pathogenesis[- ]related.*chitinase|endochitinase|class .*chitinase"
      ),
      
      hit_glucanase = detect_any(
        annotation_lc,
        "beta[- ]1,3[- ]glucanase|β[- ]1,3[- ]glucanase|pathogenesis[- ]related.*glucanase"
      ),
      
      hit_callose = detect_any(
        annotation_lc,
        "\\bpmr4\\b|\\bgsl5\\b|callose synthase"
      ),
      
      hit_defensin = detect_any(
        annotation_lc,
        "plant defensin|pathogenesis[- ]related.*defensin"
      ),
      
      ########################################################
      ## N) Beneficial-microbe / symbiosis perception
      ########################################################
      hit_NFR_NFP = detect_any(
        annotation_lc,
        "\\bnfr[15]\\b|nod factor receptor|\\bnfp\\b"
      ),
      
      hit_SYMRK = detect_any(
        annotation_lc,
        "\\bsymrk\\b|\\bdmi2\\b|symbiosis receptor[- ]like kinase"
      ),
      
      hit_EPR3 = detect_any(
        annotation_lc,
        "\\bepr3\\b|exopolysaccharide receptor 3"
      ),
      
      ########################################################
      ## O) Generic immune/pathogen annotation
      ########################################################
      explicit_immune_language = detect_any(
        annotation_lc,
        "plant immunity|immune response|pathogen response|disease resistance|defence response|defense response|microbe interaction|microbial recognition|pattern recognition receptor|pamp receptor|mamp receptor"
      ),
      
      ########################################################
      ## Primary family flags
      ########################################################
      exact_surface_PRR =
        hit_FLS2 |
        hit_EFR |
        hit_CERK1 |
        hit_LYK |
        hit_LYM_CEBiP |
        hit_PEPR |
        hit_MIK2 |
        hit_RLK7 |
        hit_RLP_immune |
        hit_LORE |
        hit_P2K_DORN |
        hit_HSL3 |
        hit_NILR1 |
        hit_BAM1 |
        hit_FER,
      
      broad_surface_candidate =
        hit_LRR_RLK |
        hit_LysM_RLK |
        hit_LRR_RLP |
        hit_LecRLK |
        hit_CRK |
        hit_WAK |
        hit_CrRLK1L,
      
      coreceptor_scaffold =
        hit_BAK1_SERK |
        hit_SOBIR1 |
        hit_BIR |
        hit_IOS1 |
        hit_QSK1,
      
      receptor_proximal =
        hit_RLCK |
        hit_BIK1 |
        hit_PBL |
        hit_PBS1 |
        hit_BSK1 |
        hit_RIPK |
        hit_PCRK |
        hit_OXI1,
      
      calcium_module =
        hit_CNGC |
        hit_CPK_CDPK |
        hit_CaM_CML |
        hit_ANN |
        hit_OSCA,
      
      ros_module =
        hit_RBOH |
        hit_NDPK2 |
        hit_DGK5 |
        hit_HPCA1,
      
      mapk_module =
        hit_MAPK |
        hit_MKK |
        hit_MAPKKK,
      
      nlr_module =
        hit_NLR_domain |
        hit_NLR_named,
      
      eti_support =
        hit_EDS1 |
        hit_PAD4 |
        hit_SAG101 |
        hit_ADR1 |
        hit_NRG1 |
        hit_NDR1 |
        hit_RIN4 |
        hit_SGT1 |
        hit_RAR1 |
        hit_HSP90,
      
      sa_nhp_module =
        hit_CBP60g |
        hit_SARD1 |
        hit_ICS1 |
        hit_EDS5 |
        hit_PBS3 |
        hit_NPR |
        hit_ALD1 |
        hit_FMO1,
      
      hormone_module =
        hit_COI1 |
        hit_JAZ |
        hit_MYC2 |
        hit_EIN2 |
        hit_EIN3_EIL |
        hit_ACS6,
      
      immune_tf_module =
        hit_WRKY22_29 |
        hit_TGA |
        hit_CBP60g |
        hit_SARD1,
      
      defence_output_module =
        hit_PR1 |
        hit_chitinase |
        hit_glucanase |
        hit_callose |
        hit_defensin,
      
      symbiosis_module =
        hit_NFR_NFP |
        hit_SYMRK |
        hit_EPR3,
      
      generic_receptor_only =
        (hit_generic_RLK | hit_generic_RLP) &
        !exact_surface_PRR &
        !broad_surface_candidate &
        !coreceptor_scaffold &
        !receptor_proximal,
      
      ########################################################
      ## Primary module
      ########################################################
      Immune_module_auto = dplyr::case_when(
        exact_surface_PRR |
          broad_surface_candidate ~
          "01_Cell-surface immune receptors",
        
        coreceptor_scaffold ~
          "02_Co-receptors / receptor scaffolds",
        
        receptor_proximal ~
          "03_Receptor-proximal RLCKs",
        
        calcium_module ~
          "04_Ca2+ signalling",
        
        ros_module ~
          "05_ROS / redox signalling",
        
        mapk_module ~
          "06_MAPK signalling",
        
        nlr_module ~
          "07_Intracellular NLR receptors",
        
        eti_support ~
          "08_ETI support / guardee machinery",
        
        sa_nhp_module ~
          "09_SA / NHP / systemic immunity",
        
        hormone_module ~
          "10_JA / ethylene defence crosstalk",
        
        immune_tf_module ~
          "11_Immune transcriptional regulators",
        
        defence_output_module ~
          "12_Defence outputs / PR / cell wall",
        
        symbiosis_module ~
          "13_Beneficial-microbe / symbiosis perception",
        
        explicit_immune_language |
          generic_receptor_only ~
          "14_Other immune-associated candidate",
        
        TRUE ~ NA_character_
      ),
      
      ########################################################
      ## Submodule / interpretable short class
      ########################################################
      Immune_submodule_auto = dplyr::case_when(
        hit_FLS2 ~ "FLS2-like flagellin PRR",
        hit_EFR ~ "EFR-like EF-Tu PRR",
        hit_CERK1 ~ "CERK1-like chitin signalling",
        hit_LYK ~ "LYK/LysM receptor kinase",
        hit_LYM_CEBiP ~ "LysM receptor protein / CEBiP-like",
        hit_PEPR ~ "PEPR DAMP receptor",
        hit_MIK2 ~ "MIK2-like receptor",
        hit_RLK7 ~ "RLK7-like peptide receptor",
        hit_RLP_immune ~ "Immune LRR-RLP-like",
        hit_LORE ~ "LORE/CARD1-like PRR",
        hit_P2K_DORN ~ "P2K1/DORN1 DAMP receptor",
        hit_HSL3 ~ "HSL3-like receptor",
        hit_NILR1 ~ "NILR1-like nematode receptor",
        hit_BAM1 ~ "BAM1-like receptor kinase",
        hit_FER ~ "FERONIA / CrRLK1L",
        
        hit_LysM_RLK ~ "LysM-RLK candidate",
        hit_LRR_RLK ~ "LRR-RLK candidate",
        hit_LRR_RLP ~ "LRR-RLP candidate",
        hit_LecRLK ~ "Lectin receptor-like kinase",
        hit_CRK ~ "Cysteine-rich RLK / CRK",
        hit_WAK ~ "Wall-associated kinase",
        hit_CrRLK1L ~ "CrRLK1L candidate",
        
        hit_BAK1_SERK ~ "BAK1/SERK co-receptor",
        hit_SOBIR1 ~ "SOBIR1 RLP co-receptor",
        hit_BIR ~ "BIR receptor regulator",
        hit_IOS1 ~ "IOS1 receptor scaffold",
        hit_QSK1 ~ "QSK1 PRR-complex regulator",
        
        hit_BIK1 ~ "BIK1 RLCK",
        hit_PBS1 ~ "PBS1 RLCK",
        hit_PBL ~ "PBL RLCK",
        hit_BSK1 ~ "BSK1 RLCK",
        hit_RIPK ~ "RIPK RLCK",
        hit_PCRK ~ "PCRK RLCK",
        hit_OXI1 ~ "OXI1 kinase",
        hit_RLCK ~ "RLCK candidate",
        
        hit_CNGC ~ "CNGC Ca2+ channel",
        hit_CPK_CDPK ~ "CPK/CDPK",
        hit_CaM_CML ~ "CaM/CML",
        hit_ANN ~ "Annexin",
        hit_OSCA ~ "OSCA Ca2+-permeable channel",
        
        hit_RBOH ~ "RBOH NADPH oxidase",
        hit_NDPK2 ~ "NDPK2",
        hit_DGK5 ~ "DGK5",
        hit_HPCA1 ~ "HPCA1 H2O2 sensor",
        
        hit_MAPKKK ~ "MAPKKK",
        hit_MKK ~ "MKK/MAPKK",
        hit_MAPK ~ "MAPK/MPK",
        
        hit_TNL ~ "TNL candidate",
        hit_CNL ~ "CNL candidate",
        hit_RNL ~ "RNL/helper-NLR candidate",
        hit_NLR_named ~ "Named NLR/R gene-like",
        hit_NLR_domain ~ "NLR / NB-ARC candidate",
        
        hit_EDS1 ~ "EDS1",
        hit_PAD4 ~ "PAD4",
        hit_SAG101 ~ "SAG101",
        hit_ADR1 ~ "ADR1 helper-NLR/hub",
        hit_NRG1 ~ "NRG1 helper-NLR/hub",
        hit_NDR1 ~ "NDR1",
        hit_RIN4 ~ "RIN4 guardee",
        hit_SGT1 ~ "SGT1",
        hit_RAR1 ~ "RAR1",
        hit_HSP90 ~ "HSP90 ETI cofactor",
        
        hit_CBP60g ~ "CBP60g",
        hit_SARD1 ~ "SARD1",
        hit_ICS1 ~ "ICS1/SID2",
        hit_EDS5 ~ "EDS5",
        hit_PBS3 ~ "PBS3",
        hit_NPR ~ "NPR SA receptor/signalling",
        hit_ALD1 ~ "ALD1",
        hit_FMO1 ~ "FMO1",
        
        hit_COI1 ~ "COI1",
        hit_JAZ ~ "JAZ",
        hit_MYC2 ~ "MYC2",
        hit_EIN2 ~ "EIN2",
        hit_EIN3_EIL ~ "EIN3/EIL",
        hit_ACS6 ~ "ACS6",
        
        hit_WRKY22_29 ~ "WRKY22/29 PTI regulator",
        hit_TGA ~ "TGA defence TF",
        
        hit_PR1 ~ "PR1",
        hit_chitinase ~ "Defence chitinase",
        hit_glucanase ~ "β-1,3-glucanase",
        hit_callose ~ "Callose synthase",
        hit_defensin ~ "Plant defensin",
        
        hit_NFR_NFP ~ "Nod-factor receptor-like",
        hit_SYMRK ~ "SYMRK/DMI2-like",
        hit_EPR3 ~ "EPR3 exopolysaccharide receptor",
        
        generic_receptor_only ~ "Generic RLK/RLP candidate",
        explicit_immune_language ~ "Immune-associated annotation",
        
        TRUE ~ NA_character_
      ),
      
      ########################################################
      ## Confidence
      ########################################################
      Evidence_level_auto = dplyr::case_when(
        exact_surface_PRR |
          coreceptor_scaffold |
          hit_BIK1 |
          hit_PBS1 |
          hit_OXI1 |
          hit_RBOH |
          hit_EDS1 |
          hit_PAD4 |
          hit_RIN4 |
          hit_SGT1 |
          hit_CBP60g |
          hit_SARD1 |
          hit_ICS1 |
          hit_ALD1 |
          hit_FMO1 |
          hit_PR1 |
          hit_WRKY22_29 |
          symbiosis_module ~
          "High",
        
        broad_surface_candidate |
          receptor_proximal |
          calcium_module |
          ros_module |
          mapk_module |
          nlr_module |
          eti_support |
          sa_nhp_module |
          hormone_module |
          defence_output_module |
          explicit_immune_language ~
          "Medium",
        
        generic_receptor_only ~
          "Low",
        
        TRUE ~
          "Not_selected"
      ),
      
      ########################################################
      ## Carrot A. dauci publication anchor
      ##
      ## Liu et al. 2024 directly examined these categories/nodes
      ## in carrot during A. dauci infection.
      ########################################################
      Carrot_Adauci_anchor =
        hit_LRR_RLK |
        hit_LysM_RLK |
        nlr_module |
        hit_CNGC |
        hit_CPK_CDPK |
        hit_CaM_CML |
        hit_RBOH |
        hit_OXI1 |
        hit_NDPK2 |
        hit_WRKY22_29 |
        hit_ACS6 |
        hit_PR1 |
        hit_EDS1 |
        hit_PBS1 |
        hit_RIN4 |
        hit_NLR_named |
        hit_SGT1 |
        hit_HSP90,
      
      ########################################################
      ## Direct heat x biotic / immunity literature priority
      ##
      ## Only explicit homolog-like names or documented families
      ## are flagged. This is a prioritisation layer, not proof that
      ## the carrot orthologue has the same function.
      ########################################################
      Direct_heat_immunity_anchor =
        hit_FLS2 |
        hit_BAM1 |
        hit_PBS1 |
        hit_RBOH |
        hit_FER |
        hit_BAK1_SERK |
        hit_SOBIR1 |
        hit_CRK |
        hit_LecRLK |
        hit_MAPK |
        hit_CBP60g |
        hit_SARD1 |
        hit_ICS1 |
        hit_EDS1 |
        hit_PAD4 |
        hit_ALD1 |
        hit_FMO1,
      
      Heat_immunity_note = dplyr::case_when(
        hit_FER ~
          "FER: direct plasma-membrane thermotolerance evidence plus established immune-receptor-complex regulation.",
        
        hit_BAM1 |
          hit_PBS1 |
          hit_RBOH ~
          "BAM1-PBS1-RBOHD: direct receptor-complex evidence for early heat signalling; components overlap with immune signalling.",
        
        hit_FLS2 ~
          "FLS2/PTI is heat-sensitive in experimental Arabidopsis systems; receptor state should not be inferred from transcript alone.",
        
        hit_BAK1_SERK ~
          "SERK/BAK1 proteins are major PRR co-receptors; crop studies link SERK-family members to high-temperature disease resistance.",
        
        hit_CRK ~
          "CRK family has experimental links to both heat responses and pathogen resistance in crops; individual homolog function requires validation.",
        
        hit_LecRLK ~
          "Lectin RLK family includes experimentally validated heat/biotic trade-off regulators; individual homolog function requires validation.",
        
        hit_SOBIR1 ~
          "SOBIR1 is central to RLP immunity and has direct high-temperature/biotic crosstalk evidence in potato.",
        
        hit_MAPK ~
          "MAPK signalling can coordinate pathogen resistance and thermotolerance; individual MAPKs differ functionally.",
        
        hit_CBP60g |
          hit_SARD1 |
          hit_ICS1 |
          hit_EDS1 |
          hit_PAD4 ~
          "SA immune module is a major temperature-sensitive defence sector.",
        
        hit_ALD1 |
          hit_FMO1 ~
          "NHP/SAR module is suppressed by warm temperature in experimental systems.",
        
        TRUE ~
          ""
      ),
      
      ########################################################
      ## Main inclusion
      ########################################################
      Include_auto = dplyr::case_when(
        is.na(Immune_module_auto) ~ FALSE,
        
        Evidence_level_auto %in% c(
          "High",
          "Medium"
        ) ~ TRUE,
        
        Evidence_level_auto == "Low" &
          INCLUDE_LOW_CONFIDENCE_GENERIC_RLK_IN_MAIN ~ TRUE,
        
        TRUE ~ FALSE
      )
    ) |>
    dplyr::mutate(
      # Optional removal of TF/HSP nodes from main figures.
      Include_auto = dplyr::if_else(
        !INCLUDE_IMMUNE_TF_HSP_NODES &
          (
            hit_WRKY22_29 |
              hit_TGA |
              hit_CBP60g |
              hit_SARD1 |
              hit_HSP90
          ),
        FALSE,
        Include_auto
      ),
      
      Evidence_reason_auto = dplyr::case_when(
        Evidence_level_auto == "High" ~
          "Explicit known immune/microbe-interaction component in supplied annotation.",
        
        Evidence_level_auto == "Medium" ~
          "Immune-relevant receptor/signalling family or pathway component; individual carrot gene requires functional validation.",
        
        Evidence_level_auto == "Low" ~
          "Generic receptor-like annotation without sufficient immune-specific/domain evidence.",
        
        TRUE ~
          "No sufficiently specific immune/microbe-interaction evidence."
      )
    )
}


##############################
## 5) Manual curation
##############################
apply_manual_overrides <- function(
    tbl,
    override_file) {
  
  template <- tbl |>
    dplyr::transmute(
      Gene,
      gene_ID,
      Immune_module_auto,
      Immune_submodule_auto,
      Evidence_level_auto,
      Include_auto,
      Module_override = "",
      Submodule_override = "",
      Short_name_override = "",
      Evidence_override = "",
      Include_override = "",
      Curator_note = ""
    )
  
  if (!file.exists(override_file)) {
    
    dir.create(
      dirname(override_file),
      recursive = TRUE,
      showWarnings = FALSE
    )
    
    readr::write_csv(
      template,
      override_file
    )
    
    message(
      "Created manual-override template: ",
      override_file
    )
  }
  
  ov <- read_csv_flexible(
    override_file
  )
  
  assert_columns(
    ov,
    c("Gene"),
    "Immune manual override file"
  )
  
  for (
    nm in c(
      "Module_override",
      "Submodule_override",
      "Short_name_override",
      "Evidence_override",
      "Include_override",
      "Curator_note"
    )
  ) {
    if (!nm %in% colnames(ov)) {
      ov[[nm]] <- ""
    }
  }
  
  ov <- ov |>
    dplyr::mutate(
      dplyr::across(
        c(
          Module_override,
          Submodule_override,
          Short_name_override,
          Evidence_override,
          Include_override,
          Curator_note
        ),
        ~ trimws(
          as.character(.x)
        )
      ),
      
      Include_override_parsed =
        parse_logical_override(
          Include_override
        )
    ) |>
    dplyr::select(
      Gene,
      Module_override,
      Submodule_override,
      Short_name_override,
      Evidence_override,
      Include_override_parsed,
      Curator_note
    )
  
  tbl |>
    dplyr::left_join(
      ov,
      by = "Gene"
    ) |>
    dplyr::mutate(
      Immune_module = dplyr::if_else(
        !is.na(Module_override) &
          nzchar(Module_override),
        Module_override,
        Immune_module_auto
      ),
      
      Immune_submodule = dplyr::if_else(
        !is.na(Submodule_override) &
          nzchar(Submodule_override),
        Submodule_override,
        Immune_submodule_auto
      ),
      
      Evidence_level = dplyr::if_else(
        !is.na(Evidence_override) &
          nzchar(Evidence_override),
        Evidence_override,
        Evidence_level_auto
      ),
      
      Short_name = dplyr::if_else(
        !is.na(Short_name_override) &
          nzchar(Short_name_override),
        Short_name_override,
        dplyr::coalesce(
          Immune_submodule,
          ""
        )
      ),
      
      Include = dplyr::if_else(
        !is.na(Include_override_parsed),
        Include_override_parsed,
        Include_auto
      ),
      
      Curator_note = dplyr::coalesce(
        Curator_note,
        ""
      )
    )
}


##############################
## 6) Heatmap functions
##############################
make_top_annotation <- function() {
  
  genotype <- c(
    rep("ROBILA", 3),
    rep("PRESTO", 3)
  )
  
  condition <- c(
    "NH", "HS-7", "HS-2",
    "NH", "HS-7", "HS-2"
  )
  
  ComplexHeatmap::HeatmapAnnotation(
    Genotype = genotype,
    Condition = condition,
    
    col = list(
      Genotype = geno_colors,
      
      Condition = c(
        "NH" = unname(
          temp_colors["T0"]
        ),
        
        "HS-7" = unname(
          temp_colors["T1"]
        ),
        
        "HS-2" = unname(
          temp_colors["T2"]
        )
      )
    ),
    
    annotation_name_gp = grid::gpar(
      fontsize = 8
    ),
    
    simple_anno_size = grid::unit(
      3.5,
      "mm"
    ),
    
    show_legend = TRUE
  )
}

draw_heatmap_object <- function(ht) {
  
  ComplexHeatmap::draw(
    ht,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    merge_legends = TRUE,
    
    padding = grid::unit(
      c(5, 5, 5, 5),
      "mm"
    )
  )
}

save_complex_heatmap <- function(
    ht,
    file_base,
    width_mm,
    height_mm,
    allow_raster = TRUE) {
  
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  
  dir.create(
    dirname(file_base),
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  render_device <- function(
    label,
    open_device) {
    
    device_before <- grDevices::dev.cur()
    
    tryCatch({
      
      open_device()
      draw_heatmap_object(ht)
      grDevices::dev.off()
      
      invisible(TRUE)
      
    }, error = function(e) {
      
      if (
        grDevices::dev.cur() !=
        device_before
      ) {
        try(
          grDevices::dev.off(),
          silent = TRUE
        )
      }
      
      warning(
        label,
        " export failed for ",
        file_base,
        ": ",
        conditionMessage(e)
      )
      
      invisible(FALSE)
    })
  }
  
  if (IMMUNE_WRITE_PDF) {
    
    render_device(
      "PDF",
      function() {
        
        grDevices::pdf(
          paste0(
            file_base,
            ".pdf"
          ),
          width = width_in,
          height = height_in,
          useDingbats = FALSE,
          onefile = TRUE
        )
      }
    )
  }
  
  if (IMMUNE_WRITE_SVG) {
    
    render_device(
      "SVG",
      function() {
        
        svglite::svglite(
          paste0(
            file_base,
            ".svg"
          ),
          width = width_in,
          height = height_in,
          bg = "white"
        )
      }
    )
  }
  
  raster_ok <-
    allow_raster &&
    height_mm <=
    IMMUNE_RASTER_MAX_HEIGHT_MM
  
  if (
    IMMUNE_WRITE_PNG &&
    raster_ok
  ) {
    
    render_device(
      "PNG",
      function() {
        
        grDevices::png(
          paste0(
            file_base,
            ".png"
          ),
          width = width_in,
          height = height_in,
          units = "in",
          res = IMMUNE_RASTER_DPI,
          type = "cairo-png",
          bg = "white"
        )
      }
    )
  }
  
  if (
    IMMUNE_WRITE_TIFF &&
    raster_ok
  ) {
    
    render_device(
      "TIFF",
      function() {
        
        grDevices::tiff(
          paste0(
            file_base,
            ".tiff"
          ),
          width = width_in,
          height = height_in,
          units = "in",
          res = IMMUNE_RASTER_DPI,
          compression = "lzw",
          type = "cairo",
          bg = "white"
        )
      }
    )
  }
  
  if (
    !raster_ok &&
    (
      IMMUNE_WRITE_PNG |
      IMMUNE_WRITE_TIFF
    )
  ) {
    
    message(
      "Raster skipped for very tall heatmap (",
      round(
        height_mm,
        1
      ),
      " mm). PDF/SVG remain complete: ",
      file_base
    )
  }
}

make_global_expression_heatmap <- function(
    zmat,
    title,
    cluster_rows = TRUE,
    show_row_names = FALSE) {
  
  if (nrow(zmat) < 2L) {
    stop(
      "At least two genes are required."
    )
  }
  
  zplot <- pmax(
    pmin(
      zmat,
      IMMUNE_GLOBAL_Z_LIMIT
    ),
    -IMMUNE_GLOBAL_Z_LIMIT
  )
  
  z_col_fun <- circlize::colorRamp2(
    c(
      -IMMUNE_GLOBAL_Z_LIMIT,
      0,
      IMMUNE_GLOBAL_Z_LIMIT
    ),
    c(
      "#2166AC",
      "#F7F7F7",
      "#B2182B"
    )
  )
  
  ComplexHeatmap::Heatmap(
    zplot,
    name = "Row z-score",
    col = z_col_fun,
    na_col = "grey90",
    
    top_annotation =
      make_top_annotation(),
    
    column_title = title,
    
    column_title_gp =
      grid::gpar(
        fontsize = 11,
        fontface = "bold"
      ),
    
    cluster_columns = FALSE,
    
    cluster_rows =
      cluster_rows,
    
    clustering_distance_rows =
      "euclidean",
    
    clustering_method_rows =
      "complete",
    
    show_row_dend = TRUE,
    
    row_dend_width =
      grid::unit(
        12,
        "mm"
      ),
    
    show_row_names =
      show_row_names,
    
    row_names_gp =
      grid::gpar(
        fontsize = 5
      ),
    
    column_names_rot = 90,
    
    column_names_gp =
      grid::gpar(
        fontsize = 8
      ),
    
    rect_gp =
      grid::gpar(
        col = NA
      ),
    
    border = FALSE,
    
    use_raster =
      nrow(zplot) > 200L,
    
    raster_quality = 3,
    
    heatmap_legend_param = list(
      title = "Expression",
      
      at = c(
        -IMMUNE_GLOBAL_Z_LIMIT,
        0,
        IMMUNE_GLOBAL_Z_LIMIT
      ),
      
      labels = c(
        paste0(
          "<= -",
          IMMUNE_GLOBAL_Z_LIMIT
        ),
        "0",
        paste0(
          ">= ",
          IMMUNE_GLOBAL_Z_LIMIT
        )
      )
    )
  )
}

get_row_order <- function(
    mat,
    cluster_rows = TRUE) {
  
  if (
    !cluster_rows ||
    nrow(mat) <= 2L
  ) {
    return(
      seq_len(
        nrow(mat)
      )
    )
  }
  
  stats::hclust(
    stats::dist(mat),
    method = "complete"
  )$order
}

make_heatmap_pair <- function(
    zmat,
    status_mat,
    labels,
    title,
    row_split = NULL,
    cluster_rows = TRUE,
    show_row_names = TRUE) {
  
  stopifnot(
    identical(
      rownames(zmat),
      rownames(status_mat)
    )
  )
  
  stopifnot(
    all(
      rownames(zmat) %in%
        names(labels)
    )
  )
  
  zplot <- pmax(
    pmin(
      zmat,
      IMMUNE_Z_LIMIT
    ),
    -IMMUNE_Z_LIMIT
  )
  
  if (is.null(row_split)) {
    
    ord <- get_row_order(
      zplot,
      cluster_rows
    )
    
    zplot <-
      zplot[
        ord,
        ,
        drop = FALSE
      ]
    
    status_mat <-
      status_mat[
        ord,
        ,
        drop = FALSE
      ]
    
    labels <-
      labels[
        rownames(zplot)
      ]
  }
  
  z_col_fun <- circlize::colorRamp2(
    c(
      -IMMUNE_Z_LIMIT,
      0,
      IMMUNE_Z_LIMIT
    ),
    c(
      "#2166AC",
      "#F7F7F7",
      "#B2182B"
    )
  )
  
  hm_expr <- ComplexHeatmap::Heatmap(
    zplot,
    name = "Row z-score",
    col = z_col_fun,
    na_col = "grey90",
    
    top_annotation =
      make_top_annotation(),
    
    column_title =
      title,
    
    column_title_gp =
      grid::gpar(
        fontsize = 11,
        fontface = "bold"
      ),
    
    cluster_columns = FALSE,
    
    cluster_rows =
      if (
        is.null(row_split)
      ) {
        FALSE
      } else {
        cluster_rows
      },
    
    cluster_row_slices = FALSE,
    
    row_split =
      row_split,
    
    row_title_gp =
      grid::gpar(
        fontsize = 8.5,
        fontface = "bold"
      ),
    
    row_title_rot = 0,
    
    row_labels =
      labels,
    
    show_row_names =
      show_row_names,
    
    row_names_gp =
      grid::gpar(
        fontsize = 6.2
      ),
    
    row_names_max_width =
      grid::unit(
        95,
        "mm"
      ),
    
    column_names_rot = 45,
    
    column_names_gp =
      grid::gpar(
        fontsize = 8
      ),
    
    column_names_centered = TRUE,
    
    rect_gp =
      grid::gpar(
        col = "white",
        lwd = 0.30
      ),
    
    border = TRUE,
    
    heatmap_legend_param = list(
      title = "Expression",
      
      at = c(
        -IMMUNE_Z_LIMIT,
        0,
        IMMUNE_Z_LIMIT
      ),
      
      labels = c(
        paste0(
          "<= -",
          IMMUNE_Z_LIMIT
        ),
        "0",
        paste0(
          ">= ",
          IMMUNE_Z_LIMIT
        )
      )
    )
  )
  
  hm_de <- ComplexHeatmap::Heatmap(
    status_mat,
    name = "DE status",
    
    col = c(
      "Down" =
        unname(
          deg_colors["Down"]
        ),
      
      "NS" =
        unname(
          deg_colors["NS"]
        ),
      
      "Up" =
        unname(
          deg_colors["Up"]
        )
    ),
    
    na_col =
      unname(
        deg_colors["NS"]
      ),
    
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    show_row_names = FALSE,
    
    column_names_rot = 45,
    
    column_names_gp =
      grid::gpar(
        fontsize = 7
      ),
    
    rect_gp =
      grid::gpar(
        col = "white",
        lwd = 0.30
      ),
    
    border = TRUE,
    
    width =
      grid::unit(
        38,
        "mm"
      ),
    
    heatmap_legend_param = list(
      title =
        "Differential\nexpression",
      
      at = c(
        "Up",
        "Down",
        "NS"
      )
    )
  )
  
  hm_expr + hm_de
}


##############################
## 7) Load and map DH13M14 annotation
##############################
if (!file.exists(IMMUNE_ANNOT_FILE)) {
  stop(
    "Annotation file not found: ",
    IMMUNE_ANNOT_FILE
  )
}

if (!file.exists(map_file)) {
  stop(
    "Gene-ID map file not found: ",
    map_file
  )
}

annot <- utils::read.delim(
  IMMUNE_ANNOT_FILE,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  quote = "",
  check.names = FALSE
)

colnames(annot) <- trimws(
  colnames(annot)
)

assert_columns(
  annot,
  c(
    "gene_ID",
    "function"
  ),
  "DH13M14 annotation"
)

annot$annotation_text <-
  build_annotation_text(
    annot
  )

map_df <- utils::read.csv(
  map_file,
  stringsAsFactors = FALSE
) |>
  dplyr::select(
    my_gene_id,
    plant2t_id
  ) |>
  dplyr::mutate(
    my_gene_id =
      trimws(
        as.character(
          my_gene_id
        )
      ),
    
    plant2t_id =
      trimws(
        as.character(
          plant2t_id
        )
      )
  ) |>
  dplyr::filter(
    !is.na(my_gene_id),
    nzchar(my_gene_id),
    !is.na(plant2t_id),
    nzchar(plant2t_id)
  ) |>
  dplyr::distinct()

annot_map_raw <- annot |>
  dplyr::transmute(
    gene_ID =
      trimws(
        as.character(
          gene_ID
        )
      ),
    
    annotation_text =
      as.character(
        annotation_text
      )
  ) |>
  dplyr::inner_join(
    map_df,
    by = c(
      "gene_ID" =
        "plant2t_id"
    )
  ) |>
  dplyr::transmute(
    Gene =
      my_gene_id,
    
    gene_ID,
    annotation_text
  )

# Detect mapping conflicts rather than silently dropping them.
id_conflicts <- annot_map_raw |>
  dplyr::count(
    Gene,
    name = "n_rows"
  ) |>
  dplyr::filter(
    n_rows > 1L
  )

if (
  nrow(id_conflicts) >
  0L
) {
  
  duplicated_identical <-
    annot_map_raw |>
    dplyr::filter(
      Gene %in%
        id_conflicts$Gene
    ) |>
    dplyr::group_by(Gene) |>
    dplyr::summarise(
      n_gene_ID =
        dplyr::n_distinct(
          gene_ID
        ),
      
      .groups = "drop"
    )
  
  serious_conflicts <-
    duplicated_identical |>
    dplyr::filter(
      n_gene_ID > 1L
    )
  
  if (
    nrow(serious_conflicts) >
    0L
  ) {
    
    stop(
      "Multiple DH13M14 gene_IDs map to the same analysis Gene for ",
      nrow(serious_conflicts),
      " genes. Resolve mapping before immune classification."
    )
  }
}

annot_map <- annot_map_raw |>
  dplyr::distinct(
    Gene,
    .keep_all = TRUE
  )

immune_catalogue_all <-
  classify_immune_annotations(
    annot_map
  )

cat(
  "Mapped genes:",
  nrow(annot_map),
  "\n"
)

cat(
  "Genes with immune/microbe-interaction annotation:",
  sum(
    !is.na(
      immune_catalogue_all$Immune_module_auto
    )
  ),
  "\n"
)


##############################
## 8) Load fitted expression object
##############################
pattern_run <- paste0(
  "^",
  SAMPLING_TIME_CHOSEN,
  "_Merged_BothGenotypes_"
)

main_dir <- find_latest_complete_run(
  base_out,
  pattern_run,
  
  required_relpaths = c(
    "01_edgeR_fit",
    "05_DEG_Contrasts/DEG_Tables"
  )
)

dge_file <- file.path(
  main_dir,
  "01_edgeR_fit",
  paste0(
    "dge_MERGED_",
    SAMPLING_TIME_CHOSEN,
    ".rds"
  )
)

if (!file.exists(dge_file)) {
  stop(
    "Merged edgeR DGE object not found: ",
    dge_file
  )
}

y <- readRDS(
  dge_file
)

logcpm <- edgeR::cpm(
  y,
  log = TRUE,
  prior.count = 2
)

meta2 <- y$samples
meta2$Sample <- rownames(meta2)

if (
  !"GenoTemp" %in%
  colnames(meta2)
) {
  
  assert_columns(
    meta2,
    c(
      "Genotype",
      "Temperature"
    ),
    "DGE sample metadata"
  )
  
  meta2$GenoTemp <- paste(
    meta2$Genotype,
    meta2$Temperature,
    sep = "_"
  )
}

missing_groups <- setdiff(
  IMMUNE_GROUP_ORDER,
  unique(
    as.character(
      meta2$GenoTemp
    )
  )
)

if (
  length(missing_groups) >
  0L
) {
  
  stop(
    "Missing required group(s): ",
    paste(
      missing_groups,
      collapse = ", "
    )
  )
}


##############################
## 9) Read ONLY the four intended DE contrasts
##############################
dir_deg <- file.path(
  main_dir,
  "05_DEG_Contrasts",
  "DEG_Tables"
)

read_one_contrast <- function(ct) {
  
  path <- file.path(
    dir_deg,
    paste0(
      "DEG_",
      ct,
      ".csv"
    )
  )
  
  if (!file.exists(path)) {
    stop(
      "Required DE contrast file not found: ",
      path
    )
  }
  
  df <- read_csv_flexible(
    path
  )
  
  if (
    !"Gene" %in%
    colnames(df)
  ) {
    
    if (
      !is.null(
        rownames(df)
      ) &&
      all(
        nzchar(
          rownames(df)
        )
      )
    ) {
      
      df$Gene <-
        rownames(df)
      
    } else {
      
      stop(
        "No Gene column in: ",
        path
      )
    }
  }
  
  assert_columns(
    df,
    c(
      "Gene",
      "logFC",
      "FDR"
    ),
    basename(path)
  )
  
  df |>
    dplyr::transmute(
      Gene =
        trimws(
          as.character(
            Gene
          )
        ),
      
      Contrast =
        ct,
      
      logFC =
        as.numeric(
          logFC
        ),
      
      FDR =
        dplyr::coalesce(
          as.numeric(
            FDR
          ),
          1
        ),
      
      Status =
        status_from_de(
          as.numeric(
            logFC
          ),
          
          dplyr::coalesce(
            as.numeric(
              FDR
            ),
            1
          )
        )
    )
}

all_deg <- lapply(
  IMMUNE_CONTRAST_ORDER,
  read_one_contrast
) |>
  dplyr::bind_rows()


##############################
## 10) Output directories + manual curation
##############################
immune_out <- file.path(
  main_dir,
  paste0(
    "09C_Immune_Microbe_HeatStress_",
    SAMPLING_TIME_CHOSEN,
    "_",
    timestamp
  )
)

dir_tables <- file.path(
  immune_out,
  "00_tables"
)

dir_global <- file.path(
  immune_out,
  "01_global_heatmaps"
)

dir_de <- file.path(
  immune_out,
  "02_DE_heatmaps"
)

dir_module <- file.path(
  immune_out,
  "03_module_heatmaps"
)

dir_qc <- file.path(
  immune_out,
  "04_QC"
)

for (
  d in c(
    immune_out,
    dir_tables,
    dir_global,
    dir_de,
    dir_module,
    dir_qc
  )
) {
  
  dir.create(
    d,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

IMMUNE_OVERRIDE_FILE <- file.path(
  dir_tables,
  "Immune_manual_overrides.csv"
)

immune_catalogue <- immune_catalogue_all |>
  dplyr::filter(
    !is.na(
      Immune_module_auto
    )
  ) |>
  apply_manual_overrides(
    IMMUNE_OVERRIDE_FILE
  )


##############################
## 11) Expression status / expressed immune catalogue
##############################
immune_catalogue <- immune_catalogue |>
  dplyr::mutate(
    In_expression_matrix =
      Gene %in%
      rownames(logcpm)
  )

immune_main <- immune_catalogue |>
  dplyr::filter(
    Include,
    In_expression_matrix
  )

immune_low_or_excluded <- immune_catalogue |>
  dplyr::filter(
    !Include |
      Evidence_level == "Low" |
      !In_expression_matrix
  )

cat(
  "Immune catalogue:",
  nrow(immune_catalogue),
  "\n"
)

cat(
  "Main expressed immune genes:",
  nrow(immune_main),
  "\n"
)


##############################
## 12) Build six-group mean expression
##############################
group_means_all <- vapply(
  IMMUNE_GROUP_ORDER,
  FUN.VALUE =
    numeric(
      nrow(logcpm)
    ),
  
  FUN = function(g) {
    
    samps <- meta2$Sample[
      as.character(
        meta2$GenoTemp
      ) == g
    ]
    
    if (
      length(samps) ==
      0L
    ) {
      
      stop(
        "No samples for group: ",
        g
      )
    }
    
    rowMeans(
      logcpm[
        ,
        samps,
        drop = FALSE
      ],
      na.rm = TRUE
    )
  }
)

rownames(group_means_all) <-
  rownames(logcpm)


##############################
## 13) GLOBAL HEATMAP — all expressed immune genes
##############################
global_immune_ids <- unique(
  immune_main$Gene
)

if (
  length(global_immune_ids) >=
  2L
) {
  
  global_immune_mean <-
    group_means_all[
      global_immune_ids,
      IMMUNE_GROUP_ORDER,
      drop = FALSE
    ]
  
  colnames(global_immune_mean) <-
    IMMUNE_GROUP_LABELS
  
  finite_rows <- apply(
    global_immune_mean,
    1L,
    function(x) {
      all(
        is.finite(x)
      )
    }
  )
  
  global_immune_mean <-
    global_immune_mean[
      finite_rows,
      ,
      drop = FALSE
    ]
  
  global_immune_z <-
    zscore_rows_safe(
      global_immune_mean
    )
  
  ht_global_immune <-
    make_global_expression_heatmap(
      global_immune_z,
      
      title = paste0(
        "Global expression panorama of plant–microbe / immune-associated genes (n=",
        nrow(
          global_immune_z
        ),
        ")"
      )
    )
  
  save_complex_heatmap(
    ht_global_immune,
    
    file.path(
      dir_global,
      "Figure_9C_GLOBAL_ALL_IMMUNE_GENES"
    ),
    
    width_mm =
      IMMUNE_GLOBAL_WIDTH_MM,
    
    height_mm = 200,
    
    allow_raster = TRUE
  )
}


##############################
## 14) GLOBAL HEATMAP — cell-surface receptor candidates
##
## Includes:
##   primary receptors + broad LRR/LysM/Lec/CRK/WAK/CrRLK1L
## Excludes:
##   RLCKs / downstream signalling
##############################
surface_ids <- immune_main |>
  dplyr::filter(
    Immune_module ==
      "01_Cell-surface immune receptors"
  ) |>
  dplyr::pull(
    Gene
  ) |>
  unique()

if (
  length(surface_ids) >=
  2L
) {
  
  surface_mean <-
    group_means_all[
      surface_ids,
      IMMUNE_GROUP_ORDER,
      drop = FALSE
    ]
  
  colnames(surface_mean) <-
    IMMUNE_GROUP_LABELS
  
  surface_z <-
    zscore_rows_safe(
      surface_mean
    )
  
  ht_surface_global <-
    make_global_expression_heatmap(
      surface_z,
      
      title = paste0(
        "Global expression of cell-surface immune receptor candidates (n=",
        nrow(surface_z),
        ")"
      )
    )
  
  save_complex_heatmap(
    ht_surface_global,
    
    file.path(
      dir_global,
      "Figure_9C_GLOBAL_CELL_SURFACE_RECEPTORS"
    ),
    
    width_mm =
      IMMUNE_GLOBAL_WIDTH_MM,
    
    height_mm = 185,
    
    allow_raster = TRUE
  )
}


##############################
## 15) GLOBAL HEATMAP — direct heat × immunity literature priorities
##############################
heat_immune_priority_ids <- immune_main |>
  dplyr::filter(
    Direct_heat_immunity_anchor
  ) |>
  dplyr::pull(
    Gene
  ) |>
  unique()

if (
  length(
    heat_immune_priority_ids
  ) >= 2L
) {
  
  heat_immune_priority_mean <-
    group_means_all[
      heat_immune_priority_ids,
      IMMUNE_GROUP_ORDER,
      drop = FALSE
    ]
  
  colnames(
    heat_immune_priority_mean
  ) <- IMMUNE_GROUP_LABELS
  
  heat_immune_priority_z <-
    zscore_rows_safe(
      heat_immune_priority_mean
    )
  
  ht_heat_immune_priority <-
    make_global_expression_heatmap(
      heat_immune_priority_z,
      
      title = paste0(
        "Heat–immunity intersection: literature-prioritised genes (n=",
        nrow(
          heat_immune_priority_z
        ),
        ")"
      )
    )
  
  save_complex_heatmap(
    ht_heat_immune_priority,
    
    file.path(
      dir_global,
      "Figure_9C_GLOBAL_DIRECT_HEAT_IMMUNITY_INTERSECTION"
    ),
    
    width_mm =
      IMMUNE_GLOBAL_WIDTH_MM,
    
    height_mm = max(
      130,
      min(
        240,
        45 +
          2.5 *
          nrow(
            heat_immune_priority_z
          )
      )
    ),
    
    allow_raster = TRUE
  )
}


##############################
## 16) Select immune genes DE in >=1 intended contrast
##############################
selected_gene_stats <- all_deg |>
  dplyr::filter(
    Gene %in%
      immune_main$Gene,
    
    Contrast %in%
      IMMUNE_CONTRAST_ORDER
  ) |>
  dplyr::group_by(
    Gene
  ) |>
  dplyr::summarise(
    n_DE_contrasts =
      sum(
        Status %in%
          c(
            "Up",
            "Down"
          ),
        na.rm = TRUE
      ),
    
    min_FDR =
      suppressWarnings(
        min(
          FDR[
            Status %in%
              c(
                "Up",
                "Down"
              )
          ],
          na.rm = TRUE
        )
      ),
    
    max_abs_logFC =
      suppressWarnings(
        max(
          abs(
            logFC[
              Status %in%
                c(
                  "Up",
                  "Down"
                )
            ]
          ),
          na.rm = TRUE
        )
      ),
    
    response_direction =
      dplyr::case_when(
        any(
          Status == "Up",
          na.rm = TRUE
        ) &
          any(
            Status == "Down",
            na.rm = TRUE
          ) ~ "Mixed",
        
        any(
          Status == "Up",
          na.rm = TRUE
        ) ~ "Up only",
        
        any(
          Status == "Down",
          na.rm = TRUE
        ) ~ "Down only",
        
        TRUE ~ "NS"
      ),
    
    .groups = "drop"
  ) |>
  dplyr::mutate(
    min_FDR =
      dplyr::if_else(
        is.infinite(
          min_FDR
        ),
        NA_real_,
        min_FDR
      ),
    
    max_abs_logFC =
      dplyr::if_else(
        is.infinite(
          max_abs_logFC
        ),
        NA_real_,
        max_abs_logFC
      )
  ) |>
  dplyr::filter(
    n_DE_contrasts >= 1L
  )

immune_de <- selected_gene_stats |>
  dplyr::inner_join(
    immune_main,
    by = "Gene"
  )

cat(
  "Immune/microbe genes DE in >=1 heat contrast:",
  nrow(immune_de),
  "\n"
)


##############################
## 17) Build DE expression and status matrices
##############################
if (
  nrow(immune_de) >
  0L
) {
  
  de_genes <-
    immune_de$Gene
  
  mean_de <-
    group_means_all[
      de_genes,
      IMMUNE_GROUP_ORDER,
      drop = FALSE
    ]
  
  colnames(mean_de) <-
    IMMUNE_GROUP_LABELS
  
  z_de <-
    zscore_rows_safe(
      mean_de
    )
  
  status_df <- all_deg |>
    dplyr::filter(
      Gene %in% de_genes,
      Contrast %in%
        IMMUNE_CONTRAST_ORDER
    ) |>
    dplyr::select(
      Gene,
      Contrast,
      Status
    ) |>
    dplyr::distinct() |>
    tidyr::complete(
      Gene = de_genes,
      Contrast =
        IMMUNE_CONTRAST_ORDER,
      fill = list(
        Status = "NS"
      )
    ) |>
    tidyr::pivot_wider(
      names_from =
        Contrast,
      values_from =
        Status,
      values_fill =
        "NS"
    ) |>
    dplyr::select(
      Gene,
      dplyr::all_of(
        IMMUNE_CONTRAST_ORDER
      )
    )
  
  status_de <-
    as.matrix(
      status_df[
        ,
        IMMUNE_CONTRAST_ORDER,
        drop = FALSE
      ]
    )
  
  rownames(status_de) <-
    status_df$Gene
  
  colnames(status_de) <-
    IMMUNE_CONTRAST_LABELS
  
  status_de[
    !status_de %in%
      c(
        "Up",
        "Down",
        "NS"
      )
  ] <- "NS"
  
  status_de <-
    status_de[
      de_genes,
      ,
      drop = FALSE
    ]
  
  stopifnot(
    identical(
      rownames(z_de),
      rownames(status_de)
    )
  )
}


##############################
## 18) Response-pattern classification
##############################
if (
  exists("status_de")
) {
  
  status_named <-
    as.data.frame(
      status_de,
      stringsAsFactors = FALSE
    ) |>
    tibble::rownames_to_column(
      "Gene"
    )
  
  response_patterns <- status_named |>
    dplyr::transmute(
      Gene,
      
      ROBILA_HS7_status =
        .data[[
          "ROBILA HS-7"
        ]],
      
      ROBILA_HS2_status =
        .data[[
          "ROBILA HS-2"
        ]],
      
      PRESTO_HS7_status =
        .data[[
          "PRESTO HS-7"
        ]],
      
      PRESTO_HS2_status =
        .data[[
          "PRESTO HS-2"
        ]],
      
      ROBILA_pattern =
        response_pattern_two(
          ROBILA_HS7_status,
          ROBILA_HS2_status
        ),
      
      PRESTO_pattern =
        response_pattern_two(
          PRESTO_HS7_status,
          PRESTO_HS2_status
        ),
      
      Cross_genotype_pattern =
        paste0(
          "ROBILA: ",
          ROBILA_pattern,
          " | PRESTO: ",
          PRESTO_pattern
        )
    )
  
  immune_de <- immune_de |>
    dplyr::left_join(
      response_patterns,
      by = "Gene"
    )
} else {
  
  response_patterns <-
    tibble::tibble()
}


##############################
## 19) Module census
##############################
module_universe <-
  immune_main |>
  dplyr::distinct(
    Gene,
    Immune_module
  )

module_de <-
  immune_de |>
  dplyr::distinct(
    Gene,
    Immune_module
  )

module_up <-
  all_deg |>
  dplyr::filter(
    Gene %in%
      immune_de$Gene,
    Status == "Up"
  ) |>
  dplyr::left_join(
    immune_main |>
      dplyr::select(
        Gene,
        Immune_module
      ),
    by = "Gene"
  ) |>
  dplyr::distinct(
    Gene,
    Immune_module
  )

module_down <-
  all_deg |>
  dplyr::filter(
    Gene %in%
      immune_de$Gene,
    Status == "Down"
  ) |>
  dplyr::left_join(
    immune_main |>
      dplyr::select(
        Gene,
        Immune_module
      ),
    by = "Gene"
  ) |>
  dplyr::distinct(
    Gene,
    Immune_module
  )

module_census <-
  module_universe |>
  dplyr::count(
    Immune_module,
    name = "n_expressed"
  ) |>
  dplyr::full_join(
    module_de |>
      dplyr::count(
        Immune_module,
        name = "n_DE"
      ),
    by = "Immune_module"
  ) |>
  dplyr::full_join(
    module_up |>
      dplyr::count(
        Immune_module,
        name = "n_Up"
      ),
    by = "Immune_module"
  ) |>
  dplyr::full_join(
    module_down |>
      dplyr::count(
        Immune_module,
        name = "n_Down"
      ),
    by = "Immune_module"
  ) |>
  dplyr::mutate(
    dplyr::across(
      dplyr::starts_with("n_"),
      ~ tidyr::replace_na(
        .x,
        0L
      )
    ),
    
    pct_DE_of_expressed =
      dplyr::if_else(
        n_expressed > 0,
        100 *
          n_DE /
          n_expressed,
        NA_real_
      ),
    
    Immune_module =
      factor(
        Immune_module,
        levels =
          IMMUNE_MODULE_ORDER
      )
  ) |>
  dplyr::arrange(
    Immune_module
  )


##############################
## 20) Module overrepresentation among heat-responsive genes
##############################
n_immune_expressed <-
  dplyr::n_distinct(
    module_universe$Gene
  )

n_immune_de <-
  dplyr::n_distinct(
    module_de$Gene
  )

module_enrichment <-
  module_census |>
  dplyr::filter(
    n_expressed > 0
  ) |>
  dplyr::rowwise() |>
  dplyr::mutate(
    a_DE_in_module =
      n_DE,
    
    b_nonDE_in_module =
      n_expressed -
      n_DE,
    
    c_DE_other =
      n_immune_de -
      n_DE,
    
    d_nonDE_other =
      (
        n_immune_expressed -
          n_expressed
      ) -
      c_DE_other,
    
    Fisher_p = {
      
      tab <- matrix(
        c(
          a_DE_in_module,
          b_nonDE_in_module,
          c_DE_other,
          d_nonDE_other
        ),
        nrow = 2L,
        byrow = TRUE
      )
      
      if (
        all(
          tab >= 0
        )
      ) {
        
        stats::fisher.test(
          tab,
          alternative =
            "greater"
        )$p.value
        
      } else {
        
        NA_real_
      }
    },
    
    Odds_ratio = {
      
      tab <- matrix(
        c(
          a_DE_in_module,
          b_nonDE_in_module,
          c_DE_other,
          d_nonDE_other
        ),
        nrow = 2L,
        byrow = TRUE
      )
      
      if (
        all(
          tab >= 0
        )
      ) {
        
        unname(
          stats::fisher.test(
            tab,
            alternative =
              "greater"
          )$estimate
        )
        
      } else {
        
        NA_real_
      }
    }
  ) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    Fisher_FDR =
      stats::p.adjust(
        Fisher_p,
        method = "BH"
      ),
    
    Enriched_FDR05 =
      !is.na(
        Fisher_FDR
      ) &
      Fisher_FDR <
      0.05
  ) |>
  dplyr::arrange(
    Fisher_FDR,
    dplyr::desc(
      Odds_ratio
    )
  )


##############################
## 21) Reference audit against Liu et al. 2024
##
## IMPORTANT:
## Their published counts were obtained using protein-domain
## HMMER + TMHMM, not free-text annotation.
## This audit tells us whether our local annotation-text approach
## is grossly incomplete for these structural classes.
##############################
observed_lrr <- immune_catalogue_all |>
  dplyr::filter(
    hit_LRR_RLK
  ) |>
  dplyr::summarise(
    n = dplyr::n_distinct(
      Gene
    )
  ) |>
  dplyr::pull(n)

observed_lysm <- immune_catalogue_all |>
  dplyr::filter(
    hit_LysM_RLK |
      hit_CERK1 |
      hit_LYK
  ) |>
  dplyr::summarise(
    n = dplyr::n_distinct(
      Gene
    )
  ) |>
  dplyr::pull(n)

observed_nlr <- immune_catalogue_all |>
  dplyr::filter(
    nlr_module
  ) |>
  dplyr::summarise(
    n = dplyr::n_distinct(
      Gene
    )
  ) |>
  dplyr::pull(n)

observed_tnl <- immune_catalogue_all |>
  dplyr::filter(
    hit_TNL
  ) |>
  dplyr::summarise(
    n = dplyr::n_distinct(
      Gene
    )
  ) |>
  dplyr::pull(n)

observed_cnl <- immune_catalogue_all |>
  dplyr::filter(
    hit_CNL
  ) |>
  dplyr::summarise(
    n = dplyr::n_distinct(
      Gene
    )
  ) |>
  dplyr::pull(n)

observed_rnl <- immune_catalogue_all |>
  dplyr::filter(
    hit_RNL
  ) |>
  dplyr::summarise(
    n = dplyr::n_distinct(
      Gene
    )
  ) |>
  dplyr::pull(n)

reference_audit <-
  LIU2024_REFERENCE_COUNTS |>
  dplyr::mutate(
    Local_annotation_text_n =
      dplyr::case_when(
        Reference_class == "LRR-RLK" ~
          observed_lrr,
        
        Reference_class == "LysM-RLK" ~
          observed_lysm,
        
        Reference_class == "NLR_total" ~
          observed_nlr,
        
        Reference_class == "CNL" ~
          observed_cnl,
        
        Reference_class == "TNL" ~
          observed_tnl,
        
        Reference_class == "RNL" ~
          observed_rnl,
        
        Reference_class == "Other_NLR" ~
          max(
            0L,
            observed_nlr -
              observed_cnl -
              observed_tnl -
              observed_rnl
          ),
        
        TRUE ~
          NA_integer_
      ),
    
    Recovery_pct =
      100 *
      Local_annotation_text_n /
      Published_n,
    
    Interpretation =
      dplyr::case_when(
        Recovery_pct >= 80 ~
          "Good recovery by annotation text.",
        
        Recovery_pct >= 50 ~
          "Partial recovery; manual/domain validation recommended.",
        
        TRUE ~
          "Low recovery: annotation-text search is incomplete; domain-based HMMER/TMHMM or Liu et al. Table S15 should be preferred for exhaustive structural census."
      )
  )

if (
  any(
    reference_audit$Recovery_pct[
      reference_audit$Reference_class %in%
      c(
        "LRR-RLK",
        "LysM-RLK",
        "NLR_total"
      )
    ] <
    50,
    na.rm = TRUE
  )
) {
  
  warning(
    "The text-annotation classifier recovered <50% of at least one ",
    "published Liu et al. 2024 structural immune class. ",
    "Use this script for expression profiling of annotated immune genes, ",
    "but do not claim an exhaustive RLK/NLR genome census without ",
    "domain-based validation or Table S15."
  )
}


##############################
## 22) COMPLETE DE immune-network heatmap
##############################
if (
  nrow(immune_de) >=
  2L
) {
  
  module_levels_present <-
    IMMUNE_MODULE_ORDER[
      IMMUNE_MODULE_ORDER %in%
        immune_de$Immune_module
    ]
  
  immune_de <- immune_de |>
    dplyr::mutate(
      Immune_module =
        factor(
          Immune_module,
          levels =
            module_levels_present
        )
    ) |>
    dplyr::arrange(
      Immune_module,
      Immune_submodule,
      Gene
    )
  
  de_genes <-
    immune_de$Gene
  
  mean_de <-
    mean_de[
      de_genes,
      ,
      drop = FALSE
    ]
  
  z_de <-
    z_de[
      de_genes,
      ,
      drop = FALSE
    ]
  
  status_de <-
    status_de[
      de_genes,
      ,
      drop = FALSE
    ]
  
  labels_all <-
    stats::setNames(
      paste0(
        immune_de$Gene,
        " | ",
        clean_label(
          immune_de$Short_name,
          52L
        )
      ),
      immune_de$Gene
    )
  
  split_modules <-
    droplevels(
      immune_de$Immune_module
    )
  
  ht_complete <-
    make_heatmap_pair(
      zmat =
        z_de,
      
      status_mat =
        status_de,
      
      labels =
        labels_all,
      
      title = paste0(
        "Figure 9C | Heat-responsive plant–microbe / immune network (n=",
        nrow(
          immune_de
        ),
        ")"
      ),
      
      row_split =
        split_modules,
      
      cluster_rows =
        IMMUNE_CLUSTER_ROWS,
      
      show_row_names =
        TRUE
    )
  
  complete_height_mm <-
    max(
      190,
      60 +
        IMMUNE_ROW_MM *
        nrow(
          immune_de
        )
    )
  
  save_complex_heatmap(
    ht_complete,
    
    file.path(
      dir_de,
      "Figure_9C_COMPLETE_DE_IMMUNE_NETWORK"
    ),
    
    width_mm =
      IMMUNE_FIG_WIDTH_MM +
      25,
    
    height_mm =
      complete_height_mm,
    
    allow_raster = TRUE
  )
}


##############################
## 23) DE cell-surface receptor heatmap
##############################
de_surface <- immune_de |>
  dplyr::filter(
    as.character(
      Immune_module
    ) ==
      "01_Cell-surface immune receptors"
  )

if (
  nrow(de_surface) >=
  2L
) {
  
  genes <-
    de_surface$Gene
  
  lab <-
    stats::setNames(
      paste0(
        genes,
        " | ",
        clean_label(
          de_surface$Short_name,
          55L
        )
      ),
      genes
    )
  
  ht <-
    make_heatmap_pair(
      zmat =
        z_de[
          genes,
          ,
          drop = FALSE
        ],
      
      status_mat =
        status_de[
          genes,
          ,
          drop = FALSE
        ],
      
      labels =
        lab,
      
      title = paste0(
        "Cell-surface immune receptor candidates — heat responsive (n=",
        length(genes),
        ")"
      ),
      
      row_split =
        NULL,
      
      cluster_rows =
        TRUE,
      
      show_row_names =
        TRUE
    )
  
  save_complex_heatmap(
    ht,
    
    file.path(
      dir_de,
      "Figure_9C_DE_CELL_SURFACE_RECEPTORS"
    ),
    
    width_mm =
      IMMUNE_FIG_WIDTH_MM,
    
    height_mm =
      max(
        110,
        45 +
          IMMUNE_ROW_MM *
          length(genes)
      ),
    
    allow_raster = TRUE
  )
}


##############################
## 24) DE NLR / ETI heatmap
##############################
de_nlr_eti <- immune_de |>
  dplyr::filter(
    as.character(
      Immune_module
    ) %in%
      c(
        "07_Intracellular NLR receptors",
        "08_ETI support / guardee machinery"
      )
  )

if (
  nrow(de_nlr_eti) >=
  2L
) {
  
  genes <-
    de_nlr_eti$Gene
  
  de_nlr_eti <-
    de_nlr_eti |>
    dplyr::mutate(
      Immune_module =
        droplevels(
          Immune_module
        )
    )
  
  lab <-
    stats::setNames(
      paste0(
        genes,
        " | ",
        clean_label(
          de_nlr_eti$Short_name,
          55L
        )
      ),
      genes
    )
  
  ht <-
    make_heatmap_pair(
      zmat =
        z_de[
          genes,
          ,
          drop = FALSE
        ],
      
      status_mat =
        status_de[
          genes,
          ,
          drop = FALSE
        ],
      
      labels =
        lab,
      
      title = paste0(
        "Intracellular NLR / ETI machinery — heat responsive (n=",
        length(genes),
        ")"
      ),
      
      row_split =
        droplevels(
          de_nlr_eti$Immune_module
        ),
      
      cluster_rows =
        TRUE,
      
      show_row_names =
        TRUE
    )
  
  save_complex_heatmap(
    ht,
    
    file.path(
      dir_de,
      "Figure_9C_DE_NLR_ETI"
    ),
    
    width_mm =
      IMMUNE_FIG_WIDTH_MM,
    
    height_mm =
      max(
        120,
        50 +
          IMMUNE_ROW_MM *
          length(genes)
      ),
    
    allow_raster = TRUE
  )
}


##############################
## 25) DE early PTI heatmap
##############################
early_pti_modules <- c(
  "01_Cell-surface immune receptors",
  "02_Co-receptors / receptor scaffolds",
  "03_Receptor-proximal RLCKs",
  "04_Ca2+ signalling",
  "05_ROS / redox signalling",
  "06_MAPK signalling"
)

de_early_pti <- immune_de |>
  dplyr::filter(
    as.character(
      Immune_module
    ) %in%
      early_pti_modules
  )

if (
  nrow(de_early_pti) >=
  2L
) {
  
  genes <-
    de_early_pti$Gene
  
  de_early_pti <-
    de_early_pti |>
    dplyr::mutate(
      Immune_module =
        factor(
          as.character(
            Immune_module
          ),
          levels =
            early_pti_modules
        )
    ) |>
    dplyr::arrange(
      Immune_module,
      Immune_submodule,
      Gene
    )
  
  genes <-
    de_early_pti$Gene
  
  lab <-
    stats::setNames(
      paste0(
        genes,
        " | ",
        clean_label(
          de_early_pti$Short_name,
          52L
        )
      ),
      genes
    )
  
  ht <-
    make_heatmap_pair(
      zmat =
        z_de[
          genes,
          ,
          drop = FALSE
        ],
      
      status_mat =
        status_de[
          genes,
          ,
          drop = FALSE
        ],
      
      labels =
        lab,
      
      title = paste0(
        "Early PTI / cell-surface signalling — heat responsive (n=",
        length(genes),
        ")"
      ),
      
      row_split =
        droplevels(
          de_early_pti$Immune_module
        ),
      
      cluster_rows =
        TRUE,
      
      show_row_names =
        TRUE
    )
  
  save_complex_heatmap(
    ht,
    
    file.path(
      dir_de,
      "Figure_9C_DE_EARLY_PTI"
    ),
    
    width_mm =
      IMMUNE_FIG_WIDTH_MM +
      10,
    
    height_mm =
      max(
        140,
        55 +
          IMMUNE_ROW_MM *
          length(genes)
      ),
    
    allow_raster = TRUE
  )
}


##############################
## 26) One DE heatmap per immune module
##############################
if (
  nrow(immune_de) >
  0L
) {
  
  modules_present <-
    unique(
      as.character(
        immune_de$Immune_module
      )
    )
  
  modules_present <-
    IMMUNE_MODULE_ORDER[
      IMMUNE_MODULE_ORDER %in%
        modules_present
    ]
  
  for (
    module_name in modules_present
  ) {
    
    module_tbl <-
      immune_de |>
      dplyr::filter(
        as.character(
          Immune_module
        ) ==
          module_name
      ) |>
      dplyr::arrange(
        Immune_submodule,
        Gene
      )
    
    genes <-
      module_tbl$Gene
    
    if (
      length(genes) <
      1L
    ) {
      next
    }
    
    lab <-
      stats::setNames(
        paste0(
          genes,
          " | ",
          clean_label(
            module_tbl$Short_name,
            55L
          )
        ),
        genes
      )
    
    # A one-gene module can still be plotted without clustering.
    ht <-
      make_heatmap_pair(
        zmat =
          z_de[
            genes,
            ,
            drop = FALSE
          ],
        
        status_mat =
          status_de[
            genes,
            ,
            drop = FALSE
          ],
        
        labels =
          lab,
        
        title = paste0(
          module_name,
          " — DE in heat stress (n=",
          length(genes),
          ")"
        ),
        
        row_split =
          NULL,
        
        cluster_rows =
          length(genes) > 2L,
        
        show_row_names =
          TRUE
      )
    
    save_complex_heatmap(
      ht,
      
      file.path(
        dir_module,
        paste0(
          "Figure_9C_",
          safe_file_stub(
            module_name
          ),
          "_DE"
        )
      ),
      
      width_mm =
        IMMUNE_FIG_WIDTH_MM,
      
      height_mm =
        max(
          95,
          42 +
            IMMUNE_ROW_MM *
            length(genes)
        ),
      
      allow_raster = TRUE
    )
  }
}


##############################
## 27) Literature-priority table
##############################
literature_rules <- tibble::tribble(
  ~Biological_layer,
  ~Genes_or_families,
  ~Why_included,
  ~Heat_intersection,
  ~Reference,
  
  "Carrot-specific immune repertoire",
  "LRR-RLK; LysM-RLK; NLR; CNGC; CDPK/CPK; CaM/CML; RBOH; OXI1; NDPK2; WRKY22; ACS6; PR1; EDS1; PBS1; RIN4; RPS2; SGT1; HSP90",
  "Directly examined in carrot during Alternaria dauci infection.",
  "These genes are profiled here under heat alone to test whether thermal stress shifts basal immune competence.",
  "Liu et al. 2024, The Plant Journal, DOI:10.1111/tpj.17049",
  
  "Flagellin perception",
  "FLS2; BAK1/SERK3",
  "Canonical bacterial-pattern receptor/co-receptor complex.",
  "Temporary heat suppresses flg22/PTI outputs; transcript changes are therefore biologically relevant but do not measure receptor activity.",
  "Sun et al. 2013 Science DOI:10.1126/science.1243825; Mühlenbeck et al. 2019 PMID:30924595",
  
  "Fungal chitin perception",
  "LYK5; CERK1; LysM-RLK/RLP",
  "Major fungal MAMP receptor system.",
  "Especially relevant to Alternaria because chitin is a fungal cell-wall MAMP.",
  "Cao et al. 2014 eLife DOI:10.7554/eLife.03766",
  
  "RLP immune complexes",
  "RLPs; SOBIR1; BAK1",
  "RLPs require receptor kinases such as SOBIR1 and often BAK1 for signalling.",
  "SOBIR1 has direct evidence for high-temperature/biotic-response reprioritisation in potato.",
  "Liebrand et al. 2013 PMID:23716655; StSOBIR1 2025 PMID:40891209",
  
  "FERONIA",
  "FER / CrRLK1L",
  "FER regulates immune receptor complex assembly and plant–microbe interactions.",
  "FER is now directly implicated in plasma-membrane thermal sensing/thermotolerance.",
  "Stegmann et al. 2017 Science DOI:10.1126/science.aal2541; Science 2026 DOI:10.1126/science.aeb1752",
  
  "Heat receptor complex",
  "BAM1; PBS1; RBOHD",
  "A plasma-membrane BAM1-PBS1 complex activates RBOHD-dependent H2O2 during heat.",
  "Strong direct molecular bridge between heat perception and components also used in immune signalling.",
  "Molecular Plant 2025 DOI:10.1016/j.molp.2025.10.021",
  
  "ROS / MAPK",
  "RBOHD/F; OXI1; MAPKKK-MKK-MPK",
  "Core early PTI outputs and signal propagation.",
  "OXI1 links MAMP-induced ROS to MAPK; crop MAPKs can coordinate heat tolerance and pathogen resistance.",
  "Ma et al. 2024 Plant Cell PMID:39566103; CaMAPK1 2024 DOI:10.3390/plants13131775",
  
  "Cysteine-rich RLKs",
  "CRK / DUF26 RLKs",
  "Large stress-responsive RLK family with immune functions.",
  "Direct crop evidence links CRKs to both heat and pathogen responses.",
  "TaCRK10 DOI:10.1111/tpj.15513; rice CRK cluster DOI:10.1111/pbi.14381",
  
  "Lectin RLKs",
  "LecRLK / G-type lectin RLK",
  "Several LecRLKs participate in microbial perception/defence.",
  "Rice PWL1 demonstrates heat-tolerance versus bacterial-resistance trade-off.",
  "PWL1 2023 DOI:10.1111/pbi.14150",
  
  "SA immunity",
  "CBP60g; SARD1; ICS1; EDS1; PAD4; PBS3; EDS5; NPR",
  "Central basal/ETI defence transcription and SA biosynthesis/signalling.",
  "Elevated temperature suppresses this defence sector upstream at CBP60g/SARD1 regulation.",
  "Kim et al. 2022 Nature DOI:10.1038/s41586-022-04902-y",
  
  "NHP / SAR",
  "ALD1; FMO1",
  "Core NHP biosynthesis / systemic acquired resistance machinery.",
  "Warm temperature suppresses SAR by limiting NHP biosynthesis.",
  "2025 PMID:40758777",
  
  "NLR / ETI",
  "NLR; EDS1; PAD4; SAG101; ADR1; NRG1; NDR1; RIN4; SGT1; RAR1; HSP90",
  "Intracellular effector recognition and ETI support.",
  "ETI outputs and several NLR-dependent resistances are temperature sensitive.",
  "Carrot Liu et al. 2024 DOI:10.1111/tpj.17049; elevated-temperature ETI study PMID:26617631",
  
  "Cell-wall/DAMP perception",
  "WAK/WAKL; PEPR1/2; MIK2; P2K1/DORN1",
  "Damage and cell-wall integrity sensing potentiates immunity.",
  "Relevant because heat can perturb membranes/cell walls and could alter danger signalling even without infection.",
  "PEPR primary studies PMID:20200150; MIK2/SCOOP18 2024 PMID:39482527"
)

readr::write_csv(
  literature_rules,
  file.path(
    dir_tables,
    "Immune_literature_rules.csv"
  )
)


##############################
## 28) Export tables
##############################
readr::write_csv(
  immune_catalogue,
  file.path(
    dir_tables,
    "Immune_gene_catalogue_complete.csv"
  )
)

readr::write_csv(
  immune_main,
  file.path(
    dir_tables,
    "Immune_gene_catalogue_main.csv"
  )
)

readr::write_csv(
  immune_low_or_excluded,
  file.path(
    dir_qc,
    "Immune_excluded_lowconfidence_or_missing.csv"
  )
)

readr::write_csv(
  immune_de,
  file.path(
    dir_tables,
    "Immune_DE_selected.csv"
  )
)

readr::write_csv(
  module_census,
  file.path(
    dir_tables,
    "Immune_module_census.csv"
  )
)

readr::write_csv(
  module_enrichment,
  file.path(
    dir_tables,
    "Immune_module_overrepresentation_Fisher.csv"
  )
)

readr::write_csv(
  response_patterns,
  file.path(
    dir_tables,
    "Immune_response_patterns.csv"
  )
)

readr::write_csv(
  reference_audit,
  file.path(
    dir_tables,
    "Immune_carrot_Liu2024_reference_audit.csv"
  )
)

# Global expression matrices.
if (
  exists(
    "global_immune_mean"
  )
) {
  
  readr::write_csv(
    tibble::rownames_to_column(
      as.data.frame(
        global_immune_mean
      ),
      "Gene"
    ) |>
      dplyr::left_join(
        immune_main |>
          dplyr::select(
            Gene,
            Immune_module,
            Immune_submodule,
            Evidence_level,
            Carrot_Adauci_anchor,
            Direct_heat_immunity_anchor
          ),
        by = "Gene"
      ),
    
    file.path(
      dir_tables,
      "Immune_GLOBAL_ALL_meanLogCPM.csv"
    )
  )
  
  readr::write_csv(
    tibble::rownames_to_column(
      as.data.frame(
        global_immune_z
      ),
      "Gene"
    ) |>
      dplyr::left_join(
        immune_main |>
          dplyr::select(
            Gene,
            Immune_module,
            Immune_submodule,
            Evidence_level,
            Carrot_Adauci_anchor,
            Direct_heat_immunity_anchor
          ),
        by = "Gene"
      ),
    
    file.path(
      dir_tables,
      "Immune_GLOBAL_ALL_zscore.csv"
    )
  )
}

if (
  exists("mean_de")
) {
  
  readr::write_csv(
    tibble::rownames_to_column(
      as.data.frame(
        mean_de
      ),
      "Gene"
    ) |>
      dplyr::left_join(
        immune_de |>
          dplyr::select(
            Gene,
            Immune_module,
            Immune_submodule,
            Evidence_level,
            n_DE_contrasts,
            min_FDR,
            max_abs_logFC,
            response_direction
          ),
        by = "Gene"
      ),
    
    file.path(
      dir_tables,
      "Immune_DE_meanLogCPM.csv"
    )
  )
  
  readr::write_csv(
    tibble::rownames_to_column(
      as.data.frame(
        z_de
      ),
      "Gene"
    ) |>
      dplyr::left_join(
        immune_de |>
          dplyr::select(
            Gene,
            Immune_module,
            Immune_submodule,
            Evidence_level,
            n_DE_contrasts,
            response_direction
          ),
        by = "Gene"
      ),
    
    file.path(
      dir_tables,
      "Immune_DE_zscore.csv"
    )
  )
  
  readr::write_csv(
    tibble::rownames_to_column(
      as.data.frame(
        status_de
      ),
      "Gene"
    ),
    
    file.path(
      dir_tables,
      "Immune_DE_status_four_contrasts.csv"
    )
  )
}


##############################
## 29) Per-contrast module counts
##############################
module_per_contrast <- all_deg |>
  dplyr::filter(
    Gene %in%
      immune_main$Gene,
    Status %in%
      c(
        "Up",
        "Down"
      )
  ) |>
  dplyr::left_join(
    immune_main |>
      dplyr::select(
        Gene,
        Immune_module,
        Immune_submodule
      ),
    by = "Gene"
  ) |>
  dplyr::distinct(
    Gene,
    Contrast,
    Immune_module,
    Status
  ) |>
  dplyr::count(
    Immune_module,
    Contrast,
    Status,
    name = "n"
  ) |>
  dplyr::left_join(
    module_census |>
      dplyr::select(
        Immune_module,
        n_expressed
      ),
    by = "Immune_module"
  ) |>
  dplyr::mutate(
    pct_of_expressed_module =
      dplyr::if_else(
        n_expressed > 0,
        100 *
          n /
          n_expressed,
        NA_real_
      )
  )

readr::write_csv(
  module_per_contrast,
  file.path(
    dir_tables,
    "Immune_module_per_contrast.csv"
  )
)


##############################
## 30) QC / reproducibility
##############################
utils::write.csv(
  meta2,
  file.path(
    dir_qc,
    "metadata_from_DGE_used.csv"
  ),
  row.names = FALSE
)

readr::write_csv(
  annot_map,
  file.path(
    dir_qc,
    "annotation_mapped_all_genes.csv"
  )
)

readr::write_csv(
  immune_catalogue_all,
  file.path(
    dir_qc,
    "immune_classifier_all_mapped_genes.csv"
  )
)

readr::write_csv(
  all_deg,
  file.path(
    dir_qc,
    "all_four_heat_contrasts.csv"
  )
)

grand_summary <- tibble::tibble(
  Sampling_time =
    SAMPLING_TIME_CHOSEN,
  
  n_contrasts =
    length(
      IMMUNE_CONTRAST_ORDER
    ),
  
  n_mapped_genes =
    nrow(
      annot_map
    ),
  
  n_immune_candidates_total =
    nrow(
      immune_catalogue
    ),
  
  n_immune_main_expressed =
    nrow(
      immune_main
    ),
  
  n_immune_DE =
    nrow(
      immune_de
    ),
  
  n_modules_expressed =
    dplyr::n_distinct(
      immune_main$Immune_module
    ),
  
  n_modules_with_DE =
    dplyr::n_distinct(
      immune_de$Immune_module
    ),
  
  n_surface_receptor_candidates_expressed =
    length(
      surface_ids
    ),
  
  n_direct_heat_immunity_priority_expressed =
    length(
      heat_immune_priority_ids
    ),
  
  Liu2024_LRR_RLK_reference =
    283L,
  
  Liu2024_LysM_RLK_reference =
    10L,
  
  Liu2024_NLR_reference =
    202L
)

readr::write_csv(
  grand_summary,
  file.path(
    dir_qc,
    "Immune_GRAND_SUMMARY.csv"
  )
)

saveRDS(
  list(
    main_dir =
      main_dir,
    
    annotation_file =
      IMMUNE_ANNOT_FILE,
    
    annotation_mapped =
      annot_map,
    
    classifier_all =
      immune_catalogue_all,
    
    immune_catalogue =
      immune_catalogue,
    
    immune_main =
      immune_main,
    
    metadata =
      meta2,
    
    DGE =
      y,
    
    logCPM =
      logcpm,
    
    group_means_all =
      group_means_all,
    
    four_contrasts =
      all_deg,
    
    immune_DE =
      immune_de,
    
    DE_mean_logCPM =
      get0(
        "mean_de",
        ifnotfound = NULL
      ),
    
    DE_zscore =
      get0(
        "z_de",
        ifnotfound = NULL
      ),
    
    DE_status =
      get0(
        "status_de",
        ifnotfound = NULL
      ),
    
    global_immune_mean =
      get0(
        "global_immune_mean",
        ifnotfound = NULL
      ),
    
    global_immune_z =
      get0(
        "global_immune_z",
        ifnotfound = NULL
      ),
    
    response_patterns =
      response_patterns,
    
    module_census =
      module_census,
    
    module_enrichment =
      module_enrichment,
    
    Liu2024_reference_audit =
      reference_audit,
    
    literature =
      literature_rules
  ),
  
  file.path(
    dir_qc,
    "Figure_9C_Immune_analysis_objects.rds"
  )
)

capture.output(
  utils::sessionInfo(),
  file =
    file.path(
      immune_out,
      "sessionInfo.txt"
    )
)


##############################
## 31) Final console summary
##############################
cat(
  "\nDONE Figure 9C immune / microbe heat-stress pipeline\n"
)

cat(
  "Input run:",
  main_dir,
  "\n"
)

cat(
  "Mapped genes:",
  nrow(
    annot_map
  ),
  "\n"
)

cat(
  "Immune candidates:",
  nrow(
    immune_catalogue
  ),
  "\n"
)

cat(
  "Main expressed immune genes:",
  nrow(
    immune_main
  ),
  "\n"
)

cat(
  "Heat-responsive immune genes:",
  nrow(
    immune_de
  ),
  "\n"
)

cat(
  "Cell-surface receptor candidates expressed:",
  length(
    surface_ids
  ),
  "\n"
)

cat(
  "Direct heat-immunity priority genes expressed:",
  length(
    heat_immune_priority_ids
  ),
  "\n"
)

cat(
  "Expression order:",
  paste(
    IMMUNE_GROUP_LABELS,
    collapse = " | "
  ),
  "\n"
)

cat(
  "Contrast order:",
  paste(
    IMMUNE_CONTRAST_LABELS,
    collapse = " | "
  ),
  "\n"
)

cat(
  "Output directory:",
  immune_out,
  "\n"
)
