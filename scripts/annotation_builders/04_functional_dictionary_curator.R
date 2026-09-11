# NOTE: companion annotation-builder reference implementation.
# Paths have been sanitized. See docs/05_functional_dictionaries.md before running.

############################################################
## 10_Functional_KEGG_family_heatmaps_LITERATURE_REFINED.R
## Differentially expressed functional families defined from
## KEGG-enhanced pathway sets
##
## Input root:
##   10_Functional_Families_KEGG_Enhanced/
##
## Expected within each family directory:
##   <family>_expression.csv
##   <family>_DE_per_contrast.csv
##   <family>_subfamily_assignment.csv
##   <family>_KEGG_provenance.csv
##
## Selection:
##   - retain a gene when Status is Up or Down in >= 1 contrast
##
## Expression panel:
##   - six precomputed mean-logCPM groups
##   - row-wise z-score
##   - fixed biological order:
##       ROBILA NH -> HS-7 -> HS-2 ->
##       PRESTO NH -> HS-7 -> HS-2
##
## Functional classification — literature-refined version:
##   - explicit gene/enzyme annotation + KEGG provenance are evaluated for
##     EVERY selected gene, not only genes previously labelled "Other"
##   - supplied subfamily assignments are retained only as a secondary
##     fallback when no stronger literature-guided rule is available
##   - unresolved genes are exported for manual curation but are NOT shown
##     as an "Other" category in the figures
##   - subfamilies are displayed in biologically meaningful pathway order
##     rather than alphabetical order
##   - manual overrides are supported
##
## Outputs:
##   - one ComplexHeatmap figure per functional family
##   - optional combined figure
##   - expression z-score + Up/Down/NS panel
##   - PDF, SVG, PNG and TIFF
##   - curation, overlap, provenance and bibliography tables
############################################################

##############################
## 0) Configuration
##############################
source("00_builder_config.R")
source("00_builder_helpers.R")

FUNCTIONAL_FAMILY_DIRS <- c(
  "01_Protein_Folding_Processing_Degradation",
  "02_Hormone_Signaling_Metabolism",
  "03_MAPK_Calcium_Stress_Signaling",
  "04_Photosynthesis_Energy_Metabolism",
  "05_Terpenoid_Secondary_Metabolism",
  "06_ROS_Redox_Glutathione",
  "07_Lipid_Membrane_Metabolism",
  "08_Carbon_Starch_Sucrose_Osmoprotection"
)

FUNCTIONAL_FAMILY_LABELS <- c(
  "01_Protein_Folding_Processing_Degradation" =
    "Protein folding, processing and degradation",
  "02_Hormone_Signaling_Metabolism" =
    "Hormone signalling and metabolism",
  "03_MAPK_Calcium_Stress_Signaling" =
    "MAPK, calcium and stress signalling",
  "04_Photosynthesis_Energy_Metabolism" =
    "Photosynthesis and energy metabolism",
  "05_Terpenoid_Secondary_Metabolism" =
    "Terpenoid and specialised metabolism",
  "06_ROS_Redox_Glutathione" =
    "ROS, redox and glutathione metabolism",
  "07_Lipid_Membrane_Metabolism" =
    "Lipid and membrane metabolism",
  "08_Carbon_Starch_Sucrose_Osmoprotection" =
    "Carbon, starch, sucrose and osmoprotection"
)

# Fixed input path supplied by the user.
# It may still be overridden from the shell:
# export FUNCTIONAL_INPUT_DIR=/another/path
FUNCTIONAL_INPUT_DIR <- FUNCTIONAL_INPUT_DIR

FUNCTIONAL_OVERRIDE_FILE <- file.path(
  FUNCTIONAL_INPUT_DIR,
  "Functional_manual_overrides.csv"
)

# Differential-expression selection.
USE_EXISTING_DE_STATUS <- TRUE
FUNCTIONAL_FDR_CUTOFF <- PADJ_CUTOFF
FUNCTIONAL_LFC_CUTOFF <- lfc_threshold

# Curation settings.
# The rule-based classifier is now applied to ALL genes.
# "Other" is never used as a plotted biological category.
UNRESOLVED_LABEL <- "Unresolved / manual curation"
INCLUDE_UNRESOLVED_IN_PLOTS <- FALSE

# Broad-family display order follows a biological story:
# signalling -> proteostasis/redox/membrane protection ->
# photosynthesis/carbon -> specialised metabolism.
FUNCTIONAL_PLOT_ORDER <- c(
  "02_Hormone_Signaling_Metabolism",
  "03_MAPK_Calcium_Stress_Signaling",
  "01_Protein_Folding_Processing_Degradation",
  "06_ROS_Redox_Glutathione",
  "07_Lipid_Membrane_Metabolism",
  "04_Photosynthesis_Energy_Metabolism",
  "08_Carbon_Starch_Sucrose_Osmoprotection",
  "05_Terpenoid_Secondary_Metabolism"
)

# Figure settings — aligned with the HSP and TF scripts.
FUNCTIONAL_Z_LIMIT       <- 2.5
FUNCTIONAL_CLUSTER_ROWS  <- TRUE
FUNCTIONAL_SPLIT_ROWS_BY_SUBFAMILY <- TRUE
FUNCTIONAL_ROW_MM        <- 3.7
FUNCTIONAL_FIG_WIDTH_MM  <- 225
FUNCTIONAL_RASTER_DPI    <- 600
FUNCTIONAL_WRITE_PNG     <- TRUE
FUNCTIONAL_WRITE_TIFF    <- TRUE
FUNCTIONAL_WRITE_PDF     <- TRUE
FUNCTIONAL_WRITE_SVG     <- TRUE

# The combined figure can be very tall because genes may occur in more than
# one broad KEGG-derived family. Individual-family figures are the main output.
FUNCTIONAL_WRITE_COMBINED <- TRUE
FUNCTIONAL_MAX_COMBINED_GENES <- 250L

set.seed(SEED)
options(stringsAsFactors = FALSE)

##############################
## 1) Packages
##############################
load_pkgs(
  c("ComplexHeatmap"),
  install_missing = INSTALL_MISSING_PKGS,
  bioc = TRUE
)
load_pkgs(
  c("dplyr", "tidyr", "tibble", "readr", "stringr",
    "circlize", "svglite"),
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
      object_name, " is missing required column(s): ",
      paste(missing_cols, collapse = ", ")
    )
  }
}

read_csv_flexible <- function(path) {
  if (!file.exists(path)) stop("Input file not found: ", path)
  x <- tryCatch(
    readr::read_csv(path, show_col_types = FALSE, progress = FALSE),
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
  as.data.frame(x, stringsAsFactors = FALSE)
}

zscore_rows_safe <- function(mat) {
  mat <- as.matrix(mat)
  storage.mode(mat) <- "numeric"
  mu <- rowMeans(mat, na.rm = TRUE)
  sds <- apply(mat, 1L, stats::sd, na.rm = TRUE)
  out <- sweep(mat, 1L, mu, "-")
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
  x <- iconv(x, to = "ASCII//TRANSLIT")
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x
}

parse_logical_override <- function(x) {
  y <- trimws(tolower(as.character(x)))
  out <- rep(NA, length(y))
  out[y %in% c("true", "t", "1", "yes", "y", "oui")] <- TRUE
  out[y %in% c("false", "f", "0", "no", "n", "non")] <- FALSE
  out
}

clean_text <- function(x) {
  x <- dplyr::coalesce(as.character(x), "")
  x <- gsub("\\s+", " ", trimws(x))
  x
}

first_informative_annotation <- function(annotation, fallback = "") {
  annotation <- clean_text(annotation)
  fallback <- clean_text(fallback)
  
  choose_one <- function(a, f) {
    if (!nzchar(a)) return(f)
    terms <- trimws(unlist(strsplit(a, ";", fixed = TRUE)))
    terms <- terms[nzchar(terms)]
    if (length(terms) == 0L) return(f)
    
    generic <- grepl(
      paste(
        "domain$|domain profile|signature$|family signature|",
        "superfamily$|motif$|atp-binding region|active-site signature|",
        "protein family$|^other$|^uncharacterized$|^hypothetical$",
        sep = ""
      ),
      terms,
      ignore.case = TRUE
    )
    
    preferred <- terms[!generic & nchar(terms) >= 4L & nchar(terms) <= 90L]
    if (length(preferred) == 0L) {
      preferred <- terms[nchar(terms) >= 4L]
    }
    if (length(preferred) == 0L) return(f)
    
    ans <- preferred[[1L]]
    ans <- gsub("[.;]+$", "", ans)
    if (nchar(ans) > 72L) ans <- paste0(substr(ans, 1L, 69L), "...")
    ans
  }
  
  vapply(
    seq_along(annotation),
    function(i) choose_one(annotation[[i]], fallback[[i]]),
    FUN.VALUE = character(1)
  )
}

normalise_subfamily <- function(x) {
  x <- clean_text(x)
  
  generic <- !nzchar(x) |
    stringr::str_to_lower(x) %in% c(
      "other",
      "others",
      "unknown",
      "unclassified",
      "unresolved",
      "na",
      "n/a"
    )
  
  x[generic] <- UNRESOLVED_LABEL
  x
}

family_prefix <- function(family_dir) family_dir

expected_family_files <- function(root, family_dir) {
  prefix <- family_prefix(family_dir)
  folder <- file.path(root, family_dir)
  
  list(
    folder = folder,
    expression = file.path(folder, paste0(prefix, "_expression.csv")),
    de = file.path(folder, paste0(prefix, "_DE_per_contrast.csv")),
    subfamily = file.path(folder, paste0(prefix, "_subfamily_assignment.csv")),
    provenance = file.path(folder, paste0(prefix, "_KEGG_provenance.csv"))
  )
}

validate_input_root <- function(root) {
  if (!dir.exists(root)) {
    stop("FUNCTIONAL_INPUT_DIR does not exist: ", root)
  }
  
  missing <- character(0)
  
  for (fam in FUNCTIONAL_FAMILY_DIRS) {
    ff <- expected_family_files(root, fam)
    
    for (nm in c("expression", "de", "subfamily", "provenance")) {
      if (!file.exists(ff[[nm]])) {
        missing <- c(missing, ff[[nm]])
      }
    }
  }
  
  if (length(missing) > 0L) {
    stop(
      "The functional input structure is incomplete. Missing file(s):\n",
      paste(missing, collapse = "\n")
    )
  }
  
  normalizePath(
    root,
    winslash = "/",
    mustWork = TRUE
  )
}


##############################
## 3) Literature-guided functional classification
##############################
## PRINCIPLE
## ------------------------------------------------------------
## classify_rule_one() uses explicit gene/enzyme names, enzyme families,
## pathway names and established pathway architecture.
##
## It does NOT infer a precise biochemical function from a generic
## "domain", "superfamily", "binding protein" or "stress protein" label.
##
## If no rule matches:
##   1. a non-generic supplied subfamily is retained as a fallback;
##   2. otherwise the gene becomes "Unresolved / manual curation".
##
## Unresolved genes remain in QC/curation outputs but are excluded from
## publication heatmaps by default.
##
## KEY PRIMARY LITERATURE USED TO BUILD THE CLASSIFICATION / ORDER
## ------------------------------------------------------------
## Proteostasis:
##   Xie et al. 2024 Mol Cell 84:3320-3335.e7
##   DOI: 10.1016/j.molcel.2024.07.033
##   Sanchez de Medina Hernandez et al. 2025 Dev Cell
##   DOI: 10.1016/j.devcel.2025.11.001
##   He et al. 2025 Plant Physiol, HsfA1a-BAG5b-autophagy
##   DOI: 10.1093/plphys/kiaf526
##
## Hormones:
##   SPL1/SPL12 -> PYL-mediated ABA thermotolerance:
##   PMID: 28400323
##   EIN3 -> ERF95/ERF97 thermotolerance:
##   PMID: 33793870
##   BES1 / BR heat responses:
##   DOI: 10.1016/j.plantsci.2020.110470
##   HSFA2 -> LOX3/OPR3 -> JA:
##   PMID: 38556243
##   MAX1/D14 strigolactone heat tolerance:
##   PMID: 37608548
##   Cytokinin pool and heat acclimation:
##   PMID: 32133021
##
## Ca2+ / MAPK / lipid signalling:
##   Finka et al. 2012 Plant Cell, CNGC thermal sensing
##   DOI: 10.1105/tpc.112.095844
##   MPK6 -> HsfA2:
##   PMID: 23638397
##   Heat -> PLD -> phosphatidic acid:
##   DOI: 10.1111/j.1365-313X.2009.03933.x
##   PLDalpha2 / FGT2 heat memory:
##   DOI: 10.1111/tpj.14927
##
## Redox:
##   Dard/Reichheld group, glutathione in heat adaptation:
##   DOI: 10.1093/jxb/erad042
##   CAT2 in long-term heat tolerance:
##   DOI: 10.1016/j.bbrc.2020.11.006
##
## Lipid remodelling:
##   Higashi et al. / Arabidopsis leaf lipid heat remodelling:
##   PMID: 32635518
##   TaMGD wheat heat tolerance:
##   PMID: 39837009
##   DGAT2/TAG accumulation under heat:
##   PMID: 38503190
##
## Carbon / osmoprotection:
##   TPS1/T6P/raffinose/carbon partitioning:
##   PMID: 37427798
##   T6P-KIN10-PIF4 thermoresponsive growth:
##   PMID: 31393060
##
## Specialised metabolism:
##   carrot phenylpropanoid protection against heat:
##   PMID: 27713760
##   recent heat multi-omics repeatedly recover phenylpropanoid/flavonoid
##   and carbon-remodelling responses (see bibliography table below).
##
## The categories below are therefore pathway modules, not claims that
## every carrot homolog has been functionally validated in heat stress.


classify_rule_one <- function(
    family_dir,
    annotation,
    pathways) {
  
  a <- stringr::str_to_lower(
    clean_text(annotation)
  )
  
  p <- stringr::str_to_lower(
    clean_text(pathways)
  )
  
  ap <- paste(a, p)
  
  # ----------------------------------------------------------
  # 01 PROTEOSTASIS / FOLDING / PROCESSING / DEGRADATION
  # ----------------------------------------------------------
  if (family_dir == "01_Protein_Folding_Processing_Degradation") {
    
    # Canonical chaperone network
    if (grepl(
      "small heat shock|hsp20|alpha[- ]crystallin|heat shock protein 1[5-9][.]?[0-9]* kda",
      ap
    )) return("Molecular chaperones — sHSP/HSP20")
    
    if (grepl(
      "\\bdnaj\\b|j[- ]domain|heat shock protein 40|hsp40",
      ap
    )) return("Molecular chaperones — HSP40/DnaJ")
    
    if (grepl(
      "chaperonin 60|cpn60|groel|hsp60",
      ap
    )) return("Molecular chaperones — HSP60/Cpn60")
    
    if (grepl(
      "chaperonin 10|cpn10|groes|hsp10",
      ap
    )) return("Molecular chaperones — HSP10/Cpn10")
    
    if (grepl(
      "hsp70|heat shock 70|heat shock cognate|hsc70|\\bbip\\b|hspa|dnak",
      ap
    )) return("Molecular chaperones — HSP70/HSC70")
    
    if (grepl(
      "hsp90|heat shock 90|endoplasmin|grp94|htpg",
      ap
    )) return("Molecular chaperones — HSP90")
    
    if (grepl(
      "hsp101|hsp100|clpb|caseinolytic peptidase b",
      ap
    )) return("Molecular chaperones — HSP100/ClpB")
    
    if (grepl(
      "chaperonin containing tcp|\\bcct[1-8]\\b|t[- ]complex protein 1|tcp[- ]1",
      ap
    )) return("CCT/TRiC chaperonin")
    
    if (grepl(
      "prefoldin|\\bpfd[1-6]\\b",
      ap
    )) return("Prefoldin")
    
    if (grepl(
      "\\bbag[1-9]\\b|bcl[- ]2[- ]associated athanogene|sti1|hop[123]?|hsp70[- ]hsp90 organizing|\\baha1\\b|p23 co[- ]chaperone|\\bsgt1\\b",
      ap
    )) return("HSP co-chaperones")
    
    # ER proteostasis / UPR
    if (grepl(
      "\\bire1[ab]?\\b|bzip60|unfolded protein response|er stress sensor",
      ap
    )) return("ER stress sensing / UPR")
    
    if (grepl(
      "calnexin|calreticulin|protein disulfide isomerase|\\bpdi\\b|ero1|peptidyl[- ]prolyl|cyclophilin|fk506[- ]binding|protein folding catalyst",
      ap
    )) return("ER folding / isomerases")
    
    if (grepl(
      "\\bhrd1\\b|sel1l|\\bos9\\b|derlin|\\bder1\\b|er[- ]associated degradation|\\berad\\b|\\bcdc48\\b|\\bufd1\\b|\\bnpl4\\b",
      ap
    )) return("ER-associated degradation (ERAD)")
    
    # Ubiquitin–proteasome system
    if (grepl(
      "ubiquitin[- ]activating enzyme|\\be1 enzyme\\b|uba[1-9]",
      ap
    )) return("Ubiquitin system — E1 activation")
    
    if (grepl(
      "ubiquitin[- ]conjugating|\\be2 enzyme\\b|\\bubc[0-9]+\\b",
      ap
    )) return("Ubiquitin system — E2 conjugation")
    
    if (grepl(
      "e3 ubiquitin|ubiquitin ligase|ring[- ]type|u[- ]box|\\bchip\\b|hect|f[- ]box",
      ap
    )) return("Ubiquitin system — E3 ligases")
    
    if (grepl(
      "deubiquitin|ubiquitin carboxyl[- ]terminal hydrolase|\\bubp[0-9]+\\b",
      ap
    )) return("Ubiquitin system — deubiquitination")
    
    if (grepl(
      "proteasome.*alpha|proteasome.*beta|20s proteasome|proteasome subunit alpha|proteasome subunit beta",
      ap
    )) return("Proteasome — 20S core")
    
    if (grepl(
      "26s protease regulatory|26s proteasome regulatory|proteasome regulatory subunit|\\brpt[1-6]\\b|\\brpn[1-9]",
      ap
    )) return("Proteasome — 19S regulatory particle")
    
    if (grepl(
      "proteasome",
      ap
    )) return("Proteasome — other defined component")
    
    # Autophagy / vacuolar proteolysis
    if (grepl(
      "\\batg1\\b|\\batg13\\b|autophagy initiation",
      ap
    )) return("Autophagy — initiation")
    
    if (grepl(
      "\\batg6\\b|beclin|vps34|vps15|phosphatidylinositol 3[- ]kinase.*autoph",
      ap
    )) return("Autophagy — PI3K/nucleation")
    
    if (grepl(
      "\\batg(3|4|5|7|8|10|12|16|18)[a-z0-9.-]*\\b|autophagosome|autophagy[- ]related protein",
      ap
    )) return("Autophagy — core conjugation/autophagosome")
    
    if (grepl(
      "vacuolar processing enzyme|vacuolar protease|cathepsin|vacuolar protein sorting|\\bvps[0-9]+\\b",
      ap
    )) return("Vacuolar degradation / trafficking")
    
    if (grepl(
      "subtilase|subtilisin|cysteine protease|cysteine proteinase|metacaspase|aspartic protease|serine carboxypeptidase|aminopeptidase|peptidase|protease",
      ap
    )) return("Proteases / peptidases")
  }
  
  
  # ----------------------------------------------------------
  # 02 HORMONE SIGNALLING + HOMEOSTASIS
  # ----------------------------------------------------------
  if (family_dir == "02_Hormone_Signaling_Metabolism") {
    
    # ABA
    if (grepl(
      "9[- ]cis[- ]epoxycarotenoid dioxygenase|\\bnced[0-9]*\\b|zeaxanthin epoxidase|\\baba1\\b|\\baba2\\b|short[- ]chain dehydrogenase.*aba|abscisic aldehyde oxidase|\\baao3\\b",
      ap
    )) return("ABA — biosynthesis")
    
    if (grepl(
      "cyp707a|aba 8['′-]?-hydroxylase|aba glucosyltransferase|\\bbg1\\b|\\bbg2\\b|aba beta[- ]glucosidase",
      ap
    )) return("ABA — catabolism / conjugation")
    
    if (grepl(
      "abcg25|abcg40|ait1|npf4[.]6|aba transporter",
      ap
    )) return("ABA — transport")
    
    if (grepl(
      "pyrabactin resistance|\\bpyr1\\b|\\bpyl[0-9]+\\b|\\brcar[0-9]+\\b",
      ap
    )) return("ABA — PYR/PYL perception")
    
    if (grepl(
      "protein phosphatase 2c|\\babi1\\b|\\babi2\\b|\\bhab1\\b|\\bhai[123]\\b|\\bpp2c",
      ap
    )) return("ABA — PP2C/SnRK2 core signalling")
    
    if (grepl(
      "\\bsnrk2[.]?[0-9]*\\b|open stomata 1|\\bost1\\b|aba[- ]activated protein kinase",
      ap
    )) return("ABA — PP2C/SnRK2 core signalling")
    
    if (grepl(
      "\\bareb[0-9]*\\b|\\babf[1-4]\\b|\\babi5\\b|aba[- ]responsive element[- ]binding",
      ap
    )) return("ABA — transcriptional output")
    
    # Auxin
    if (grepl(
      "tryptophan aminotransferase.*arabidopsis|\\btaa1\\b|\\btar[123]\\b|yucca|\\byuc[0-9]+\\b",
      ap
    )) return("Auxin — biosynthesis")
    
    if (grepl(
      "\\bgh3[.-]?[0-9]*\\b|dioxygenase for auxin oxidation|\\bdao[12]\\b|iaa[- ]amido synthetase|auxin conjugat",
      ap
    )) return("Auxin — homeostasis / conjugation")
    
    if (grepl(
      "pin[- ]formed|\\bpin[1-9]\\b|aux1|lax[1-3]|abcb[0-9]+|pgp[0-9]+|auxin transporter|auxin efflux|auxin influx",
      ap
    )) return("Auxin — transport")
    
    if (grepl(
      "\\btir1\\b|\\bafb[1-5]\\b|auxin receptor|aux/iaa|auxin[- ]responsive protein iaa|auxin response factor|\\barf[0-9]+\\b",
      ap
    )) return("Auxin — perception / signalling")
    
    # Ethylene
    if (grepl(
      "1[- ]aminocyclopropane[- ]1[- ]carboxylate synthase|\\bacs[0-9]+\\b|1[- ]aminocyclopropane[- ]1[- ]carboxylate oxidase|\\baco[0-9]+\\b",
      ap
    )) return("Ethylene — biosynthesis")
    
    if (grepl(
      "ethylene receptor|\\betr[12]\\b|\\bers[12]\\b|\\bein4\\b|\\bctr1\\b|\\bein2\\b|\\bein3\\b|ein3[- ]like|\\beil[1-5]\\b|\\bebf[12]\\b|ethylene response factor|ethylene[- ]responsive factor",
      ap
    )) return("Ethylene — perception / signalling")
    
    # Cytokinin
    if (grepl(
      "isopentenyltransferase|adenylate isopentenyltransferase|\\bipt[0-9]+\\b|lonely guy|\\blog[0-9]+\\b",
      ap
    )) return("Cytokinin — biosynthesis / activation")
    
    if (grepl(
      "cytokinin oxidase|cytokinin dehydrogenase|\\bckx[0-9]+\\b|cytokinin glucosyltransferase",
      ap
    )) return("Cytokinin — catabolism / conjugation")
    
    if (grepl(
      "abcg14|purine permease|\\bpup[0-9]+\\b|equilibrative nucleoside transporter|cytokinin transporter",
      ap
    )) return("Cytokinin — transport")
    
    if (grepl(
      "\\bahk[234]\\b|cytokinin receptor|histidine phosphotransfer protein|\\bahp[1-6]\\b|response regulator.*arr|\\barr[0-9]+\\b",
      ap
    )) return("Cytokinin — two-component signalling")
    
    # Gibberellin
    if (grepl(
      "ent[- ]copalyl diphosphate synthase|ent[- ]kaurene synthase|ent[- ]kaurene oxidase|kaurenoic acid oxidase|ga20[- ]oxidase|ga3[- ]oxidase|\\bga20ox|\\bga3ox",
      ap
    )) return("Gibberellin — biosynthesis")
    
    if (grepl(
      "ga2[- ]oxidase|\\bga2ox|gibberellin catabol",
      ap
    )) return("Gibberellin — catabolism")
    
    if (grepl(
      "gibberellin receptor|\\bgid1[abc]?\\b|\\bdella\\b|\\bgai\\b|\\brga\\b|\\brgl[123]\\b|\\bsly1\\b",
      ap
    )) return("Gibberellin — GID1/DELLA signalling")
    
    # Brassinosteroid
    if (grepl(
      "\\bdwf4\\b|\\bcpd\\b|cyp90|cyp85|det2|brassinosteroid biosynth",
      ap
    )) return("Brassinosteroid — biosynthesis")
    
    if (grepl(
      "\\bbri1\\b|\\bbak1\\b|serk3|\\bbsk[0-9]+\\b|\\bbsu1\\b|\\bbin2\\b|\\bbzr1\\b|\\bbes1\\b|brassinosteroid receptor|brassinosteroid signalling|brassinosteroid signaling",
      ap
    )) return("Brassinosteroid — signalling")
    
    # Jasmonate / oxylipin
    if (grepl(
      "13[- ]lipoxygenase|\\blo(x|x[0-9]+)\\b|allene oxide synthase|\\baos\\b|allene oxide cyclase|\\baoc[0-9]*\\b|12[- ]oxophytodienoate reductase|\\bopr3\\b|opcl1|jasmonic acid biosynth|alpha[- ]linolenic acid metabolism",
      ap
    )) return("Jasmonate — biosynthesis / oxylipin")
    
    if (grepl(
      "jasmonate catabol|cyp94|jasmonate sulfotransferase|\\bst2a\\b|jasmonate conjugat",
      ap
    )) return("Jasmonate — catabolism / conjugation")
    
    if (grepl(
      "\\bcoi1\\b|jasmonate receptor|jasmonate zim|\\bjaz[0-9]+\\b|\\bmyc2\\b|\\bmed25\\b",
      ap
    )) return("Jasmonate — COI1/JAZ signalling")
    
    # Salicylate
    if (grepl(
      "isochorismate synthase|\\bics1\\b|\\bsid2\\b|avrpphb susceptible 3|\\bpbs3\\b|enhanced pseudomonas susceptibility 1|\\beps1\\b|salicylic acid biosynth",
      ap
    )) return("Salicylate — biosynthesis")
    
    if (grepl(
      "salicylate hydroxylase|\\bs3h\\b|\\bs5h\\b|benzoic acid.*methyltransferase|\\bbsmt\\b|salicylic acid catabol|salicylate glucosyltransferase",
      ap
    )) return("Salicylate — catabolism / conjugation")
    
    if (grepl(
      "\\bnpr1\\b|\\bnpr3\\b|\\bnpr4\\b|nonexpressor of pathogenesis|\\btga[1-7]\\b|salicylic acid signaling|salicylic acid signalling",
      ap
    )) return("Salicylate — NPR/TGA signalling")
    
    # Strigolactone / karrikin
    if (grepl(
      "\\bd27\\b|\\bccd7\\b|\\bmax3\\b|\\bccd8\\b|\\bmax4\\b|\\bmax1\\b|strigolactone biosynth",
      ap
    )) return("Strigolactone — biosynthesis")
    
    if (grepl(
      "\\bd14\\b|\\bmax2\\b|\\bd3\\b|\\bsmxl[0-9]+\\b|\\bd53\\b|strigolactone receptor|strigolactone signaling|strigolactone signalling",
      ap
    )) return("Strigolactone — D14/MAX2 signalling")
    
    if (grepl(
      "\\bkai2\\b|karrikin|\\bsmax1\\b",
      ap
    )) return("Karrikin/KAI2 — signalling")
  }
  
  
  # ----------------------------------------------------------
  # 03 MAPK / CALCIUM / EARLY STRESS SIGNALLING
  # ----------------------------------------------------------
  if (family_dir == "03_MAPK_Calcium_Stress_Signaling") {
    
    # Receptors before generic kinase calls
    if (grepl(
      "leucine[- ]rich repeat receptor[- ]like kinase|lrr[- ]rlk|receptor[- ]like kinase|receptor kinase",
      ap
    )) return("Cell-surface receptor kinases")
    
    if (grepl(
      "receptor[- ]like cytoplasmic kinase|\\brlck\\b|\\bbik1\\b|\\bpbs1\\b|pbs1[- ]like",
      ap
    )) return("Receptor-proximal RLCKs")
    
    # Ca2+ influx / transport
    if (grepl(
      "cyclic nucleotide[- ]gated channel|\\bcngc[0-9]+\\b",
      ap
    )) return("Ca2+ influx — CNGC")
    
    if (grepl(
      "glutamate receptor[- ]like|\\bglr[0-9.]+\\b",
      ap
    )) return("Ca2+ influx — GLR")
    
    if (grepl(
      "\\bosca[0-9.]+\\b|hyperosmolality[- ]gated calcium|mechanosensitive channel|\\bmca[12]\\b|two[- ]pore channel|\\btpc1\\b|annexin",
      ap
    )) return("Ca2+ influx — other channels")
    
    if (grepl(
      "calcium[- ]transporting p[- ]type atpase|autoinhibited ca2\\+[- ]atpase|\\baca[0-9]+\\b|\\beca[0-9]+\\b|cation calcium exchanger|\\bcax[0-9]+\\b|calcium exchanger|calcium pump",
      ap
    )) return("Ca2+ efflux / sequestration")
    
    # Ca2+ decoding
    if (grepl(
      "calcium[- ]dependent protein kinase|\\bcdpk\\b|\\bcpk[0-9]+\\b",
      ap
    )) return("Ca2+ decoding — CPK/CDPK")
    
    if (grepl(
      "calcineurin b[- ]like|\\bcbl[0-9]+\\b|cbl[- ]interacting protein kinase|\\bcipk[0-9]+\\b",
      ap
    )) return("Ca2+ decoding — CBL/CIPK")
    
    if (grepl(
      "calmodulin[- ]like|calmodulin protein|\\bcml[0-9]+\\b|\\bcam[0-9]+\\b",
      ap
    )) return("Ca2+ decoding — CaM/CML")
    
    # ROS/lipid early signalling
    if (grepl(
      "respiratory burst oxidase homolog|\\brboh[a-j]\\b|nadph oxidase.*respiratory burst",
      ap
    )) return("ROS signalling — RBOH")
    
    if (grepl(
      "phospholipase d|\\bpld[a-z0-9]+\\b",
      ap
    )) return("Lipid signalling — PLD/phosphatidic acid")
    
    if (grepl(
      "phosphoinositide phospholipase c|phospholipase c|\\bplc[0-9]+\\b|non[- ]specific phospholipase c|\\bnpc[0-9]+\\b",
      ap
    )) return("Lipid signalling — PLC")
    
    if (grepl(
      "diacylglycerol kinase|\\bdgk[0-9]+\\b",
      ap
    )) return("Lipid signalling — DGK/phosphatidic acid")
    
    # MAPK tiers — order from upstream to downstream
    if (grepl(
      "map kinase kinase kinase|mapkkk|mekk|map3k",
      ap
    )) return("MAPK cascade — MAPKKK")
    
    if (grepl(
      "map kinase kinase|mapkk|\\bmkk[0-9]+\\b|map2k",
      ap
    )) return("MAPK cascade — MKK/MAPKK")
    
    if (grepl(
      "mitogen[- ]activated protein kinase|map kinase|\\bmpk[0-9]+\\b|\\bmapk[0-9]+\\b",
      ap
    )) return("MAPK cascade — MPK/MAPK")
    
    # Negative/adapter regulators
    if (grepl(
      "map kinase phosphatase|dual specificity phosphatase|protein phosphatase 2c|\\bpp2c\\b",
      ap
    )) return("Stress signalling phosphatases")
    
    if (grepl(
      "14[- ]3[- ]3 protein|general regulatory factor|\\bgrf[0-9]+\\b",
      ap
    )) return("14-3-3 signalling adaptors")
    
    if (grepl(
      "histidine kinase|histidine phosphotransfer|response regulator|two[- ]component",
      ap
    )) return("Two-component signalling")
  }
  
  
  # ----------------------------------------------------------
  # 04 PHOTOSYNTHESIS / ENERGY
  # ----------------------------------------------------------
  if (family_dir == "04_Photosynthesis_Energy_Metabolism") {
    
    # PSII antenna -> PSII -> OEC
    if (grepl(
      "light[- ]harvesting complex ii|\\blhcb[1-7]\\b|chlorophyll a[- ]b binding.*ii",
      ap
    )) return("PSII — LHCII antenna")
    
    if (grepl(
      "photosystem ii.*d1|\\bpsba\\b|photosystem ii.*d2|\\bpsbd\\b|\\bpsbb\\b|\\bpsbc\\b|photosystem ii reaction cent|photosystem ii core",
      ap
    )) return("PSII — reaction centre/core")
    
    if (grepl(
      "oxygen[- ]evolving|oxygen evolving|\\bpsbo\\b|\\bpsbp\\b|\\bpsbq\\b|\\bpsbr\\b|\\bpsbu\\b",
      ap
    )) return("PSII — oxygen-evolving complex")
    
    if (grepl(
      "\\bpsbs\\b|nonphotochemical quenching|npq[14]|violaxanthin de[- ]epoxidase|\\bvde\\b|xanthophyll cycle",
      ap
    )) return("Photoprotection / NPQ")
    
    # Inter-system chain
    if (grepl(
      "cytochrome b6|cytochrome b6f|\\bpeta\\b|\\bpetb\\b|\\bpetc\\b|\\bpetd\\b|cytochrome f",
      ap
    )) return("Cytochrome b6f complex")
    
    if (grepl(
      "plastocyanin|\\bpete\\b",
      ap
    )) return("Plastocyanin / intersystem transfer")
    
    # PSI
    if (grepl(
      "light[- ]harvesting complex i|\\blhca[1-6]\\b",
      ap
    )) return("PSI — LHCI antenna")
    
    if (grepl(
      "photosystem i|\\bpsa[a-l0-9]*\\b",
      ap
    )) return("PSI — reaction centre/core")
    
    if (grepl(
      "ferredoxin[- ]nadp|ferredoxin nadp reductase|\\bfnr[0-9]*\\b",
      ap
    )) return("Ferredoxin–FNR")
    
    if (grepl(
      "ferredoxin",
      ap
    )) return("Ferredoxin / electron carriers")
    
    if (grepl(
      "\\bpgr5\\b|\\bpgrl1\\b|ndh complex|nad\\(p\\)h dehydrogenase.*chloroplast|cyclic electron flow",
      ap
    )) return("Cyclic electron flow / chloroplast NDH")
    
    if (grepl(
      "chloroplast atp synthase|photosynthetic atp synthase|atp synthase.*chloroplast|cf0|cf1",
      ap
    )) return("Chloroplast ATP synthase")
    
    # Carbon fixation
    if (grepl(
      "ribulose[- ]1,5[- ]bisphosphate carboxylase|rubisco large|rubisco small|\\brbcl\\b|\\brbcs\\b",
      ap
    )) return("Calvin cycle — Rubisco")
    
    if (grepl(
      "rubisco activase|\\brca\\b",
      ap
    )) return("Calvin cycle — Rubisco activase")
    
    if (grepl(
      "phosphoribulokinase|sedoheptulose[- ]1,7[- ]bisphosphatase|fructose[- ]1,6[- ]bisphosphatase.*chloroplast|transketolase|ribose[- ]5[- ]phosphate isomerase|ribulose[- ]phosphate 3[- ]epimerase|glyceraldehyde[- ]3[- ]phosphate dehydrogenase.*chloroplast|calvin[- ]benson|carbon fixation in photosynthetic",
      ap
    )) return("Calvin cycle — regeneration/reduction")
    
    # Photorespiration
    if (grepl(
      "photorespirat|glycolate oxidase|glycine dehydrogenase|glycine decarboxylase|serine hydroxymethyltransferase|hydroxypyruvate reductase|glycerate kinase",
      ap
    )) return("Photorespiration")
    
    # Pigments
    if (grepl(
      "chlorophyll synthase|protochlorophyllide oxidoreductase|magnesium chelatase|chlorophyllide a oxygenase|chlorophyll biosynth|chlorophyll catabol|stay[- ]green",
      ap
    )) return("Chlorophyll metabolism")
    
    # Mitochondrial energy
    if (grepl(
      "complex i|nadh dehydrogenase.*mitochond|succinate dehydrogenase|complex ii|cytochrome bc1|complex iii|cytochrome c oxidase|complex iv|oxidative phosphorylation",
      ap
    )) return("Mitochondrial electron transport")
    
    if (grepl(
      "alternative oxidase|\\baox[0-9]*\\b|uncoupling protein",
      ap
    )) return("Alternative respiration")
    
    if (grepl(
      "atp synthase|f1f0 atp synthase",
      ap
    )) return("ATP synthase — other")
  }
  
  
  # ----------------------------------------------------------
  # 05 SPECIALISED METABOLISM
  # ----------------------------------------------------------
  if (family_dir == "05_Terpenoid_Secondary_Metabolism") {
    
    # Aromatic precursor / phenylpropanoid branches
    if (grepl(
      "shikimate pathway|3[- ]dehydroquinate|shikimate dehydrogenase|chorismate synthase",
      ap
    )) return("Shikimate / aromatic precursor pathway")
    
    if (grepl(
      "phenylalanine ammonia[- ]lyase|\\bpal[0-9]*\\b|cinnamate 4[- ]hydroxylase|\\bc4h\\b|4[- ]coumarate.*coa ligase|\\b4cl[0-9]*\\b",
      ap
    )) return("Phenylpropanoid — gateway PAL/C4H/4CL")
    
    if (grepl(
      "hydroxycinnamoyl.*transferase|\\bhct\\b|coumarate 3[- ]hydroxylase|\\bc3h\\b|caffeoyl shikimate esterase|\\bcse\\b",
      ap
    )) return("Phenylpropanoid — hydroxycinnamate branch")
    
    if (grepl(
      "cinnamoyl[- ]coa reductase|\\bccr[0-9]*\\b|cinnamyl alcohol dehydrogenase|\\bcad[0-9]*\\b|caffeoyl[- ]coa o[- ]methyltransferase|ccoaomt|ferulate 5[- ]hydroxylase|\\bf5h\\b|laccase|monolignol|lignin",
      ap
    )) return("Phenylpropanoid — monolignol/lignin")
    
    if (grepl(
      "chalcone synthase|\\bchs\\b|chalcone isomerase|\\bchi\\b|flavanone 3[- ]hydroxylase|\\bf3h\\b|flavonoid biosynth",
      ap
    )) return("Flavonoid — core pathway")
    
    if (grepl(
      "flavonol synthase|\\bfls[0-9]*\\b|flavonol",
      ap
    )) return("Flavonoid — flavonol")
    
    if (grepl(
      "dihydroflavonol 4[- ]reductase|\\bdfr\\b|anthocyanidin synthase|\\bans\\b|udp[- ]glucose.*flavonoid 3[- ]o[- ]glucosyltransferase|ufgt|anthocyan",
      ap
    )) return("Flavonoid — anthocyanin")
    
    if (grepl(
      "leucoanthocyanidin reductase|\\blar\\b|anthocyanidin reductase|\\banr\\b|proanthocyanidin|condensed tannin",
      ap
    )) return("Flavonoid — proanthocyanidin")
    
    if (grepl(
      "isoflavone synthase|isoflavonoid",
      ap
    )) return("Isoflavonoid")
    
    if (grepl(
      "coumarin|scopoletin|feruloyl[- ]coa 6['′-]?-hydroxylase|f6['′]?h|scopoletin 8[- ]hydroxylase",
      ap
    )) return("Coumarin / scopoletin")
    
    # Carrot-specific defence metabolites when explicitly annotated
    if (grepl(
      "6[- ]methoxymellein|methoxymellein|mellein|polyketide synthase.*mellein",
      ap
    )) return("Carrot defence metabolites — 6-methoxymellein/polyketide")
    
    if (grepl(
      "polyacetylene|falcarinol|falcarindiol|acetylenase",
      ap
    )) return("Carrot defence metabolites — polyacetylenes")
    
    # Isoprenoid backbone
    if (grepl(
      "acetyl[- ]coa acetyltransferase|acetoacetyl[- ]coa thiolase|hydroxymethylglutaryl[- ]coa synthase|\\bhmgs\\b|hydroxymethylglutaryl[- ]coa reductase|\\bhmgr\\b|mevalonate kinase|phosphomevalonate kinase|mevalonate diphosphate decarboxylase",
      ap
    )) return("Isoprenoid backbone — mevalonate (MVA)")
    
    if (grepl(
      "1[- ]deoxy[- ]d[- ]xylulose 5[- ]phosphate synthase|\\bdxs\\b|1[- ]deoxy[- ]d[- ]xylulose 5[- ]phosphate reductoisomerase|\\bdxr\\b|methylerythritol|mep pathway|\\bmct\\b|\\bcmk\\b|\\bmds\\b|\\bhds\\b|\\bhdr\\b",
      ap
    )) return("Isoprenoid backbone — MEP")
    
    if (grepl(
      "geranyl diphosphate synthase|\\bgpps\\b",
      ap
    )) return("Prenyl diphosphate synthesis — GPP")
    
    if (grepl(
      "farnesyl diphosphate synthase|\\bfpps\\b",
      ap
    )) return("Prenyl diphosphate synthesis — FPP")
    
    if (grepl(
      "geranylgeranyl diphosphate synthase|\\bggpps\\b",
      ap
    )) return("Prenyl diphosphate synthesis — GGPP")
    
    if (grepl(
      "monoterpene synthase|monoterpen|limonene synthase|pinene synthase|geraniol synthase",
      ap
    )) return("Terpenoid — monoterpene")
    
    if (grepl(
      "sesquiterpene synthase|sesquiterpen|farnesene synthase",
      ap
    )) return("Terpenoid — sesquiterpene")
    
    if (grepl(
      "diterpene synthase|diterpen|ent[- ]kaurene synthase",
      ap
    )) return("Terpenoid — diterpene")
    
    if (grepl(
      "squalene synthase|squalene epoxidase|oxidosqualene cyclase|triterpen|sterol biosynth",
      ap
    )) return("Terpenoid — triterpene/sterol")
    
    # Carotenoid and prenylquinone
    if (grepl(
      "phytoene synthase|\\bpsy\\b|phytoene desaturase|\\bpds\\b|zeta[- ]carotene desaturase|\\bzds\\b|carotene isomerase|crtiso|lycopene beta[- ]cyclase|lcyb|lycopene epsilon[- ]cyclase|lcye|beta[- ]carotene hydroxylase|carotenoid biosynth",
      ap
    )) return("Carotenoid biosynthesis")
    
    if (grepl(
      "tocopherol|homogentisate phytyltransferase|\\bvte[1-6]\\b|plastoquinone|prenylquinone",
      ap
    )) return("Tocopherol / prenylquinone")
    
    if (grepl(
      "alkaloid biosynth|alkaloid",
      ap
    )) return("Alkaloid metabolism")
  }
  
  
  # ----------------------------------------------------------
  # 06 ROS / REDOX / GLUTATHIONE
  # ----------------------------------------------------------
  if (family_dir == "06_ROS_Redox_Glutathione") {
    
    if (grepl(
      "respiratory burst oxidase homolog|\\brboh[a-j]\\b|nadph oxidase",
      ap
    )) return("ROS production — RBOH/NADPH oxidase")
    
    if (grepl(
      "superoxide dismutase|\\bsod\\b|copper zinc superoxide|manganese superoxide|iron superoxide",
      ap
    )) return("ROS scavenging — SOD")
    
    if (grepl(
      "catalase|\\bcat[123]\\b",
      ap
    )) return("ROS scavenging — catalase")
    
    if (grepl(
      "ascorbate peroxidase|\\bapx[0-9]*\\b",
      ap
    )) return("Ascorbate–glutathione cycle — APX")
    
    if (grepl(
      "monodehydroascorbate reductase|\\bmdhar\\b",
      ap
    )) return("Ascorbate–glutathione cycle — MDHAR")
    
    if (grepl(
      "dehydroascorbate reductase|\\bdhar\\b",
      ap
    )) return("Ascorbate–glutathione cycle — DHAR")
    
    if (grepl(
      "glutathione reductase|\\bgr[12]\\b",
      ap
    )) return("Ascorbate–glutathione cycle — GR")
    
    if (grepl(
      "l[- ]galactose pathway|vtc[1245]|gdp[- ]mannose.*ascorbate|ascorbate biosynth",
      ap
    )) return("Ascorbate biosynthesis")
    
    if (grepl(
      "gamma[- ]glutamylcysteine synthetase|glutamate[- ]cysteine ligase|\\bgsh1\\b|glutathione synthetase|\\bgsh2\\b",
      ap
    )) return("Glutathione biosynthesis")
    
    if (grepl(
      "glutathione s[- ]transferase|glutathione transferase|\\bgst[a-z0-9]+\\b",
      ap
    )) return("Glutathione transferases (GST)")
    
    if (grepl(
      "glutathione peroxidase|\\bgpx[0-9]*\\b|gpx[- ]like",
      ap
    )) return("Glutathione peroxidase / GPXL")
    
    if (grepl(
      "thioredoxin reductase|nadph[- ]thioredoxin reductase|\\bntra\\b|\\bntrb\\b|\\bntrc\\b",
      ap
    )) return("Thioredoxin system — reductases")
    
    if (grepl(
      "thioredoxin|\\btrx[- ]?[a-z0-9]+\\b",
      ap
    )) return("Thioredoxin system — TRX")
    
    if (grepl(
      "glutaredoxin|\\bgrx[a-z0-9-]*\\b",
      ap
    )) return("Glutaredoxin system")
    
    if (grepl(
      "peroxiredoxin|peroxiredoxin|2[- ]cys peroxiredoxin|\\bprx[a-z0-9-]*\\b",
      ap
    )) return("Peroxiredoxin")
    
    if (grepl(
      "class iii peroxidase|secretory peroxidase|guaiacol peroxidase|peroxidase [0-9]+",
      ap
    )) return("Class III peroxidases")
    
    if (grepl(
      "methionine sulfoxide reductase|\\bmsr[ab]",
      ap
    )) return("Protein redox repair — MSR")
    
    if (grepl(
      "glyoxalase i|lactoylglutathione lyase|glyoxalase ii|hydroxyacylglutathione hydrolase|methylglyoxal",
      ap
    )) return("Methylglyoxal detoxification / glyoxalase")
    
    if (grepl(
      "glucose[- ]6[- ]phosphate dehydrogenase|6[- ]phosphogluconate dehydrogenase|nadp[- ]dependent isocitrate dehydrogenase|nadp[- ]malic enzyme",
      ap
    )) return("NADPH supply for redox homeostasis")
    
    if (grepl(
      "alternative oxidase|\\baox[0-9]*\\b",
      ap
    )) return("Mitochondrial redox — AOX")
    
    if (grepl(
      "peroxisom|glycolate oxidase",
      ap
    )) return("Peroxisomal redox")
    
    if (grepl(
      "nitric oxide|gsnor|s[- ]nitrosoglutathione reductase",
      ap
    )) return("NO / S-nitrosoglutathione redox")
  }
  
  
  # ----------------------------------------------------------
  # 07 LIPID / MEMBRANE METABOLISM
  # ----------------------------------------------------------
  if (family_dir == "07_Lipid_Membrane_Metabolism") {
    
    if (grepl(
      "acetyl[- ]coa carboxylase|ketoacyl[- ]acp synthase|\\bkas[123]\\b|acyl carrier protein|malonyl[- ]coa.*acp|fatty acid synth",
      ap
    )) return("Fatty acid synthesis")
    
    if (grepl(
      "fatty acid desaturase|\\bfad[2-8]\\b|stearoyl[- ]acp desaturase|desaturase",
      ap
    )) return("Fatty acid desaturation")
    
    if (grepl(
      "3[- ]ketoacyl[- ]coa synthase|\\bkcs[0-9]+\\b|very[- ]long[- ]chain fatty|fatty acid elongase|elongation of very long",
      ap
    )) return("Very-long-chain FA elongation")
    
    if (grepl(
      "long[- ]chain acyl[- ]coa synthetase|\\blacs[0-9]+\\b|acyl[- ]coa synthetase|acyl[- ]acp thioesterase",
      ap
    )) return("Acyl activation / acyl editing")
    
    if (grepl(
      "monogalactosyldiacylglycerol synthase|\\bmgd[123]\\b",
      ap
    )) return("Thylakoid lipids — MGDG")
    
    if (grepl(
      "digalactosyldiacylglycerol synthase|\\bdgd[12]\\b",
      ap
    )) return("Thylakoid lipids — DGDG")
    
    if (grepl(
      "sulfoquinovosyldiacylglycerol|\\bsqd[12]\\b|sulfolipid",
      ap
    )) return("Thylakoid lipids — SQDG")
    
    if (grepl(
      "glycerol[- ]3[- ]phosphate acyltransferase|\\bgpat[0-9]+\\b|lysophosphatidic acid acyltransferase|\\blpaat\\b",
      ap
    )) return("Glycerolipid backbone synthesis")
    
    if (grepl(
      "diacylglycerol acyltransferase|\\bdgat[123]\\b|phospholipid:diacylglycerol acyltransferase|\\bpdat\\b|triacylglycerol",
      ap
    )) return("Triacylglycerol synthesis / storage")
    
    if (grepl(
      "phospholipase d|\\bpld[a-z0-9]+\\b",
      ap
    )) return("Phospholipid signalling — PLD/PA")
    
    if (grepl(
      "phospholipase c|\\bplc[0-9]+\\b|non[- ]specific phospholipase c|\\bnpc[0-9]+\\b",
      ap
    )) return("Phospholipid signalling — PLC")
    
    if (grepl(
      "phospholipase a|\\bpla[12]\\b|patatin[- ]like phospholipase",
      ap
    )) return("Phospholipid remodelling — PLA/lipases")
    
    if (grepl(
      "diacylglycerol kinase|\\bdgk[0-9]+\\b|phosphatidic acid phosphatase|lipin",
      ap
    )) return("Phosphatidic acid / DAG interconversion")
    
    if (grepl(
      "phosphatidylcholine|phosphatidylethanolamine|phosphatidylserine|phosphatidylinositol|cdp[- ]diacylglycerol|phospholipid biosynth",
      ap
    )) return("Structural phospholipid metabolism")
    
    if (grepl(
      "sphingolipid|ceramide synthase|long[- ]chain base|sphingosine|glycosylceramide|gipc",
      ap
    )) return("Sphingolipid / ceramide metabolism")
    
    if (grepl(
      "sterol biosynth|cycloartenol|sitosterol|campesterol|sterol methyltransferase|sterol glucosyltransferase",
      ap
    )) return("Sterol metabolism")
    
    if (grepl(
      "cutin|suberin|wax biosynth|eceriferum|\\bcer[0-9]+\\b|fatty acyl[- ]coa reductase|wax synthase",
      ap
    )) return("Cuticle / wax / suberin")
    
    if (grepl(
      "lipoxygenase|\\blox[0-9]+\\b|allene oxide synthase|allene oxide cyclase|oxylipin|alpha[- ]linolenic acid metabolism",
      ap
    )) return("Oxylipin / JA-precursor lipid metabolism")
    
    if (grepl(
      "acyl[- ]coa oxidase|\\bacx[0-9]+\\b|multifunctional protein.*beta[- ]oxidation|3[- ]ketoacyl[- ]coa thiolase|fatty acid degradation|beta[- ]oxidation|β[- ]oxidation",
      ap
    )) return("Fatty acid beta-oxidation")
    
    if (grepl(
      "lipid transfer protein|nsltp|lipid transfer",
      ap
    )) return("Lipid transfer / trafficking")
  }
  
  
  # ----------------------------------------------------------
  # 08 CARBON / SUGAR / OSMOPROTECTION
  # ----------------------------------------------------------
  if (family_dir == "08_Carbon_Starch_Sucrose_Osmoprotection") {
    
    if (grepl(
      "hexokinase|\\bhxk[123]\\b|fructokinase|\\bfrk[1-7]\\b",
      ap
    )) return("Hexose phosphorylation / sugar sensing")
    
    if (grepl(
      "phosphofructokinase|fructose[- ]bisphosphate aldolase|triosephosphate isomerase|phosphoglycerate kinase|enolase|pyruvate kinase|glycolysis",
      ap
    )) return("Glycolysis")
    
    if (grepl(
      "pyruvate dehydrogenase|pyruvate decarboxylase|acetyl[- ]coa.*pyruvate",
      ap
    )) return("Pyruvate metabolism")
    
    if (grepl(
      "citrate synthase|aconitase|isocitrate dehydrogenase|2[- ]oxoglutarate dehydrogenase|succinate dehydrogenase|fumarase|malate dehydrogenase|tricarboxylic acid|tca cycle|citrate cycle",
      ap
    )) return("TCA cycle")
    
    if (grepl(
      "isocitrate lyase|malate synthase|glyoxylate cycle",
      ap
    )) return("Glyoxylate cycle")
    
    if (grepl(
      "glucose[- ]6[- ]phosphate dehydrogenase|6[- ]phosphogluconate dehydrogenase|transketolase|transaldolase|pentose phosphate",
      ap
    )) return("Pentose phosphate pathway")
    
    # Starch
    if (grepl(
      "adp[- ]glucose pyrophosphorylase|\\bagpase\\b|granule[- ]bound starch synthase|\\bgbss\\b|starch synthase|starch branching enzyme|\\bsbe[123]\\b",
      ap
    )) return("Starch — synthesis")
    
    if (grepl(
      "beta[- ]amylase|\\bbam[1-9]\\b|alpha[- ]amylase|\\bamy[1-3]\\b|starch debranching|isoamylase|\\bisa[1-3]\\b|glucan water dikinase|\\bgwd\\b|phosphoglucan water dikinase|\\bpwd\\b",
      ap
    )) return("Starch — degradation/remobilisation")
    
    # Sucrose
    if (grepl(
      "sucrose[- ]phosphate synthase|\\bsps\\b|sucrose[- ]phosphate phosphatase|\\bspp\\b",
      ap
    )) return("Sucrose — synthesis")
    
    if (grepl(
      "sucrose synthase|\\bsus[1-6]\\b|invertase|beta[- ]fructofuranosidase",
      ap
    )) return("Sucrose — cleavage/remobilisation")
    
    if (grepl(
      "sucrose transporter|sucrose carrier|\\bsut[1-9]\\b|\\bsuc[1-9]\\b|sugar will eventually be exported transporter|\\bsweet[0-9]+\\b|monosaccharide transporter",
      ap
    )) return("Sugar transport")
    
    # T6P / trehalose
    if (grepl(
      "trehalose[- ]6[- ]phosphate synthase|\\btps[0-9]+\\b",
      ap
    )) return("Trehalose/T6P — TPS")
    
    if (grepl(
      "trehalose[- ]6[- ]phosphate phosphatase|\\btpp[a-z0-9]+\\b",
      ap
    )) return("Trehalose/T6P — TPP")
    
    if (grepl(
      "trehalase|trehalose metabolism",
      ap
    )) return("Trehalose — turnover")
    
    # Raffinose family
    if (grepl(
      "galactinol synthase|\\bgols[0-9]+\\b",
      ap
    )) return("Raffinose-family oligosaccharides — galactinol")
    
    if (grepl(
      "raffinose synthase|stachyose synthase|raffinose family|raffinose metabolism",
      ap
    )) return("Raffinose-family oligosaccharides — synthesis")
    
    # Compatible solutes
    if (grepl(
      "delta[- ]1[- ]pyrroline[- ]5[- ]carboxylate synthase|pyrroline[- ]5[- ]carboxylate synthetase|\\bp5cs[12]?\\b|pyrroline[- ]5[- ]carboxylate reductase|\\bp5cr\\b",
      ap
    )) return("Compatible solutes — proline synthesis")
    
    if (grepl(
      "proline dehydrogenase|\\bprodh[12]?\\b|pyrroline[- ]5[- ]carboxylate dehydrogenase|proline catabol",
      ap
    )) return("Compatible solutes — proline catabolism")
    
    if (grepl(
      "choline monooxygenase|betaine aldehyde dehydrogenase|glycine betaine|betaine biosynth",
      ap
    )) return("Compatible solutes — glycine betaine")
    
    if (grepl(
      "myo[- ]inositol|inositol oxygenase|sorbitol|mannitol|polyol",
      ap
    )) return("Compatible solutes — polyols/inositol")
    
    if (grepl(
      "glutamate decarboxylase|\\bgad[0-9]+\\b|gaba transaminase|gamma[- ]aminobutyric|gaba shunt",
      ap
    )) return("GABA shunt")
    
    # Cell wall carbohydrate allocation
    if (grepl(
      "cellulose synthase|\\bcesa[0-9]+\\b|cellulose synthase[- ]like|\\bcsl[a-h]",
      ap
    )) return("Cell wall carbon — cellulose/hemicellulose")
    
    if (grepl(
      "xyloglucan endotransglucosylase|xyloglucan endotransglycosylase|\\bxth[0-9]+\\b|pectin methylesterase|polygalacturonase|pectate lyase|pectin",
      ap
    )) return("Cell wall carbon — matrix polysaccharides")
  }
  
  UNRESOLVED_LABEL
}


classify_rule_vector <- function(
    family_dir,
    annotation,
    pathways) {
  
  mapply(
    FUN = function(ann, pth) {
      classify_rule_one(
        family_dir,
        ann,
        pth
      )
    },
    annotation,
    pathways,
    USE.NAMES = FALSE
  )
}


##############################
## 3B) Biologically meaningful subfamily order
##############################
subfamily_order_for_family <- function(family_dir) {
  
  orders <- list(
    
    "01_Protein_Folding_Processing_Degradation" = c(
      "Molecular chaperones — sHSP/HSP20",
      "Molecular chaperones — HSP40/DnaJ",
      "Molecular chaperones — HSP60/Cpn60",
      "Molecular chaperones — HSP10/Cpn10",
      "Molecular chaperones — HSP70/HSC70",
      "Molecular chaperones — HSP90",
      "Molecular chaperones — HSP100/ClpB",
      "CCT/TRiC chaperonin",
      "Prefoldin",
      "HSP co-chaperones",
      "ER stress sensing / UPR",
      "ER folding / isomerases",
      "ER-associated degradation (ERAD)",
      "Ubiquitin system — E1 activation",
      "Ubiquitin system — E2 conjugation",
      "Ubiquitin system — E3 ligases",
      "Ubiquitin system — deubiquitination",
      "Proteasome — 20S core",
      "Proteasome — 19S regulatory particle",
      "Proteasome — other defined component",
      "Autophagy — initiation",
      "Autophagy — PI3K/nucleation",
      "Autophagy — core conjugation/autophagosome",
      "Vacuolar degradation / trafficking",
      "Proteases / peptidases"
    ),
    
    "02_Hormone_Signaling_Metabolism" = c(
      "ABA — biosynthesis",
      "ABA — catabolism / conjugation",
      "ABA — transport",
      "ABA — PYR/PYL perception",
      "ABA — PP2C/SnRK2 core signalling",
      "ABA — transcriptional output",
      "Auxin — biosynthesis",
      "Auxin — homeostasis / conjugation",
      "Auxin — transport",
      "Auxin — perception / signalling",
      "Ethylene — biosynthesis",
      "Ethylene — perception / signalling",
      "Cytokinin — biosynthesis / activation",
      "Cytokinin — catabolism / conjugation",
      "Cytokinin — transport",
      "Cytokinin — two-component signalling",
      "Gibberellin — biosynthesis",
      "Gibberellin — catabolism",
      "Gibberellin — GID1/DELLA signalling",
      "Brassinosteroid — biosynthesis",
      "Brassinosteroid — signalling",
      "Jasmonate — biosynthesis / oxylipin",
      "Jasmonate — catabolism / conjugation",
      "Jasmonate — COI1/JAZ signalling",
      "Salicylate — biosynthesis",
      "Salicylate — catabolism / conjugation",
      "Salicylate — NPR/TGA signalling",
      "Strigolactone — biosynthesis",
      "Strigolactone — D14/MAX2 signalling",
      "Karrikin/KAI2 — signalling"
    ),
    
    "03_MAPK_Calcium_Stress_Signaling" = c(
      "Cell-surface receptor kinases",
      "Receptor-proximal RLCKs",
      "Ca2+ influx — CNGC",
      "Ca2+ influx — GLR",
      "Ca2+ influx — other channels",
      "Ca2+ efflux / sequestration",
      "Ca2+ decoding — CaM/CML",
      "Ca2+ decoding — CPK/CDPK",
      "Ca2+ decoding — CBL/CIPK",
      "ROS signalling — RBOH",
      "Lipid signalling — PLD/phosphatidic acid",
      "Lipid signalling — PLC",
      "Lipid signalling — DGK/phosphatidic acid",
      "MAPK cascade — MAPKKK",
      "MAPK cascade — MKK/MAPKK",
      "MAPK cascade — MPK/MAPK",
      "Stress signalling phosphatases",
      "14-3-3 signalling adaptors",
      "Two-component signalling"
    ),
    
    "04_Photosynthesis_Energy_Metabolism" = c(
      "PSII — LHCII antenna",
      "PSII — reaction centre/core",
      "PSII — oxygen-evolving complex",
      "Photoprotection / NPQ",
      "Cytochrome b6f complex",
      "Plastocyanin / intersystem transfer",
      "PSI — LHCI antenna",
      "PSI — reaction centre/core",
      "Ferredoxin / electron carriers",
      "Ferredoxin–FNR",
      "Cyclic electron flow / chloroplast NDH",
      "Chloroplast ATP synthase",
      "Calvin cycle — Rubisco",
      "Calvin cycle — Rubisco activase",
      "Calvin cycle — regeneration/reduction",
      "Photorespiration",
      "Chlorophyll metabolism",
      "Mitochondrial electron transport",
      "Alternative respiration",
      "ATP synthase — other"
    ),
    
    "05_Terpenoid_Secondary_Metabolism" = c(
      "Shikimate / aromatic precursor pathway",
      "Phenylpropanoid — gateway PAL/C4H/4CL",
      "Phenylpropanoid — hydroxycinnamate branch",
      "Phenylpropanoid — monolignol/lignin",
      "Flavonoid — core pathway",
      "Flavonoid — flavonol",
      "Flavonoid — anthocyanin",
      "Flavonoid — proanthocyanidin",
      "Isoflavonoid",
      "Coumarin / scopoletin",
      "Carrot defence metabolites — 6-methoxymellein/polyketide",
      "Carrot defence metabolites — polyacetylenes",
      "Isoprenoid backbone — mevalonate (MVA)",
      "Isoprenoid backbone — MEP",
      "Prenyl diphosphate synthesis — GPP",
      "Prenyl diphosphate synthesis — FPP",
      "Prenyl diphosphate synthesis — GGPP",
      "Terpenoid — monoterpene",
      "Terpenoid — sesquiterpene",
      "Terpenoid — diterpene",
      "Terpenoid — triterpene/sterol",
      "Carotenoid biosynthesis",
      "Tocopherol / prenylquinone",
      "Alkaloid metabolism"
    ),
    
    "06_ROS_Redox_Glutathione" = c(
      "ROS production — RBOH/NADPH oxidase",
      "ROS scavenging — SOD",
      "ROS scavenging — catalase",
      "Ascorbate biosynthesis",
      "Ascorbate–glutathione cycle — APX",
      "Ascorbate–glutathione cycle — MDHAR",
      "Ascorbate–glutathione cycle — DHAR",
      "Ascorbate–glutathione cycle — GR",
      "Glutathione biosynthesis",
      "Glutathione transferases (GST)",
      "Glutathione peroxidase / GPXL",
      "Thioredoxin system — reductases",
      "Thioredoxin system — TRX",
      "Glutaredoxin system",
      "Peroxiredoxin",
      "Class III peroxidases",
      "Protein redox repair — MSR",
      "Methylglyoxal detoxification / glyoxalase",
      "NADPH supply for redox homeostasis",
      "Mitochondrial redox — AOX",
      "Peroxisomal redox",
      "NO / S-nitrosoglutathione redox"
    ),
    
    "07_Lipid_Membrane_Metabolism" = c(
      "Fatty acid synthesis",
      "Fatty acid desaturation",
      "Very-long-chain FA elongation",
      "Acyl activation / acyl editing",
      "Glycerolipid backbone synthesis",
      "Thylakoid lipids — MGDG",
      "Thylakoid lipids — DGDG",
      "Thylakoid lipids — SQDG",
      "Structural phospholipid metabolism",
      "Phospholipid remodelling — PLA/lipases",
      "Phospholipid signalling — PLC",
      "Phospholipid signalling — PLD/PA",
      "Phosphatidic acid / DAG interconversion",
      "Triacylglycerol synthesis / storage",
      "Sphingolipid / ceramide metabolism",
      "Sterol metabolism",
      "Cuticle / wax / suberin",
      "Oxylipin / JA-precursor lipid metabolism",
      "Fatty acid beta-oxidation",
      "Lipid transfer / trafficking"
    ),
    
    "08_Carbon_Starch_Sucrose_Osmoprotection" = c(
      "Hexose phosphorylation / sugar sensing",
      "Glycolysis",
      "Pyruvate metabolism",
      "TCA cycle",
      "Glyoxylate cycle",
      "Pentose phosphate pathway",
      "Starch — synthesis",
      "Starch — degradation/remobilisation",
      "Sucrose — synthesis",
      "Sucrose — cleavage/remobilisation",
      "Sugar transport",
      "Trehalose/T6P — TPS",
      "Trehalose/T6P — TPP",
      "Trehalose — turnover",
      "Raffinose-family oligosaccharides — galactinol",
      "Raffinose-family oligosaccharides — synthesis",
      "Compatible solutes — proline synthesis",
      "Compatible solutes — proline catabolism",
      "Compatible solutes — glycine betaine",
      "Compatible solutes — polyols/inositol",
      "GABA shunt",
      "Cell wall carbon — cellulose/hemicellulose",
      "Cell wall carbon — matrix polysaccharides"
    )
  )
  
  orders[[family_dir]] %||% character(0)
}


order_subfamilies_for_plot <- function(
    family_dir,
    observed) {
  
  observed <- unique(
    as.character(observed)
  )
  
  desired <- subfamily_order_for_family(
    family_dir
  )
  
  # Literature-defined categories first in biological pathway order.
  # Any valid supplied fallback categories are placed afterwards,
  # alphabetically. Unresolved is never added to plotted levels.
  c(
    desired[desired %in% observed],
    sort(
      setdiff(
        observed,
        c(
          desired,
          UNRESOLVED_LABEL
        )
      )
    )
  )
}


##############################
## 4) Manual overrides
##############################
apply_manual_overrides <- function(tbl, override_file) {
  template <- tbl |>
    dplyr::transmute(
      Family_ID,
      Family,
      Gene,
      gene_ID,
      Subfamily_source,
      Subfamily_auto,
      Short_name_auto,
      Include_auto,
      Subfamily_override = "",
      Short_name_override = "",
      Include_override = "",
      Curator_note = ""
    )
  
  if (!file.exists(override_file)) {
    readr::write_csv(template, override_file)
    message("Created manual-override template: ", override_file)
  }
  
  ov <- read_csv_flexible(override_file)
  assert_columns(ov, c("Family_ID", "Gene"), "Functional override file")
  
  for (nm in c(
    "Subfamily_override", "Short_name_override",
    "Include_override", "Curator_note"
  )) {
    if (!nm %in% colnames(ov)) ov[[nm]] <- ""
  }
  
  ov <- ov |>
    dplyr::mutate(
      dplyr::across(
        c(
          Subfamily_override,
          Short_name_override,
          Include_override,
          Curator_note
        ),
        ~ trimws(as.character(.x))
      ),
      Include_override_parsed = parse_logical_override(Include_override)
    ) |>
    dplyr::select(
      Family_ID,
      Gene,
      Subfamily_override,
      Short_name_override,
      Include_override_parsed,
      Curator_note
    )
  
  tbl |>
    dplyr::left_join(ov, by = c("Family_ID", "Gene")) |>
    dplyr::mutate(
      Subfamily = dplyr::if_else(
        !is.na(Subfamily_override) & nzchar(Subfamily_override),
        Subfamily_override,
        Subfamily_auto
      ),
      Short_name = dplyr::if_else(
        !is.na(Short_name_override) & nzchar(Short_name_override),
        Short_name_override,
        Short_name_auto
      ),
      Include = dplyr::if_else(
        !is.na(Include_override_parsed),
        Include_override_parsed,
        Include_auto
      ),
      Curator_note = dplyr::coalesce(Curator_note, ""),
      Label = dplyr::if_else(
        !is.na(Short_name) & nzchar(Short_name),
        paste0(Gene, " | ", Short_name),
        Gene
      )
    )
}

##############################
## 5) Figure functions
##############################
draw_heatmap_object <- function(ht) {
  ComplexHeatmap::draw(
    ht,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    merge_legends = TRUE,
    padding = grid::unit(c(5, 5, 5, 5), "mm")
  )
}

save_complex_heatmap <- function(ht, file_base, width_mm, height_mm) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  dir.create(dirname(file_base), recursive = TRUE, showWarnings = FALSE)
  
  render_device <- function(label, open_device) {
    device_before <- grDevices::dev.cur()
    tryCatch(
      {
        open_device()
        draw_heatmap_object(ht)
        grDevices::dev.off()
        invisible(TRUE)
      },
      error = function(e) {
        if (grDevices::dev.cur() != device_before) {
          try(grDevices::dev.off(), silent = TRUE)
        }
        warning(
          label, " export failed for ", file_base, ": ",
          conditionMessage(e)
        )
        invisible(FALSE)
      }
    )
  }
  
  if (FUNCTIONAL_WRITE_PDF) {
    render_device("PDF", function() {
      grDevices::pdf(
        paste0(file_base, ".pdf"),
        width = width_in,
        height = height_in,
        useDingbats = FALSE,
        onefile = TRUE
      )
    })
  }
  
  if (FUNCTIONAL_WRITE_SVG) {
    render_device("SVG", function() {
      svglite::svglite(
        paste0(file_base, ".svg"),
        width = width_in,
        height = height_in,
        bg = "white"
      )
    })
  }
  
  if (FUNCTIONAL_WRITE_PNG) {
    render_device("PNG", function() {
      grDevices::png(
        paste0(file_base, ".png"),
        width = width_in,
        height = height_in,
        units = "in",
        res = FUNCTIONAL_RASTER_DPI,
        type = "cairo-png",
        bg = "white"
      )
    })
  }
  
  if (FUNCTIONAL_WRITE_TIFF) {
    render_device("TIFF", function() {
      grDevices::tiff(
        paste0(file_base, ".tiff"),
        width = width_in,
        height = height_in,
        units = "in",
        res = FUNCTIONAL_RASTER_DPI,
        compression = "lzw",
        type = "cairo",
        bg = "white"
      )
    })
  }
}

make_top_annotation <- function() {
  genotype <- c(rep("ROBILA", 3), rep("PRESTO", 3))
  condition <- c("NH", "HS-7", "HS-2", "NH", "HS-7", "HS-2")
  
  ComplexHeatmap::HeatmapAnnotation(
    Genotype = genotype,
    Condition = condition,
    col = list(
      Genotype = geno_colors,
      Condition = c(
        "NH" = unname(temp_colors["T0"]),
        "HS-7" = unname(temp_colors["T1"]),
        "HS-2" = unname(temp_colors["T2"])
      )
    ),
    annotation_name_gp = grid::gpar(fontsize = 8),
    simple_anno_size = grid::unit(3.5, "mm"),
    show_legend = TRUE
  )
}

make_heatmap_pair <- function(
    zmat,
    status_mat,
    labels,
    title,
    row_split = NULL,
    cluster_rows = TRUE
) {
  stopifnot(identical(rownames(zmat), rownames(status_mat)))
  stopifnot(all(rownames(zmat) %in% names(labels)))
  
  zplot <- pmax(
    pmin(zmat, FUNCTIONAL_Z_LIMIT),
    -FUNCTIONAL_Z_LIMIT
  )
  
  z_col_fun <- circlize::colorRamp2(
    c(-FUNCTIONAL_Z_LIMIT, 0, FUNCTIONAL_Z_LIMIT),
    c("#2166AC", "#F7F7F7", "#B2182B")
  )
  
  hm_expr <- ComplexHeatmap::Heatmap(
    zplot,
    name = "Row z-score",
    col = z_col_fun,
    na_col = "grey90",
    top_annotation = make_top_annotation(),
    column_title = title,
    column_title_gp = grid::gpar(fontsize = 11, fontface = "bold"),
    cluster_columns = FALSE,
    cluster_rows = cluster_rows,
    cluster_row_slices = FALSE,
    row_split = row_split,
    row_title_gp = grid::gpar(fontsize = 8.5, fontface = "bold"),
    row_title_rot = 0,
    row_labels = labels,
    show_row_names = TRUE,
    row_names_gp = grid::gpar(fontsize = 7.0),
    row_names_max_width = grid::unit(82, "mm"),
    column_names_rot = 45,
    column_names_gp = grid::gpar(fontsize = 8),
    column_names_centered = TRUE,
    rect_gp = grid::gpar(col = "white", lwd = 0.35),
    border = TRUE,
    heatmap_legend_param = list(
      title = "Expression",
      at = c(-FUNCTIONAL_Z_LIMIT, 0, FUNCTIONAL_Z_LIMIT),
      labels = c(
        paste0("<= -", FUNCTIONAL_Z_LIMIT),
        "0",
        paste0(">= ", FUNCTIONAL_Z_LIMIT)
      )
    )
  )
  
  hm_de <- ComplexHeatmap::Heatmap(
    status_mat,
    name = "DE status",
    col = c(
      "Down" = unname(deg_colors["Down"]),
      "NS" = unname(deg_colors["NS"]),
      "Up" = unname(deg_colors["Up"])
    ),
    na_col = unname(deg_colors["NS"]),
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    show_row_names = FALSE,
    column_names_rot = 45,
    column_names_gp = grid::gpar(fontsize = 7),
    rect_gp = grid::gpar(col = "white", lwd = 0.35),
    border = TRUE,
    width = grid::unit(35, "mm"),
    heatmap_legend_param = list(
      title = "Differential\nexpression",
      at = c("Up", "Down", "NS")
    )
  )
  
  hm_expr + hm_de
}

##############################
## 6) Input validation and output directories
##############################
FUNCTIONAL_INPUT_DIR <- validate_input_root(FUNCTIONAL_INPUT_DIR)

functional_out <- file.path(
  base_out,
  paste0(
    "Functional_KEGG_family_heatmaps_",
    SAMPLING_TIME_CHOSEN,
    "_",
    timestamp
  )
)

dir_tables <- file.path(functional_out, "00_tables")
dir_families <- file.path(functional_out, "01_family_heatmaps")
dir_combined <- file.path(functional_out, "02_combined_figure")
dir_qc <- file.path(functional_out, "03_QC")

for (d in c(
  functional_out,
  dir_tables,
  dir_families,
  dir_combined,
  dir_qc
)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

log_file <- file.path(functional_out, "Functional_KEGG_pipeline.log")
if (file.exists(log_file)) file.remove(log_file)

log_line <- function(...) {
  txt <- paste0(..., collapse = "")
  cat(txt, "\n")
  cat(txt, "\n", file = log_file, append = TRUE)
}

log_line("START functional KEGG-family pipeline")
log_line("Input root: ", FUNCTIONAL_INPUT_DIR)
log_line("Output root: ", functional_out)

##############################
## 7) Read and harmonise all family datasets
##############################
expression_columns <- c(
  "ROBILA_T0", "ROBILA_T1", "ROBILA_T2",
  "PRESTO_T0", "PRESTO_T1", "PRESTO_T2"
)

plot_column_order <- c(
  "ROBILA_T0", "ROBILA_T1", "ROBILA_T2",
  "PRESTO_T0", "PRESTO_T1", "PRESTO_T2"
)

plot_column_labels <- c(
  "ROBILA NH", "ROBILA HS-7", "ROBILA HS-2",
  "PRESTO NH", "PRESTO HS-7", "PRESTO HS-2"
)

contrast_order <- c(
  "ROBILA_T1_vs_T0", "ROBILA_T2_vs_T0",
  "PRESTO_T1_vs_T0", "PRESTO_T2_vs_T0"
)

contrast_labels <- c(
  "ROBILA HS-7", "ROBILA HS-2",
  "PRESTO HS-7", "PRESTO HS-2"
)

family_objects <- list()
curation_list <- list()
selected_de_list <- list()

for (family_dir in FUNCTIONAL_FAMILY_DIRS) {
  ff <- expected_family_files(FUNCTIONAL_INPUT_DIR, family_dir)
  family_label <- unname(FUNCTIONAL_FAMILY_LABELS[family_dir])
  
  log_line("Reading: ", family_label)
  
  expr <- read_csv_flexible(ff$expression)
  de <- read_csv_flexible(ff$de)
  subfam <- read_csv_flexible(ff$subfamily)
  provenance <- read_csv_flexible(ff$provenance)
  
  assert_columns(
    expr,
    c("my_gene_id", "gene_ID", "func_clean", expression_columns),
    paste0(family_dir, " expression table")
  )
  assert_columns(
    de,
    c("Gene", "Contrast", "logFC", "FDR"),
    paste0(family_dir, " DE table")
  )
  assert_columns(
    subfam,
    c("Gene", "Subfamily"),
    paste0(family_dir, " subfamily table")
  )
  assert_columns(
    provenance,
    c("my_gene_id", "Pathway", "pathway_num", "func_clean"),
    paste0(family_dir, " KEGG provenance table")
  )
  
  expr <- expr |>
    dplyr::mutate(
      my_gene_id = trimws(as.character(my_gene_id)),
      gene_ID = trimws(as.character(gene_ID)),
      func_clean = clean_text(func_clean)
    ) |>
    dplyr::distinct(my_gene_id, .keep_all = TRUE)
  
  de <- de |>
    dplyr::mutate(
      Gene = trimws(as.character(Gene)),
      Contrast = trimws(as.character(Contrast)),
      logFC = as.numeric(logFC),
      FDR = as.numeric(FDR)
    )
  
  if (!"Status" %in% colnames(de) || !USE_EXISTING_DE_STATUS) {
    de <- de |>
      dplyr::mutate(
        Status = dplyr::case_when(
          is.finite(FDR) &
            FDR <= FUNCTIONAL_FDR_CUTOFF &
            logFC >= FUNCTIONAL_LFC_CUTOFF ~ "Up",
          is.finite(FDR) &
            FDR <= FUNCTIONAL_FDR_CUTOFF &
            logFC <= -FUNCTIONAL_LFC_CUTOFF ~ "Down",
          TRUE ~ "NS"
        )
      )
  } else {
    de <- de |>
      dplyr::mutate(
        Status = stringr::str_to_title(trimws(as.character(Status))),
        Status = dplyr::if_else(
          Status %in% c("Up", "Down"),
          Status,
          "NS"
        )
      )
  }
  
  missing_contrasts <- setdiff(contrast_order, unique(de$Contrast))
  if (length(missing_contrasts) > 0L) {
    stop(
      family_label,
      " is missing expected contrast(s): ",
      paste(missing_contrasts, collapse = ", ")
    )
  }
  
  selected_stats <- de |>
    # Selection must be based ONLY on the four intended heat-vs-NH contrasts.
    dplyr::filter(Contrast %in% contrast_order) |>
    dplyr::group_by(Gene) |>
    dplyr::summarise(
      n_DE_contrasts = sum(Status %in% c("Up", "Down"), na.rm = TRUE),
      min_FDR = suppressWarnings(
        min(FDR[Status %in% c("Up", "Down")], na.rm = TRUE)
      ),
      max_abs_logFC = suppressWarnings(
        max(abs(logFC[Status %in% c("Up", "Down")]), na.rm = TRUE)
      ),
      response_direction = dplyr::case_when(
        any(Status == "Up", na.rm = TRUE) &
          any(Status == "Down", na.rm = TRUE) ~ "Mixed",
        any(Status == "Up", na.rm = TRUE) ~ "Up only",
        any(Status == "Down", na.rm = TRUE) ~ "Down only",
        TRUE ~ "NS"
      ),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      min_FDR = dplyr::if_else(
        is.infinite(min_FDR),
        NA_real_,
        min_FDR
      ),
      max_abs_logFC = dplyr::if_else(
        is.infinite(max_abs_logFC),
        NA_real_,
        max_abs_logFC
      )
    ) |>
    dplyr::filter(n_DE_contrasts >= 1L)
  
  selected_genes <- selected_stats$Gene
  
  subfam2 <- subfam |>
    dplyr::transmute(
      Gene = trimws(as.character(Gene)),
      Subfamily_source = normalise_subfamily(Subfamily)
    ) |>
    dplyr::distinct(Gene, .keep_all = TRUE)
  
  provenance2 <- provenance |>
    dplyr::mutate(
      my_gene_id = trimws(as.character(my_gene_id)),
      Pathway = trimws(as.character(Pathway)),
      pathway_num = trimws(as.character(pathway_num)),
      provenance_func_clean = clean_text(func_clean)
    ) |>
    dplyr::group_by(my_gene_id) |>
    dplyr::summarise(
      KEGG_Pathways = paste(
        sort(unique(Pathway[nzchar(Pathway)])),
        collapse = ";"
      ),
      KEGG_pathway_num = paste(
        sort(unique(pathway_num[nzchar(pathway_num)])),
        collapse = ";"
      ),
      provenance_func_clean = {
        vals <- provenance_func_clean[nzchar(provenance_func_clean)]
        if (length(vals) == 0L) "" else vals[[1L]]
      },
      .groups = "drop"
    )
  
  cur <- selected_stats |>
    dplyr::left_join(
      expr |>
        dplyr::select(
          my_gene_id,
          gene_ID,
          func_clean,
          dplyr::all_of(expression_columns)
        ),
      by = c("Gene" = "my_gene_id")
    ) |>
    dplyr::left_join(subfam2, by = "Gene") |>
    dplyr::left_join(
      provenance2,
      by = c("Gene" = "my_gene_id")
    ) |>
    dplyr::mutate(
      Family_ID = family_dir,
      Family = family_label,
      func_clean = dplyr::coalesce(func_clean, ""),
      provenance_func_clean = dplyr::coalesce(provenance_func_clean, ""),
      func_clean = dplyr::if_else(
        nzchar(func_clean),
        func_clean,
        provenance_func_clean
      ),
      Subfamily_source = normalise_subfamily(Subfamily_source),
      
      # Primary classification from explicit annotation + KEGG provenance.
      Subfamily_rule = classify_rule_vector(
        family_dir,
        func_clean,
        KEGG_Pathways
      ),
      
      # Source assignment is only a fallback when no stronger rule matches.
      Subfamily_auto = dplyr::case_when(
        Subfamily_rule != UNRESOLVED_LABEL ~ Subfamily_rule,
        Subfamily_source != UNRESOLVED_LABEL ~ Subfamily_source,
        TRUE ~ UNRESOLVED_LABEL
      ),
      
      Short_name_auto = first_informative_annotation(
        func_clean,
        fallback = Subfamily_auto
      ),
      
      In_expression_matrix = Gene %in% expr$my_gene_id,
      
      Include_auto = In_expression_matrix &
        (
          INCLUDE_UNRESOLVED_IN_PLOTS |
            Subfamily_auto != UNRESOLVED_LABEL
        ),
      
      Classification_evidence = dplyr::case_when(
        Subfamily_rule != UNRESOLVED_LABEL ~
          "Literature-guided explicit gene/enzyme/KEGG rule.",
        Subfamily_source != UNRESOLVED_LABEL ~
          "No stronger rule matched; supplied KEGG-enhanced subfamily retained as secondary evidence.",
        TRUE ~
          "Unresolved after extended literature-guided classification; excluded from plots pending manual curation."
      ),
      
      Classification_confidence = dplyr::case_when(
        Subfamily_rule != UNRESOLVED_LABEL ~ "High/explicit annotation rule",
        Subfamily_source != UNRESOLVED_LABEL ~ "Medium/supplied pathway assignment",
        TRUE ~ "Low/unresolved"
      )
    )
  
  expr_selected <- cur |>
    dplyr::filter(In_expression_matrix) |>
    dplyr::select(Gene, dplyr::all_of(expression_columns))
  
  mean_mat <- as.matrix(
    expr_selected[, plot_column_order, drop = FALSE]
  )
  storage.mode(mean_mat) <- "numeric"
  rownames(mean_mat) <- expr_selected$Gene
  colnames(mean_mat) <- plot_column_labels
  zmat <- zscore_rows_safe(mean_mat)
  
  status_df <- de |>
    dplyr::filter(
      Gene %in% selected_genes,
      Contrast %in% contrast_order
    ) |>
    dplyr::select(Gene, Contrast, Status) |>
    dplyr::distinct() |>
    tidyr::complete(
      Gene = selected_genes,
      Contrast = contrast_order,
      fill = list(Status = "NS")
    ) |>
    tidyr::pivot_wider(
      names_from = Contrast,
      values_from = Status,
      values_fill = "NS"
    ) |>
    dplyr::select(Gene, dplyr::all_of(contrast_order))
  
  status_mat <- as.matrix(
    status_df[, contrast_order, drop = FALSE]
  )
  rownames(status_mat) <- status_df$Gene
  colnames(status_mat) <- contrast_labels
  status_mat[!status_mat %in% c("Up", "Down", "NS")] <- "NS"
  
  family_objects[[family_dir]] <- list(
    Family_ID = family_dir,
    Family = family_label,
    mean_logCPM = mean_mat,
    zscore = zmat,
    status = status_mat,
    DE = de,
    provenance = provenance2
  )
  
  curation_list[[family_dir]] <- cur
  
  selected_de_list[[family_dir]] <- selected_stats |>
    dplyr::mutate(Family_ID = family_dir, Family = family_label) |>
    dplyr::left_join(
      de |>
        dplyr::filter(Contrast %in% contrast_order) |>
        dplyr::select(Gene, Contrast, logFC, FDR, Status) |>
        tidyr::pivot_wider(
          names_from = Contrast,
          values_from = c(logFC, FDR, Status),
          names_glue = "{Contrast}_{.value}"
        ),
      by = "Gene"
    )
  
  log_line(
    "  DE genes: ", length(selected_genes),
    " | expression matches: ", sum(cur$In_expression_matrix)
  )
}

curation_all <- dplyr::bind_rows(curation_list) |>
  dplyr::arrange(
    match(Family_ID, FUNCTIONAL_PLOT_ORDER),
    Subfamily_auto,
    Gene
  )

curation_all <- apply_manual_overrides(
  curation_all,
  FUNCTIONAL_OVERRIDE_FILE
)

##############################
## 8) Export curation, DE, overlap and bibliography tables
##############################
readr::write_csv(
  curation_all |>
    dplyr::select(
      Family_ID,
      Family,
      Gene,
      gene_ID,
      func_clean,
      KEGG_Pathways,
      KEGG_pathway_num,
      n_DE_contrasts,
      min_FDR,
      max_abs_logFC,
      response_direction,
      Subfamily_source,
      Subfamily_rule,
      Subfamily_auto,
      Subfamily,
      Classification_confidence,
      Short_name_auto,
      Short_name,
      Label,
      Classification_evidence,
      Include_auto,
      Include,
      Curator_note,
      In_expression_matrix
    ),
  file.path(dir_tables, "Functional_family_curated.csv")
)

readr::write_csv(
  dplyr::bind_rows(selected_de_list),
  file.path(dir_tables, "Functional_DE_selected.csv")
)

# Explicit curation queue. These genes NEVER appear as an "Other" block
# in publication heatmaps unless the user manually overrides them.
readr::write_csv(
  curation_all |>
    dplyr::filter(
      Subfamily == UNRESOLVED_LABEL |
        !Include
    ) |>
    dplyr::select(
      Family_ID,
      Family,
      Gene,
      gene_ID,
      func_clean,
      KEGG_Pathways,
      KEGG_pathway_num,
      Subfamily_source,
      Subfamily_rule,
      Subfamily_auto,
      Subfamily,
      Classification_confidence,
      Classification_evidence,
      Include,
      Curator_note
    ),
  file.path(
    dir_tables,
    "Functional_UNRESOLVED_excluded_for_manual_curation.csv"
  )
)

readr::write_csv(
  curation_all |>
    dplyr::filter(Include) |>
    dplyr::count(Family_ID, Family, Subfamily, name = "n_genes") |>
    dplyr::arrange(
      match(Family_ID, FUNCTIONAL_PLOT_ORDER),
      Subfamily
    ),
  file.path(dir_tables, "Functional_subfamily_summary.csv")
)


# How much of the previous "Other"/unresolved pool was rescued?
classification_audit <- curation_all |>
  dplyr::mutate(
    Source_was_unresolved = Subfamily_source == UNRESOLVED_LABEL,
    Rule_resolved = Subfamily_rule != UNRESOLVED_LABEL,
    Final_resolved = Subfamily != UNRESOLVED_LABEL
  ) |>
  dplyr::group_by(Family_ID, Family) |>
  dplyr::summarise(
    n_selected_DE = dplyr::n(),
    n_source_unresolved = sum(Source_was_unresolved, na.rm = TRUE),
    n_rescued_by_rule = sum(Source_was_unresolved & Rule_resolved, na.rm = TRUE),
    n_final_unresolved = sum(!Final_resolved, na.rm = TRUE),
    pct_source_unresolved_rescued = dplyr::if_else(
      n_source_unresolved > 0,
      100 * n_rescued_by_rule / n_source_unresolved,
      NA_real_
    ),
    .groups = "drop"
  ) |>
  dplyr::arrange(
    match(Family_ID, FUNCTIONAL_PLOT_ORDER)
  )

readr::write_csv(
  classification_audit,
  file.path(
    dir_qc,
    "Functional_classification_audit_other_rescue.csv"
  )
)

# KEGG-derived categories are intentionally non-exclusive. This table reports
# genes occurring in more than one broad family so that overlaps are explicit.
overlap_table <- curation_all |>
  dplyr::filter(Include) |>
  dplyr::group_by(Gene) |>
  dplyr::summarise(
    n_families = dplyr::n_distinct(Family_ID),
    Families = paste(unique(Family), collapse = "; "),
    Family_IDs = paste(unique(Family_ID), collapse = ";"),
    .groups = "drop"
  ) |>
  dplyr::filter(n_families > 1L) |>
  dplyr::arrange(dplyr::desc(n_families), Gene)

readr::write_csv(
  overlap_table,
  file.path(dir_tables, "Functional_gene_overlap_across_families.csv")
)

bibliography <- tibble::tribble(
  ~Family_ID, ~Module, ~Classification_basis, ~Reference, ~Identifier,
  
  "General_KEGG",
  "Pathway provenance",
  "KEGG provenance is supporting evidence; broad KEGG pathways are non-exclusive and are not used alone to assign a precise enzyme function.",
  "Kanehisa et al. — KEGG pathway framework",
  "KEGG",
  
  "01_Protein_Folding_Processing_Degradation",
  "Proteasome / stress granules",
  "26S proteasome components are separated from chaperones and autophagy because proteasomes directly participate in heat-stress-granule turnover.",
  "Xie et al. 2024, Molecular Cell",
  "DOI:10.1016/j.molcel.2024.07.033",
  
  "01_Protein_Folding_Processing_Degradation",
  "Selective autophagy",
  "ATG machinery and selective autophagy are treated as a distinct proteostasis arm for heat-damaged protein/aggregate clearance.",
  "Sanchez de Medina Hernandez et al. 2025, Developmental Cell",
  "DOI:10.1016/j.devcel.2025.11.001",
  
  "01_Protein_Folding_Processing_Degradation",
  "Hsf–BAG–autophagy",
  "BAG co-chaperones and ATG components are separated because BAG5b links HsfA1a to autophagy-mediated thermotolerance.",
  "He et al. 2025, Plant Physiology",
  "DOI:10.1093/plphys/kiaf526",
  
  "01_Protein_Folding_Processing_Degradation",
  "UPR / translation",
  "ER-UPR, ER folding and ERAD are separated from cytosolic HSP families.",
  "Dannfald et al. 2025, Plant Physiology",
  "PMID:39688875",
  
  "02_Hormone_Signaling_Metabolism",
  "ABA core signalling",
  "ABA is separated into biosynthesis/homeostasis/transport/PYR-PYL perception/PP2C-SnRK2/transcriptional output.",
  "SPL1/SPL12-PYL thermotolerance study",
  "PMID:28400323",
  
  "02_Hormone_Signaling_Metabolism",
  "Ethylene",
  "ACS/ACO biosynthesis is separated from receptor–CTR1–EIN2–EIN3/ERF signalling.",
  "ERF95/ERF97 heat-response study",
  "PMID:33793870",
  
  "02_Hormone_Signaling_Metabolism",
  "Brassinosteroid",
  "BR biosynthesis is separated from BRI1/BAK1-BIN2-BES1/BZR signalling.",
  "Setsungnern et al. 2020, Plant Science",
  "DOI:10.1016/j.plantsci.2020.110470",
  
  "02_Hormone_Signaling_Metabolism",
  "Jasmonate",
  "LOX/AOS/AOC/OPR biosynthesis and COI1/JAZ signalling are distinct axes.",
  "HSFA2-LOX3/OPR3-JA heat study",
  "PMID:38556243",
  
  "02_Hormone_Signaling_Metabolism",
  "Cytokinin",
  "IPT/LOG biosynthesis, CKX turnover, transport and AHK/AHP/ARR two-component signalling are separated.",
  "Prerostova et al. 2020",
  "PMID:32133021",
  
  "02_Hormone_Signaling_Metabolism",
  "Strigolactone",
  "D27/CCD7/CCD8/MAX1 biosynthesis is separated from D14/MAX2/SMXL signalling.",
  "AtMYBS1-MAX1-D14 heat study",
  "PMID:37608548",
  
  "03_MAPK_Calcium_Stress_Signaling",
  "Thermal Ca2+ influx",
  "CNGC channels are separated from CaM/CML, CPK and CBL/CIPK decoding modules.",
  "Finka et al. 2012, The Plant Cell",
  "DOI:10.1105/tpc.112.095844",
  
  "03_MAPK_Calcium_Stress_Signaling",
  "MAPK",
  "MAPKKK, MKK and MPK tiers are separated to preserve signalling direction.",
  "MPK6-HsfA2 heat-response study",
  "PMID:23638397",
  
  "03_MAPK_Calcium_Stress_Signaling",
  "Phospholipid signalling",
  "PLD, PLC and DGK are separated because heat rapidly remodels PA/PIP2 signalling.",
  "Mishkind et al. 2009, The Plant Journal",
  "DOI:10.1111/j.1365-313X.2009.03933.x",
  
  "03_MAPK_Calcium_Stress_Signaling",
  "Heat memory / phosphatidic acid",
  "PLDalpha2 and PP2C FGT2 connect membrane-lipid signalling to heat memory.",
  "Urrea Castellanos et al. 2020, The Plant Journal",
  "DOI:10.1111/tpj.14927",
  
  "04_Photosynthesis_Energy_Metabolism",
  "Photosynthetic electron transport",
  "LHCII, PSII core, OEC, cytochrome b6f, PSI, Fd/FNR, CEF and ATP synthase are separated rather than pooled as photosynthesis.",
  "Daily heat stress / cytochrome b6f study",
  "PMID:42418814",
  
  "04_Photosynthesis_Energy_Metabolism",
  "Thylakoid membrane dependency",
  "Photosynthetic function is interpreted together with galactolipid status under heat.",
  "Arabidopsis DGDG heat study",
  "PMID:21458884",
  
  "05_Terpenoid_Secondary_Metabolism",
  "Phenylpropanoid/flavonoid",
  "PAL/C4H/4CL, hydroxycinnamate, lignin, flavonoid, flavonol and anthocyanin branches are separated.",
  "Carrot phenylpropanoid heat-protection study",
  "PMID:27713760",
  
  "05_Terpenoid_Secondary_Metabolism",
  "Heat multi-omics",
  "Phenylpropanoid and flavonoid branches can move differently under heat and should not be merged into one secondary-metabolism class.",
  "Rosa hybrida heat transcriptome/metabolome study",
  "PMID:39304829",
  
  "05_Terpenoid_Secondary_Metabolism",
  "Carotenoid",
  "Carotenoid biosynthesis is kept distinct from MVA/MEP and terpene-synthase branches.",
  "ABA-mediated carotenoid reprogramming under heat",
  "PMID:40935072",
  
  "06_ROS_Redox_Glutathione",
  "Glutathione",
  "GSH biosynthesis, GST, GPX and the AsA-GSH cycle are separated because glutathione status is functionally required for high-temperature adaptation.",
  "Glutathione-mediated heat responses in Arabidopsis",
  "DOI:10.1093/jxb/erad042",
  
  "06_ROS_Redox_Glutathione",
  "Catalase",
  "Catalase is separated from other ROS enzymes because CAT2 has a specific role in long-duration heat.",
  "CAT2 long-term heat study",
  "DOI:10.1016/j.bbrc.2020.11.006",
  
  "06_ROS_Redox_Glutathione",
  "Ascorbate–glutathione cycle",
  "APX, MDHAR, DHAR and GR are represented as distinct biochemical steps.",
  "Heat-shock compartmental AsA-GSH study",
  "PMID:19236663",
  
  "07_Lipid_Membrane_Metabolism",
  "Global lipid remodelling",
  "FA synthesis/desaturation, phospholipids, galactolipids, TAG, sphingolipids and sterols are separated because they respond differently to heat.",
  "Arabidopsis leaf lipid heat-remodelling study",
  "PMID:32635518",
  
  "07_Lipid_Membrane_Metabolism",
  "Thylakoid galactolipids",
  "MGDG, DGDG and SQDG are split because thylakoid lipid composition directly affects heat performance.",
  "TaMGD heat-tolerance study",
  "PMID:39837009",
  
  "07_Lipid_Membrane_Metabolism",
  "TAG storage",
  "DGAT/PDAT/TAG is separated from structural membrane lipids; TAG accumulation can sequester heat-generated lipid intermediates.",
  "DGAT2/TAG heat study",
  "PMID:38503190",
  
  "08_Carbon_Starch_Sucrose_Osmoprotection",
  "T6P / carbon allocation",
  "TPS/TPP/trehalose are separated from sucrose and starch because T6P is a carbon-signalling node in temperature responses.",
  "TPS1 / raffinose / thermotolerance study",
  "PMID:37427798",
  
  "08_Carbon_Starch_Sucrose_Osmoprotection",
  "Thermoresponsive growth",
  "T6P is treated as a signalling/metabolic module rather than simply an osmolyte.",
  "T6P-KIN10-PIF4 study",
  "PMID:31393060",
  
  "08_Carbon_Starch_Sucrose_Osmoprotection",
  "Compatible solutes",
  "Proline, raffinose-family oligosaccharides and polyols are separated from central carbon metabolism.",
  "Tomato compatible-solute heat study",
  "PMID:25747289"
)

readr::write_csv(
  bibliography,
  file.path(dir_tables, "Functional_classification_literature_rules.csv")
)

##############################
## 9) One figure per broad functional family
##############################
combined_z <- list()
combined_status <- list()
combined_meta <- list()

for (family_dir in FUNCTIONAL_FAMILY_DIRS) {
  obj <- family_objects[[family_dir]]
  family_label <- obj$Family
  
  fam_tbl <- curation_all |>
    dplyr::filter(
      Family_ID == family_dir,
      Include,
      Subfamily != UNRESOLVED_LABEL
    ) |>
    dplyr::mutate(
      Subfamily = factor(
        Subfamily,
        levels = order_subfamilies_for_plot(
          family_dir,
          Subfamily
        )
      )
    ) |>
    dplyr::arrange(Subfamily, Gene)
  
  if (nrow(fam_tbl) == 0L) {
    warning("No genes retained for family: ", family_label)
    next
  }
  
  genes <- fam_tbl$Gene
  genes <- genes[
    genes %in% rownames(obj$zscore) &
      genes %in% rownames(obj$status)
  ]
  fam_tbl <- fam_tbl |>
    dplyr::filter(Gene %in% genes) |>
    dplyr::arrange(Subfamily, Gene)
  genes <- fam_tbl$Gene
  
  z_fam <- obj$zscore[genes, , drop = FALSE]
  st_fam <- obj$status[genes, , drop = FALSE]
  
  # Biological pathway order, not alphabetical order.
  subfamily_levels <- order_subfamilies_for_plot(
    family_dir,
    fam_tbl$Subfamily
  )
  
  row_split <- NULL
  if (FUNCTIONAL_SPLIT_ROWS_BY_SUBFAMILY) {
    row_split <- factor(
      fam_tbl$Subfamily,
      levels = subfamily_levels
    )
  }
  
  labels <- stats::setNames(fam_tbl$Label, fam_tbl$Gene)
  
  ht <- make_heatmap_pair(
    zmat = z_fam,
    status_mat = st_fam,
    labels = labels,
    title = paste0(
      family_label,
      " — DE in at least one heat-stress contrast"
    ),
    row_split = row_split,
    cluster_rows = FUNCTIONAL_CLUSTER_ROWS
  )
  
  n <- length(genes)
  n_slices <- if (is.null(row_split)) 1L else nlevels(droplevels(row_split))
  height_mm <- max(
    105,
    48 + FUNCTIONAL_ROW_MM * n + 5 * n_slices
  )
  
  file_base <- file.path(
    dir_families,
    paste0(
      family_dir,
      "_",
      safe_file_stub(family_label),
      "_heatmap"
    )
  )
  
  save_complex_heatmap(
    ht,
    file_base,
    width_mm = FUNCTIONAL_FIG_WIDTH_MM,
    height_mm = height_mm
  )
  
  readr::write_csv(
    tibble::rownames_to_column(as.data.frame(obj$mean_logCPM[genes, , drop = FALSE]), "Gene") |>
      dplyr::left_join(
        fam_tbl |>
          dplyr::select(Gene, gene_ID, Subfamily, Short_name, Label),
        by = "Gene"
      ),
    file.path(
      dir_tables,
      paste0(family_dir, "_meanLogCPM_selected.csv")
    )
  )
  
  readr::write_csv(
    tibble::rownames_to_column(as.data.frame(z_fam), "Gene") |>
      dplyr::left_join(
        fam_tbl |>
          dplyr::select(Gene, gene_ID, Subfamily, Short_name, Label),
        by = "Gene"
      ),
    file.path(
      dir_tables,
      paste0(family_dir, "_zscore_selected.csv")
    )
  )
  
  # Unique row keys are required because KEGG-derived families overlap.
  row_keys <- paste0(family_dir, "::", genes)
  rownames(z_fam) <- row_keys
  rownames(st_fam) <- row_keys
  
  combined_z[[family_dir]] <- z_fam
  combined_status[[family_dir]] <- st_fam
  combined_meta[[family_dir]] <- fam_tbl |>
    dplyr::mutate(Row_key = row_keys)
  
  log_line("Saved family: ", family_label, " | genes: ", n)
}

##############################
## 10) Optional combined figure
##############################
if (FUNCTIONAL_WRITE_COMBINED && length(combined_z) > 0L) {
  combined_tbl <- dplyr::bind_rows(combined_meta) |>
    dplyr::arrange(
      match(Family_ID, FUNCTIONAL_PLOT_ORDER),
      Subfamily,
      Gene
    )
  
  if (nrow(combined_tbl) <= FUNCTIONAL_MAX_COMBINED_GENES) {
    z_combined <- do.call(rbind, combined_z)
    st_combined <- do.call(rbind, combined_status)
    
    keys <- combined_tbl$Row_key
    z_combined <- z_combined[keys, , drop = FALSE]
    st_combined <- st_combined[keys, , drop = FALSE]
    
    split_labels <- paste(
      combined_tbl$Family,
      combined_tbl$Subfamily,
      sep = " — "
    )
    row_split <- factor(
      split_labels,
      levels = unique(split_labels)
    )
    
    labels <- stats::setNames(combined_tbl$Label, combined_tbl$Row_key)
    
    ht_combined <- make_heatmap_pair(
      zmat = z_combined,
      status_mat = st_combined,
      labels = labels,
      title = "Heat-responsive functional modules in carrot leaves",
      row_split = row_split,
      cluster_rows = FUNCTIONAL_CLUSTER_ROWS
    )
    
    n_slices <- nlevels(droplevels(row_split))
    combined_height_mm <- max(
      220,
      58 + FUNCTIONAL_ROW_MM * nrow(combined_tbl) + 4 * n_slices
    )
    
    save_complex_heatmap(
      ht_combined,
      file.path(
        dir_combined,
        "Functional_KEGG_families_combined"
      ),
      width_mm = FUNCTIONAL_FIG_WIDTH_MM + 12,
      height_mm = combined_height_mm
    )
    
    log_line("Saved combined figure | rows: ", nrow(combined_tbl))
  } else {
    warning(
      "Combined figure skipped because ", nrow(combined_tbl),
      " rows exceed FUNCTIONAL_MAX_COMBINED_GENES = ",
      FUNCTIONAL_MAX_COMBINED_GENES,
      ". Individual-family figures were still generated."
    )
    log_line(
      "Combined figure skipped | rows: ", nrow(combined_tbl),
      " | maximum: ", FUNCTIONAL_MAX_COMBINED_GENES
    )
  }
}

##############################
## 11) Reproducibility record
##############################
saveRDS(
  list(
    input_root = FUNCTIONAL_INPUT_DIR,
    curation = curation_all,
    family_objects = family_objects,
    overlap = overlap_table,
    bibliography = bibliography
  ),
  file.path(dir_qc, "Functional_KEGG_analysis_objects.rds")
)

capture.output(
  utils::sessionInfo(),
  file = file.path(functional_out, "sessionInfo.txt")
)

cat("\nDONE functional KEGG-family pipeline\n")
cat("Input directory:", FUNCTIONAL_INPUT_DIR, "\n")
cat("Output directory:", functional_out, "\n")
cat("Formats: PDF, SVG, PNG and TIFF\n")
cat("Expression order: ROBILA NH | HS-7 | HS-2 | PRESTO NH | HS-7 | HS-2\n")
cat("Unresolved category shown in figures: NO\n")
cat("Manual override file:", FUNCTIONAL_OVERRIDE_FILE, "\n")
