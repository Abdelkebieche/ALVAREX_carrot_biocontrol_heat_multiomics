# NOTE: companion annotation-builder reference implementation.
# Paths have been sanitized. See docs/05_functional_dictionaries.md before running.

############################################################
## 09A_HSP_family_heatmaps.R
## Figure 9A — Differentially expressed chaperones/HSPs
##
## Input model:
##   - same P1 merged data structure as 01_load_and_prepare_data_P1_merged.R
##   - one DGEList containing ROBILA + PRESTO
##   - cell-means design: ~ 0 + GenoTemp (+ Block)
##
## Selection:
##   - gene is retained when Status is Up or Down in >= 1 of the four
##     genotype-specific heat-vs-NH contrasts only
##
## Expression heatmap:
##   - complete TMM logCPM matrix rebuilt from counts
##   - means calculated for the six Genotype x Temperature groups
##   - row-wise z-score used for colour
##   - fixed biological column order; no column clustering
##
## HSP classification:
##   - conservative annotation/domain keyword rules
##   - known false-positive domain matches are excluded
##   - ambiguous cases are exported for manual curation
##
## Outputs:
##   - one heatmap per HSP family
##   - combined Figure 9A
##   - PDF, SVG, PNG and TIFF
##   - curation, exclusion, expression and summary tables
##   - fixed plotting order:
##       ROBILA NH -> ROBILA HS-7 -> ROBILA HS-2 ->
##       PRESTO NH -> PRESTO HS-7 -> PRESTO HS-2
############################################################

##############################
## 0) Configuration
##############################
source("00_builder_config.R")
source("00_builder_helpers.R")

# Directory containing the two HSP input tables.
# It can be overridden without editing this script:
# export HSP_INPUT_DIR=/path/to/HSP_Figure_9A_inputs
# Dossier contenant les fichiers HSP
HSP_INPUT_DIR <- dirname(HSP_ANNOT_FILE)

# Fichiers d’entrée
HSP_ANNOTATION_FILE <- file.path(
  HSP_INPUT_DIR,
  "HSP_all_genes.csv"
)

HSP_DE_FILE <- file.path(
  HSP_INPUT_DIR,
  "HSP_DE_per_contrast.csv"
)

# Fichier de correction manuelle
HSP_OVERRIDE_FILE <- file.path(
  HSP_INPUT_DIR,
  "HSP_manual_overrides.csv"
)
# Optional manual override file. The script creates a template on first run.  
HSP_OVERRIDE_FILE <- file.path(HSP_INPUT_DIR, "HSP_manual_overrides.csv")

# Safety parameter: if more than one Treatment remains after P1 filtering,
# specify the treatment to plot, e.g. "EAU". Leave NULL only when one level exists.
HSP_TREATMENT_CHOSEN <- NULL

# Selection follows the already-computed DE Status column by default.
# Set FALSE only if Status must be reconstructed from FDR and logFC.
USE_EXISTING_DE_STATUS <- TRUE
HSP_FDR_CUTOFF <- PADJ_CUTOFF
HSP_LFC_CUTOFF <- lfc_threshold

# Figure settings
HSP_Z_LIMIT       <- 2.5
HSP_CLUSTER_ROWS  <- TRUE
HSP_ROW_MM        <- 3.7
HSP_FIG_WIDTH_MM  <- 215
HSP_RASTER_DPI    <- 600
HSP_WRITE_PNG     <- TRUE
HSP_WRITE_TIFF    <- TRUE
HSP_WRITE_PDF     <- TRUE
HSP_WRITE_SVG     <- TRUE

# Combined-figure composition:
# TRUE  = include explicitly annotated HSP70/HSP90 co-chaperones in the combined figure.
# FALSE = restrict the combined figure to canonical HSP/chaperone families.
# Co-chaperones are always retained in the curation/QC outputs.
INCLUDE_COCHAPERONES_IN_COMBINED <- TRUE

set.seed(SEED)
options(stringsAsFactors = FALSE)

##############################
## 1) Packages
##############################
load_pkgs(c("edgeR", "ComplexHeatmap"),
          install_missing = INSTALL_MISSING_PKGS, bioc = TRUE)
load_pkgs(c("dplyr", "tidyr", "tibble", "readr", "stringr",
            "circlize", "svglite"),
          install_missing = INSTALL_MISSING_PKGS, bioc = FALSE)

##############################
## 2) Utility functions
##############################
`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || all(is.na(x))) y else x
}

assert_columns <- function(x, required, object_name) {
  missing_cols <- setdiff(required, colnames(x))
  if (length(missing_cols) > 0L) {
    stop(object_name, " is missing required column(s): ",
         paste(missing_cols, collapse = ", "))
  }
}

zscore_rows_safe <- function(mat) {
  mat <- as.matrix(mat)
  mu  <- rowMeans(mat, na.rm = TRUE)
  sds <- apply(mat, 1L, stats::sd, na.rm = TRUE)
  out <- sweep(mat, 1L, mu, "-")
  valid <- is.finite(sds) & sds > 0
  out[valid, ] <- sweep(out[valid, , drop = FALSE], 1L, sds[valid], "/")
  out[!valid, ] <- 0
  out[!is.finite(out)] <- 0
  out
}

first_existing <- function(paths) {
  ok <- paths[file.exists(paths)]
  if (length(ok) == 0L) return(NA_character_)
  normalizePath(ok[[1L]], winslash = "/", mustWork = TRUE)
}

read_csv_flexible <- function(path) {
  if (!file.exists(path)) stop("Input file not found: ", path)
  x <- tryCatch(
    readr::read_csv(path, show_col_types = FALSE, progress = FALSE),
    error = function(e) NULL
  )
  if (is.null(x) || ncol(x) <= 1L) {
    x <- readr::read_delim(path, delim = ";", show_col_types = FALSE,
                           progress = FALSE)
  }
  as.data.frame(x, stringsAsFactors = FALSE)
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

##############################
## 3) HSP classification rules
##############################
# Literature basis for the conservative rules below:
# - sHSPs: Hsp20/alpha-crystallin domain is necessary but can also occur in
#   non-canonical ACD proteins; catalytic multidomain hits are therefore flagged.
#   Waters et al. 2008, DOI: 10.1007/s12192-008-0023-7
#   Bondino et al. 2012, DOI: 10.1007/s00425-011-1575-9
# - DnaJ/J-proteins: DNAJA-like = J + zinc-finger + C-terminal domains;
#   DNAJB-like = J + C-terminal domains; DNAJC-like = J domain with other architecture.
#   DOI: 10.1111/nph.14827; DOI: 10.1139/gen-2017-0206
# - HSP70: cytosolic, organellar and ER/BiP members are separated when explicit.
#   Sung et al. 2001, DOI: 10.1104/pp.126.2.789
# - HSP90: cytosolic, plastid, mitochondrial and ER/endoplasmin isoforms.
#   Krishna & Gloor 2001, PMID: 11599565
# - Chaperonins: type I Cpn60/Cpn10 and type II CCT/TCP-1 are distinct systems.
#   Hill & Hemmingsen 2001, PMID: 11599560
# - HSP100/ClpB: cytosolic HSP101 and organellar ClpB forms.
#   Myung et al. 2006, DOI: 10.1111/j.1365-313X.2006.02940.x

classify_hsp_annotations <- function(tbl) {
  assert_columns(tbl, c("Gene", "func_clean"), "HSP annotation table")
  
  tbl |>
    dplyr::mutate(
      annotation_raw = dplyr::coalesce(as.character(func_clean), ""),
      annotation_lc  = stringr::str_to_lower(annotation_raw),
      
      # Known misleading domain-name matches: these are not sufficient to call an HSP.
      fp_groes_enzyme = stringr::str_detect(
        annotation_lc,
        "alcohol dehydrogenase groes-like|cinnamyl alcohol dehydrogenase|zinc-binding dehydrogenase|zinc-containing alcohol dehydrogenases|quinone oxidoreductase"
      ),
      fp_hsp70_nbd = stringr::str_detect(annotation_lc, "nbd_sugar-kinase_hsp70_actin") &
        !stringr::str_detect(
          annotation_lc,
          "hsp70 protein|heat shock 70|heat shock protein 70|70kda heat shock|chaperone protein dnak|molecular chaperone dnak|hspa[0-9]"
        ),
      fp_hsp90_hatpase = stringr::str_detect(
        annotation_lc,
        "histidine kinase|ethylene receptor|phytochrome|pyruvate dehydrogenase.*kinase|branched-chain alpha-ketoacid dehydrogenase kinase"
      ) & !stringr::str_detect(
        annotation_lc,
        "hsp90 protein|heat shock protein 90|heat shock 90|90kda heat shock protein signature|heat shock hsp90 proteins family|endoplasmin|grp94|chaperone protein htpg"
      ),
      ambiguous_acd_enzyme = stringr::str_detect(
        annotation_lc,
        "hsp20/alpha crystallin|small heat shock protein.*domain|acd_"
      ) & stringr::str_detect(
        annotation_lc,
        "udp-gluc|udp-glycosyl|glycosyltransferase|lipase|acylhydrolase|polysaccharide lyase"
      ),
      
      is_cochaperone = stringr::str_detect(
        annotation_lc,
        "activator of hsp90 atpase|\\baha1\\b|stress-induced phosphoprotein|hsp70/hsp90-organizing|hsp70-hsp90 organizing|\\bsti1\\b|\\bhop\\b|hsp70-binding protein|nucleotide exchange factor fes1|\\bfes1\\b|hspbp"
      ),
      is_hsp100 = stringr::str_detect(
        annotation_lc,
        "atp-dependent chaperone clpb|chaperone[_ ]clpb|chaperone protein clpb|chaperone clpb|clpa/b signature|\\bhsp100\\b|\\bhsp101\\b"
      ),
      is_hsp90 = stringr::str_detect(
        annotation_lc,
        "hsp90 protein|heat shock protein 90|heat shock 90|90kda heat shock protein signature|heat shock hsp90 proteins family|endoplasmin|grp94|chaperone protein htpg|hatpase_hsp90-like"
      ),
      is_hsp70 = stringr::str_detect(
        annotation_lc,
        "hsp70 protein|heat shock 70|heat shock protein 70|70kda heat shock protein signature|heat shock hsp70 proteins family|chaperone protein dnak|molecular chaperone dnak|luminal-binding protein|stromal 70 kda|hspa[0-9]|hsca homolog"
      ),
      is_cct = stringr::str_detect(
        annotation_lc,
        "t-complex protein|tailless complex polypeptide|tcp1_|chaperonins tcp-1|chap_cct|\\bcct\\b"
      ),
      is_hsp10 = stringr::str_detect(
        annotation_lc,
        "\\bcpn10\\b|10kda chaperonin|10 kda chaperonin|chaperonin 10"
      ),
      is_hsp60 = stringr::str_detect(
        annotation_lc,
        "chaperonin groel|\\bgroel\\b|60 kda chaperonin|60kda chaperonin|chaperonins cpn60|\\bcpn60\\b"
      ) & !is_cct,
      is_dnaj = stringr::str_detect(
        annotation_lc,
        "dnaj domain|dnaj domain signature|dnaj domain profile|\\bdnaj\\b|molecular chaperone dnaj|chaperone protein dnaj"
      ),
      is_shsp = stringr::str_detect(
        annotation_lc,
        "hsp20/alpha crystallin|small heat shock protein|small hsp|17\\.[0-9]+ kda class|26\\.[0-9]+ kda heat shock|acd_shsp|acd_schsp"
      ),
      
      Exclusion_reason = dplyr::case_when(
        fp_groes_enzyme ~ "GroES-like fold occurs in an alcohol/dehydrogenase enzyme; not evidence for HSP10.",
        fp_hsp70_nbd ~ "Generic sugar-kinase/HSP70/actin NBD match without a canonical HSP70 annotation.",
        fp_hsp90_hatpase ~ "Generic HATPase domain in a kinase/receptor/phytochrome; not evidence for HSP90.",
        ambiguous_acd_enzyme ~ "Alpha-crystallin/Hsp20 hit is combined with an unrelated catalytic-domain annotation; sequence-level validation required.",
        TRUE ~ NA_character_
      ),
      
      HSP_family_auto = dplyr::case_when(
        !is.na(Exclusion_reason) ~ "Excluded / ambiguous HSP-like",
        is_cochaperone ~ "Co-chaperones",
        is_hsp100 ~ "HSP100 / ClpB",
        is_hsp90 ~ "HSP90",
        is_hsp70 ~ "HSP70 / HSC70",
        is_cct ~ "CCT / TCP-1",
        is_hsp10 ~ "HSP10 / Cpn10",
        is_hsp60 ~ "HSP60 / Cpn60",
        is_dnaj ~ "HSP40 / DnaJ",
        is_shsp ~ "HSP20 / sHSP",
        TRUE ~ "Unclassified"
      ),
      
      HSP_subfamily_auto = dplyr::case_when(
        !is.na(Exclusion_reason) ~ "Requires validation",
        
        is_cochaperone & stringr::str_detect(annotation_lc, "aha1|activator of hsp90") ~ "Aha1",
        is_cochaperone & stringr::str_detect(annotation_lc, "sti1|stress-induced phosphoprotein|organizing|\\bhop\\b") ~ "HOP / STI1",
        is_cochaperone & stringr::str_detect(annotation_lc, "fes1|hsp70-binding|hspbp") ~ "HSP70 nucleotide-exchange factor",
        is_cochaperone ~ "Other co-chaperone",
        # HSP100/ClpB: localisation is assigned only when supported by annotation.
        # Generic ClpA/B signatures remain lower-confidence HSP100-like candidates.
        is_hsp100 & stringr::str_detect(annotation_lc, "clpb3|chloroplast|chloroplastic|plastid") ~ "ClpB3-like, chloroplast",
        is_hsp100 & stringr::str_detect(annotation_lc, "mitochond") ~ "ClpB-like, mitochondrial",
        is_hsp100 & stringr::str_detect(
          annotation_lc,
          "hsp101|atp-dependent chaperone clpb|chaperone[_ ]clpb|chaperone protein clpb"
        ) ~ "ClpB / HSP101-like",
        is_hsp100 ~ "HSP100/ClpA-B-like, localisation unresolved",
        
        # HSP90: do not infer cytosolic localisation when organellar evidence is absent.
        is_hsp90 & stringr::str_detect(annotation_lc, "endoplasmin|grp94") ~ "ER HSP90 / endoplasmin",
        is_hsp90 & stringr::str_detect(annotation_lc, "hsp90-5|chloroplast|chloroplastic|plastid") ~ "Chloroplast HSP90",
        is_hsp90 & stringr::str_detect(annotation_lc, "hsp90-6|mitochond") ~ "Mitochondrial HSP90",
        is_hsp90 & stringr::str_detect(annotation_lc, "cytosol|cytosolic") ~ "Cytosolic HSP90-like",
        is_hsp90 ~ "HSP90-like, localisation unresolved",
        
        # HSP70: DnaK/HscA homology alone is not used to infer mitochondria.
        # Organellar localisation requires explicit localisation evidence.
        is_hsp70 & stringr::str_detect(annotation_lc, "luminal-binding|hspa5|\\bbip\\b") ~ "ER BiP / HSPA5-like",
        is_hsp70 & stringr::str_detect(annotation_lc, "stromal|chloroplast|chloroplastic|plastid") ~ "Chloroplast HSP70",
        is_hsp70 & stringr::str_detect(annotation_lc, "hspa9|heat shock protein 9|mitochond") ~ "Mitochondrial HSP70-like",
        is_hsp70 & stringr::str_detect(annotation_lc, "hspa4|heat shock 70 kda protein 4") ~ "HSP110 / HSPH-like",
        is_hsp70 & stringr::str_detect(annotation_lc, "heat shock cognate|hsc70") ~ "HSC70-like",
        is_hsp70 & stringr::str_detect(annotation_lc, "cytosol|cytosolic") ~ "Cytosolic HSP70-like",
        is_hsp70 ~ "HSP70-like, localisation unresolved",
        
        is_cct & stringr::str_detect(annotation_lc, "delta") ~ "CCT delta subunit",
        is_cct ~ "CCT / TCP-1 subunit",
        
        is_hsp10 & stringr::str_detect(annotation_lc, "mitochond") ~ "Mitochondrial Cpn10",
        is_hsp10 & stringr::str_detect(annotation_lc, "chloroplast|chloroplastic|plastid") ~ "Chloroplast Cpn10",
        is_hsp10 ~ "Cpn10",
        
        is_hsp60 & stringr::str_detect(annotation_lc, "chloroplast|chloroplastic|plastid") ~ "Chloroplast Cpn60",
        is_hsp60 & stringr::str_detect(annotation_lc, "mitochond") ~ "Mitochondrial Cpn60",
        is_hsp60 ~ "Cpn60 / GroEL-like",
        
        is_dnaj & stringr::str_detect(annotation_lc, "dnaj_zf|zinc finger cr-type") &
          stringr::str_detect(annotation_lc, "dnaj_c|dnaj c terminal") ~ "DnaJ-A-like",
        is_dnaj & !stringr::str_detect(annotation_lc, "dnaj_zf|zinc finger cr-type") &
          stringr::str_detect(annotation_lc, "dnaj_c|dnaj c terminal") ~ "DnaJ-B-like",
        is_dnaj & stringr::str_detect(annotation_lc, "dnaj domain|dnaj domain signature|dnaj domain profile") ~ "DnaJ-C-like",
        is_dnaj ~ "DnaJ-related, domain incomplete",
        
        is_shsp & stringr::str_detect(annotation_lc, "class ii") ~ "Cytosolic class II sHSP",
        is_shsp & stringr::str_detect(annotation_lc, "class i") ~ "Cytosolic class I sHSP",
        is_shsp & stringr::str_detect(annotation_lc, "mitochond") ~ "Mitochondrial sHSP",
        is_shsp & stringr::str_detect(annotation_lc, "chloroplast|chloroplastic|plastid") ~ "Chloroplast sHSP",
        is_shsp ~ "sHSP / ACD protein",
        
        TRUE ~ "Unclassified"
      ),
      
      Evidence_auto = dplyr::case_when(
        !is.na(Exclusion_reason) ~ Exclusion_reason,
        is_hsp100 ~ "ClpB name plus AAA+/ClpA-B signatures.",
        is_hsp90 ~ "HSP90 name/signature or Hsp90-specific HATPase annotation.",
        is_hsp70 ~ "HSP70/DnaK name plus HSP70 family signatures.",
        is_cct ~ "Type-II chaperonin CCT/TCP-1 subunit annotation.",
        is_hsp10 ~ "Cpn10/10-kDa cochaperonin signature.",
        is_hsp60 ~ "GroEL/Cpn60/60-kDa chaperonin signatures.",
        is_dnaj ~ "J-domain/DnaJ annotation; A/B/C-like class inferred from domain architecture.",
        is_shsp ~ "Hsp20/alpha-crystallin/small-HSP domain annotation.",
        is_cochaperone ~ "Explicit HSP70/HSP90 co-chaperone annotation.",
        TRUE ~ "No sufficiently specific HSP evidence in the supplied annotation."
      ),
      
      Confidence_auto = dplyr::case_when(
        !is.na(Exclusion_reason) ~ "Low",
        is_hsp100 & stringr::str_detect(
          annotation_lc,
          "hsp101|atp-dependent chaperone clpb|chaperone[_ ]clpb|chaperone protein clpb"
        ) ~ "High",
        is_hsp100 & stringr::str_detect(annotation_lc, "clpa/b signature") ~ "Medium",
        is_hsp90 & stringr::str_detect(annotation_lc, "90kda heat shock protein signature|hsp90 protein|endoplasmin") ~ "High",
        is_hsp70 & stringr::str_detect(annotation_lc, "70kda heat shock protein signature|heat shock hsp70 proteins family|dnak") ~ "High",
        is_cct | is_hsp10 | is_hsp60 ~ "High",
        is_dnaj & stringr::str_detect(annotation_lc, "dnaj domain|dnaj domain signature|dnaj domain profile") ~ "High",
        is_dnaj ~ "Medium",
        is_shsp & stringr::str_detect(annotation_lc, "small heat shock protein.*domain|17\\.[0-9]+ kda class|26\\.[0-9]+ kda heat shock") ~ "High",
        is_shsp ~ "Medium",
        is_cochaperone ~ "High",
        TRUE ~ "Low"
      ),
      
      Include_auto = HSP_family_auto %in% c(
        "HSP20 / sHSP", "HSP40 / DnaJ", "HSP60 / Cpn60",
        "HSP10 / Cpn10", "CCT / TCP-1", "HSP70 / HSC70",
        "HSP90", "HSP100 / ClpB", "Co-chaperones"
      ),
      
      Short_name_auto = paste0(HSP_subfamily_auto, " | ", Gene)
    )
}

apply_manual_overrides <- function(tbl, override_file) {
  template <- tbl |>
    dplyr::transmute(
      Gene,
      Family_auto = HSP_family_auto,
      Subfamily_auto = HSP_subfamily_auto,
      Confidence_auto,
      Include_auto,
      Family_override = "",
      Subfamily_override = "",
      Short_name_override = "",
      Include_override = "",
      Curator_note = ""
    )
  
  if (!file.exists(override_file)) {
    dir.create(dirname(override_file), recursive = TRUE, showWarnings = FALSE)
    readr::write_csv(template, override_file)
    message("Created manual-override template: ", override_file)
  }
  
  ov <- read_csv_flexible(override_file)
  assert_columns(ov, c("Gene"), "HSP manual override file")
  
  for (nm in c("Family_override", "Subfamily_override", "Short_name_override",
               "Include_override", "Curator_note")) {
    if (!nm %in% colnames(ov)) ov[[nm]] <- ""
  }
  
  ov <- ov |>
    dplyr::mutate(
      dplyr::across(c(Family_override, Subfamily_override,
                      Short_name_override, Include_override, Curator_note),
                    ~ trimws(as.character(.x))),
      Include_override_parsed = parse_logical_override(Include_override)
    ) |>
    dplyr::select(Gene, Family_override, Subfamily_override,
                  Short_name_override, Include_override_parsed, Curator_note)
  
  tbl |>
    dplyr::left_join(ov, by = "Gene") |>
    dplyr::mutate(
      HSP_family = dplyr::if_else(
        !is.na(Family_override) & nzchar(Family_override),
        Family_override, HSP_family_auto
      ),
      HSP_subfamily = dplyr::if_else(
        !is.na(Subfamily_override) & nzchar(Subfamily_override),
        Subfamily_override, HSP_subfamily_auto
      ),
      Short_name = dplyr::if_else(
        !is.na(Short_name_override) & nzchar(Short_name_override),
        Short_name_override, paste0(HSP_subfamily, " | ", Gene)
      ),
      Include = dplyr::if_else(
        !is.na(Include_override_parsed),
        Include_override_parsed, Include_auto
      ),
      Curator_note = dplyr::coalesce(Curator_note, "")
    )
}

##############################
## 4) Figure device functions
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
  width_in  <- width_mm / 25.4
  height_in <- height_mm / 25.4
  dir.create(dirname(file_base), recursive = TRUE, showWarnings = FALSE)
  
  render_device <- function(label, open_device) {
    device_before <- grDevices::dev.cur()
    tryCatch({
      open_device()
      draw_heatmap_object(ht)
      grDevices::dev.off()
      invisible(TRUE)
    }, error = function(e) {
      if (grDevices::dev.cur() != device_before) {
        try(grDevices::dev.off(), silent = TRUE)
      }
      warning(label, " export failed for ", file_base, ": ",
              conditionMessage(e))
      invisible(FALSE)
    })
  }
  
  if (HSP_WRITE_PDF) {
    render_device("PDF", function() {
      grDevices::pdf(paste0(file_base, ".pdf"), width = width_in,
                     height = height_in, useDingbats = FALSE, onefile = TRUE)
    })
  }
  
  if (HSP_WRITE_SVG) {
    render_device("SVG", function() {
      svglite::svglite(paste0(file_base, ".svg"), width = width_in,
                       height = height_in, bg = "white")
    })
  }
  
  if (HSP_WRITE_PNG) {
    render_device("PNG", function() {
      grDevices::png(paste0(file_base, ".png"), width = width_in,
                     height = height_in, units = "in", res = HSP_RASTER_DPI,
                     type = "cairo-png", bg = "white")
    })
  }
  
  if (HSP_WRITE_TIFF) {
    render_device("TIFF", function() {
      grDevices::tiff(paste0(file_base, ".tiff"), width = width_in,
                      height = height_in, units = "in", res = HSP_RASTER_DPI,
                      compression = "lzw", type = "cairo", bg = "white")
    })
  }
}

get_row_order <- function(mat, cluster_rows = TRUE) {
  if (!cluster_rows || nrow(mat) <= 2L) return(seq_len(nrow(mat)))
  stats::hclust(stats::dist(mat), method = "complete")$order
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
        "HS-2" = unname(temp_colors["T2"]),
        "HS-7" = unname(temp_colors["T1"])
      )
    ),
    annotation_name_gp = grid::gpar(fontsize = 8),
    simple_anno_size = grid::unit(3.5, "mm"),
    show_legend = TRUE
  )
}

make_heatmap_pair <- function(zmat, status_mat, labels, title,
                              row_split = NULL, cluster_rows = TRUE) {
  stopifnot(identical(rownames(zmat), rownames(status_mat)))
  stopifnot(all(rownames(zmat) %in% names(labels)))
  
  zplot <- pmax(pmin(zmat, HSP_Z_LIMIT), -HSP_Z_LIMIT)
  
  if (is.null(row_split)) {
    ord <- get_row_order(zplot, cluster_rows)
    zplot <- zplot[ord, , drop = FALSE]
    status_mat <- status_mat[ord, , drop = FALSE]
    labels <- labels[rownames(zplot)]
  }
  
  z_col_fun <- circlize::colorRamp2(
    c(-HSP_Z_LIMIT, 0, HSP_Z_LIMIT),
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
    cluster_rows = if (is.null(row_split)) FALSE else cluster_rows,
    cluster_row_slices = FALSE,
    row_split = row_split,
    row_title_gp = grid::gpar(fontsize = 9, fontface = "bold"),
    row_title_rot = 0,
    row_labels = labels,
    show_row_names = TRUE,
    row_names_gp = grid::gpar(fontsize = 7.2),
    row_names_max_width = grid::unit(75, "mm"),
    column_names_rot = 45,
    column_names_gp = grid::gpar(fontsize = 8),
    column_names_centered = TRUE,
    rect_gp = grid::gpar(col = "white", lwd = 0.35),
    border = TRUE,
    heatmap_legend_param = list(
      title = "Expression",
      at = c(-HSP_Z_LIMIT, 0, HSP_Z_LIMIT),
      labels = c(paste0("<= -", HSP_Z_LIMIT), "0",
                 paste0(">= ", HSP_Z_LIMIT))
    )
  )
  
  hm_de <- ComplexHeatmap::Heatmap(
    status_mat,
    name = "DE status",
    col = c("Down" = unname(deg_colors["Down"]),
            "NS" = unname(deg_colors["NS"]),
            "Up" = unname(deg_colors["Up"])),
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
## 5) Check inputs and output directories
##############################
# Permit a convenient fallback when the files are beside the R scripts.
HSP_ANNOTATION_FILE <- first_existing(c(
  HSP_ANNOTATION_FILE,
  file.path(getwd(), "HSP_all_genes.csv")
))
HSP_DE_FILE <- first_existing(c(
  HSP_DE_FILE,
  file.path(getwd(), "HSP_DE_per_contrast.csv")
))

if (is.na(HSP_ANNOTATION_FILE)) {
  stop("HSP_all_genes.csv was not found. Expected under: ", HSP_INPUT_DIR,
       " or the working directory.")
}
if (is.na(HSP_DE_FILE)) {
  stop("HSP_DE_per_contrast.csv was not found. Expected under: ", HSP_INPUT_DIR,
       " or the working directory.")
}

hsp_out <- file.path(
  base_out,
  paste0("Figure_9A_HSP_families_", SAMPLING_TIME_CHOSEN, "_", timestamp)
)
dir_tables   <- file.path(hsp_out, "00_tables")
dir_family   <- file.path(hsp_out, "01_family_heatmaps")
dir_combined <- file.path(hsp_out, "02_combined_figure")
dir_qc       <- file.path(hsp_out, "03_QC")
for (d in c(hsp_out, dir_tables, dir_family, dir_combined, dir_qc)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

log_file <- file.path(hsp_out, "Figure_9A_HSP_pipeline.log")
log_line <- function(...) {
  txt <- paste0(..., collapse = "")
  cat(txt, "\n")
  cat(txt, "\n", file = log_file, append = TRUE)
}
if (file.exists(log_file)) file.remove(log_file)

log_line("START Figure 9A HSP pipeline")
log_line("HSP annotation: ", HSP_ANNOTATION_FILE)
log_line("HSP DE table: ", HSP_DE_FILE)
log_line("Output: ", hsp_out)

##############################
## 6) Rebuild the complete TMM logCPM matrix
##    using the same logic as script 01
##############################
stopifnot(file.exists(counts_file), file.exists(meta_file))
count_mat <- read_featurecounts_matrix(counts_file)

meta <- utils::read.csv(meta_file, sep = ";", stringsAsFactors = FALSE) |>
  dplyr::mutate(
    Sample        = sprintf("%04d", as.integer(as.character(Sample))),
    Genotype      = trimws(Genotype),
    Sampling_time = trimws(Sampling_time),
    Block         = trimws(Block),
    Temperature   = normalize_temperature(Temperature),
    Treatment     = normalize_treatment(Treatment)
  ) |>
  dplyr::mutate(
    Genotype      = factor(Genotype, levels = c("ROBILA", "PRESTO")),
    Sampling_time = factor(Sampling_time),
    Block         = factor(Block),
    Temperature   = factor(Temperature, levels = c("T0", "T1", "T2")),
    Treatment     = factor(Treatment, levels = c("EAU", "SDP"))
  ) |>
  dplyr::filter(Sampling_time == SAMPLING_TIME_CHOSEN) |>
  droplevels()

if (nrow(meta) == 0L) {
  stop("No samples remain for Sampling_time = ", SAMPLING_TIME_CHOSEN)
}

remaining_treatments <- levels(droplevels(meta$Treatment))
if (length(remaining_treatments) > 1L) {
  if (is.null(HSP_TREATMENT_CHOSEN) || !nzchar(HSP_TREATMENT_CHOSEN)) {
    stop(
      "More than one Treatment is present after P1 filtering: ",
      paste(remaining_treatments, collapse = ", "), ".\n",
      "Set HSP_TREATMENT_CHOSEN (for example, 'EAU') to avoid silently pooling treatments."
    )
  }
  meta <- meta |>
    dplyr::filter(as.character(Treatment) == HSP_TREATMENT_CHOSEN) |>
    droplevels()
  if (nrow(meta) == 0L) stop("No samples remain for Treatment = ", HSP_TREATMENT_CHOSEN)
}

common <- intersect(meta$Sample, colnames(count_mat))
if (length(common) == 0L) stop("No sample overlap between metadata and count matrix.")
if (length(common) < nrow(meta)) {
  warning(nrow(meta) - length(common), " metadata samples are absent from the count matrix.")
}

meta <- meta |>
  dplyr::filter(Sample %in% common) |>
  dplyr::arrange(match(Sample, colnames(count_mat)))
count_px <- count_mat[, meta$Sample, drop = FALSE]
stopifnot(identical(colnames(count_px), meta$Sample))

meta$GenoTemp <- factor(
  paste(meta$Genotype, meta$Temperature, sep = "_"),
  levels = HEATMAP_GROUPS_ORDER
)

if (any(table(meta$GenoTemp) == 0L)) {
  stop("At least one required Genotype x Temperature group has zero samples: ",
       paste(names(table(meta$GenoTemp))[table(meta$GenoTemp) == 0L],
             collapse = ", "))
}

use_block <- nlevels(meta$Block) > 1L
if (use_block) {
  design <- stats::model.matrix(~ 0 + GenoTemp + Block, data = meta)
} else {
  design <- stats::model.matrix(~ 0 + GenoTemp, data = meta)
}
rownames(design) <- meta$Sample

if (qr(design)$rank < ncol(design)) {
  stop("Design matrix is not full rank: ", qr(design)$rank, "/", ncol(design))
}

y <- edgeR::DGEList(counts = round(count_px), samples = meta)
keep <- edgeR::filterByExpr(y, design = design)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- edgeR::calcNormFactors(y, method = "TMM")
logcpm <- edgeR::cpm(y, log = TRUE, prior.count = 2)

cat("Samples retained:", nrow(meta), "\n")
print(table(meta$Genotype, meta$Temperature))
cat("Genes retained by filterByExpr:", nrow(logcpm), "/", nrow(count_px), "\n")

# Group means in the biological order defined in the configuration.
mean_logcpm <- vapply(
  levels(meta$GenoTemp),
  FUN.VALUE = numeric(nrow(logcpm)),
  FUN = function(grp) {
    idx <- which(meta$GenoTemp == grp)
    rowMeans(logcpm[, idx, drop = FALSE], na.rm = TRUE)
  }
)
rownames(mean_logcpm) <- rownames(logcpm)

# Explicit biological order requested for the final figure:
# ROBILA: NH -> HS-7 -> HS-2, then PRESTO: NH -> HS-7 -> HS-2.
plot_group_order <- c(
  "ROBILA_T0", "ROBILA_T1", "ROBILA_T2",
  "PRESTO_T0", "PRESTO_T1", "PRESTO_T2"
)
mean_logcpm <- mean_logcpm[, plot_group_order, drop = FALSE]
colnames(mean_logcpm) <- c(
  "ROBILA NH", "ROBILA HS-7", "ROBILA HS-2",
  "PRESTO NH", "PRESTO HS-7", "PRESTO HS-2"
)
zscore_mat <- zscore_rows_safe(mean_logcpm)

##############################
## 7) Read DE calls and annotations
##############################
hsp_de  <- read_csv_flexible(HSP_DE_FILE)
hsp_ann <- read_csv_flexible(HSP_ANNOTATION_FILE)

assert_columns(hsp_de, c("Gene", "Contrast", "logFC", "FDR"), "HSP DE table")
assert_columns(hsp_ann, c("my_gene_id", "gene_ID", "func_clean"),
               "HSP annotation table")

hsp_de <- hsp_de |>
  dplyr::mutate(
    Gene = trimws(as.character(Gene)),
    Contrast = trimws(as.character(Contrast)),
    logFC = as.numeric(logFC),
    FDR = as.numeric(FDR)
  )

if (!"Status" %in% colnames(hsp_de) || !USE_EXISTING_DE_STATUS) {
  hsp_de <- hsp_de |>
    dplyr::mutate(
      Status = dplyr::case_when(
        is.finite(FDR) & FDR <= HSP_FDR_CUTOFF & logFC >= HSP_LFC_CUTOFF ~ "Up",
        is.finite(FDR) & FDR <= HSP_FDR_CUTOFF & logFC <= -HSP_LFC_CUTOFF ~ "Down",
        TRUE ~ "NS"
      )
    )
} else {
  hsp_de <- hsp_de |>
    dplyr::mutate(
      Status = stringr::str_to_title(trimws(as.character(Status))),
      Status = dplyr::if_else(Status %in% c("Up", "Down"), Status, "NS")
    )
}

contrast_order <- c(
  "ROBILA_T1_vs_T0", "ROBILA_T2_vs_T0",
  "PRESTO_T1_vs_T0", "PRESTO_T2_vs_T0"
)
contrast_labels <- c(
  "ROBILA HS-7", "ROBILA HS-2",
  "PRESTO HS-7", "PRESTO HS-2"
)

missing_contrasts <- setdiff(contrast_order, unique(hsp_de$Contrast))
if (length(missing_contrasts) > 0L) {
  stop("Missing expected contrast(s) in HSP DE table: ",
       paste(missing_contrasts, collapse = ", "))
}

selected_gene_stats <- hsp_de |>
  # IMPORTANT: selection is restricted to the four heat-vs-NH contrasts
  # displayed/tested in this figure. This prevents a gene from entering the
  # figure only because it is DE in an unrelated contrast present in the input.
  dplyr::filter(Contrast %in% contrast_order) |>
  dplyr::group_by(Gene) |>
  dplyr::summarise(
    n_DE_contrasts = sum(Status %in% c("Up", "Down"), na.rm = TRUE),
    min_FDR = suppressWarnings(min(FDR[Status %in% c("Up", "Down")], na.rm = TRUE)),
    max_abs_logFC = suppressWarnings(max(abs(logFC[Status %in% c("Up", "Down")]), na.rm = TRUE)),
    response_direction = dplyr::case_when(
      any(Status == "Up", na.rm = TRUE) & any(Status == "Down", na.rm = TRUE) ~ "Mixed",
      any(Status == "Up", na.rm = TRUE) ~ "Up only",
      any(Status == "Down", na.rm = TRUE) ~ "Down only",
      TRUE ~ "NS"
    ),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    min_FDR = dplyr::if_else(is.infinite(min_FDR), NA_real_, min_FDR),
    max_abs_logFC = dplyr::if_else(is.infinite(max_abs_logFC), NA_real_, max_abs_logFC)
  ) |>
  dplyr::filter(n_DE_contrasts >= 1L)

selected_genes <- selected_gene_stats$Gene
cat("Genes DE in >=1 contrast:", length(selected_genes), "\n")

hsp_ann2 <- hsp_ann |>
  dplyr::transmute(
    Gene = trimws(as.character(my_gene_id)),
    gene_ID = trimws(as.character(gene_ID)),
    func_clean = as.character(func_clean)
  ) |>
  dplyr::distinct(Gene, .keep_all = TRUE)

missing_ann <- setdiff(selected_genes, hsp_ann2$Gene)
if (length(missing_ann) > 0L) {
  warning(length(missing_ann), " selected DE genes lack annotation. They will be exported as unclassified.")
}

curation <- selected_gene_stats |>
  dplyr::left_join(hsp_ann2, by = "Gene") |>
  dplyr::mutate(
    gene_ID = dplyr::coalesce(gene_ID, NA_character_),
    func_clean = dplyr::coalesce(func_clean, "")
  ) |>
  classify_hsp_annotations() |>
  apply_manual_overrides(HSP_OVERRIDE_FILE)

##############################
## 8) Status matrix
##############################
status_df <- hsp_de |>
  dplyr::filter(Gene %in% selected_genes, Contrast %in% contrast_order) |>
  dplyr::select(Gene, Contrast, Status) |>
  dplyr::distinct() |>
  tidyr::complete(Gene = selected_genes,
                  Contrast = contrast_order,
                  fill = list(Status = "NS")) |>
  tidyr::pivot_wider(names_from = Contrast, values_from = Status,
                     values_fill = "NS") |>
  dplyr::select(Gene, dplyr::all_of(contrast_order))

status_mat <- as.matrix(status_df[, contrast_order, drop = FALSE])
rownames(status_mat) <- status_df$Gene
colnames(status_mat) <- contrast_labels
status_mat[!status_mat %in% c("Up", "Down", "NS")] <- "NS"

##############################
## 9) Match selected genes to expression matrix
##############################
# The featureCounts matrix may use either my_gene_id (DC_Chr_...) or the
# corresponding DH13M14 gene_ID. Select the mapping that recovers the most
# selected HSP genes, without guessing from a single identifier.
direct_matches <- sum(selected_genes %in% rownames(mean_logcpm))
geneid_matches <- sum(hsp_ann2$gene_ID[hsp_ann2$Gene %in% selected_genes] %in%
                        rownames(mean_logcpm), na.rm = TRUE)

if (geneid_matches > direct_matches) {
  geneid_to_gene <- stats::setNames(hsp_ann2$Gene, hsp_ann2$gene_ID)
  old_ids <- rownames(mean_logcpm)
  mapped_ids <- unname(geneid_to_gene[old_ids])
  new_ids <- ifelse(!is.na(mapped_ids) & nzchar(mapped_ids), mapped_ids, old_ids)
  
  if (anyDuplicated(new_ids)) {
    stop("Mapping featureCounts gene_ID to my_gene_id creates duplicated row names. ",
         "Resolve duplicate gene mappings before plotting.")
  }
  rownames(mean_logcpm) <- new_ids
  rownames(zscore_mat)  <- new_ids
  rownames(logcpm)      <- new_ids
  cat("Expression identifier mapping: gene_ID -> my_gene_id
")
} else {
  cat("Expression identifier mapping: direct my_gene_id match
")
}

curation$In_expression_matrix <- curation$Gene %in% rownames(mean_logcpm)

missing_expr <- curation |>
  dplyr::filter(Include, !In_expression_matrix)
if (nrow(missing_expr) > 0L) {
  warning(nrow(missing_expr),
          " included genes are absent after filterByExpr and cannot be plotted. See QC table.")
}

plot_tbl <- curation |>
  dplyr::filter(Include, In_expression_matrix) |>
  dplyr::mutate(
    HSP_family = factor(
      HSP_family,
      levels = c(
        "HSP20 / sHSP", "HSP40 / DnaJ", "HSP60 / Cpn60",
        "HSP10 / Cpn10", "CCT / TCP-1", "HSP70 / HSC70",
        "HSP90", "HSP100 / ClpB", "Co-chaperones"
      )
    )
  ) |>
  dplyr::filter(!is.na(HSP_family)) |>
  dplyr::arrange(HSP_family, HSP_subfamily, Gene)

if (nrow(plot_tbl) == 0L) {
  stop("No genes remain after DE selection, HSP curation and expression matching.")
}

plot_genes <- plot_tbl$Gene
mean_plot <- mean_logcpm[plot_genes, , drop = FALSE]
z_plot    <- zscore_mat[plot_genes, , drop = FALSE]
status_plot <- status_mat[plot_genes, , drop = FALSE]
labels_all <- stats::setNames(plot_tbl$Short_name, plot_tbl$Gene)

##############################
## 10) Export tables and QC
##############################
readr::write_csv(
  curation |>
    dplyr::select(
      Gene, gene_ID, func_clean,
      n_DE_contrasts, min_FDR, max_abs_logFC, response_direction,
      HSP_family_auto, HSP_subfamily_auto, Confidence_auto,
      Evidence_auto, Include_auto, Exclusion_reason,
      HSP_family, HSP_subfamily, Short_name, Include, Curator_note,
      In_expression_matrix
    ),
  file.path(dir_tables, "HSP_family_curated.csv")
)

readr::write_csv(
  curation |>
    dplyr::filter(!Include | !is.na(Exclusion_reason) |
                    HSP_family %in% c("Excluded / ambiguous HSP-like", "Unclassified")) |>
    dplyr::select(Gene, gene_ID, func_clean, HSP_family_auto,
                  HSP_subfamily_auto, Confidence_auto,
                  Exclusion_reason, Evidence_auto, Curator_note),
  file.path(dir_tables, "HSP_excluded_or_ambiguous_candidates.csv")
)

readr::write_csv(
  curation |>
    dplyr::filter(Include, !In_expression_matrix) |>
    dplyr::select(Gene, gene_ID, HSP_family, HSP_subfamily,
                  n_DE_contrasts, func_clean),
  file.path(dir_qc, "HSP_DE_missing_after_filterByExpr.csv")
)

readr::write_csv(
  selected_gene_stats |>
    dplyr::left_join(
      hsp_de |>
        dplyr::filter(Contrast %in% contrast_order) |>
        dplyr::select(Gene, Contrast, logFC, FDR, Status) |>
        tidyr::pivot_wider(
          names_from = Contrast,
          values_from = c(logFC, FDR, Status),
          names_glue = "{Contrast}_{.value}"
        ),
      by = "Gene"
    ),
  file.path(dir_tables, "HSP_DE_selected.csv")
)

readr::write_csv(
  tibble::rownames_to_column(as.data.frame(mean_plot), "Gene") |>
    dplyr::left_join(plot_tbl |>
                       dplyr::select(Gene, HSP_family, HSP_subfamily, Short_name),
                     by = "Gene"),
  file.path(dir_tables, "HSP_meanLogCPM_six_groups.csv")
)

readr::write_csv(
  tibble::rownames_to_column(as.data.frame(z_plot), "Gene") |>
    dplyr::left_join(plot_tbl |>
                       dplyr::select(Gene, HSP_family, HSP_subfamily, Short_name),
                     by = "Gene"),
  file.path(dir_tables, "HSP_zscore_six_groups.csv")
)

family_summary <- plot_tbl |>
  dplyr::count(HSP_family, HSP_subfamily, name = "n_genes") |>
  dplyr::arrange(HSP_family, dplyr::desc(n_genes), HSP_subfamily)
readr::write_csv(family_summary,
                 file.path(dir_tables, "HSP_family_summary.csv"))

status_summary <- hsp_de |>
  dplyr::filter(Gene %in% plot_genes, Contrast %in% contrast_order,
                Status %in% c("Up", "Down")) |>
  dplyr::left_join(plot_tbl |>
                     dplyr::select(Gene, HSP_family, HSP_subfamily),
                   by = "Gene") |>
  dplyr::count(HSP_family, Contrast, Status, name = "n_genes") |>
  dplyr::arrange(HSP_family, match(Contrast, contrast_order), Status)
readr::write_csv(status_summary,
                 file.path(dir_tables, "HSP_DE_counts_by_family_and_contrast.csv"))

literature_rules <- tibble::tribble(
  ~Family, ~Operational_evidence, ~Curation_warning, ~Reference,
  "HSP20 / sHSP", "Hsp20/alpha-crystallin domain plus explicit small-HSP annotation", "An ACD alone does not prove canonical sHSP function", "Waters et al. 2008; Bondino et al. 2012",
  "HSP40 / DnaJ", "J-domain; A/B/C-like subclass inferred from zinc-finger and C-terminal domains", "Exact subclass/localisation may require sequence-domain confirmation", "Rajan & D'Silva 2017; Zhang et al. 2018",
  "HSP60 / Cpn60", "GroEL/Cpn60/60-kDa chaperonin signatures", "Separate from type-II CCT/TCP-1", "Hill & Hemmingsen 2001",
  "HSP10 / Cpn10", "Cpn10/10-kDa cochaperonin signature", "GroES-like enzyme folds are excluded", "Hill & Hemmingsen 2001",
  "HSP70 / HSC70", "HSP70/DnaK name and family signatures", "Generic sugar-kinase/HSP70/actin NBD is excluded; organellar localisation is not inferred from DnaK alone", "Sung et al. 2001",
  "HSP90", "HSP90-specific name/signature/HATPase annotation", "Generic HATPase domains in receptors and kinases are excluded; localisation is left unresolved unless supported", "Krishna & Gloor 2001",
  "HSP100 / ClpB", "Explicit ClpB/HSP101 evidence is preferred; generic ClpA/B signatures are lower confidence", "ClpA/B-like ATPases without explicit ClpB/HSP101 evidence require validation", "Myung et al. 2006",
  "Co-chaperones", "Explicit Aha1, HOP/STI1 or HSP70 nucleotide-exchange-factor annotation", "Displayed separately from canonical HSP families", "Fellerer et al. 2011 and family-specific literature"
)
readr::write_csv(literature_rules,
                 file.path(dir_tables, "HSP_classification_literature_rules.csv"))

##############################
## 11) One figure per family
##############################
families_present <- levels(droplevels(plot_tbl$HSP_family))

for (fam in families_present) {
  fam_tbl <- plot_tbl |>
    dplyr::filter(as.character(HSP_family) == fam) |>
    dplyr::arrange(HSP_subfamily, Gene)
  
  genes <- fam_tbl$Gene
  z_fam <- zscore_mat[genes, , drop = FALSE]
  st_fam <- status_mat[genes, , drop = FALSE]
  lab_fam <- stats::setNames(fam_tbl$Short_name, fam_tbl$Gene)
  
  ht_fam <- make_heatmap_pair(
    zmat = z_fam,
    status_mat = st_fam,
    labels = lab_fam,
    title = paste0(fam, " — DE in at least one heat-stress contrast"),
    row_split = NULL,
    cluster_rows = HSP_CLUSTER_ROWS
  )
  
  n <- length(genes)
  height_mm <- max(95, 42 + HSP_ROW_MM * n)
  file_base <- file.path(dir_family,
                         paste0("Figure_9A_", safe_file_stub(fam), "_heatmap"))
  save_complex_heatmap(ht_fam, file_base,
                       width_mm = HSP_FIG_WIDTH_MM,
                       height_mm = height_mm)
  cat("Saved family:", fam, "—", n, "genes\n")
}

##############################
## 12) Combined Figure 9A
##############################
combined_tbl <- plot_tbl
if (!INCLUDE_COCHAPERONES_IN_COMBINED) {
  combined_tbl <- combined_tbl |>
    dplyr::filter(as.character(HSP_family) != "Co-chaperones")
}

combined_tbl <- combined_tbl |>
  dplyr::arrange(HSP_family, HSP_subfamily, Gene)
combined_genes <- combined_tbl$Gene
combined_split <- droplevels(combined_tbl$HSP_family)

ht_combined <- make_heatmap_pair(
  zmat = zscore_mat[combined_genes, , drop = FALSE],
  status_mat = status_mat[combined_genes, , drop = FALSE],
  labels = stats::setNames(combined_tbl$Short_name, combined_tbl$Gene),
  title = "Figure 9A | Heat-responsive HSP and co-chaperone families",
  row_split = combined_split,
  cluster_rows = HSP_CLUSTER_ROWS
)

combined_height_mm <- max(180, 55 + HSP_ROW_MM * nrow(combined_tbl))
save_complex_heatmap(
  ht_combined,
  file.path(dir_combined, "Figure_9A_HSP_families_combined"),
  width_mm = HSP_FIG_WIDTH_MM + 10,
  height_mm = combined_height_mm
)

##############################
## 13) Reproducibility record
##############################
utils::write.csv(meta,
                 file.path(dir_qc, "metadata_used_for_HSP_heatmaps.csv"),
                 row.names = FALSE)
utils::write.csv(as.data.frame(design),
                 file.path(dir_qc, "design_matrix_used.csv"),
                 row.names = TRUE)
saveRDS(
  list(
    metadata = meta,
    design = design,
    dge = y,
    logCPM = logcpm,
    mean_logCPM = mean_logcpm,
    zscore = zscore_mat,
    curation = curation,
    plot_table = plot_tbl,
    status_matrix = status_mat
  ),
  file.path(dir_qc, "Figure_9A_HSP_analysis_objects.rds")
)

capture.output(utils::sessionInfo(),
               file = file.path(hsp_out, "sessionInfo.txt"))

cat("\nDONE Figure 9A HSP pipeline\n")
cat("Included and plotted genes:", nrow(plot_tbl), "\n")
cat("Families plotted:", paste(families_present, collapse = "; "), "\n")
cat("Combined figure formats: PDF, SVG, PNG, TIFF\n")
cat("Output directory:", hsp_out, "\n")
