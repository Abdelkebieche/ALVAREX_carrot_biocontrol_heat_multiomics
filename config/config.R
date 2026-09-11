# Central configuration for the publication repository.
#
# No institutional/HPC path is stored in this file. Override paths with
# environment variables or copy this file locally and edit the copy.

CARROT_CONFIG_LOADED <- TRUE

find_repo_root <- function(start = getwd()) {
  d <- normalizePath(start, mustWork = FALSE)
  repeat {
    if (file.exists(file.path(d, 'README.md')) &&
        file.exists(file.path(d, 'config', 'config.R'))) return(d)
    parent <- dirname(d)
    if (identical(parent, d)) break
    d <- parent
  }
  normalizePath(start, mustWork = FALSE)
}

repo_root_default <- find_repo_root()
REPO_ROOT <- Sys.getenv('CARROT_REPO_ROOT', unset = repo_root_default)
DATA_ROOT <- Sys.getenv('CARROT_DATA_ROOT', unset = file.path(REPO_ROOT, 'data'))
OUTPUT_ROOT <- Sys.getenv('CARROT_OUTPUT_ROOT', unset = file.path(REPO_ROOT, 'outputs'))

# RNA-seq
RNA_COUNTS_FILE <- Sys.getenv(
  'CARROT_RNA_COUNTS',
  unset = file.path(DATA_ROOT, 'transcriptomics', 'gene_counts_featureCounts.txt')
)
RNA_METADATA_FILE <- Sys.getenv(
  'CARROT_RNA_METADATA',
  unset = file.path(REPO_ROOT, 'metadata', 'rnaseq_metadata.csv')
)
RNA_WGCNA_ROOT <- Sys.getenv(
  'CARROT_RNA_WGCNA_ROOT',
  unset = file.path(OUTPUT_ROOT, 'transcriptomics', 'wgcna')
)
MASIGPRO_ROOT <- Sys.getenv(
  'CARROT_MASIGPRO_ROOT',
  unset = file.path(OUTPUT_ROOT, 'transcriptomics', 'masigpro')
)
BIOLOGICAL_HEATMAP_ROOT <- Sys.getenv(
  'CARROT_BIOLOGICAL_HEATMAP_ROOT',
  unset = file.path(OUTPUT_ROOT, 'transcriptomics', 'biological_heatmaps')
)

# Annotation resources
GENE_ID_MAP_FILE <- Sys.getenv(
  'CARROT_GENE_ID_MAP',
  unset = file.path(DATA_ROOT, 'annotation', 'GeneID_mapping_my_to_PlanT2T.csv')
)
DH13_ANNOT_FILE <- Sys.getenv(
  'CARROT_DH13_ANNOT',
  unset = file.path(DATA_ROOT, 'annotation', 'Daucus_carota_DH13M14.txt')
)
HSP_ANNOT_FILE <- Sys.getenv(
  'CARROT_HSP_ANNOT',
  unset = file.path(DATA_ROOT, 'annotation', 'HSP_all_genes.csv')
)
FUNCTIONAL_INPUT_DIR <- Sys.getenv(
  'CARROT_FUNCTIONAL_INPUT_DIR',
  unset = file.path(DATA_ROOT, 'annotation', 'functional_families')
)
ANNOTATION_BUILDER_OUTPUT <- Sys.getenv(
  'CARROT_ANNOTATION_OUTPUT',
  unset = file.path(OUTPUT_ROOT, 'annotation_builders')
)

ORGDB_URL <- 'https://biobigdata.nju.edu.cn/plant2t/orgdb/org.Daucus.carota.DH13M14.eg.db_1.0.tar.gz'
ORGDB_PKG <- 'org.Daucus.carota.DH13M14.eg.db'

# LC-HRMS
LCMS_FILE <- Sys.getenv(
  'CARROT_LCMS_FEATURES',
  unset = file.path(DATA_ROOT, 'metabolomics', 'lcms_features.csv')
)
LCMS_METADATA_FILE <- Sys.getenv(
  'CARROT_LCMS_METADATA',
  unset = file.path(DATA_ROOT, 'metabolomics', 'metadata.csv')
)
LCMS_WGCNA_ROOT <- Sys.getenv(
  'CARROT_LCMS_WGCNA_ROOT',
  unset = file.path(OUTPUT_ROOT, 'metabolomics', 'wgcna')
)
LCMS_FIGURE_ROOT <- Sys.getenv(
  'CARROT_LCMS_FIGURE_ROOT',
  unset = file.path(OUTPUT_ROOT, 'metabolomics', 'figures')
)
LCMS_REFERENCE_RUN <- Sys.getenv('CARROT_LCMS_REFERENCE_RUN', unset = '')

# Phenotype
PHENOTYPE_FILE <- Sys.getenv(
  'CARROT_PHENOTYPE_FILE',
  unset = file.path(DATA_ROOT, 'phenotype', 'Data_phenotype.xlsx')
)
PHENOTYPE_OUTPUT <- Sys.getenv(
  'CARROT_PHENOTYPE_OUTPUT',
  unset = file.path(OUTPUT_ROOT, 'phenotype')
)

# Publication analysis constants
RNA_SOFT_POWER <- c(PRESTO = 12, ROBILA = 11)
LCMS_SOFT_POWER <- 12
KME_SCREEN <- 0.70
KME_HIGH_CONF <- 0.85
FDR_CUTOFF <- 0.05

for (d in c(OUTPUT_ROOT, RNA_WGCNA_ROOT, MASIGPRO_ROOT,
            BIOLOGICAL_HEATMAP_ROOT, LCMS_WGCNA_ROOT,
            LCMS_FIGURE_ROOT, PHENOTYPE_OUTPUT, ANNOTATION_BUILDER_OUTPUT)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
