# Compatibility configuration for the functional dictionary reference builders.
# This file intentionally contains no institutional or HPC-specific paths.

if (!exists('CARROT_CONFIG_LOADED')) source('../../config/config.R')

SAMPLING_TIME_CHOSEN <- Sys.getenv('CARROT_DICTIONARY_SAMPLING_TIME', unset = 'P1')
SEED <- 123
alpha <- 0.05
PADJ_CUTOFF <- 0.05
lfc_threshold <- 0.5
INSTALL_MISSING_PKGS <- FALSE

timestamp <- format(Sys.time(), '%Y%m%d_%H%M%S')
base_out <- ANNOTATION_BUILDER_OUTPUT
counts_file <- RNA_COUNTS_FILE
meta_file <- RNA_METADATA_FILE
map_file <- GENE_ID_MAP_FILE

# Labels and palettes are retained only for reproducing historical figures.
temp_colors <- c('T0'='#2166AC', 'T1'='#B2182B', 'T2'='#F46D43')
trt_colors  <- c('EAU'='#4575B4', 'SDP'='#D73027')
geno_colors <- c('ROBILA'='#1B7837', 'PRESTO'='#762A83')
deg_colors  <- c('NS'='grey80', 'Up'='#D73027', 'Down'='#4575B4')
temp_display <- c('T0'='NH', 'T1'='HS-7', 'T2'='HS-2')
HEATMAP_GROUPS_ORDER <- c('ROBILA_T0','ROBILA_T1','ROBILA_T2',
                          'PRESTO_T0','PRESTO_T1','PRESTO_T2')

set.seed(SEED)
options(stringsAsFactors = FALSE)
