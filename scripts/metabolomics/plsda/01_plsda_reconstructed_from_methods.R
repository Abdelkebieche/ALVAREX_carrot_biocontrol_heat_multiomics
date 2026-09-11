# Reconstructed from the manuscript method specification.
# This was NOT present in the supplied RStudio export. Validate against the
# original analysis script before using for a final archival release.

if (!exists('CARROT_CONFIG_LOADED')) source('../../../config/config.R')
source('../wgcna/00_helpers.R')

required <- c('dplyr','mixOmics','ggplot2')
for (p in required) if (!requireNamespace(p, quietly=TRUE)) stop('Missing package: ', p)

lcms <- read_lcms_matrix(LCMS_FILE)
X_all <- lcms$mat
meta <- read.csv(LCMS_METADATA_FILE, sep=';', stringsAsFactors=FALSE)
meta$Sample <- sprintf('%04d', as.integer(as.character(meta$Sample)))
meta$Temperature <- normalize_temperature(meta$Temperature)
meta$Treatment <- normalize_treatment(meta$Treatment)
meta$Sampling_time <- trimws(meta$Sampling_time)

out_dir <- file.path(LCMS_FIGURE_ROOT, 'PLSDA')
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

run_one <- function(genotype, stage) {
  m <- meta[meta$Genotype == genotype & meta$Sampling_time == stage, , drop=FALSE]
  samples <- intersect(m$Sample, colnames(X_all))
  m <- m[match(samples, m$Sample), , drop=FALSE]
  X <- t(X_all[, samples, drop=FALSE])
  X <- log2(X + 1)
  keep <- apply(X, 2, var, na.rm=TRUE) > 0
  X <- X[, keep, drop=FALSE]
  Y <- if (stage == 'P1') factor(m$Temperature) else interaction(m$Temperature, m$Treatment, drop=TRUE)
  fit <- mixOmics::plsda(X, Y, ncomp=2, scale=TRUE)
  perf <- mixOmics::perf(fit, validation='Mfold', folds=4, nrepeat=50, progressBar=FALSE)
  saveRDS(fit, file.path(out_dir, paste0(genotype,'_',stage,'_plsda.rds')))
  saveRDS(perf, file.path(out_dir, paste0(genotype,'_',stage,'_cv.rds')))
  grDevices::pdf(file.path(out_dir, paste0(genotype,'_',stage,'_PLSDA.pdf')), width=7, height=6)
  mixOmics::plotIndiv(fit, comp=c(1,2), group=Y, ind.names=FALSE,
                      legend=TRUE, title=paste(genotype, stage))
  grDevices::dev.off()
}

for (g in c('PRESTO','ROBILA')) for (s in c('P1','P3','P4')) run_one(g,s)
