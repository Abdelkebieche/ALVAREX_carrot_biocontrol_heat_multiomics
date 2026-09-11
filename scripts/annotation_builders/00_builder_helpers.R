# Minimal shared helpers required by the historical annotation-builder scripts.

load_pkgs <- function(pkgs, install_missing = FALSE, bioc = FALSE) {
  if (bioc && !requireNamespace('BiocManager', quietly = TRUE)) {
    if (install_missing) install.packages('BiocManager') else stop('Missing BiocManager')
  }
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      if (!install_missing) stop('Missing package: ', p)
      if (bioc) BiocManager::install(p, ask = FALSE, update = FALSE) else install.packages(p)
    }
    suppressPackageStartupMessages(library(p, character.only = TRUE))
  }
  invisible(TRUE)
}

clean_sample_names <- function(x) {
  x <- gsub('^X', '', x)
  m <- regmatches(x, regexpr('\\d{4}', x))
  ifelse(m == '' | is.na(m), x, m)
}

normalize_temperature <- function(x) {
  x_chr <- trimws(as.character(x))
  if (all(grepl('^T[0-9]+$', x_chr))) return(x_chr)
  if (all(x_chr %in% c('0','1','2'))) return(paste0('T', x_chr))
  suppressWarnings(x_num <- as.numeric(x_chr))
  if (!all(is.na(x_num)) && all(x_num %in% c(0,1,2))) return(paste0('T', x_num))
  x_chr
}

normalize_treatment <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- gsub('\\s+', '', x)
  ifelse(x %in% c('EAU','WATER','H2O','CTRL','CONTROL'), 'EAU',
         ifelse(x %in% c('SDP'), 'SDP', x))
}

read_featurecounts_matrix <- function(path) {
  fc <- utils::read.table(path, header=TRUE, sep='\t', comment.char='#', quote='',
                          check.names=FALSE, stringsAsFactors=FALSE)
  stopifnot('Geneid' %in% colnames(fc))
  drop_cols <- intersect(colnames(fc), c('Chr','Start','End','Strand','Length'))
  mat <- fc[, setdiff(colnames(fc), drop_cols), drop=FALSE]
  rownames(mat) <- mat$Geneid
  mat$Geneid <- NULL
  mat <- as.matrix(mat)
  mode(mat) <- 'numeric'
  colnames(mat) <- make.unique(clean_sample_names(colnames(mat)))
  mat
}

zscore_rows <- function(m) {
  z <- t(scale(t(as.matrix(m))))
  z[is.na(z)] <- 0
  z
}

make_status <- function(df, alpha, lfc_threshold) {
  df$FDR[is.na(df$FDR)] <- 1
  df$Status <- ifelse(df$FDR < alpha & df$logFC >= lfc_threshold, 'Up',
                      ifelse(df$FDR < alpha & df$logFC <= -lfc_threshold, 'Down', 'NS'))
  df
}
