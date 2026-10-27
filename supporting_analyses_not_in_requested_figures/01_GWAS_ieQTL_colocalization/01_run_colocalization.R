#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
})

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(7, 8)) {
  stop(paste(
    "Usage: gwas_ieqtl_coloc.R <gwas_base_dir> <ieqtl_file>",
    "<tissue> <cell_type> <n_ieqtl> <out_dir> <min_shared_snps> [trait]"
  ))
}

gwas_base_dir <- args[[1]]
ieqtl_file <- args[[2]]
tissue <- args[[3]]
cell_type <- args[[4]]
n_ieqtl <- as.integer(args[[5]])
out_dir <- args[[6]]
min_shared_snps <- as.integer(args[[7]])
trait_filter <- if (length(args) == 8) args[[8]] else NULL

required_ieqtl <- c(
  "pheno_id", "variant_id", "af", "beta_g2", "beta_se_g2", "pval_g2"
)

ieqtl <- fread(ieqtl_file)
missing_ieqtl <- setdiff(required_ieqtl, names(ieqtl))
if (length(missing_ieqtl) > 0) {
  stop("Missing ieQTL columns: ", paste(missing_ieqtl, collapse = ", "))
}
if (!is.finite(n_ieqtl) || n_ieqtl <= 0) stop("Invalid ieQTL sample size")
if (!is.finite(min_shared_snps) || min_shared_snps < 2) stop('min_shared_snps must be >= 2')

ieqtl[, chrom := sub("_.*$", "", variant_id)]
ieqtl <- ieqtl[
  is.finite(af) & af > 0 & af < 1 &
    is.finite(beta_g2) & is.finite(beta_se_g2) & beta_se_g2 > 0 &
    is.finite(pval_g2) & pval_g2 > 0 & pval_g2 <= 1
]
setkey(ieqtl, variant_id)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
trait_dirs <- list.dirs(gwas_base_dir, recursive = FALSE, full.names = TRUE)
if (!is.null(trait_filter)) {
  trait_dirs <- trait_dirs[basename(trait_dirs) == trait_filter]
}
if (!length(trait_dirs)) stop('No matching trait directories found')
# Preflight ALL outputs before starting: never erase a previous analysis.
planned <- file.path(out_dir, paste0('GWAS_', basename(trait_dirs), '_ieQTL.pph4'))
if (any(file.exists(planned))) stop('Output exists; select a fresh output directory')
error_file <- file.path(out_dir, 'coloc_errors.tsv')
if (file.exists(error_file)) stop('Existing error log; select a fresh output directory')
writeLines('trait\tchrom\tgene\terror', error_file)

for (trait_dir in trait_dirs) {
  trait <- basename(trait_dir)
  outfile <- file.path(out_dir, paste0("GWAS_", trait, "_ieQTL.pph4"))
  wrote_header <- FALSE

  for (chr_num in 1:26) {
    gwas_file <- file.path(trait_dir, paste0("chr", chr_num, ".tsv"))
    if (!file.exists(gwas_file) || file.info(gwas_file)$size <= 30) next

    gwas <- fread(gwas_file)
    required_gwas <- c("SNP", "Freq", "b", "se", "p")
    missing_gwas <- setdiff(required_gwas, names(gwas))
    if (length(missing_gwas) > 0) {
      stop("Missing GWAS columns in ", gwas_file, ": ",
           paste(missing_gwas, collapse = ", "))
    }

    gwas <- unique(gwas[
      is.finite(Freq) & Freq > 0 & Freq < 1 &
        is.finite(b) & is.finite(se) & se > 0 & is.finite(p) & p > 0 & p <= 1
    ], by = "SNP")
    if (nrow(gwas) == 0) next

    gwas[, maf := pmin(Freq, 1 - Freq)]
    gwas[, n_est := 1 / (se^2 * 2 * maf * (1 - maf))]
    n_gwas <- as.integer(round(median(gwas$n_est[is.finite(gwas$n_est)], na.rm = TRUE)))
    if (!is.finite(n_gwas) || n_gwas <= 0) next

    iq <- ieqtl[chrom == as.character(chr_num)]
    if (nrow(iq) == 0) next
    merged <- merge(gwas, iq, by.x = "SNP", by.y = "variant_id")
    if (nrow(merged) == 0) next

    for (gene in unique(merged$pheno_id)) {
      x <- merged[pheno_id == gene]
      x <- unique(x, by = "SNP")
      if (nrow(x) < min_shared_snps) next

      d1 <- list(
        beta = x$b,
        varbeta = x$se^2,
        snp = x$SNP,
        type = "quant",
        N = rep(n_gwas, nrow(x)),
        MAF = x$maf
      )
      d2 <- list(
        beta = x$beta_g2,
        varbeta = x$beta_se_g2^2,
        snp = x$SNP,
        type = "quant",
        N = rep(n_ieqtl, nrow(x)),
        MAF = pmin(x$af, 1 - x$af)
      )

      result <- tryCatch(coloc.abf(d1, d2), error = function(e) {
        fwrite(data.table(trait = trait, chrom = chr_num, gene = gene,
                          error = conditionMessage(e)), error_file, sep = '\t',
               append = TRUE, col.names = FALSE)
        warning('coloc failed for ', trait, '/', chr_num, '/', gene, ': ', conditionMessage(e))
        NULL
      })
      if (is.null(result)) next

      out <- as.data.table(as.list(result$summary))
      out[, `:=`(
        SNP = x$SNP[which.min(x$p)],
        phenotype = gene,
        chrom = paste0("chr", chr_num),
        tissue = tissue,
        cell_type = cell_type,
        trait = trait,
        n_gwas = n_gwas,
        n_ieqtl = n_ieqtl
      )]
      fwrite(out, outfile, sep = "\t", append = wrote_header,
             col.names = !wrote_header)
      wrote_header <- TRUE
    }
  }
}
