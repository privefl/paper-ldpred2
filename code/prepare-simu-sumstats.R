library(bigsnpr)
ukb <- snp_attach("data/UKBB_imp_HM3.rds")
G <- ukb$genotypes
CHR <- as.integer(ukb$map$chromosome)
POS <- ukb$map$physical.pos

bigassertr::assert_dir("data/sumstats-simu")

set.seed(1)
ind.val <- sort(sample(nrow(G), 10e3))
ind.gwas <- sort(sample(setdiff(rows_along(G), ind.val), 300e3))

pheno_files <- list.files("data/pheno-simu")
all_pheno <- lapply(pheno_files, function(file) {
  readRDS(file.path("data/pheno-simu", file))
})

library(future.batchtools)
plan(batchtools_slurm(resources = list(
  t = "12:00:00", c = "16", mem = "128g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

future.apply::future_lapply(1:22, function(chr) {

  ind.chr <- which(CHR == chr)

  purrr::walk2(pheno_files, all_pheno, ~ {

    gwas_file <- file.path("data/sumstats-simu",
                           sub("\\.rds$", paste0("_chr", chr, ".rds"), .x))

    if (!file.exists(gwas_file)) {
      print(gwas_file)
      gwas <- big_univLogReg(G, .y[ind.gwas], ind.train = ind.gwas,
                             ind.col = ind.chr, ncores = 16)
      saveRDS(gwas, gwas_file)
    }
  })
})
