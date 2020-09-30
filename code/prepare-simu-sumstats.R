library(bigsnpr)
ukb <- snp_attach("data/UKBB_imp_HM3.rds")
G <- ukb$genotypes
CHR <- as.integer(ukb$map$chromosome)
POS <- ukb$map$physical.pos

load("data/ind_gwas_val_test.RData")

#### Simus varying disease architecture ####

pheno_files <- list.files("data/pheno-simu")

library(future.batchtools)
plan(batchtools_slurm(resources = list(
  t = "12:00:00", c = "16", mem = "128g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

bigassertr::assert_dir("data/sumstats-simu")

future.apply::future_lapply(1:22, function(chr) {

  ind.chr <- which(CHR == chr)

  purrr::walk(pheno_files, ~ {

    gwas_file <- file.path("data/sumstats-simu",
                           sub("\\.rds$", paste0("_chr", chr, ".rds"), .x))

    if (!file.exists(gwas_file)) {
      print(gwas_file)
      y <- readRDS(file.path("data/pheno-simu", .x))
      gwas <- big_univLogReg(G, y[ind.gwas], ind.train = ind.gwas,
                             ind.col = ind.chr, ncores = 16)
      saveRDS(gwas, gwas_file)
    }
  })
})

#### Simus varying sample size ####

grid <- expand.grid(
  pheno_file = pheno_files[grepl("all_40_3000_15", pheno_files)],
  N = c(10e3, 20e3, 50e3, 120e3),
  chr = 1:22
)

library(future.batchtools)
plan(batchtools_slurm(resources = list(
  t = "12:00:00", c = "16", mem = "128g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

furrr::future_pmap(grid, function(pheno_file, N, chr) {

  ind.chr <- which(CHR == chr)

  gwas_file <- file.path(
    "data/sumstats-simu",
    sub("\\.rds$", paste0("_chr", chr, "_", N, ".rds"), pheno_file))

  if (!file.exists(gwas_file)) {
    print(gwas_file)
    y <- readRDS(file.path("data/pheno-simu", pheno_file))
    ind.gwas.sub <- sort(sample(ind.gwas, size = N))
    gwas <- big_univLogReg(G, y[ind.gwas.sub], ind.train = ind.gwas.sub,
                           ind.col = ind.chr, ncores = 16)
    saveRDS(gwas, gwas_file)
  }
})
