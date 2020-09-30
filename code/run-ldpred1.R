library(bigsnpr)
ukb <- snp_attach("data/UKBB_imp_HM3.rds")
G <- ukb$genotypes
ukb$map$chromosome <- as.integer(ukb$map$chromosome)

load("data/ind_val_test.RData")

library(dplyr)
files <- tibble(basename = list.files("data/sumstats"),
                res_file = file.path("results/ldpred1", basename))
bigassertr::assert_dir("results/ldpred1")

files_sub <- files %>%
  filter(!file.exists(res_file)) %>%
  print()

library(future.batchtools)
plan(batchtools_slurm(resources = list(
  t = "60:00:00", c = "4", mem = "32g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

furrr::future_pmap(files_sub, function(basename, res_file) {

  y <- readRDS(file.path("data/pheno", basename))
  sumstats <- readRDS(file.path("data/sumstats", basename))

  tmp <- tempfile(tmpdir = "tmp-data")
  file_sumstats <- paste0(tmp, ".txt")
  cbind(ukb$map[sumstats$`_NUM_ID_`, ], beta = sumstats$beta, pval = sumstats$p) %>%
    bigreadr::fwrite2(file_sumstats, sep = "\t")
  # readLines(file_sumstats, n = 5)
  (Neff <- round(median(sumstats$n_eff)))

  file_hdf5 <- paste0(tmp, ".hdf5")
  ldpred <- "ldpred/LDpred.py"
  # system(glue::glue("python3 {ldpred} coord --help"))
  system(glue::glue(
    "python3 {ldpred} coord",
    " --gf data/UKBB_imp_HM3_val",
    " --ssf {file_sumstats}",
    " --skip-coordination",
    " --rs marker.ID --A1 allele1 --A2 allele2 --pos physical.pos --chr chromosome",
    " --pval pval --eff beta --beta",
    " --N {Neff}",
    " --out {file_hdf5}"
  ))

  # system(glue::glue("python3 {ldpred} gibbs --help"))
  system(glue::glue(
    "python3 {ldpred} gibbs",
    " --cf {file_hdf5}",
    " --ldr {round(nrow(ukb$map) / 3000)}",
    " --N {Neff}",
    " --ldf {tmp}",
    " --out {tmp}"
  ))

  ext <- c(sprintf("_LDpred_p%.4e.txt", c(1, 0.3, 0.1, 0.03, 0.01, 0.003, 0.001)),
           "_LDpred-inf.txt")
  files_ldpred <- paste0(tmp, ext)
  beta_ldpred <- sapply(files_ldpred, function(file) {
    res_ldpred <- bigreadr::fread2(file, select = c(3, 7))
    beta_ldpred <- rep(0, nrow(ukb$map))
    beta_ldpred[match(res_ldpred$sid, ukb$map$marker.ID)] <- res_ldpred[[2]]
    beta_ldpred
  })

  unlink(paste0(tmp, "*"))
  saveRDS(beta_ldpred, res_file)
})
