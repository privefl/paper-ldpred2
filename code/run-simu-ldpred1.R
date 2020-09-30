# system("git clone  https://github.com/bvilhjal/ldpred.git --branch 1.0.0")

library(bigsnpr)
ukb <- snp_attach("data/UKBB_imp_HM3.rds")
G <- ukb$genotypes

load("data/ind_gwas_val_test.RData")

library(dplyr)
files <- tibble(
  basename = list.files("data/sumstats-simu"),
  chr = as.integer(sub(".*_chr([0-9]+).*\\.rds$", "\\1", basename)),
  N = sub(".*_chr[0-9]+(.*)\\.rds$", "\\1", basename),
  pheno_file = file.path("data/pheno-simu", sub("_chr[0-9]+.*(\\.rds)$", "\\1", basename)),
  gwas_file = file.path("data/sumstats-simu", basename),
  res_file = file.path("results/simu-ldpred1", sub("_chr[0-9]+", "", basename))
) %>%
  mutate(N = ifelse(N == "", 300000L, as.integer(sub("_", "", N)))) %>%
  group_by(pheno_file, N, res_file) %>%
  arrange(chr) %>%
  summarize(gwas_files = list(c(gwas_file))) %>%
  print()
# files <- tibble(basename = list.files("data/pheno-simu"),
#                 res_file = file.path("results/simu-ldpred1", basename))
dim(files)
table(files$N)
bigassertr::assert_dir("results/simu-ldpred1")

files_sub <- files %>%
  filter(
    # sub("(.*_)([0-9]+)(\\.rds)$", "\\2", basename) %in% 1:3, #iter
    !file.exists(res_file)) %>%
  print()


library(future.batchtools)
plan(batchtools_slurm(resources = list(
  t = "3-00:00", c = "4", mem = "32g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

furrr::future_pmap(files_sub, function(pheno_file, N, res_file, gwas_files) {

  # num <- 1
  # pheno_file <- files_sub$pheno_file[num]
  # N <- files_sub$N[num]
  # res_file <- files_sub$res_file[num]
  # gwas_files <- files_sub$gwas_files[[num]]

  y <- readRDS(pheno_file)
  gwas <- do.call("rbind", lapply(gwas_files, readRDS))

  tmp <- tempfile(tmpdir = "tmp-data")
  file_sumstats <- paste0(tmp, ".txt")
  cbind(ukb$map, beta = gwas$estim, pval = predict(gwas, log10 = FALSE)) %>%
    na.omit() %>%
    bigreadr::fwrite2(file_sumstats, sep = "\t") %>%
    readLines(n = 5)
  (Neff <- round(4 / (1 / sum(y[ind.gwas] == 0) + 1 / sum(y[ind.gwas] == 1)) * N / 300e3))

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
