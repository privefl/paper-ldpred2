# zip <- runonce::download_file(
#   "https://cnsgenomics.com/software/gctb/download/gctb_2.02_Linux.zip",
#   dir = "tmp-data")
# unzip(zip, exdir = "tmp-data")
gctb <- "tmp-data/gctb_2.02_Linux/gctb"
# bigsnpr:::make_executable(gctb)
plink <- bigsnpr::download_plink("tmp-data")

## Compute LD

library(future.batchtools)
plan(batchtools_slurm(resources = list(
  t = "2-00:00", c = 2, mem = "125g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))
future.apply::future_lapply(1:22, function(chr) {
  if (!file.exists(paste0("tmp-data/ldm_", chr, ".ldm.shrunk.bin"))) {
    system(glue::glue(
      "{plink} --bfile data/UKBB_imp_HM3_val",
      " --chr {chr} --maf 0.01 --make-bed",
      " --out tmp-data/myplink_{chr}"
    ))
    system.time(
      system(glue::glue(
        "{gctb} --bfile tmp-data/myplink_{chr}",
        " --make-shrunk-ldm",
        # " --gen-map tmp-data/chr{chr}.OMNI.interpolated_genetic_map", # too few variants matched
        " --out tmp-data/ldm_{chr}"
      ))
    )
    # 19 min for chr 22 (16410 SNPs)
    # 1h30 for chr 14 (37031 SNPs)
    unlink(paste0("tmp-data/myplink_", chr, "*"))
  }
}, future.globals = c("gctb", "plink"))


## Compute SBayesR

library(dplyr)
files <- tidyr::expand_grid(
  basename = list.files("data/sumstats"),
  chr = 1:22
) %>%
  mutate(
    res_file = paste0("results/sbayesr/", sub("\\.rds$", "", basename), "_chr", chr),
    gwas_file = file.path("data/sumstats", basename)) %>%
  print()

bigassertr::assert_dir("results/sbayesr")

files_sub <- files %>%
  arrange(basename, chr) %>%
  filter(!file.exists(paste0(res_file, ".snpRes"))) %>%
  select(-basename) %>%
  print()

obj.bed <- bigsnpr::bed("data/UKBB_imp_HM3_val.bed")
af <- bigsnpr::bed_MAF(obj.bed, ncores = 15)$af

rsid <- bigreadr::fread2("data/UKBB_imp_HM3_val.bim")[[2]]

library(future.batchtools)
plan(batchtools_slurm(resources = list(
  t = "12:00:00", c = 4, mem = "64g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

furrr::future_pmap(files_sub, function(chr, res_file, gwas_file) {

  tmp <- tempfile(tmpdir = "tmp-data", fileext = ".ma")
  library(dplyr)
  sumstats <- readRDS(gwas_file)
  sumstats %>%
    mutate(SNP = rsid[`_NUM_ID_`], freq = af[`_NUM_ID_`], N = round(n_eff)) %>%
    select(SNP, A1 = a1, A2 = a0, freq, b = beta, se = beta_se, p, N) %>%
    bigreadr::fwrite2(tmp, sep = " ") %>%
    readLines(n = 5) %>%
    writeLines()

  system(glue::glue(
    gctb,
    " --ldm tmp-data/ldm_{chr}.ldm.shrunk",
    " --sbayes R --pi 0.95,0.02,0.02,0.01 --gamma 0.0,0.01,0.1,1",
    " --gwas-summary {tmp}",
    " --chain-length 10000 --burn-in 2000",
    " --out {res_file} --out-freq 100"
  ))

  file.remove(tmp)
})
