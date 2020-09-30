## reference data
# download.file("https://www.dropbox.com/s/p9aqanhxvxaqv8k/ldblk_1kg_eur.tar.gz?dl=1",
#               destfile = "data/ldblk_1kg_eur.tar.gz")
# untar("data/ldblk_1kg_eur.tar.gz")

## software
# system("git clone https://github.com/getian107/PRScs.git")
# system("python --version")  # 3.7.4

library(dplyr)
files <- tidyr::expand_grid(
  basename = list.files("data/sumstats"),
  chr = 1:22,
  phi = c(10^-(6:0), NA)
) %>%
  mutate(
    res_file = sprintf("results/prscs/%s_pst_eff_a1_b0.5_phi%.e_chr%d.txt",
                       sub("\\.rds$", "", basename), phi, chr),
    gwas_file = file.path("data/sumstats", basename)) %>%
  print()

bigassertr::assert_dir("results/prscs")

files_sub <- files %>%
  arrange(basename, chr) %>%
  filter(
    # grepl("Asthma", basename), #trait
    !file.exists(res_file)) %>%
  select(-basename) %>%
  print()


rsid <- bigreadr::fread2("data/UKBB_imp_HM3_val.bim")[[2]]

library(future.batchtools)
NCORES <- 8
plan(batchtools_slurm(resources = list(
  t = "4:00:00", c = NCORES, mem = "64g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))


furrr::future_pmap(files_sub, function(chr, phi, res_file, gwas_file) {

  tmp <- tempfile(tmpdir = "tmp-data", fileext = ".txt")
  library(dplyr)
  sumstats <- readRDS(gwas_file)
  sumstats %>%
    mutate(SNP = rsid[`_NUM_ID_`]) %>%
    select(SNP, A1 = a1, A2 = a0, BETA = beta, P = p) %>%
    bigreadr::fwrite2(tmp, sep = "\t") %>%
    readLines(n = 5)

  prscs <- "PRScs/PRScs.py"
  library(glue)
  # system(glue("{prscs} --help"))

  system(glue(
    "OMP_NUM_THREADS=", NCORES,
    " python {prscs}",
    " --ref_dir=ldblk_1kg_eur",
    " --bim_prefix=data/UKBB_imp_HM3_val",
    " --sst_file={tmp}",
    " --n_gwas={round(median(sumstats$n_eff))}",
    if (is.na(phi)) "" else " --phi={phi}",
    " --chrom={chr}",
    " --out_dir={sub('_pst.+$', '', res_file)}"
  ))

  file.remove(tmp)
})
