library(bigsnpr)
ukb <- snp_attach("data/UKBB_imp_HM3.rds")
G <- ukb$genotypes

set.seed(1)
ind.val <- sort(sample(nrow(G), 10e3))
ind.test <- setdiff(rows_along(G), ind.val)


library(dplyr)
files <- tibble(basename = list.files("data/sumstats"),
                res_file = file.path("results/lassosum", basename))
bigassertr::assert_dir("results/lassosum")

files_sub <- files %>%
  filter(!file.exists(res_file)) %>%
  print()

AUC_no_NA <- function(pred, target) {
  is_not_na <- which(!is.na(target))
  AUC(pred[is_not_na], target[is_not_na])
}

rsid <- bigreadr::fread2("data/UKBB_imp_HM3_val.bim")[[2]]

library(future.batchtools)
NCORES <- 16
plan(batchtools_slurm(resources = list(
  t = "4:00:00", c = NCORES, mem = "128g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

furrr::future_pmap(files_sub, function(basename, res_file) {

  y <- readRDS(file.path("data/pheno", basename))
  sumstats <- readRDS(file.path("data/sumstats", basename))

  library(lassosum)
  doParallel::registerDoParallel(cl <- parallel::makeCluster(NCORES))
  system.time(
    out <- lassosum.pipeline(
      cor = with(sumstats, beta / beta_se / sqrt(n_eff)),
      snp = rsid[sumstats$`_NUM_ID_`],
      A1 = sumstats$a1,
      A2 = sumstats$a0,
      exclude.ambiguous = FALSE,
      test.bfile = "data/UKBB_imp_HM3_val",
      LDblocks = "EUR.hg19",
      cluster = cl,
      destandardize = TRUE
    )
  ) # 25 min
  parallel::stopCluster(cl)

  auc.val <- apply(-do.call("cbind", out$pgs), 2, AUC_no_NA, target = y[ind.val])
  beta_lassosum <- do.call("cbind", out$beta)[, which.max(auc.val)]

  ind <- which(beta_lassosum != 0)
  pred <- big_prodVec(G, beta_lassosum[ind], ind.row = ind.test,
                      ind.col = sumstats$`_NUM_ID_`[ind])

  saveRDS(-pred, res_file)
})
