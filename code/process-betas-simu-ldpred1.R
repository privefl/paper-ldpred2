library(bigsnpr)
ukb <- snp_attach("data/UKBB_imp_HM3.rds")
G <- ukb$genotypes

load("data/ind_gwas_val_test.RData")

library(dplyr)
pred_ldpred1 <-
  tibble(
    res_file = list.files("results/simu-ldpred1", full.names = TRUE),
    betas = lapply(res_file, function(.x) {
      as.list(setNames(as_tibble(readRDS(.x)), c(rep("grid", 7), "inf")))
    })
  ) %>%
  tidyr::unnest_longer("betas", indices_to = "method") %>%
  print()

str(all_betas <- do.call("cbind", pred_ldpred1$betas))

(NCORES <- as.integer(Sys.getenv("SLURM_JOB_CPUS_PER_NODE")) - 1L)
bigparallelr::set_blas_ncores(NCORES)
pred <- runonce::save_run(big_prodMat(G, all_betas),
                          file = "tmp-data/all_pred_simu_ldpred1.rds")

pred_ldpred1$pred <- as.list(as.data.frame(pred))
library(future.apply)
plan(multisession, workers = NCORES)
pheno_files <- sub("results/simu-ldpred1/(.*_15_[0-9]+).*",
                   "data/pheno-simu/\\1.rds", pred_ldpred1$res_file)
# all(file.exists(pheno_files))
pred_ldpred1$pheno <- lapply(pheno_files, readRDS)

res_ldpred1 <- pred_ldpred1 %>%
  mutate(auc_val = purrr::map2_dbl(pred, pheno, ~ AUC(-.x[ind.val], .y[ind.val]))) %>%
  group_by(res_file, method) %>%
  arrange(desc(auc_val)) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(params = sub("^(.+_15)_([0-9]+)(.*)\\.rds$", "\\1\\3", basename(res_file)),
         iter = sub("^(.+_15)_([0-9]+)(.*)\\.rds$", "\\2", basename(res_file))) %>%
  transmute(params, iter,
            Method = ifelse(method == "grid", "best_ldpred1_grid", "ldpred1_inf"),
            AUC = purrr::map2_dbl(pred, pheno, ~ AUC(-.x[ind.test], .y[ind.test]))) %>%
  print()

saveRDS(res_ldpred1, "results/all-simu-ldpred1.rds")
