library(bigsnpr)
ukb <- snp_attach("data/UKBB_imp_HM3.rds")
G <- ukb$genotypes

set.seed(1)
ind.val <- sort(sample(nrow(G), 10e3))
ind.gwas <- sort(sample(setdiff(rows_along(G), ind.val), 300e3))
ind.test <- setdiff(rows_along(G), c(ind.gwas, ind.val))

library(tidyverse)
pred_ldpred1 <-
  tibble(res_file = list.files("results/simu-ldpred1", full.names = TRUE),
       betas = map(res_file, ~ as.list(setNames(as_tibble(readRDS(.x)),
                                                c(rep("grid", 7), "inf"))))) %>%
  unnest_longer("betas", indices_to = "method") %>%
  print()

str(all_betas <- do.call("cbind", pred_ldpred1$betas))

(NCORES <- as.integer(Sys.getenv("SLURM_JOB_CPUS_PER_NODE")) - 1L)
bigparallelr::set_blas_ncores(NCORES)
pred <- big_prodMat(G, all_betas)

pred_ldpred1$pred <- as.list(as.data.frame(pred))
pred_ldpred1$pheno <- map(pred_ldpred1$res_file, ~ {
  readRDS(sub("results/simu-ldpred1", "data/pheno-simu", .x, fixed = TRUE))
})
res_ldpred1 <- pred_ldpred1 %>%
  mutate(auc_val = map2_dbl(pred, pheno, ~ AUC(-.x[ind.val], .y[ind.val]))) %>%
  group_by(res_file, method) %>%
  arrange(desc(auc_val)) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(basename = gsubfn::gsubfn("(.+)_([0-9]+)\\.rds$", "\\1 \\2", basename(res_file))) %>%
  separate("basename", c("params", "iter"), sep = " ", convert = TRUE) %>%
  transmute(params, iter,
            Method = ifelse(method == "grid", "best_ldpred1_grid", "ldpred1_inf"),
            AUC = map2_dbl(pred, pheno, ~ AUC(-.x[ind.test], .y[ind.test]))) %>%
  print()

saveRDS(res_ldpred1, "results/all-simu-ldpred1.rds")
