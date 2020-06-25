library(bigsnpr)
ukb <- snp_attach("data/UKBB_imp_HM3.rds")
G <- ukb$genotypes

set.seed(1)
ind.val <- sort(sample(nrow(G), 10e3))
ind.test <- setdiff(rows_along(G), ind.val)

library(tidyverse)
pred_ldpred1 <-
  tibble(res_file = list.files("results/ldpred1", full.names = TRUE),
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
  readRDS(sub("results/ldpred1", "data/pheno", .x, fixed = TRUE))
})

AUC_no_NA <- function(pred, target) {
  is_not_na <- which(!is.na(target))
  AUC(pred[is_not_na], target[is_not_na])
}
AUCBoot_no_NA <- function(pred, target) {
  is_not_na <- which(!is.na(target))
  AUCBoot(pred[is_not_na], target[is_not_na])
}

library(future)
plan(multisession, workers = NCORES)
res_ldpred1 <- pred_ldpred1 %>%
  mutate(auc_val = map2_dbl(pred, pheno, ~ AUC_no_NA(-.x[ind.val], .y[ind.val]))) %>%
  group_by(res_file, method) %>%
  arrange(desc(auc_val)) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(Trait = gsubfn::strapply(
    res_file, "results/ldpred1/(.+)\\.rds$", simplify = 'c')) %>%
  transmute(Trait,
            Method = ifelse(method == "grid", "best_ldpred1_grid", "ldpred1_inf"),
            AUCBoot = furrr::future_map2(pred, pheno, ~ {
              AUCBoot_no_NA(-.x[ind.test], .y[ind.test])
            })) %>%
  print()

saveRDS(res_ldpred1, "results/all-ldpred1.rds")
