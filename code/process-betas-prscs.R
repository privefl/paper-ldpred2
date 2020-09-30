library(bigsnpr)
ukb <- snp_attach("data/UKBB_imp_HM3.rds")
G <- ukb$genotypes

load("data/ind_val_test.RData")

library(dplyr)
beta_prscs <-
  tibble(res_file = list.files("results/prscs", "\\.txt$", full.names = TRUE)) %>%
  mutate(Trait = sub("(.+/)([[:alnum:]]+)(.+_chr[0-9]+\\.txt)$", "\\2", res_file),
         chr = as.integer(sub("(.+_chr)([0-9]+)(\\.txt)$", "\\2", res_file)),
         phi = sub("(.+_phi)(.+)(_chr[0-9]+\\.txt)$", "\\2", res_file)) %>%
  group_by(Trait, phi) %>%
  arrange(chr) %>%
  summarize(beta = list(bigreadr::fread2(res_file))) %>%
  print()

betas <- matrix(0, ncol(G), nrow(beta_prscs))
for (k in cols_along(betas)) {
  tmp <- beta_prscs$beta[[k]]
  betas[match(tmp$V2, ukb$map$rsid), k] <- tmp$V6
}

(NCORES <- as.integer(Sys.getenv("SLURM_JOB_CPUS_PER_NODE")) - 1L)
# bigparallelr::set_blas_ncores(NCORES)
# pred <- big_prodMat(G, betas)
pred <- big_parallelize(G, function(X, ind, betas) {
  bigstatsr::big_prodMat(X, betas[ind, ], ind.col = ind)
}, p.combine = plus, ncores = NCORES, betas = betas)

beta_prscs$pred <- as.list(as.data.frame(pred))
beta_prscs$pheno <-
  lapply(paste0("data/pheno/", beta_prscs$Trait, ".rds"), readRDS)

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

res_prscs <- beta_prscs %>%
  filter(phi != "auto") %>%
  mutate(auc_val = furrr::future_map2_dbl(
    pred, pheno, ~ AUC_no_NA(-.x[ind.val], .y[ind.val]))) %>%
  group_by(Trait) %>%
  arrange(desc(auc_val)) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(Trait, Method = "PRS-CS", AUCBoot = furrr::future_map2(
    pred, pheno, ~ AUCBoot_no_NA(-.x[ind.test], .y[ind.test]))) %>%
  print()

res_prscs_auto <- beta_prscs %>%
  filter(phi == "auto") %>%
  transmute(Trait, Method = "PRS-CS-auto", AUCBoot = furrr::future_map2(
    pred, pheno, ~ AUCBoot_no_NA(-.x[ind.test], .y[ind.test]))) %>%
  print()

saveRDS(dplyr::bind_rows(res_prscs, res_prscs_auto), "results/all-prscs.rds")
