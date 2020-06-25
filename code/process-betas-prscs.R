library(bigsnpr)
ukb <- snp_attach("data/UKBB_imp_HM3.rds")
G <- ukb$genotypes

set.seed(1)
ind.val <- sort(sample(nrow(G), 10e3))
ind.test <- setdiff(rows_along(G), ind.val)

library(tidyverse)
beta_prscs <-
  tibble(res_file = list.files("results/prscs", "\\.txt$", full.names = TRUE)) %>%
  mutate(Trait = sub("(.+/)([[:alnum:]]+)(.+_chr[0-9]+\\.txt)$", "\\2", res_file),
         chr = as.integer(sub("(.+_chr)([0-9]+)(\\.txt)$", "\\2", res_file)),
         phi = as.double(sub("(.+_phi)(.+)(_chr[0-9]+\\.txt)$", "\\2", res_file))) %>%
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
bigparallelr::set_blas_ncores(NCORES)
pred <- big_prodMat(G, betas)

beta_prscs$pred <- as.list(as.data.frame(pred))
beta_prscs$pheno <- map(beta_prscs$Trait, ~ {
  readRDS(paste0("data/pheno/", ., ".rds"))
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
res_prscs <- beta_prscs %>%
  mutate(auc_val = map2_dbl(pred, pheno, ~ AUC_no_NA(-.x[ind.val], .y[ind.val]))) %>%
  group_by(Trait) %>%
  arrange(desc(auc_val)) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(Trait,
            Method = "PRS-CS",
            AUCBoot = furrr::future_map2(pred, pheno, ~ {
              AUCBoot_no_NA(-.x[ind.test], .y[ind.test])
            })) %>%
  print()

saveRDS(res_prscs, "results/all-prscs.rds")
