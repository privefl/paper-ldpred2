library(bigsnpr)
ukb <- snp_attach("data/UKBB_imp_HM3.rds")
G <- ukb$genotypes

load("data/ind_val_test.RData")

library(dplyr)
res_sbayesr <-
  tibble(res_file = list.files("results/sbayesr", "\\.snpRes$", full.names = TRUE)) %>%
  mutate(Trait = sub(".+/([[:alnum:]]+)_chr[0-9]+\\.snpRes$", "\\1", res_file),
         chr = as.integer(sub(".+/[[:alnum:]]+_chr([0-9]+)\\.snpRes", "\\1", res_file))) %>%
  group_by(Trait) %>%
  arrange(chr) %>%
  summarize(beta = list(bigreadr::fread2(res_file))) %>%
  print()

betas <- matrix(0, ncol(G), nrow(res_sbayesr))
for (k in cols_along(betas)) {
  tmp <- res_sbayesr$beta[[k]]
  ind <- match(tmp$Name, ukb$map$marker.ID)
  stopifnot(!anyNA(ind))
  stopifnot(all.equal(ukb$map$allele1[ind], tmp$A2))
  stopifnot(all.equal(ukb$map$allele2[ind], tmp$A1))
  betas[ind, k] <- tmp$A1Effect
}

(NCORES <- as.integer(Sys.getenv("SLURM_JOB_CPUS_PER_NODE")) - 1L)
# bigparallelr::set_blas_ncores(NCORES)
pred <- big_parallelize(G, function(X, ind, betas, ind.test) {
  bigstatsr::big_prodMat(X, betas[ind, ], ind.row = ind.test, ind.col = ind)
}, p.combine = plus, ncores = NCORES, betas = betas, ind.test = ind.test)
res_sbayesr$pred <- as.list(as.data.frame(pred))
res_sbayesr$beta <- NULL

AUCBoot_no_NA <- function(pred, target) {
  is_not_na <- which(!is.na(target))
  AUCBoot(pred[is_not_na], target[is_not_na])
}

res_sbayesr$AUCBoot <- purrr::pmap(res_sbayesr, function(Trait, pred) {
  y <- readRDS(paste0("data/pheno/", Trait, ".rds"))
  AUCBoot_no_NA(pred, y[ind.test])
})
res_sbayesr$pred <- NULL
res_sbayesr$Method <- "SBayesR"

saveRDS(res_sbayesr, "results/all-sbayesr.rds")
