library(bigsnpr)
ukb <- snp_attach("data/UKBB_imp_HM3.rds")
G <- ukb$genotypes
CHR <- as.integer(ukb$map$chromosome)

set.seed(1)
ind.val <- sort(sample(nrow(G), 10e3))
ind.test <- setdiff(rows_along(G), ind.val)


library(dplyr)
files <- tidyr::expand_grid(
  basename = list.files("data/sumstats"),
  chr = 1:22
) %>%
  mutate(
    res_file = file.path("results/ldpred2", stringr::str_replace(
      basename, "\\.rds$", paste0("_chr", chr, ".rds"))),
    pheno_file = file.path("data/pheno",  basename),
    gwas_file = file.path("data/sumstats", basename)) %>%
  print()
table(files$chr)
all(file.exists(files$pheno_file)) & all(file.exists(files$gwas_file))
bigassertr::assert_dir("results/ldpred2")

files_sub <- files %>%
  arrange(chr) %>%
  filter(
    # grepl("BRCA", basename), #trait
    !file.exists(res_file)) %>%
  print()

library(future.batchtools)
NCORES <- 8
plan(batchtools_slurm(resources = list(
  t = "6:00:00", c = NCORES, mem = "64g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

furrr::future_pmap(files_sub[-1], function(res_file, chr, pheno_file, gwas_file) {

  y <- readRDS(pheno_file)
  sumstats <- readRDS(gwas_file)

  ind.val2 <- ind.val[!is.na(y[ind.val])]

  ## indices in 'sumstats'
  ind.chr <- which(sumstats$chr == chr)
  (df_beta <- sumstats[ind.chr, c("beta", "beta_se", "n_eff")])
  ## indices in 'G'
  ind.chr2 <- sumstats$`_NUM_ID_`[ind.chr]
  ## indices in 'corr'
  ind.chr3 <- match(ind.chr2, which(CHR == chr))
  corr <- readRDS(paste0("data/corr/chr", chr, ".rds"))[ind.chr3, ind.chr3]

  (ldsc <- snp_ldsc2(corr, df_beta))
  h2_est <- ldsc[["h2"]]

  if (h2_est < 1e-4) {

    rep(list(rep(0, length(ind.test))), 4) %>%
      setNames(c("beta_inf", "best_beta_grid_nosp", "best_beta_grid_sp", "beta_auto")) %>%
      saveRDS(res_file)

  } else {

    tmp <- tempfile(tmpdir = "tmp-data")
    corr <- bigsparser::as_SFBM(as(corr, "dgCMatrix"), tmp)
    on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)

    # LDpred2-inf
    beta_inf <- snp_ldpred2_inf(corr, df_beta, h2_est)

    # LDpred2-grid
    (h2_seq <- round(h2_est * c(0.7, 1, 1.4), 4))
    (p_seq <- signif(seq_log(1e-5, 1, length.out = 21), 2))
    (params <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(FALSE, TRUE)))

    beta_grid <- snp_ldpred2_grid(corr, df_beta, params, ncores = NCORES)

    pred_grid <- big_prodMat(G, beta_grid, ind.row = ind.val2,
                             ind.col = ind.chr2)
    params$score <- big_univLogReg(as_FBM(pred_grid), y[ind.val2])$score

    library(dplyr)
    best_beta_grid_nosp <- params %>%
      mutate(id = row_number()) %>%
      filter(!sparse) %>%
      arrange(desc(score)) %>%
      slice(1) %>%
      pull(id) %>%
      beta_grid[, .]

    best_beta_grid_sp <- params %>%
      mutate(id = row_number()) %>%
      filter(sparse) %>%
      arrange(desc(score)) %>%
      slice(1) %>%
      pull(id) %>%
      beta_grid[, .]

    # LDpred2-auto
    multi_auto <- snp_ldpred2_auto(corr, df_beta, h2_init = h2_est,
                                   vec_p_init = seq_log(1e-4, 0.9, 30),
                                   ncores = NCORES)
    beta_auto <- sapply(multi_auto, function(auto) auto$beta_est)
    pred_auto <- big_prodMat(G, beta_auto, ind.row = ind.val, ind.col = ind.chr2)
    sc <- apply(pred_auto, 2, sd)
    final_beta_auto <- rowMeans(beta_auto[, abs(sc - median(sc)) < 3 * mad(sc)])

    # best from all previous 4 + all auto
    all_beta <- cbind(beta_inf, best_beta_grid_nosp, best_beta_grid_sp,
                      final_beta_auto, beta_auto)
    pred_val <- big_prodMat(G, all_beta, ind.row = ind.val2, ind.col = ind.chr2)
    score <- big_univLogReg(as_FBM(pred_val), y[ind.val2])$score
    best_beta <- all_beta[, which.max(score)]

    # compute predictions for test set
    pred_test <- tibble(beta_inf, best_beta_grid_nosp, best_beta_grid_sp,
                        beta_auto = final_beta_auto, best_beta) %>%
      lapply(function(beta) big_prodVec(G, beta, ind.row = ind.test,
                                        ind.col = ind.chr2))

    # save results
    saveRDS(pred_test, res_file)

  }
})
