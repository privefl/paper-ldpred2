library(bigsnpr)
ukb <- snp_attach("data/UKBB_imp_HM3.rds")
G <- ukb$genotypes
CHR <- as.integer(ukb$map$chromosome)

load("data/ind_gwas_val_test.RData")

library(dplyr)
files <- tibble(
  basename = list.files("data/sumstats-simu", pattern = "all_40_3000_15"),
  chr = as.integer(sub(".*_chr([0-9]+).*\\.rds$", "\\1", basename)),
  N = sub(".*_chr[0-9]+(.*)\\.rds$", "\\1", basename),
  pheno_file = file.path("data/pheno-simu", sub("_chr[0-9]+.*(\\.rds)$", "\\1", basename)),
  gwas_file = file.path("data/sumstats-simu", basename)
) %>%
  mutate(N = ifelse(N == "", 300000L, as.integer(sub("_", "", N))),
         res_file = file.path("results/simu-ldpred2-gwide",
                              paste0(sub("_chr[0-9]+.*", "", basename), "_", N, ".rds"))) %>%
  group_by(pheno_file, N, res_file) %>%
  arrange(chr) %>%
  summarize(gwas_files = list(c(gwas_file))) %>%
  print()
table(files$N)
bigassertr::assert_dir("results/simu-ldpred2-gwide")

files_sub <- files %>%
  filter(!file.exists(res_file)) %>%
  # group_by(N) %>% slice(1:3) %>% ungroup() %>%
  print()

library(future.batchtools)
NCORES <- 30
plan(batchtools_slurm(resources = list(
  t = "2-00:00", c = NCORES + 2, mem = "250g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

furrr::future_pmap(files_sub, function(pheno_file, gwas_files, N, res_file) {

  # (pheno_file <- files_sub$pheno_file[7])
  # (gwas_files <- files_sub$gwas_files[[7]])
  # (N <- files_sub$N[7])
  # (res_file <- files_sub$res_file[7])

  y <- readRDS(pheno_file)
  (Neff <- 4 / (1 / sum(y[ind.gwas] == 0) + 1 / sum(y[ind.gwas] == 1)) * N / 300e3)
  sumstats <- do.call("rbind", lapply(gwas_files, readRDS)) %>%
    transmute(beta = estim, beta_se = std.err, n_eff = Neff,
              "_NUM_ID_" = row_number(), chr = CHR) %>%
    na.omit()

  ind.val2 <- ind.val[!is.na(y[ind.val])]

  tmp <- tempfile(tmpdir = "tmp-data")
  on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)

  for (chr in 1:22) {

    # print(chr)

    ## indices in 'sumstats'
    ind.chr <- which(sumstats$chr == chr)
    ## indices in 'G'
    ind.chr2 <- sumstats$`_NUM_ID_`[ind.chr]
    ## indices in 'corr'
    ind.chr3 <- match(ind.chr2, which(CHR == chr))

    corr0 <- readRDS(paste0("data/corr/chr", chr, ".rds"))[ind.chr3, ind.chr3]

    if (chr == 1) {
      ld <- Matrix::colSums(corr0^2)
      corr <- as_SFBM(corr0, tmp)
    } else {
      ld <- c(ld, Matrix::colSums(corr0^2))
      corr$add_columns(corr0, nrow(corr))
    }
  }

  df_beta <- sumstats
  (ldsc <- with(df_beta, snp_ldsc(ld, length(ld), chi2 = (beta / beta_se)^2,
                                  sample_size = n_eff, blocks = NULL)))
  h2_est <- ldsc[["h2"]]

  # LDpred2-inf
  beta_inf <- snp_ldpred2_inf(corr, df_beta, h2_est)

  # LDpred2-grid
  (h2_seq <- round(h2_est * c(0.7, 1, 1.4), 4))
  (p_seq <- signif(seq_log(1e-5, 1, length.out = 21), 2))
  (params <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(FALSE, TRUE)))

  beta_grid <- snp_ldpred2_grid(corr, df_beta, params, ncores = NCORES)
  params$sparsity <- colMeans(beta_grid == 0)

  bigparallelr::set_blas_ncores(NCORES)
  pred_grid <- big_prodMat(G, beta_grid, ind.row = ind.val2,
                           ind.col = df_beta[["_NUM_ID_"]])
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
  pred_auto <- big_prodMat(G, beta_auto, ind.row = ind.val,
                           ind.col = df_beta[["_NUM_ID_"]], ncores = NCORES)
  sc <- apply(pred_auto, 2, sd)
  keep <- abs(sc - median(sc)) < 3 * mad(sc)
  final_beta_auto <- rowMeans(beta_auto[, keep])

  # compute predictions for test set
  betas <- cbind(beta_inf, best_beta_grid_nosp, best_beta_grid_sp,
                 beta_auto = final_beta_auto)
  pred_test <- big_prodMat(G, betas, ind.row = ind.test,
                           ind.col = df_beta[["_NUM_ID_"]], ncores = NCORES)

  # save results
  res <- list(pred = setNames(as.data.frame(pred_test), colnames(betas)),
              params = params, auto = multi_auto[keep], h2_ldsc = h2_est)
  saveRDS(res, res_file)
})
