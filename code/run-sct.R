library(bigsnpr)
ukb <- snp_attach("data/UKBB_imp_HM3.rds")
G <- ukb$genotypes
CHR <- as.integer(ukb$map$chromosome)
POS <- ukb$map$physical.pos

set.seed(1)
ind.val <- sort(sample(nrow(G), 10e3))
ind.test <- setdiff(rows_along(G), ind.val)


library(dplyr)
files <- tibble(basename = list.files("data/sumstats"),
                res_file = file.path("results/sct", basename))
bigassertr::assert_dir("results/sct")

files_sub <- files %>%
  filter(!file.exists(res_file)) %>%
  print()

library(future.batchtools)
NCORES <- 16
plan(batchtools_slurm(resources = list(
  t = "12:00:00", c = NCORES, mem = "128g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

furrr::future_pmap(files_sub, function(basename, res_file) {

  # basename <- "CAD.rds"
  y <- readRDS(file.path("data/pheno", basename))
  sumstats <- readRDS(file.path("data/sumstats", basename))

  ind.val2 <- ind.val[!is.na(y[ind.val])]

  ### Clumping

  chi2 <- with(sumstats, (beta / beta_se)^2)
  lpval <- rep(NA, ncol(G))
  lpval[sumstats$`_NUM_ID_`] <-
    -pchisq(chi2, df = 1, lower.tail = FALSE, log.p = TRUE) / log(10)

  system.time(
    all_keep <- snp_grid_clumping(G, CHR, POS, lpS = lpval, ind.row = ind.val,
                                  exclude = which(is.na(lpval)), ncores = NCORES)
  ) # 3420 sec

  beta <- rep(0, ncol(G))
  beta[sumstats$`_NUM_ID_`] <- sumstats$beta

  ### Thresholding

  system.time(
    multi_PRS <- snp_grid_PRS(
      G, all_keep, betas = beta, lpS = lpval, ind.row = ind.val2,
      grid.lpS.thr = seq_log(0.1, 0.9999 * max(lpval, na.rm = TRUE), 50),
      backingfile = sprintf("tmp-data/%s_scores", sub("\\.rds$", "", basename)),
      ncores = NCORES
    )
  ) # 1965 sec

  ### Stacking

  system.time(
    final_mod <- snp_grid_stacking(multi_PRS, y[ind.val2],
                                   K = 5, ncores = NCORES)
  ) # 80 sec
  # mod <- final_mod$mod
  # plot(mod)
  # summary(mod)

  new_beta <- final_mod$beta.G

  length(ind <- which(new_beta != 0))  # 491,240
  # summary(new_beta)
  # summary(new_beta[which(sign(new_beta * beta) < 0)])

  pred_sct <- final_mod$intercept +
    big_prodVec(G, new_beta[ind], ind.row = ind.test, ind.col = ind)

  ### Best C+T

  library(tidyverse)
  grid2 <- attr(all_keep, "grid") %>%
    mutate(thr.lp = list(attr(multi_PRS, "grid.lpS.thr")), id = row_number()) %>%
    unnest(cols = "thr.lp")
  s <- nrow(grid2)
  grid2$auc <- big_apply(multi_PRS, a.FUN = function(X, ind, s, y.val) {
    # Sum over all chromosomes, for the same C+T parameters
    single_PRS <- rowSums(X[, ind + s * (0:21)])
    bigstatsr::AUC(single_PRS, y.val)
  }, ind = 1:s, s = s, y.val = y[ind.val2],
  a.combine = 'c', block.size = 1, ncores = NCORES)

  max_prs <- grid2 %>% arrange(desc(auc)) %>% slice(1:10) %>% print() %>% slice(1)
  ind.keep <- unlist(map(all_keep, max_prs$id))
  pred_clumping <- snp_PRS(G, beta[ind.keep], ind.test = ind.test, ind.keep = ind.keep,
                           lpS.keep = lpval[ind.keep], thr.list = max_prs$thr.lp)

  # Save results
  saveRDS(list(SCT = pred_sct, `C+T` = pred_clumping), res_file)
})
