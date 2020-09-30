library(bigsnpr)
library(bigreadr)

NCORES <- 30  # TO MODIFY

## Information for the variants provided in the LD reference
map_ldref <- readRDS("ld-ref/map.rds")

## Breast cancer summary statistics

# download.file("ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/MichailidouK_29059683_GCST004988/oncoarray_bcac_public_release_oct17%20(1).txt.gz",
#               destfile = "tmp-data/sumstats_BRCA.txt.gz")
# R.utils::gunzip("tmp-data/sumstats_BRCA.txt.gz")
sumstats <- fread2("tmp-data/sumstats_BRCA.txt", na.strings = "NULL",
                   select = c("chr", "position_b37", "a0", "a1",
                              "bcac_onco_icogs_gwas_beta",
                              "bcac_onco_icogs_gwas_se"),
                   col.names = c("chr", "pos", "a0", "a1", "beta", "beta_se"))
sumstats$n_eff <- 4 / (1 / 137045 + 1 / 119078)

info_snp <- snp_match(sumstats, map_ldref)
# 11,792,542 variants to be matched.
# 0 ambiguous SNPs have been removed.
# 1,054,287 variants have been matched; 0 were flipped and 0 were reversed.
(info_snp <- tidyr::drop_na(tibble::as_tibble(info_snp)))

sd_ldref <- with(info_snp, sqrt(2 * af_UKBB * (1 - af_UKBB)))
sd_ss <- with(info_snp, 2 / sqrt(n_eff * beta_se^2))

is_bad <-
  sd_ss < (0.5 * sd_ldref) | sd_ss > (sd_ldref + 0.1) | sd_ss < 0.1 | sd_ldref < 0.05

library(ggplot2)
qplot(sd_ldref, sd_ss, color = is_bad, alpha = I(0.5)) +
  theme_bigstatsr() +
  coord_equal() +
  scale_color_viridis_d(direction = -1) +
  geom_abline(linetype = 2, color = "red") +
  labs(x = "Standard deviations derived from allele frequencies of the LD reference",
       y = "Standard deviations derived from the summary statistics",
       color = "Removed?")

df_beta <- info_snp[!is_bad, ]

# Here, you also want to restrict to the variants present
# in your test data as well. For this, you can use something like
in_test <- vctrs::vec_in(df_beta[, c("chr", "pos")], map_test[, c("chr", "pos")])
df_beta <- df_beta[in_test, ]

tmp <- tempfile(tmpdir = "tmp-data")

for (chr in 1:22) {

  cat(chr, ".. ", sep = "")

  ## indices in 'df_beta'
  ind.chr <- which(df_beta$chr == chr)
  ## indices in 'map_ldref'
  ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]
  ## indices in 'corr_chr'
  ind.chr3 <- match(ind.chr2, which(map_ldref$chr == chr))

  corr_chr <- readRDS(paste0("ld-ref/LD_chr", chr, ".rds"))[ind.chr3, ind.chr3]

  if (chr == 1) {
    corr <- as_SFBM(corr_chr, tmp)
  } else {
    corr$add_columns(corr_chr, nrow(corr))
  }
}

# Heritability estimation of LD score regression
# to be used as a starting value in LDpred2-auto
(ldsc <- with(df_beta, snp_ldsc(ld, ld_size = nrow(map_ldref),
                                chi2 = (beta / beta_se)^2,
                                sample_size = n_eff,
                                ncores = NCORES)))
#        int     int_se         h2      h2_se
# 1.09045635 0.01332333 0.12732974 0.01104159
h2_est <- ldsc[["h2"]]


# LDpred2-auto
multi_auto <- snp_ldpred2_auto(corr, df_beta, h2_init = h2_est,
                               vec_p_init = seq_log(1e-4, 0.9, 30),
                               ncores = NCORES)  # 13 hours
beta_auto <- sapply(multi_auto, function(auto) auto$beta_est)

# Use your favorite tool (e.g. PLINK) to get the predictions 'pred_auto'
# corresponding to the vectors of effects in 'beta_auto', e.g. here I can do
library(bigsnpr)
ukb <- snp_attach("data/UKBB_imp_HM3.rds")
G <- ukb$genotypes
map <- dplyr::transmute(ukb$map,
                        chr = as.integer(chromosome), pos = physical.pos,
                        a0 = allele1, a1 = allele2)  # reversed somehow..
map_pgs <- df_beta[1:4]; map_pgs$beta <- 1
map_pgs2 <- snp_match(map_pgs, map)

load("data/ind_val_test.RData")
pred_auto <- big_prodMat(G, sweep(beta_auto, 1, map_pgs2$beta, '*'),
                         ind.row = ind.test,
                         ind.col = map_pgs2[["_NUM_ID_"]],
                         ncores = NCORES)

# Filter outlier predictions and average remaining ones -> your polygenic score
sc <- apply(pred_auto, 2, sd)
keep <- abs(sc - median(sc)) < 3 * mad(sc)
final_pred_auto <- rowMeans(pred_auto[, keep])

# Testing of final PGS
AUCBoot_no_NA <- function(pred, target) {
  is_not_na <- which(!is.na(target))
  AUCBoot(pred[is_not_na], target[is_not_na])
}
y <- readRDS("data/pheno/BRCA.rds")
AUCBoot_no_NA(final_pred_auto, y[ind.test])
#        Mean        2.5%       97.5%          Sd
# 0.655530387 0.650819295 0.660281045 0.002426852

# cleanup
file.remove(paste0(tmp, ".sbk"))
