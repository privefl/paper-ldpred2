#### Simulations ####

library(bigsnpr)
ukb <- snp_attach("data/UKBB_imp_HM3.rds")
G <- ukb$genotypes

load("data/ind_gwas_val_test.RData")

sd <- readRDS("data/sd.rds")

sumstats_simu <- do.call("rbind", lapply(1:22, function(chr) {
  readRDS(paste0("data/sumstats-simu/all_40_300_15_1_chr", chr, ".rds"))
}))
y <- readRDS("data/pheno-simu/all_40_300_15_1.rds")
Neff <- 4 / (1 / sum(y[ind.gwas] == 0) + 1 / sum(y[ind.gwas] == 1))

sd_val <- sd
sd_ss <- 2 / sqrt(Neff * sumstats_simu$std.err^2)

qplot(sd_val, sd_ss, alpha = I(0.5)) +
  theme_bigstatsr() +
  coord_equal() +
  geom_abline(linetype = 2, color = "red") +
  labs(x = "Standard deviations in the validation set",
       y = "Standard deviations derived from the summary statistics")

# ggsave("figures/sd-approx-simu.png", width = 7, height = 7)


#### Real data ####

library(bigsnpr)
library(bigreadr)

## Information for the variants provided in the LD reference
info <- readRDS(url("https://github.com/privefl/bigsnpr/raw/master/data-raw/hm3_variants.rds"))
str(info)

#### Coronary artery disease (CAD) ####

# download.file("ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/NikpayM_26343387_GCST003116/cad.add.160614.website.txt",
#               destfile = "tmp-data/sumstats_CAD.txt")
sumstats <- fread2("tmp-data/sumstats_CAD.txt",
                   select = c("chr", "bp_hg19", "noneffect_allele",
                              "effect_allele", "beta", "se_dgc", "p_dgc"),
                   col.names = c("chr", "pos", "a0", "a1", "beta", "beta_se", "p"))
sumstats$n_eff <- 4 / (1 / 60801 + 1 / 123504)

info_snp <- snp_match(sumstats, info)
# 9,455,778 variants to be matched.
# 0 ambiguous SNPs have been removed.
# 1,116,009 variants have been matched; 0 were flipped and 1,011,595 were reversed.
(info_snp <- tidyr::drop_na(tibble::as_tibble(info_snp)))

sd_ldref <- with(info_snp, sqrt(2 * maf * (1 - maf)))
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
# ggsave("figures/sd-approx-CAD.png", width = 10, height = 7)


#### Asthma ####

# download.file("ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/DemenaisF_29273806_GCST006862/TAGC_meta-analyses_results_for_asthma_risk.zip",
#               destfile = "tmp-data/sumstats_asthma.zip")
# unzip("tmp-data/sumstats_asthma.zip", exdir = "tmp-data")
sumstats <- fread2(
  "tmp-data/TAGC meta-analyses results for asthma risk/TAGC_Multiancestry_and_European-Ancestry_Meta-analyses_Results.tsv",
  select = c(1, 3:5, 18:20),
  col.names = c("chr", "pos", "a0", "a1", "beta", "beta_se", "p")
)
sumstats$n_eff <- 4 / (1 / 19954 + 1 / 107715)

info_snp <- snp_match(sumstats, info)
# 2,001,280 variants to be matched.
# 0 ambiguous SNPs have been removed.
# 986,459 variants have been matched; 0 were flipped and 298,357 were reversed.
(info_snp <- tidyr::drop_na(tibble::as_tibble(info_snp)))

sd_ldref <- with(info_snp, sqrt(2 * maf * (1 - maf)))
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
# ggsave("figures/sd-approx-asthma.png", width = 10, height = 7)

