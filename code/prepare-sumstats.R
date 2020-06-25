library(bigsnpr)
library(bigreadr)
library(dplyr)
library(ggplot2)

ukb <- snp_attach("data/UKBB_imp_HM3.rds")
G <- ukb$genotypes
map <- transmute(ukb$map,
                 chr = as.integer(chromosome), pos = physical.pos,
                 a0 = allele1, a1 = allele2)  # reversed somehow..

bigassertr::assert_dir("data/sumstats")

set.seed(1); ind.val <- sort(sample(nrow(G), 10e3))
(NCORES <- as.integer(Sys.getenv("SLURM_JOB_CPUS_PER_NODE")) - 1L)
sd <- runonce::save_run(
  big_parallelize(G, function(X, ind, ind.val) {
    bigstatsr::big_scale()(X, ind.val, ind)$scale
  }, p.combine = "c", ncores = NCORES, ind.val = ind.val),
  file = "data/sd.rds"
)

#### Breast cancer (BRCA) ####

# download.file("ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/MichailidouK_29059683_GCST004988/oncoarray_bcac_public_release_oct17%20(1).txt.gz",
#               destfile = "tmp-data/sumstats_BRCA.txt.gz")
# R.utils::gunzip("tmp-data/sumstats_BRCA.txt.gz")
sumstats <- fread2("tmp-data/sumstats_BRCA.txt", na.strings = "NULL",
                   select = c("chr", "position_b37", "a0", "a1",
                              "bcac_onco_icogs_gwas_eaf_controls",
                              "bcac_onco_icogs_gwas_beta",
                              "bcac_onco_icogs_gwas_se",
                              "bcac_onco_icogs_gwas_P1df"),
                   col.names = c("chr", "pos", "a0", "a1", "freq",
                                 "beta", "beta_se", "p"))
sumstats$n_eff <- 4 / (1 / 137045 + 1 / 119078)

info_snp <- snp_match(sumstats, map)
# 11,792,542 variants to be matched.
# 1,633,946 ambiguous SNPs have been removed.
# 1,096,979 variants have been matched; 0 were flipped and 0 were reversed.
info_snp <- bigsnpr::snp_match(sumstats, map, strand_flip = FALSE)
# 11,792,542 variants to be matched.
# 1,188,767 variants have been matched; 0 were flipped and 0 were reversed.
(info_snp <- tidyr::drop_na(as_tibble(info_snp)))

sd_val <- sd[info_snp$`_NUM_ID_`]
sd_ss <- with(info_snp, 1 / sqrt(n_eff / 4 * beta_se^2))

is_bad <-
  sd_ss < (0.5 * sd_val) | sd_ss > (sd_val + 0.1) | sd_ss < 0.1 | sd_val < 0.05
qplot(sd_val, sd_ss, color = is_bad, alpha = I(0.5)) +
  theme_bigstatsr() +
  coord_equal() +
  scale_color_viridis_d(direction = -1) +
  geom_abline(linetype = 2, color = "red") +
  labs(x = "Standard deviations in the validation set",
       y = "Standard deviations derived from the summary statistics",
       color = "Removed?")

# ggsave("figures/sd-approx-BRCA.png", width = 10, height = 7)

info_snp[is_bad, ] %>%
  arrange(p) %>%
  head(20) %>%
  mutate(freq2 = big_scale()(G, ind.row = ind.val, ind.col = `_NUM_ID_`)$center / 2) %>%
  select(-`_NUM_ID_.ss`, -`_NUM_ID_`)

saveRDS(info_snp[!is_bad, ], "data/sumstats/BRCA.rds")


#### Rheumatoid arthritis (RA) ####

# download.file("ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/OkadaY_24390342_GCST002318/RA_GWASmeta_TransEthnic_v2.txt.gz",
#               destfile = "tmp-data/sumstats_RA.txt.gz")
# R.utils::gunzip("tmp-data/sumstats_RA.txt.gz")
sumstats <- fread2("tmp-data/sumstats_RA.txt", select = 2:7,
                   col.names = c("chr", "pos", "a1", "a0", "or", "p"))
sumstats <- sumstats %>%
  mutate(beta = log(or),
         chi2 = qchisq(p, df = 1, lower.tail = FALSE),
         beta_se = ifelse(chi2 > 1e-4, abs(beta) / sqrt(chi2), NA),
         n_eff = 4 / (1 / 29880 + 1 / 73758))

info_snp <- snp_match(sumstats, map)
# 9,739,303 variants to be matched.
# 1,498,996 ambiguous SNPs have been removed.
# 1,090,905 variants have been matched; 0 were flipped and 530,577 were reversed.
info_snp <- bigsnpr::snp_match(sumstats, map, strand_flip = FALSE)
# 9,739,303 variants to be matched.
# 1,182,239 variants have been matched; 0 were flipped and 576,124 were reversed.
(info_snp <- tidyr::drop_na(as_tibble(info_snp)))

sd_val <- sd[info_snp$`_NUM_ID_`]
sd_ss <- with(info_snp, 1 / sqrt(n_eff / 4 * beta_se^2))

is_bad <-
  sd_ss < (0.5 * sd_val) | sd_ss > (sd_val + 0.1) | sd_ss < 0.1 | sd_val < 0.05
qplot(sd_val, sd_ss, color = is_bad, alpha = I(0.5)) +
  theme_bigstatsr() +
  scale_color_viridis_d(direction = -1) +
  geom_abline(linetype = 2)

info_snp[is_bad, ] %>%
  arrange(p) %>%
  head(20) %>%
  mutate(freq2 = big_scale()(G, ind.row = ind.val, ind.col = `_NUM_ID_`)$center / 2) %>%
  select(-`_NUM_ID_.ss`, -`_NUM_ID_`)

saveRDS(info_snp[!is_bad, ], "data/sumstats/RA.rds")


#### Type 1 diabetes (T1D) ####

# To request a download link, go to
# https://datadryad.org/stash/dataset/doi:10.5061/dryad.ns8q3
# download.file("http://merritt.cdlib.org/cloudcontainer/mrtstore2/35227889.tar.gz",
#               destfile = "tmp-data/sumstats_T1D.tar.gz")
# untar("tmp-data/sumstats_T1D.tar.gz", exdir = "tmp-data")
sumstats <- bigreadr::fread2(
  paste0("tmp-data/meta_chr_", 1:22),
  select = c("chromosome", "position", "a0", "a1", "beta.meta", "se.meta",
             "p.meta", "qc.check", "EUR_MAF_1kG"),
  col.names = c("chr", "pos", "a0", "a1", "beta", "beta_se", "p", "qc", "maf"))

sumstats2 <- sumstats %>%
  filter(qc == "PASS") %>%
  tidyr::drop_na() %>%
  mutate(n_eff = 4 / (1 / 5913 + 1 / 8828), qc = NULL) %>%
  filter(maf > (1 / sqrt(n_eff))) %>%
  as_tibble() %>%
  print()

info_snp <- snp_match(sumstats2, map)
# 4,187,590 variants to be matched.
# 596,370 ambiguous SNPs have been removed.
# 502,114 variants have been matched; 0 were flipped and 0 were reversed.
info_snp <- bigsnpr::snp_match(sumstats2, map, strand_flip = FALSE)
# 4,187,590 variants to be matched.
# 548,659 variants have been matched; 0 were flipped and 0 were reversed.
(info_snp <- as_tibble(info_snp))
table(info_snp$chr)

sd_val <- sd[info_snp$`_NUM_ID_`]
sd_ss <- with(info_snp, 1 / sqrt(n_eff / 4 * beta_se^2))

is_bad <-
  sd_ss < (0.5 * sd_val) | sd_ss > (sd_val + 0.1) | sd_ss < 0.1 | sd_val < 0.05
qplot(sd_val, sd_ss, color = is_bad, alpha = I(0.5)) +
  theme_bigstatsr() +
  scale_color_viridis_d(direction = -1) +
  geom_abline(linetype = 2)

info_snp[is_bad, ] %>%
  arrange(p) %>%
  head(20) %>%
  mutate(freq2 = big_scale()(G, ind.row = ind.val, ind.col = `_NUM_ID_`)$center / 2) %>%
  select(-`_NUM_ID_.ss`, -`_NUM_ID_`)

saveRDS(info_snp[!is_bad, ], "data/sumstats/T1D.rds")


#### Type 2 diabetes (T2D) ####

# Download from http://diagram-consortium.org/downloads.html
# DIAGRAM 1000G GWAS meta-analysis Stage 1 Summary statistics
# Published in Scott et al (2017)
# unzip("tmp-data/METAANALYSIS_DIAGRAM_SE1.zip", exdir = "tmp-data")
sumstats <- fread2("tmp-data/METAANALYSIS_DIAGRAM_SE1.txt")
sumstats <- tidyr::separate(sumstats, "Chr:Position", c("chr", "pos"), convert = TRUE)
names(sumstats) <- c("chr", "pos", "a1", "a0", "beta", "beta_se", "p", "N")
sumstats$n_eff <- 4 / (1 / 26676 + 1 / 132532)

info_snp <- snp_match(sumstats, map)
# 12,056,346 variants to be matched.
# 1,835,114 ambiguous SNPs have been removed.
# 1,099,430 variants have been matched; 213 were flipped and 534,429 were reversed.
info_snp <- bigsnpr::snp_match(sumstats, map, strand_flip = FALSE)
# 12,056,346 variants to be matched.
# 1,191,209 variants have been matched; 0 were flipped and 580,358 were reversed.
(info_snp <- as_tibble(info_snp))

sd_val <- sd[info_snp$`_NUM_ID_`]
sd_ss <- with(info_snp, 1 / sqrt(n_eff / 4 * beta_se^2))

is_bad <-
  sd_ss < (0.5 * sd_val) | sd_ss > (sd_val + 0.1) | sd_ss < 0.1 | sd_val < 0.05
qplot(sd_val, sd_ss, color = is_bad, alpha = I(0.5)) +
  theme_bigstatsr() +
  scale_color_viridis_d(direction = -1) +
  geom_abline(linetype = 2)

info_snp[is_bad, ] %>%
  arrange(p) %>%
  head(20) %>%
  mutate(freq2 = big_scale()(G, ind.row = ind.val, ind.col = `_NUM_ID_`)$center / 2) %>%
  select(-`_NUM_ID_.ss`, -`_NUM_ID_`)

saveRDS(info_snp[!is_bad, ], "data/sumstats/T2D.rds")


#### Prostate cancer (PRCA) ####

# download.file("http://practical.icr.ac.uk/blog/wp-content/uploads/uploadedfiles/oncoarray/MetaSummaryData/meta_v3_onco_euro_overall_ChrAll_1_release.zip",
#               destfile = "tmp-data/sumstats_PRCA.zip")
# unzip("tmp-data/sumstats_PRCA.zip", exdir = "tmp-data")
sumstats <- fread2(
  "tmp-data/meta_v3_onco_euro_overall_ChrAll_1_release.txt",
  select = c("Chr", "position", "Allele1", "Allele2", "Effect", "StdErr",
             "Pvalue", "Freq1"),
  col.names = c("chr", "pos", "a1", "a0", "beta", "beta_se", "p", "freq")
)
sumstats2 <- sumstats %>%
  mutate(a0 = toupper(a0), a1 = toupper(a1),
         n_eff = 4 / (1 / 79148 + 1 / 61106)) %>%
  filter(pmin(freq, 1 - freq) > (1 / sqrt(n_eff)))

info_snp <- snp_match(sumstats2, map)
# 13,611,076 variants to be matched.
# 1,889,205 ambiguous SNPs have been removed.
# 1,101,409 variants have been matched; 454 were flipped and 535,358 were reversed.
(info_snp <- as_tibble(info_snp))

sd_val <- sd[info_snp$`_NUM_ID_`]
sd_ss <- with(info_snp, 1 / sqrt(n_eff / 4 * beta_se^2))

is_bad <-
  sd_ss < (0.5 * sd_val) | sd_ss > (sd_val + 0.1) | sd_ss < 0.1 | sd_val < 0.05
qplot(sd_val, sd_ss, color = is_bad, alpha = I(0.5)) +
  theme_bigstatsr() +
  scale_color_viridis_d(direction = -1) +
  geom_abline(linetype = 2)

info_snp[is_bad, ] %>%
  arrange(p) %>%
  head(20) %>%
  mutate(freq2 = big_scale()(G, ind.row = ind.val, ind.col = `_NUM_ID_`)$center / 2) %>%
  select(-`_NUM_ID_.ss`, -`_NUM_ID_`)

saveRDS(info_snp[!is_bad, ], "data/sumstats/PRCA.rds")


#### Depression (MDD) ####

# download.file("ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/WrayNR_29700475_GCST005839/MDD2018_ex23andMe.gz",
#               destfile = "tmp-data/sumstats_MDD.txt.gz")
# R.utils::gunzip("tmp-data/sumstats_MDD.txt.gz")
sumstats <- fread2("tmp-data/sumstats_MDD.txt", fill = TRUE,
                   select = c("CHR", "BP", "A1", "A2", "OR", "SE", "P", "Nca", "Nco"),
                   col.names = c("chr", "pos", "a1", "a0", "or", "beta_se", "p", "Nca", "Nco"))
sumstats <- sumstats %>%
  as_tibble() %>%
  mutate(beta = log(or), or = NULL, chr = as.integer(chr),
         n_eff = 4 / (1 / Nca + 1 / Nco), Nca = NULL, Nco = NULL) %>%
  tidyr::drop_na() %>%
  filter(n_eff > (0.5 * max(n_eff))) %>%
  print()

info_snp <- snp_match(sumstats, map)
# 8,511,960 variants to be matched.
# 1,223,174 ambiguous SNPs have been removed.
# 1,080,828 variants have been matched; 0 were flipped and 525,366 were reversed.
info_snp <- snp_match(sumstats, map, strand_flip = FALSE)
# 8,511,960 variants to be matched.
# 1,171,119 variants have been matched; 0 were flipped and 570,586 were reversed.
(info_snp <- as_tibble(info_snp))

sd_val <- sd[info_snp$`_NUM_ID_`]
sd_ss <- with(info_snp, 1 / sqrt(n_eff / 4 * beta_se^2))

is_bad <-
  sd_ss < (0.5 * sd_val) | sd_ss > (sd_val + 0.1) | sd_ss < 0.1 | sd_val < 0.05
qplot(sd_val, sd_ss, color = is_bad, alpha = I(0.5)) +
  theme_bigstatsr() +
  scale_color_viridis_d(direction = -1) +
  geom_abline(linetype = 2)

info_snp[is_bad, ] %>%
  arrange(p) %>%
  head(20) %>%
  mutate(freq2 = big_scale()(G, ind.row = ind.val, ind.col = `_NUM_ID_`)$center / 2) %>%
  select(-`_NUM_ID_.ss`, -`_NUM_ID_`)

saveRDS(info_snp[!is_bad, ], "data/sumstats/MDD.rds")


#### Coronary artery disease (CAD) ####

# download.file("ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/NikpayM_26343387_GCST003116/cad.add.160614.website.txt",
#               destfile = "tmp-data/sumstats_CAD.txt")
sumstats <- fread2("tmp-data/sumstats_CAD.txt",
                   select = c("chr", "bp_hg19", "noneffect_allele",
                              "effect_allele", "beta", "se_dgc", "p_dgc"),
                   col.names = c("chr", "pos", "a0", "a1", "beta", "beta_se", "p"))
sumstats$n_eff <- 4 / (1 / 60801 + 1 / 123504)

info_snp <- snp_match(sumstats, map)
# 9,455,778 variants to be matched.
# 1,332,725 ambiguous SNPs have been removed.
# 1,112,939 variants have been matched; 0 were flipped and 771,046 were reversed.
info_snp <- snp_match(sumstats, map, strand_flip = FALSE)
# 9,455,778 variants to be matched.
# 1,112,939 variants have been matched; 0 were flipped and 771,046 were reversed.
(info_snp <- as_tibble(info_snp))

sd_val <- sd[info_snp$`_NUM_ID_`]
sd_ss <- with(info_snp, 1 / sqrt(n_eff / 4 * beta_se^2))

is_bad <-
  sd_ss < (0.5 * sd_val) | sd_ss > (sd_val + 0.1) | sd_ss < 0.1 | sd_val < 0.05
qplot(sd_val, sd_ss, color = is_bad, alpha = I(0.5)) +
  theme_bigstatsr() +
  scale_color_viridis_d(direction = -1) +
  geom_abline(linetype = 2)

info_snp[is_bad, ] %>%
  arrange(p) %>%
  head(20) %>%
  mutate(freq2 = big_scale()(G, ind.row = ind.val, ind.col = `_NUM_ID_`)$center / 2) %>%
  select(-`_NUM_ID_.ss`, -`_NUM_ID_`)

saveRDS(info_snp[!is_bad, ], "data/sumstats/CAD.rds")


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

info_snp <- snp_match(sumstats, map)
# 2,001,280 variants to be matched.
# 0 ambiguous SNPs have been removed.
# 983,838 variants have been matched; 0 were flipped and 0 were reversed.
(info_snp <- as_tibble(info_snp))

sd_val <- sd[info_snp$`_NUM_ID_`]
sd_ss <- with(info_snp, 1 / sqrt(n_eff / 4 * beta_se^2))

is_bad <-
  sd_ss < (0.5 * sd_val) | sd_ss > (sd_val + 0.1) | sd_ss < 0.1 | sd_val < 0.05
qplot(sd_val, sd_ss, color = is_bad, alpha = I(0.5)) +
  theme_bigstatsr() +
  scale_color_viridis_d(direction = -1) +
  geom_abline(linetype = 2)

info_snp[is_bad, ] %>%
  arrange(p) %>%
  head(20) %>%
  mutate(freq2 = big_scale()(G, ind.row = ind.val, ind.col = `_NUM_ID_`)$center / 2) %>%
  select(-`_NUM_ID_.ss`, -`_NUM_ID_`)

saveRDS(info_snp[!is_bad, ], "data/sumstats/Asthma.rds")
