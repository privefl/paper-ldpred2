sumstats <- bigreadr::fread2("https://data.broadinstitute.org/alkesgroup/UKBB/UKBB_409K/cov_SMOKING_STATUS.sumstats.gz")

hm3_info <- readRDS(url("https://github.com/privefl/bigsnpr/raw/master/data-raw/hm3_variants.rds"))

match_info <- bigsnpr::snp_match(
  sumstats = dplyr::rename(sumstats, chr = CHR, pos = POS, a1 = A1, a0 = A2, beta = Beta),
  info_snp = hm3_info
)
# 11,904,913 variants to be matched.
# 1,790,164 ambiguous SNPs have been removed.
# 1,119,564 variants have been matched; 0 were flipped and 783,941 were reversed.

# https://github.com/privefl/paper-ldpred2/blob/master/paper/paper-ldpred2-supp.pdf
df <- dplyr::mutate(match_info,
                    sd_af = sqrt(2 * EAF * (1 - EAF)),
                    sd_ss = 2 / (se * sqrt(N)))

# Reduce points for plotting (but try to keep outliers)
hist(ratio <- with(df, sd_ss / sd_af))
hist(out <- abs(ratio - median(ratio)) / mad(ratio))
ind <- sample(nrow(df), 50e3, prob = out)

library(ggplot2)
theme_set(bigstatsr::theme_bigstatsr())

ggplot(df[ind, ]) +
  geom_point(aes(sd_af, sd_ss))
# ggsave("tmp-figures/qcplot1.png", width = 7, height = 6)

ggplot(df[ind, ]) +
  geom_point(aes(sd_af, sd_ss, color = INFO)) +
  scale_color_viridis_c()
# ggsave("tmp-figures/qcplot2.png", width = 7, height = 6)

ggplot(df[ind, ]) +
  geom_point(aes(sd_af, sd_ss / sqrt(INFO), color = INFO)) +
  scale_color_viridis_c()
# ggsave("tmp-figures/qcplot3.png", width = 7, height = 6)

ggplot(df[ind, ]) +
  geom_point(aes(sd_af, sd_ss / sqrt(INFO),
                 color = ifelse(chr %in% c(6, 8, 11), chr, "Other"))) +
  labs(color = "Chr")
# ggsave("tmp-figures/qcplot4.png", width = 7, height = 6)
