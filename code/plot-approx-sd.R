library(bigsnpr)
ukb <- snp_attach("data/UKBB_imp_HM3.rds")
G <- ukb$genotypes

load("data/ind_gwas_val_test.RData")

sd <- readRDS("data/sd.rds")

#### Simulations ####

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

## See prepare-sumstats.R for BRCA
