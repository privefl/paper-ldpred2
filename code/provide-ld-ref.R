library(bigsnpr)
library(bigreadr)
library(dplyr)
library(ggplot2)

#### Filter data from allele frequency errors ####

ukb <- snp_attach("data/UKBB_imp_HM3.rds")
G <- ukb$genotypes
map <- transmute(ukb$map,
                 chr = as.integer(chromosome), pos = physical.pos,
                 a0 = allele1, a1 = allele2, rsid)  # reversed somehow..

(NCORES <- as.integer(Sys.getenv("SLURM_JOB_CPUS_PER_NODE")) - 1L)
af <- runonce::save_run(
  big_colstats(G, ncores = NCORES)$sum / (2 * nrow(G)),
  file = "data/af.rds"
)

obj.bed <- bed(download_1000G("tmp-data"))
map2 <- dplyr::transmute(obj.bed$map, chr = chromosome, pos = physical.pos,
                         a1 = allele1, a0 = allele2, beta = 1)
info_snp <- snp_match(map2, map)
# 1,664,852 variants to be matched.
# 0 ambiguous SNPs have been removed.
# 1,058,123 variants have been matched; 0 were flipped and 0 were reversed.

fam2 <- bigreadr::fread2(sub_bed(obj.bed$bedfile, ".fam2"))
ind_eur <- which(fam2$`Super Population` == "EUR")
af2 <- bed_MAF(obj.bed, ind.row = ind_eur, ind.col = info_snp$`_NUM_ID_.ss`,
               ncores = NCORES)

af_UKBB <- af[info_snp$`_NUM_ID_`]
af_1000G <- af2$af

N <- nrow(G)
af_avg <- (af_UKBB * N + af_1000G * 503) / (N + 503)
chi2 <- (af_UKBB - af_1000G)^2 / (af_avg * (1 - af_avg) * (1 / N + 1 / 503))
hist(pval <- pchisq(chi2, df = 1, lower.tail = FALSE), breaks = 0:50 / 50)

is_bad <- pval < 1e-5 | af_1000G < 0.01 | af_1000G > 0.99 | af_UKBB < 0.005 | af_UKBB > 0.995
qplot(af_UKBB, af_1000G, color = is_bad) +
  theme_bigstatsr() +
  geom_abline(color = "red") +
  coord_equal() +
  scale_color_viridis_d(direction = -1) +
  labs(x = "Allele frequency in UKBB", y = "Allele frequency in 1000G",
       color = "Removed?")

# ggsave("figures/af-outliers.png", width = 9, height = 7)

ind_keep <- info_snp$`_NUM_ID_`[which(!is_bad)]
map_keep <- map[ind_keep, ]
map_keep$af_UKBB <- af_UKBB[which(!is_bad)]
vctrs::vec_duplicate_any(map_keep[, 1:2])  # FALSE

# add positions in different builds
liftOver <- runonce::download_file(
  "http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver", "tmp-data")
bigsnpr:::make_executable(liftOver)
map_keep$pos_hg19 <- snp_modifyBuild(map_keep, liftOver, to = "hg19")$pos
map_keep$pos_hg17 <- snp_modifyBuild(map_keep, liftOver, to = "hg17")$pos

bigassertr::assert_dir("ld-ref")

#### Compute LD matrices ####

CHR <- map_keep$chr
POS2 <- snp_asGeneticPos(CHR, map_keep$pos, dir = "tmp-data", ncores = NCORES)

library(future.batchtools)
NCORES <- 15
plan(batchtools_slurm(resources = list(
  t = "12:00:00", c = NCORES, mem = "125g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

future.apply::future_lapply(1:22, function(chr) {

  ind.chr <- which(CHR == chr)

  corr <- snp_cor(
    G, ind.col = ind_keep[ind.chr], infos.pos = POS2[ind.chr],
    size = 3 / 1000, ncores = NCORES
  )

  saveRDS(corr, file = paste0("ld-ref/LD_chr", chr, ".rds"), version = 2)
})
# 10 min for chr 22 with 15,449 variants

sum(sapply(list.files("ld-ref", full.names = TRUE), file.size)) / 1024^3
# 8.5 GB

# Compute LD scores
map_keep$ld <- do.call('c', lapply(1:22, function(chr) {
  cat(chr, ".. ", sep = "")
  corr_chr <- readRDS(paste0("ld-ref/LD_chr", chr, ".rds"))
  Matrix::colSums(corr_chr^2)
}))

saveRDS(map_keep, "ld-ref/map.rds", version = 2)
