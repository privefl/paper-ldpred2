library(bigsnpr)
ukb <- snp_attach("data/UKBB_imp_HM3.rds")
G <- ukb$genotypes
CHR <- as.integer(ukb$map$chromosome)
POS <- ukb$map$physical.pos
(NCORES <- as.integer(Sys.getenv("SLURM_JOB_CPUS_PER_NODE")) - 1L)

bigassertr::assert_dir("data/pheno-simu")

#### Anywhere on the genome ####

params <- expand.grid(
  h2   = 0.4,
  M    = c(300, 3e3, 30e3, 300e3),
  K    = 0.15,
  iter = 1:10
)

purrr::pwalk(params, function(h2, M, K, iter) {

  pheno_file <- sprintf("data/pheno-simu/all_%d_%d_%d_%d.rds",
                        100 * h2, M, 100 * K, iter)

  if (!file.exists(pheno_file)) {
    print(pheno_file)
    simu <- snp_simuPheno(G, h2, M, K, ncores = NCORES)
    saveRDS(simu$pheno, pheno_file)
  }
})


#### In HLA ####

ind.HLA <- snp_indLRLDR(CHR, POS, LD.wiki34[12, ])

params2 <- expand.grid(
  h2   = 0.3,
  M    = c(300, 3000),
  K    = 0.15,
  iter = 1:10
)

purrr::pwalk(params2, function(h2, M, K, iter) {

  pheno_file <- sprintf("data/pheno-simu/HLA_%d_%d_%d_%d.rds",
                        100 * h2, M, 100 * K, iter)

  if (!file.exists(pheno_file)) {
    print(pheno_file)
    simu <- snp_simuPheno(G, h2, M, K, ind.possible = ind.HLA, ncores = NCORES)
    saveRDS(simu$pheno, pheno_file)
  }
})


#### In both ####

h2 <- 0.4
M_HLA <- 300
M_all <- 10e3
K <- 0.15

for (iter in 1:10) {

  pheno_file <- sprintf("data/pheno-simu/both_%d_%d_%d.rds",
                        100 * h2, 100 * K, iter)

  if (!file.exists(pheno_file)) {
    print(pheno_file)
    simu1 <- snp_simuPheno(G, h2, M_HLA, ind.possible = ind.HLA)
    simu2 <- snp_simuPheno(G, h2, M_all, ncores = NCORES)
    simu_liab <- (simu1$pheno + simu2$pheno) / sqrt(2)
    simu_pheno <- (simu_liab > qnorm(K, lower.tail = FALSE)) + 0L
    saveRDS(simu_pheno, pheno_file)
  }
}
