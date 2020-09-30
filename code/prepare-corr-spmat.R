library(bigsnpr)
ukb <- snp_attach("data/UKBB_imp_HM3.rds")
G <- ukb$genotypes
CHR <- as.integer(ukb$map$chromosome)
POS <- ukb$map$physical.pos

load("data/ind_val_test.RData")

(NCORES <- as.integer(Sys.getenv("SLURM_JOB_CPUS_PER_NODE")) - 1L)
POS2 <- snp_asGeneticPos(CHR, POS, dir = "tmp-data", ncores = NCORES)

bigassertr::assert_dir("data/corr")

for (chr in 1:22) {

  ind.chr <- which(CHR == chr)

  runonce::save_run(
    snp_cor(
      G, ind.row = ind.val, ind.col = ind.chr,
      alpha = 1, infos.pos = POS2[ind.chr], size = 3 / 1000,
      ncores = NCORES
    ),
    file = paste0("data/corr/chr", chr, ".rds")
  )
}
## On a node with 16 cores and 128 GB RAM:
#     user   system  elapsed
# 6737.315   31.980  650.724
# 6929.821   40.563  665.379
# 7321.282   39.027  644.371
# 5373.981   32.860  486.469
# 5638.935 1368.178  632.483
# 7840.772   52.543  693.595
# 4477.675   34.834  424.483
# 4875.977   37.007  448.571
# 3611.675   26.642  325.730
# 4434.387  542.683  445.706
# 4250.445   32.365  383.651
# 4057.731  801.579  424.581
# 3038.082   19.844  265.982
# 2445.123   17.932  220.284
# 2182.027   15.725  193.493
# 1982.660   16.003  182.821
# 1632.852   15.099  147.933
# 2177.030  372.233  216.424
#  857.160    7.769   80.343
# 1803.857   13.320  154.808
#  788.966    6.565   70.550
#  861.304    6.443   72.679
# -> 2.2H in total

sum(file.size(paste0("data/corr/chr", 1:22, ".rds"))) / 1024^3
# -> 9.5 GB total
