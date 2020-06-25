# file.symlink("../../UKBB/", "UKBB")

map_hapmap3 <- bigreadr::fread2("ldblk_1kg_eur/snpinfo_1kg_hm3")

library(doParallel)
cl <- makeCluster(22)
parallel::clusterExport(cl, "map_hapmap3")
list_snp_id <- parLapply(cl, 1:22, function(chr) {
  mfi <- paste0("UKBB/mfi/ukb_mfi_chr", chr, "_v3.txt")
  infos_chr <- bigreadr::fread2(mfi, showProgress = FALSE)
  joined <- dplyr::inner_join(cbind(chr = chr, infos_chr), map_hapmap3[1:2],
                              by = c("chr" = "CHR", "V2" = "SNP"))
  with(joined[!vctrs::vec_duplicate_detect(joined$V2), ],
       paste(chr, V3, V4, V5, sep = "_"))
})
stopCluster(cl)

sum(lengths(list_snp_id))  # 1,117,493


sample <- bigreadr::fread2("UKBB/ukb58024_imp_chr1_v3_s487296.sample")
str(sample)
sample <- sample[-1, ]


csv <- "UKBB/ukb41181.csv"
df0 <- bigreadr::fread2(
  csv,
  select = c("eid", "22020-0.0", paste0("22009-0.", 1:16)),
  col.names = c("eid", "used_in_pca", paste0("PC", 1:16))
)

# not removed & quality controled / unrelated & white-British
ind.indiv <- match(df0$eid, sample$ID_2)
sub <- which(!is.na(ind.indiv) & df0$used_in_pca)

# Genetically homogeneous
dist <- bigutilsr::dist_ogk(as.matrix(df0[sub, -(1:2)]))
sub2 <- sub[log(dist) < 5]
length(sub2) # 362,320
saveRDS(sub2, "data/ind_sub_csv.rds")


(NCORES <- as.integer(Sys.getenv("SLURM_JOB_CPUS_PER_NODE")) - 1L)

system.time(
  rds <- bigsnpr::snp_readBGEN(
    bgenfiles   = glue::glue("UKBB/bgen/ukb_imp_chr{chr}_v3.bgen", chr = 1:22),
    list_snp_id = list_snp_id,
    backingfile = "data/UKBB_imp_HM3",
    ind_row     = ind.indiv[sub2],
    ncores      = NCORES
  )
) # 43 min with 23 cores

library(bigsnpr)
ukb <- snp_attach("data/UKBB_imp_HM3.rds")
G <- ukb$genotypes
CHR <- as.numeric(ukb$map$chromosome)
POS <- ukb$map$physical.pos
dim(G) # 362,320 x 1,117,493
file.size(G$backingfile) / 1024^3  # 377 GB

# Write bed file
ukb$map <- dplyr::mutate(ukb$map, chromosome = as.integer(chromosome),
                         marked.ID = rsid, genetic.dist = 0)
ukb$fam <- snp_fake(n = nrow(G), m = 1)$fam
ukb$fam$sample.ID <- sample$ID_2[ind.indiv[sub2]]

set.seed(1); ind.val <- sort(sample(nrow(G), 10e3))
snp_writeBed(ukb, bedfile = "data/UKBB_imp_HM3_val.bed", ind.row = ind.val)
