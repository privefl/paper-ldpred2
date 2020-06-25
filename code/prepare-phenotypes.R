library(tidyverse)
library(bigreadr)

bigassertr::assert_dir("data/pheno")

csv <- "UKBB/ukb41181.csv"
sub <- readRDS("data/ind_sub_csv.rds")

# Cancers
df_cancer0 <- fread2(csv, select = paste0("2453-", 0:2, ".0"))
df_cancer1 <- fread2(csv, select = c(paste0("20001-0.", 0:5),
                                     paste0("20001-1.", 0:5),
                                     paste0("20001-2.", 0:5)))
df_cancer2 <- fread2(csv, select = c(paste0("40001-", 0:2, ".0"),
                                     paste0("40002-0.", 0:13),
                                     paste0("40002-1.", 0:13),
                                     paste0("40002-2.", 0:13),
                                     paste0("40006-", 0:31, ".0"),
                                     paste0("41202-0.", 0:379),
                                     paste0("41204-0.", 0:434)))

# Other illnesses
df_illness <- fread2(csv, select = c(paste0("20002-0.", 0:28),
                                     paste0("20002-1.", 0:28),
                                     paste0("20002-2.", 0:28)))
df_ICD10 <- fread2(csv, select = c(paste0("40001-", 0:2, ".0"),
                                   paste0("40002-0.", 0:13),
                                   paste0("40002-1.", 0:13),
                                   paste0("40002-2.", 0:13),
                                   paste0("41202-0.", 0:379),
                                   paste0("41204-0.", 0:434)))

#### Breast cancer (BRCA) ####

ind_BRCA <- sort(unique(unlist(c(
  lapply(df_cancer1, function(x) which(x == 1002)),
  lapply(df_cancer2, function(x) which(substr(x, 1, 3) %in% c("C50", "D05")))
))))

y <- rep(NA, nrow(df_cancer0))
y[rowMeans(df_cancer0 == 0, na.rm = TRUE) == 1] <- 0
y[ind_BRCA] <- 1

sex <- fread2(csv, select = "22001-0.0")[[1]]
y[sex == 1] <- NA  ## exclude men for breast cancer

table(y[sub], exclude = NULL)
#      0      1   <NA>
# 169807  13932 178581

saveRDS(y[sub], "data/pheno/BRCA.rds")


#### Rheumatoid arthritis (RA) ####

ind_RA <- sort(unique(unlist(c(
  lapply(df_illness, function(x) which(x == 1464)),
  lapply(df_ICD10,   function(x) which(substr(x, 1, 3) %in% c("M05", "M06")))
))))
ind_muscu <- sort(unique(unlist(c(
  lapply(df_illness, function(x) which(x %in% c(1295, 1464:1467, 1477, 1538))),
  lapply(df_ICD10,   function(x) which(substr(x, 1, 1) == "M"))
))))
y <- rep(0, nrow(df_illness)); y[ind_muscu] <- NA; y[ind_RA] <- 1

table(y[sub], exclude = NULL)
#      0      1   <NA>
# 231938   6819 123563

saveRDS(y[sub], "data/pheno/RA.rds")


#### Type 1 diabetes (T1D) ####

ind_diabetes <- sort(unique(unlist(c(
  lapply(df_illness, function(x) which(x %in% 1220:1223)),
  lapply(df_ICD10, function(x) which(substr(x, 1, 3) %in% paste0("E", 10:14)))
))))
ind_TD1 <- sort(unique(unlist(c(
  lapply(df_illness, function(x) which(x == 1222)),
  lapply(df_ICD10, function(x) which(substr(x, 1, 3) == "E10"))
))))
ind_TD2 <- sort(unique(unlist(c(
  lapply(df_illness, function(x) which(x == 1223)),
  lapply(df_ICD10, function(x) which(substr(x, 1, 3) == "E11"))
))))
y <- rep(0, nrow(df_illness))
y[ind_diabetes] <- NA
y[ind_TD1] <- 1
y[ind_TD2] <- NA

table(y[sub], exclude = NULL)
#      0      1   <NA>
# 337087    823  24410

saveRDS(y[sub], "data/pheno/T1D.rds")


#### Type 2 diabetes (T2D) ####

y <- rep(0, nrow(df_illness))
y[ind_diabetes] <- NA
y[ind_TD2] <- 1
y[ind_TD1] <- NA

table(y[sub], exclude = NULL)
#      0      1   <NA>
# 337087  18182   7051

saveRDS(y[sub], "data/pheno/T2D.rds")


#### Prostate cancer (PRCA) ####

ind_PRCA <- sort(unique(unlist(c(
  lapply(df_cancer1, function(x) which(x == 1044)),
  lapply(df_cancer2, function(x) which(x %in% c("C61", "D075")))
))))
y <- rep(NA, nrow(df_cancer0))
y[rowMeans(df_cancer0 == 0, na.rm = TRUE) == 1] <- 0
y[ind_PRCA] <- 1

sex <- fread2(csv, select = "22001-0.0")[[1]]
y[sex == 0] <- NA  ## exclude women for prostate cancer

table(y[sub], exclude = NULL)
#      0      1   <NA>
# 149974   8660 203686

saveRDS(y[sub], "data/pheno/PRCA.rds")


#### Depression (MDD) ####

ind_MDD <- sort(unique(unlist(c(
  lapply(df_illness, function(x) which(x == 1286)),
  lapply(df_ICD10, function(x) which(substr(x, 1, 3) %in% c("F32", "F33")))
))))
ind_psy <- sort(unique(unlist(c(
  lapply(df_illness, function(x) which(x %in% 1286:1291)),
  lapply(df_ICD10, function(x) which(substr(x, 1, 1) == "F"))
))))
y <- rep(0, nrow(df_illness)); y[ind_psy] <- NA; y[ind_MDD] <- 1

batch <- fread2(csv, select = "22000-0.0")[[1]]
y[batch < 0] <- NA  ## exclude pilot data as it was included in sumstats

table(y[sub], exclude = NULL)
#      0      1   <NA>
# 275585  27037  59698

saveRDS(y[sub], "data/pheno/MDD.rds")


#### Coronary artery disease (CAD) ####

df_heart <- fread2(csv, select = c(paste0("6150-0.", 0:3),
                                   paste0("6150-1.", 0:3),
                                   paste0("6150-2.", 0:3)))
ind_CAD <- sort(unique(unlist(c(
  lapply(df_heart,   function(x) which(x == 1)),
  lapply(df_illness, function(x) which(x == 1075)),
  lapply(df_ICD10,   function(x) which(substr(x, 1, 3) %in% paste0("I", 21:23))),
  lapply(df_ICD10,   function(x) which(x == "I252"))
))))
ind_heart <- sort(unique(unlist(c(
  lapply(df_heart,   function(x) which(x %in% 1:3)),
  lapply(df_illness, function(x) which(x %in% 1074:1080)),
  lapply(df_ICD10,   function(x) which(substr(x, 1, 1) == "I"))
))))
y <- rep(0, nrow(df_heart)); y[ind_heart] <- NA; y[ind_CAD] <- 1

table(y[sub], exclude = NULL)
#      0      1   <NA>
# 228915  15099 118306

saveRDS(y[sub], "data/pheno/CAD.rds")


#### Asthma ####

ind_asthma <- sort(unique(unlist(c(
  lapply(df_illness, function(x) which(x == 1111)),
  lapply(df_ICD10,   function(x) which(substr(x, 1, 3) == "J45"))
))))
ind_respi <- sort(unique(unlist(c(
  lapply(df_illness, function(x) which(x %in% 1111:1125)),
  lapply(df_ICD10,   function(x) which(substr(x, 1, 1) == "J"))
))))
y <- rep(0, nrow(df_illness)); y[ind_respi] <- NA; y[ind_asthma] <- 1

table(y[sub], exclude = NULL)
#      0      1   <NA>
# 275392  48688  38240

saveRDS(y[sub], "data/pheno/Asthma.rds")
