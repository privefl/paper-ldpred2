library(bigsnpr)
ukb <- snp_attach("data/UKBB_imp_HM3.rds")
G <- ukb$genotypes

set.seed(1)
ind.val <- sort(sample(nrow(G), 10e3))
ind.test <- setdiff(rows_along(G), ind.val)

library(tidyverse)
library(future)
(NCORES <- as.integer(Sys.getenv("SLURM_JOB_CPUS_PER_NODE")) - 1L)
plan(multisession, workers = NCORES)

pred_ldpred2 <-
  tibble(res_file = list.files("results/ldpred2", "\\.rds$", full.names = TRUE)) %>%
  mutate(Trait = sub("(.*/)([[:alnum:]]+)(_chr[0-9]+\\.rds)$", "\\2", res_file),
         chr = as.integer(sub("(.*_chr)([0-9]+)(\\.rds)$", "\\2", res_file)),
         pred = furrr::future_map(res_file, readRDS), res_file = NULL) %>%
  unnest_longer(col = "pred", indices_to = "Method") %>%
  filter(Method %in% c("best_beta", "beta_auto")) %>%
  mutate(Method = case_when(Method == "beta_auto" ~ "LDpred2_auto",
                            Method == "best_beta" ~ "LDpred2_best")) %>%
  group_by_at(-(2:3)) %>%
  summarise(pred = list(do.call(plus, map(pred, replace_na, replace = 0)))) %>%
  ungroup() %>%
  print()

pred_sct <-
  tibble(res_file = list.files("results/sct", "\\.rds$", full.names = TRUE)) %>%
  mutate(Trait = sub("\\.rds$", "", basename(res_file)),
         pred = furrr::future_map(res_file, readRDS), res_file = NULL) %>%
  unnest_longer(col = "pred", indices_to = "Method") %>%
  print()

pred_lassosum <-
  tibble(Method = "lassosum",
         res_file = list.files("results/lassosum", "\\.rds$", full.names = TRUE)) %>%
  mutate(Trait = sub("(results/.+/)(.+)(\\.rds$)", "\\2", res_file),
         pred = furrr::future_map(res_file, readRDS), res_file = NULL) %>%
  print()

AUCBoot_no_NA <- function(pred, target) {
  is_not_na <- which(!is.na(target))
  AUCBoot(pred[is_not_na], target[is_not_na])
}

res_boot <- bind_rows(pred_ldpred2, pred_sct, pred_lassosum) %>%
  mutate(AUCBoot = furrr::future_map2(pred, Trait, ~ {
           pheno <- readRDS(paste0("data/pheno/", .y, ".rds"))
           AUCBoot_no_NA(.x, pheno[ind.test])
         }),
         pred = NULL) %>%
  print()

auc_real <- res_boot %>%
  # bind_rows(readRDS("results/all-ldpred1.rds")) %>%
  bind_rows(readRDS("results/all-prscs.rds")) %>%
  unnest_wider("AUCBoot") %>%
  print()

ggplot(auc_real, aes(Trait, Mean, fill = Method, color = Method)) +
  bigstatsr::theme_bigstatsr() +
  geom_hline(yintercept = 0.5, linetype = 2) +
  geom_col(position = position_dodge(), alpha = 0.6, color = "black", size = 1) +
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`),
                position = position_dodge(width = 0.9),
                color = "black", width = 0.2, size = 1) +
  scale_y_continuous(limits = c(0.5, 0.8), minor_breaks = 0:50 / 50,
                     oob = scales::rescale_none) +
  labs(y = "AUC") +
  theme(legend.position = c(0.12, 0.79))

# ggsave("figures/AUC-all.pdf", width = 11, height = 6.5)

auc_real %>%
  mutate_at(3:5, ~ round(. * 100, 1)) %>%
  mutate(AUC = sprintf("%.1f [%.1f-%.1f]", Mean, `2.5%`, `97.5%`)) %>%
  select(-(3:6)) %>%
  pivot_wider(names_from = "Method", values_from = "AUC") %>%
  print(width = Inf) %>%
  xtable::xtable(align = "|l|c|c|c|c|c|c|c|") %>%
  print(include.rownames = FALSE)

#   Trait  LDpred2_auto     LDpred2_best     SCT
# 1 Asthma 58.2 [57.9-58.5] 58.6 [58.3-58.9] 58.3 [58.0-58.6]
# 2 BRCA   65.6 [65.2-66.1] 64.5 [64.0-65.0] 62.3 [61.8-62.8]
# 3 CAD    62.1 [61.6-62.6] 61.5 [61.0-62.0] 60.6 [60.1-61.1]
# 4 MDD    58.9 [58.5-59.2] 57.9 [57.5-58.2] 57.4 [57.1-57.8]
# 5 PRCA   70.2 [69.6-70.7] 68.2 [67.6-68.8] 66.5 [65.9-67.1]
# 6 RA     59.7 [58.9-60.4] 60.0 [59.3-60.7] 57.4 [56.6-58.1]
# 7 T1D    76.6 [74.8-78.5] 77.0 [75.2-78.8] 72.4 [70.5-74.3]
# 8 T2D    63.9 [63.5-64.3] 63.3 [62.9-63.7] 62.4 [62.0-62.9]
#   `C+T`            lassosum         `PRS-CS`
# 1 56.7 [56.4-57.0] 57.6 [57.3-57.8] 57.1 [56.8-57.4]
# 2 62.9 [62.5-63.4] 65.2 [64.7-65.6] 63.3 [62.8-63.7]
# 3 61.6 [61.1-62.1] 62.5 [62.1-63.0] 61.8 [61.4-62.3]
# 4 58.5 [58.1-58.9] 59.1 [58.8-59.5] 57.5 [57.2-57.9]
# 5 67.3 [66.7-67.9] 69.3 [68.7-69.9] 67.5 [66.9-68.0]
# 6 59.1 [58.4-59.8] 59.2 [58.6-60.0] 59.2 [58.5-59.9]
# 7 74.4 [72.5-76.2] 74.4 [72.6-76.2] 74.2 [72.4-76.0]
# 8 59.9 [59.4-60.3] 62.7 [62.3-63.2] 62.2 [61.7-62.6]

# % latex table generated in R 3.5.1 by xtable 1.8-4 package
# % Sun Jun 21 12:54:17 2020
# \begin{table}[ht]
# \centering
# \begin{tabular}{|c|c|c|c|c|c|c|}
# \hline
# Trait & LDpred2\_auto & LDpred2\_best & SCT & C+T & lassosum & PRS-CS \\
# \hline
# Asthma & 58.2 [57.9-58.5] & 58.6 [58.3-58.9] & 58.3 [58.0-58.6] & 56.7 [56.4-57.0] & 57.6 [57.3-57.8] & 57.1 [56.8-57.4] \\
# BRCA & 65.6 [65.2-66.1] & 64.5 [64.0-65.0] & 62.3 [61.8-62.8] & 62.9 [62.5-63.4] & 65.2 [64.7-65.6] & 63.3 [62.8-63.7] \\
# CAD & 62.1 [61.6-62.6] & 61.5 [61.0-62.0] & 60.6 [60.1-61.1] & 61.6 [61.1-62.1] & 62.5 [62.1-63.0] & 61.8 [61.4-62.3] \\
# MDD & 58.9 [58.5-59.2] & 57.9 [57.5-58.2] & 57.4 [57.1-57.8] & 58.5 [58.1-58.9] & 59.1 [58.8-59.5] & 57.5 [57.2-57.9] \\
# PRCA & 70.2 [69.6-70.7] & 68.2 [67.6-68.8] & 66.5 [65.9-67.1] & 67.3 [66.7-67.9] & 69.3 [68.7-69.9] & 67.5 [66.9-68.0] \\
# RA & 59.7 [58.9-60.4] & 60.0 [59.3-60.7] & 57.4 [56.6-58.1] & 59.1 [58.4-59.8] & 59.2 [58.6-60.0] & 59.2 [58.5-59.9] \\
# T1D & 76.6 [74.8-78.5] & 77.0 [75.2-78.8] & 72.4 [70.5-74.3] & 74.4 [72.5-76.2] & 74.4 [72.6-76.2] & 74.2 [72.4-76.0] \\
# T2D & 63.9 [63.5-64.3] & 63.3 [62.9-63.7] & 62.4 [62.0-62.9] & 59.9 [59.4-60.3] & 62.7 [62.3-63.2] & 62.2 [61.7-62.6] \\
# \hline
# \end{tabular}
# \end{table}
