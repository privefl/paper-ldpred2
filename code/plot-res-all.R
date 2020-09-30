library(bigsnpr)
ukb <- snp_attach("data/UKBB_imp_HM3.rds")
G <- ukb$genotypes

load("data/ind_val_test.RData")

library(future)
(NCORES <- as.integer(Sys.getenv("SLURM_JOB_CPUS_PER_NODE")) - 1L)
plan(multisession, workers = NCORES)

library(dplyr)
library(tidyr)
pred_ldpred2_gwide <-
  tibble(res_file = list.files("results/ldpred2-gwide", "\\.rds$", full.names = TRUE)) %>%
  mutate(Trait = sub("\\.rds$", "", basename(res_file)),
         pred = furrr::future_map(res_file, ~ as.list(readRDS(.)$pred)),
         res_file = NULL) %>%
  unnest_longer(col = "pred", indices_to = "Method") %>%
  mutate(Method = paste0(sub("beta", "ldpred2", Method), "_gwide")) %>%
  print()

pred_sct <-
  tibble(res_file = list.files("results/sct", "\\.rds$", full.names = TRUE)) %>%
  mutate(Trait = sub("\\.rds$", "", basename(res_file)),
         pred = furrr::future_map(res_file, readRDS), res_file = NULL) %>%
  unnest_longer(col = "pred", indices_to = "Method") %>%
  print()

pred_lassosum <-
  tibble(res_file = list.files("results/lassosum", "\\.rds$", full.names = TRUE),
         Trait = sub("(results/.+/)(.+)(\\.rds$)", "\\2", res_file),
         pred = furrr::future_map(res_file, readRDS), res_file = NULL) %>%
  unnest_longer(col = "pred", indices_to = "Method") %>%
  select(-res_file) %>%
  print()

AUCBoot_no_NA <- function(pred, target) {
  is_not_na <- which(!is.na(target))
  AUCBoot(pred[is_not_na], target[is_not_na])
}

res_boot <- bind_rows(pred_ldpred2_gwide, pred_sct, pred_lassosum) %>%
  mutate(AUCBoot = furrr::future_map2(pred, Trait, ~ {
           pheno <- readRDS(paste0("data/pheno/", .y, ".rds"))
           AUCBoot_no_NA(.x, pheno[ind.test])
         }),
         pred = NULL) %>%
  print()

auc_real <- res_boot %>%
  bind_rows(readRDS("results/all-prscs.rds")) %>%
  bind_rows(readRDS("results/all-sbayesr.rds")) %>%
  mutate(Method = sub("^best_", "", Method),
         Method = stringr::str_replace_all(Method, "_", "-"),
         Method = sub("ldpred", "LDpred", Method)) %>%
  unnest_wider("AUCBoot") %>%
  print()


library(ggplot2)
# http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/#a-colorblind-friendly-palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
               "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

ggplot({
  auc_real %>%
    mutate(Method = factor(Method, levels = c(
      "LDpred2-inf-gwide", "LDpred2-grid-nosp-gwide", "LDpred2-auto-gwide",
      "SCT", "C+T", "lassosum", "PRS-CS", "SBayesR"))) %>%
    na.omit()
}, aes(Trait, Mean, fill = Method, color = Method)) +
  bigstatsr::theme_bigstatsr(0.9) +
  geom_hline(yintercept = 0.5, linetype = 2) +
  geom_col(position = position_dodge(), alpha = 0.6, color = "black", size = 1) +
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`),
                position = position_dodge(width = 0.9),
                color = "black", width = 0.2, size = 1) +
  scale_y_continuous(limits = c(0.5, 0.8), minor_breaks = 0:50 / 50,
                     oob = scales::rescale_none) +
  labs(y = "AUC") +
  theme(legend.position = c(0.25, 0.8)) +
  guides(fill = guide_legend(ncol = 2)) +
  scale_fill_manual(values = cbPalette)

# ggsave("figures/AUC-all.pdf", width = 12, height = 6.5)


ggplot({
  auc_real %>%
    mutate(Method = factor(Method, levels = c(
      "LDpred2-grid-nosp-gwide", "LDpred2-auto-gwide",
      "lassosum", "lassosum-auto", "PRS-CS", "PRS-CS-auto", "SBayesR"))) %>%
    na.omit()
}, aes(Trait, Mean, fill = Method, color = Method)) +
  bigstatsr::theme_bigstatsr() +
  geom_hline(yintercept = 0.5, linetype = 2) +
  geom_col(position = position_dodge(), alpha = 0.6, color = "black", size = 1) +
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`),
                position = position_dodge(width = 0.9),
                color = "black", width = 0.2, size = 1) +
  scale_y_continuous(limits = c(0.5, 0.8), minor_breaks = 0:50 / 50,
                     oob = scales::rescale_none) +
  labs(y = "AUC") +
  theme(legend.position = c(0.25, 0.8)) +
  guides(fill = guide_legend(ncol = 2)) +
  scale_fill_manual(values = cbPalette[c(2, 5, 6, 3, 7, 8, 4)])

# ggsave("figures/AUC-auto.pdf", width = 12, height = 6.5)


auc_real %>%
  mutate_at(3:5, ~ round(. * 100, 1)) %>%
  mutate(AUC = sprintf("%.1f [%.1f-%.1f]", Mean, `2.5%`, `97.5%`)) %>%
  select(-(3:6)) %>%
  pivot_wider(names_from = "Trait", values_from = "AUC") %>%
  print(width = Inf) %>%
  xtable::xtable(align = "|l|c|c|c|c|c|c|c|c|c|") %>%
  print(include.rownames = FALSE)

#   Method                  Asthma           BRCA             CAD
#   <chr>                   <chr>            <chr>            <chr>
# 1 LDpred2-inf-gwide       57.2 [56.9-57.4] 61.6 [61.1-62.1] 61.6 [61.1-62.0]
# 2 LDpred2-grid-nosp-gwide 59.3 [59.0-59.6] 65.5 [65.0-66.0] 63.6 [63.2-64.1]
# 3 LDpred2-grid-sp-gwide   59.4 [59.1-59.7] 65.7 [65.3-66.2] 63.5 [63.1-64.0]
# 4 LDpred2-auto-gwide      58.4 [58.2-58.7] 65.6 [65.2-66.1] 61.8 [61.3-62.3]
# 5 SCT                     58.3 [58.0-58.6] 62.3 [61.8-62.8] 60.6 [60.1-61.1]
# 6 C+T                     56.7 [56.4-57.0] 62.9 [62.5-63.4] 61.6 [61.1-62.1]
# 7 lassosum                57.6 [57.3-57.8] 65.2 [64.7-65.6] 62.5 [62.1-63.0]
# 8 lassosum-auto           57.6 [57.3-57.8] 65.1 [64.7-65.6] 61.4 [61.0-61.9]
# 9 PRS-CS                  57.1 [56.8-57.4] 63.3 [62.8-63.7] 61.8 [61.4-62.3]
# 10 PRS-CS-auto             56.6 [56.3-56.8] 63.4 [62.9-63.8] 60.9 [60.4-61.3]
# 11 SBayesR                 57.6 [57.4-57.9] 65.7 [65.2-66.2] 62.2 [61.7-62.6]
#   MDD              PRCA             RA               T1D              T2D
#   <chr>            <chr>            <chr>            <chr>            <chr>
# 1 59.0 [58.6-59.3] 64.4 [63.8-65.0] 59.3 [58.6-60.0] 71.0 [69.1-72.9] 61.6 [61.1-62.0]
# 2 59.1 [58.7-59.4] 70.2 [69.6-70.8] 60.3 [59.6-61.0] 78.4 [76.7-80.1] 64.3 [63.9-64.7]
# 3 59.0 [58.7-59.4] 70.2 [69.6-70.7] 60.5 [59.8-61.3] 78.2 [76.5-80.0] 64.0 [63.6-64.5]
# 4 59.0 [58.6-59.4] 70.1 [69.5-70.7] 59.7 [59.0-60.5] 77.7 [75.9-79.4] 63.8 [63.3-64.2]
# 5 57.4 [57.1-57.8] 66.5 [65.9-67.1] 57.4 [56.7-58.1] 72.4 [70.5-74.2] 62.5 [62.0-62.9]
# 6 58.5 [58.1-58.9] 67.3 [66.7-67.9] 59.1 [58.4-59.8] 74.4 [72.6-76.2] 59.9 [59.4-60.3]
# 7 59.1 [58.8-59.5] 69.3 [68.7-69.9] 59.3 [58.5-59.9] 74.3 [72.5-76.1] 62.7 [62.3-63.2]
# 8 54.9 [54.6-55.3] 69.3 [68.7-69.8] 58.2 [57.5-58.9] 75.8 [74.1-77.5] 62.4 [62.0-62.8]
# 9 57.5 [57.2-57.9] 67.5 [66.9-68.1] 59.2 [58.5-59.9] 74.2 [72.4-76.0] 62.2 [61.7-62.6]
# 10 53.9 [53.5-54.3] 67.2 [66.6-67.7] 58.6 [57.9-59.3] 73.9 [72.1-75.7] 62.4 [62.0-62.9]
# 11 58.8 [58.5-59.2] 69.6 [69.1-70.2] 56.2 [55.5-56.9] 58.1 [56.1-60.2] 64.1 [63.7-64.5]
