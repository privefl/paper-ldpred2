library(bigsnpr)
ukb <- snp_attach("data/UKBB_imp_HM3.rds")
G <- ukb$genotypes

load("data/ind_val_test.RData")

library(dplyr)
library(tidyr)
library(future)
(NCORES <- as.integer(Sys.getenv("SLURM_JOB_CPUS_PER_NODE")) - 1L)
plan(multisession, workers = NCORES)

pred_ldpred2 <-
  tibble(res_file = list.files("results/ldpred2", "\\.rds$", full.names = TRUE)) %>%
  mutate(Trait = sub("(.*/)([[:alnum:]]+)(_chr[0-9]+\\.rds)$", "\\2", res_file),
         chr = as.integer(sub("(.*_chr)([0-9]+)(\\.rds)$", "\\2", res_file)),
         pred = furrr::future_map(res_file, readRDS), res_file = NULL) %>%
  unnest_longer(col = "pred", indices_to = "Method") %>%
  filter(Method != "best_beta") %>%
  group_by_at(-(2:3)) %>%
  summarise(pred = list(do.call(plus, lapply(pred, replace_na, replace = 0)))) %>%
  ungroup() %>%
  mutate(Method = paste0(sub("beta", "ldpred2", Method), "_perchr")) %>%
  print()

pred_ldpred2_gwide <-
  tibble(res_file = list.files("results/ldpred2-gwide", "\\.rds$", full.names = TRUE)) %>%
  mutate(Trait = sub("\\.rds$", "", basename(res_file)),
         pred = furrr::future_map(res_file, ~ as.list(readRDS(.)$pred)),
         res_file = NULL) %>%
  unnest_longer(col = "pred", indices_to = "Method") %>%
  mutate(Method = paste0(sub("beta", "ldpred2", Method), "_gwide")) %>%
  print()


AUCBoot_no_NA <- function(pred, target) {
  is_not_na <- which(!is.na(target))
  AUCBoot(pred[is_not_na], target[is_not_na])
}

auc_real <- bind_rows(pred_ldpred2, pred_ldpred2_gwide) %>%
  mutate(AUCBoot = furrr::future_map2(pred, Trait, ~ {
    pheno <- readRDS(paste0("data/pheno/", .y, ".rds"))
    AUCBoot_no_NA(.x, pheno[ind.test])
  }), pred = NULL) %>%
  bind_rows(readRDS("results/all-ldpred1.rds")) %>%
  unnest_wider("AUCBoot") %>%
  mutate(Method = sub("^best_", "", Method),
         Method = stringr::str_replace_all(Method, "_", "-"),
         Method = sub("ldpred", "LDpred", Method)) %>%
  print()

library(ggplot2)
# http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/#a-colorblind-friendly-palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
               "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


## LDpred1 vs LDpred2

ggplot({
  auc_real %>%
    mutate(Method = factor(Method, levels = c(
      "LDpred1-inf", "LDpred2-inf-gwide", "LDpred1-grid",
      "LDpred2-grid-nosp-gwide", "LDpred2-grid-sp-gwide", "LDpred2-auto-gwide"))) %>%
    na.omit()
}, aes(Trait, Mean, fill = Method, color = Method)) +
  bigstatsr::theme_bigstatsr() +
  geom_hline(yintercept = 0.5, linetype = 2) +
  geom_col(position = position_dodge(), alpha = 0.6, color = "black", size = 1) +
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`),
                position = position_dodge(width = 0.9),
                color = "black", width = 0.2, size = 1) +
  scale_y_continuous(limits = c(0.5, 0.82), minor_breaks = 0:50 / 50,
                     oob = scales::rescale_none) +
  labs(y = "AUC") +
  theme(legend.position = c(0.2, 0.78)) +
  scale_fill_manual(values = cbPalette[c(6, 7, 4, 2, 5, 8)])

# ggsave("figures/AUC-real.pdf", width = 11.5, height = 6.5)


## LDpred2 per chromosome vs genome-wide

ggplot({
  auc_real %>%
    mutate(Method = factor(Method, levels = c(
      "LDpred2-inf-perchr", "LDpred2-inf-gwide",
      "LDpred2-grid-nosp-perchr", "LDpred2-grid-nosp-gwide",
      "LDpred2-auto-perchr", "LDpred2-auto-gwide"))) %>%
    na.omit()
}, aes(Trait, Mean, fill = Method, color = Method)) +
  bigstatsr::theme_bigstatsr() +
  geom_hline(yintercept = 0.5, linetype = 2) +
  geom_col(position = position_dodge(), alpha = 0.6, color = "black", size = 1) +
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`),
                position = position_dodge(width = 0.9),
                color = "black", width = 0.2, size = 1) +
  scale_y_continuous(limits = c(0.5, 0.82), minor_breaks = 0:50 / 50,
                     oob = scales::rescale_none) +
  labs(y = "AUC") +
  theme(legend.position = c(0.16, 0.76)) +
  scale_fill_manual(values = cbPalette[c(2, 5, 6, 3, 4, 8)])

# ggsave("figures/AUC-real-ldpred2.pdf", width = 12, height = 6)


## Tables of all results

auc_real %>%
  mutate_at(3:5, ~ round(. * 100, 1)) %>%
  mutate(AUC = sprintf("%.1f [%.1f-%.1f]", Mean, `2.5%`, `97.5%`)) %>%
  select(-(3:6)) %>%
  pivot_wider(names_from = "Trait", values_from = "AUC") %>%
  print(width = Inf) %>%
  xtable::xtable(align = "|l|c|c|c|c|c|c|c|c|c|") %>%
  print(include.rownames = FALSE)

#   Method                   Asthma           BRCA             CAD
#   <chr>                    <chr>            <chr>            <chr>
# 1 LDpred2-grid-nosp-perchr 58.8 [58.5-59.1] 64.6 [64.1-65.1] 61.5 [61.0-62.0]
# 2 LDpred2-grid-sp-perchr   58.5 [58.2-58.8] 64.6 [64.1-65.0] 61.2 [60.8-61.7]
# 3 LDpred2-auto-perchr      58.2 [57.9-58.5] 65.6 [65.2-66.1] 62.1 [61.6-62.6]
# 4 LDpred2-inf-perchr       57.1 [56.8-57.4] 61.9 [61.4-62.4] 60.8 [60.4-61.3]
# 5 LDpred2-inf-gwide        57.2 [56.9-57.4] 61.6 [61.2-62.1] 61.6 [61.1-62.1]
# 6 LDpred2-grid-nosp-gwide  59.3 [59.0-59.6] 65.5 [65.0-66.0] 63.6 [63.2-64.1]
# 7 LDpred2-grid-sp-gwide    59.4 [59.1-59.7] 65.7 [65.3-66.2] 63.5 [63.1-64.0]
# 8 LDpred2-auto-gwide       58.4 [58.2-58.7] 65.6 [65.2-66.1] 61.8 [61.3-62.3]
# 9 LDpred1-grid             58.6 [58.3-58.8] 58.9 [58.4-59.4] 60.3 [59.9-60.8]
# 10 LDpred1-inf             56.9 [56.6-57.1] 58.8 [58.3-59.3] 59.3 [58.9-59.8]
#   MDD              PRCA             RA               T1D              T2D
#   <chr>            <chr>            <chr>            <chr>            <chr>
# 1 57.9 [57.5-58.2] 68.2 [67.6-68.7] 59.9 [59.2-60.6] 77.0 [75.2-78.8] 63.3 [62.9-63.7]
# 2 58.2 [57.8-58.6] 67.8 [67.3-68.4] 60.0 [59.3-60.7] 75.9 [74.0-77.7] 63.2 [62.7-63.6]
# 3 58.9 [58.5-59.2] 70.2 [69.6-70.7] 59.7 [58.9-60.4] 76.6 [74.7-78.5] 63.9 [63.5-64.3]
# 4 58.9 [58.5-59.2] 65.0 [64.4-65.6] 59.6 [58.9-60.3] 74.5 [72.6-76.3] 61.5 [61.1-61.9]
# 5 59.0 [58.6-59.3] 64.4 [63.8-65.0] 59.3 [58.6-60.0] 71.0 [69.1-72.8] 61.6 [61.1-62.0]
# 6 59.1 [58.7-59.4] 70.2 [69.6-70.7] 60.3 [59.6-61.0] 78.4 [76.7-80.2] 64.3 [63.8-64.7]
# 7 59.0 [58.7-59.4] 70.2 [69.6-70.7] 60.5 [59.8-61.3] 78.3 [76.5-80.0] 64.0 [63.6-64.5]
# 8 59.0 [58.6-59.4] 70.1 [69.5-70.7] 59.7 [59.0-60.4] 77.7 [75.9-79.4] 63.8 [63.3-64.2]
# 9 56.4 [56.1-56.8] 63.5 [62.9-64.1] 59.1 [58.4-59.8] 57.4 [55.5-59.2] 50.3 [49.8-50.7]
# 10 56.5 [56.1-56.8] 61.6 [61.0-62.3] 59.5 [58.8-60.2] 57.7 [55.8-59.5] 51.8 [51.4-52.2]
