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
  print()

AUCBoot_no_NA <- function(pred, target) {
  is_not_na <- which(!is.na(target))
  AUCBoot(pred[is_not_na], target[is_not_na])
}

res_ldpred2 <- pred_ldpred2 %>%
  group_by_at(-(2:3)) %>%
  summarise(pred = list(do.call(plus, map(pred, replace_na, replace = 0)))) %>%
  ungroup() %>%
  mutate(Method = sub("beta", "ldpred2", Method),
         AUCBoot = furrr::future_map2(pred, Trait, ~ {
           pheno <- readRDS(paste0("data/pheno/", .y, ".rds"))
           AUCBoot_no_NA(.x, pheno[ind.test])
         }),
         pred = NULL) %>%
  print()

auc_real <- res_ldpred2 %>%
  filter(Method != "best_ldpred2") %>%
  bind_rows(readRDS("results/all-ldpred1.rds")) %>%
  mutate(Method = factor(Method, levels = c(
    "ldpred1_inf", "ldpred2_inf", "best_ldpred1_grid", "best_ldpred2_grid_nosp",
    "best_ldpred2_grid_sp", "ldpred2_auto"))) %>%
  unnest_wider("AUCBoot") %>%
  print()

cols <- rep(NA, 6)
cols[c(1, 3)] <- viridisLite::inferno(10)[c(5, 7)]
cols[-c(1, 3)] <- viridisLite::viridis(8)[c(3:5, 7)]

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
  theme(legend.position = c(0.15, 0.79)) +
  scale_fill_manual(values = cols)

# ggsave("figures/AUC-real.pdf", width = 11, height = 6.5)

auc_real %>%
  mutate_at(3:5, ~ round(. * 100, 1)) %>%
  mutate(AUC = sprintf("%.1f [%.1f-%.1f]", Mean, `2.5%`, `97.5%`)) %>%
  select(-(3:6)) %>%
  pivot_wider(names_from = "Method", values_from = "AUC") %>%
  print(width = Inf) %>%
  xtable::xtable(align = "|l|c|c|c|c|c|c|c|") %>%
  print(include.rownames = FALSE)

#   Trait  best_ldpred2_grid_nosp best_ldpred2_grid_sp ldpred2_auto
#   <chr>  <chr>                  <chr>                <chr>
# 1 Asthma 58.8 [58.5-59.1]       58.5 [58.2-58.8]     58.2 [57.9-58.5]
# 2 BRCA   64.6 [64.1-65.1]       64.5 [64.1-65.0]     65.6 [65.2-66.1]
# 3 CAD    61.5 [61.0-62.0]       61.2 [60.8-61.7]     62.1 [61.6-62.6]
# 4 MDD    57.9 [57.5-58.2]       58.2 [57.8-58.6]     58.9 [58.5-59.2]
# 5 PRCA   68.2 [67.6-68.7]       67.8 [67.2-68.4]     70.2 [69.6-70.8]
# 6 RA     59.9 [59.2-60.6]       60.0 [59.3-60.7]     59.7 [59.0-60.4]
# 7 T1D    77.0 [75.2-78.8]       75.9 [74.0-77.7]     76.6 [74.8-78.5]
# 8 T2D    63.3 [62.9-63.7]       63.2 [62.8-63.6]     63.9 [63.5-64.3]
#   ldpred2_inf      best_ldpred1_grid ldpred1_inf
#   <chr>            <chr>             <chr>
# 1 57.1 [56.8-57.4] 58.6 [58.3-58.8]  56.9 [56.6-57.1]
# 2 61.9 [61.4-62.4] 58.9 [58.4-59.4]  58.8 [58.3-59.3]
# 3 60.8 [60.4-61.3] 60.3 [59.9-60.8]  59.3 [58.9-59.8]
# 4 58.9 [58.5-59.2] 56.4 [56.1-56.8]  56.5 [56.1-56.8]
# 5 65.1 [64.5-65.7] 63.5 [62.9-64.1]  61.6 [61.0-62.3]
# 6 59.6 [58.8-60.3] 59.1 [58.4-59.8]  59.5 [58.8-60.2]
# 7 74.5 [72.6-76.3] 57.4 [55.5-59.2]  57.7 [55.8-59.5]
# 8 61.5 [61.1-61.9] 50.3 [49.8-50.7]  51.8 [51.4-52.2]

# % latex table generated in R 3.5.1 by xtable 1.8-4 package
# % Thu May 28 15:38:13 2020
# \begin{table}[ht]
# \centering
# \begin{tabular}{|c|c|c|c|c|c|c|}
# \hline
# Trait & best\_ldpred2\_grid\_nosp & best\_ldpred2\_grid\_sp & ldpred2\_auto & ldpred2\_inf & best\_ldpred1\_grid & ldpred1\_inf \\
# \hline
# Asthma & 58.8 [58.5-59.1] & 58.5 [58.2-58.8] & 58.2 [57.9-58.5] & 57.1 [56.8-57.4] & 58.6 [58.3-58.8] & 56.9 [56.6-57.1] \\
# BRCA & 64.6 [64.1-65.1] & 64.5 [64.1-65.0] & 65.6 [65.2-66.1] & 61.9 [61.4-62.4] & 58.9 [58.4-59.4] & 58.8 [58.3-59.3] \\
# CAD & 61.5 [61.0-62.0] & 61.2 [60.8-61.7] & 62.1 [61.6-62.6] & 60.8 [60.4-61.3] & 60.3 [59.9-60.8] & 59.3 [58.9-59.8] \\
# MDD & 57.9 [57.5-58.2] & 58.2 [57.8-58.6] & 58.9 [58.5-59.2] & 58.9 [58.5-59.2] & 56.4 [56.1-56.8] & 56.5 [56.1-56.8] \\
# PRCA & 68.2 [67.6-68.7] & 67.8 [67.2-68.4] & 70.2 [69.6-70.8] & 65.1 [64.5-65.7] & 63.5 [62.9-64.1] & 61.6 [61.0-62.3] \\
# RA & 59.9 [59.2-60.6] & 60.0 [59.3-60.7] & 59.7 [59.0-60.4] & 59.6 [58.8-60.3] & 59.1 [58.4-59.8] & 59.5 [58.8-60.2] \\
# T1D & 77.0 [75.2-78.8] & 75.9 [74.0-77.7] & 76.6 [74.8-78.5] & 74.5 [72.6-76.3] & 57.4 [55.5-59.2] & 57.7 [55.8-59.5] \\
# T2D & 63.3 [62.9-63.7] & 63.2 [62.8-63.6] & 63.9 [63.5-64.3] & 61.5 [61.1-61.9] & 50.3 [49.8-50.7] & 51.8 [51.4-52.2] \\
# \hline
# \end{tabular}
# \end{table}
