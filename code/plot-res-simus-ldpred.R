library(bigsnpr)
ukb <- snp_attach("data/UKBB_imp_HM3.rds")
G <- ukb$genotypes

set.seed(1)
ind.val <- sort(sample(nrow(G), 10e3))
ind.gwas <- sort(sample(setdiff(rows_along(G), ind.val), 300e3))
ind.test <- setdiff(rows_along(G), c(ind.gwas, ind.val))

library(tidyverse)
library(foreach)
pred_ldpred2 <- tibble(basename = list.files("data/pheno-simu")) %>%
  mutate(
    pheno_file = file.path("data/pheno-simu", basename),
    res_file = map(file.path("results/simu-ldpred2", basename), ~ {
      foreach(chr = 1:22, .combine = 'c') %do% {
        sub("\\.rds$", paste0("_chr", chr, ".rds"), .x)
      }
    }))

res_ldpred2 <- pred_ldpred2 %>%
  mutate(basename = gsubfn::gsubfn("(.+)_([0-9]+)\\.rds$", "\\1 \\2", basename)) %>%
  separate("basename", c("params", "iter"), sep = " ", convert = TRUE) %>%
  unnest_longer(col = "res_file", indices_to = "chr") %>%
  filter(file.exists(res_file)) %>%
  mutate(pred = map(res_file, readRDS), res_file = NULL) %>%
  unnest_longer(col = "pred", indices_to = "Method") %>%
  group_by_at(-(4:5)) %>%
  summarise(pred = list(do.call(plus, pred))) %>%
  ungroup() %>%
  mutate(Method = sub("beta", "ldpred2", Method),
         pheno = map(pheno_file, readRDS), pheno_file = NULL,
         AUC = map2_dbl(pred, pheno, ~ AUC(.x, .y[ind.test])),
         pred = NULL, pheno = NULL) %>%
  print()

auc_simu <- res_ldpred2 %>%
  bind_rows(readRDS("results/all-simu-ldpred1.rds")) %>%
  group_by(params, Method) %>%
  summarise(aucs = {
    boot <- replicate(1e4, mean(sample(AUC, replace = TRUE)))
    q <- quantile(boot, c(0.025, 0.975))
    list(c(mean = mean(boot), inf = q[[1]], sup = q[[2]]))
  }) %>%
  ungroup() %>%
  unnest_wider("aucs") %>%
  print()

cols <- rep(NA, 6)
cols[c(1, 3)] <- viridisLite::inferno(10)[c(5, 7)]
cols[-c(1, 3)] <- viridisLite::viridis(8)[c(3:5, 7)]

auc_simu %>%
  mutate(Method = factor(Method, levels = c(
    "ldpred1_inf", "ldpred2_inf", "best_ldpred1_grid", "best_ldpred2_grid_nosp",
    "best_ldpred2_grid_sp", "ldpred2_auto"))) %>%
  ggplot(aes(params, mean, fill = Method, color = Method)) +
  bigstatsr::theme_bigstatsr() +
  geom_hline(yintercept = 0.5, linetype = 2) +
  geom_col(position = position_dodge(), alpha = 0.6, color = "black", size = 1) +
  geom_errorbar(aes(ymin = inf, ymax = sup),
                position = position_dodge(width = 0.9), color = "black", width = 0.2, size = 1) +
  scale_y_continuous(limits = c(0.5, 0.9), minor_breaks = 0:50 / 50,
                     oob = scales::rescale_none) +
  labs(x = "Simulation", y = "AUC") +
  theme(legend.position = c(0.51, 0.8)) +
  scale_fill_manual(values = cols)

# ggsave("figures/AUC-simu.pdf", width = 11, height = 6.5)

auc_simu %>%
  mutate_at(3:5, ~ round(. * 100, 1)) %>%
  mutate(AUC = sprintf("%.1f [%.1f-%.1f]", mean, inf, sup)) %>%
  select(-(3:5)) %>%
  pivot_wider(names_from = "Method", values_from = "AUC") %>%
  print(width = Inf) %>%
  xtable::xtable(align = "|l|c|c|c|c|c|c|c|") %>%
  print(include.rownames = FALSE)

#   params           best_ldpred1_grid best_ldpred2_grid_nosp
#   <chr>            <chr>             <chr>
# 1 all_40_300_15    73.6 [72.0-75.4]  81.7 [81.6-81.9]
# 2 all_40_3000_15   70.4 [69.8-70.9]  78.9 [78.7-79.0]
# 3 all_40_30000_15  69.1 [68.8-69.3]  70.6 [70.4-70.8]
# 4 all_40_300000_15 69.2 [69.0-69.4]  69.3 [69.0-69.6]
# 5 HLA_30_300_15    54.7 [53.9-55.6]  76.3 [75.7-77.0]
# 6 HLA_30_3000_15   54.5 [53.8-55.1]  76.4 [76.1-76.6]
#   best_ldpred2_grid_sp ldpred1_inf      ldpred2_auto     ldpred2_inf
#   <chr>                <chr>            <chr>            <chr>
# 1 81.7 [81.6-81.9]     68.9 [68.4-69.3] 82.0 [81.8-82.2] 70.3 [69.9-70.6]
# 2 78.8 [78.7-78.9]     69.3 [69.0-69.6] 79.3 [79.2-79.4] 70.0 [69.7-70.2]
# 3 70.5 [70.3-70.7]     69.5 [69.2-69.7] 70.8 [70.6-71.0] 69.9 [69.8-70.1]
# 4 69.4 [69.1-69.6]     69.5 [69.3-69.6] 69.8 [69.6-70.0] 69.9 [69.7-70.0]
# 5 76.5 [75.8-77.1]     54.5 [53.0-55.9] 73.1 [72.0-74.1] 75.6 [74.8-76.3]
# 6 76.5 [76.3-76.7]     54.5 [53.6-55.5] 72.9 [71.6-74.0] 75.7 [75.4-75.9]

# % latex table generated in R 3.5.1 by xtable 1.8-4 package
# % Fri May 29 16:53:24 2020
# \begin{table}[ht]
# \centering
# \begin{tabular}{|c|c|c|c|c|c|c|}
# \hline
# params & best\_ldpred1\_grid & best\_ldpred2\_grid\_nosp & best\_ldpred2\_grid\_sp & ldpred1\_inf & ldpred2\_auto & ldpred2\_inf \\
# \hline
# all\_40\_300\_15 & 73.6 [72.0-75.4] & 81.7 [81.6-81.9] & 81.7 [81.6-81.9] & 68.9 [68.4-69.3] & 82.0 [81.8-82.2] & 70.3 [69.9-70.6] \\
# all\_40\_3000\_15 & 70.4 [69.8-70.9] & 78.9 [78.7-79.0] & 78.8 [78.7-78.9] & 69.3 [69.0-69.6] & 79.3 [79.2-79.4] & 70.0 [69.7-70.2] \\
# all\_40\_30000\_15 & 69.1 [68.8-69.3] & 70.6 [70.4-70.8] & 70.5 [70.3-70.7] & 69.5 [69.2-69.7] & 70.8 [70.6-71.0] & 69.9 [69.8-70.1] \\
# all\_40\_300000\_15 & 69.2 [69.0-69.4] & 69.3 [69.0-69.6] & 69.4 [69.1-69.6] & 69.5 [69.3-69.6] & 69.8 [69.6-70.0] & 69.9 [69.7-70.0] \\
# HLA\_30\_300\_15 & 54.7 [53.9-55.6] & 76.3 [75.7-77.0] & 76.5 [75.8-77.1] & 54.5 [53.0-55.9] & 73.1 [72.0-74.1] & 75.6 [74.8-76.3] \\
# HLA\_30\_3000\_15 & 54.5 [53.8-55.1] & 76.4 [76.1-76.6] & 76.5 [76.3-76.7] & 54.5 [53.6-55.5] & 72.9 [71.6-74.0] & 75.7 [75.4-75.9] \\
# \hline
# \end{tabular}
# \end{table}
