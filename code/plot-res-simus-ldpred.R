library(bigsnpr)
ukb <- snp_attach("data/UKBB_imp_HM3.rds")
G <- ukb$genotypes

load("data/ind_gwas_val_test.RData")

library(dplyr)
library(tidyr)

pred_ldpred2 <- tibble(basename = list.files("data/pheno-simu")) %>%
  mutate(
    pheno_file = file.path("data/pheno-simu", basename),
    res_file = lapply(file.path("results/simu-ldpred2", basename), function(.x) {
      sapply(1:22, function(chr) sub("\\.rds$", paste0("_chr", chr, ".rds"), .x))
    })) %>%
  print()

res_ldpred2 <- pred_ldpred2 %>%
  mutate(basename = gsubfn::gsubfn("(.+)_([0-9]+)\\.rds$", "\\1 \\2", basename)) %>%
  separate("basename", c("params", "iter"), sep = " ", convert = TRUE) %>%
  unnest_longer(col = "res_file", indices_to = "chr") %>%
  mutate(pred = lapply(res_file, readRDS), res_file = NULL) %>%
  unnest_longer(col = "pred", indices_to = "Method") %>%
  group_by_at(-(4:5)) %>%
  summarise(pred = list(do.call(plus, pred))) %>%
  ungroup() %>%
  mutate(Method = paste0(sub("beta", "ldpred2", Method), "_perchr"),
         pheno = lapply(pheno_file, readRDS), pheno_file = NULL,
         AUC = purrr::map2_dbl(pred, pheno, ~ AUC(.x, .y[ind.test])),
         pred = NULL, pheno = NULL) %>%
  print()

res_ldpred1 <- filter(readRDS("results/all-simu-ldpred1.rds"),
                      params %in% res_ldpred2$params)

auc_simu <-
  rbind(res_ldpred2, res_ldpred1) %>%
  mutate(Method = sub("^best_", "", Method),
         Method = stringr::str_replace_all(Method, "_", "-"),
         Method = sub("ldpred", "LDpred", Method)) %>%
  group_by(params, Method) %>%
  summarise(aucs = {
    boot <- replicate(1e4, mean(sample(AUC, replace = TRUE)))
    q <- quantile(boot, c(0.025, 0.975))
    list(c(mean = mean(boot), inf = q[[1]], sup = q[[2]]))
  }) %>%
  ungroup() %>%
  unnest_wider("aucs") %>%
  mutate(Method = factor(Method, ordered = TRUE, levels = c(
    "LDpred1-inf", "LDpred2-inf-perchr", "LDpred1-grid",
    "LDpred2-grid-nosp-perchr", "LDpred2-grid-sp-perchr", "LDpred2-auto-perchr")),
    params = sub("_15$", "", params)) %>%
  print()


library(ggplot2)
# http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/#a-colorblind-friendly-palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
               "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

ggplot(auc_simu, aes(params, mean, fill = Method, color = Method)) +
  bigstatsr::theme_bigstatsr() +
  geom_hline(yintercept = 0.5, linetype = 2) +
  geom_col(position = position_dodge(), alpha = 0.6, color = "black", size = 1) +
  geom_errorbar(aes(ymin = inf, ymax = sup),
                position = position_dodge(width = 0.9), color = "black", width = 0.2, size = 1) +
  scale_y_continuous(limits = c(0.5, 0.92), minor_breaks = 0:50 / 50,
                     oob = scales::rescale_none) +
  labs(x = "Simulation", y = "AUC") +
  theme(legend.position = c(0.45, 0.8)) +
  scale_fill_manual(values = cbPalette[c(6, 7, 4, 2, 5, 8)])

# ggsave("figures/AUC-simu.pdf", width = 12.5, height = 6.5)

auc_simu %>%
  mutate_at(3:5, ~ round(. * 100, 1)) %>%
  mutate(AUC = sprintf("%.1f [%.1f-%.1f]", mean, inf, sup)) %>%
  select(-(3:5)) %>%
  pivot_wider(names_from = "params", values_from = "AUC") %>%
  arrange(Method) %>%
  print(width = Inf) %>%
  xtable::xtable(align = "|l|c|c|c|c|c|c|c|c|") %>%
  print(include.rownames = FALSE)

#   Method                   all_40_300       all_40_3000      all_40_30000
#   <ord>                    <chr>            <chr>            <chr>
# 1 LDpred1-inf              68.9 [68.4-69.3] 69.3 [68.9-69.6] 69.5 [69.2-69.7]
# 2 LDpred2-inf-perchr       70.3 [69.9-70.6] 70.0 [69.7-70.2] 69.9 [69.7-70.1]
# 3 LDpred1-grid             73.6 [72.0-75.4] 70.4 [69.8-70.9] 69.1 [68.8-69.3]
# 4 LDpred2-grid-nosp-perchr 81.7 [81.6-81.9] 78.9 [78.7-79.0] 70.6 [70.4-70.8]
# 5 LDpred2-grid-sp-perchr   81.7 [81.6-81.9] 78.8 [78.7-78.9] 70.5 [70.3-70.7]
# 6 LDpred2-auto-perchr      82.0 [81.8-82.2] 79.3 [79.2-79.4] 70.8 [70.6-71.0]
#   all_40_300000    both_40          HLA_30_300       HLA_30_3000
#   <chr>            <chr>            <chr>            <chr>
# 1 69.5 [69.3-69.6] 62.9 [61.6-64.3] 54.4 [52.9-55.9] 54.5 [53.7-55.5]
# 2 69.9 [69.7-70.0] 74.2 [73.8-74.5] 75.6 [74.8-76.3] 75.7 [75.4-75.9]
# 3 69.2 [69.0-69.4] 61.7 [60.4-63.1] 54.7 [53.9-55.5] 54.5 [53.8-55.1]
# 4 69.3 [69.0-69.6] 72.8 [67.3-75.7] 76.4 [75.7-77.0] 76.4 [76.1-76.6]
# 5 69.4 [69.1-69.6] 72.8 [67.3-75.7] 76.5 [75.8-77.1] 76.5 [76.3-76.7]
# 6 69.8 [69.6-70.0] 73.7 [73.3-74.2] 73.1 [72.0-74.1] 72.9 [71.6-74.0]
