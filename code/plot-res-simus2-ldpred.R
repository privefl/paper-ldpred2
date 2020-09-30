library(bigsnpr)
ukb <- snp_attach("data/UKBB_imp_HM3.rds")
G <- ukb$genotypes

load("data/ind_gwas_val_test.RData")

library(dplyr)
library(tidyr)

pred_ldpred2 <-
  tibble(basename = list.files("results/simu-ldpred2", "all_40_3000_15")) %>%
  mutate(
    pheno_file = file.path("data/pheno-simu", sub("_chr[0-9]+.*", ".rds", basename)),
    res_file = file.path("results/simu-ldpred2", basename)) %>%
  print()

res_ldpred2 <- pred_ldpred2 %>%
  mutate(N = sub(".*_chr[0-9]+(.*)\\.rds$", "\\1", basename),
         N = ifelse(N == "", 300000L, as.integer(sub("_", "", N))),
         chr = as.integer(sub(".*_chr([0-9]+).*\\.rds$", "\\1", basename)),
         basename = NULL) %>%
  mutate(pred = lapply(res_file, readRDS), res_file = NULL) %>%
  unnest_longer(col = "pred", indices_to = "Method") %>%
  group_by(pheno_file, N, Method) %>%
  summarise(pred = list(do.call(plus, pred))) %>%
  ungroup() %>%
  mutate(Method = paste0(sub("beta", "ldpred2", Method), "_perchr"),
         pheno = lapply(pheno_file, readRDS), pheno_file = NULL,
         AUC = purrr::map2_dbl(pred, pheno, ~ AUC(.x, .y[ind.test])),
         pred = NULL, pheno = NULL) %>%
  print()

res_ldpred2_gwide <-
  tibble(res_file = list.files("results/simu-ldpred2-gwide/", full.names = TRUE),
         pheno_file = file.path("data/pheno-simu",
                                sub("_[0-9]+\\.rds", ".rds", basename(res_file))),
         N = as.integer(sub(".+_([0-9]+)\\.rds", "\\1", res_file))) %>%
  mutate(pred = lapply(res_file, function(file) as.list(readRDS(file)$pred)),
         res_file = NULL) %>%
  unnest_longer(col = "pred", indices_to = "Method") %>%
  mutate(Method = paste0(sub("beta", "ldpred2", Method), "_gwide"),
         pheno = lapply(pheno_file, readRDS), pheno_file = NULL,
         AUC = purrr::map2_dbl(pred, pheno, ~ AUC(.x, .y[ind.test])),
         pred = NULL, pheno = NULL) %>%
  print()

res_ldpred1 <- readRDS("results/all-simu-ldpred1.rds") %>%
  filter(grepl("all_40_3000_15", params)) %>%
  mutate(N = sub("all_40_3000_15(.*)$", "\\1", params),
         N = ifelse(N == "", 300000L, as.integer(sub("_", "", N))),
         params = NULL, iter = NULL)

auc_simu <-
  bind_rows(res_ldpred2, res_ldpred2_gwide, res_ldpred1) %>%
  mutate(Method = sub("^best_", "", Method),
         Method = stringr::str_replace_all(Method, "_", "-"),
         Method = sub("ldpred", "LDpred", Method)) %>%
  group_by(N, Method) %>%
  summarise(aucs = {
    boot <- replicate(1e4, mean(sample(AUC, replace = TRUE)))
    q <- quantile(boot, c(0.025, 0.975))
    list(c(mean = mean(boot), inf = q[[1]], sup = q[[2]]))
  }) %>%
  ungroup() %>%
  unnest_wider("aucs") %>%
  print()


library(ggplot2)
# http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/#a-colorblind-friendly-palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
               "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


## Grid models

ggplot({
  auc_simu %>%
    mutate(Method = factor(Method, levels = c(
      "LDpred1-grid", "LDpred2-grid-nosp-perchr", "LDpred2-grid-sp-perchr",
      "LDpred2-grid-nosp-gwide", "LDpred2-grid-sp-gwide"))) %>%
    na.omit() %>%
    mutate(N = N / 1000)
}, aes(as.factor(N), mean, fill = Method, color = Method)) +
  bigstatsr::theme_bigstatsr() +
  geom_hline(yintercept = 0.5, linetype = 2) +
  geom_col(position = position_dodge(), alpha = 0.6, color = "black", size = 1) +
  geom_errorbar(aes(ymin = inf, ymax = sup),
                position = position_dodge(width = 0.9), color = "black", width = 0.2, size = 1) +
  scale_y_continuous(limits = c(0.5, 0.82), minor_breaks = 0:50 / 50,
                     oob = scales::rescale_none) +
  labs(x = "Sample size (x1000)", y = "AUC") +
  scale_fill_manual(values = cbPalette[c(1, 2, 5, 6, 3)]) +
  theme(legend.position = c(0.23, 0.78))

# ggsave("figures/AUC-simu2-grid.pdf", width = 9, height = 5.5)


## LDpred2 models

ggplot({
  auc_simu %>%
    mutate(Method = factor(Method, levels = c(
      "LDpred2-inf-perchr", "LDpred2-inf-gwide",
      "LDpred2-grid-nosp-perchr", "LDpred2-grid-nosp-gwide",
      "LDpred2-auto-perchr", "LDpred2-auto-gwide"))) %>%
    na.omit() %>%
    mutate(N = N / 1000)
}, aes(as.factor(N), mean, fill = Method, color = Method)) +
  bigstatsr::theme_bigstatsr() +
  geom_hline(yintercept = 0.5, linetype = 2) +
  geom_col(position = position_dodge(), alpha = 0.6, color = "black", size = 1) +
  geom_errorbar(aes(ymin = inf, ymax = sup),
                position = position_dodge(width = 0.9), color = "black", width = 0.2, size = 1) +
  scale_y_continuous(limits = c(0.5, 0.82), minor_breaks = 0:50 / 50,
                     oob = scales::rescale_none) +
  labs(x = "Sample size (x1000)", y = "AUC") +
  scale_fill_manual(values = cbPalette[c(2, 5, 6, 3, 4, 8)]) +
  theme(legend.position = c(0.22, 0.75))

# ggsave("figures/AUC-simu2-ldpred2.pdf", width = 10, height = 5.5)


## Tables of all results

auc_simu %>%
  mutate_at(3:5, ~ round(. * 100, 1)) %>%
  mutate(AUC = sprintf("%.1f [%.1f-%.1f]", mean, inf, sup)) %>%
  select(-(3:5)) %>%
  pivot_wider(names_from = "N", values_from = "AUC") %>%
  print(width = Inf) %>%
  xtable::xtable(align = "|l|c|c|c|c|c|c|") %>%
  print(include.rownames = FALSE)

#   Method                   `10000`          `20000`          `50000`
#   <chr>                    <chr>            <chr>            <chr>
# 1 LDpred1-grid             56.1 [55.8-56.4] 59.3 [58.9-59.8] 66.3 [65.2-67.1]
# 2 LDpred1-inf              55.7 [55.4-56.0] 57.7 [57.4-58.0] 61.2 [60.9-61.5]
# 3 LDpred2-auto-gwide       55.9 [55.5-56.2] 59.3 [58.9-59.8] 67.5 [67.3-67.8]
# 4 LDpred2-auto-perchr      54.1 [53.7-54.5] 57.3 [56.8-57.8] 66.6 [66.3-66.9]
# 5 LDpred2-grid-nosp-gwide  55.7 [55.2-56.0] 59.6 [59.3-60.0] 67.4 [67.2-67.7]
# 6 LDpred2-grid-nosp-perchr 54.4 [54.0-54.8] 58.6 [58.1-59.1] 66.8 [66.5-67.1]
# 7 LDpred2-grid-sp-gwide    56.1 [55.7-56.5] 59.6 [59.3-60.0] 67.4 [67.2-67.7]
# 8 LDpred2-grid-sp-perchr   54.5 [54.2-54.9] 58.4 [57.9-58.9] 66.8 [66.6-67.1]
# 9 LDpred2-inf-gwide        55.5 [55.2-55.8] 57.7 [57.5-58.0] 61.4 [61.1-61.6]
# 10 LDpred2-inf-perchr       54.6 [54.3-55.0] 57.2 [56.8-57.5] 61.3 [61.0-61.5]
#   `120000`         `300000`
#   <chr>            <chr>
# 1 70.7 [69.2-71.7] 70.4 [69.8-70.9]
# 2 65.0 [64.7-65.3] 69.3 [69.0-69.7]
# 3 74.7 [74.6-74.9] 79.3 [79.2-79.4]
# 4 74.7 [74.5-74.8] 79.3 [79.2-79.4]
# 5 74.6 [74.5-74.8] 79.1 [79.0-79.3]
# 6 74.2 [73.9-74.5] 78.9 [78.7-79.0]
# 7 74.6 [74.4-74.8] 79.2 [79.1-79.3]
# 8 74.2 [73.9-74.4] 78.8 [78.7-78.9]
# 9 65.4 [65.1-65.6] 69.9 [69.7-70.1]
# 10 65.4 [65.2-65.6] 70.0 [69.7-70.2]
