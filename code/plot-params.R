library(dplyr)
library(purrr)

## LDpred2-gwide params
files <- list.files("results/ldpred2-gwide", full.names = TRUE)
res_ldpred2_gwide <- tibble(res = lapply(files, readRDS)) %>%
  tidyr::unnest_wider("res")

res_ldpred2_gwide %>%
  transmute(
    Trait = sub("\\.rds$", "", basename(files)),
    h2_ldsc_obs = map_dbl(params, ~ median(.x$h2)),
    sparsity_grid = map_dbl(params, ~ {
      .x %>%
        filter(sparse) %>%
        arrange(desc(score)) %>%
        slice(1) %>%
        pull("sparsity")
    }),
    p_auto  = map_dbl(auto, ~ median(map_dbl(.x, "p_est"))),
    h2_auto_obs = map_dbl(auto, ~ median(map_dbl(.x, "h2_est"))),
    prev_pop = map_dbl(Trait, ~ {
      mean(readRDS(paste0("data/pheno/", .x, ".rds")), na.rm = TRUE)
    }),
    h2_auto_liab = map2_dbl(h2_auto_obs, prev_pop, ~ .x * bigsnpr::coef_to_liab(.y))
  ) %>%
  print() %>%
  xtable::xtable(align = "|l|c|c|c|c|c|c|c|") %>%
  print(include.rownames = FALSE)
#   Trait  h2_ldsc_obs sparsity_grid   p_auto h2_auto_obs prev_pop h2_auto_liab
#   <chr>        <dbl>         <dbl>    <dbl>       <dbl>    <dbl>        <dbl>
# 1 Asthma      0.104          0.544 0.000738      0.0534  0.150         0.0640
# 2 BRCA        0.127          0.414 0.00474       0.141   0.0758        0.136
# 3 CAD         0.0819         0.441 0.000978      0.0428  0.0619        0.0387
# 4 MDD         0.0923         0.377 0.0675        0.0877  0.0893        0.0890
# 5 PRCA        0.188          0.516 0.00181       0.189   0.0546        0.164
# 6 RA          0.416          0.193 0.0774        0.562   0.0286        0.406
# 7 T1D         0.955          0.498 0.00482       1.40    0.00244       0.575
# 8 T2D         0.148          0.498 0.00310       0.131   0.0512        0.112


## Look at CAD more specifically
res <- readRDS("results/ldpred2-gwide/CAD.rds")
library(ggplot2)
theme_set(bigstatsr::theme_bigstatsr())
ggplot(res$params, aes(x = p, y = score, color = as.factor(h2))) +
  geom_point() +
  geom_line() +
  scale_x_log10(breaks = 10^(-5:0), minor_breaks = res$params$p) +
  facet_wrap(~ sparse, labeller = label_both) +
  labs(y = "GLM Z-Score", color = "h2") +
  theme(legend.position = "top", panel.spacing = unit(1, "lines"))

# ggsave("figures/CAD-grid.pdf", width = 11, height = 5.5)

length(multi_auto <- res$auto) # 23 (out of 30)
map_dbl(multi_auto, "p_est")
#  [1] 0.0009807511 0.0009516860 0.0009784240 0.0010072454 0.0009774177 0.0009749286
#  [7] 0.0009947648 0.0010186770 0.0009347747 0.0009413393 0.0009877178 0.0009754323
# [13] 0.0010086016 0.0010132343 0.0010032294 0.0009281907 0.0009426621 0.0009643458
# [19] 0.0009659263 0.0010186284 0.0009637553 0.0010112630 0.0010289278
map_dbl(multi_auto, "h2_est")
#  [1] 0.04280892 0.04244302 0.04279115 0.04322751 0.04270675 0.04264875 0.04308405
#  [8] 0.04325921 0.04226913 0.04209307 0.04304932 0.04280672 0.04314088 0.04322056
# [15] 0.04309024 0.04205662 0.04241823 0.04262959 0.04265346 0.04323514 0.04273634
# [22] 0.04315293 0.04335818

auto <- multi_auto[[1]]
auto <- multi_auto[[23]]
bigstatsr::plot_grid(
  qplot(y = auto$path_p_est) +
    geom_hline(yintercept = auto$p_est, col = "blue") +
    scale_y_log10() +
    labs(y = "p"),
  qplot(y = auto$path_h2_est) +
    geom_hline(yintercept = auto$h2_est, col = "blue") +
    labs(y = "h2"),
  ncol = 1, align = "hv"
)

# ggsave("figures/CAD-auto.pdf", width = 8, height = 5)
# ggsave("figures/CAD-auto2.pdf", width = 8, height = 5)


## Look at GWAS sumstats
sumstats <- readRDS("data/sumstats/CAD.rds")
fun.pred <- function(xtr) {
  stats::pchisq(xtr, df = 1, lower.tail = FALSE, log.p = TRUE) / log(10)
}
gwas <- structure(data.frame(score = with(sumstats, (beta / beta_se)^2)),
                  class = c("mhtest", "data.frame"),
                  transfo = identity,
                  predict = fun.pred)
bigsnpr::snp_manhattan(gwas, sumstats$chr, sumstats$pos, npoint = 50e3)
