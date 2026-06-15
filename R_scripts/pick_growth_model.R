library(brms)
library(pcvr)
library(ggplot2)
setwd("~/maize-abiotic-stress-paper-2026/")

df <- read.csv("datasets/2024-11-13_VIS_SV_MG001_MG002_dead-plant-filtered_calcualted-values_D-DH-H-C.csv")

# Pick subset of genotypes to test model fits on

set.seed(456)
samp_genos <- sample(unique(df$genotype), 9)
sub <- df[df$genotype %in% samp_genos, ]

# Visual guess at starting priors

(p <- ggplot(sub, aes(x = day, y = area_cm2, group = barcode, color = treatment)) +
  facet_wrap(~genotype) +
  geom_line() +
  theme_bw()
)
p +
  geom_line(
    data = growthSim("power law", params = list("A" = 4, "B" = 1), t = max(sub$day)),
    aes(x = time, y = y, group = id), color = "black"
  ) +
  labs(title = "PL Check")

p +
  geom_line(
    data = growthSim("logistic", params = list("A" = 50, "B" = 8, "C" = 3), t = max(sub$day)),
    aes(x = time, y = y, group = id), color = "black"
  ) +
  labs(title = "Logistic Check")

p +
  geom_line(
    data = growthSim("gompertz", params = list("A" = 50, "B" = 8, "C" = 0.25), t = max(sub$day)),
    aes(x = time, y = y, group = id), color = "black"
  ) +
  labs(title = "Gompertz Check")

# Set up skeletons of models (avoids recompilation time)

if (!file.exists("datasets/picking_growth_model.rdata")) {

  pl_ss <- growthSS(model = "power law", form = area_cm2 ~ day | barcode / treatment,
    df = df[df$genotype == "G1", ], start = list("A" = 4, "B" = 1), type = "brms", sigma = "spline")
  pl_skel <- fitGrowth(pl_ss, iter = 2000, backend = "cmdstanr", cores = 4, chains = 4)

  lg_ss <- growthSS(model = "logistic", form = area_cm2 ~ day | barcode / treatment,
    df = df[df$genotype == "G1", ], start = list("A" = 50, "B" = 8, "C" = 3),
    type = "brms", sigma = "spline")
  lg_skel <- fitGrowth(lg_ss, iter = 2000, backend = "cmdstanr", cores = 4, chains = 4)

  gom_ss <- growthSS(model = "gompertz", form = area_cm2 ~ day | barcode / treatment,
    df = df[df$genotype == "G1", ], start = list("A" = 50, "B" = 8, "C" = 0.25),
    type = "brms", sigma = "spline")
  gom_skel <- fitGrowth(gom_ss, iter = 2000, backend = "cmdstanr", cores = 4, chains = 4)

  lin_ss <- growthSS(model = "linear", form = area_cm2 ~ day | barcode / treatment,
    df = df[df$genotype == "G1", ], start = list("A" = 10),
    type = "brms", sigma = "spline")
  lin_skel <- fitGrowth(lin_ss, iter = 2000, backend = "cmdstanr", cores = 4, chains = 4)

  # Run models for each of the sampled genotypes

  res_list <- lapply(samp_genos, function(geno){
    d <- df[df$genotype == geno, ]
    if (geno != samp_genos[1]) {
      pl_mod <- update(pl_skel, newdata = d, cores = 4, chains = 4, iter = 2000, backend = "cmdstanr")
      lg_mod <- update(lg_skel, newdata = d, cores = 4, chains = 4, iter = 2000, backend = "cmdstanr")
      gom_mod <- update(gom_skel, newdata = d, cores = 4, chains = 4, iter = 2000, backend = "cmdstanr")
      lin_mod <- update(lin_skel, newdata = d, cores = 4, chains = 4, iter = 2000, backend = "cmdstanr")
    } else {
      pl_mod <- pl_skel
      lg_mod <- lg_skel
      gom_mod <- gom_skel
      lin_mod <- lin_skel 
    }
    pl_mod <- brms::add_criterion(pl_mod, c("loo", "waic"))
    lg_mod <- brms::add_criterion(lg_mod, c("loo", "waic"))
    gom_mod <- brms::add_criterion(gom_mod, c("loo", "waic"))
    lin_mod <- brms::add_criterion(lin_mod, c("loo", "waic"))
    return(list("p" = pl_mod, "lg" = lg_mod, "g" = gom_mod, "ln" = lin_mod))
  })

  # Save result (this takes non-trivial minutes to run)

  save(res_list, file = "datasets/picking_growth_model.rdata")
} else {
  print(load("datasets/picking_growth_model.rdata"))
}

# Compare LOO and WAIC between models

res <- do.call(rbind, lapply(seq_along(res_list), function(i) {
  res <- res_list[[i]]
  names(res)[1] <- "x"
  geno <- samp_genos[i]
  lc <- as.data.frame(do.call(loo_compare,
    append(res, list("criterion" = "loo",
      "model_names" = c("power_law", "logistic", "gompertz", "linear")))))
  lc$mod <- rownames(lc)
  lc$preference <- seq_len(nrow(lc))
  lc$geno <- geno
  lc$criterion <- "LOO IC"
  #lc <- as.data.frame(tidyr::pivot_longer(lc[, c(7:9)], cols = c(1:2)))
  lcsub <- lc[, c(1, 7:12)]
  lcsub$sig_improvement <- abs(lcsub$elpd_diff) > (3 * lcsub$se_looic)
  colnames(lcsub)[2:3] <- c("val", "se_val")

  waicc <- as.data.frame(do.call(loo_compare,
    append(res, list("criterion" = "waic",
      "model_names" = c("power_law", "logistic", "gompertz", "linear")))))
  waicc$mod <- rownames(waicc)
  waicc$preference <- seq_len(nrow(waicc))
  waicc$geno <- geno
  waicc$criterion <- "WAIC"
  #lc <- as.data.frame(tidyr::pivot_longer(lc[, c(7:9)], cols = c(1:2)))
  waiccsub <- waicc[, c(1, 7:12)]
  waiccsub$sig_improvement <- abs(waiccsub$elpd_diff) > (3 * waiccsub$se_waic)
  colnames(waiccsub)[2:3] <- c("val", "se_val")

  return(rbind(lcsub, waiccsub))
}))

head(res)

ggplot(res[res$criterion == "LOO IC", ], aes(x = mod, y = val, fill = mod)) +
  facet_wrap(~geno) +
  geom_col() +
  geom_errorbar(aes(ymin = val - 3 * se_val, ymax = val + 3 * se_val),
    width = 0.25, linewidth = 0.5) +
  geom_point(data = res[res$criterion == "LOO IC" & res$sig_improvement, ]) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(title = "LOO IC",
    subtitle = "Significant differences marked with a point",
    x = "Model", y = "LOO IC"
  )
ggsave("supplement/looic.png", width = 8, height = 6, dpi = 300, bg = "#ffffff")

ggplot(res[res$criterion == "WAIC", ], aes(x = mod, y = val, fill = mod)) +
  facet_wrap(~geno) +
  geom_col() +
  geom_errorbar(aes(ymin = val - 3 * se_val, ymax = val + 3 * se_val),
    width = 0.25, linewidth = 0.5) +
  geom_point(data = res[res$criterion == "WAIC" & res$sig_improvement, ]) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(title = "WAIC",
    subtitle = "Significant differences marked with a point",
    x = "Model", y = "WAIC"
  )
ggsave("supplement/waic.png", width = 8, height = 6, dpi = 300, bg = "#ffffff")
