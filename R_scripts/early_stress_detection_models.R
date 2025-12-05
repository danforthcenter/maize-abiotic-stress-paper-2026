setwd("~/scripts/gehan_lab/annaCasto/adaptive_design_maize_phenotyper/")

library(ggplot2)
library(data.table)
library(pcvr)
library(patchwork)
library(brms)
testing_model_building <- FALSE
theme_set(pcv_theme() + theme(axis.text.x.bottom = element_text(hjust = 0.5)))
df <- as.data.frame(fread("2024-11-13_VIS_SV_MG001_MG002_dead-plant-filtered_calcualted-values_D-DH-H-C.csv"))
dim(df)
head(df)

#* ***** `Checking data` *****
#* just taking a look at the data to see if there are any weird looking points

ggplot(df[df$genotype %in% paste0("G", 1:12), ],
       aes(x = day, y = area_cm2, group = barcode, color = treatment)) +
  facet_wrap(~genotype) +
  geom_line() +
  scale_x_continuous(breaks = seq(3, 12, 3))

ggplot(df[df$genotype %in% paste0("G", 13:24), ],
       aes(x = day, y = area_cm2, group = barcode, color = treatment)) +
  facet_wrap(~genotype) +
  geom_line() +
  scale_x_continuous(breaks = seq(3, 12, 3))

ggplot(df[df$genotype %in% paste0("G", 25:36), ],
       aes(x = day, y = area_cm2, group = barcode, color = treatment)) +
  facet_wrap(~genotype) +
  geom_line() +
  scale_x_continuous(breaks = seq(3, 12, 3))

ggplot(df[df$genotype %in% paste0("G", 37:47), ],
       aes(x = day, y = area_cm2, group = barcode, color = treatment)) +
  facet_wrap(~genotype) +
  geom_line() +
  scale_x_continuous(breaks = seq(3, 12, 3))
#* everything looks totally reasonable to me for area at least.

#* ***** `Model Building` *****
#* I think we've been using power law models for these.
#* One of the awkward bits about power law is that the parameters
#* are correlated since both are acting directly on time:
#* `y ~ A * x ^ B`
#* One way we might get around that is by only estimating B per group.
#* We could technically do that with A instead, but B has more interesting
#* interpretation I think. If we were to hold B constant then I'd think we should use
#* an exponential model:
#* `y ~ A * e ^ (B * x)`
#* instead.
#* 
#* It will almost certainly be a worse fit if we only estimate one parameter instead of 2
#* by group but it might be worth it for interpretation.
#* I'll do a bit and check the posterior predictive intervals.
#* 
#* Before I can do that I'm going to push the days by 1 so that I don't have a bunch of
#* stuff getting multiplied by 0. We could add an intercept term to handle it but I'd rather
#* treat time as a little more consistent.

df$day_1 <- df$day + 1

if (testing_model_building) {
  ggplot() +
    geom_line(data = growthSim("power law", params = list("A" = 3, "B" = 1), t = max(df$day), n = 10),
              aes(x = time, y = y, group = id), color = "blue") +
    geom_line(data = growthSim("power law", params = list("A" = 6, "B" = 0.9), t = max(df$day), n = 10),
              aes(x = time, y = y, group = id), color = "red") +
    geom_line(data = df[df$genotype == "G1", ],
              aes(x = day_1, y = area_cm2, group = barcode), color = "gray40") +
    facet_wrap(~treatment) +
    labs(title = "Considering global priors")
  
  
  ss1 <- growthSS("power law", area_cm2 ~ day_1 | barcode / treatment,
                  sigma = "logistic", df = df[df$genotype == "G1", ],
                  pars = "B", start = list("A" = 3, "B" = 1,
                                           "sigmaA" = 5, "sigmaB" = 5, "sigmaC" = 3))
  ss1 #* ah, forgot that pars doesn't work for brms models. I'll change that I think.
  #* for now that means I'm going to use brms.
  
  
  bf1 <- bf(
    area_cm2 ~ A * day_1 ^ B,
    A ~ 1,
    B ~ 0 + treatment,
    nlf(sigma ~ sigmaA/(1 + exp((sigmaB - day_1)/sigmaC))),
    sigmaA + sigmaB + sigmaC ~ 0 + treatment,
    nu ~ 1,
    autocor = ~arma(~day_1 | barcode:treatment, 1, 1),
    nl = TRUE, family = "student"
  )
  prior1 <- ss1$prior
  m1 <- brm(bf1, data = df[df$genotype == "G1", ], prior = prior1,
            cores = 4, chains = 4, iter = 2000,
            backend = "cmdstanr", control = list(adapt_delta = 0.99, max_treedepth = 20))
  m1
  growthPlot(m1, form = ss1$pcvrForm, df = ss1$df)
  
  bf2 <- bf(
    area_cm2 ~ A * day_1 ^ B,
    A ~ 0 + treatment,
    B ~ 0 + treatment,
    nlf(sigma ~ sigmaA/(1 + exp((sigmaB - day_1)/sigmaC))),
    sigmaA + sigmaB + sigmaC ~ 0 + treatment,
    nu ~ 1,
    autocor = ~arma(~day_1 | barcode:treatment, 1, 1),
    nl = TRUE, family = "student"
  )
  
  m2 <- brm(bf2, data = df[df$genotype == "G1", ], prior = prior1,
            cores = 4, chains = 4, iter = 2000,
            backend = "cmdstanr", control = list(adapt_delta = 0.99, max_treedepth = 20))
  m2
  growthPlot(m2, form = ss1$pcvrForm, df = ss1$df)
  
  bf3 <- bf(
    area_cm2 ~ A * day_1 ^ B,
    A ~ 0 + treatment,
    B ~ 0 + treatment,
    sigma ~ s(day_1, by = treatment),
    nu ~ 1,
    autocor = ~arma(~day_1 | barcode:treatment, 1, 1),
    nl = TRUE, family = "student"
  )
  default_prior <- get_prior(bf3, data = df)
  prior3 <- rbind(prior1[1:3, ], default_prior[default_prior$dpar == "sigma", ])
  m3 <- brm(bf3, data = df[df$genotype == "G1", ], prior = prior3,
            cores = 4, chains = 4, iter = 2000,
            backend = "cmdstanr", control = list(adapt_delta = 0.99, max_treedepth = 20))
  m3
  
  bf4 <- bf(
    area_cm2 ~ A * day_1 ^ B,
    A ~ 1,
    B ~ 0 + treatment,
    sigma ~ s(day_1, by = treatment),
    nu ~ 1,
    autocor = ~arma(~day_1 | barcode:treatment, 1, 1),
    nl = TRUE, family = "student"
  )
  m4 <- brm(bf4, data = df[df$genotype == "G1", ], prior = prior3,
            cores = 4, chains = 4, iter = 2000,
            backend = "cmdstanr", control = list(adapt_delta = 0.99, max_treedepth = 20))
  m4
  
  bf5 <- bf( # note I don't like the idea of doing this, but we'll try it.
    area_cm2 ~ A * day_1 ^ B,
    A ~ 0 + treatment,
    B ~ 1,
    sigma ~ s(day_1, by = treatment),
    nu ~ 1,
    autocor = ~arma(~day_1 | barcode:treatment, 1, 1),
    nl = TRUE, family = "student"
  )
  m5 <- brm(bf4, data = df[df$genotype == "G1", ], prior = prior3,
            cores = 4, chains = 4, iter = 2000,
            backend = "cmdstanr", control = list(adapt_delta = 0.99, max_treedepth = 20))
  m5
  
  growthPlot(m1, form = ss1$pcvrForm, df = ss1$df) +
    labs(title = "m1")
  growthPlot(m2, form = ss1$pcvrForm, df = ss1$df) +
    labs(title = "m2")
  growthPlot(m3, form = ss1$pcvrForm, df = ss1$df) +
    labs(title = "m3")
  growthPlot(m4, form = ss1$pcvrForm, df = ss1$df) +
    labs(title = "m4")
  growthPlot(m5, form = ss1$pcvrForm, df = ss1$df) +
    labs(title = "m5")
  
  m1 <- add_criterion(m1, "loo")
  m2 <- add_criterion(m2, "loo")
  m3 <- add_criterion(m3, "loo")
  m4 <- add_criterion(m4, "loo")
  m5 <- add_criterion(m5, "loo")
  
  loo_compare(m1, m2, m3, m4, m5) # m2 preferred, makes sense. m1 is essentially the same.
}

#* ***** `Model to Use` *****
#* 
#* I think that m3 looks the best in terms of posterior predictive.
#* m2 is preferred by loo
#* m1 is almost identical in loo to m2 and has simpler hypothesis test interpretation if A is constant.
#*
#* Given how those look I'm settling on m1 to start with. If it fails dramatically with some group
#* then I might adjust, but it would be nice to keep the hypothesis testing simpler.
#* The variance submodel can also change, the reward for loo is that it's way simpler to fit
#* even a 3 parameter model compared to something like a spline. I might play with the prior for it
#* a little and try to push it to be a bit tighter, but the data is definitely getting off of that
#* prior anyway.

bf1 <- bf(
  area_cm2 ~ A * day_1 ^ B,
  A ~ 1,
  B ~ 0 + treatment,
  nlf(sigma ~ sigmaA/(1 + exp((sigmaB - day_1)/sigmaC))),
  sigmaA + sigmaB + sigmaC ~ 0 + treatment,
  nu ~ 1,
  autocor = ~arma(~day_1 | barcode:treatment, 1, 1),
  nl = TRUE, family = "student"
)
prior1 <- set_prior("normal(2.7, 0.8)", dpar = "nu", class = "Intercept") +
  set_prior("lognormal(log(3), 0.25)", nlpar = "A", lb = 0) + # will be estimated per geno
  set_prior("normal(1, 0.25)", nlpar = "B", lb = 0) + # assuming linear is average
  set_prior("lognormal(log(3), 0.25)", nlpar = "sigmaA", lb = 0) + # encouraging low asymptote
  set_prior("lognormal(log(5), 0.25)", nlpar = "sigmaB", lb = 0) + # encouraging reaching asymptote fast
  set_prior("lognormal(log(5), 0.25)", nlpar = "sigmaC", lb = 0) # encouraging reaching asymptote fast

#* ***** `Fit model to all Genotypes` *****

if (!file.exists("all_mods.rdata")) {
  options(cmdstanr_write_stan_file_dir = getwd())
  
  mod_skeleton <- brm(bf1, data = df[df$genotype == "G1" & df$day <= 4, ],
                      prior = prior1,
                      cores = 4, chains = 4, iter = 2000,
                      backend = "cmdstanr", control = list(adapt_delta = 0.99,
                                                           max_treedepth = 20),
                      file = "mod_skeleton")
  
  all_mods <- lapply(unique(df$genotype), function(GENO) {
    cat(paste0("\n\n\nGeno: ", GENO, "\n\n\n"))
    models <- lapply(seq(4, 12, 2), function(tm) {
      cat(paste0("\n\nGeno: ", GENO, " - ", tm , "\n\n"))
      m_iter <- update(mod_skeleton,
                       newdata = df[df$genotype == GENO & df$day <= tm, ],
                       cores = 4,
                       backend = "cmdstanr",
                       control = list(adapt_delta = 0.99,
                                      max_treedepth = 20))
      return(m_iter)
    })
    out <- setNames(list(models), GENO)
    return(out)
  })
  names(all_mods) <- unique(df$genotype)
  save(all_mods, file = "all_mods.rdata")
} else {
  print(load("all_mods.rdata"))
}

#* ***** `Visualize all Models' Posterior Predictive Distributions` *****

patch_des <- c(
  area(1, 2, 1, 3),
  area(1, 4, 1, 5),
  area(2, 1, 2, 2),
  area(2, 3, 2, 4),
  area(2, 5, 2, 6)
)
# plot(patch_des)

patches <- lapply(names(all_mods), function(nm) {
  # grab just the models of this genotype
  l <- all_mods[[nm]][[1]]
  # make pcvr::growthPlot of each (this should really be a geom)
  plots <- lapply(l, function(m) {
    growthPlot(m,
               form = area_cm2 ~ day_1 | barcode/treatment,
               df = df[df$genotype == nm, ])
  }
  )
  # get limits of final plot
  x_lims <- ggplot2::layer_scales(plots[[length(plots)]])$x$range$range
  y_lims <- ggplot2::layer_scales(plots[[length(plots)]])$y$range$range
  # apply limits to each plot
  plots <- lapply(plots, function(p) {
    p_mod <- p +
      coord_cartesian(xlim = x_lims, ylim = y_lims) +
      labs(title = paste0("modeled to day ", max(p$data$day_1)),
           x = "Day")
    return(p_mod)
  })
  # assemble a patchwork of the plots
  patch <- Reduce("+", plots) +
    plot_annotation(title = nm) +
    plot_layout(
      axes = "collect",
      axis_titles = "collect",
      guides = "collect",
      design = patch_des) &
    theme(legend.position = "right")
  ggsave(
    paste0("model_plots/", nm, "_patch.png"),
    patch,
    width = 11, height = 7, dpi = 300,
    bg = "#ffffff")
  return(patch)
})

#* ***** `Visualize all Models posterior Distributions` *****
#* ***** `Univariate Density`
#* This is what pcvr::distPlot is for, but it's a relatively untested function so we'll see if it's
#* flexible enough for what I want here.
devtools::load_all("~/pcvr")

for (i in seq_along(unique(df$genotype))) {
  geno <- unique(df$genotype)[i]
  p <- distributionPlot(all_mods[[i]][[geno]],
                   form = area_cm2 ~ day_1 | barcode / treatment,
                   df = df[df$genotype == geno, ],
                   params = "B",
                   virOptions = c(
                     "plasma", "mako", "viridis", "cividis"
                   ))
  p[[1]] <- p[[1]] + ggplot2::labs(x = "Days", y = ~~Area~cm^2)
  p[[3]] <- p[[3]] + ggplot2::labs(x = "Days", y = ~~Area~cm^2)
  p[[5]] <- p[[5]] + ggplot2::labs(x = "Days", y = ~~Area~cm^2)
  p[[7]] <- p[[7]] + ggplot2::labs(x = "Days", y = ~~Area~cm^2)
  p <- p +
    plot_layout(axes = "collect", axis_titles = "collect") +
    plot_annotation(title = paste0(geno, " posteriors of B over time")
                    #* PENDING:
                    #* I'll want to add some information about the hypothesis testing here.
                    #* I'm thinking "first time X hypothesis is true: Y day"
                    #* or something like that.
                    )
  ggsave(
    paste0("distribution_plots/", geno, "_distPlot.png"),
    p,
    width = 8, height = 6, dpi = 300, bg = "#ffffff"
  )
}

#* ***** `Bivariate Density`
#* This will be a plot of A and B like I first did for Dom's power law stuff to show the correlation
#* differences between stressors.
#* 
#* Not sure that I need this given that A is not estimated by group? It might still be interesting.

for (mod_list in all_mods) {
  geno <- names(mod_list)
  ddfl <- do.call(rbind, lapply(mod_list[[geno]], function(fit) {
    ddf <- as.data.frame(fit)
    ddf <- ddf[, grepl("^b_(A|B)", colnames(ddf))]
    ndraws <- nrow(ddf)
    ddf$id <- seq_len(ndraws)
    ddf <- reshape2::melt(ddf, id.vars = "id")
    ddf$time <- max(fit$data[["day_1"]])
    ddf$variable <- gsub("^b_", "", ddf$variable)
    ddf$variable <- gsub("Drought_Heat", "DroughtHeat", ddf$variable)
    unique(ddf$variable)
    ddf$par <- gsub("_(treatment)*.*", "", ddf$variable)
    ddf$trt <- gsub(".*_(treatment)*", "", ddf$variable)
    par_dfs <- split(ddf, ddf$par)
    par_dfs$A <- do.call(rbind, lapply(unique(par_dfs$B$trt), function(trt) {
      d <- par_dfs$A
      d$trt <- trt
      d$variable <- paste0("A_treatment", trt)
      return(d)
    }))
    ddf <- do.call(rbind, par_dfs)
    ddfl <- reshape2::dcast(ddf[,-2], ...  ~ par, value.var = "value")
  }))
  p <- ggplot(ddfl,
         aes(x = A, y = B, color = interaction(trt, time))) +
    geom_point(data = ddfl[ddfl$id %in% sample(seq_len(ndraws), 500),],
               size = 0.25, alpha = 0.25) +
    geom_density2d(aes(group = interaction(trt, time)), linewidth = 0.5) +
    scale_color_manual(values = rbind(
      viridis::plasma(5, 1, 0.3, 0.9),
      viridis::mako(5, 1, 0.3, 0.9),
      viridis::cividis(5, 1, 0.3, 0.9),
      viridis::viridis(5, 1, 0.3, 0.9)
    )) +
    guides(color = guide_legend(override.aes = list(linewidth = 5), title = "Treatment +\nTime")) +
    theme_bw() +
    theme(legend.position = "bottom") +
    labs(x = "A Draw (estimated at Genotype level)", y = "B Draw",
         title = paste0("Joint Posterior Distribution of ", geno))
  ggsave(
    paste0("bivariate_posterior_plots/", geno, "_biPostPlot.png"),
    p,
    width = 9, height = 7, dpi = 300, bg = "#ffffff"
  )
}

#* ***** `Hypothesis Testing`
#* 
#* Here I'm not totally sure what I want to test.
#* Intuitively I think that something like "X is less than control"
#* makes sense to me since higher values of B mean faster and probably healthier plants.
#* So in that case I'd want to:
#* For every genotype
#*   On each day:
#*     Test B_trt_x * EFFECT_SIZE < B_control
#*     report posterior probability and estimated effect size
#*   Find first time (X) that post.prob > CUTOFF
#*   report hypothesis test from Time_X
#* Check number of Genotypes for which subsequent days added "useful information"
#*   Barplot or something?
#* Check number of genotypes where conclusion could be made early
#*   Barplot or something?
#* Check distribution of conclusion times (semi-continuous)   
#*   histogram/freqpoly plot?
#* Check distribution of effect sizes over time
#*   histogram/freqpoly plot?

i <- 1
geno <- "G1"

h <- do.call(rbind, lapply(seq_along(all_mods), function(i) {
  geno <- names(all_mods[[i]])
  fits <- all_mods[[i]][[geno]]
  hyps_df <- do.call(rbind, lapply(fits, function(fit) {
    conditions <- c("Drought_Heat", "Drought", "Heat")
    do.call(rbind, lapply(conditions, function(trt) {
      hyp_df <- brms::hypothesis(fit,
                                 paste0("B_treatment", trt, " < B_treatmentControl")
                                 )$hypothesis
      hyp_df$trt <- trt
      hyp_df$time <- max(fit$data$day_1)
      hyp_df
    }))
  }))
  hyps_df$geno <- geno
  return(hyps_df)
}))

head(h)
dim(h)

d <- split(h, interaction(h$trt, h$geno))[[1]]
h2 <- do.call(rbind, lapply(split(h, interaction(h$trt, h$geno)), function(d) {
  d$sig <- FALSE
  d$firstSig <- NA
  d$firstEst <- NA
  d$estDiff <- NA
  if (any(d$Post.Prob > 0.95)) {
    d$sig[which(d$Post.Prob > 0.95)] <- TRUE
    d$firstSig <- d$time[min(which(d$Post.Prob > 0.95))]
    d$firstEst <- d$Estimate[min(which(d$Post.Prob > 0.95))]
    d$estDiff <- abs(d$firstEst - d$Estimate[nrow(d)])
  }
  d
}))
h3 <- h2[which(h2$time == h2$firstSig), ]

#* Check number of Genotypes for which subsequent days added "useful information" in estimate
#*   Barplot or something?
#*   Check difference in estimate from first significance to end?
#*   Correlation of time of first difference vs magnitude of final difference is Malia's interest

ggplot(h3,
       aes(x = as.factor(firstSig), y = estDiff)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() +
  facet_wrap(~trt) +
  labs(x = "Time of First Significant Hypothesis Test",
       y = "Difference in Mean Estimate between\nfirst significant and final hypothesis test")
ggsave("early_significant_tests_effect_size_boxplot.png", width = 7, height = 6, dpi = 300, 
       bg = "#ffffff")
head(h3)

#* Check number of genotypes where conclusion could be made early
#*   Barplot or something?

ggplot(h3, aes(x = firstSig)) +
  facet_wrap(~trt) +
  geom_bar() +
  geom_text(aes(label = after_stat(count)), vjust = -0.5, stat = "count") +
  scale_x_continuous(breaks = seq(5, 13, 2)) +
  geom_vline(xintercept = 12, color = "red", linetype = 5) +
  labs(x = "First time with significant hypothesis test relative to control",
       y = "Count")
ggsave("early_significant_tests_histogram.png", width = 7, height = 6, dpi = 300, 
       bg = "#ffffff")

write.csv(h3, "h3_tests.csv", row.names = FALSE)


#* I could also not plot those but do some error propagation and include bars around the points
#* and/or make a model of the correlation.

#* process is:
#* for each row of my hypotheses:
#*   Grab the relevant model (geno and time)
#*   Grab the relevant treatment's draws and controls draws from model
#*   Calculate mean of each set of draws
#*   Calculate error for each draw from mean
#*   Calculate mean difference
#*   Calculate propagated error of mean difference
#*   return mean difference with propagated error and relevant metadata
h4 <- do.call(rbind, lapply(seq_len(nrow(h3)), function(i) {
  geno <- h3[i, "geno"]
  trt <- h3[i, "trt"]
  firstSig <- h3[i, "firstSig"]
  if (firstSig == 13) {
    return(NULL) # because it's already the last day.
  }
  final_h <- h2[
    h2$geno == geno &
      h2$trt == trt &
      h2$time == 13,
  ]
  est_diff <- abs(final_h$Estimate - h3[i, "Estimate"])
  est_diff_error <- sqrt(final_h$Est.Error ^ 2 + h3[i, "Est.Error"] ^ 2)
  out <- data.frame(
    est_diff = est_diff,
    est_diff_error = est_diff_error,
    firstSig = firstSig,
    geno = geno,
    trt = trt
  )
  return(out)
}))


ggplot(h4, aes(x = firstSig, y = est_diff)) +
  facet_wrap(~trt) +
  geom_hline(yintercept = 0, color = "red", linetype = 5) +
  geom_pointrange(aes(ymin = est_diff - est_diff_error,
                      ymax = est_diff + est_diff_error),
                  size = 0.25, linewidth = 0.25,
                  position = position_jitter(height = 0, width = 0.35)) +
  scale_x_continuous(breaks = seq(5, 13, 2)) +
  labs(x = "Day of First Significant Hypothesis Test",
       y = "Difference in First vs Final Estimate\n(with propagated error)",
       title = "Early testing of Rate Parameter in Power-Law Growth")
ggsave("early_estimate_differences.png", width = 7, height = 6, dpi = 300, 
       bg = "#ffffff")

if (FALSE) {
  mf1 <- bf(est_diff | se(est_diff_error, sigma = TRUE) ~ firstSig + (1 | trt),
            sigma ~ s(firstSig, by = trt, k = 4), family = "student")
  mprior1 <- set_prior("normal(0, 1)", coef = "firstSig")
  mm1 <- brm(mf1, data = h4, prior = mprior1,
             iter = 2000, cores = 4, chains = 4,
             backend = "cmdstanr",
             control = list(adapt_delta = 0.99, max_treedepth = 20))
}

mf2 <- bf(est_diff | se(est_diff_error) ~ firstSig + (1 | trt), family = "student")
mprior2 <- set_prior("normal(0, 1)", coef = "firstSig")
mm2 <- brm(mf2, data = h4, prior = mprior2,
           iter = 2000, cores = 4, chains = 4,
           backend = "cmdstanr",
           control = list(adapt_delta = 0.99, max_treedepth = 20))

if (FALSE) {
  mf3 <- bf(est_diff ~ firstSig + (1 | trt), family = "student")
  mprior3 <- set_prior("normal(0, 1)", coef = "firstSig")
  mm3 <- brm(mf3, data = h4, prior = mprior3,
             iter = 2000, cores = 4, chains = 4,
             backend = "cmdstanr",
             control = list(adapt_delta = 0.99, max_treedepth = 20))
  
  mm1
  mm2
  mm3
  
  mm1 <- add_criterion(mm1, "loo")
  mm2 <- add_criterion(mm2, "loo")
  mm3 <- add_criterion(mm3, "loo")
  loo_compare(mm1, mm2, mm3)
  
  conditional_effects(mm2)
  
  mcmc_plot(mm1, type = "trace")
  mcmc_plot(mm2, type = "trace")
  mcmc_plot(mm3, type = "trace")
}


sub <- expand.grid(firstSig = seq(5, 11, 2),
                   trt = unique(h4$trt),
                   est_diff_error = mean(h4$est_diff_error))
test <- cbind(
  sub,
  predict(mm2, newdata = sub,
          probs = seq(0.01, 0.99, 0.02))[, -c(1:2)]
)

max_prime <- 0.99
min_prime <- 0.01
max_obs <- 49
min_obs <- 1
c1 <- (max_prime - min_prime) / (max_obs - min_obs)
longPreds <- do.call(rbind, lapply(seq_len(nrow(test)), function(r) {
  sub <- test[r, ]
  lp <- do.call(rbind, lapply(seq(1, 49, 2), function(i) {
    min <- paste0("Q", i)
    max <- paste0("Q", 100 - i)
    iter <- sub[, c("firstSig", "trt")]
    iter$q <- round(1 - (c1 * (i - max_obs) + max_prime), 2)
    iter$min <- sub[[min]]
    iter$max <- sub[[max]]
    return(iter)
  }))
  return(lp)
}))

mm2_hyp <- hypothesis(mm2, "firstSig < 0")$hyp
mm2_est <- mm2_hyp$Estimate
mm2_post_prob <- mm2_hyp$Post.Prob

ggplot(h4, aes(x = firstSig, y = est_diff)) +
  facet_wrap(~trt) +
  lapply(unique(longPreds$q), function(q) {
    ribbon_plot <- ggplot2::geom_ribbon(
      data = longPreds[longPreds$q == q, ],
      ggplot2::aes(
        x = .data[["firstSig"]],
        ymin = .data[["min"]],
        ymax = .data[["max"]],
        group = interaction(.data[["trt"]]),
        fill = .data[["q"]]
      ), alpha = 0.5
    )
    return(ribbon_plot)
  }) +
  scale_fill_viridis(option = "plasma", direction = -1) +
  geom_hline(yintercept = 0, color = "black", linetype = 5) +
  geom_pointrange(data = h4,
    aes(ymin = est_diff - est_diff_error,
                      ymax = est_diff + est_diff_error),
                  size = 0.25, linewidth = 0.25,
                  position = position_jitter(height = 0, width = 0.25)) +
  scale_x_continuous(breaks = seq(5, 13, 2)) +
  labs(x = "Day of First Significant Hypothesis Test",
       y = "Difference in First vs Final Estimate\n(with propagated error)",
       title = "Meta Analysis of Estimate Improvement given Additional Data",
       subtitle = paste0("Posterior probability of decreasing error with more data: ",
                         round(mm2_post_prob * 100, 2), "%\n",
                         "Estimated slope: ", round(mm2_est, 3)),
       fill = "CI")
ggsave("meta_model.png", width = 7, height = 6, dpi = 300, 
       bg = "#ffffff")















