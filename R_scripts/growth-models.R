library(tidyverse)
library(pcvr)
library(brms)
#library(cmdstanr)
library(patchwork)

setwd("")

dat = read.csv(file = "data-file.csv")

# prior testing

test_priors <- list("A" = 4, "B" = 2.5)
plotPrior(test_priors, "power law")[[1]]

test_priors <- list("A" = 20, "B" = 6, "C" = 2)
plotPrior(test_priors, "logistic")[[1]]

test_priors <- list("A" = 3)
plotPrior(test_priors, "linear")[[1]]

################
# pcvr growth modeling

ss <- growthSS(model = "power law", form = area_cm2 ~ day | replicate / treatment,
               df = Geno, start = list("A" = 4, "B" = 2.5), type = "brms", sigma = "spline")

##
Geno = dat %>% filter(genotype == "G1")
genotype <- unique(Geno$genotype_name)

ss <- growthSS(model = "power law", form = area_cm2 ~ day | replicate / treatment,
               df = Geno, start = list("A" = 4, "B" = 2.5), type = "brms", sigma = "spline")

fit <- fitGrowth(ss)

saveRDS(fit, paste("./growth-model/fit_objects/", genotype, "_fit.RDS", sep = ""))

pdf(file = paste("./growth-model/growth_plots/", genotype, "_growth_plot.pdf", sep = ""), width = 8, height = 3)
growthPlot(fit = fit, form = ss$pcvrForm, df = ss$df) +
  facet_wrap("treatment", nrow = 1)+
  scale_x_continuous(breaks = c(0,1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12))
dev.off()

##
Geno = dat %>% filter(genotype == "G2")
genotype <- unique(Geno$genotype_name)

ss <- growthSS(model = "power law", form = area_cm2 ~ day | replicate / treatment,
               df = Geno, start = list("A" = 4, "B" = 2.5), type = "brms", sigma = "spline")

fit <- fitGrowth(ss)

summary(fit)

saveRDS(fit, paste("./growth-model/fit_objects/", genotype, "_fit.RDS", sep = ""))

pdf(file = paste("./growth-model/growth_plots/", genotype, "_growth_plot.pdf", sep = ""), width = 8, height = 3)
growthPlot(fit = fit, form = ss$pcvrForm, df = ss$df) +
  facet_wrap("treatment", nrow = 1)+
  scale_x_continuous(breaks = c(0,1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12))
dev.off()

##
Geno = dat %>% filter(genotype == "G3")
genotype <- unique(Geno$genotype_name)

ss <- growthSS(model = "power law", form = area_cm2 ~ day | replicate / treatment,
               df = Geno, start = list("A" = 4, "B" = 2.5), type = "brms", sigma = "spline")

fit <- fitGrowth(ss)
summary(fit)

saveRDS(fit, paste("./growth-model/fit_objects/", genotype, "_fit.RDS", sep = ""))

pdf(file = paste("./growth-model/growth_plots/", genotype, "_growth_plot.pdf", sep = ""), width = 8, height = 3)
growthPlot(fit = fit, form = ss$pcvrForm, df = ss$df) +
  facet_wrap("treatment", nrow = 1)+
  scale_x_continuous(breaks = c(0,1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12))
dev.off()

##
Geno = dat %>% filter(genotype == "G4")
genotype <- unique(Geno$genotype_name)

ss <- growthSS(model = "power law", form = area_cm2 ~ day | replicate / treatment,
               df = Geno, start = list("A" = 4, "B" = 2.5), type = "brms", sigma = "spline")

fit <- fitGrowth(ss)
summary(fit)

saveRDS(fit, paste("./growth-model/fit_objects/", genotype, "_fit.RDS", sep = ""))

pdf(file = paste("./growth-model/growth_plots/", genotype, "_growth_plot.pdf", sep = ""), width = 8, height = 3)
growthPlot(fit = fit, form = ss$pcvrForm, df = ss$df) +
  facet_wrap("treatment", nrow = 1)+
  scale_x_continuous(breaks = c(0,1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12))
dev.off()

##
Geno = dat %>% filter(genotype == "G5")
genotype <- unique(Geno$genotype_name)

ss <- growthSS(model = "power law", form = area_cm2 ~ day | replicate / treatment,
               df = Geno, start = list("A" = 4, "B" = 2.5), type = "brms", sigma = "spline")

fit <- fitGrowth(ss)
summary(fit)

saveRDS(fit, paste("./growth-model/fit_objects/", genotype, "_fit.RDS", sep = ""))

pdf(file = paste("./growth-model/growth_plots/", genotype, "_growth_plot.pdf", sep = ""), width = 8, height = 3)
growthPlot(fit = fit, form = ss$pcvrForm, df = ss$df) +
  facet_wrap("treatment", nrow = 1)+
  scale_x_continuous(breaks = c(0,1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12))
dev.off()

##
Geno = dat %>% filter(genotype == "G6")
genotype <- unique(Geno$genotype_name)

ss <- growthSS(model = "power law", form = area_cm2 ~ day | replicate / treatment,
               df = Geno, start = list("A" = 4, "B" = 2.5), type = "brms", sigma = "spline")

fit <- fitGrowth(ss)

saveRDS(fit, paste("./growth-model/fit_objects/", genotype, "_fit.RDS", sep = ""))

pdf(file = paste("./growth-model/growth_plots/", genotype, "_growth_plot.pdf", sep = ""), width = 8, height = 3)
growthPlot(fit = fit, form = ss$pcvrForm, df = ss$df) +
  facet_wrap("treatment", nrow = 1)+
  scale_x_continuous(breaks = c(0,1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12))
dev.off()

##
Geno = dat %>% filter(genotype == "G42")
genotype <- unique(Geno$genotype_name)

ss <- growthSS(model = "power law", form = area_cm ~ day | replicate / treatment,
               df = Geno, start = list("A" = 4, "B" = 2.5), type = "brms", sigma = "spline")

fit <- fitGrowth(ss)
summary(fit)

saveRDS(fit, paste("./growth-model/fit_objects/", genotype, "_fit.RDS", sep = ""))

fit <- readRDS("./growth-model/fit_objects/b73_fit.RDS")

#pdf(file = paste("./growth-model/growth_plots/", genotype, "_growth_plot.pdf", sep = ""), width = 8, height = 3)
growthPlot(fit = fit, form = ss$pcvrForm, df = ss$df) +
  facet_wrap("treatment", nrow = 1)+
  scale_x_continuous(breaks = c(0,1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13))
dev.off()

Geno <- subdat %>% filter(genotype == genolist[1])


##
# loop through the genotypes left

subdat <- dat %>% filter(genotype == "G1"| genotype == "G2"| genotype == "G3"| genotype == "G4"| genotype == "G5"|
                            genotype == "G6"| genotype == "G7"| genotype == "G8"| genotype == "G9"| genotype == "G10"|
                            genotype == "G11"| genotype == "G12"| genotype == "G13"| genotype == "G14")
genolist <- unique(subdat$genotype)

for (i in genolist) {
  Geno <- subdat %>% filter(genotype == i)
  genotype <- unique(Geno$genotype_name)
  ss <- growthSS(model = "power law", form = area_cm ~ day | replicate / treatment,
                 df = Geno, start = list("A" = 4, "B" = 2.5), type = "brms", sigma = "spline")
  fit <- fitGrowth(ss)
  saveRDS(fit, paste("./growth-model/fit_objects/", genotype, "_fit.RDS", sep = ""))
}

##
# loading fit objects and plotting growth curves

setwd("/growth-model/fit_objects/")
files <- list.files("./")
dfs <- lapply(files, readRDS)

n <- as.data.frame(files)
n <- separate(n, files, into = c("names", "ext"), sep = "_", remove = F)
nms <- n$names
names(dfs) <- nms
form = area_cm ~ day | replicate/treatment

##

b108 <- growthPlot(fit = dfs$b108, form = form, df = dat[dat$genotype_name=="b108",]) + 
  facet_wrap("treatment", nrow = 1) +
  labs(title = "B108")

b73 <- growthPlot(fit = dfs$b73, form = form, df = dat[dat$genotype_name=="b73",]) + 
  facet_wrap("treatment", nrow = 1) +
  labs(title = "B73")

b84 <- growthPlot(fit = dfs$b84, form = form, df = dat[dat$genotype_name=="b84",]) + 
  facet_wrap("treatment", nrow = 1) +
  labs(title = "B84")

b97 <- growthPlot(fit = dfs$b97, form = form, df = dat[dat$genotype_name=="b97",]) + 
  facet_wrap("treatment", nrow = 1) +
  labs(title = "B97")

pdf(file = "./growth-model/growth_plots/plots1.pdf", width = 8, height = 11)
b108 / b73 / b84 / b97
dev.off()

##
cml103 <- growthPlot(fit = dfs$cml103, form = form, df = dat[dat$genotype_name=="cml103",]) + 
  facet_wrap("treatment", nrow = 1) +
  labs(title = "cml103")

cml228 <- growthPlot(fit = dfs$cml228, form = form, df = dat[dat$genotype_name=="cml228",]) + 
  facet_wrap("treatment", nrow = 1) +
  labs(title = "cml228")

cml247 <- growthPlot(fit = dfs$cml247, form = form, df = dat[dat$genotype_name=="cml247",]) + 
  facet_wrap("treatment", nrow = 1) +
  labs(title = "cml247")

cml277 <- growthPlot(fit = dfs$cml277, form = form, df = dat[dat$genotype_name=="cml277",]) + 
  facet_wrap("treatment", nrow = 1) +
  labs(title = "cml277")

pdf(file = "./growth-model/growth_plots/plots2.pdf", width = 8, height = 11)
cml103 / cml228 / cml247 / cml277
dev.off()

##
cml322 <- growthPlot(fit = dfs$cml322, form = form, df = dat[dat$genotype_name=="cml322",]) + 
  facet_wrap("treatment", nrow = 1) +
  labs(title = "cml322")

cml333 <- growthPlot(fit = dfs$cml333, form = form, df = dat[dat$genotype_name=="cml333",]) + 
  facet_wrap("treatment", nrow = 1) +
  labs(title = "cml333")

cml52 <- growthPlot(fit = dfs$cml52, form = form, df = dat[dat$genotype_name=="cml52",]) + 
  facet_wrap("treatment", nrow = 1) +
  labs(title = "cml52")

cml69 <- growthPlot(fit = dfs$cml69, form = form, df = dat[dat$genotype_name=="cml69",]) + 
  facet_wrap("treatment", nrow = 1) +
  labs(title = "cml69")

pdf(file = "./growth-model/growth_plots/plots3.pdf", width = 8, height = 11)
cml322 / cml333 / cml52 / cml69
dev.off()

##
cr1ht <- growthPlot(fit = dfs$cr1ht, form = form, df = dat[dat$genotype_name=="cr1ht",]) + 
  facet_wrap("treatment", nrow = 1) +
  labs(title = "cr1ht")

hp301 <- growthPlot(fit = dfs$hp301, form = form, df = dat[dat$genotype_name=="hp301",]) + 
  facet_wrap("treatment", nrow = 1) +
  labs(title = "hp301")

il14h <- growthPlot(fit = dfs$il14h, form = form, df = dat[dat$genotype_name=="il14h",]) + 
  facet_wrap("treatment", nrow = 1) +
  labs(title = "il14h")

ki11 <- growthPlot(fit = dfs$ki11, form = form, df = dat[dat$genotype_name=="ki11",]) + 
  facet_wrap("treatment", nrow = 1) +
  labs(title = "ki11")

pdf(file = "./growth-model/growth_plots/plots4.pdf", width = 8, height = 11)
cr1ht / hp301 / il14h / ki11
dev.off()

##
ki3 <- growthPlot(fit = dfs$ki3, form = form, df = dat[dat$genotype_name=="ki3",]) + 
  facet_wrap("treatment", nrow = 1) +
  labs(title = "ki3")

ky21 <- growthPlot(fit = dfs$ky21, form = form, df = dat[dat$genotype_name=="ky21",]) + 
  facet_wrap("treatment", nrow = 1) +
  labs(title = "ky21")

lh145 <- growthPlot(fit = dfs$lh145, form = form, df = dat[dat$genotype_name=="lh145",]) + 
  facet_wrap("treatment", nrow = 1) +
  labs(title = "lh145")

lh149 <- growthPlot(fit = dfs$lh149, form = form, df = dat[dat$genotype_name=="lh149",]) + 
  facet_wrap("treatment", nrow = 1) +
  labs(title = "lh149")

pdf(file = "./growth-model/growth_plots/plots5.pdf", width = 8, height = 11)
ki3 / ky21 / lh145 / lh149
dev.off()

##
lh197 <- growthPlot(fit = dfs$lh197, form = form, df = dat[dat$genotype_name=="lh197",]) + 
  facet_wrap("treatment", nrow = 1) +
  labs(title = "lh197")

m162w <- growthPlot(fit = dfs$m162w, form = form, df = dat[dat$genotype_name=="m162w",]) + 
  facet_wrap("treatment", nrow = 1) +
  labs(title = "m162w")

m37w <- growthPlot(fit = dfs$m37w, form = form, df = dat[dat$genotype_name=="m37w",]) + 
  facet_wrap("treatment", nrow = 1) +
  labs(title = "m37w")

minimaize <- growthPlot(fit = dfs$minimaize, form = form, df = dat[dat$genotype_name=="minimaize",]) + 
  facet_wrap("treatment", nrow = 1) +
  labs(title = "minimaize")

pdf(file = "./growth-model/growth_plots/plots6.pdf", width = 8, height = 11)
lh197 / m162w / m37w / minimaize
dev.off()

##
mo17 <- growthPlot(fit = dfs$mo17, form = form, df = dat[dat$genotype_name=="mo17",]) + 
  facet_wrap("treatment", nrow = 1) +
  labs(title = "mo17")

mo18w <- growthPlot(fit = dfs$mo18w, form = form, df = dat[dat$genotype_name=="mo18w",]) + 
  facet_wrap("treatment", nrow = 1) +
  labs(title = "mo18w")

ms71 <- growthPlot(fit = dfs$ms71, form = form, df = dat[dat$genotype_name=="ms71",]) + 
  facet_wrap("treatment", nrow = 1) +
  labs(title = "ms71")

nc350 <- growthPlot(fit = dfs$nc350, form = form, df = dat[dat$genotype_name=="nc350",]) + 
  facet_wrap("treatment", nrow = 1) +
  labs(title = "nc350")

pdf(file = "./growth-model/growth_plots/plots7.pdf", width = 8, height = 11)
mo17 / mo18w / ms71 / nc350
dev.off()

##
nc358 <- growthPlot(fit = dfs$nc358, form = form, df = dat[dat$genotype_name=="nc358",]) + 
  facet_wrap("treatment", nrow = 1) +
  labs(title = "nc358")

nkh8431 <- growthPlot(fit = dfs$nkh8431, form = form, df = dat[dat$genotype_name=="nkh8431",]) + 
  facet_wrap("treatment", nrow = 1) +
  labs(title = "nkh8431")

oh33 <- growthPlot(fit = dfs$oh33, form = form, df = dat[dat$genotype_name=="oh33",]) + 
  facet_wrap("treatment", nrow = 1) +
  labs(title = "oh33")

oh43 <- growthPlot(fit = dfs$oh43, form = form, df = dat[dat$genotype_name=="oh43",]) + 
  facet_wrap("treatment", nrow = 1) +
  labs(title = "oh43")

pdf(file = "./growth-model/growth_plots/plots8.pdf", width = 8, height = 11)
nc358 / nkh8431 / oh33 / oh43
dev.off()

##
oh7b <- growthPlot(fit = dfs$oh7b, form = form, df = dat[dat$genotype_name=="oh7b",]) + 
  facet_wrap("treatment", nrow = 1) +
  labs(title = "oh7b")

p39 <- growthPlot(fit = dfs$p39, form = form, df = dat[dat$genotype_name=="p39",]) + 
  facet_wrap("treatment", nrow = 1) +
  labs(title = "p39")

ph207 <- growthPlot(fit = dfs$ph207, form = form, df = dat[dat$genotype_name=="ph207",]) + 
  facet_wrap("treatment", nrow = 1) +
  labs(title = "ph207")

phb47 <- growthPlot(fit = dfs$phb47, form = form, df = dat[dat$genotype_name=="phb47",]) + 
  facet_wrap("treatment", nrow = 1) +
  labs(title = "phb47")

pdf(file = "./growth-model/growth_plots/plots9.pdf", width = 8, height = 11)
oh7b / p39 / ph207 / phb47
dev.off()

##
phg29 <- growthPlot(fit = dfs$phg29, form = form, df = dat[dat$genotype_name=="phg29",]) + 
  facet_wrap("treatment", nrow = 1) +
  labs(title = "phg29")

phg47 <- growthPlot(fit = dfs$phg47, form = form, df = dat[dat$genotype_name=="phg47",]) + 
  facet_wrap("treatment", nrow = 1) +
  labs(title = "phg47")

phg50 <- growthPlot(fit = dfs$phg50, form = form, df = dat[dat$genotype_name=="phg50",]) + 
  facet_wrap("treatment", nrow = 1) +
  labs(title = "phg50")

phg72 <- growthPlot(fit = dfs$phg72, form = form, df = dat[dat$genotype_name=="phg72",]) + 
  facet_wrap("treatment", nrow = 1) +
  labs(title = "phg72")

pdf(file = "./growth-model/growth_plots/plots10.pdf", width = 8, height = 11)
phg29 / phg47 / phg50 / phg72
dev.off()

##
phg80 <- growthPlot(fit = dfs$phg80, form = form, df = dat[dat$genotype_name=="phg80",]) + 
  facet_wrap("treatment", nrow = 1) +
  labs(title = "phg80")

phj40 <- growthPlot(fit = dfs$phj40, form = form, df = dat[dat$genotype_name=="phj40",]) + 
  facet_wrap("treatment", nrow = 1) +
  labs(title = "phj40")

phw65 <- growthPlot(fit = dfs$phw65, form = form, df = dat[dat$genotype_name=="phw65",]) + 
  facet_wrap("treatment", nrow = 1) +
  labs(title = "phw65")

tx303 <- growthPlot(fit = dfs$tx303, form = form, df = dat[dat$genotype_name=="tx303",]) + 
  facet_wrap("treatment", nrow = 1) +
  labs(title = "tx303")

pdf(file = "./growth-model/growth_plots/plots11.pdf", width = 8, height = 11)
phg80 / phj40 / phw65 / tx303
dev.off()

##
tzi8 <- growthPlot(fit = dfs$tzi8, form = form, df = dat[dat$genotype_name=="tzi8",]) + 
  facet_wrap("treatment", nrow = 1) +
  labs(title = "tzi8")

w604s <- growthPlot(fit = dfs$w604s, form = form, df = dat[dat$genotype_name=="w604s",]) + 
  facet_wrap("treatment", nrow = 1) +
  labs(title = "w604s")

w605s <- growthPlot(fit = dfs$w605s, form = form, df = dat[dat$genotype_name=="w605s",]) + 
  facet_wrap("treatment", nrow = 1) +
  labs(title = "w605s")

pdf(file = "./growth-model/growth_plots/plots12.pdf", width = 8, height = 8.25)
tzi8 / w604s / w605s
dev.off()

#####
setwd("growth-model/fit_objects/")
files <- list.files("./")
dfs <- lapply(files, readRDS)

n <- as.data.frame(files)
n <- separate(n, files, into = c("names", "ext"), sep = "_", remove = F)
nms <- n$names
names(dfs) <- nms

genolist <- unique(dat$genotype_name)
genolist <- genolist[2:47]

genotype <- "b73"
s <- summary(dfs[[genotype]])
y <- s$fixed
y <- rownames_to_column(y, var = "Parameter")
y$genotype <- genotype
df <- y

for (i in genolist) {
  genotype <- i
  x <- summary(dfs[[genotype]])
  z <- x$fixed
  z <- rownames_to_column(z, var = "Parameter")
  z$genotype <- genotype
  df <- rbind(df,z)
}

write.csv(df, file = "growth-model/parameter-estimates.csv", row.names = F)

