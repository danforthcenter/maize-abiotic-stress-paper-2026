library(brms)
library(tidyverse)

##
# hypothesis testing
setwd("./")
dat = read.csv(file = "./pheno_data/11222024_C_D_H_DH_univar_traits_outlier_free.csv")

setwd("growth-model/fit_objects/")
files <- list.files("./")
dfs <- lapply(files, readRDS)
n <- as.data.frame(files)
n <- separate(n, files, into = c("names", "ext"), sep = "_", remove = F)
nms <- n$names
names(dfs) <- nms


##
#### Test if B < 1
temp_df <- data.frame(Hypothesis = NA,
                      Estimate = NA, 
                      Est.Error = NA, 
                      CI.Lower = NA,
                      CI.Upper = NA,
                      Evid.Ratio = NA, 
                      Post.Prob = NA, 
                      Star = NA)
genolist <- unique(dat$genotype_name)
genolist <- genolist[2:47]

genotype <- "b73"
ss <- growthSS(model = "power law", form = area_cm ~ day | replicate / treatment,
               df = dat[dat$genotype_name == genotype,], start = list("A" = 4, "B" = 2.5), type = "brms", sigma = "spline")
w <- testGrowth(ss, dfs[[genotype]], "B_treatmentControl > 1")
x <- testGrowth(ss, dfs[[genotype]], "B_treatmentDrought < 1")
y <- testGrowth(ss, dfs[[genotype]], "B_treatmentHeat < 1")
z <- testGrowth(ss, dfs[[genotype]], "B_treatmentDrought_Heat < 1")
temp_df[1,] <- w$hypothesis
temp_df[2,] <- x$hypothesis
temp_df[3,] <- y$hypothesis
temp_df[4,] <- z$hypothesis
temp_df$Genotype <- genotype
hypothesis_tests <- temp_df

for (i in genolist) {
  genotype <- i
  ss <- growthSS(model = "power law", form = area_cm ~ day | replicate / treatment,
                 df = dat[dat$genotype_name == genotype,], start = list("A" = 4, "B" = 2.5), type = "brms", sigma = "spline")
  w <- testGrowth(ss, dfs[[genotype]], "B_treatmentControl < 1")
  x <- testGrowth(ss, dfs[[genotype]], "B_treatmentDrought < 1")
  y <- testGrowth(ss, dfs[[genotype]], "B_treatmentHeat < 1")
  z <- testGrowth(ss, dfs[[genotype]], "B_treatmentDrought_Heat < 1")
  temp_df[1,] <- w$hypothesis
  temp_df[2,] <- x$hypothesis
  temp_df[3,] <- y$hypothesis
  temp_df[4,] <- z$hypothesis
  temp_df$Genotype <- genotype
  hypothesis_tests <- rbind(hypothesis_tests, temp_df)
}


write.csv(hypothesis_tests, "growth-model/B_less_1_test.csv", row.names = F)


#### Test if B > 1
temp_df <- data.frame(Hypothesis = NA,
                      Estimate = NA, 
                      Est.Error = NA, 
                      CI.Lower = NA,
                      CI.Upper = NA,
                      Evid.Ratio = NA, 
                      Post.Prob = NA, 
                      Star = NA)
genolist <- unique(dat$genotype_name)
genolist <- genolist[2:47]

genotype <- "b73"
ss <- growthSS(model = "power law", form = area_cm ~ day | replicate / treatment,
               df = dat[dat$genotype_name == genotype,], start = list("A" = 4, "B" = 2.5), type = "brms", sigma = "spline")
w <- testGrowth(ss, dfs[[genotype]], "B_treatmentControl > 1")
x <- testGrowth(ss, dfs[[genotype]], "B_treatmentDrought > 1")
y <- testGrowth(ss, dfs[[genotype]], "B_treatmentHeat > 1")
z <- testGrowth(ss, dfs[[genotype]], "B_treatmentDrought_Heat > 1")
temp_df[1,] <- w$hypothesis
temp_df[2,] <- x$hypothesis
temp_df[3,] <- y$hypothesis
temp_df[4,] <- z$hypothesis
temp_df$Genotype <- genotype
hypothesis_tests <- temp_df

for (i in genolist) {
  genotype <- i
  ss <- growthSS(model = "power law", form = area_cm ~ day | replicate / treatment,
                 df = dat[dat$genotype_name == genotype,], start = list("A" = 4, "B" = 2.5), type = "brms", sigma = "spline")
  w <- testGrowth(ss, dfs[[genotype]], "B_treatmentControl > 1")
  x <- testGrowth(ss, dfs[[genotype]], "B_treatmentDrought > 1")
  y <- testGrowth(ss, dfs[[genotype]], "B_treatmentHeat > 1")
  z <- testGrowth(ss, dfs[[genotype]], "B_treatmentDrought_Heat > 1")
  temp_df[1,] <- w$hypothesis
  temp_df[2,] <- x$hypothesis
  temp_df[3,] <- y$hypothesis
  temp_df[4,] <- z$hypothesis
  temp_df$Genotype <- genotype
  hypothesis_tests <- rbind(hypothesis_tests, temp_df)
}


write.csv(hypothesis_tests, "growth-model/B_greater_1_test.csv", row.names = F)

