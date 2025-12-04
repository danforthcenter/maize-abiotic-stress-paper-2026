library(pcvr) 
library(tidyverse)
library(ggpubr)
library(ggridges)
library(patchwork)

setwd("./")

# read in watering data with pcvr
metadata1 = bw.water(file = "./water-use-efficiency/MG001_E_060722_metadata.json", envKey = "environment")
metadata2 = bw.water(file = "./water-use-efficiency/MG002_E_062222_metadata.json", envKey = "environment")

watering = rbind(metadata1, metadata2)
watering = filter(watering, barcode != "Cc00AB000000" & barcode != "Cc00AA000000")

barcodes1 <- read.csv("/Users/acasto/Documents/Projects/USDA_maize_stress/maize_abiotic_stress/pheno_data/run1_barcodes_merged.csv")
barcodes2 <- read.csv("/Users/acasto/Documents/Projects/USDA_maize_stress/maize_abiotic_stress/pheno_data/run2_barcodes_merged.csv")

barcodes <- rbind(barcodes1, barcodes2)
colnames(barcodes)[1] <- "barcode"

watering <- full_join(barcodes, watering, by = "barcode")
watering <- watering %>% filter(weight_after != -1)

# 
write.csv(watering, "./water-use-efficiency/MG001_MG002_watering_data.csv", row.names = F)


###

###
## pWUE calculation

watering = read.csv("./water-use-efficiency/MG001_MG002_watering_data.csv")
watering$timestamp <- as.POSIXct(watering$timestamp)

dat = read.csv(file = "./pheno_data/2024-11-13_VIS_SV_MG001_MG002_dead-plant-filtered_calcualted-values_D-DH-H-C.csv")

dat = separate(dat, col = timestamp, into = c("date", "time"), remove = T, sep = "T")
dat$time = str_trunc(dat$time, 8, side = "right", ellipsis = "")
dat = unite(dat, col = "timestamp", "date", "time", sep = " ")
dat$timestamp = as.POSIXct(dat$timestamp)


wue = pwue(df = dat, w = watering, pheno = "area_cm2", time = "timestamp", 
           id = "barcode", offset = 0, method = "rate")

sum(is.infinite(wue$pWUE))

wue <- wue %>% filter(barcode != "Ea041AA116281" & barcode != "Ea045AA115240")
wue$replicate <- as.factor(wue$replicate)


ggplot(wue, aes(x = day, y = pWUE, group = treatment))+
  geom_smooth(aes(colour = treatment), se=F, method = "loess")+
  facet_wrap("genotype_name")+
  theme_light() + theme(strip.text = element_text(color = "black"))


## plot of average pWUE 
wue <- read.csv("pheno_data/water-use-efficiency/2025-01-15_pcvr_wue_filtered_total_water.csv")
plotdf <- wue  %>% group_by(genotype_name, treatment, day) %>% summarise(avg_pWUE = mean(pWUE),
                                                                      sd_pWUE = sd(pWUE))

ggplot(plotdf, aes(x = day, y = avg_pWUE, group = treatment))+
  geom_line(aes(colour = treatment))+
  geom_errorbar(aes(ymax = avg_pWUE + sd_pWUE, ymin = avg_pWUE-sd_pWUE, color = treatment), linewidth = 0.2, width = 0.2)+
  facet_wrap("genotype_name")+
  scale_x_continuous(breaks = c(0,1,2,3,4,5,6,7,8,9,10,11,12))+
  #labs(title = "pWUE")+
  theme_light() + theme(strip.text = element_text(color = "black"))


## plots of pWUE, total_water, and pheno_diff in each treatment with all replicates for diagnostic purposes

wue$replicate <- as.factor(wue$replicate)

trt <- "Control"
plotdf <- wue %>% filter(treatment == trt) %>% filter(total_water >1) 

ggplot(plotdf, aes(x = day, y = pWUE, group = replicate))+
  geom_line(aes(colour = replicate))+
  facet_wrap("genotype_name")+
  scale_x_continuous(breaks = c(0,1,2,3,4,5,6,7,8,9,10,11,12))+
  labs(title = paste("pWUE of ", trt ," replicates"))+
  theme_light() + theme(strip.text = element_text(color = "black"))

ggplot(plotdf, aes(x = day, y = total_water, group = replicate))+
  geom_line(aes(colour = replicate))+
  facet_wrap("genotype_name")+
  scale_x_continuous(breaks = c(0,1,2,3,4,5,6,7,8,9,10,11,12))+
  labs(title = "water used of Control replicates")+
  theme_light() + theme(strip.text = element_text(color = "black"))

ggplot(plotdf, aes(x = day, y = pheno_diff, group = replicate))+
  geom_line(aes(colour = replicate))+
  facet_wrap("genotype_name")+
  scale_x_continuous(breaks = c(0,1,2,3,4,5,6,7,8,9,10,11,12))+
  labs(title = "area difference of Drought_Heat replicates")+
  theme_light() + theme(strip.text = element_text(color = "black"))

 ggplot(wue, aes(x = day, y = total_water, group = treatment))+
  geom_smooth(aes(colour = treatment), se=F)+
  #geom_point(aes(color = treatment))+
  scale_x_continuous(breaks = c(0,1,2,3,4,5,6,7,8,9,10,11,12))+
  facet_wrap("genotype_name")+
  theme_light() + theme(strip.text = element_text(color = "black"))

wue = wue[,2:ncol(wue)]

write.csv(wue, "./water-use-efficiency/2025-01-15_pcvr_wue_filtered_total_water.csv", row.names = F)


## plots of water used and change in area. 
wue <- read.csv("pheno_data/water-use-efficiency/2025-01-15_pcvr_wue_filtered_total_water.csv")
plotdf <- wue %>% group_by(treatment, day) %>% summarise(avg_water_used = mean(total_water), sd_water_used = sd(total_water)/sqrt(n()),
                                                                        avg_pheno_diff = mean(pheno_diff), sd_pheno_diff = sd(pheno_diff)/sqrt(n())) #%>% filter(treatment == "Drought") 
water <- ggplot(plotdf, aes(x = day, y = avg_water_used))+
  geom_line(aes(colour = treatment))+
  geom_errorbar(aes(ymax = avg_water_used + sd_water_used, ymin = avg_water_used - sd_water_used, color = treatment), linewidth = 0.2, width = 0.2)+
  #facet_wrap("genotype_name")+
  scale_color_manual(values = c("Control" = "#7CAE00", "Drought"="#00BFC4", "Heat"="#F8766D", "Drought_Heat"="#C77CFF"))+
  scale_x_continuous(breaks = c(0,1,2,3,4,5,6,7,8,9,10,11,12))+
  labs(y = "Water used (g)", x = "Day")+
  theme_light() + theme(strip.text = element_text(color = "black"), 
                        legend.position = "none")

pheno <- ggplot(plotdf, aes(x = day, y = avg_pheno_diff, group = treatment))+
  geom_line(aes(colour = treatment))+
  geom_errorbar(aes(ymax = avg_pheno_diff + sd_pheno_diff, ymin = avg_pheno_diff - sd_pheno_diff, color = treatment), linewidth = 0.2, width = 0.2)+
  #facet_wrap("genotype_name")+
  scale_color_manual(values = c("Control" = "#7CAE00", "Drought"="#00BFC4", "Heat"="#F8766D", "Drought_Heat"="#C77CFF"))+
  scale_x_continuous(breaks = c(0,1,2,3,4,5,6,7,8,9,10,11,12))+
  labs(y = "Change in plant area (cm^2)", x = "Day")+
  theme_light() + theme(strip.text = element_text(color = "black"))

water + pheno
