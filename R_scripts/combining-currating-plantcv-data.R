library(tidyverse)
library(ggpubr)
library(pheatmap)
library(patchwork)
library(ggridges)
library(pcvr)


setwd("./")

#read in raw data
dat1 <- read.csv("07252022_VIS_SV_MG001_results-single-value-traits.csv")
dat2 <- read.csv("07282022_VIS_SV_MG002_results-single-value-traits.csv")

#read in multi value file and pull out center of mass
dat1.m <- read.csv("raw-data/07252022_VIS_SV_MG001_results-multi-value-traits.csv.gz")
dat1.m <- dat1.m %>% filter(trait == "center_of_mass")
dat1.mw <- dat1.m %>% pivot_wider(id_cols = c("rotation", "timestamp", "barcode", "image"),
                                   names_from = "label",
                                   values_from = "value")
dat1.mw <- dat1.mw[,1:6]
names(dat1.mw)[5] <- "center_of_mass_x"
names(dat1.mw)[6] <- "center_of_mass_y"

dat1 <- left_join(dat1, dat1.mw, by = c("rotation", "timestamp", "barcode", "image"))
rm(dat1.m)
rm(dat1.mw)

dat2.m <- read.csv("raw-data/07282022_VIS_SV_MG002_results-multi-value-traits.csv.gz")
dat2.m <- dat2.m %>% filter(trait == "center_of_mass")
dat2.mw <- dat2.m %>% pivot_wider(id_cols = c("rotation", "timestamp", "barcode", "image"),
                                  names_from = "label",
                                  values_from = "value")
dat2.mw <- dat2.mw[,1:6]
names(dat2.mw)[5] <- "center_of_mass_x"
names(dat2.mw)[6] <- "center_of_mass_y"

dat2 <- left_join(dat2, dat2.mw, by = c("rotation", "timestamp", "barcode", "image"))
rm(dat2.m)
rm(dat2.mw)

###
# read in barcodes file and combine gentoype and treatment information
barcodes1 <- read.csv("metadata-files/run1_barcodes_merged.csv")
names(barcodes1)[1] <- "barcode"
barcodes2 <- read.csv("metadata-files/run2_barcodes_merged.csv")
names(barcodes2)[1] <- "barcode"

dat1 <- left_join(barcodes1, dat1, by = "barcode")
dat2 <- left_join(barcodes2, dat2, by = "barcode")

# calculate experiment days
dat1$date <- str_trunc(dat1$timestamp, 10, c("right"), ellipsis = "")
dat1$date <- as.Date(dat1$date, format = '%Y-%m-%d')
day1 = as.Date('2022-06-08')
dat1$day <- difftime(dat1$date, day1, units = c("days"))
dat1$day <- as.integer(dat1$day)

dat2$date <- str_trunc(dat2$timestamp, 10, c("right"), ellipsis = "")
dat2$date <- as.Date(dat2$date, format = '%Y-%m-%d')
day1 = as.Date('2022-06-23')
dat2$day <- difftime(dat2$date, day1, units = c("days"))
dat2$day <- as.integer(dat2$day)

dat <- rbind(dat1, dat2)

# add subpoulatino data
subpop <- read.csv("metadata-files/genotype-names.csv")

dat.full <- left_join(subpop, dat, by= "Genotype")

dat.full <- dat.full[,-c(24,29,31,33,35,37,39,41,49:51,55,56,59,61,63)]
dat.full <- dat.full[,c(1:15,51,52,16:ncol(dat.full))]

write.csv(dat.full, "raw-data/2025-08-27_full-dataset-no-filtering-annotations-com.csv", row.names = F)

# calculating center of mass from the reference line instead of from the top of the image
dat.full <- read.csv("raw-data/2025-08-27_full-dataset-no-filtering-annotations-com.csv")

reference_line_px <- 1325

dat.full <- dat.full %>% mutate(comy_intuitive_px = reference_line_px - center_of_mass_y, 
                                comy_percent = comy_intuitive_px/height_above_reference_pixels)


dat.full <- drop_na(dat.full)

# converting pixels to cm values
dat.cm <- dat.full %>%
  group_by(barcode, Genotype, Subpopulation, Genotype_name, Treatment, Replicate, day, timestamp)  %>% summarise(
                                                                                                  height_above_reference_cm = max(height_above_reference_pixels*1.2)/39,
                                                                                                  height_above_reference_pixels = max(height_above_reference_pixels),
                                                                                                  area_cm2 = mean(area_pixels*1.44)/1456,
                                                                                                  area_pixels = mean(area_pixels),
                                                                                                  width_pixels = max(width_pixels),
                                                                                                  width_cm = max(width_pixels*1.2)/39,
                                                                                                  height_pixels = max(height_pixels), 
                                                                                                  height_cm = max(height_pixels*1.2)/39,
                                                                                                  hue_circular_mean_degrees = mean(hue_circular_mean_degrees),
                                                                                                  comy_pixels = mean(center_of_mass_y),
                                                                                                  comy_cm = mean(center_of_mass_y*1.2)/39,
                                                                                                  comy_intuitive_pixels = mean(comy_intuitive_px),
                                                                                                  comy_intuitive_cm = mean(comy_intuitive_px*1.2)/39)




write.csv(dat.cm, "raw-data/2025-08-27_no-filtering-cm-converted.csv", row.names = F)

dat_hdh <- dat.cm %>% filter(Treatment == "Control" | Treatment == "Drought" | Treatment == "Heat" | Treatment == "Drought_Heat")
write.csv(dat_hdh, "raw-data/2025-08-27_no-filtering_C-D-H-DH.csv", row.names = F)
