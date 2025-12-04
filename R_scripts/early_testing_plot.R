library(tidyverse)

setwd("./")
theme_set(theme_light() + theme(axis.text = element_text(size = 9),
                                axis.title = element_text(size = 10),
                                strip.text = element_text(color = "black")))

h3 <- read.csv("early_testing/h3_tests.csv")

write.csv(h3, "early_testing/supplemental_table_5.csv", row.names = F)


ggplot(h3, aes(x = firstSig)) +
  facet_wrap(~trt) +
  geom_bar() +
  geom_text(aes(label = after_stat(count)), vjust = -0.5, stat = "count", size = 3) +
  scale_x_continuous(breaks = seq(5, 13, 2)) +
  labs(x = "Day of significant hypothesis test relative to control",
       y = "Count")
