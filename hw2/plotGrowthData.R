# Bio185d
# Joe Wirth, September 2022
# Cevi Bainton and Kenneth Mitchel September 2022
# Run on local RStudio

# import packages and helper code
library(Rmisc)
library(ggplot2)
library(dplyr)
source("/root/my_files/plotHelper.R")

# input file names
raw.fn <- "/root/my_files/data/absorbance_raw.txt"
key.fn <- "/root/my_files/data/plate_key.txt"

# output file names (CHANGE THE PATH IF YOURS IS DIFFERENT!)
pos.fn <- "/root/my_files/plots/positive_control.png"
neg.fn <- "/root/my_files/plots/negative_control.png"
exp.fn <- "/root/my_files/plots/experiment.png"

# import the plate key and the raw data as data.frames

absorbances <- read.delim("./data/absorbance_raw.txt")
plate.key <- read.delim("./data/plate_key.txt")


# remove the unused wells from both data.frames

plate.key <- plate.key %>% 
  filter(!is.na(strain)) # remove all strains where is na

absorbances <- absorbances %>% 
  filter(!is.na(val)) # remove all rows where val is na

# keep only the measurements from every other (even) hours

absorbances <- absorbances %>%
  filter(time %% 2 == 0) # filter: keep all times where mod 2 of time is 0


# add strain and condition data from key data.frame to the raw data.frame

absorbances <-   merge(x = absorbances, y = plate.key, by = "well")


# convert the strain column to a factor with the following levels:
##### wt, mut1, mut2, mut3, mut4


absorbances$strain <- factor(absorbances$strain, c("wt", "mut1", "mut2", "mut3", "mut4"))

# get three separate statistical summaries, one for each condition

summarySE(absorbances, "val", "condition")

# TODO: generate three separate plots, one for each statistical summary
#### x-axis: time (hrs)
#### y-axis: A630

# absorbances <- absorbances %>%
#   mutate(log_val = log10(val))

absorbance_means <-absorbances  %>% 
  group_by(time, strain, condition) %>%
  # summarise(mean_log_val = mean(log_val), se_log_val = sd(log_val) / sqrt(n()))
  summarise(mean_val = mean(val), se_val = sd(val) / sqrt(n()))

  
# scale scale_y_log10

# absorbance_means %>%
#   ggplot(aes(x = time, y = mean_val, colour = strain)) +
#   geom_point() +
#   geom_line() +
#   geom_errorbar(aes(ymin = mean_val - se_val, ymax = mean_val + se_val)) + 
#   facet_grid(`strain` ~ `condition`) + 
#   scale_y_log10() +
#   xlab("Time (hours)") + ylab("A630 (AU)")
# ggsave("Double faceted absorbance.png")

absorbance_means %>% 
  filter(`condition` == 'exp') %>%
  ggplot(aes(x = time, y = mean_val, colour = strain)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = mean_val - se_val, ymax = mean_val + se_val)) + 
  scale_y_log10() +
  xlab("Time (hours)") + ylab("Log 10 A630 (AU)") + ggtitle("Experimental Absorbance")
ggsave("plots/Experimental_absorbance.png")

absorbance_means %>% 
  filter(`condition` == 'pos') %>%
  ggplot(aes(x = time, y = mean_val, colour = strain)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = mean_val - se_val, ymax = mean_val + se_val)) + 
  scale_y_log10() +
  xlab("Time (hours)") + ylab("Log 10 A630 (AU)") + ggtitle("Positive Control Absorbance")
ggsave("plots/Positive_absorbance.png")

absorbance_means %>% 
  filter(`condition` == 'neg') %>%
  ggplot(aes(x = time, y = mean_val, colour = strain)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = mean_val - se_val, ymax = mean_val + se_val)) + 
  scale_y_log10() +
  xlab("Time (hours)") + ylab("Log 10 A630 (AU)") + ggtitle("Negative Control Absorbance")
ggsave("plots/Negative_absorbance.png")
# average and get se for each data point

# save the plots (positive control plot named 'pos.plot' in example below)
#### example: savePlotAsPng(pos.plot, pos.fn)


# extra credit: choose better colors, shapes, and linetypes for the plots