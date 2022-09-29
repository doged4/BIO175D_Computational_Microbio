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

# TODO: delete this
# absorbances %>% 
#   group_by(condition) %>%
#   summarize("count" = n(),  "mean" = mean(val), "sd" = sd(val), "se" = sd(val)/sqrt(n()))

summarySE(absorbances, "val", "condition")

# TODO: generate three separate plots, one for each statistical summary
#### x-axis: time (hrs)
#### y-axis: A630

# use facet wrapping?

absorbances %>%
  # group_by(strain) %>%
  ggplot() +
  geom_point(aes(x = time, y = val, colour = strain)) +
  geom_line(aes(x = time, y = )) + 
  facet_wrap(`strain` ~ `condition`)
 

# save the plots (positive control plot named 'pos.plot' in example below)
#### example: savePlotAsPng(pos.plot, pos.fn)


# extra credit: choose better colors, shapes, and linetypes for the plots