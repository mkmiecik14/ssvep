# Script for preprocessing the visual eprime data files
# Matt Kmiecik
# Started 21 FEB 2020

# Load libraries
library(tidyverse)

# Load data
vis_edata <- 
  read_delim(
    "../data/0-global-vars/visual-eprime-data.txt", # raw data
    delim = "\t", # tab delimited
    guess_max = 2000 # this is necessary due to parsing failures 
  )

# Selecting and cleaning up all the visual eprime data
vis_behav_data <- 
  vis_edata %>%
  select(
    ss = Subject, 
    session = Session, 
    date = SessionDate, 
    time = SessionTime,
    order = Block,
    stim = BlockList,
    rating_stim1 = RateUnplblock1.RESP,
    rating_stim2 = RateUnpl1block2.RESP,
    rating_stim3 = RateUnplblock3.RESP,
    rating_stim4 = RateUnplblock4.RESP,
    rating_stim5 = RateUnplblock5.RESP,
    rt_block1 = RateUnplblock1.RT,
    rt_block2 = RateUnpl1block2.RT,
    rt_block3 = RateUnplblock3.RT,
    rt_block4 = RateUnplblock4.RT,
    rt_block5 = RateUnplblock5.RT
    ) %>%
  mutate(date = as.Date(date, format = "%m-%d-%Y")) %>% # converts to date
  gather(meas, value, -ss, -session, -date, -time, -order, -stim) %>% # long format
  separate(meas, into = c("dv", "when")) %>% # separates out column
  filter(complete.cases(value)) %>% # removes "missing" data
  distinct(.) %>% # removes repeated rows
  arrange(ss, order, dv) # orders for better viewing

# Data Editing ----
# The following edits are to raw data that was either entered incorrectly or for
# some other reason
# Subject 354 rating the first stim block an 11 not a 1 (see notes)
vis_behav_data[vis_behav_data$ss == 354 & vis_behav_data$order == 1,]$value <- 11

# Saves out data
save(vis_behav_data, file = "../output/vis-behav-data.RData")
write_csv(vis_behav_data, path = "../output/vis-behav-data.csv")

# Removes R objects created from script
rm(vis_behav_data, vis_edata)

