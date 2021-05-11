# Script for preparing subject information
# Matt Kmiecik
# Started 13 March 2020

# Libraries
library(tidyverse); library(readxl)

# Subject group codes
subj_codes_path <- "../data/0-global-vars/codes_better_v8.xlsx" # path

ss_codes <- 
  read_excel(path = subj_codes_path) %>%
  select(
    ss = `Subject ID`, 
    group = FinalGroup, 
    endo = Endo, 
    notes = `reassignment notes`
  )

save(ss_codes, file = "../output/ss-codes.RData") # saves out

# Notes from EEG
eeg_info_path <- "../data/0-global-vars/vis-subj-info.xlsx"

vis_subj_info <- 
  read_excel(path = eeg_info_path, sheet = "all") %>%
  rename(ss = ssid) %>%
  mutate(ss = as.numeric(ss))

save(vis_subj_info, file = "../output/vis-subj-info.RData") # saves out
write_csv(vis_subj_info, path = "../output/vis-subj-info.csv") # saves out

# cleans up script objects
rm(subj_codes_path, ss_codes, eeg_info_path, vis_subj_info)
