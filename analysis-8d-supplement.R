# Analysis 8-d supplement
# CPP Development 1 year later
#
# Matt Kmiecik
# Started 22 April 2021

# Libraries ----
library(tidyverse); library(broom); library(RColorBrewer); library(patchwork)
library(extrafont) # font_import() loadfonts() loadfonts(device = "postscript")
library(lmSupport)

# Visualization tools ----
rdgy_pal <- brewer.pal(n = 11, name = "RdGy") # color palette
source("topo_tools.R") # topo plotting tools

# theme for figures font size and type
fig_theme <- theme(text = element_text(size = 10, family = "ArialMT"))

# Reads in data ----
load("../output/spectral-data.RData")      # spectral data (in dB)
load("../output/ss-codes.RData")           # subject codes
load("../output/vis-subj-info.RData")      # eeg notes
load("../output/vis-behav-data.RData")     # behavioral data
load("../output/med-assess.RData")         # medical data from assessment visit
load("../output/med-screen.RData")         # medical data from screen visit

# These are the ss that have HC, DYS, DYSB classification
ss_group_incl <- ss_codes %>% filter(group %in% c("HC", "DYS", "DYSB"))

# A dataframe to use for filtering EEG subjects to keep
eeg_ss_to_keep <- 
  vis_subj_info %>% 
  select(ss, eeg_status) %>% 
  filter(eeg_status %in% "keep")

# Spectral data
spec_data <- 
  stim_bl_long %>%
  filter(
    freq == 25, # 25 Hz spectral power
    ss %in% eeg_ss_to_keep$ss, # Good EEG data
    ss %in% ss_group_incl$ss,  # DYS, DYSB, and HC diagnosis
    stim > 0
  ) %>% # 25Hz, no baseline
  mutate(stim = as.numeric(stim))

# Spectral data with unpleas ratings
spec_ratings_data <- 
  spec_data %>% 
  left_join(., vis_behav_data %>% filter(dv == "rating"), by = c("ss", "stim")) %>% # add unpleas ratings
  group_by(elec, ss) %>%
  mutate(
    dB_mc = scale(dB, scale = FALSE), 
    stim_mc = scale(stim, scale = FALSE)
  ) %>%
  ungroup()

final_subs <- 
  spec_ratings_data %>% 
  select(ss) %>% 
  distinct() %>%
  left_join(., ss_codes, by = "ss")

final_subs %>% count(group)

# CPP Criteria data
cpp_data <- 
  read_excel(path = "../data/0-global-vars/cpp-dev-1-year.xlsx") %>%
  rename(ss = subject_id) %>% 
  select(-finalgroup) 

# CPP criteria for EEG participants
cpp_data_eeg <- 
  cpp_data %>%
  filter(
    ss %in% eeg_ss_to_keep$ss, # Good EEG data
    ss %in% ss_group_incl$ss   # DYS, DYSB, and HC diagnosis
  ) %>%
  left_join(., ss_codes, by = "ss")

# Results
cpp_data_eeg %>% count(group, aibs)
cpp_data_eeg %>% count(group, aicsi)
cpp_data_eeg %>% count(group, ibsorbps)

# chi-square
cpp_mat <- 
  cpp_data_eeg %>% 
  count(group, ibsorbps) %>%
  pivot_wider(id_cols = group, names_from = ibsorbps, values_from = n) %>%
  rename(no = `0`, yes = `1`)

cpp_mat_data <- 
  cpp_mat %>% 
  select(-group) %>% 
  as.matrix(.)

rownames(cpp_mat_data) <- cpp_mat$group
fisher.test(cpp_mat_data) # testing CPP development 1 year later