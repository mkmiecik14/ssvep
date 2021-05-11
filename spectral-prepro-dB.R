# Preprocessing script to organize the spectral data into R
# Matt Kmiecik
# Started 24 FEB 2020

# Load packages ----
library(tidyverse); library(R.matlab); library(readxl)

# Loading various data ----

# Loading in spectral results from matlab (.mat)
spec_files <- list.files("../data/5-spec-res/", pattern = "*spectral-res.mat")

# Loading in channel locations
chan_locs <- read_xlsx("../data/0-global-vars/elec_pos.xlsx")

# Gathering subject numbers
subjs <- as.numeric(gsub("[^[:digit:].]", "", spec_files))

# Reading in and unpacking spectral results ----
spec_res <- 
  spec_files %>%
  map(~readMat(file.path("../data/5-spec-res", .x))) %>%
  map("spec.res") %>%
  map(~as.matrix(.x))

# Pulling out various spectral results and trial numbers ----

# Baseline
bl_spectra  <- 
  spec_res %>% 
  map(1) %>% 
  map_dfr(~as_tibble(.x)) %>%
  mutate(
    elec = rep(chan_locs$labels, times = nrow(.)/length(chan_locs$labels)), # electrodes
    ss = rep(subjs, each = 32) # subject numbers
  )
  
bl_freqs    <- spec_res %>% map(2) %>% map_dfr(~as.data.frame(.x))
bl_trials   <- spec_res %>% map(3) %>% map_dfr(~as.data.frame(.x))

# Stimulation intensities
stim_spectra  <- spec_res %>% map(4)

s1_stim <- 
  stim_spectra %>% 
  map(~.x[,,1]) %>% 
  map_dfr(~as_tibble(.x)) %>%
  mutate(
    elec = rep(chan_locs$labels, times = nrow(.)/length(chan_locs$labels)), # electrodes
    ss = rep(subjs, each = 32) # subject numbers
    )

s2_stim <- 
  stim_spectra %>% 
  map(~.x[,,2]) %>% 
  map_dfr(~as_tibble(.x)) %>%
  mutate(
    elec = rep(chan_locs$labels, times = nrow(.)/length(chan_locs$labels)), # electrodes
    ss = rep(subjs, each = 32) # subject numbers
  )

s3_stim <- 
  stim_spectra %>% 
  map(~.x[,,3]) %>% 
  map_dfr(~as_tibble(.x)) %>%
  mutate(
    elec = rep(chan_locs$labels, times = nrow(.)/length(chan_locs$labels)), # electrodes
    ss = rep(subjs, each = 32) # subject numbers
  )

s4_stim <- 
  stim_spectra %>% 
  map(~.x[,,4]) %>% 
  map_dfr(~as_tibble(.x)) %>%
  mutate(
    elec = rep(chan_locs$labels, times = nrow(.)/length(chan_locs$labels)), # electrodes
    ss = rep(subjs, each = 32) # subject numbers
  )

s5_stim <- 
  stim_spectra %>% 
  map(~.x[,,5]) %>% 
  map_dfr(~as_tibble(.x)) %>%
  mutate(
    elec = rep(chan_locs$labels, times = nrow(.)/length(chan_locs$labels)), # electrodes
    ss = rep(subjs, each = 32) # subject numbers
  )
  
# Table that helps convert to freq bins
freq_values <- 
  tibble(
    freq = paste0("V", 1:257), 
    new = seq(0, (length(freq)-1)/2, .5)
  )

# Processing baseline spectral results
bl_spectra_long <- 
  bl_spectra %>% 
  gather(freq, dB, -ss, -elec) %>% 
  as_tibble(.) %>%
  left_join(., freq_values, by = "freq") %>%
  select(-freq) %>% rename(freq = new) %>%
  mutate(stim = "0") # to make later combining easier

# Combines stimulations
stim_comb <- 
  bind_rows(s1_stim, s2_stim, s3_stim, s4_stim, s5_stim, .id = "stim")

# Converts stimulations to long
stim_spectra_long <- 
  stim_comb %>% 
  gather(freq, dB, -stim, -ss, -elec) %>% 
  left_join(., freq_values, by = "freq") %>%
  select(-freq) %>% 
  rename(freq = new)

# Combines baseline with stimulations (THE data)
stim_bl_long <- bind_rows(bl_spectra_long, stim_spectra_long)

# Organizing accepted trial counts ----

# Baseline trial acceptance
bl_trials_data <- 
  bl_trials %>% 
  mutate(trials = V1, subjs = subjs) %>% 
  select(subjs, trials, -V1) %>%
  mutate(stim = "0")

# Stimulation accepted trials
stim_trials <- 
  spec_res %>% 
  map(6) %>% 
  map_dfr(~as.data.frame(.x)) %>%
  mutate(subjs = subjs) %>%
  rename(s1 = V1, s2 = V2, s3 = V3, s4 = V4, s5 = V5) %>%
  select(subjs, s1:s5) %>%
  as_tibble(.)

# Converts to long
stim_trials_long <- 
  stim_trials %>% 
  gather(stim, trials, -subjs) %>%
  mutate(stim = gsub("s", "", stim))

# Combines baseline with stimulation accepted trial counts
trial_counts <-
  bind_rows(bl_trials_data, stim_trials_long) %>% 
  select(ss = subjs, stim, trials) %>% # reorders for readability
  as_tibble(.)

# Saving out data ----
save(trial_counts, file = "../output/trial-counts.RData")
write_csv(trial_counts, path = "../output/trial-counts.csv")

save(stim_bl_long, file = "../output/spectral-data.RData")
write_csv(stim_bl_long, path = "../output/spectral-data.csv")

# Cleans workspace ----
rm(
  bl_freqs, 
  bl_spectra, 
  bl_spectra_long, 
  bl_trials, 
  bl_trials_data, 
  chan_locs, 
  freq_values, 
  s1_stim,
  s2_stim,
  s3_stim,
  s4_stim,
  s5_stim,
  spec_files,
  spec_res,
  stim_bl_long,
  stim_comb,
  stim_spectra,
  stim_spectra_long, 
  stim_trials, 
  stim_trials_long, 
  subjs, 
  trial_counts
  )
