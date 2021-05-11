# Analysis 8-d:
#
# This analysis models predicts unpleasantness ratings as a function of
# stimulation brightness and spectral power density (PSD)
#
# And then assesses whether these effects are moderated by pain ratings
# 
# Level 1:
# unpleasantness ~ 1 + stim + power 
#
# Level 2:
# 1 ~ 1 + mens_pain + fu_pain + bsi_total
# stim ~ 1 + mens_pain + fu_pain + bsi_total
# power ~ 1 + mens_pain + fu_pain + bsi_total
#
# Matt Kmiecik
# Started 8 May 2020

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

# Taking a look at the broadband spectral results first ----

# Subject-wise broadband data
bb_ss <- 
  stim_bl_long %>%
  filter(
    ss %in% eeg_ss_to_keep$ss, # Good EEG data
    ss %in% ss_group_incl$ss,  # DYS, DYSB, and HC diagnosis
    stim > 0
    )

# Summarised broadband data
bb_sum <- 
  bb_ss %>% 
  group_by(elec, freq, stim) %>% 
  summarise(m = mean(dB), sd = sd(dB), n = n(), sem = sd/sqrt(n)) %>%
  ungroup()

# Summarised across stim
bb_ss_stim <-
  bb_ss %>%
  group_by(elec, freq, ss) %>%
  summarise(m = mean(dB), n = n()) %>%
  ungroup()

# Summarised across subjects
# CI calculation from here: https://www.cyclismo.org/tutorial/R/confidence.html#calculating-a-confidence-interval-from-a-t-distribution
bb_sum_sum <-
  bb_ss_stim %>%
  group_by(elec, freq) %>%
  summarise(
    M = mean(m), 
    SD = sd(m), 
    N = n(), 
    SEM = SD/sqrt(N),
    error = qt(.975, df = N-1)*SD/sqrt(N),
    lower = M - error,
    upper = M + error
    ) %>%
  ungroup()

# Version to play around with
ggplot(bb_sum_sum %>% filter(elec == "Oz"), aes(freq, M)) +
  #geom_ribbon(aes(ymin = M-SEM, ymax = M+SEM), fill = "grey", alpha = 2/3) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey", alpha = 2/3) +
  geom_line(size = 1) + 
  coord_cartesian(xlim = c(0, 50), ylim = c(-45, -20)) +
  scale_x_continuous(breaks = seq(0, 50, 5), minor_breaks = NULL) +
  scale_y_continuous(minor_breaks = NULL) +
  labs(x = "Frequency (Hz)", y = "Log Power Spectral Density", caption = "Shading is SEM") +
  theme_classic()

# Final version (SHADING IS 95% CI)
bb_spect_fig <- 
  ggplot(bb_sum_sum %>% filter(elec == "Oz"), aes(freq, M)) +
  geom_vline(size = .6, xintercept = 25, linetype = 2, color = rdgy_pal[3]) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = rdgy_pal[7], alpha = 1) +
  geom_line(size = 1) + 
  coord_cartesian(xlim = c(0, 30), ylim = c(-45, -20)) +
  scale_x_continuous(breaks = seq(0, 30, 5), minor_breaks = NULL) +
  scale_y_continuous(minor_breaks = NULL) +
  labs(x = "Frequency (Hz)", y = "Log Power Spectral Density (dB) \n") +
  theme_classic() + 
  fig_theme
# uncomment to save out image
# ggsave(
#   plot = bb_spect_fig, 
#   path = "figs/publication/", 
#   filename = "bb-spect-fig.svg", 
#   width = 4, 
#   height = 3, 
#   unit = "in"
#   )

# Main data structure to work with ----
spec_data <- 
  stim_bl_long %>%
  filter(
    freq == 25, # 25 Hz spectral power
    ss %in% eeg_ss_to_keep$ss, # Good EEG data
    ss %in% ss_group_incl$ss,  # DYS, DYSB, and HC diagnosis
    stim > 0
  ) %>% # 25Hz, no baseline
  mutate(stim = as.numeric(stim))

# Quick proof that cortical excitability increases with stimulation, esp at Oz ----

# Prepares lvl1 data
lvl1_data_spec <- 
  spec_data %>%
  group_by(elec, ss) %>%
  mutate(stim_mc = as.numeric(scale(stim, scale = FALSE))) %>% # mean centering
  ungroup()

# lvl1 modeling
lvl1_mod_spec <- 
  lvl1_data_spec %>%
  nest_by(elec, ss) %>%
  mutate(modA = list(lm(dB ~ 1 + stim_mc, data = data)))

# lvl1 estimates
lvl1_est_spec <- 
  lvl1_mod_spec %>%
  summarise(broom::tidy(modA)) %>%
  ungroup() %>%
  mutate(term = gsub("[\\(\\)]", "", term))

# lvl2 modeling
lvl2_mod_spec <- 
  lvl1_est_spec %>%
  nest_by(elec, term) %>%
  mutate(
    modA = list(lm(estimate ~ 1, data = data)), 
    modC = list(lm(estimate ~ 0, data = data))
  )

# level 2 estimates organized for visualization
lvl2_est_spec_modA <-
  lvl2_mod_spec %>%
  summarise(broom::tidy(modA, conf.int = TRUE, conf.lvl = .95)) %>%
  ungroup() %>%
  mutate(term = gsub("[\\(\\)]", "", term)) %>%
  mutate(source = rep(lvl2_mod_spec$term, each = 1)) %>%
  select(elec, source, term:conf.high) %>%
  group_by(source, term) %>%
  mutate(
    p.fdr = p.adjust(p.value, method = "fdr"), 
    sig = p.value < .05,
    sig.fdr = p.fdr < .05
    ) %>%
  ungroup() %>%
  left_join(., elec_locs, by = c("elec" = "labels"))

# for publication table
# calculating sums of squares
lvl2_mod_spec_SS <- 
  lvl2_mod_spec %>%
  ungroup() %>%
  group_by(elec, term) %>%
  do(as.data.frame(t(unlist(modelCompare(.$modC[[1]], .$modA[[1]]))))) %>%
  ungroup() %>%
  rename(source = term) # helps with join

# joins sums of squares with estimates to for table for publication
spec_res_table <- 
  left_join(lvl2_est_spec_modA, lvl2_mod_spec_SS, by = c("elec", "source")) %>%
  mutate(
    SS = sseC - sseA, 
    MS = SS/nDF, 
    MSE = sseA/dDF,
    source = ifelse(source == "stim_mc", "Brightness", source)
    ) %>%
  select(
    Electrode = elec, 
    Source = source, 
    Term = term, 
    b = estimate,
    ll = conf.low,
    ul = conf.high,
    SE = std.error,
    F = Fstat,
    p,
    p.fdr,
    SS, 
    SSE = sseA, 
    nDF, 
    dDF, 
    MS, 
    MSE,
    PRE
    )
# write_csv(spec_res_table, path = "doc/viz/spec-res-table.csv") # writes to csv


# Estimates point plot
ggplot(lvl2_est_spec_modA, aes(estimate, elec, color = interaction(sig, sig.fdr))) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_point() +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = .3) +
  theme_minimal() +
  labs(x = "Estimate", y = "Electrode", caption = "95% CI error bars.") +
  facet_wrap(~source, scales = "free") +
  theme(legend.position = "bottom")

# Interpolates the MLM estimates for topo plots
mlm_interp_spec <- 
  lvl2_est_spec_modA %>%
  split(interaction(.$source, .$term, sep = "-")) %>%
  map_dfr(~ topo_interp(data = .x, dv = "estimate", gridRes = 67), .id = "band") %>%
  as_tibble(.) %>%
  separate(band, into = c("source", "term"), sep = "-")


################################################################################
# Settings for topo plots to make them uniform 
maskRing <- circleFun(diameter = 1.42) # creates the mask
contour_alpha <- 1/3
contour_color <- "black"
headshape_size <- .25 # used to be .5
electrode_size <- 1.25
nose_size <- .25 # used to be .5
bwidth <- .5 # width of colorbar
bheight <- .1 # height of colorbar
################################################################################

# INTERCEPT - INTERCEPT TOPO PLOT
terms <- c("Intercept")
this_interp <- mlm_interp_spec %>% filter(source == "Intercept", term %in% terms)  
this_orig <- lvl2_est_spec_modA %>% filter(source == "Intercept", term %in% terms) 
this_min <- min(this_orig$estimate) # minumum for legend
this_max <- max(this_orig$estimate) # maximum for legend
minmax <- max(abs(this_min), abs(this_max))

spec_topo_int <-
  ggplot(this_interp, aes(x = x, y = y, fill = estimate)) +
  coord_equal() + # equalizes coordantes
  geom_raster(interpolate = TRUE) + # basis of the topo
  stat_contour(aes(z = estimate), colour = contour_color, alpha = contour_alpha) + # contour lines
  theme_topo() + # topo theme is added (white background etc.)
  # plots headshape
  geom_path(data = headShape, aes(x, y, z = NULL, fill = NULL), size = headshape_size) +
  geom_point(data = this_orig, aes(x, y, shape = sig.fdr), size = electrode_size) + # plots elecs
  scale_shape_manual(values = c(19)) + # all were sig
  # plots nose
  geom_path(data = nose, aes(x, y, z = NULL, fill = NULL), size = nose_size) +
  # creates a mask to remove points outside defined circle
  geom_path(
    data = maskRing,
    aes(x, y, z = NULL, fill = NULL),
    colour = "white",
    size = 6
  ) +
  # colors here
  # note: oob = squish forces everything outside the colour limits to equal
  # nearest colour boundary (i.e., below min colours = min colour)
  scale_fill_gradientn(
    colours = brewer.pal(n = 9, "Purples"), # jet_colors(10)
    # may want to hard code these for comparison across groups
    limits = c(-36, -24), # these should be determined from the uninterpolated data (i think)
    breaks = c(-36, -30, -24), labels = c(-36, -30, -24),
    guide = "colourbar",
    oob = squish,
    name = "Int"
  ) + 
  guides(
    shape = FALSE, 
    fill = guide_colourbar(
      title.position = "top", 
      title.hjust = 0.5, 
      frame.colour = "black", 
      ticks.colour = "black", 
      barwidth = unit(bwidth, "in"),
      barheight = unit(bheight, "in")
      )
    ) +
  theme(legend.position = "bottom") +
  fig_theme
spec_topo_int

# STIM - INTERCEPT TOPO PLOT
terms <- c("Intercept")
this_interp <- mlm_interp_spec %>% filter(source == "stim_mc", term %in% terms)  
this_orig <- lvl2_est_spec_modA %>% filter(source == "stim_mc", term %in% terms) 
this_min <- min(this_orig$estimate) # minumum for legend
this_max <- max(this_orig$estimate) # maximum for legend
minmax <- max(abs(this_min), abs(this_max))

spec_topo_stim <-
  ggplot(this_interp, aes(x = x, y = y, fill = estimate)) +
  coord_equal() + # equalizes coordantes
  geom_raster(interpolate = TRUE) + # basis of the topo
  stat_contour(aes(z = estimate), colour = "black", alpha = contour_alpha) + # contour lines
  theme_topo() + # topo theme is added (white background etc.)
  # plots headshape
  geom_path(data = headShape, aes(x, y, z = NULL, fill = NULL), size = headshape_size) +
  geom_point(data = this_orig, aes(x, y, shape = sig.fdr), size = electrode_size) + # plots elecs
  scale_shape_manual(values = c(19)) + # all were sig
  # plots nose
  geom_path(data = nose, aes(x, y, z = NULL, fill = NULL), size = nose_size) +
  # creates a mask to remove points outside defined circle
  geom_path(
    data = maskRing,
    aes(x, y, z = NULL, fill = NULL),
    colour = "white",
    size = 6
  ) +
  # colors here
  # note: oob = squish forces everything outside the colour limits to equal
  # nearest colour boundary (i.e., below min colours = min colour)
  scale_fill_gradientn(
    colours = rev(brewer.pal(n = 11, "RdGy")), # jet_colors(10)
    # may want to hard code these for comparison across groups
    limits = c(-2, +2), # customized after looking at min and max defined above
    breaks = c(-2, 0, +2), labels = c(-2, 0, +2),
    guide = "colourbar",
    oob = squish,
    name = "Stim"
  ) + 
  guides(
    shape = FALSE, 
    fill = guide_colourbar(
      title.position = "top", 
      title.hjust = 0.5, 
      frame.colour = "black", 
      ticks.colour = "black", 
      barwidth = unit(bwidth, "in"),
      barheight = unit(bheight, "in")
    )
  ) +
  theme(legend.position = "bottom") +
  fig_theme

# Puts plots together using patchwork
spec_stim_plot <- bb_spect_fig | spec_topo_int | spec_topo_stim + plot_layout(heights = 1)
# #Standard PDF save
# ggsave(
#   filename = "spec-stim-plot-2.pdf",
#   plot = spec_stim_plot,
#   path = "doc/viz/",
#   width = 6.5, height = 5, units = "in"
#   )
# # embeds fonts
# embed_fonts(
#   "doc/viz/spec-stim-plot-2.pdf",
#   outfile = "doc/viz/spec-stim-plot-2-embed.pdf"
#   )

# MAIN ANALYSIS ----

# Preparing data for modeling
spec_ratings_data <- 
  spec_data %>% 
  left_join(., vis_behav_data %>% filter(dv == "rating"), by = c("ss", "stim")) %>% # add unpleas ratings
  group_by(elec, ss) %>%
  mutate(
    dB_mc = scale(dB, scale = FALSE), 
    stim_mc = scale(stim, scale = FALSE)
  ) %>%
  ungroup()

# Level 1 mods
lvl1_mod <- 
  spec_ratings_data %>%
  nest_by(elec, ss) %>%
  mutate(mod = list(lm(value ~ 1 + stim_mc + dB_mc, data = data)))

# Level 1 estimates
lvl1_est <- 
  lvl1_mod %>%
  summarise(broom::tidy(mod)) %>%
  ungroup() %>%
  mutate(
    term = case_when(
      term == "(Intercept)" ~ "Intercept",
      TRUE ~ as.character(term)
    )
  )

# Level 2 data preparation

# Simplifies assessment measures
subs <- unique(lvl1_est$ss) # subjects to extract from med data
med_assess_simple <- 
  med_assess %>% 
  select(ss = subject_id, fupain = bt10b_fupain, bsi1:bsi6) %>%
  mutate(bsi_total = bsi1+bsi2+bsi3+bsi4+bsi5+bsi6) %>%
  filter(ss %in% subs)

# Simplifies screen data measures
med_screen_simple <- 
  med_screen %>% 
  select(ss = subject_id, mens_pain = mh23, mens_pain_nsaid = mh23a) %>%
  filter(ss %in% subs) %>%
  group_by(ss) %>%
  summarise(
    mens_pain = mean(mens_pain), 
    #mens_pain_nsaid = mean(mens_pain_nsaid),
    n = n()
  ) %>% # ss 139 had 2 scores for some reason (moved to arm2)
  ungroup()

# Simple medical data frame
med_data <- 
  left_join(med_assess_simple, med_screen_simple, by = "ss") %>%
  select(ss, fupain, bsi_total, mens_pain) %>%
  mutate(
    fupain_mc = as.numeric(scale(fupain, scale = FALSE)), 
    bsi_total_mc = as.numeric(scale(bsi_total, scale = FALSE)),
    mens_pain_mc = as.numeric(scale(mens_pain, scale= FALSE))
  )

# Level 2 data
lvl2_data <- lvl1_est %>% left_join(., med_data, by = "ss")

# Level 2 mod
lvl2_mod <- 
  lvl2_data %>% 
  nest_by(elec, term) %>% 
  # models all model comparisons to be compared to the augmented model (mod)
  mutate(
    mod = list(lm(estimate ~ 1 + fupain_mc + mens_pain_mc + bsi_total_mc, data = data)),
    mod_int = list(lm(estimate ~ 0 + fupain_mc + mens_pain_mc + bsi_total_mc, data = data)),
    mod_fupain = list(lm(estimate ~ 1 + mens_pain_mc + bsi_total_mc, data = data)),
    mod_mens = list(lm(estimate ~ 1 + fupain_mc + bsi_total_mc, data = data)),
    mod_bsi = list(lm(estimate ~ 1 + fupain_mc + mens_pain_mc, data = data))
    )

# POWER ANALYSIS
# Looking at predicting PSD ~ brightness
# using partial eta squared of .78 as effect size
modelPower(pc = 1, pa = 2, N = 147, alpha = .05, peta2 = .78) # estimates power
modelPower(pc = 1, pa = 2, alpha = .05, power = .8, peta2 = .051) # estimates sample 

# Level 2 estimates
full_lvl2_est <- 
  lvl2_mod %>%
  summarise(broom::tidy(mod, conf.int = TRUE, conf.level = .95)) %>%
  ungroup() %>%
  mutate(
    term = case_when(
      term == "(Intercept)" ~ "Intercept",
      TRUE ~ as.character(term)
    ),
    sig = p.value < .05,
    source = rep(lvl2_mod$term, each = 4),
    term = fct_relevel(term, c("Intercept", "mens_pain_mc", "fupain_mc", "bsi_total_mc"))
  ) %>%
  select(source, term, elec, estimate:sig) %>%
  group_by(source, term) %>%
  mutate(
    p.fdr = p.adjust(p.value, method = "fdr"), 
    sig.fdr = p.fdr < .05,
    p_cor = interaction(sig, sig.fdr),
    p_cor = case_when(
      p_cor == "TRUE.TRUE" ~ "p < .05 (FDR)",
      p_cor == "FALSE.FALSE" ~ "p > .05",
      p_cor == "TRUE.FALSE" ~ "p < .05 (uncorrected)",
      TRUE ~ as.character(p_cor)
      )
    ) %>%
  ungroup() %>%
  left_join(., elec_locs, by = c("elec" = "labels"))

# For publication table ----
# Calculating sums of squares for all the terms

# Intercept
int_SS <- 
  lvl2_mod %>%
  ungroup() %>%
  group_by(elec, term) %>%
  do(as.data.frame(t(unlist(modelCompare(.$mod_int[[1]], .$mod[[1]]))))) %>%
  ungroup() %>%
  rename(source = term) %>%
  mutate(term = "Intercept") # helps with join

# fupain
fupain_SS <- 
  lvl2_mod %>%
  ungroup() %>%
  group_by(elec, term) %>%
  do(as.data.frame(t(unlist(modelCompare(.$mod_fupain[[1]], .$mod[[1]]))))) %>%
  ungroup() %>%
  rename(source = term) %>%
  mutate(term = "fupain_mc") # helps with join

# menspain 
menspain_SS <- 
  lvl2_mod %>%
  ungroup() %>%
  group_by(elec, term) %>%
  do(as.data.frame(t(unlist(modelCompare(.$mod_mens[[1]], .$mod[[1]]))))) %>%
  ungroup() %>%
  rename(source = term) %>%
  mutate(term = "mens_pain_mc") # helps with join

# bsi
bsi_SS <- 
  lvl2_mod %>%
  ungroup() %>%
  group_by(elec, term) %>%
  do(as.data.frame(t(unlist(modelCompare(.$mod_bsi[[1]], .$mod[[1]]))))) %>%
  ungroup() %>%
  rename(source = term) %>%
  mutate(term = "bsi_total_mc") # helps with join

# combining all SS tables
all_SS <- 
  bind_rows(int_SS, fupain_SS, menspain_SS, bsi_SS) %>%
  select(elec, source, term, sseC, sseA, nDF, dDF, Fstat, p, PRE, DeltaR2)

# combining SS tables with estimates
# Intercept    fupain_mc    mens_pain_mc bsi_total_mc
mlm_res <-
  left_join(full_lvl2_est, all_SS, by = c("elec", "source", "term")) %>%
  mutate(
    SS = sseC - sseA, 
    MS = SS/nDF, 
    MSE = sseA/dDF,
    source = case_when(
      source == "dB_mc" ~ "Power",
      source == "stim_mc" ~ "Brightness",
      TRUE ~ as.character(source)
    ),
    term = case_when(
      term == "fupain_mc" ~ "First Urge Pain",
      term == "mens_pain_mc" ~ "Menstrual Pain",
      term == "bsi_total_mc" ~ "BSI",
      TRUE ~ as.character(term)
    )
    ) %>%
  select(
    Electrode = elec, 
    Source = source, 
    Term = term, 
    b = estimate,
    ll = conf.low,
    ul = conf.high,
    SE = std.error,
    F = Fstat,
    p,
    p.fdr,
    SS, 
    SSE = sseA, 
    nDF, 
    dDF, 
    MS, 
    MSE,
    PRE
    )
# write_csv(mlm_res, path = "doc/viz/mlm-res-table.csv") # writes to csv
  

# Visualizing Intercept results ----
ggplot(
  full_lvl2_est %>% filter(source == "Intercept"), 
  aes(estimate, elec, color = interaction(sig, sig.fdr))
) +
  geom_point() +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 1) +
  geom_vline(xintercept = 0, linetype = 2, color = "black") +
  scale_color_manual(values = c(rdgy_pal[7], rdgy_pal[2])) +
  labs(
    x = "Estimate", 
    y = "Electrode", 
    title = "Intercept Models",
    caption = "95% CI error bars."
  ) +
  theme_minimal() +
  facet_wrap(~term, scales = "free", nrow = 1) +
  theme(legend.position = "bottom")

# Scatter plots of this
this_elec <- "Oz"
this_term <- "Intercept"
lvl2_int_data <- lvl2_data %>% filter(elec == this_elec, term == this_term)
pj <- position_jitter(width = .1)

int_plot1 <- 
  ggplot(
  full_lvl2_est %>% 
    filter(elec == this_elec, source == this_term, term == "Intercept"),
  aes(term, estimate)
  ) +
  geom_point(
    data = lvl2_int_data, 
    aes(y = estimate), 
    shape = 1, 
    position = pj,
    alpha = 1/3
    ) +
  geom_point() +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = .1) +
  labs(x = "Term", y = "Estimate") +
  coord_cartesian(ylim = c(0, 20)) +
  theme_minimal()

int_plot2 <-
  ggplot(
  lvl2_int_data %>% pivot_longer(cols = fupain:mens_pain),
  aes(value, estimate)
  ) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  labs(x = "Participant Rating", y = "Estimate") +
  theme_minimal() +
  facet_wrap(~name, scales = "free")

# Together on same graph
int_plot1 / int_plot2

# Manuscript figure:
bsi_scatter <- 
  ggplot(lvl2_int_data, aes(bsi_total, estimate)) +
  geom_point(alpha = 1/2) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  labs(x = "\n BSI Total", y = "") +
  theme_classic() +
  fig_theme

menspain_scatter <- 
  ggplot(lvl2_int_data, aes(mens_pain, estimate)) +
  geom_point(alpha = 1/2) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  #geom_rug(position = "jitter", color = rdgy_pal[10], length = unit(0.015, "npc")) +
  labs(
    x = "\n Menstrual Pain", 
    y = "Mean Unpleasantness Rating \n"
    ) +
  theme_classic() +
  fig_theme

# fupain_scatter <- 
#   ggplot(lvl2_int_data, aes(fupain, estimate)) +
#   geom_point(alpha = 1/2) +
#   geom_smooth(method = "lm", se = FALSE, color = rdgy_pal[3]) +
#   #geom_rug(position = "jitter", color = rdgy_pal[10], length = unit(0.015, "npc")) +
#   labs(
#     x = "\n Pain at First Urge", 
#     y = ""
#     ) +
#   theme_classic() +
#   fig_theme



# Redoing the fu_pain_scatter as this is not the correct
# partial regression plot
library(car) # use this for the avPlots() function
this_mod <- lvl2_mod %>% filter(elec == "Oz", term == "Intercept") # get the model
part_reg_data <- avPlots(this_mod$mod[[1]])
fupain_resids <- as_tibble(part_reg_data$fupain_mc)

fupain_scatter <- 
  ggplot(fupain_resids, aes(fupain_mc, estimate)) +
  geom_point(alpha = 1/2) +
  geom_smooth(method = "lm", se = FALSE, color = rdgy_pal[3]) +
    labs(
    x = "\n Resid(Pain at First Urge)", 
    y = ""
  ) +
  theme_classic() +
  fig_theme

scatter_plots <- menspain_scatter + fupain_scatter + bsi_scatter
scatter_plots

# Visualizing Stim results ----
ggplot(
  full_lvl2_est %>% filter(source == "stim_mc"), 
  aes(estimate, elec, color = interaction(sig, sig.fdr))
) +
  geom_point() +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 1) +
  geom_vline(xintercept = 0, linetype = 2, color = "black") +
  scale_color_manual(values = c(rdgy_pal[7], rdgy_pal[4], rdgy_pal[2])) +
  labs(
    x = "Estimate", 
    y = "Electrode", 
    title = "Stim Models",
    caption = "95% CI error bars."
  ) +
  theme_minimal() +
  facet_wrap(~term, scales = "free", nrow = 1) +
  theme(legend.position = "bottom")

# Interpolate data for topos
mlm_interp_unpleas <- 
  full_lvl2_est %>%
  split(interaction(.$source, .$term, sep = "-")) %>%
  map_dfr(~ topo_interp(data = .x, dv = "estimate", gridRes = 67), .id = "band") %>%
  as_tibble(.) %>%
  separate(band, into = c("source", "term"), sep = "-")

# STIM - INTERCEPT TOPO PLOT
terms <- c("Intercept")
this_source <- "stim_mc"
this_interp <- mlm_interp_unpleas %>% filter(source == this_source, term %in% terms)  
this_orig <- full_lvl2_est %>% filter(source == this_source, term %in% terms) 
this_min <- min(this_orig$estimate) # minumum for legend
this_max <- max(this_orig$estimate) # maximum for legend
minmax <- max(abs(this_min), abs(this_max))

stim_int_plot <-  
  ggplot(this_interp, aes(x = x, y = y, fill = estimate)) +
  coord_equal() + # equalizes coordantes
  geom_raster(interpolate = TRUE) + # basis of the topo
  stat_contour(aes(z = estimate), colour = "black", alpha = contour_alpha) + # contour lines
  theme_topo() + # topo theme is added (white background etc.)
  # plots headshape
  geom_path(data = headShape, aes(x, y, z = NULL, fill = NULL), size = headshape_size) +
  geom_point(data = this_orig, aes(x, y, shape = interaction(sig, sig.fdr)), size = electrode_size) + # plots elecs
  scale_shape_manual(values = c(19)) + # ALL ARE SIGNIFICANT 
  # plots nose
  geom_path(data = nose, aes(x, y, z = NULL, fill = NULL), size = nose_size) +
  # creates a mask to remove points outside defined circle
  geom_path(
    data = maskRing,
    aes(x, y, z = NULL, fill = NULL),
    colour = "white",
    size = 6
  ) +
  # colors here
  # note: oob = squish forces everything outside the colour limits to equal
  # nearest colour boundary (i.e., below min colours = min colour)
  scale_fill_gradientn(
    colours = rev(brewer.pal(n = 11, "RdGy")), # jet_colors(10)
    # may want to hard code these for comparison across groups
    limits = c(-.70, +.70), # customized after looking at min and max defined above
    breaks = c(-.70, 0, +.70), labels = c(-.70, 0, +.70),
    guide = "colourbar",
    oob = squish,
    name = "INT"
  ) + 
  guides(
    shape = FALSE, 
    fill = guide_colourbar(
      title.position = "top", 
      title.hjust = 0.5, 
      frame.colour = "black", 
      ticks.colour = "black", 
      barwidth = unit(bwidth, "in"),
      barheight = unit(bheight, "in")
      )
  ) +
  theme(legend.position = "bottom") +
  fig_theme
stim_int_plot

# STIM - MENSTRUAL PAIN
terms <- c("mens_pain_mc") # "fupain_mc", "bsi_total_mc"
this_source <- "stim_mc"
this_interp <- mlm_interp_unpleas %>% filter(source == this_source, term %in% terms)  
this_orig <- full_lvl2_est %>% filter(source == this_source, term %in% terms) 
this_min <- min(this_orig$estimate) # minumum for legend
this_max <- max(this_orig$estimate) # maximum for legend
minmax <- max(abs(this_min), abs(this_max))

stim_mens_plot <-  
  ggplot(this_interp, aes(x = x, y = y, fill = estimate)) +
  coord_equal() + # equalizes coordantes
  geom_raster(interpolate = TRUE) + # basis of the topo
  stat_contour(aes(z = estimate), colour = "black", alpha = contour_alpha) + # contour lines
  theme_topo() + # topo theme is added (white background etc.)
  # plots headshape
  geom_path(data = headShape, aes(x, y, z = NULL, fill = NULL), size = headshape_size) +
  geom_point(data = this_orig, aes(x, y, shape = interaction(sig, sig.fdr)), size = electrode_size) + # plots elecs
  scale_shape_manual(values = c(1, 17)) +
  # plots nose
  geom_path(data = nose, aes(x, y, z = NULL, fill = NULL), size = nose_size) +
  # creates a mask to remove points outside defined circle
  geom_path(
    data = maskRing,
    aes(x, y, z = NULL, fill = NULL),
    colour = "white",
    size = 6
  ) +
  # colors here
  # note: oob = squish forces everything outside the colour limits to equal
  # nearest colour boundary (i.e., below min colours = min colour)
  scale_fill_gradientn(
    colours = brewer.pal(n = 11, "BrBG"), # jet_colors(10)
    # may want to hard code these for comparison across groups
    limits = c(-.01, +.01), # customized after looking at min and max defined above
    breaks = c(-.01, 0, +.01), labels = c(-.01, 0, +.01),
    guide = "colourbar",
    oob = squish,
    name = "MENS"
  ) + 
  guides(
    shape = FALSE, 
    fill = guide_colourbar(
      title.position = "top", 
      title.hjust = 0.5, 
      frame.colour = "black", 
      ticks.colour = "black", 
      barwidth = unit(bwidth, "in"),
      barheight = unit(bheight, "in")
      )
  ) +
  theme(legend.position = "bottom") +
  fig_theme
stim_mens_plot

# STIM - FU PAIN
terms <- c("fupain_mc") # "bsi_total_mc"
this_source <- "stim_mc"
this_interp <- mlm_interp_unpleas %>% filter(source == this_source, term %in% terms)  
this_orig <- full_lvl2_est %>% filter(source == this_source, term %in% terms) 
this_min <- min(this_orig$estimate) # minumum for legend
this_max <- max(this_orig$estimate) # maximum for legend
minmax <- max(abs(this_min), abs(this_max))

stim_fupain_plot <-  
  ggplot(this_interp, aes(x = x, y = y, fill = estimate)) +
  coord_equal() + # equalizes coordantes
  geom_raster(interpolate = TRUE) + # basis of the topo
  stat_contour(aes(z = estimate), colour = "black", alpha = contour_alpha) + # contour lines
  theme_topo() + # topo theme is added (white background etc.)
  # plots headshape
  geom_path(data = headShape, aes(x, y, z = NULL, fill = NULL), size = headshape_size) +
  geom_point(data = this_orig, aes(x, y, shape = interaction(sig, sig.fdr)), size = electrode_size) + # plots elecs
  scale_shape_manual(values = c(1)) + # All non sig
  # plots nose
  geom_path(data = nose, aes(x, y, z = NULL, fill = NULL), size = nose_size) +
  # creates a mask to remove points outside defined circle
  geom_path(
    data = maskRing,
    aes(x, y, z = NULL, fill = NULL),
    colour = "white",
    size = 6
  ) +
  # colors here
  # note: oob = squish forces everything outside the colour limits to equal
  # nearest colour boundary (i.e., below min colours = min colour)
  scale_fill_gradientn(
    colours = brewer.pal(n = 11, "PuOr"), # jet_colors(10)
    # may want to hard code these for comparison across groups
    limits = c(-.01, +.01), # customized after looking at min and max defined above
    breaks = c(-.01, 0, +.01), labels = c(-.01, 0, +.01),
    guide = "colourbar",
    oob = squish,
    name = "FUP"
  ) + 
  guides(
    shape = FALSE, 
    fill = guide_colourbar(
      title.position = "top", 
      title.hjust = 0.5, 
      frame.colour = "black", 
      ticks.colour = "black", 
      barwidth = unit(bwidth, "in"),
      barheight = unit(bheight, "in")
      )
  ) +
  theme(legend.position = "bottom") +
  fig_theme
stim_fupain_plot

# STIM - BSI
terms <- c("bsi_total_mc") 
this_source <- "stim_mc"
this_interp <- mlm_interp_unpleas %>% filter(source == this_source, term %in% terms)  
this_orig <- full_lvl2_est %>% filter(source == this_source, term %in% terms) 
this_min <- min(this_orig$estimate) # minumum for legend
this_max <- max(this_orig$estimate) # maximum for legend
minmax <- max(abs(this_min), abs(this_max))

stim_bsi_plot <-  
  ggplot(this_interp, aes(x = x, y = y, fill = estimate)) +
  coord_equal() + # equalizes coordantes
  geom_raster(interpolate = TRUE) + # basis of the topo
  stat_contour(aes(z = estimate), colour = "black", alpha = contour_alpha) + # contour lines
  theme_topo() + # topo theme is added (white background etc.)
  # plots headshape
  geom_path(data = headShape, aes(x, y, z = NULL, fill = NULL), size = headshape_size) +
  geom_point(data = this_orig, aes(x, y, shape = interaction(sig, sig.fdr)), size = electrode_size) + # plots elecs
  scale_shape_manual(values = c(1, 17, 19)) +
  # plots nose
  geom_path(data = nose, aes(x, y, z = NULL, fill = NULL), size = nose_size) +
  # creates a mask to remove points outside defined circle
  geom_path(
    data = maskRing,
    aes(x, y, z = NULL, fill = NULL),
    colour = "white",
    size = 6
  ) +
  # colors here
  # note: oob = squish forces everything outside the colour limits to equal
  # nearest colour boundary (i.e., below min colours = min colour)
  scale_fill_gradientn(
    colours = rev(brewer.pal(n = 11, "RdBu")), # jet_colors(10)
    # may want to hard code these for comparison across groups
    limits = c(-.2, +.2), # customized after looking at min and max defined above
    breaks = c(-.2, 0, +.2), labels = c(-.2, 0, +.2),
    guide = "colourbar",
    oob = squish,
    name = "BSI"
  ) + 
  guides(
    shape = FALSE, 
    fill = guide_colourbar(
      title.position = "top", 
      title.hjust = 0.5, 
      frame.colour = "black", 
      ticks.colour = "black", 
      barwidth = unit(bwidth, "in"),
      barheight = unit(bheight, "in")
      )
  ) +
  theme(legend.position = "bottom") +
  fig_theme
stim_bsi_plot

# Putting stim topos together with patchwork
stim_topos <- 
  stim_int_plot + 
  stim_mens_plot + 
  stim_fupain_plot + 
  stim_bsi_plot + 
  plot_layout(width = 1, nrow = 1)
stim_topos

# Visualizing power results ----
ggplot(
  full_lvl2_est %>% filter(source == "dB_mc"), 
  aes(estimate, elec, color = interaction(sig, sig.fdr))
) +
  geom_point() +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 1) +
  geom_vline(xintercept = 0, linetype = 2, color = "black") +
  scale_color_manual(values = c(rdgy_pal[7], rdgy_pal[4], rdgy_pal[2])) +
  labs(
    x = "Estimate", 
    y = "Electrode", 
    title = "Power Models",
    caption = "95% CI error bars."
  ) +
  theme_minimal() +
  facet_wrap(~term, scales = "free", nrow = 1) +
  theme(legend.position = "bottom")

# POWER - INTERCEPT TOPO PLOT
terms <- c("Intercept")
this_source <- "dB_mc"
this_interp <- mlm_interp_unpleas %>% filter(source == this_source, term %in% terms)  
this_orig <- full_lvl2_est %>% filter(source == this_source, term %in% terms) 
this_min <- min(this_orig$estimate) # minumum for legend
this_max <- max(this_orig$estimate) # maximum for legend
minmax <- max(abs(this_min), abs(this_max))

power_int_plot <-  
  ggplot(this_interp, aes(x = x, y = y, fill = estimate)) +
  coord_equal() + # equalizes coordantes
  geom_raster(interpolate = TRUE) + # basis of the topo
  stat_contour(aes(z = estimate), colour = "black", alpha = contour_alpha) + # contour lines
  theme_topo() + # topo theme is added (white background etc.)
  # plots headshape
  geom_path(data = headShape, aes(x, y, z = NULL, fill = NULL), size = headshape_size) +
  geom_point(data = this_orig, aes(x, y, shape = interaction(sig, sig.fdr)), size = electrode_size) + # plots elecs
  scale_shape_manual(values = c(1, 17, 19)) + # all were sig
  # plots nose
  geom_path(data = nose, aes(x, y, z = NULL, fill = NULL), size = nose_size) +
  # creates a mask to remove points outside defined circle
  geom_path(
    data = maskRing,
    aes(x, y, z = NULL, fill = NULL),
    colour = "white",
    size = 6
  ) +
  # colors here
  # note: oob = squish forces everything outside the colour limits to equal
  # nearest colour boundary (i.e., below min colours = min colour)
  scale_fill_gradientn(
    colours = rev(brewer.pal(n = 11, "RdGy")), # jet_colors(10)
    # may want to hard code these for comparison across groups
    limits = c(-.3, +.3), # customized after looking at min and max defined above
    breaks = c(-.3, 0, +.3), labels = c(-.3, 0, +.3),
    guide = "colourbar",
    oob = squish,
    name = "INT"
  ) + 
  guides(
    shape = FALSE, 
    fill = guide_colourbar(
      title.position = "top", 
      title.hjust = 0.5, 
      frame.colour = "black", 
      ticks.colour = "black",
      barwidth = unit(bwidth, "in"),
      barheight = unit(bheight, "in")
      )
  ) +
  theme(legend.position = "bottom") +
  fig_theme
power_int_plot

# POWER - Pain Scores TOPO PLOTs
terms <- c("mens_pain_mc") # "fupain_mc", "bsi_total_mc"
this_source <- "dB_mc"
this_interp <- mlm_interp_unpleas %>% filter(source == this_source, term %in% terms)  
this_orig <- full_lvl2_est %>% filter(source == this_source, term %in% terms) 
this_min <- min(this_orig$estimate) # minumum for legend
this_max <- max(this_orig$estimate) # maximum for legend
minmax <- max(abs(this_min), abs(this_max))

power_menspain_plot <- 
  ggplot(this_interp, aes(x = x, y = y, fill = estimate)) +
  coord_equal() + # equalizes coordantes
  geom_raster(interpolate = TRUE) + # basis of the topo
  stat_contour(aes(z = estimate), colour = "black", alpha = contour_alpha) + # contour lines
  theme_topo() + # topo theme is added (white background etc.)
  # plots headshape
  geom_path(data = headShape, aes(x, y, z = NULL, fill = NULL), size = headshape_size) +
  geom_point(data = this_orig, aes(x, y, shape = p_cor), size = electrode_size) + # plots elecs
  scale_shape_manual(values = c(1, 3, 19)) + # all were sig
  # plots nose
  geom_path(data = nose, aes(x, y, z = NULL, fill = NULL), size = nose_size) +
  # creates a mask to remove points outside defined circle
  geom_path(
    data = maskRing,
    aes(x, y, z = NULL, fill = NULL),
    colour = "white",
    size = 6
  ) +
  # colors here
  # note: oob = squish forces everything outside the colour limits to equal
  # nearest colour boundary (i.e., below min colours = min colour)
  scale_fill_gradientn(
    colours = brewer.pal(n = 11, "BrBG"), # RdGy
    # may want to hard code these for comparison across groups
    limits = c(-.01, +.01), # customized after looking at min and max defined above
    breaks = c(-.01, 0, +.01), labels = c(-.01, 0, +.01),
    guide = "colourbar",
    oob = squish,
    name = "MENS"
  ) + 
  guides(
    shape = FALSE,
    # shape = guide_legend(
    #   title = "Significance",
    #   title.position = "top",
    #   title.hjust = .5
    #   ), 
    fill = guide_colourbar(
      title.position = "top", 
      title.hjust = 0.5, 
      frame.colour = "black", 
      ticks.colour = "black", 
      barwidth = unit(bwidth, "in"),
      barheight = unit(bheight, "in")
    )
  ) +
  theme(legend.position = "bottom") +
  fig_theme
power_menspain_plot  

# POWER - Pain Scores TOPO PLOTs
terms <- c("fupain_mc") #  "bsi_total_mc" , "mens_pain_mc"
this_source <- "dB_mc"
this_interp <- mlm_interp_unpleas %>% filter(source == this_source, term %in% terms)  
this_orig <- full_lvl2_est %>% filter(source == this_source, term %in% terms) 
this_min <- min(this_orig$estimate) # minumum for legend
this_max <- max(this_orig$estimate) # maximum for legend
minmax <- max(abs(this_min), abs(this_max))

power_fupain_plot <- 
  ggplot(this_interp, aes(x = x, y = y, fill = estimate)) +
  coord_equal() + # equalizes coordantes
  geom_raster(interpolate = TRUE) + # basis of the topo
  stat_contour(aes(z = estimate), colour = "black", alpha = contour_alpha) + # contour lines
  theme_topo() + # topo theme is added (white background etc.)
  # plots headshape
  geom_path(data = headShape, aes(x, y, z = NULL, fill = NULL), size = headshape_size) +
  geom_point(data = this_orig, aes(x, y, shape = p_cor), size = electrode_size) + # plots elecs
  scale_shape_manual(values = c(17, 1)) +
  # plots nose
  geom_path(data = nose, aes(x, y, z = NULL, fill = NULL), size = nose_size) +
  # creates a mask to remove points outside defined circle
  geom_path(
    data = maskRing,
    aes(x, y, z = NULL, fill = NULL),
    colour = "white",
    size = 6
  ) +
  # colors here
  # note: oob = squish forces everything outside the colour limits to equal
  # nearest colour boundary (i.e., below min colours = min colour)
  scale_fill_gradientn(
    colours = brewer.pal(n = 11, "PuOr"), # jet_colors(10)
    # may want to hard code these for comparison across groups
    limits = c(-.02, +.02), # customized after looking at min and max defined above
    breaks = c(-.02, 0, +.02), labels = c(-.02, 0, +.02),
    guide = "colourbar",
    oob = squish,
    name = "FUP"
  ) + 
  guides(
    shape = FALSE,
    # shape = guide_legend(
    #   title = "Significance",
    #   title.position = "top",
    #   title.hjust = .5
    # ), 
    fill = guide_colourbar(
      title.position = "top", 
      title.hjust = 0.5, 
      frame.colour = "black", 
      ticks.colour = "black", 
      barwidth = unit(bwidth, "in"),
      barheight = unit(bheight, "in")
    )
  ) +
  theme(legend.position = "bottom") +
  fig_theme
power_fupain_plot  

# POWER - Pain Scores TOPO PLOTs
terms <- c("bsi_total_mc") #   , "mens_pain_mc" "fupain_mc"
this_source <- "dB_mc"
this_interp <- mlm_interp_unpleas %>% filter(source == this_source, term %in% terms)  
this_orig <- full_lvl2_est %>% filter(source == this_source, term %in% terms) 
this_min <- min(this_orig$estimate) # minumum for legend
this_max <- max(this_orig$estimate) # maximum for legend
minmax <- max(abs(this_min), abs(this_max))

power_bsi_plot <- 
  ggplot(this_interp, aes(x = x, y = y, fill = estimate)) +
  coord_equal() + # equalizes coordantes
  geom_raster(interpolate = TRUE) + # basis of the topo
  stat_contour(aes(z = estimate), colour = "black", alpha = contour_alpha) + # contour lines
  theme_topo() + # topo theme is added (white background etc.)
  # plots headshape
  geom_path(data = headShape, aes(x, y, z = NULL, fill = NULL), size = headshape_size) +
  geom_point(data = this_orig, aes(x, y, shape = p_cor), size = electrode_size) + # plots elecs
  scale_shape_manual(values = c(19, 1)) +
  # plots nose
  geom_path(data = nose, aes(x, y, z = NULL, fill = NULL), size = nose_size) +
  # creates a mask to remove points outside defined circle
  geom_path(
    data = maskRing,
    aes(x, y, z = NULL, fill = NULL),
    colour = "white",
    size = 6
  ) +
  # colors here
  # note: oob = squish forces everything outside the colour limits to equal
  # nearest colour boundary (i.e., below min colours = min colour)
  scale_fill_gradientn(
    colours = rev(brewer.pal(n = 11, "RdBu")), # jet_colors(10)
    # may want to hard code these for comparison across groups
    limits = c(-.2, +.2), # customized after looking at min and max defined above
    breaks = c(-.2, 0, +.2), labels = c(-.2, 0, +.2),
    guide = "colourbar",
    oob = squish,
    name = "BSI"
  ) + 
  guides(
    shape = FALSE,
    # shape = guide_legend(
    #   title = "Significance",
    #   title.position = "top",
    #   title.hjust = .5
    # ), 
    fill = guide_colourbar(
      title.position = "top", 
      title.hjust = 0.5, 
      frame.colour = "black", 
      ticks.colour = "black", 
      barwidth = unit(bwidth, "in"),
      barheight = unit(bheight, "in")
    )
  ) +
  theme(legend.position = "bottom") +
  fig_theme
power_bsi_plot

# Scatterplots of Oz unpleas and stim
spec_ratings_oz <- spec_ratings_data %>% filter(elec == "Oz") # simps data

# STIM
oz_stim_slopes <- 
  ggplot(spec_ratings_oz, aes(stim, value, group = ss)) +
  geom_line(
    stat = "smooth",
    method = "lm", 
    formula = y ~ x,
    se = FALSE, 
    color = rdgy_pal[8], 
    size = .5,
    position = position_dodge(width = .2, preserve = "total"),
    alpha = 1/2
    ) +
  geom_smooth(aes(group = 1), method = "lm", color = rdgy_pal[3], se = FALSE) +
  scale_y_continuous(breaks = seq(0, 20, 5), minor_breaks = NULL) +
  labs(x = "\n Brightness Level", y = "Unpleasantness Rating \n") +
  theme_classic() +
  fig_theme
oz_stim_slopes

# POWER
oz_power_slopes <-
  ggplot(spec_ratings_oz, aes(dB, value, group = ss)) +
  geom_line(
    stat = "smooth",
    method = "lm", 
    formula = y ~ x,
    se = FALSE, 
    color = rdgy_pal[8], 
    size = .5,
    position = position_dodge(width = .2, preserve = "total"),
    alpha = 1/2
    ) +
  geom_smooth(aes(group = 1), method = "lm", color = rdgy_pal[3], se = FALSE) +
  scale_y_continuous(breaks = seq(0, 20, 5), minor_breaks = NULL) +
  labs(x = "\n PSD", y = "Unpleasantness Rating \n") +
  theme_classic() +
  fig_theme
oz_power_slopes
oz_stim_slopes / oz_power_slopes

# Putting power topos together with patchwork
power_topos <- 
  power_int_plot + 
  power_menspain_plot + 
  power_fupain_plot + 
  power_bsi_plot + 
  plot_layout(width = 1, nrow = 1)

# Putting all plots together
# main_figure <- scatter_plots / stim_topos / power_topos
# main_figure

# New top panel
top_pan <- (fupain_scatter/plot_spacer()) | plot_spacer() | (oz_stim_slopes / oz_power_slopes)
# # SVG save
# ggsave(
#   filename = "top-panel-v2.svg",
#   plot = top_pan,
#   path = "doc/viz/",
#   width = 6.5, height = 4, units = "in"
#   )
# Bottom panel
bottom_pan <-
  scatter_plots /
  #(plot_spacer() + plot_spacer() + plot_spacer()) /
  stim_topos /
  power_topos


# # Standard PDF save
# ggsave(
#   filename = "main-figure.pdf",
#   plot = main_figure,
#   path = "doc/viz/",
#   width = 6.5, height = 7, units = "in"
#   )
# # embeds fonts
# embed_fonts(
#   "doc/viz/main-figure.pdf",
#   outfile = "doc/viz/main-figure-embed.pdf"
#   )

# # Standard PDF save
# ggsave(
#   filename = "bottom-panel.pdf",
#   plot = bottom_pan,
#   path = "doc/viz/",
#   width = 6.5, height = 7, units = "in"
#   )
# # embeds fonts
# embed_fonts(
#   "doc/viz/bottom-panel.pdf",
#   outfile = "doc/viz/bottom-panel-embed.pdf"
#   )

# Checking the trial counts for paper ----
load("../output/trial-counts.RData")

tcount_data <- 
  trial_counts %>% 
  # filters out unused subjects and stims
  filter(ss %in% unique(lvl1_data_spec$ss), stim > 0)

tcount_sum <- 
  tcount_data %>%
  group_by(stim) %>%
  summarise(sum = sum(trials), n = n(), poss = n*20, perc = sum/poss)

# Percentage of epochs excluded (< 1%)
100 - (sum(tcount_sum$sum)/sum(tcount_sum$poss))*100


# DEMOGRAPHICS TABLE ----
# med_assess
# med_screen

# mh20 - 20. How many days do you have menstrual pelvic pain in an average month?
# mh21 - 21. How many days of work, school, or your usual daily activities have you missed because of your painful periods IN THE LAST 3 MONTHS?
# mh27 - 27. What is the average number of days you bleed with each period (please count even very light bleeding)?

# for people that have menstrual pain, do they have period pain?
# mh19 - 19. Did this start at the time of menarche (first menses)?
# mh19a - 19a. If not, how old were you when your painful periods started?
# mh19b - 19b. What is the total number of years since menarche (first period) that you went without a period? For example, due to continuous birth control use, lactation, pregnancy, or other reasons.

# report the % of people that had pain since menarche
# and y% had menstrual that began after menarche for median 25-75 CI # of years

demo_data <- 
  med_screen %>%
  select(
    ss = subject_id, 
    age = mh2_age, 
    race_am_in = mh3_race___1, # American Indian or Alaskan Native
    race_asian = mh3_race___2, # Asian
    race_pac_isl = mh3_race___3, # Native Hawaiian or Pacific Islander
    race_aa = mh3_race___4, # Black or African American
    race_white = mh3_race___5, # White
    ethn = mh4_ethnicity,
    mh18_painfulperiodsyn,
    mh20, 
    mh21, 
    mh27, 
    mh19, 
    mh19a, 
    mh19b
    ) %>%
  filter(ss %in% unique(lvl1_data_spec$ss)) %>%
  # ss139 had duplicate entries with one of the entries with NA
  mutate(mh27 = ifelse(ss == 139 & is.na(mh27), 5, mh27)) %>%
  # ss24 has NA for mh18_painfulperiodsys but did complete mh19 qs
  mutate(
    mh18_painfulperiodsyn = ifelse(
      ss == 24 & is.na(mh18_painfulperiodsyn), 
      1, 
      mh18_painfulperiodsyn
      )
    ) %>%
  distinct() %>% # removes the duplicate entries for ss139
  # Participant 24 had missing age and it was confirmed via paperwork that her
  # age was 19 years old; therefore, this is corrected here:
  mutate(age = ifelse(is.na(age) & ss == 24, 19, age)) %>%
  # Participant reported bleeding for 5 days (mh27) in paper data; corrected here:
  mutate(mh27 = ifelse(is.na(mh27) & ss == 24, 5, mh27))

# Age
demo_data %>% 
  select(age) %>% 
  summarise(
    m = mean(age, na.rm = TRUE), 
    sd = sd(age, na.rm = TRUE),
    n = n(),
    n.na = sum(is.na(age)),
    min = min(age, na.rm = TRUE),
    max = max(age, na.rm = TRUE)
    )

race_data <- 
  demo_data %>% 
  select(ss, starts_with("race")) %>%
  mutate(mult = race_am_in + race_asian + race_pac_isl + race_aa + race_white)

# Identified with a single race
race_data %>% 
  filter(mult == 1) %>%
  select(-mult) %>%
  pivot_longer(cols = -ss) %>%
  filter(value > 0) %>%
  count(name)

# Identified with multiple races (13 subjects)
race_data %>% filter(mult > 1)

# Did not provide a race
race_data %>% filter(mult < 1)

# Ethnicity
# Hispanic or Latino = 1
# Not Hispanic or Latino = 2
demo_data %>% select(ss, ethn) %>% count(ethn)

# Do you have painful periods (mh18_painfulperiodsys)
demo_data %>% 
  select(ss, mh18_painfulperiodsyn) %>%
  count(mh18_painfulperiodsyn)

demo_data %>% filter(is.na(mh18_painfulperiodsyn)) # HERE 

# Temporal questions (mh20, mh21, mh27)
demo_data %>% 
  select(ss, mh20:mh27) %>%
  pivot_longer(-ss) %>%
  group_by(name) %>%
  summarise(
    m = mean(value, na.rm = TRUE), 
    sd = sd(value, na.rm = TRUE),
    n = n(),
    n.na = sum(is.na(value)),
    min = min(value, na.rm = TRUE),
    max = max(value, na.rm = TRUE)
  )

# Menarche questions
demo_data %>% select(ss, mh19:mh19b)

demo_data %>% select(ss, mh19) %>% count(mh19)

demo_data %>% 
  select(ss, mh19a) %>%
  summarise(
    m = mean(mh19a, na.rm = TRUE), 
    sd = sd(mh19a, na.rm = TRUE),
    n = n(),
    n.na = sum(is.na(mh19a)),
    min = min(mh19a, na.rm = TRUE),
    max = max(mh19a, na.rm = TRUE)
  )

# mh19b
# 999 Never, always had regular period
# 0 Less than 1 year
# 1 1 year
# 2 2 years
# 11 More than 10 years
demo_data %>%
  select(ss, mh19b) %>%
  count(mh19b)

# # 31 ss consistently did not respond to mh questions
# test <- demo_data %>% filter(is.na(mh19))
# write_csv(x = test, path = "missing-mh.csv")

# Partial regression plots plots for fupain and mens pain at Oz ----

# Model to compute resids from:
db_model <- lvl2_mod %>% filter(elec == "Oz", term == "dB_mc")
stim_model <- lvl2_mod %>% filter(elec == "Oz", term == "stim_mc")

# computing partial regression residuals for plotting
db_resids <- avPlots(db_model$mod[[1]])
stim_resids <- avPlots(stim_model$mod[[1]])

# gathering specific resids for paper
db_resids_fupain <- as_tibble(db_resids$fupain_mc)
stim_resids_menspain <- as_tibble(stim_resids$mens_pain_mc)

# plotting partial regressions (added variable plots)

# PSD slope ~ pain at first urge 
db_fupain <- 
  ggplot(db_resids_fupain, aes(fupain_mc, estimate)) +
  geom_point(alpha = 1/2) +
  geom_smooth(method = "lm", se = FALSE, color = rdgy_pal[3]) +
  labs(
    x = "\n Residuals (fup)", 
    y = "Residuals (psdSlope ~ fup)"
  ) +
  theme_classic() +
  fig_theme

# Brightness slope ~ pain at first urge
stim_menspain <- 
  ggplot(stim_resids_menspain, aes(mens_pain_mc, estimate)) +
  geom_point(alpha = 1/2) +
  geom_smooth(method = "lm", se = FALSE, color = rdgy_pal[3]) +
  labs(
    x = "\n Residuals (menspain)", 
    y = "Residuals (brightSlope ~ menspain)"
  ) +
  theme_classic() +
  fig_theme

# stacks and saves for paper
av_plots_paper <- stim_menspain / db_fupain
# ggsave(
#   plot = av_plots_paper,
#   path = "doc/viz/",
#   filename = "av-plots-paper-v2.svg",
#   width = 2,
#   height = 4,
#   unit = "in"
#   )

# Density plots for Menstrual pain, BSI, and First urge pain
density_data <- 
  filter(lvl2_data, elec == "Oz", term == "Intercept") %>%  # just to get sample data
  select(ss, fupain, bsi_total, mens_pain) #%>%
  #pivot_longer(-ss)

density_data_long <- 
  filter(lvl2_data, elec == "Oz", term == "Intercept") %>%  # just to get sample data
  select(ss, fupain, bsi_total, mens_pain) %>%
  pivot_longer(-ss)

# Summary stats for these moderating variables
density_data_long %>% 
  group_by(name) %>%
  summarise(
    m = mean(value), 
    sd = sd(value), 
    n = n(), 
    sem = sd/sqrt(n), 
    min = min(value),
    max = max(value)
    ) %>%
  ungroup()

cor_res <-
  psych::corr.test(
  density_data %>% select(-ss), 
  method = "pearson", 
  adjust = "none"
  )

# plotting
dark2_pal <- brewer.pal(8, "Dark2")

mp_density <- 
  ggplot(density_data, aes(mens_pain)) +
  geom_density(color = dark2_pal[3], fill = dark2_pal[3], alpha = 1/3) +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, .06)) +
  labs(x = "Menstrual Pain (VAS)", y = "Density") +
  scale_y_continuous(minor_breaks = NULL) +
  scale_x_continuous(minor_breaks = NULL) +
  theme_minimal() +
  fig_theme

fup_density <- 
  ggplot(density_data, aes(fupain)) +
  geom_density(color = dark2_pal[2], fill = dark2_pal[2], alpha = 1/3) +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, .06)) +
  labs(x = "First Urge Pain (VAS)", y = "Density") +
  scale_y_continuous(minor_breaks = NULL) +
  scale_x_continuous(minor_breaks = NULL) +
  theme_minimal() +
  fig_theme

bsi_density <- 
  ggplot(density_data, aes(bsi_total)) +
  geom_density(color = dark2_pal[1], fill = dark2_pal[1], alpha = 1/3) +
  coord_cartesian(xlim = c(0, 15), ylim = c(0, .25)) +
  labs(x = "BSI Total Score", y = "Density") +
  scale_y_continuous(minor_breaks = NULL) +
  scale_x_continuous(minor_breaks = NULL) +
  theme_minimal() +
  fig_theme

other_density <-
  ggplot(
  density_data_long %>% filter(name %in% c("fupain", "mens_pain")), 
  aes(value, group = name, color = name, fill = name)
  ) +
  geom_density(alpha = 1/3) +
  scale_color_manual(values = c(dark2_pal[2], dark2_pal[3])) +
  scale_fill_manual(values = c(dark2_pal[2], dark2_pal[3])) +
  coord_cartesian(ylim = c(0, .06)) +
  labs(x = "Pain Rating (VAS)", y = "Density") +
  scale_y_continuous(minor_breaks = NULL) +
  scale_x_continuous(minor_breaks = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  fig_theme

#density_plot <- bsi_density + other_density # patches the two plots together
# ggsave(
#   plot = density_plot,
#   path = "doc/viz/",
#   filename = "density-plot-v1.svg",
#   width = 3,
#   height = 2.5,
#   unit = "in"
#   )

density_plot_2 <- mp_density + fup_density + bsi_density
# ggsave(
#   plot = density_plot_2,
#   path = "doc/viz/",
#   filename = "density-plot-2-v1.svg",
#   width = 5,
#   height = 3,
#   unit = "in"
#   )