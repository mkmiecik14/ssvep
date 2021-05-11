# Preprocessing some medical data
# Matt Kmiecik
# Started 24 April 2020

# Load libraries ----
library(tidyverse)

# Load data ----
# arm1 <- 
#   read_csv(file = "data/0-global-vars/crampp-arm1-medical.csv") %>%
#   filter(complete.cases(subject_id), redcap_event_name == "screen_visit_arm_1") %>%
#   select(-record_number) %>% # removes record number as is diff across arms
#   mutate_at(vars(subject_id, mh37a), as.numeric) %>% # converts to numeric
#   mutate_at(vars(mh29_acu_mos), as.character) # converts to character
# 
# arm2 <- 
#   read_csv(file = "data/0-global-vars/crampp-arm2-medical.csv") %>%
#   filter(complete.cases(subject_id), redcap_event_name == "screen_visit_arm_1") %>%
#   select(-record_number) %>% # removes record number as is diff across arms
#   mutate_at(vars(subject_id, mh37a), as.numeric) %>% # converts to numeric
#   mutate_at(vars(mh29_acu_mos), as.character) # converts to character

# reads in data
arm1 <- read_csv(file = "../data/0-global-vars/crampp-arm1-medical.csv")
arm2 <- read_csv(file = "../data/0-global-vars/crampp-arm2-medical.csv")

# creating a key to link record and subject id for later
arm1_key <-
  arm1 %>% 
  select(record_number, subject_id) %>% 
  filter(complete.cases(.))

arm2_key <-
  arm2 %>% 
  select(record_number, subject_id) %>% 
  filter(complete.cases(.)) %>%
  distinct()

# Separating out arms into screen and assessment visit
arm1_screen <- 
  arm1 %>% 
  filter(redcap_event_name == "screen_visit_arm_1") %>% 
  select(-record_number) %>%
  mutate_at(vars(subject_id, mh37a), as.numeric) %>% # converts to numeric
  mutate_at(vars(mh29_acu_mos), as.character) # converts to character

arm1_assess <- 
  arm1 %>% 
  filter(redcap_event_name == "assessment_visit_1_arm_1") %>%
  select(-subject_id) %>%
  left_join(., arm1_key, by = "record_number") %>%
  mutate_at(vars(bt1_oztoday, subject_id), as.numeric) %>%
  select(-record_number)

arm2_screen <- 
  arm2 %>% 
  filter(redcap_event_name == "screen_visit_arm_1") %>% 
  select(-record_number) %>%
  mutate_at(vars(subject_id, mh37a), as.numeric) %>% # converts to numeric
  mutate_at(vars(mh29_acu_mos), as.character) # converts to character

arm2_assess <- 
  arm2 %>% 
  filter(redcap_event_name == "assessment_visit_1_arm_1") %>%
  mutate_at(vars(bt1_oztoday), as.numeric) %>%
  select(-record_number)

# Removing columns that have all NAs
arm1_screen_clean <- arm1_screen[,colSums(is.na(arm1_screen)) < nrow(arm1_screen)]
arm2_screen_clean <- arm2_screen[,colSums(is.na(arm2_screen)) < nrow(arm2_screen)]
arm1_assess_clean <- arm1_assess[,colSums(is.na(arm1_assess)) < nrow(arm1_assess)]
arm2_assess_clean <- arm2_assess[,colSums(is.na(arm2_assess)) < nrow(arm2_assess)]

# Combining arm 1 and arm2 data
med_screen <- bind_rows(arm1_screen_clean, arm2_screen_clean)
med_assess <- 
  bind_rows(arm1_assess_clean, arm2_assess_clean) %>%
  select(redcap_event_name, subject_id, bt0a_pretest:bsi18_based_questionnaire_complete)

# Fix the scores for a few individuals
# These folks skipped a branching logic question thats why the data is missing

# Subject 24
# Notes from Kevin:
# Participant 24 -- she phone screened as 6, and every year after that has 
# always been mh23: 67,71, 73 and mh23a:43,50,50
# so she is very stable. I would use the average of those numbers to compute 
# mh23 and mh23a
# mh23: 67,71, 73 and mh23a: 43,50,50
mh23_ss24   <- mean(c(67, 71, 73))
mh23a_ss24  <- mean(c(43, 50, 50))

# Subject 289
# She phone screened as a 5. Here diary data suggests that when 
# its bad (without nsaids) she is a 6. With NSAIDS a 5. So mh23:60 mh23a:50
mh23_ss289  <- 60 
mh23a_ss289 <- 50

# Fixes these subjects' data
med_screen <- 
  med_screen %>% 
  mutate(
    mh23 = replace(mh23, subject_id == 24, mh23_ss24),
    mh23a = replace(mh23a, subject_id == 24, mh23a_ss24)
    ) %>%
  mutate(
    mh23 = replace(mh23, subject_id == 289, mh23_ss289),
    mh23a = replace(mh23a, subject_id == 289, mh23a_ss289)
  )

# Saving out data ----

# Medical screen data
save(med_screen, file = "../output/med-screen.RData")    # RDATA
write.csv(med_screen, file = "../output/med-screen.csv") # CSV
# Medical assessment visit data
save(med_assess, file = "../output/med-assess.RData")    # RDATA
write.csv(med_assess, file = "../output/med-assess.csv") # CSV

# Cleans up script objects ----
rm(
  arm1,
  arm1_assess,
  arm1_assess_clean,
  arm1_key,
  arm1_screen,
  arm1_screen_clean,
  arm2,
  arm2_assess,
  arm2_assess_clean,
  arm2_key,
  arm2_screen,
  arm2_screen_clean,
  med_assess,
  med_screen,
  mh23_ss24,
  mh23_ss289,
  mh23a_ss24,
  mh23a_ss289
  )
