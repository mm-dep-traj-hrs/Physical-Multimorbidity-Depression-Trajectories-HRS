# HEADER ===
# Title:        Prevalent Physical Multimorbidity Patterns and Depressive Symptom 
#               Trajectories among Community Dwelling U.S. Older Adults
# Author:       Nicholas Bishop, University of Arizona
# Last updated: 05-04-2026
# Description:  Data preparation, multimorbidity classification, and inverse
#               probability of dropout weight (IPW) construction for analyses
#               of depression trajectories among older U.S. adults using HRS
#               waves 11-15 (2012-2020).
# Data sources:
#   - RAND HRS Longitudinal File 2020 (v2): randhrs1992_2020v2.sav
#   - HRS Childhood Health and Family file: childhoodhealthfamily.sav
#   Both files require registration and approval at https://hrsdata.isr.umich.edu/.
#
# R version:    [e.g., 4.5.2]
# Mplus:        Output .dat files are formatted for Mplus via MplusAutomation.


#  R setup ----------------------------------------------------------------

# Install all packages
install.packages("haven")
install.packages("tidyverse")
install.packages("janitor")
install.packages("psych")
install.packages("rio")
install.packages("car")
install.packages("summarytools")
install.packages("Hmisc")
install.packages("labelled")
install.packages("vtable")
install.packages("dineq")
install.packages("clintools")
# install stable version of MplusAutomation
install.packages("C:/Downloads/MplusAutomation_1.1.0.tar.gz", repos = NULL, type = "source")
install.packages("WeightIt")
install.packages("cobalt")
install.packages("sjPlot")
install.packages("purrr")

# Load in all packages
library(haven)
library(tidyverse)
library(janitor)
library(psych)
library(rio)
library(car)
library(summarytools)
library(Hmisc)
library(labelled)
library(vtable)
library(dineq)
library(MplusAutomation)
library(clintools)
library(missForest)
library(doParallel)
library(WeightIt)
library(cobalt)
library(sjPlot)
library(purrr)

# open data ---------------------------------------------------------------
# set to your local data directory
setwd("C:/YOURPATH/")

# read rand hrs data
rand <- read_spss("randhrs1992_2022v1.sav")

# Filter data by non-missing CESD, non-missing chronic conditions, age greater than 65, 
# excluding proxy responses and missing/zero weights
rand2 <- rand %>% 
  dplyr::filter(!is.na(R11CESD) & R11HIBP %in% c(0,1) & R11DIAB %in% c(0,1) & 
                  R11CANCR %in% c(0,1) & R11LUNG %in% c(0,1) & 
                  R11HEART %in% c(0,1) & R11STROK %in% c(0,1) & R11ARTHR %in% c(0,1) &
                  R11AGEY_B >= 65 & R11PROXY == 0 & R11WTRESP > 0) %>% 
  clean_names(., "snake")

nrow(rand2)

# select variables
rand3 <- rand2  %>%
  dplyr::select(
    # constant variables 
    hhidpn, hhid, pn, raracem, rahispan, raedyrs, rabplace, ragender, r12mstat,
    
    # survey design variables
    r11wtresp, raestrat, hhid,
    
    # indicators of response status
    r11iwstat, r12iwstat, r13iwstat, r14iwstat, r15iwstat,
    
    # cesd
    r11cesd, r12cesd, r13cesd, r14cesd, r15cesd, 
    
    #r11age
    r11agey_b, r12agey_b, r13agey_b, r14agey_b, r15agey_b,
    
    # 2012 vars
    r11mstat,
    h11itot, h11atotb, r11bmi, r11cesd, r11hibp, r11diab, 
    r11cancr,r11lung, r11heart, r11strok, r11arthr, r11drink, 
    r11drinkd, r11drinkn, r11vgactx,
    r11smoken, r11doctor, r11doctim, 
  ) %>%
  remove_val_labels()

# childhood and health data
chf <- read_spss("AGGCHLDFH2016A_R.sav") %>%
  clean_names(., "snake") %>%
  dplyr::select(hhid, pn, moeduc, faeduc, rthlthch,
                famfin, movfin, fmfinh, faunem, fjob)

# merge the datasets
df1 <- left_join(rand3, chf, by=c("hhid","pn")) 

sumtable(df1)

# recoding variables ---------------------------------------
df2 <- df1 %>%
  dplyr::mutate(
    # sex/gender
    female = case_when(ragender == 1 ~ 0,
                       ragender == 2 ~ 1),
    # race/eth
    white = case_when(
      raracem == 1 & rahispan == 0 ~ 1,
      raracem %in% c(2, 3) | rahispan == 1 ~ 0,
      TRUE ~ NA_real_),
    
    black = case_when(
      raracem == 2 & rahispan == 0 ~ 1,
      raracem %in% c(1, 3) | rahispan == 1 ~ 0,
      TRUE ~ NA_real_),
    
    other = case_when(
      raracem == 3 & rahispan == 0 ~ 1,
      raracem %in% c(1, 2) | rahispan == 1 ~ 0,
      TRUE ~ NA_real_),
    
    latinx = case_when(
      rahispan == 1 ~ 1,
      rahispan == 0 ~ 0,
      TRUE ~ NA_real_),
    
    # marital status
    mar_12 = case_when (r11mstat >= 1 & r11mstat <= 3 ~ 1,
                        r11mstat >= 4 & r11mstat <= 8 ~ 0,
                        TRUE ~ NA_real_),
    
    mar_f_12 = as.factor(mar_12),
    
    # education
    edu_c = as.factor(case_when(
      raedyrs < 12 ~ 1,
      raedyrs == 12 ~ 2,
      raedyrs >= 13 ~ 3,
      TRUE ~ NA_real_)),
    
    educ_yr_f = as.factor(raedyrs),
    
    educ_c_1 = case_when (edu_c == 1 ~ 1, edu_c %in% c(2,3) ~ 0, TRUE ~ NA_real_),
    educ_c_2 = case_when (edu_c == 2 ~ 1, edu_c %in% c(1,3) ~ 0, TRUE ~ NA_real_),
    educ_c_3 = case_when (edu_c == 3 ~ 1, edu_c %in% c(1,2) ~ 0, TRUE ~ NA_real_),
    
    # weighted household income quartiles
    inc_q_12 = as.factor(ntiles.wtd(h11itot, 4, weights = r11wtresp)),
    
    inc_q1_12 = case_when (inc_q_12 == "1" ~ 1, inc_q_12 %in% c("2","3","4") ~ 0, TRUE ~ NA_real_),
    inc_q2_12 = case_when (inc_q_12 == "2" ~ 1, inc_q_12 %in% c("1","3","4") ~ 0, TRUE ~ NA_real_),
    inc_q3_12 = case_when (inc_q_12 == "3" ~ 1, inc_q_12 %in% c("1","2","4") ~ 0, TRUE ~ NA_real_),
    inc_q4_12 = case_when (inc_q_12 == "4" ~ 1, inc_q_12 %in% c("1","2","3") ~ 0, TRUE ~ NA_real_),
    
    # weighted household wealth quartiles
    wlt_q_12 = as.factor(ntiles.wtd(h11atotb , 4, weights = r11wtresp)),
    
    wlt_q1_12 = case_when (wlt_q_12 == "1" ~ 1, wlt_q_12 %in% c("2","3","4") ~ 0, TRUE ~ NA_real_),
    wlt_q2_12 = case_when (wlt_q_12 == "2" ~ 1, wlt_q_12 %in% c("1","3","4") ~ 0, TRUE ~ NA_real_),
    wlt_q3_12 = case_when (wlt_q_12 == "3" ~ 1, wlt_q_12 %in% c("1","2","4") ~ 0, TRUE ~ NA_real_),
    wlt_q4_12 = case_when (wlt_q_12 == "4" ~ 1, wlt_q_12 %in% c("1","2","3") ~ 0, TRUE ~ NA_real_),
    
    # nativity
    nonus_b = case_when(
      rabplace %in% c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10) ~ 0,
      rabplace == 11 ~ 1,
      TRUE ~ NA_real_),
    
    # bmi dummy variables
    bmi_12un = case_when(
      r11bmi < 18.5 ~ 1,
      r11bmi >= 18.5 ~ 0,
      TRUE ~ NA_real_),
    
    bmi_12n = case_when(
      r11bmi < 18.5 ~ 0,
      between(r11bmi, 18.5, 24.99) ~ 1,
      r11bmi > 24.99 ~ 0,
      TRUE ~ NA_real_),
    
    bmi_12ow = case_when(
      r11bmi < 25 ~ 0,
      between(r11bmi, 25, 29.99) ~ 1,
      r11bmi > 29.99 ~ 0,
      TRUE ~ NA_real_),
    
    bmi_12ob = case_when(
      r11bmi < 30 ~ 0,
      r11bmi >= 30 ~ 1,
      TRUE ~ NA_real_),
    
    # activity levels (keeping only vigorous) 
    vact_12c = case_when(
      r11vgactx %in% c(1,2) ~ 1,
      r11vgactx %in% c(3,4,5) ~ 0,
      TRUE ~ NA_real_),
    
    r11vgactx_f = as.factor(r11vgactx),
    
    # number of doctor's visits (weighted quartiles)
    r11doctim_q = ntiles.wtd(r11doctim, 4, weights = r11wtresp),
    
    r11doctim_q1 = case_when(r11doctim_q == 1 ~ 1,
                             r11doctim_q %in% c(2,3,4) ~ 0,
                             TRUE ~ NA_real_),
    r11doctim_q2 = case_when(r11doctim_q == 2 ~ 1,
                             r11doctim_q %in% c(1,3,4) ~ 0,
                             TRUE ~ NA_real_),
    r11doctim_q3 = case_when(r11doctim_q == 3 ~ 1,
                             r11doctim_q %in% c(1,2,4) ~ 0,
                             TRUE ~ NA_real_),
    r11doctim_q4 = case_when(r11doctim_q == 4 ~ 1,
                             r11doctim_q %in% c(1,2,3) ~ 0,
                             TRUE ~ NA_real_),
    # chronic conditions 
    hbp12_c = case_when(r11hibp == 1 ~ 1,
                        r11hibp == 0 ~ 0,
                        TRUE ~ NA_real_), 
    
    dia12_c = case_when(r11diab == 1 ~ 1,
                        r11diab == 0 ~ 0,
                        TRUE ~ NA_real_),
    
    can12_c = case_when(r11cancr == 1 ~ 1,
                        r11cancr == 0 ~ 0,
                        TRUE ~ NA_real_),
    
    lung12_c = case_when(r11lung == 1 ~ 1,
                         r11lung == 0 ~ 0,
                         TRUE ~ NA_real_),
    
    hrt12_c = case_when(r11heart == 1 ~ 1,
                        r11heart == 0 ~ 0,
                        TRUE ~ NA_real_),  
    
    str12_c = case_when(r11strok == 1 ~ 1,
                        r11strok == 0 ~ 0,
                        TRUE ~ NA_real_),
    
    art12_c = case_when(r11arthr == 1 ~ 1,
                        r11arthr == 0 ~ 0,
                        TRUE ~ NA_real_),
    
    # condition count 
    mcc_12 = rowSums(data.frame(hbp12_c,dia12_c,can12_c,lung12_c,hrt12_c,str12_c,art12_c), na.rm = FALSE),
    
    # childhood adversity
    # parent's education
    # coding missing as 1 following Montez and Hayward 2014
    ch_fedlt8 = case_when(
      faeduc %in% c(0,1,2,3,4,5,6,7,7.5) ~ 1,
      faeduc %in% c(8,8.5,9,10,11,12,13,14,15,16,17) ~ 0,
      TRUE ~ 1),
    
    # coding missing as 1 following Montez and Hayward 2014
    ch_medlt8 = case_when(
      moeduc %in% c(0,1,2,3,4,5,6,7,7.5) ~ 1,
      moeduc %in% c(8,8.5,9,10,11,12,13,14,15,16,17) ~ 0,
      TRUE ~ 1),
    
    # financial status
    ch_poor = case_when(
      famfin == 5 ~ 1,
      famfin %in% c(1,3) ~ 0,
      TRUE ~ NA_real_),
    
    # move due to financial status
    ch_move = case_when(
      movfin == 1 ~ 1,
      movfin == 5 ~ 0,
      TRUE ~ NA_real_),
    
    # family given financial help
    ch_fhelp = case_when(
      fmfinh == 1 ~ 1,
      fmfinh == 5 ~ 0,
      TRUE ~ NA_real_),
    
    # never lived with father
    ch_nlf = case_when(
      faunem == 7 ~ 1,
      faunem %in% c(1,5,6) ~ 0,
      TRUE ~ NA_real_),
    
    # father occupational status
    # armed services coded as missing
    ch_fblcl = case_when(
      fjob == 5 ~ 1,
      fjob %in% c(1,2,3,4) ~ 0,
      TRUE ~ NA_real_),
    
    # index of childhood circumstances
    # imputing missing for parent's educ as Montez Hayward 2014
    ch_circum_n = rowSums(data.frame(is.na(ch_fedlt8), is.na(ch_medlt8), is.na(ch_poor), is.na(ch_move), 
                                     is.na(ch_fhelp), is.na(ch_nlf), is.na(ch_fblcl)), na.rm = TRUE),
    ch_circum_sum = rowSums(data.frame(ch_fedlt8, ch_medlt8, ch_poor, ch_move, ch_fhelp, ch_nlf, ch_fblcl), na.rm = TRUE),
    
    # creating count variable with missing for cases with more than 1 missing item
    ch_circum = case_when(ch_circum_sum %in% c(0,1,2,3,4,5) ~ ch_circum_sum,
                          ch_circum_sum == 6 ~ 5,
                          ch_circum_n <= 5 ~ NA_real_, 
                          TRUE ~ NA_real_),
    
    # create categorical dummy variables
    ch_cir_0 = case_when (ch_circum == 0 ~ 1, ch_circum %in% c(1,2,3,4,5) ~ 0, TRUE ~ NA_real_),
    ch_cir_1 = case_when (ch_circum == 1 ~ 1, ch_circum %in% c(0,2,3,4,5) ~ 0, TRUE ~ NA_real_),
    ch_cir_2 = case_when (ch_circum == 2 ~ 1, ch_circum %in% c(0,1,3,4,5) ~ 0, TRUE ~ NA_real_),
    ch_cir_3 = case_when (ch_circum == 3 ~ 1, ch_circum %in% c(0,1,2,4,5) ~ 0, TRUE ~ NA_real_),
    ch_cir_4 = case_when (ch_circum == 4 ~ 1, ch_circum %in% c(0,1,2,3,5) ~ 0, TRUE ~ NA_real_),
    ch_cir_5 = case_when (ch_circum == 5 ~ 1, ch_circum %in% c(0,1,2,3,4) ~ 0, TRUE ~ NA_real_),
    
    # childhood health
    ch_hlth_p = case_when (rthlthch %in% c(1,2,3) ~ 0,
                           rthlthch %in% c(4,5) ~ 1, 
                           TRUE ~ NA_real_),
    
    # current smoker 
    smoke_12 = r11smoken,
    
    # alcohol use
    r11drnkwk=r11drinkn*r11drinkd,
    drnk_12 = case_when(r11drnkwk == 0 ~ 0, 
                        female == 0 & (r11drnkwk >= 1 & r11drnkwk <= 14) ~ 1,
                        female == 1 & (r11drnkwk >= 1 & r11drnkwk <= 7) ~ 1,
                        female == 0 & (r11drnkwk > 14) ~ 2,
                        female == 1 & (r11drnkwk > 7) ~ 2,
                        TRUE ~ NA_real_), 
    drnk_12_0 = case_when(drnk_12 == 0 ~ 1, drnk_12 %in% c(1,2) ~ 0, TRUE ~ NA_real_),
    drnk_12_1 = case_when(drnk_12 == 1 ~ 1, drnk_12 %in% c(0,2) ~ 0, TRUE ~ NA_real_),
    drnk_12_2 = case_when(drnk_12 == 2 ~ 1, drnk_12 %in% c(0,1) ~ 0, TRUE ~ NA_real_),
    
  ) %>% dplyr::select(r11wtresp, raestrat, hhid, pn, hhidpn,
                      r11cesd, r12cesd, r13cesd, r14cesd, r15cesd,
                      r11agey_b, r12agey_b, r13agey_b, r14agey_b, r15agey_b,
                      r11iwstat, r12iwstat, r13iwstat, r14iwstat, r15iwstat,
                      r11hibp, r11diab, r11cancr, r11lung, r11heart, r11strok, r11arthr,
                      hbp12_c, dia12_c, can12_c, lung12_c, hrt12_c, str12_c, art12_c, mcc_12,   
                      ragender, female, r11doctim, 
                      r11doctim_q, r11doctim_q1, r11doctim_q2, r11doctim_q3, r11doctim_q4, 
                      r11agey_b, white, black, other, latinx,
                      r12mstat, mar_12, mar_f_12,
                      educ_yr_f, edu_c, educ_c_1, educ_c_2, educ_c_3,
                      h11itot, inc_q_12, inc_q1_12, inc_q2_12, inc_q3_12, inc_q4_12, 
                      h11atotb, wlt_q_12, wlt_q1_12, wlt_q2_12, wlt_q3_12, wlt_q4_12,
                      nonus_b, 
                      bmi_12un, bmi_12n, bmi_12ow, bmi_12ob,
                      vact_12c,
                      ch_circum, ch_cir_0, ch_cir_1, ch_cir_2, ch_cir_3, ch_cir_4, ch_cir_5,
                      ch_hlth_p,
                      drnk_12, drnk_12_0, drnk_12_1, drnk_12_2,
                      smoke_12)

# identifying most common physical disease combinations --------
# select only those with one or more condition
mcc_combo_1 <- df2 %>% 
  dplyr::filter (mcc_12 >= 1) %>% 
  dplyr::select(hbp12_c,dia12_c,can12_c,lung12_c,hrt12_c,str12_c,art12_c)

sumtable(mcc_combo_1)

# identify most common combinations for those with >= 1 condtion
mcc_combo_2 <- mcc_combo_1  %>% 
  group_by(hbp12_c,dia12_c,can12_c,lung12_c,hrt12_c,str12_c,art12_c) %>% 
  dplyr::summarise(n = n()) %>%
  arrange(desc(n)) %>% 
  mutate (pct = n/9004)

View(mcc_combo_2)

# create categorical variables identifying physical multimorbidity pattern membership
df3 <- df2 %>% 
  mutate(combos = case_when (mcc_12 == 0 ~ 0,
                             mcc_12 == 1 ~ 1,
                             hbp12_c == 1 & dia12_c == 0 & can12_c == 0 & lung12_c == 0 & hrt12_c == 0 & str12_c == 0 & art12_c == 1 ~ 2,
                             hbp12_c == 1 & dia12_c == 0 & can12_c == 0 & lung12_c == 0 & hrt12_c == 1 & str12_c == 0 & art12_c == 1 ~ 3,
                             hbp12_c == 1 & dia12_c == 1 & can12_c == 0 & lung12_c == 0 & hrt12_c == 0 & str12_c == 0 & art12_c == 1 ~ 4,
                             hbp12_c == 1 & dia12_c == 1 & can12_c == 0 & lung12_c == 0 & hrt12_c == 1 & str12_c == 0 & art12_c == 1 ~ 5,
                             hbp12_c == 1 & dia12_c == 0 & can12_c == 1 & lung12_c == 0 & hrt12_c == 0 & str12_c == 0 & art12_c == 1 ~ 6,
                             hbp12_c == 1 & dia12_c == 1 & can12_c == 0 & lung12_c == 0 & hrt12_c == 0 & str12_c == 0 & art12_c == 0 ~ 7,
                             hbp12_c == 0 & dia12_c == 0 & can12_c == 0 & lung12_c == 0 & hrt12_c == 1 & str12_c == 0 & art12_c == 1 ~ 8,
                             hbp12_c == 0 & dia12_c == 0 & can12_c == 1 & lung12_c == 0 & hrt12_c == 0 & str12_c == 0 & art12_c == 1 ~ 9,
                             hbp12_c == 1 & dia12_c == 0 & can12_c == 0 & lung12_c == 1 & hrt12_c == 0 & str12_c == 0 & art12_c == 1 ~ 10,
                             TRUE ~ 11),
         no_cond_12 = case_when (mcc_12 == 0 ~ 1, TRUE ~ 0), 
         hbp_1_12 = case_when (mcc_12 == 1 & hbp12_c == 1 ~ 1, TRUE ~ 0), 
         diab_1_12 = case_when (mcc_12 == 1 & dia12_c == 1 ~ 1, TRUE ~ 0), 
         can_1_12 = case_when (mcc_12 == 1 & can12_c == 1 ~ 1, TRUE ~ 0), 
         lung_1_12 = case_when (mcc_12 == 1 & lung12_c == 1 ~ 1, TRUE ~ 0),
         hrt_1_12 = case_when (mcc_12 == 1 & hrt12_c == 1 ~ 1, TRUE ~ 0),
         strk_1_12 = case_when (mcc_12 == 1 & str12_c == 1 ~ 1, TRUE ~ 0),
         art_1_12 = case_when (mcc_12 == 1 & art12_c == 1 ~ 1, TRUE ~ 0), 
         hbp_art = case_when (combos == 2 ~ 1, TRUE ~ 0), 
         hbp_hrt_art = case_when (combos == 3 ~ 1, TRUE ~ 0),
         hbp_dia_art = case_when (combos == 4 ~ 1, TRUE ~ 0),
         hbp_dia_hrt_art = case_when (combos == 5 ~ 1, TRUE ~ 0), 
         hbp_can_art = case_when (combos == 6 ~ 1, TRUE ~ 0), 
         hbp_dia = case_when (combos == 7 ~ 1, TRUE ~ 0),
         hrt_art = case_when (combos == 8 ~ 1, TRUE ~ 0), 
         can_art = case_when (combos == 9 ~ 1, TRUE ~ 0), 
         hbp_lng_art = case_when (combos == 10 ~ 1, TRUE ~ 0),
         oth_comb = case_when (combos == 11 ~ 1, TRUE ~ 0),
         oth_sin_cond_12 = case_when(can_1_12 == 1 | hrt_1_12  == 1 | diab_1_12 == 1 | lung_1_12  == 1 | strk_1_12 == 1 ~ 1, TRUE ~ 0),
         
         combos_nm_reduc = case_when (no_cond_12 == 1 ~ "No conditions", 
                                      art_1_12 == 1 ~ "Arthritis only",  
                                      hbp_1_12 == 1 ~ "HBP only",  
                                      can_1_12 == 1 | hrt_1_12  == 1 | diab_1_12 == 1 | lung_1_12  == 1 | strk_1_12 == 1 ~ "Other single condition", 
                                      hbp_art == 1 ~ "ART+HBP",
                                      hbp_hrt_art == 1 ~ "ART+HBP+HRT",
                                      hbp_dia_art == 1 ~ "ART+HBP+DIAB",
                                      hbp_dia_hrt_art == 1 ~ "ART+HBP+DIA+HRT",
                                      hbp_can_art == 1 ~ "ART+HBP+CAN",
                                      hbp_dia == 1 ~ "HBP+DIA",
                                      hrt_art == 1 ~ "ART+HRT",
                                      can_art == 1 ~ "ART+CAN",
                                      hbp_lng_art == 1 ~ "ART+HBP+LNG",
                                      oth_comb == 1 ~ "Other condition combination"), 
         combos_nm_reduc = factor(combos_nm_reduc, levels=c("No conditions", 
                                                            "Arthritis only",  
                                                            "HBP only",  
                                                            "Other single condition", 
                                                            "ART+HBP",
                                                            "HBP+DIA",
                                                            "ART+HRT",
                                                            "ART+CAN",
                                                            "ART+HBP+HRT",
                                                            "ART+HBP+DIAB",
                                                            "ART+HBP+CAN",
                                                            "ART+HBP+LNG",
                                                            "ART+HBP+DIA+HRT",
                                                            "Other condition combination"))
  )

sumtable(df3, digits = 4)

# creating mplus file for analysis  ---------------------------------------
formplus_1 <- as_tibble(as.data.frame(df3))
prepareMplusData(formplus_1,"mplus_file_df3.dat")


# creating inverse probability weight -------------------------------------

# filter data by non-missing cesd, non-missing chronic conditions, age greater than 65, 
# excluding proxy responses and missing/zero weights

wide_ipw <- read_spss("randhrs1992_2020v2.sav") %>% 
  dplyr::filter(!is.na(R11CESD) & R11AGEY_B >= 65 & R11PROXY == 0 & R11WTRESP > 0 & 
                  R11HIBP %in% c(0,1) & R11DIAB %in% c(0,1) & R11CANCR %in% c(0,1) & R11LUNG %in% c(0,1) & 
                  R11HEART %in% c(0,1) & R11STROK %in% c(0,1) & R11ARTHR %in% c(0,1)) %>% 
  clean_names(., "snake") %>% 
  dplyr::select(
    hhidpn, hhid, pn,
    # time-stable demographics (wave 11 baseline)
    ragender, raracem,                        
    # respondent weight (wave 11)
    r11wtresp,
    # age
    r11agey_b, r12agey_b, r13agey_b, r14agey_b, r15agey_b,
    # marital status
    r11mstat, r12mstat, r13mstat, r14mstat, r15mstat,
    # household income
    h11itot, h12itot, h13itot, h14itot, h15itot,
    # household wealth
    h11atotb, h12atotb, h13atotb, h14atotb, h15atotb,
    # BMI
    r11bmi, r12bmi, r13bmi, r14bmi, r15bmi,
    # CES-D depression score
    r11cesd, r12cesd, r13cesd, r14cesd, r15cesd,
    # hypertension ever
    r11hibp, r12hibp, r13hibp, r14hibp, r15hibp,
    # diabetes ever
    r11diab, r12diab, r13diab, r14diab, r15diab,
    # cancer ever
    r11cancr, r12cancr, r13cancr, r14cancr, r15cancr,
    # lung disease ever
    r11lung, r12lung, r13lung, r14lung, r15lung,
    # heart disease ever
    r11heart, r12heart, r13heart, r14heart, r15heart,
    # stroke ever
    r11strok, r12strok, r13strok, r14strok, r15strok,
    # arthritis ever
    r11arthr, r12arthr, r13arthr, r14arthr, r15arthr,
    # drinks alcohol (0/1)
    r11drink, r12drink, r13drink, r14drink, r15drink,
    # drinking days per week
    r11drinkd, r12drinkd, r13drinkd, r14drinkd, r15drinkd,
    # drinks per drinking day
    r11drinkn, r12drinkn, r13drinkn, r14drinkn, r15drinkn,
    # vigorous physical activity
    r11vgactx, r12vgactx, r13vgactx, r14vgactx, r15vgactx,
    # total cognition score
    r11cogtot, r12cogtot, r13cogtot, r14cogtotp, r15cogtotp, 
    # ADL difficulties (0-5)
    r11adl5a, r12adl5a, r13adl5a, r14adl5a, r15adl5a,
    # smokes now
    r11smoken, r12smoken, r13smoken, r14smoken, r15smoken,
    # saw doctor
    r11doctor, r12doctor, r13doctor, r14doctor, r15doctor,
    # number of doctor visits
    r11doctim, r12doctim, r13doctim, r14doctim, r15doctim,
    # interview status (for dropout classification)
    r11iwstat, r12iwstat, r13iwstat, r14iwstat, r15iwstat
  ) %>%
  remove_val_labels()

# merge in already coded static variables
df3_r <- df3 %>% 
  dplyr::select(hhid, pn, female, mar_12, educ_c_1, educ_c_2, educ_c_3,
                black, other, latinx, nonus_b,
                ch_cir_1, ch_cir_2, ch_cir_3, ch_cir_4, ch_cir_5, ch_hlth_p)

# merge the datasets
ipw_df1 <- left_join(wide_ipw, df3_r, by=c("hhid","pn"))

sumtable(ipw_df1)

# create cross-wave summary variables
ipw_df2 <- ipw_df1 %>% 
  mutate(across(-c(hhidpn, female, mar_12, black, other, latinx, nonus_b,
                   ch_cir_1, ch_cir_2, ch_cir_3, ch_cir_4, ch_cir_5, ch_hlth_p, educ_c_1, educ_c_2,
                   educ_c_3), as.numeric)) %>%
  
  # --- continuous: row-median across waves ---
  mutate(
    mdn_age    = rowMedians(data.frame(r11agey_b, r12agey_b, r13agey_b, r14agey_b, r15agey_b), na.rm = TRUE),
    mdn_income = rowMedians(data.frame(h11itot,   h12itot,   h13itot,   h14itot,   h15itot),   na.rm = TRUE),
    mdn_wealth = rowMedians(data.frame(h11atotb,  h12atotb,  h13atotb,  h14atotb,  h15atotb),  na.rm = TRUE),
    mdn_bmi    = rowMedians(data.frame(r11bmi,    r12bmi,    r13bmi,    r14bmi,    r15bmi),     na.rm = TRUE),
    mdn_cesd   = rowMedians(data.frame(r11cesd,   r12cesd,   r13cesd,   r14cesd,   r15cesd),   na.rm = TRUE),
    mdn_adl    = rowMedians(data.frame(r11adl5a,  r12adl5a,  r13adl5a,  r14adl5a,  r15adl5a),  na.rm = TRUE),
    mdn_cogtot = rowMedians(data.frame(r11cogtot, r12cogtot, r13cogtot, r14cogtotp, r15cogtotp), na.rm = TRUE),
    mdn_drinkd = rowMedians(data.frame(r11drinkd, r12drinkd, r13drinkd, r14drinkd, r15drinkd), na.rm = TRUE),
    mdn_drinkn = rowMedians(data.frame(r11drinkn, r12drinkn, r13drinkn, r14drinkn, r15drinkn), na.rm = TRUE),
    mdn_doctim = rowMedians(data.frame(r11doctim, r12doctim, r13doctim, r14doctim, r15doctim), na.rm = TRUE),
    mdn_vigact = rowMedians(data.frame(r11vgactx, r12vgactx, r13vgactx, r14vgactx, r15vgactx), na.rm = TRUE),
    mdn_mstat = rowMedians(data.frame(r11mstat, r12mstat, r13mstat, r14mstat, r15mstat), na.rm = TRUE),
    
  ) %>%
  # --- binary ever-occurred: row-max across waves ---
  mutate(
    across(c(r11hibp, r12hibp, r13hibp, r14hibp, r15hibp,
             r11diab, r12diab, r13diab, r14diab, r15diab,
             r11cancr, r12cancr, r13cancr, r14cancr, r15cancr,
             r11lung,  r12lung,  r13lung,  r14lung,  r15lung,
             r11heart, r12heart, r13heart, r14heart, r15heart,
             r11strok, r12strok, r13strok, r14strok, r15strok,
             r11arthr, r12arthr, r13arthr, r14arthr, r15arthr),
           ~ case_when(. == 1 ~ 1, . == 0 ~ 0, TRUE ~ NA_real_)),
    
    ever_hibp  = pmax(r11hibp,  r12hibp,  r13hibp,  r14hibp,  r15hibp,  na.rm = TRUE),
    ever_diab  = pmax(r11diab,  r12diab,  r13diab,  r14diab,  r15diab,  na.rm = TRUE),
    ever_cancr = pmax(r11cancr, r12cancr, r13cancr, r14cancr, r15cancr, na.rm = TRUE),
    ever_lung  = pmax(r11lung,  r12lung,  r13lung,  r14lung,  r15lung,  na.rm = TRUE),
    ever_heart = pmax(r11heart, r12heart, r13heart, r14heart, r15heart, na.rm = TRUE),
    ever_strok = pmax(r11strok, r12strok, r13strok, r14strok, r15strok, na.rm = TRUE),
    ever_arthr = pmax(r11arthr, r12arthr, r13arthr, r14arthr, r15arthr, na.rm = TRUE),
    ever_smok  = pmax(r11smoken, r12smoken, r13smoken, r14smoken, r15smoken, na.rm = TRUE),
    ever_drink = pmax(r11drink,  r12drink,  r13drink,  r14drink,  r15drink,  na.rm = TRUE),
    ever_doc   = pmax(r11doctor, r12doctor, r13doctor, r14doctor, r15doctor, na.rm = TRUE),
    # clean up -Inf from pmax on all-NA rows
    ever_hibp  = case_when(is.infinite(ever_hibp)  ~ NA_real_, TRUE ~ ever_hibp),
    ever_diab  = case_when(is.infinite(ever_diab)  ~ NA_real_, TRUE ~ ever_diab),
    ever_cancr = case_when(is.infinite(ever_cancr) ~ NA_real_, TRUE ~ ever_cancr),
    ever_lung  = case_when(is.infinite(ever_lung)  ~ NA_real_, TRUE ~ ever_lung),
    ever_heart = case_when(is.infinite(ever_heart) ~ NA_real_, TRUE ~ ever_heart),
    ever_strok = case_when(is.infinite(ever_strok) ~ NA_real_, TRUE ~ ever_strok),
    ever_arthr = case_when(is.infinite(ever_arthr) ~ NA_real_, TRUE ~ ever_arthr),
    ever_smok  = case_when(is.infinite(ever_smok)  ~ NA_real_, TRUE ~ ever_smok),
    ever_drink = case_when(is.infinite(ever_drink) ~ NA_real_, TRUE ~ ever_drink),
    ever_doc   = case_when(is.infinite(ever_doc)   ~ NA_real_, TRUE ~ ever_doc)
  ) %>%
  # --- derived categorical variables ---
  mutate(
    # BMI categories (ref = normal weight)
    bmi_cat = as.factor(case_when(
      mdn_bmi < 18.5                  ~ 1,   # underweight
      mdn_bmi >= 18.5 & mdn_bmi < 25 ~ 2,   # normal (reference)
      mdn_bmi >= 25   & mdn_bmi < 30 ~ 3,   # overweight
      mdn_bmi >= 30                   ~ 4,   # obese
      TRUE ~ NA_real_
    )),
    bmi_cat = factor(bmi_cat, levels = c("2", "1", "3", "4")),
    
    vact_c = case_when(
      mdn_vigact %in% c(1,2) ~ 1,
      mdn_vigact %in% c(3,4,5) ~ 0,
      TRUE ~ NA_real_),
    
    mdn_mar_stat = case_when (mdn_mstat >= 1 & mdn_mstat <= 3 ~ 1,
                              mdn_mstat >= 4 & mdn_mstat <= 8 ~ 0,
                              TRUE ~ NA_real_),
    
    mdn_drnkwk=mdn_drinkn*mdn_drinkd,
    mdn_drnk = case_when(mdn_drnkwk == 0 ~ 0, 
                         female == 0 & (mdn_drnkwk >= 1 & mdn_drnkwk <= 14) ~ 1,
                         female == 1 & (mdn_drnkwk >= 1 & mdn_drnkwk <= 7) ~ 1,
                         female == 0 & (mdn_drnkwk > 14) ~ 2,
                         female == 1 & (mdn_drnkwk > 7) ~ 2, 
                         TRUE ~ NA_real_), 
    adl_any= case_when (mdn_adl == 0 ~ 0,
                        mdn_adl >= 1 ~ 1, 
                        TRUE ~ NA_real_)
  ) 

# Recode r##iwstat to dropout categories
# 1 = Resp alive (observed)
# 4 = NR, alive (sporadic/refusal - came back possible)
# 5 = NR, died this wave | these two are
# 6 = NR, died prev wave | permanent exits
# 7 = NR, dropped from sample (permanent)
# 0 = Inapplicable (not yet in sample)

ipw_df3 <- ipw_df2 %>%
  dplyr::mutate(
    # recode iwstat
    across(c(r11iwstat, r12iwstat, r13iwstat, r14iwstat, r15iwstat),
           ~ case_when(. == 1         ~ "resp",
                       . == 4         ~ "nr_alive",
                       . %in% c(5, 6) ~ "died",
                       . == 7         ~ "dropped",
                       . == 0         ~ "inap",
                       TRUE           ~ NA_character_)),
    
    # flag any death across waves
    has_death = (r11iwstat == "died" | r12iwstat == "died" | 
                   r13iwstat == "died" | r14iwstat == "died" | 
                   r15iwstat == "died"),
    has_death = case_when(is.na(has_death) ~ FALSE, TRUE ~ has_death),
    
    # flag any drop across waves
    has_dropped = (r11iwstat == "dropped" | r12iwstat == "dropped" | 
                     r13iwstat == "dropped" | r14iwstat == "dropped" | 
                     r15iwstat == "dropped"),
    has_dropped = case_when(is.na(has_dropped) ~ FALSE, TRUE ~ has_dropped),
    
    # last non missing wave status 
    last_obs = case_when(!is.na(r15iwstat) ~ r15iwstat,
                         !is.na(r14iwstat) ~ r14iwstat,
                         !is.na(r13iwstat) ~ r13iwstat,
                         !is.na(r12iwstat) ~ r12iwstat,
                         !is.na(r11iwstat) ~ r11iwstat,
                         TRUE              ~ NA_character_),
    
    # classify dropout type using flags above
    # trailing_nr for those whose last recorded wave was non-response
    dropout_type = case_when(has_death              ~ "permanent_death",
                             has_dropped            ~ "permanent_dropped",
                             last_obs == "nr_alive" ~ "trailing_nr",
                             TRUE                   ~ "retained"),
    
    # Binary flags derived from dropout_type
    permanent_exit  = dropout_type %in% c("permanent_death", "permanent_dropped"),
    trailing_nr     = dropout_type == "trailing_nr",
    
    # dropout flag for IPW
    continuous_dropout = permanent_exit | trailing_nr,  
    
    # for IPW, this is outcome
    no_dropout_d = as.factor(as.integer(!continuous_dropout))
  ) 

sumtable(ipw_df3, digits = 4)

# missing data imputation for dropout ipw covariates
# use doParallel to register cores, machine specific
registerDoParallel(cores = 15)

# select and type-cast variables
ipw_df3_imp <- ipw_df3 %>%
  dplyr::select(
    # continuous
    hhidpn, r11wtresp, mdn_age, 
    mdn_cesd, mdn_cogtot, 
    # demographics
    female, mdn_mar_stat,
    black, other, latinx, nonus_b,
    educ_c_1, educ_c_2, educ_c_3,
    ch_cir_1, ch_cir_2, ch_cir_3, ch_cir_4, ch_cir_5, ch_hlth_p,
    # continuous to be quartiles later
    mdn_income, mdn_wealth, mdn_doctim,
    # BMI dummies
    bmi_cat,
    # behaviors/health utilization
    ever_hibp, ever_diab, ever_cancr, ever_lung, ever_heart,
    ever_strok, ever_arthr, ever_smok, mdn_drnk, vact_c, adl_any,
    no_dropout_d
  ) %>%
  dplyr::mutate(
    across(c(hhidpn, r11wtresp, mdn_age, 
             mdn_cesd, mdn_cogtot,
             mdn_income, mdn_wealth, mdn_doctim), as.numeric),
    across(c(female, mdn_mar_stat,
             black, other, latinx, nonus_b,
             educ_c_1, educ_c_2, educ_c_3,
             ch_cir_1, ch_cir_2, ch_cir_3, ch_cir_4, ch_cir_5, ch_hlth_p,
             bmi_cat,
             ever_hibp, ever_diab, ever_cancr, ever_lung, ever_heart,
             ever_strok, ever_arthr, ever_smok, mdn_drnk, vact_c, adl_any,
             no_dropout_d),
           as.factor)
  )

vtable(ipw_df3_imp)
sumtable(ipw_df3_imp)

# Run missForest (drop ID before imputing)
vars_dropout_df <- ipw_df3_imp %>%
  dplyr::select(-hhidpn) %>%
  as.data.frame()

set.seed(123)
vars_dropout_imp <- missForest(
  vars_dropout_df,
  verbose = TRUE,
  ntree = 100,
  variablewise = TRUE,
  parallelize = "forests"   
)

# Check OOB imputation error
vars_dropout_imp$OOBerror

data.frame(
  variable = names(vars_dropout_imp$OOBerror),
  error_type = ifelse(names(vars_dropout_imp$OOBerror) == "MSE", "MSE", "PFC"),
  oob_error = vars_dropout_imp$OOBerror,
  col_name = names(vars_dropout_df)
) %>% arrange(desc(oob_error))

# extract and re-attach id
imputed_dropout <- vars_dropout_imp$ximp

# re-identify income wealth and doctor visit quartiles
imputed_dropout_final <- bind_cols(
  ipw_df3_imp %>% dplyr::select(hhidpn),
  imputed_dropout
) %>% 
  dplyr::mutate(
    inc_q  = as.factor(ntile(mdn_income, 4)),
    wlth_q = as.factor(ntile(mdn_wealth, 4)),
    doctim_q = as.factor(ntile(mdn_doctim, 4))
  )

vtable(imputed_dropout_final)

sumtable(imputed_dropout_final, digits = 4, group.weights = "r11wtresp")

save(imputed_dropout_final, file = "imputed_dropout_final.Rdata")

# Check VIF before weighting
vif_dropout <- glm(no_dropout_d ~   # your dropout indicator (0/1)
                     mdn_age + 
                     mdn_cesd + adl_any + mdn_cogtot +
                     inc_q + wlth_q + doctim_q +
                     female + mdn_mar_stat +
                     black + other + latinx + nonus_b +
                     educ_c_1 + educ_c_3 +
                     ch_cir_1 + ch_cir_2 + ch_cir_3 + ch_cir_4 + ch_cir_5 + ch_hlth_p +
                     bmi_cat +
                     ever_hibp + ever_diab + ever_cancr + ever_lung + ever_heart +
                     ever_strok + ever_arthr + ever_smok + mdn_drnk + vact_c,
                   data = imputed_dropout_final,       
                   family = "binomial")

car::vif(vif_dropout)

W_dropout <- weightit(
  no_dropout_d ~
    mdn_age + 
    mdn_cesd + adl_any + mdn_cogtot +
    inc_q + wlth_q + doctim_q +
    female + mdn_mar_stat +
    black + other + latinx + nonus_b +
    educ_c_1 + educ_c_3 +
    ch_cir_1 + ch_cir_2 + ch_cir_3 + ch_cir_4 + ch_cir_5 + ch_hlth_p +
    bmi_cat +
    ever_hibp + ever_diab + ever_cancr + ever_lung + ever_heart +
    ever_strok + ever_arthr + ever_smok + mdn_drnk + vact_c,
  data = imputed_dropout_final,
  method = "ps",
  estimand = "ATE",
  s.weights = 'r11wtresp')

summary(W_dropout)

# check balance
bal.tab(W_dropout, thresholds = c(m = .05))

# extract weights and bind back
wght_dropout <- tibble(ipdw_wghtit = W_dropout$weights)

vtable(imputed_dropout_final)

imputed_dropout_final_2 <- imputed_dropout_final %>% 
  dplyr::mutate(no_dropout_d_n = case_when (no_dropout_d == "0" ~ 0, 
                                            no_dropout_d == "1" ~ 1),
                r11wtresp_c = r11wtresp / mean(r11wtresp, na.rm = TRUE)
  )

vtable(imputed_dropout_final_2)

# calculate manually to cross-check (mirrors your IPW approach)
glm_dropout <- glm(
  no_dropout_d_n ~
    mdn_age + 
    mdn_cesd + adl_any + mdn_cogtot +
    inc_q + wlth_q + doctim_q +
    female + mdn_mar_stat +
    black + other + latinx + nonus_b +
    educ_c_1 + educ_c_3 +
    ch_cir_1 + ch_cir_2 + ch_cir_3 + ch_cir_4 + ch_cir_5 + ch_hlth_p +
    bmi_cat +
    ever_hibp + ever_diab + ever_cancr + ever_lung + ever_heart +
    ever_strok + ever_arthr + ever_smok + mdn_drnk + vact_c,
  data = imputed_dropout_final_2,
  family = "binomial", weights = r11wtresp_c
)

summary(glm_dropout)
tab_model(glm_dropout)
glm_probs_dropout <- data.frame(probs = predict(glm_dropout, type = "response"))

# calculate manual IPW
df_ipdw <- cbind(imputed_dropout_final_2, glm_probs_dropout, wght_dropout) %>%
  dplyr::mutate(
    ipdw = ((1 / probs) * no_dropout_d_n) + 
      ((1 / (1 - probs)) * (1 - no_dropout_d_n))
  )
# both approaches match

# trim at 98th percentile, calculate weighted survey weight
df_ipdw <- df_ipdw %>%
  dplyr::mutate(ipdw_pct = percent_rank(ipdw_wghtit),
                svy_ipw_wgt = ipdw_wghtit*r11wtresp)

nrow(df_ipdw)

df_ipdw_trim <- df_ipdw %>%
  dplyr::filter(ipdw_pct < .98)

nrow(df_ipdw_trim)

# combine with full data, filter out cases without IPW
sumtable(df_ipdw_trim)

df_ipdw_trim_2 <- df_ipdw_trim %>% 
  dplyr::select(hhidpn, svy_ipw_wgt)

df3_ipw <- left_join(df3, df_ipdw_trim_2, by = "hhidpn") %>% 
  dplyr::filter(!is.na(svy_ipw_wgt))

sumtable(df3_ipw)

# creating an mplus file for analysis 
formplus_2 <- as_tibble(as.data.frame(df3_ipw))

getwd()

prepareMplusData(formplus_2,"mplus_file_df3_ipw.dat")

# complete remaining analyses in Mplus




