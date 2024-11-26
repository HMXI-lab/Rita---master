library(tidyverse) 
library(dbplyr) 
library(readxl)

setwd("~/ITQB/Plate readings")

Exp3to25 <- read_excel("ODs_all_TPEN_exp_days_R.xlsx")
# View(Exp3to25)

Exp3to25_v1 <- Exp3to25 %>%
  pivot_longer(names_to = "Condition", values_to = "OD_578", cols = -Time_h)
# Exp3to25_v1
# View(Exp3to25_v1)

# Exp3to25_v1 <- Exp3to25_v1 %>% filter(Time_h < "25")

Exp3to25_v2 <- separate(Exp3to25_v1, Condition, 
                    c("Exp_number",
                      "Strain",
                      "Transfer_medium",
                      "Zn_condition", 
                      "Tech_Repl",
                      "Biol_Repl",
                      "TPEN_conc",
                      "Well_number"), 
                    sep = "_", remove = FALSE,
                    convert = FALSE)
# View(Exp3to25_v2)

Exp3to25_v2$OD_578 <- as.numeric (Exp3to25_v2$OD_578) #Because OD_578 values are character

# table(Exp3to25_v2$TPEN_conc)

Exp3to25_v2$TPEN_conc <- factor(Exp3to25_v2$TPEN_conc, levels = c("ctrTPENcells",
                                                                  "ctrTPENwZn",
                                                                  "ctrTPEN",
                                                                  "1.2",
                                                                  "2.4",
                                                                  "3.6",
                                                                  "4.8",
                                                                  "6",
                                                                  "7.2",
                                                                  "8.4",
                                                                  "9.6",
                                                                  "10.8",
                                                                  "12",
                                                                  "14",
                                                                  "16",
                                                                  "18",
                                                                  "20",
                                                                  "22",
                                                                  "24",
                                                                  "26",
                                                                  "28"))

# table(Exp3to25_v2$TPEN_conc)

Exp3to25_v2 <- Exp3to25_v2 %>% mutate(TPEN_conc = recode (TPEN_conc,
                                                  "ctrTPENcells" = "blank" ,
                                                  "ctrTPENwZn" = "17.4 uM Zn",
                                                  "ctrTPEN" = "0 uM",
                                                  "1.2" = "6 uM",
                                                  "2.4" = "12 uM",
                                                  "3.6"= "18 uM",
                                                  "4.8"= "24 uM",
                                                  "6" = "30 uM",
                                                  "7.2" = "36 uM",
                                                  "8.4" = "42 uM",
                                                  "9.6" = "48 uM",
                                                  "10.8" = "54 uM",
                                                  "12" = "60 uM",
                                                  "14" = "70 uM",
                                                  "16" = "80 uM",
                                                  "18" = "90 uM",
                                                  "20" = "100 uM",
                                                  "22" = "110 uM",
                                                  "24" = "120 uM",
                                                  "26" = "130 uM",
                                                  "28" = "140 uM"))

# table(Exp3to25_v2$TPEN_conc)


Exp3to25_v2$Strain <- factor(Exp3to25_v2$Strain, levels = c("5001",
                                                          "5002",
                                                          "5003",
                                                          "5004",
                                                          "5006",
                                                          "5007",
                                                          "5011",
                                                          "5019",
                                                          "5026",
                                                          "5032",
                                                          "5037",
                                                          "5038",
                                                          "5041",
                                                          "5046",
                                                          "5048",
                                                          "5071",
                                                          "5076",
                                                          "5078",
                                                          "12004",
                                                          "12010",
                                                          "12016",
                                                          "12097",
                                                          "12194",
                                                          "12195",
                                                          "14002",
                                                          "14071",
                                                          "14073",
                                                          "14076",
                                                          "14115",
                                                          "23006",
                                                          "HM02",
                                                          "HM04"))

# table(Exp3to25_v2$Strain)

Exp3to25_v2$Exp_number <- factor(Exp3to25_v2$Exp_number, levels = c("Exp3",
                                                            "Exp4",
                                                            "Exp5",
                                                            "Exp6",
                                                            "Exp7",
                                                            "Exp8a",
                                                            "Exp8b",
                                                            "Exp9",
                                                            "Exp10",
                                                            "Exp11",
                                                            "Exp12",
                                                            "Exp13",
                                                            "Exp14",
                                                            "Exp15",
                                                            "Exp16",
                                                            "Exp17",
                                                            "Exp18",
                                                            "Exp19",
                                                            "Exp20",
                                                            "Exp21",
                                                            "Exp22",
                                                            "Exp23",
                                                            "Exp24",
                                                            "Exp25"))

# table(Exp3to25_v2$Exp_number)

# Exp3to25_v2 %>%
#   ggplot(aes(x= Time_h, y=OD_578, color = Tech_Repl)) +
#   geom_line() +
#   facet_grid (Strain  ~ TPEN_conc)

Exp3to25_v3<- Exp3to25_v2 %>% group_by(Strain, Transfer_medium, Zn_condition, Tech_Repl, Biol_Repl, Well_number) %>%
  mutate(blank = Time_h ==0)

# View(Exp3to25_v3)
# table(Exp3to25_v3$TPEN_conc)
# table(Exp3to25_v3$Exp_number)

Exp3to25_v4 <- Exp3to25_v3 %>%
  group_by(Strain, Transfer_medium, Zn_condition, Tech_Repl, Biol_Repl, Well_number) %>%
  mutate(OD = OD_578 - OD_578[blank]) %>%
  ungroup()


# View(Exp3to25_v4)
# table(Exp3to25_v4$TPEN_conc)
# table(Exp3to25_v4$Exp_number)

# stocks conc / 0.001(mM --> uM) / d.f (=2)--> conc TPEN (now in the table)

Exp3to25_v4_M8 <- Exp3to25_v4 %>% filter(Transfer_medium == "M8")
# table(Exp3to25_v4_M8$Exp_number)
Exp3to25_v4_mGAM <- Exp3to25_v4 %>% filter(Transfer_medium == "mGAM")
# table(Exp3to25_v4_mGAM$Exp_number)

SAVE_PLOTS <- T

plotStrain <- function(STRAIN1) {
  p <- Exp3to25_v4_M8 %>% 
    filter(Strain == STRAIN1, Time_h < 30) %>% 
    ggplot(mapping=aes(x = Time_h,
                       y = OD_578,
                       color = Transfer_medium)) +
    geom_point() + 
    geom_line() +
    facet_grid (TPEN_conc ~ Exp_number + Biol_Repl + Tech_Repl) +
    theme_bw() +
    theme(legend.position = "bottom", 
          strip.background = element_rect(fill = "lightgrey", colour = "lightgrey"),
          panel.grid.minor = element_blank(),
          axis.text = element_text(size = 5),  # Decrease the size of the numbers on the axes
          strip.text.x = element_text(size = 5),  # Decrease facet label size
          strip.text.y = element_text(size = 5)
    ) + 
    ggtitle(STRAIN1)
  print(p)
}


if (SAVE_PLOTS) {
  pdf("AllDataM8_3_25.pdf", width = 8, height = 8)
  
  Exp3to25_v4_M8 %>% filter(Time_h < 30) %>% pull(Strain) %>% unique() %>%
    lapply(plotStrain) %>% invisible()
  dev.off()
}

plotStrain <- function(STRAIN1) {
  p <- Exp3to25_v4_mGAM %>% 
    filter(Strain == STRAIN1, Time_h <30) %>% 
    ggplot(mapping=aes(x = Time_h,
                       y = OD_578,
                       color = Transfer_medium)) +
    geom_point() + 
    geom_line() +
    facet_grid (TPEN_conc ~ Exp_number + Biol_Repl + Tech_Repl) +
    theme_bw() +
    theme(legend.position = "bottom", 
          strip.background = element_rect(fill = "lightgrey", colour = "lightgrey"),
          panel.grid.minor = element_blank(),
          axis.text = element_text(size = 5),  # Decrease the size of the numbers on the axes
          strip.text.x = element_text(size = 5),  # Decrease facet label size
          strip.text.y = element_text(size = 5)
    ) + 
    ggtitle(STRAIN1)
  print(p)
}


if (SAVE_PLOTS) {
  pdf("AllDataMGAM_3_25.pdf", width = 8, height = 8)
  
  Exp3to25_v4_mGAM %>% filter(Time_h < 30) %>% pull(Strain) %>% unique() %>%
    lapply(plotStrain) %>% invisible()
  dev.off()
}

######Clean the files

#Model for Cleaning the files


#-------------mGAM--------------------------
All_R_mGAM_clean <- Exp3to25_v4_mGAM  %>%  filter (
    !(Exp_number %in% c("Exp3","Exp4", "Exp5", "Exp6", "Exp7", "Exp9")), #Exp3:
    #Exp4: this Exp day had problems, because all strains used had growth curves which should be discarded thanks to inconsistency related to other Exps
    #Exp5:
    #Exp6:
    #Exp7: this Exp day had problems, because all strains used had growth curves which should be discarded thanks to inconsistency related to other Exps
    #Exp9: this Exp day had problems, because all strains used had growth curves which should be discarded thanks to inconsistency related to other Exps
    !(Exp_number == "Exp11" &  Tech_Repl == "TechRep2" & TPEN_conc == "100 uM"), #comparing to other BRs, this conc results were unexpected (probable cause: pipetting error)
    !(Exp_number == "Exp23" &  Tech_Repl == "TechRep3" & TPEN_conc == "blank"), #contamination w/ tips
    !(Strain %in% c("5041", "14115")), #variable growth curves #doesn't grow
    !(Strain == "5001" &  Exp_number == "Exp19" & Tech_Repl == "TechRep1" & Time_h > 20), #because growth curves presented strange results after 20h --> possible cause: solution went to the membr. when we put more plates in the plate reader
    !(Strain == "5001" & Exp_number == "Exp10" & Tech_Repl == "TechRep3" &  TPEN_conc == "90 uM"), #comparing to other BRs, this was unexpected (probable cause: pipetting error)
    !(Strain == "5003" & Exp_number == "Exp10" & Tech_Repl == "TechRep3" &  TPEN_conc == "90 uM"), #comparing to other BRs, this was unexpected (probable cause: pipetting error)
    !(Strain == "5006" &  Exp_number == "Exp10"), #this Exp day had problems, because all strains used had growth curves which should be discarded thanks to inconsistency related to other Exps
    !(Strain == "5006" &  Exp_number == "Exp18"), #strange curves + other better BRs
    !(Strain == "5006" &  Exp_number == "Exp21" & Tech_Repl == "TechRep3"),
    !(Strain == "5006" & Exp_number == "Exp19" & Tech_Repl == "TechRep2" &  TPEN_conc == "blank"), #I added strains to the control
    !(Strain == "5007" & Exp_number == "Exp19" & Tech_Repl == "TechRep2" &  TPEN_conc == "blank"), #I added strains to the control
    !(Strain == "5007" & Exp_number == "Exp16" & Tech_Repl == "TechRep3" &  TPEN_conc == "70 uM"), #
    !(Strain == "5007" & Exp_number == "Exp17" & Tech_Repl == "TechRep2" &  TPEN_conc == "80 uM"), #
    !(Strain == "5007" & Exp_number == "Exp19"), #
    !(Strain == "5011" &  Exp_number == "Exp19"), #
    !(Strain == "5019" & Exp_number == "Exp11" & Tech_Repl == "TechRep3"), #comparing to other BRs, this was unexpected (probable cause: pipetting error)
    !(Strain == "5026" &  Exp_number == "Exp12"), #
    !(Strain == "5026" &  Exp_number == "Exp11" & Tech_Repl == "TechRep2" & TPEN_conc == "100 µM"), # The only one TechRep which presented growth. Other BioRep never presented growth at this concentration.
    !(Strain == "5026" & Exp_number == "Exp19" & Tech_Repl == "TechRep1" &  TPEN_conc == "130 uM"), #
    !(Strain == "5026" & Exp_number == "Exp11" &  TPEN_conc == "90 uM"), #comparing to other BRs, these conc results was unexpected (probable cause: pipetting error)
    !(Strain == "5026" & Exp_number == "Exp11" &  TPEN_conc == "100 uM"), #comparing to other BRs, these conc results was unexpected (probable cause: pipetting error)
    !(Strain == "5032" &  Exp_number == "Exp17" & Tech_Repl == "TechRep3" ), #
    !(Strain == "5032" & Exp_number == "Exp24" &  Tech_Repl == "TechRep1" & TPEN_conc == "80 uM"), #something happened at this conc (affected all strains)(probably cause: pipetting errors)
    !(Strain == "5038" &  Exp_number == "Exp15"), #lost of 2 first reads + 3 better BRs
    !(Strain == "5038" &  Exp_number == "Exp17" & Tech_Repl == "TechRep3" & TPEN_conc == "80 µM"), #There was a contamination w/ 14002 of TechRep3
    !(Strain == "5038" & Exp_number == "Exp15" & TPEN_conc == "blank"),
    !(Strain == "5038" & Exp_number == "Exp19" & Tech_Repl == "TechRep1" &  TPEN_conc == "130 uM"), #  
    !(Strain == "5046" & Exp_number == "Exp14" & Tech_Repl == "TechRep1" &  TPEN_conc == "120 uM"), #
    !(Strain == "5046" & Exp_number == "Exp14" & Tech_Repl == "TechRep1" &  TPEN_conc == "130 uM"), #
    !(Strain == "5048" & Exp_number == "Exp14" & Tech_Repl == "TechRep1" &  TPEN_conc == "120 uM"), #
    !(Strain == "5048" & Exp_number == "Exp14" & Tech_Repl == "TechRep1" &  TPEN_conc == "130 uM"), #
    !(Strain == "5071" & Exp_number == "Exp19" & Tech_Repl == "TechRep1" &  TPEN_conc == "130 uM"), #
    !(Strain == "5071" & Exp_number == "Exp19" & Tech_Repl == "TechRep2" &  TPEN_conc == "130 uM"), #
    !(Strain == "5071" &  Exp_number == "Exp18"), #possible 5078 contamination
    !(Strain == "5076" &  Exp_number == "Exp18"), #
    !(Strain == "12016" &  Exp_number == "Exp14"), #Aparentely I already deleted all conc curves except blank--> delete Exp14 (completely)
    !(Strain == "12194" & Exp_number == "Exp19" & Tech_Repl == "TechRep1" &  TPEN_conc == "130 uM"), #
    !(Strain == "14071" &  Exp_number == "Exp17"), #
    !(Strain == "14071" &  Exp_number == "Exp25"), #thanks mistake in medium quantity calculations --> bubbles --> pipetting errors
    !(Strain == "14071" & Exp_number == "Exp16" & Tech_Repl == "TechRep2" &  TPEN_conc == "0 uM"), #
    !(Strain == "14073" & Exp_number == "Exp17" & Tech_Repl == "TechRep2" &  TPEN_conc == "90 uM"), # The only one TechRep which presented growth. Other BioRep never presented growth at this concentration.
    !(Strain == "HM04" &  Exp_number == "Exp11")) #
  
 
SAVE_PLOTS <- T

plotStrain <- function(STRAIN1) {
  p <- All_R_mGAM_clean %>% 
    filter(Strain == STRAIN1, Time_h < 22) %>% 
    ggplot(mapping=aes(x = Time_h,
                       y = OD_578,
                       color = Transfer_medium)) +
    geom_point() + 
    geom_line() +
    facet_grid (TPEN_conc ~ Exp_number + Biol_Repl + Tech_Repl) +
    theme_bw() +
    theme(legend.position = "bottom", 
          strip.background = element_rect(fill = "lightgrey", colour = "lightgrey"),
          panel.grid.minor = element_blank()
    ) + 
    ggtitle(STRAIN1)
  print(p)
}


if (SAVE_PLOTS) {
  pdf("AllDataMGAM_3_25_clean_Rita.pdf", width = 8, height = 8)
  
  All_R_mGAM_clean %>% filter(Time_h < 30) %>% pull(Strain) %>% unique() %>%
    lapply(plotStrain) %>% invisible()
  dev.off()
}

#Model for calculating the OD and sd. Here, you have to group by with your variables!!!
#Calculate the median OD per Strain, Transfer_medium, Zn_condition, TPEN_conc
All_R_mGAM_clean_summary <- All_R_mGAM_clean  %>% group_by(Time_h, Strain,Transfer_medium, Zn_condition,  TPEN_conc) %>%
  summarise(ODmean = mean(OD, na.rm = TRUE),ODsd = sd (OD, na.rm = TRUE)) 

All_R_mGAM_clean_summary %>% 
  filter(Time_h < 24, TPEN_conc != "blank") %>% 
  ggplot(mapping=aes(x = Time_h,
                     y = ODmean)) +
  #geom_point() + 
  geom_line()+
  geom_ribbon(aes(x=Time_h, ymax = ODmean + ODsd, ymin = ODmean - ODsd), 
              alpha = 0.5, colour = NA)+
  facet_grid (Strain ~ TPEN_conc) +
  theme_bw() +
  theme(legend.position = "bottom", 
        strip.background = element_rect(fill = "lightgrey", colour = "lightgrey"),
        panel.grid.minor = element_blank())

SAVE_PLOTS <- T

plotStrain <- function(STRAIN1) {
  p <- All_R_mGAM_clean_summary %>% 
    filter(Strain == STRAIN1, Time_h < 22) %>% 
    ggplot(mapping=aes(x = Time_h,
                       y = ODmean,
                       color = TPEN_conc)) +
    geom_point() + 
    geom_line() +
    geom_ribbon(aes(x=Time_h, ymax = ODmean + ODsd, ymin = ODmean - ODsd, fill = TPEN_conc),
                alpha = 0.5, colour = NA) + 
    facet_wrap (TPEN_conc ~ .) +
    theme_bw() +
    theme(legend.position = "bottom", 
          strip.background = element_rect(fill = "lightgrey", colour = "lightgrey"),
          panel.grid.minor = element_blank()
    ) + 
    ggtitle(STRAIN1)
  print(p)
}


if (SAVE_PLOTS) {
  pdf("AllDataMGAM_3_25_clean_Rita_summary.pdf", width = 8, height = 8)
  
  All_R_mGAM_clean_summary %>% filter(Time_h < 30) %>% pull(Strain) %>% unique() %>%
    lapply(plotStrain) %>% invisible()
  dev.off()
}

#-------------M8------------------------------
All_R_M8_clean <- Exp3to25_v4_M8  %>%  filter (
  !(Exp_number %in% c("Exp3", "Exp4", "Exp5", "Exp8a", "Exp8b")),  #Exp3: 1st exps --> possible contaminations + many pipetting errors
  #Exp4: 1st exps --> solution in the membrane + pipetting errors + possible contamination (control)
  #Exp5: 1st exps --> TPEN concentrations (2nd TPEN stocks preparation) wasn't ok + pipetting errors
  #Exp8a: 2x strains (some wells)
  #Exp8b: unexpected growth curves of 5037 and 5038 + cross contamination (HM04 - higher concentrations - With 5078?)
  !(Exp_number == "Exp24" & Tech_Repl == "TechRep3"), #solution in the membr
  !(Strain %in% c("5019", "5071", "5011", "14115", "5041", "5007")), #don't grow (5011: sometimes don't grow)
  !(Strain == "5001"), #slow growth in M8 --> new approach (transfer w/o Zn vs Zn & ZnCl2 vs TPEN chelator)
  !(Strain == "5002"), #slow growth in M8 --> new approach (transfer w/o Zn vs Zn & ZnCl2 vs TPEN chelator)
  !(Strain == "5003"), #slow growth in M8 --> new approach (transfer w/o Zn vs Zn & ZnCl2 vs TPEN chelator)
  !(Strain == "5004"), #slow growth in M8 --> new approach (transfer w/o Zn vs Zn & ZnCl2 vs TPEN chelator)
  !(Strain == "5026" &  Exp_number == "Exp15"), # strange growth curves (from another strain?)
  !(Strain == "5026" &  Exp_number == "Exp18" & TPEN_conc == "54 uM"), #grew only in this BR --> delete
  !(Strain == "5032" &  Exp_number == "Exp17"), #Exp17 growth curves seems different. However BioReps aren't very consistent.
  !(Strain == "5032" &  Exp_number == "Exp20" & TPEN_conc == "36 uM"), #solution in the membr
  !(Strain == "5032" &  Exp_number == "Exp20" & Tech_Repl == "TechRep2" & TPEN_conc == "42 uM"), #solution in the membr
  !(Strain == "5038" &  Exp_number == "Exp24" & Tech_Repl == "TechRep1" & TPEN_conc == "17.4 uM Zn"), #only TR when didn't grow
  !(Strain == "5076" &  Exp_number == "Exp18"), # remain: 2 TechReps
  # !(Strain == "12004" &  Exp_number == "Exp10"),
  !(Strain == "12016" &  Exp_number == "Exp13"), #2x strains (some wells) and lost of reads (thanks to stacker challenges)
  !(Strain == "12097" &  Exp_number == "Exp13"), #2x strains (some wells) and lost of reads (thanks to stacker challenges)
  !(Strain == "12194" &  Exp_number == "Exp17"), #unexpected growth curve (different from other BioReps)
  !(Strain == "HM04" &  Exp_number == "Exp12")) #Exp12 growth curves seems different.


plotStrain <- function(STRAIN1) {
  p <- All_R_M8_clean %>% 
    filter(Strain == STRAIN1, Time_h < 30) %>% 
    ggplot(mapping=aes(x = Time_h,
                       y = OD_578,
                       color = Transfer_medium)) +
    geom_point() + 
    geom_line() +
    facet_grid (TPEN_conc ~ Exp_number + Biol_Repl + Tech_Repl) +
    theme_bw() +
    theme(legend.position = "bottom", 
          strip.background = element_rect(fill = "lightgrey", colour = "lightgrey"),
          panel.grid.minor = element_blank()
    ) + 
    ggtitle(STRAIN1)
  print(p)
}


if (SAVE_PLOTS) {
  pdf("AllDataM8_3_25_clean_Rita.pdf", width = 8, height = 8)
  
  All_R_M8_clean %>% filter(Time_h < 30) %>% pull(Strain) %>% unique() %>%
    lapply(plotStrain) %>% invisible()
  dev.off()
}


 #Model for calculating the OD and sd. Here, you have to group by with your variables!!!
#Calculate the median OD per Strain, Transfer_medium, Zn_condition, TPEN_conc

All_R_M8_clean_summary <- All_R_M8_clean  %>% group_by(Time_h,Strain, Transfer_medium, Zn_condition, TPEN_conc) %>%
  summarise(ODmean = mean(OD, na.rm = TRUE),ODsd = sd (OD, na.rm = TRUE)) 


All_R_M8_clean_summary %>% 
  filter(Time_h < 24, TPEN_conc != "blank") %>% 
  ggplot(mapping=aes(x = Time_h,
                     y = ODmean)) +
  #geom_point() + 
  geom_line()+
  geom_ribbon(aes(x=Time_h, ymax = ODmean + ODsd, ymin = ODmean - ODsd), 
              alpha = 0.5, colour = NA)+
  facet_grid (Strain ~ TPEN_conc) +
  theme_bw() +
  theme(legend.position = "bottom", 
        strip.background = element_rect(fill = "lightgrey", colour = "lightgrey"),
        panel.grid.minor = element_blank())

SAVE_PLOTS <- T

plotStrain <- function(STRAIN1) {
  p <- All_R_M8_clean_summary %>% 
    filter(Strain == STRAIN1, Time_h < 30) %>% 
    ggplot(mapping=aes(x = Time_h,
                       y = ODmean,
                       color = TPEN_conc)) +
    geom_point() + 
    geom_line() +
    geom_ribbon(aes(x=Time_h, ymax = ODmean + ODsd, ymin = ODmean - ODsd, fill = TPEN_conc),
                alpha = 0.5, color = NA) + 
    facet_wrap (TPEN_conc ~ .) +
    theme_bw() +
    theme(legend.position = "bottom", 
          strip.background = element_rect(fill = "lightgrey", colour = "lightgrey"),
          panel.grid.minor = element_blank()
    ) + 
    ggtitle(STRAIN1)
  print(p)
}


if (SAVE_PLOTS) {
  pdf("AllDataM8_3_25_clean_Rita_summary.pdf", width = 8, height = 8)
  
  All_R_M8_clean_summary %>% filter(Time_h < 30) %>% pull(Strain) %>% unique() %>%
    lapply(plotStrain) %>% invisible()
  dev.off()
}

######------ M8 ------------------------

#1. Calculate the AUC
#2. Do the representation with the AUC
#3. Generate heatmap



gc_M8 <- All_R_M8_clean %>% filter (Time_h <22)


colnames(gc_M8) <- c("time", "sample_id", "experiment_number", "strain", "overnight", "Zn_condition", "tech_repl", "biol_repl", "conc",  "well", "OD", "blank", "ODblank")

gc_M8 <- gc_M8 %>% mutate(control = conc == "0 uM", ODc0 = ODblank)
#print(class(gc_M8$time))
# print(unique(gc_M8$time))

cumulmin <- function(v) {
  N <- length(v)
  if (N==0) return()
  sapply(1:N, function(i)
    min(v[i:N])
  )
}

determineTimeCutoffs <- function(GCs, decrease = 0.02) {
  GCs_controls <- GCs %>% filter(control) %>% mutate(time = round(time, 0)) %>%
    group_by(strain, experiment_number, time) %>%
    summarise(ODc0 = median(ODc0)) %>% arrange(time) %>%
    group_by(strain, experiment_number) %>% mutate(ODc0_adj = ODc0 + ODc0[n()]*decrease*(time[n()]-time))
  
  GCs_controls %>% group_by(strain, experiment_number) %>% filter(ODc0_adj == max(ODc0_adj)) %>%
    ungroup() %>% select(strain, experiment_number, time) %>% rename(cutoff_time = time)
}

gc_M8 <- gc_M8 %>% mutate(f_experiment_number = as.factor(experiment_number))

time_cutoffs <- determineTimeCutoffs(gc_M8) %>%
  mutate(f_experiment_number = as.factor(experiment_number))

 print(time_cutoffs)

gc_M8_control <- gc_M8 %>% 
  filter(control) %>% 
  group_by(f_experiment_number, experiment_number, strain, time) %>%
  summarise(OD = median(OD))


markArtefacts <- function(GCs) {
  
  GCsf_M8 <- GCs %>% group_by(experiment_number, biol_repl) %>% 
    mutate(ODc01 = ODc0 / median(ODc0[ control & time == max(time) ])) %>% arrange(time)
  
  GCsf_M8 <- GCsf_M8 %>% mutate(
    min_cutoff_time = min(time),
    is_start = time == min_cutoff_time & !is.na(ODc01) & !is.infinite(ODc01) & ODc01 < 1) %>%
    mutate(mean_start_OD = mean(ODc01[is_start]), sd_start_OD = sd(ODc01[is_start])) %>% 
    group_by(sample_id) %>%
    mutate(
      # only check for increased growth in first half, decrease later is ok
      increased = time < 4 & ODc01 > 4*sd_start_OD+cumulmin(ODc01),
      discard_conc = (sum(increased) > 3) | !(min(ODc01[time < min_cutoff_time + 1.5]) < 4 * sd_start_OD) )
  
  GCsf_M8 %>% ungroup() %>% select(-mean_start_OD, -sd_start_OD, -min_cutoff_time)
}

# gc_M8 %>% filter(experiment_number=="Exp107", plate_number == "plate4") %>% count(strain)

gc_M8 <- gc_M8 %>% left_join(time_cutoffs) %>% filter(time <= cutoff_time) %>% 
  select(-cutoff_time)

GCsf_M8 <- gc_M8 %>% group_by(strain) %>% do(markArtefacts(.))  



#GCsf <- GCsf %>% filter(!increased, !discard_conc)

calcAUC <- function(OD, time) {
  OD <- OD - min(OD)
  N <- length(OD)
  if (N<2) return(NA)
  s1 <- 1:(N-1)
  s2 <- s1 + 1
  sum( (OD[s1] + OD[s2])/2 * (time[s2]-time[s1]) )
}

species <- read_excel("Strains.xlsx")

species$strain <- as.character(species$strain)

GCsf_M8$strain <- as.character(GCsf_M8$strain)

GCsf_M8 <- GCsf_M8 %>% left_join(species, by = "strain")


AUCs_M8 <- GCsf_M8 %>% group_by(Phylum, Class, Order, Family, Genus,strain, name, short, experiment_number,biol_repl, tech_repl,
                                overnight, Zn_condition, conc, control, well) %>% 
  summarise(AUC = calcAUC(ODc01, time), maxOD = max(ODc01), finalOD = ODc01[n()])


AUCs_M8 <-  AUCs_M8 %>% group_by(Phylum, Class, Order, Family, Genus, strain, name, short, experiment_number, biol_repl, tech_repl) %>% 
  mutate(AUC = AUC / median(AUC[control]))


AUCs_M8_mean <-  AUCs_M8 %>% group_by(Phylum, Class, Order, Family, Genus, strain, name, short, conc) %>% 
  summarise (AUCsmean = mean(AUC), sdAUC = sd (AUC), maxODmean = mean (maxOD), sdmax = sd (maxOD), finalODmean = mean(finalOD), sdfinal = sd(finalOD)) %>% ungroup()



graphics.off()

AUCs_M8_mean %>% 
  filter (conc != "blank") %>%
  ggplot(aes(conc, short, fill = AUCsmean)) +
  geom_raster() +
  #geom_text(data = pvs %>% filter(qv_rescue < 0.01), label = "*") +
  scale_fill_gradient2(low = "black", high = "blueviolet") +
  #geom_text_repel(aes(label = Well_n)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 6.5),
    axis.text.y = element_text(size = 10)) +
  ylab("strains") +
  xlab ("TPEN_conc")
  


#-------------------To group M8 results by class in heatmap-------------------
AUCs_M8_mean_2  <- AUCs_M8_mean  %>% group_by(Class,  conc) %>%
  summarise(AUCsmean = mean(AUCsmean, na.rm = TRUE)) 


AUCs_M8_mean_2$Class<- factor(AUCs_M8_mean_2$Class, levels = c("Gammaproteobacteria", "Bacilli", "Erysipelotrichia", "Clostridia" ))


AUCs_M8_mean_2 %>% 
  filter (conc != "blank") %>%
  ggplot(aes(conc, Class, fill = AUCsmean)) +
  geom_raster() +
  #geom_text(data = pvs %>% filter(qv_rescue < 0.01), label = "*") +
  scale_fill_gradient2(low = "black", high = "blueviolet") +
  #geom_text_repel(aes(label = Well_n)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 12)) +
  xlab ("TPEN_conc")



AUCs_M8_mean %>% write_csv2("AUCs_M8mean.tsv")

AUCs_M8_select <-  AUCs_M8_mean %>% filter (conc != "blank") %>% select(short, conc, AUCsmean) %>% ungroup ()

AUCs_M8_select_spread <-  AUCs_M8_select  %>% pivot_wider (names_from = short, values_from = AUCsmean) 
AUCs_M8_select_spread

AUCs_M8_select_spread_m <- as.matrix( AUCs_M8_select_spread[, -1], rownames.force = NA)
AUCs_M8_select_spread_m

rownames( AUCs_M8_select_spread_m) <-  AUCs_M8_select_spread$conc
AUCs_M8_select_spread_m

heatmap( AUCs_M8_select_spread_m, scale = "none")


heatmap.2( AUCs_M8_select_spread_m, scale = "none", col = bluered(100), 
           trace = "none", density.info = "none",
           margins = c(8,11),
           #ColSideColors=col_colors,
           key = TRUE,
           cexCol = 1.2,
           cexRow = 1.2)


heatmap.2( AUCs_M8_select_spread_m, 
           scale = "none", 
           col = colorRampPalette(c("ivory", "darkgoldenrod1", "red"))(100),
           trace = "none", 
           margins = c(10, 11),
           key = TRUE,
           cexCol = 1.2,
           cexRow = 1.2,
           breaks = seq(0, 2, length.out = 101))


######------ mGAM ------------------------

#1. Calculate the AUC
#2. Do the representation with the AUC
#3. Generate heatmap



gc_mGAM <- All_R_mGAM_clean %>% filter (Time_h <22)


colnames(gc_mGAM) <- c("time", "sample_id", "experiment_number", "strain", "overnight", "Zn_condition", "tech_repl", "biol_repl", "conc",  "well", "OD", "blank", "ODblank")

gc_mGAM <- gc_mGAM %>% mutate(control = conc == "0 uM", ODc0 = ODblank)

cumulmin <- function(v) {
  N <- length(v)
  if (N==0) return()
  sapply(1:N, function(i)
    min(v[i:N])
  )
}

determineTimeCutoffs <- function(GCs, decrease = 0.02) {
  GCs_controls <- GCs %>% filter(control) %>% mutate(time = round(time, 0)) %>%
    group_by(strain, experiment_number, time) %>%
    summarise(ODc0 = median(ODc0)) %>% arrange(time) %>%
    group_by(strain, experiment_number) %>% mutate(ODc0_adj = ODc0 + ODc0[n()]*decrease*(time[n()]-time))
  
  GCs_controls %>% group_by(strain, experiment_number) %>% filter(ODc0_adj == max(ODc0_adj)) %>%
    ungroup() %>% select(strain, experiment_number, time) %>% rename(cutoff_time = time)
}

gc_mGAM <- gc_mGAM %>% mutate(f_experiment_number = as.factor(experiment_number))

time_cutoffs <- determineTimeCutoffs(gc_mGAM) %>% mutate(f_experiment_number = as.factor(experiment_number))

gc_mGAM_control <- gc_mGAM %>% filter(control) %>% group_by(f_experiment_number, experiment_number, strain, time) %>%
  summarise(OD = median(OD))


markArtefacts <- function(GCs) {
  
  GCsf_mGAM <- GCs %>% group_by(experiment_number, biol_repl) %>% 
    mutate(ODc01 = ODc0 / median(ODc0[ control & time == max(time) ])) %>% arrange(time)
  
  GCsf_mGAM <- GCsf_mGAM %>% mutate(
    min_cutoff_time = min(time),
    is_start = time == min_cutoff_time & !is.na(ODc01) & !is.infinite(ODc01) & ODc01 < 1) %>%
    mutate(mean_start_OD = mean(ODc01[is_start]), sd_start_OD = sd(ODc01[is_start])) %>% 
    group_by(sample_id) %>%
    mutate(
      # only check for increased growth in first half, decrease later is ok
      increased = time < 4 & ODc01 > 4*sd_start_OD+cumulmin(ODc01),
      discard_conc = (sum(increased) > 3) | !(min(ODc01[time < min_cutoff_time + 1.5]) < 4 * sd_start_OD) )
  
  GCsf_mGAM %>% ungroup() %>% select(-mean_start_OD, -sd_start_OD, -min_cutoff_time)
}

# gc_mGAM %>% filter(experiment_number=="Exp107", plate_number == "plate4") %>% count(strain)

gc_mGAM <- gc_mGAM %>% left_join(time_cutoffs) %>% filter(time <= cutoff_time) %>% 
  select(-cutoff_time)

GCsf_mGAM <- gc_mGAM %>% group_by(strain) %>% do(markArtefacts(.))  



#GCsf <- GCsf %>% filter(!increased, !discard_conc)

calcAUC <- function(OD, time) {
  OD <- OD - min(OD)
  N <- length(OD)
  if (N<2) return(NA)
  s1 <- 1:(N-1)
  s2 <- s1 + 1
  sum( (OD[s1] + OD[s2])/2 * (time[s2]-time[s1]) )
}

species <- read_excel("Strains.xlsx")

species$strain <- as.character(species$strain)

GCsf_mGAM$strain <- as.character(GCsf_mGAM$strain)

GCsf_mGAM <- GCsf_mGAM %>% left_join(species, by = "strain")


AUCs_mGAM <- GCsf_mGAM %>% group_by(Phylum, Class, Order, Family, Genus, strain, name, short, experiment_number,biol_repl, tech_repl,
                                    overnight, Zn_condition, conc, control, well) %>% 
  summarise(AUC = calcAUC(ODc01, time), maxOD = max(ODc01), finalOD = ODc01[n()])


AUCs_mGAM <-  AUCs_mGAM %>% group_by(Phylum, Class, Order, Family, Genus, strain, name, short, experiment_number, biol_repl, tech_repl) %>% 
  mutate(AUC = AUC / median(AUC[control]))


AUCs_mGAM_mean <-  AUCs_mGAM %>% group_by(Phylum, Class, Order, Family, Genus, strain, name, short, conc) %>% 
  summarise (AUCsmean = mean(AUC), sdAUC = sd (AUC), maxODmean = mean (maxOD), sdmax = sd (maxOD), finalODmean = mean(finalOD), sdfinal = sd(finalOD)) %>% ungroup()



graphics.off()

AUCs_mGAM_mean %>% 
  filter (conc != "blank") %>%
  ggplot(aes(conc, short, fill = AUCsmean)) +
  geom_raster() +
  #geom_text(data = pvs %>% filter(qv_rescue < 0.01), label = "*") +
  scale_fill_gradient2(low = "black", high = "blueviolet") +
  #geom_text_repel(aes(label = Well_n)) +
  theme_bw() +
  ylab("strains") +
  xlab ("TPEN_conc")
  



AUCs_mGAM_mean_2  <- AUCs_mGAM_mean  %>% group_by(Class,  conc) %>%
  summarise(AUCsmean = mean(AUCsmean, na.rm = TRUE)) 


AUCs_mGAM_mean_2$Class<- factor(AUCs_mGAM_mean_2$Class, levels = c("Gammaproteobacteria", "Clostridia", "Bacilli", "Bacteroidia", "Erysipelotrichia"))


AUCs_mGAM_mean_2 %>% 
  filter (conc != "blank") %>%
  ggplot(aes(conc, Class, fill = AUCsmean)) +
  geom_raster() +
  #geom_text(data = pvs %>% filter(qv_rescue < 0.01), label = "*") +
  scale_fill_gradient2(low = "black", high = "blueviolet") +
  #geom_text_repel(aes(label = Well_n)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 12)) +
  xlab ("TPEN_conc")


AUCs_mGAM_mean %>% write_csv2("AUCs_mGAMmean.tsv")

AUCs_mGAM_select <-  AUCs_mGAM_mean %>% filter (conc != "blank") %>% select(short, conc, AUCsmean) %>% ungroup ()

AUCs_mGAM_select_spread <-  AUCs_mGAM_select  %>% pivot_wider (names_from = short, values_from = AUCsmean) 
AUCs_mGAM_select_spread

AUCs_mGAM_select_spread_m <- as.matrix( AUCs_mGAM_select_spread[, -1], rownames.force = NA)
AUCs_mGAM_select_spread_m

rownames( AUCs_mGAM_select_spread_m) <-  AUCs_mGAM_select_spread$conc
AUCs_mGAM_select_spread_m

heatmap( AUCs_mGAM_select_spread_m, scale = "none")


heatmap.2( AUCs_mGAM_select_spread_m, scale = "none", col = bluered(100), 
           trace = "none", density.info = "none",
           margins = c(8,11),
           #ColSideColors=col_colors,
           key = TRUE,
           cexCol = 1.2,
           cexRow = 1.2)


heatmap.2( AUCs_mGAM_select_spread_m, 
           scale = "none", 
           col = colorRampPalette(c("ivory", "darkgoldenrod1", "red"))(100),
           trace = "none", 
           margins = c(10, 11),
           key = TRUE,
           cexCol = 1.2,
           cexRow = 1.2,
           breaks = seq(0, 2, length.out = 101))

