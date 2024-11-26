setwd("~/ITQB/Plate readings")

library(tidyverse)
library(dbplyr) 
library(readxl)

Exps <- Exps_26_27_28  <- read_excel("ODs_control_TPENvsZnCl2_exp_days__R.xlsx")
# View(Exps) 

Exps_v1 <- Exps %>%
  pivot_longer(names_to = "Condition", values_to = "OD_578", cols = -Time_h)
# Exps_v1
# View(Exps_v1)


Exps_v2 <- separate(Exps_v1, Condition, 
                     c("Exp_day",
                       "Strain",
                       "Transfer_medium",
                       "Zn_condition", 
                       "Tech_Repl",
                       "Biol_Repl",
                       "drug_conc",
                       "Well_number"), 
                     sep = "_", remove = FALSE,
                     convert = FALSE)
# View(Exps_v2)

Exps_v2$OD_578 <- as.numeric (Exps_v2$OD_578) #Because OD_578 values are character

Exps_v2$drug_conc <- factor(Exps_v2$drug_conc, levels = c("ctr",
                                                          "20T",
                                                          "40Zn",
                                                          "20T40Zn"))

Exps_v2 <- Exps_v2 %>% mutate(drug_conc = recode (drug_conc,
                                                             "ctr" = "untreated",
                                                             "20T" = "100 uM TPEN",
                                                             "40Zn" = "200 uM ZnCl2",
                                                             "20T40Zn" = "100 uM TPEN + 200 uM ZnCl2"))

Exps_v2$Strain <- factor(Exps_v2$Strain, levels = c("5001",
                                                    "5002",
                                                    "5003",
                                                    "5004",
                                                    "5006",
                                                    "5026",
                                                    "5032",
                                                    "5037",
                                                    "5038",
                                                    "5046",
                                                    "5048",
                                                    "5076",
                                                    "5078",
                                                    "23006",
                                                    "HM02",
                                                    "HM04"))

Exps_v2$Exp_day <- factor(Exps_v2$Exp_day, levels = c("Exp26",
                                                      "Exp27",
                                                      "Exp28"))

Exps_v2 %>%
  ggplot(aes(x= Time_h, y=OD_578, color = drug_conc)) +
  geom_line() +
  facet_grid (Strain ~ Exp_day + Tech_Repl)
#5002 doesn't appear in Exp27, because it doesn't grow in mGAM --> I couldn't use it 

Exps_v3<- Exps_v2 %>% group_by(Strain, Transfer_medium, Zn_condition, Tech_Repl, Biol_Repl, Well_number) %>%
  mutate(blank = Time_h ==0)

#View(Exps_v3)

Exps_v4 <- Exps_v3 %>%
  group_by(Strain, Transfer_medium, Zn_condition, Tech_Repl, Biol_Repl, Well_number) %>%
  mutate(OD = OD_578 - OD_578[blank]) %>%
  ungroup()

#To filter the strains didn't grow in one or more Exp day
Exps_v4 <- Exps_v4  %>%  filter (
  !(Exp_day == "Exp26" & Strain == "5048"), #doesn't grow
  !(Exp_day == "Exp26" & Strain == "5006"), #doesn't grow
  !(Strain %in% c("5076", "23006"))) #didn't grow

#To check whether filters were ok
Exps_v4 %>%
  ggplot(aes(x= Time_h, y=OD_578, color = drug_conc)) +
  geom_line() +
  facet_grid (Strain ~ Exp_day + Tech_Repl)

# View(Exps_v4)

#Obtain the mean and sd
Exps_v5 <- Exps_v4 %>% group_by(Time_h, Strain, Transfer_medium, drug_conc) %>%
  summarise(ODmean = mean(OD), sdOD = sd(OD), n())

Exps_v5$Strain <- factor(Exps_v5$Strain, levels = c("5001",
                                                      "5002",
                                                      "5003",
                                                      "5004",
                                                      "5006",
                                                       "5026",
                                                      "5032",
                                                      "5037",
                                                      "5038",
                                                      "5046",
                                                      "5048",
                                                      "5076",
                                                      "5078",
                                                      "23006",
                                                      "HM02",
                                                      "HM04"))


Exps_v5 <- Exps_v5 %>% mutate(drug_conc = recode (drug_conc,
                                                             "ctr" = "untreated",
                                                             "20T" = "100 uM TPEN",
                                                             "40Zn" = "200 uM ZnCl2",
                                                             "20T40Zn" = "100 uM TPEN + 200 uM ZnCl2"))



Exps_v5  %>% ggplot(aes(x = Time_h,
             y = ODmean,
             color = drug_conc,
             fill = drug_conc
  )) +
  geom_point() +
  geom_line() +
  geom_ribbon(aes(ymin=ODmean-sdOD, ymax=ODmean+sdOD),
              alpha = 0.3, colour = NA) +
  facet_wrap (Strain ~ .) +
  theme_bw()

# Assuming Exps_v5 contains your cleaned data
Exps_v6 <- Exps_v4 %>% mutate(control = drug_conc == "200 uM ZnCl2", ODc0 = OD)

cumulmin <- function(v) {
  N <- length(v)
  if (N==0) return()
  sapply(1:N, function(i)
    min(v[i:N])
  )
}

determineTimeCutoffs <- function(GCs, decrease = 0.02) {
  GCs_controls <- GCs %>% filter(control) %>% mutate(Time_h = round(Time_h, 0)) %>%
    group_by(Strain, Exp_day, Time_h) %>%
    summarise(ODc0 = median(ODc0)) %>% arrange(Time_h) %>%
    group_by(Strain, Exp_day) %>% mutate(ODc0_adj = ODc0 + ODc0[n()]*decrease*(Time_h[n()]-Time_h))
  
  GCs_controls %>% group_by(Strain, Exp_day) %>% filter(ODc0_adj == max(ODc0_adj)) %>%
    ungroup() %>% select(Strain, Exp_day, Time_h) %>% rename(cutoff_time = Time_h)
}

Exps_v6 <- Exps_v6 %>% mutate(f_experiment_number = as.factor(Exp_day))

time_cutoffs <- determineTimeCutoffs(Exps_v6) %>%
  mutate(f_experiment_number = as.factor(Exp_day))

# print(time_cutoffs)

Exps_v6_control <- Exps_v6 %>% 
  filter(control) %>% 
  group_by(f_experiment_number, Exp_day, Strain, Time_h) %>%
  summarise(OD = median(OD))

Exps_v6 <- Exps_v6 %>% left_join(time_cutoffs) %>% filter(Time_h <= cutoff_time) %>% 
  select(-cutoff_time)

calcAUC <- function(OD, time) {
  OD <- OD - min(OD)
  N <- length(OD)
  if (N<2) return(NA)
  s1 <- 1:(N-1)
  s2 <- s1 + 1
  sum( (OD[s1] + OD[s2])/2 * (time[s2]-time[s1]) )
}

AUCs <- Exps_v6 %>% group_by(Strain,Exp_day,Biol_Repl, Tech_Repl,
                                                 drug_conc, control, Well_number) %>% 
  summarise(AUC = calcAUC(ODc0, Time_h), maxOD = max(ODc0), finalOD = ODc0[n()])


AUCs <-  AUCs %>% group_by(Strain,Exp_day) %>% 
  mutate(AUC = AUC / median(AUC[control], na.rm = TRUE))


AUCs_mean <-  AUCs %>%  group_by(Strain, drug_conc,Exp_day) %>% 
  summarise (AUCsmean = mean(AUC), sdAUC = sd (AUC), maxODmean = mean (maxOD), sdmax = sd (maxOD), finalODmean = mean(finalOD), sdfinal = sd(finalOD)) %>% ungroup()

AUCs_mean <- AUCs_mean %>% mutate(Strain = recode (Strain,
                                               "5001" = "B. vulgatus",
                                               "5002" = "B. uniformis",
                                               "5003" = "B. fragilis",
                                               "5004" = "B. thetaiotaomicron",
                                               "5006" = "T. ramosa",
                                               "5026" = " E. bolteae",
                                               "5032" = "C. perfringens",
                                               "5037" = "L. saccharolytica",
                                               "5038" = "S. salivarius",
                                               "5046" = "M. gnavus",
                                               "5048" = "C. comes",
                                               "5076" = "D. formicigenerans",
                                               "5078" = "E. coli ED1a",
                                               "23006" = "C. difficile",
                                               "HM02" = "K. pneumoniae",
                                               "HM04" = "E. faecalis"))


AUCs_mean$drug_conc <- factor (AUCs_mean$drug_conc, levels = c(
"untreated", "100 uM TPEN", "200 uM ZnCl2", "100 uM TPEN + 200 uM ZnCl2"))

# Create the heatmap
AUCs_mean %>% filter (Exp_day != "Exp28") %>%
  ggplot(aes(x = drug_conc, y = Strain, fill = AUCsmean)) +
  geom_raster() +
  scale_fill_gradient2(low = "black", high = "blueviolet") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 9, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    legend.text = element_text(size = 10)
  ) +
  xlab("drug_conc") +
  ylab("Specie/Strain")
