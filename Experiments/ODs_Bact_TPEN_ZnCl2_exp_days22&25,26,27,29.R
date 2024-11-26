setwd("~/ITQB/Plate readings")

library(tidyverse)
library(dbplyr) 
library(readxl)

Exps_Bact <- Exps_Bacteroides  <- read_excel("ODs_Bact_TPEN_ZnCl2_exp_days22&25,26,27,29.xlsx")
View(Exps_Bact)

Exps_Bact_v1 <- Exps_Bact %>%
  pivot_longer(names_to = "Condition", values_to = "OD_578", cols = -Time_h)
# Exps_Bact_v1
View(Exps_Bact_v1)


Exps_Bact_v2 <- separate(Exps_Bact_v1, Condition, 
                     c("Exp_day",
                       "Strain",
                       "Transfer_medium",
                       "Zn_condition", 
                       "Tech_Repl",
                       "Biol_Repl",
                       "TPEN_or_ZnCl2_conc",
                       "Well_number"), 
                     sep = "_", remove = FALSE,
                     convert = FALSE)
print(Exps_Bact_v2)
View(Exps_Bact_v2)

Exps_Bact_v2$OD_578 <- as.numeric (Exps_Bact_v2$OD_578) #Because OD_578 values are character

# table(Exps_Bact_v2$TPEN_or_ZnCl2_conc)

Exps_Bact_v2$TPEN_or_ZnCl2_conc <- factor(Exps_Bact_v2$TPEN_or_ZnCl2_conc, levels = c("ctrZnCl2",
                                                                                      "10",
                                                                                      "25",
                                                                                      "50",
                                                                                      "75",
                                                                                      "100",
                                                                                      "ctrTPEN",
                                                                                      "1.8",
                                                                                      "2.4",
                                                                                      "3.0",
                                                                                      "3.6",
                                                                                      "4.2",
                                                                                      "4.8",
                                                                                      "5.4",
                                                                                      "6.0"))


# table(Exps_Bact_v2$TPEN_or_ZnCl2_conc)

Exps_Bact_v2 <- Exps_Bact_v2 %>% mutate(TPEN_or_ZnCl2_conc = recode (TPEN_or_ZnCl2_conc,
                                                                     "ctrZnCl2" = "0.1M HCl",
                                                                     "10" = "100 uM ZnCl2",
                                                                     "25" = "250 uM ZnCl2",
                                                                     "50"= "500 uM ZnCl2",
                                                                     "75"= "750 uM ZnCl2",
                                                                     "100" = "1 mM ZnCl2",
                                                                     "ctrTPEN" = "DMSO",
                                                                     "1.8" = "18 uM TPEN",
                                                                     "2.4" = "24 uM TPEN",
                                                                     "3.0" = "30 uM TPEN",
                                                                     "3.6" = "36 uM TPEN",
                                                                     "4.2" = "42 uM TPEN",
                                                                     "4.8" = "48 uM TPEN",
                                                                     "5.4" = "54 uM TPEN",
                                                                     "6.0" = "60 uM TPEN"))

Exps_Bact_v2$Strain <- factor(Exps_Bact_v2$Strain, levels = c("5001",
                                                              "5002",
                                                              "5003",
                                                              "5004"))

# Exps_Bact_v2 %>%
#   ggplot(aes(x= Time_h, y=OD_578, color = TPEN_or_ZnCl2_conc)) +
#   geom_line() +
#   facet_grid (Strain ~ Tech_Repl + Exp_day)

Exps_Bact_v3<- Exps_Bact_v2 %>% group_by(Strain, Zn_condition, Tech_Repl, Biol_Repl, Well_number) %>%
  mutate(blank = Time_h ==0)

#View(Exps_Bact_v3)

Exps_Bact_v4 <- Exps_Bact_v3 %>%
  group_by(Strain, Transfer_medium, Zn_condition, Tech_Repl, Biol_Repl, Well_number) %>%
  mutate(OD = OD_578 - OD_578[blank]) %>%
  ungroup()

View(Exps_Bact_v4)


#Obtain the mean and sd
Exps_Bact_v5 <- Exps_Bact_v4 %>% group_by(Time_h, Strain, Zn_condition, TPEN_or_ZnCl2_conc) %>%
  summarise(ODmean = mean(OD), sdOD = sd(OD), n())

Exps_Bact_v5$Strain <- factor(Exps_Bact_v5$Strain, levels = c("5001",
                                                              "5002",
                                                              "5003",
                                                              "5004"))


Exps_Bact_v5 <- Exps_Bact_v5 %>% mutate(TPEN_or_ZnCl2_conc = recode (TPEN_or_ZnCl2_conc,
                                                                     "ctrZnCl2" = "0.1M HCl",
                                                                     "10" = "100 uM ZnCl2",
                                                                     "25" = "250 uM ZnCl2",
                                                                     "50"= "500 uM ZnCl2",
                                                                     "75"= "750 uM ZnCl2",
                                                                     "100" = "1 mM ZnCl2",
                                                                     "ctrTPEN" = "DMSO",
                                                                     "1.8" = "18 uM TPEN",
                                                                     "2.4" = "24 uM TPEN",
                                                                     "3.0" = "30 uM TPEN",
                                                                     "3.6" = "36 uM TPEN",
                                                                     "4.2" = "42 uM TPEN",
                                                                     "4.8" = "48 uM TPEN",
                                                                     "5.4" = "54 uM TPEN",
                                                                     "6.0" = "60 uM TPEN"))


# stocks conc / 0.001(mM --> uM) / d.f (=2)--> conc TPEN (now in the table)


Exps_Bact_v5  %>% ggplot(aes(x = Time_h,
             y = ODmean,
             color = TPEN_or_ZnCl2_conc,
             fill = TPEN_or_ZnCl2_conc
  )) +
  geom_point() +
  geom_line() +
  geom_ribbon(aes(ymin=ODmean-sdOD, ymax=ODmean+sdOD),
              alpha = 0.3, colour = NA) +
  facet_wrap (Strain ~ Zn_condition) +
  theme_bw()

# Assuming Exps_Bact_v4 or Exps_Bact_v5 contains your cleaned Bacteroides data
# Filter the data to exclude blanks and unnecessary data points
Exps_Bact_TPEN <- Exps_Bact_v5 
#-------------Zn_condition == "woZn"-----------
# Ensure OD values are numeric and calculate AUC
Exps_Bact_TPEN$ODmean <- as.numeric(Exps_Bact_TPEN$ODmean)

Exps_Bact_v4_wo_Zn <- Exps_Bact_v4 %>% filter(Zn_condition == "woZn") %>% mutate(control = TPEN_or_ZnCl2_conc == "DMSO", ODc0 = OD)
table(Exps_Bact_v4_wo_Zn$Biol_Repl)

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

Exps_Bact_v4_wo_Zn <- Exps_Bact_v4_wo_Zn %>% mutate(f_experiment_number = as.factor(Exp_day))

time_cutoffs <- determineTimeCutoffs(Exps_Bact_v4_wo_Zn) %>%
  mutate(f_experiment_number = as.factor(Exp_day))

# print(time_cutoffs)

Exps_Bact_v4_wo_Zn_control <- Exps_Bact_v4_wo_Zn %>% 
  filter(control) %>% 
  group_by(f_experiment_number, Exp_day, Strain, Time_h) %>%
  summarise(OD = median(OD))

Exps_Bact_v4_wo_Zn <- Exps_Bact_v4_wo_Zn %>% left_join(time_cutoffs) %>% filter(Time_h <= cutoff_time) %>% 
  select(-cutoff_time)

calcAUC <- function(OD, time) {
  OD <- OD - min(OD)
  N <- length(OD)
  if (N<2) return(NA)
  s1 <- 1:(N-1)
  s2 <- s1 + 1
  sum( (OD[s1] + OD[s2])/2 * (time[s2]-time[s1]) )
}

AUCs_M8_wo_Zn <- Exps_Bact_v4_wo_Zn %>% group_by(Strain,Exp_day,Biol_Repl, Tech_Repl,
                                        TPEN_or_ZnCl2_conc, control, Well_number) %>% 
  summarise(AUC = calcAUC(ODc0, Time_h), maxOD = max(ODc0), finalOD = ODc0[n()])


AUCs_M8_wo_Zn <-  AUCs_M8_wo_Zn %>% group_by(Strain,Exp_day) %>% 
  mutate(AUC = AUC / median(AUC[control], na.rm = TRUE))


AUCs_M8_mean_wo_Zn <-  AUCs_M8_wo_Zn %>%  group_by(Strain, TPEN_or_ZnCl2_conc) %>% 
  summarise (AUCsmean = mean(AUC), sdAUC = sd (AUC), maxODmean = mean (maxOD), sdmax = sd (maxOD), finalODmean = mean(finalOD), sdfinal = sd(finalOD)) %>% ungroup()

AUCs_M8_mean_wo_Zn <- AUCs_M8_mean_wo_Zn %>% mutate (Strain = recode (Strain,
                                                                      "5001" = "B. vulgatus",
                                                                      "5002" = "B. uniformis",
                                                                      "5003" = "B. fragilis",
                                                                      "5004" = "B. thetaiotaomicron"))

AUCs_M8_mean_wo_Zn$TPEN_or_ZnCl2_conc <- factor (AUCs_M8_mean_wo_Zn$TPEN_or_ZnCl2_conc, levels = c(
  "1 mM ZnCl2", "750 uM ZnCl2", "500 uM ZnCl2", "250 uM ZnCl2", "100 uM ZnCl2", "0.1M HCl", "DMSO", "18 uM TPEN", "24 uM TPEN", "30 uM TPEN", "36 uM TPEN", "42 uM TPEN", "48 uM TPEN", "54 uM TPEN", "60 uM TPEN"))

# Create the heatmap
AUCs_M8_mean_wo_Zn %>% filter (Strain != "B. fragilis") %>% #B. fragilis presented not consistent results
  ggplot(aes(x = TPEN_or_ZnCl2_conc, y = Strain, fill = AUCsmean)) +
  geom_raster() +
  scale_fill_gradient2(low = "black", high = "blueviolet") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12),
    legend.position = "bottom",
    legend.text = element_text(size = 8)
  ) +
  xlab("drug_conc") +
  ylab("Specie")


#-------------Zn_condition == "Zn"-----------
# Ensure OD values are numeric and calculate AUC
Exps_Bact_TPEN$ODmean <- as.numeric(Exps_Bact_TPEN$ODmean)

Exps_Bact_v4_Zn <- Exps_Bact_v4 %>% filter(Zn_condition == "Zn") %>% mutate(control = TPEN_or_ZnCl2_conc == "DMSO", ODc0 = OD)
table(Exps_Bact_v4_Zn$Biol_Repl)

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

Exps_Bact_v4_Zn <- Exps_Bact_v4_Zn %>% mutate(f_experiment_number = as.factor(Exp_day))

time_cutoffs <- determineTimeCutoffs(Exps_Bact_v4_Zn) %>%
  mutate(f_experiment_number = as.factor(Exp_day))

# print(time_cutoffs)

Exps_Bact_v4_Zn_control <- Exps_Bact_v4_Zn %>% 
  filter(control) %>% 
  group_by(f_experiment_number, Exp_day, Strain, Time_h) %>%
  summarise(OD = median(OD))

Exps_Bact_v4_Zn <- Exps_Bact_v4_Zn %>% left_join(time_cutoffs) %>% filter(Time_h <= cutoff_time) %>% 
  select(-cutoff_time)

calcAUC <- function(OD, time) {
  OD <- OD - min(OD)
  N <- length(OD)
  if (N<2) return(NA)
  s1 <- 1:(N-1)
  s2 <- s1 + 1
  sum( (OD[s1] + OD[s2])/2 * (time[s2]-time[s1]) )
}

AUCs_M8_Zn <- Exps_Bact_v4_Zn %>% group_by(Strain,Exp_day,Biol_Repl, Tech_Repl,
                                                 TPEN_or_ZnCl2_conc, control, Well_number) %>% 
  summarise(AUC = calcAUC(ODc0, Time_h), maxOD = max(ODc0), finalOD = ODc0[n()])


AUCs_M8_Zn <-  AUCs_M8_Zn %>% group_by(Strain,Exp_day) %>% 
  mutate(AUC = AUC / median(AUC[control], na.rm = TRUE))


AUCs_M8_mean_Zn <-  AUCs_M8_Zn %>%  group_by(Strain, TPEN_or_ZnCl2_conc) %>% 
  summarise (AUCsmean = mean(AUC), sdAUC = sd (AUC), maxODmean = mean (maxOD), sdmax = sd (maxOD), finalODmean = mean(finalOD), sdfinal = sd(finalOD)) %>% ungroup()

AUCs_M8_mean_Zn <- AUCs_M8_mean_Zn %>% mutate (Strain = recode (Strain,
                                                                "5001" = "B. vulgatus",
                                                                "5002" = "B. uniformis",
                                                                "5003" = "B. fragilis",
                                                                "5004" = "B. thetaiotaomicron"))

AUCs_M8_mean_Zn$TPEN_or_ZnCl2_conc <- factor (AUCs_M8_mean_Zn$TPEN_or_ZnCl2_conc, levels = c(
  "1 mM ZnCl2", "750 uM ZnCl2", "500 uM ZnCl2", "250 uM ZnCl2", "100 uM ZnCl2", "0.1M HCl", "DMSO", "18 uM TPEN", "24 uM TPEN", "30 uM TPEN", "36 uM TPEN", "42 uM TPEN", "48 uM TPEN", "54 uM TPEN", "60 uM TPEN"))


# Create the heatmap
AUCs_M8_mean_Zn %>% filter (Strain != "B. fragilis") %>% #B. fragilis presented not consistent results
  ggplot(aes(x = TPEN_or_ZnCl2_conc, y = Strain, fill = AUCsmean)) +
  geom_raster() +
  scale_fill_gradient2(low = "black", high = "blueviolet") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12),
    legend.position = "bottom",
    legend.text = element_text(size = 10)
  ) +
  xlab("drug_conc") +
  ylab("Specie")
