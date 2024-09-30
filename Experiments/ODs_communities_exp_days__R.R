setwd("~/ITQB/Plate readings")

library(tidyverse)
library(dbplyr) 
library(readxl)

Exps <- Exps_30_42  <- read_excel("ODs_communities_exp_days.xlsx")
# View(Exps)

Exps_v1 <- Exps %>%
  pivot_longer(names_to = "Condition", values_to = "OD_578", cols = -Time_h)
# Exps_v1
# View(Exps_v1)


Exps_v2 <- separate(Exps_v1, Condition, 
                     c("Exp_day",
                       "Community",
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


Exps_v2$Community <- factor(Exps_v2$Community, levels = c("control",
                                                            "C1",
                                                            "C2",
                                                            "C3",
                                                            "C4",
                                                            "C5"))

Exps_v2$drug_conc <- factor(Exps_v2$drug_conc, levels = c("nodrug",
                                                          "1.74",
                                                          "10",
                                                          "0",
                                                          "1.8",
                                                          "2.4",
                                                          "3.0",
                                                          "3.6",
                                                          "4.2",
                                                          "4.8",
                                                          "5.4",
                                                          "6.0",
                                                          "7.0",
                                                          "8.0"))

Exps_v2 <- Exps_v2 %>% mutate(drug_conc =  recode (drug_conc,
                                                   "nodrug" = "blank",
                                                   "1.74" = "17.4 uM Zn",
                                                   "10" = "50 uM Zn",
                                                   "0" = "0 uM TPEN",
                                                   "1.8" = "9 uM TPEN",
                                                   "2.4" = "12 uM TPEN",
                                                   "3.0" = "15 uM TPEN",
                                                   "3.6" = "18 uM TPEN",
                                                   "4.2" = "21 uM TPEN",
                                                   "4.8" = "24 uM TPEN",
                                                   "5.4" = "27 uM TPEN",
                                                   "6.0" = "30 uM TPEN",
                                                   "7.0" = "35 uM TPEN",
                                                   "8.0" = "40 uM TPEN"))

Exps_v2 <- Exps_v2 %>% mutate(Community =  recode (Community,
                                                   "control" = "blank"))

Exps_v2 <- Exps_v2 %>% mutate(Tech_Repl =  recode (Tech_Repl,
                                                   "TechRep1" = "TR1",
                                                   "TechRep2" = "TR2",
                                                   "TechRep3" = "TR3"))

Exps_v2$Exp_day <- factor(Exps_v2$Exp_day, levels = c("Exp30",
                                                      "Exp31",
                                                      "Exp32",
                                                      "Exp33",
                                                      "Exp34",
                                                      "Exp35",
                                                      "Exp36",
                                                      "Exp37",
                                                      "Exp38",
                                                      "Exp39",
                                                      "Exp40",
                                                      "Exp41",
                                                      "Exp42"))



filtered_data <- Exps_v2 %>%
  filter(!(Exp_day == "Exp38" & Community == "C2" & Tech_Repl == "TechRep2" & drug_conc == "15uM"),
         !(Exp_day %in% c("Exp30", "Exp31", "Exp32", "Exp33", "Exp34", "Exp35", "Exp36", "Exp37", "Exp38", "Exp39")),
         !(Exp_day == "Exp40" & drug_conc == "9uM"))

# Exp30 presented very strange results, growth curves didn't make sense
# Exp31 presented very high growth probably thanks to contamination of communities
# Exp36 presented growth curves very similar on the beginning (different from other experiments, so we should discard these results)
# Exp40 would present different growth curves than expected because some wells of 96-well plates didn't have the correct volume (A3 and E3), volume wasn't enough), so I'm filtering this concentration
# I selected only the last Exp days, because results seemed to be kind of consistent
## (Communities grew better under Zn conditions than when Zn was chelated by TPEN)

Exps_v2 %>%
  ggplot(aes(x= Time_h, y=OD_578, color = drug_conc)) +
  geom_line() +
  facet_grid (Community ~ Exp_day + Tech_Repl) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 3),  # Decrease the size of the numbers on the axes
    axis.text.y = element_text(size = 5),
    strip.text.x = element_text(size = 7)
)

filtered_data %>%
  ggplot(aes(x= Time_h, y=OD_578, color = drug_conc)) +
  geom_line() +
  facet_grid (Community ~ Exp_day + Tech_Repl) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6))

Exps_v3 <- Exps_v2 %>% group_by(Community, Transfer_medium, Zn_condition, Tech_Repl, Biol_Repl, Well_number) %>%
  mutate(blank = Time_h ==0)

Exps_v3 <- filtered_data %>% group_by(Community, Transfer_medium, Zn_condition, Tech_Repl, Biol_Repl, Well_number) %>%
  mutate(blank = Time_h ==0)

# View(Exps_v3)

Exps_v4 <- Exps_v3 %>%
  group_by(Community, Transfer_medium, Zn_condition, Tech_Repl, Biol_Repl, Well_number) %>%
  mutate(OD = OD_578 - OD_578[blank]) %>%
  ungroup()


# View(Exps_v4)

#Obtain the mean and sd
Exps_v5 <- Exps_v4 %>% group_by(Time_h, Community, Transfer_medium, drug_conc) %>%
  summarise(ODmean = mean(OD), sdOD = sd(OD), n())

Exps_v5 <- Exps_v5 %>% mutate(drug_conc =  recode (drug_conc,
                                                   "17.4 uM" = "17.4 uM ZnCl2",
                                                   "50uM" = "50 uM ZnCl2",
                                                   "0uM" = "0 uM TPEN",
                                                   "9uM" = "9 uM TPEN",
                                                   "12uM" = "12 uM TPEN",
                                                   "15uM" = "15 uM TPEN",
                                                   "18uM" = "18 uM TPEN",
                                                   "21uM" = "21 uM TPEN",
                                                   "24uM" = "24 uM TPEN",
                                                   "27uM" = "27 uM TPEN",
                                                   "30uM" = "30 uM TPEN",
                                                   "35uM" = "35 uM TPEN",
                                                   "40uM" = "40 uM TPEN"))

Exps_v5  %>% ggplot(aes(x = Time_h,
             y = ODmean,
             color = drug_conc,
             fill = drug_conc
  )) +
  geom_point() +
  geom_line() +
  geom_ribbon(aes(ymin=ODmean-sdOD, ymax=ODmean+sdOD),
              alpha = 0.1, 
              colour = NA) +
  facet_wrap (Community ~ .) +
  theme_bw()

# Assuming Exps_v4 or Exps_v5 contains your cleaned commmunity's data
# Filter the data to exclude blanks and unnecessary data points

# Ensure OD values are numeric and calculate AUC
Exps_v5$ODmean <- as.numeric(Exps_v5$ODmean)

#To normalize according to 8.7 uM Zn condition
Exps_v4 <- Exps_v4 %>% filter (Community != "blank", drug_conc != "blank") %>% mutate(control = drug_conc == "0 uM TPEN", ODc0 = OD)

# table(Exps_v4$Biol_Repl)

#To determine the function
cumulmin <- function(v) {
  N <- length(v)
  if (N==0) return()
  sapply(1:N, function(i)
    min(v[i:N])
  )
}

#To AUC calculation (including normalization)
determineTimeCutoffs <- function(GCs, decrease = 0.02) {
  GCs_controls <- GCs %>% filter(control) %>% mutate(Time_h = round(Time_h, 0)) %>%
    group_by(Community, Exp_day, Time_h) %>%
    summarise(ODc0 = median(ODc0)) %>% arrange(Time_h) %>%
    group_by(Community, Exp_day) %>% mutate(ODc0_adj = ODc0 + ODc0[n()]*decrease*(Time_h[n()]-Time_h))
  
  GCs_controls %>% group_by(Community, Exp_day) %>% filter(ODc0_adj == max(ODc0_adj)) %>%
    ungroup() %>% select(Community, Exp_day, Time_h) %>% rename(cutoff_time = Time_h)
}

Exps_v4 <- Exps_v4 %>% mutate(f_experiment_number = as.factor(Exp_day))

# time_cutoffs <- determineTimeCutoffs(Exps_v4) %>%
#   mutate(f_experiment_number = as.factor(Exp_day))
# print(time_cutoffs)

Exps_v4_control <- Exps_v4 %>% 
  filter(control) %>% 
  group_by(f_experiment_number, Exp_day, Community, Time_h) %>%
  summarise(OD = median(OD))

# Exps_v4 <- Exps_v4 %>% left_join(time_cutoffs) %>% filter(Time_h <= cutoff_time) %>% 
#   select(-cutoff_time)

Exps_v4 <- Exps_v4 %>% mutate (cutoff_time = 20) %>% filter (Time_h <= cutoff_time) %>% select(-cutoff_time)

calcAUC <- function(OD, time) {
  OD <- OD - min(OD)
  N <- length(OD)
  if (N<2) return(NA)
  s1 <- 1:(N-1)
  s2 <- s1 + 1
  sum( (OD[s1] + OD[s2])/2 * (time[s2]-time[s1]) )
}

AUCs <- Exps_v4 %>% group_by(Community,Exp_day,Biol_Repl, Tech_Repl,
                                                 drug_conc, control, Well_number) %>% 
  summarise(AUC = calcAUC(ODc0, Time_h), maxOD = max(ODc0), finalOD = ODc0[n()])


AUCs <-  AUCs %>% group_by(Community,Exp_day) %>% 
  mutate(AUC = AUC / median(AUC[control], na.rm = TRUE))


AUCs_mean <-  AUCs %>%  group_by(Community, drug_conc) %>% 
  summarise (AUCsmean = mean(AUC), sdAUC = sd (AUC), maxODmean = mean (maxOD), sdmax = sd (maxOD), finalODmean = mean(finalOD), sdfinal = sd(finalOD)) %>% ungroup()

AUCs_mean$drug_conc <- factor (AUCs_mean$drug_conc, levels = c(
  "blank", "17.4 uM Zn", "50 uM Zn", "0 uM TPEN", "9 uM TPEN", "12 uM TPEN", "15 uM TPEN", "18 uM TPEN",
  "21 uM TPEN", "24 uM TPEN", "27 uM TPEN", "30 uM TPEN", "35 uM TPEN", "40 uM TPEN"))

AUCs_mean$Community <- factor (AUCs_mean$Community, levels = c(
  "C1", "C2", "C3", "C4", "C5"))
  
# Create the heatmap
AUCs_mean %>%
  ggplot(aes(x = drug_conc, y = Community, fill = AUCsmean)) +
  geom_raster() +
  scale_fill_gradient2(low = "black", high = "blueviolet") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12),
    legend.position = "bottom",
    legend.text = element_text(size = 8)
  )

AUCs_mean   %>%
  ggplot(aes(drug_conc, Community, fill = AUCsmean)) +
  geom_raster() +
  geom_text(aes(label = sprintf("%.2f", AUCsmean)), color = "black", size = 3) +
  scale_fill_gradientn(colors = c("ivory", "darkgoldenrod1", "blue"), values = c(0, 0.4, 1.4)) +
  theme_bw() +
    theme(
      axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 12),
      legend.text = element_text(size = 8)
  )
#Comment
#Comment 

