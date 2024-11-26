setwd("~/ITQB/Plate readings/Exp26, 27, 28 - 28022024, 28022024, 01032024")

library(tidyverse)

library(readxl)
Exp26 <- Exp26_5plates_28022024 <- read_excel("Exp26_5plates_28022024_Bact_final_final.xlsx")
View(Exp26) 

Exp26_v1 <- Exp26 %>%
  pivot_longer(names_to = "Condition", values_to = "OD_578", cols = -Time_h)
Exp26_v1
View(Exp26_v1)


Exp26_v2 <- separate(Exp26_v1, Condition, 
                     c("Strain",
                       "Transfer_medium",
                       "Zn_condition", 
                       "Tech_Repl",
                       "Biol_Repl",
                       "TPEN_or_ZnCl2_conc",
                       "Well_number"), 
                     sep = "_", remove = FALSE,
                     convert = FALSE)
View(Exp26_v2)

Exp26_v2$OD_578 <- as.numeric (Exp26_v2$OD_578) #Because OD_578 values are character

Exp26_v2$TPEN_or_ZnCl2_conc <- factor(Exp26_v2$TPEN_or_ZnCl2_conc, levels = c("ctrZnCl2",
                                                            "10",
                                                            "25",
                                                            "ctrTPEN",
                                                            "1.8",
                                                            "2.4",
                                                            "3.0",
                                                            "3.6",
                                                            "4.2",
                                                            "4.8",
                                                            "5.4",
                                                            "6.0"))

Exp26_v2 <- Exp26_v2 %>% mutate(TPEN_or_ZnCl2_conc = recode (TPEN_or_ZnCl2_conc,
                                                    "ctrZnCl2" = "0.1HCl",
                                                    "10" = "100 uM ZnCl2",
                                                    "25" = "250 uM ZnCl2",
                                                    "ctrTPEN" = "DMSO",
                                                    "1.8" = "18 uM TPEN",
                                                    "2.4" = "24 uM TPEN",
                                                    "3.0" = "30 uM TPEN",
                                                    "3.6" = " 36 uM TPEN",
                                                    "4.2" = "42 uM TPEN",
                                                    "4.8" = "48 uM TPEN",
                                                    "5.4" = "54 uM TPEN",
                                                    "6.0" = "60 uM TPEN"))

Exp26_v2$Strain <- factor(Exp26_v2$Strain, levels = c("5001",
                                                    "5002",
                                                    "5003",
                                                    "5004"))

# Exp26_v2 %>%
#   ggplot(aes(x= Time_h, y=OD_578, color = Tech_Repl)) +
#   geom_line() +
#   facet_grid (Strain ~ TPEN_or_ZnCl2_conc)

#I don't need this plot because at this experiment I only had M8 medium.
# Exp26_v2 %>% filter (Transfer_medium == "mGAM") %>%
#   ggplot(aes(x= Time_h, y=OD_578, color = Tech_Repl)) +
#   geom_line() +
#   facet_grid (Strain ~ TPEN_or_ZnCl2_conc)


Exp26_v2 %>% filter (Transfer_medium == "M8") %>%
  ggplot(aes(x= Time_h, y=OD_578, color = Tech_Repl)) +
  geom_line() +
  facet_grid (Strain ~ TPEN_or_ZnCl2_conc)


Exp26_v3<- Exp26_v2 %>% group_by(Strain, Transfer_medium, Zn_condition, Tech_Repl, Biol_Repl, Well_number) %>%
  mutate(blank = Time_h ==0)

View(Exp26_v3)

Exp26_v4 <- Exp26_v3 %>%
  group_by(Strain, Transfer_medium, Zn_condition, Tech_Repl, Biol_Repl, Well_number) %>%
  mutate(OD = OD_578 - OD_578[blank]) %>%
  ungroup()


View(Exp26_v4)

#Obtain the mean and sd
Exp26_v5 <- Exp26_v4 %>% group_by(Time_h, Strain, Transfer_medium, TPEN_or_ZnCl2_conc) %>%
  summarise(ODmean = mean(OD), sdOD = sd(OD), n())

Exp26_v5$Strain <- factor(Exp26_v5$Strain, levels = c("5001",
                                                      "5002",
                                                      "5003",
                                                      "5004"))


Exp26_v5 <- Exp26_v5 %>% mutate(TPEN_or_ZnCl2_conc =  recode (TPEN_or_ZnCl2_conc,
                                                              "ctrZnCl2" = "0.1HCl",
                                                              "10" = "100 uM ZnCl2",
                                                              "25" = "250 uM ZnCl2",
                                                              "ctrTPEN" = "DMSO",
                                                              "1.8" = "18 uM TPEN",
                                                              "2.4" = "24 uM TPEN",
                                                              "3.0" = "30 uM TPEN",
                                                              "3.6" = " 36 uM TPEN",
                                                              "4.2" = "42 uM TPEN",
                                                              "4.8" = "48 uM TPEN",
                                                              "5.4" = "54 uM TPEN",
                                                              "6.0" = "60 uM TPEN"))



# stocks conc / 0.001(mM --> uM) / d.f (=2)--> conc TPEN (now in the table)


Exp26_v5  %>% ggplot(aes(x = Time_h,
             y = ODmean,
             color = Transfer_medium
  )) +
  geom_point() +
  geom_line() +
  geom_ribbon(aes(ymin=ODmean-sdOD, ymax=ODmean+sdOD, fill = Transfer_medium),
              alpha = 0.3, colour = NA) +
  facet_grid (Strain ~ TPEN_or_ZnCl2_conc) +
  theme_bw()

#Filter to have just 1 strain in the plot
filtered_strain <- Exp26_v5 %>% filter (Strain == "5001")

filtered_strain <- Exp26_v5 %>% filter (Strain == "5002")

filtered_strain <- Exp26_v5 %>% filter (Strain == "5003")

filtered_strain <- Exp26_v5 %>% filter (Strain == "5004")

#Filter to compare the blanks (req_conc vs blank from the other stocks - 100 uM (ZnCl2 stock) vs TPEN blank) with a required concentration from a certain strain in the plot
filtered_data <- filtered_strain %>%
  filter(TPEN_or_ZnCl2_conc %in% c("100 uM", "blank"))

filtered_data <- filtered_strain %>%
  filter(TPEN_or_ZnCl2_conc %in% c("250 uM", "blank"))

filtered_data <- filtered_strain %>%
  filter(TPEN_or_ZnCl2_conc %in% c("18 uM", "0.1HCl"))

filtered_data <- filtered_strain %>%
  filter(TPEN_or_ZnCl2_conc %in% c("24 uM", "0.1HCl"))

filtered_data <- filtered_strain %>%
  filter(TPEN_or_ZnCl2_conc %in% c("30 uM", "0.1HCl"))

filtered_data <- filtered_strain %>%
  filter(TPEN_or_ZnCl2_conc %in% c("36 uM", "0.1HCl"))

filtered_data <- filtered_strain %>%
  filter(TPEN_or_ZnCl2_conc %in% c("42 uM", "0.1HCl"))

filtered_data <- filtered_strain %>%
  filter(TPEN_or_ZnCl2_conc %in% c("48 uM", "0.1HCl"))

filtered_data <- filtered_strain %>%
  filter(TPEN_or_ZnCl2_conc %in% c("54 uM", "0.1HCl"))

filtered_data <- filtered_strain %>%
  filter(TPEN_or_ZnCl2_conc %in% c("60 uM", "0.1HCl"))


ggplot(filtered_data, aes(x = Time_h, y = ODmean, color = Transfer_medium)) +
  geom_point() +
  geom_line() +
  geom_ribbon(aes(ymin = ODmean - sdOD, ymax = ODmean + sdOD, fill = Transfer_medium),
              alpha = 0.3, color = NA) +
  facet_wrap(~ Strain + TPEN_or_ZnCl2_conc) +
  theme_bw()
