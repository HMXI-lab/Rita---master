setwd("~/ITQB/Plate readings/Exp26, 27, 28 - 28022024, 28022024, 01032024")

library(tidyverse)

library(readxl)
Exp27 <- Exp27_5plates_29022024 <- read_excel("Exp27_5plates_29022024_TPEN_ZnCl2_final_final.xlsx")
View(Exp27) 

Exp27_v1 <- Exp27 %>%
  pivot_longer(names_to = "Condition", values_to = "OD_578", cols = -Time_h)
Exp27_v1
View(Exp27_v1)


Exp27_v2 <- separate(Exp27_v1, Condition, 
                     c("Strain",
                       "Transfer_medium",
                       "Zn_condition", 
                       "Tech_Repl",
                       "Biol_Repl",
                       "drug_conc",
                       "Well_number"), 
                     sep = "_", remove = FALSE,
                     convert = FALSE)
View(Exp27_v2)

Exp27_v2$OD_578 <- as.numeric (Exp27_v2$OD_578) #Because OD_578 values are character

Exp27_v2$drug_conc <- factor(Exp27_v2$drug_conc, levels = c("ctr",
                                                            "20T",
                                                            "40Zn",
                                                            "20T40Zn"))

Exp27_v2 <- Exp27_v2 %>% mutate(drug_conc = recode (drug_conc,
                                                             "ctr" = "untreated",
                                                             "20T" = "100 uM TPEN",
                                                             "40Zn" = "200 uM ZnCl2",
                                                             "20T40Zn" = "100 uM TPEN + 200 uM ZnCl2"))

Exp27_v2$Strain <- factor(Exp27_v2$Strain, levels = c("5001",
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

Exp27_v2 %>%
  ggplot(aes(x= Time_h, y=OD_578, color = drug_conc)) +
  geom_line() +
  facet_grid (Strain ~ Tech_Repl)

Exp27_v3<- Exp27_v2 %>% group_by(Strain, Transfer_medium, Zn_condition, Tech_Repl, Biol_Repl, Well_number) %>%
  mutate(blank = Time_h ==0)

View(Exp27_v3)

Exp27_v4 <- Exp27_v3 %>%
  group_by(Strain, Transfer_medium, Zn_condition, Tech_Repl, Biol_Repl, Well_number) %>%
  mutate(OD = OD_578 - OD_578[blank]) %>%
  ungroup()


View(Exp27_v4)

#Obtain the mean and sd
Exp27_v5 <- Exp27_v4 %>% group_by(Time_h, Strain, Transfer_medium, drug_conc) %>%
  summarise(ODmean = mean(OD), sdOD = sd(OD), n())

Exp27_v5$Strain <- factor(Exp27_v5$Strain, levels = c("5001",
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


Exp27_v5 <- Exp27_v5 %>% mutate(drug_conc = recode (drug_conc,
                                                             "ctr" = "untreated",
                                                             "20T" = "100 uM TPEN",
                                                             "40Zn" = "200 uM ZnCl2",
                                                             "20T40Zn" = "100 uM TPEN + 200 uM ZnCl2"))


# stocks conc / 0.001(mM --> uM) / d.f (=2)--> conc TPEN (now in the table)


Exp27_v5  %>% ggplot(aes(x = Time_h,
             y = ODmean,
             color = drug_conc
  )) +
  geom_point() +
  geom_line() +
  geom_ribbon(aes(ymin=ODmean-sdOD, ymax=ODmean+sdOD),
              alpha = 0.3, colour = NA) +
  facet_grid (Strain ~.) +
  theme_bw()

#Filter to have just 1 strain in the plot
filtered_strain <- Exp27_v5 %>% filter (Strain == "5001")

filtered_strain <- Exp27_v5 %>% filter (Strain == "5002")

filtered_strain <- Exp27_v5 %>% filter (Strain == "5003")

filtered_strain <- Exp27_v5 %>% filter (Strain == "5004")

#Filter to compare the blanks (req_conc vs blank from the other stocks - 100 uM (ZnCl2 stock) vs TPEN blank) with a required concentration from a certain strain in the plot
filtered_data <- filtered_strain %>%
  filter(drug_conc %in% c("100 uM", "blank"))

filtered_data <- filtered_strain %>%
  filter(drug_conc %in% c("250 uM", "blank"))

filtered_data <- filtered_strain %>%
  filter(drug_conc %in% c("18 uM", "0.1HCl"))

filtered_data <- filtered_strain %>%
  filter(drug_conc %in% c("24 uM", "0.1HCl"))

filtered_data <- filtered_strain %>%
  filter(drug_conc %in% c("30 uM", "0.1HCl"))

filtered_data <- filtered_strain %>%
  filter(drug_conc %in% c("36 uM", "0.1HCl"))

filtered_data <- filtered_strain %>%
  filter(drug_conc %in% c("42 uM", "0.1HCl"))

filtered_data <- filtered_strain %>%
  filter(drug_conc %in% c("48 uM", "0.1HCl"))

filtered_data <- filtered_strain %>%
  filter(drug_conc %in% c("54 uM", "0.1HCl"))

filtered_data <- filtered_strain %>%
  filter(drug_conc %in% c("60 uM", "0.1HCl"))


ggplot(filtered_data, aes(x = Time_h, y = ODmean, color = Transfer_medium)) +
  geom_point() +
  geom_line() +
  geom_ribbon(aes(ymin = ODmean - sdOD, ymax = ODmean + sdOD, fill = Transfer_medium),
              alpha = 0.3, color = NA) +
  facet_wrap(~ Strain + drug_conc) +
  theme_bw()
