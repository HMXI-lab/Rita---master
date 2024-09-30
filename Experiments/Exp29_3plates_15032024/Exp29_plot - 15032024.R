setwd("~/ITQB/Plate readings/Exp29_3plates_15032024")

library(tidyverse)

library(readxl)
Exp29 <- Exp29_3plates_15032024 <- read_excel("Exp29_3plates_15032024_final_final.xlsx")
View(Exp29) 

Exp29_v1 <- Exp29 %>%
  pivot_longer(names_to = "Condition", values_to = "OD_578", cols = -Time_h)
Exp29_v1
View(Exp29_v1)


Exp29_v2 <- separate(Exp29_v1, Condition, 
                     c("Strain",
                       "Transfer_medium",
                       "Zn_condition", 
                       "Tech_Repl",
                       "Biol_Repl",
                       "TPEN_or_ZnCl2_conc",
                       "Well_number"), 
                     sep = "_", remove = FALSE,
                     convert = FALSE)
View(Exp29_v2)

Exp29_v2$OD_578 <- as.numeric (Exp29_v2$OD_578) #Because OD_578 values are character

Exp29_v2$TPEN_or_ZnCl2_conc <- factor(Exp29_v2$TPEN_or_ZnCl2_conc, levels = c("ctrZnCl2",
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

Exp29_v2 <- Exp29_v2 %>% mutate(TPEN_or_ZnCl2_conc = recode (TPEN_or_ZnCl2_conc,
                                                    "ctrZnCl2" = "0.1HCl",
                                                    "10" = "100uM",
                                                    "25" = "250uM",
                                                    "ctrTPEN" = "DMSO",
                                                    "1.8" = "18uM",
                                                    "2.4" = "24uM",
                                                    "3.0" = "30uM",
                                                    "3.6" = "36uM",
                                                    "4.2" = "42uM",
                                                    "4.8" = "48uM",
                                                    "5.4" = "54uM",
                                                    "6.0" = "60uM"))

Exp29_v2$Strain <- factor(Exp29_v2$Strain, levels = c("5001",
                                                    "5002",
                                                    "5003",
                                                    "5004"))

# Exp29_v2 %>%
#   ggplot(aes(x= Time_h, y=OD_578, color = Tech_Repl)) +
#   geom_line() +
#   facet_grid (Strain ~ TPEN_or_ZnCl2_conc)

#I don't need this plot because at this experiment I only had M8 medium.
# Exp29_v2 %>% filter (Transfer_medium == "mGAM") %>%
#   ggplot(aes(x= Time_h, y=OD_578, color = Tech_Repl)) +
#   geom_line() +
#   facet_grid (Strain ~ TPEN_or_ZnCl2_conc)

Exp29_v2 %>% filter (Transfer_medium == "M8") %>%
  ggplot(aes(x= Time_h, y=OD_578, color = Tech_Repl)) +
  geom_line() +
  facet_grid (Strain ~ TPEN_or_ZnCl2_conc + Zn_condition)


Exp29_v3<- Exp29_v2 %>% group_by(Strain, Transfer_medium, Zn_condition, Tech_Repl, Biol_Repl, Well_number) %>%
  mutate(blank = Time_h ==0)

View(Exp29_v3)

Exp29_v4 <- Exp29_v3 %>%
  group_by(Strain, Transfer_medium, Zn_condition, Tech_Repl, Biol_Repl, Well_number) %>%
  mutate(OD = OD_578 - OD_578[blank]) %>%
  ungroup()


View(Exp29_v4)

#Obtain the mean and sd
Exp29_v5 <- Exp29_v4 %>% group_by(Time_h, Strain, TPEN_or_ZnCl2_conc, Zn_condition) %>%
  summarise(ODmean = mean(OD), sdOD = sd(OD), n())

Exp29_v5$Strain <- factor(Exp29_v5$Strain, levels = c("5001",
                                                      "5002",
                                                      "5003",
                                                      "5004"))


Exp29_v5 <- Exp29_v5 %>% mutate(TPEN_or_ZnCl2_conc =  recode (TPEN_or_ZnCl2_conc,
                                                              "100uM" = "100 uM ZnCl2",
                                                              "250uM" = "250 uM ZnCl2",
                                                              "18uM" = "18 uM TPEN",
                                                              "24uM" = "24 uM TPEN",
                                                              "30uM" = "30 uM TPEN",
                                                              "36uM" = "36 uM TPEN",
                                                              "42uM" = "42 uM TPEN",
                                                              "48uM" = "48 uM TPEN",
                                                              "54uM" = "54 uM TPEN",
                                                              "60uM" = "60 uM TPEN"))





# stocks conc / 0.001(mM --> uM) / d.f (=2)--> conc TPEN (now in the table)

Exp29_v5  %>% ggplot(aes(x = Time_h,
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

# #Filter to have just 1 strain in the plot
# filtered_strain <- Exp29_v5 %>% filter (Strain == "5001")
# 
# filtered_strain <- Exp29_v5 %>% filter (Strain == "5002")
# 
# filtered_strain <- Exp29_v5 %>% filter (Strain == "5003")
# 
# filtered_strain <- Exp29_v5 %>% filter (Strain == "5004")

##Filter to compare the blanks (req_conc vs blank from the other stocks - 100 uM (ZnCl2 stock) vs TPEN blank) with a required concentration from a certain strain in the plot
# filtered_data <- filtered_strain %>%
#   filter(TPEN_or_ZnCl2_conc %in% c("100 uM", "blank"))
# 
# filtered_data <- filtered_strain %>%
#   filter(TPEN_or_ZnCl2_conc %in% c("250 uM", "blank"))
# 
# filtered_data <- filtered_strain %>%
#   filter(TPEN_or_ZnCl2_conc %in% c("18 uM", "0.1HCl"))
# 
# filtered_data <- filtered_strain %>%
#   filter(TPEN_or_ZnCl2_conc %in% c("24 uM", "0.1HCl"))
# 
# filtered_data <- filtered_strain %>%
#   filter(TPEN_or_ZnCl2_conc %in% c("30 uM", "0.1HCl"))
# 
# filtered_data <- filtered_strain %>%
#   filter(TPEN_or_ZnCl2_conc %in% c("36 uM", "0.1HCl"))
# 
# filtered_data <- filtered_strain %>%
#   filter(TPEN_or_ZnCl2_conc %in% c("42 uM", "0.1HCl"))
# 
# filtered_data <- filtered_strain %>%
#   filter(TPEN_or_ZnCl2_conc %in% c("48 uM", "0.1HCl"))
# 
# filtered_data <- filtered_strain %>%
#   filter(TPEN_or_ZnCl2_conc %in% c("54 uM", "0.1HCl"))
# 
# filtered_data <- filtered_strain %>%
#   filter(TPEN_or_ZnCl2_conc %in% c("60 uM", "0.1HCl"))
# 
# 
# ggplot(filtered_data, aes(x = Time_h, y = ODmean, color = Transfer_medium)) +
#   geom_point() +
#   geom_line() +
#   geom_ribbon(aes(ymin = ODmean - sdOD, ymax = ODmean + sdOD, fill = Transfer_medium),
#               alpha = 0.3, color = NA) +
#   facet_wrap(~ Strain + TPEN_or_ZnCl2_conc) +
#   theme_bw()
