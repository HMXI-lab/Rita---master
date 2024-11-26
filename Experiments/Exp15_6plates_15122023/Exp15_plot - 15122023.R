setwd("C:/Users/sergi/Downloads/ITQB/Plate readings/Exp15_6plates_15122023")

library(tidyverse)

library(readxl)
Exp15 <- Exp15_6plates_05122023 <- read_excel("Exp15_6plates_15122023_final_final.xlsx")
View(Exp15) 

Exp15_v1 <- Exp15 %>%
  pivot_longer(names_to = "Condition", values_to = "OD_578", cols = -Time_h)
Exp15_v1
View(Exp15_v1)


Exp15_v2 <- separate(Exp15_v1, Condition, 
                     c("Strain",
                       "Transfer_medium",
                       "Zn_condition", 
                       "Tech_Repl",
                       "Biol_Repl",
                       "TPEN_conc",
                       "Well_number"), 
                     sep = "_", remove = FALSE,
                     convert = FALSE)
View(Exp15_v2)

Exp15_v2$OD_578 <- as.numeric (Exp15_v2$OD_578) #Because OD_578 values are character

Exp15_v2$TPEN_conc <- factor(Exp15_v2$TPEN_conc, levels = c("ctrTPENcells",
                                                            "ctrTPENwZn",
                                                            "ctrTPEN",
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

Exp15_v2 <- Exp15_v2 %>% mutate(TPEN_conc = recode (TPEN_conc,
                                                    "ctrTPENcells" = "blank",
                                                    "ctrTPENwZn" = "blankZn",
                                                    "ctrTPEN" = "0 uM",
                                                    "3.6"= "18 uM",
                                                    "4.8"= "24 uM",
                                                    "6" = "30 uM",
                                                    "7.2" = "36 uM",
                                                    "8.4" = "42 uM",
                                                    "9.6" = "48 uM",
                                                    "10.8" = "54 uM",
                                                    "12" = " 60 uM",
                                                    "14" = "70 uM",
                                                    "16" = "80 uM",
                                                    "18" = "90 uM",
                                                    "20" = "100 uM",
                                                    "22" = "110 uM",
                                                    "24" = "120 uM",
                                                    "26" = "130 uM",
                                                    "28" = "140 uM"))

Exp15_v2$Strain <- factor(Exp15_v2$Strain, levels = c("5026",
                                                    "5038",
                                                    "5046",
                                                    "5048",
                                                    "5076",
                                                    "5078",
                                                    "12016",
                                                    "12097",
                                                    "12194",
                                                    "14071",
                                                    "14073",
                                                    "14076"))

# Exp15_v2 %>%
#   ggplot(aes(x= Time_h, y=OD_578, color = Tech_Repl)) +
#   geom_line() +
#   facet_grid (Strain ~ TPEN_conc)


Exp15_v2 %>% filter (Transfer_medium == "mGAM") %>%
  ggplot(aes(x= Time_h, y=OD_578, color = Tech_Repl)) +
  geom_line() +
  facet_grid (Strain ~ TPEN_conc)

Exp15_v2 %>% filter (Transfer_medium == "M8") %>%
  ggplot(aes(x= Time_h, y=OD_578, color = Tech_Repl)) +
  geom_line() +
  facet_grid (Strain ~ TPEN_conc)


Exp15_v3<- Exp15_v2 %>% group_by(Strain, Transfer_medium, Zn_condition, Tech_Repl, Biol_Repl, Well_number) %>%
  mutate(blank = Time_h ==0)

View(Exp15_v3)

Exp15_v4 <- Exp15_v3 %>%
  group_by(Strain, Transfer_medium, Zn_condition, Tech_Repl, Biol_Repl, Well_number) %>%
  mutate(OD = OD_578 - OD_578[blank]) %>%
  ungroup()


View(Exp15_v4)

#Obtain the mean and sd
Exp15_v5 <- Exp15_v4 %>% group_by(Time_h, Strain, Transfer_medium, TPEN_conc) %>%
  summarise(ODmean = mean(OD), sdOD = sd(OD), n())

Exp15_v5$Strain <- factor(Exp15_v5$Strain, levels = c("5026",
                                                      "5038",
                                                      "5046",
                                                      "5048",
                                                      "5076",
                                                      "5078",
                                                      "12016",
                                                      "12097",
                                                      "12194",
                                                      "14071",
                                                      "14073",
                                                      "14076"))


Exp15_v5 <- Exp15_v5 %>% mutate(TPEN_conc = recode (TPEN_conc,
                                                    "ctrTPENcells" = "blank",
                                                    "ctrTPENwZn" = "blankZn",
                                                    "ctrTPEN" = "0 uM",
                                                    "3.6"= "18 uM",
                                                    "4.8"= "24 uM",
                                                    "6" = "30 uM",
                                                    "7.2" = "36 uM",
                                                    "8.4" = "42 uM",
                                                    "9.6" = "48 uM",
                                                    "10.8" = "54 uM",
                                                    "12" = " 60 uM",
                                                    "14" = "70 uM",
                                                    "16" = "80 uM",
                                                    "18" = "90 uM",
                                                    "20" = "100 uM",
                                                    "22" = "110 uM",
                                                    "24" = "120 uM",
                                                    "26" = "130 uM",
                                                    "28" = "140 uM"))


# stocks conc / 0.001(mM --> uM) / d.f (=2)--> conc TPEN (now in the table)


Exp15_v5  %>% filter (Transfer_medium == "M8") %>%
  ggplot(aes(x = Time_h,
             y = ODmean,
             color = Transfer_medium
  )) +
  geom_point() + 
  geom_line() +
  geom_ribbon(aes(ymin=ODmean-sdOD, ymax=ODmean+sdOD, fill = Transfer_medium), 
              alpha = 0.3, colour = NA) + 
  facet_grid (Strain ~ TPEN_conc) +
  theme_bw()


Exp15_v5  %>% filter (Transfer_medium == "mGAM") %>%
  ggplot(aes(x = Time_h,
             y = ODmean,
             color = Transfer_medium
  )) +
  geom_point() + 
  geom_line() +
  geom_ribbon(aes(ymin=ODmean-sdOD, ymax=ODmean+sdOD, fill = Transfer_medium), 
              alpha = 0.3, colour = NA) + 
  facet_grid (Strain ~ TPEN_conc) +
  theme_bw()
