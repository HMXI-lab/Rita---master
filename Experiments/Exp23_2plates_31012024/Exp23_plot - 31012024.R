setwd("C:/Users/sergi/Downloads/ITQB/Plate readings/Exp23_2plates_31012024")

library(tidyverse)
library(readxl)

Exp23 <- Exp23_4plates_24012024 <- read_excel("Exp23_2plates_31012024_final_final.xlsx")
# View(Exp23)

Exp23_v1 <- Exp23 %>%
  pivot_longer(names_to = "Condition", values_to = "OD_578", cols = -Time_h)
# Exp23_v1
# View(Exp23_v1)


Exp23_v2 <- separate(Exp23_v1, Condition, 
                     c("Strain",
                       "Transfer_medium",
                       "Zn_condition", 
                       "Tech_Repl",
                       "Biol_Repl",
                       "TPEN_conc",
                       "Well_number"), 
                     sep = "_", remove = FALSE,
                     convert = FALSE)
# View(Exp23_v2)

Exp23_v2$OD_578 <- as.numeric (Exp23_v2$OD_578) #Because OD_578 values are character

Exp23_v2$Strain <- factor(Exp23_v2$Strain, levels = c("5006",
                                                      "5007",
                                                      "23006",
                                                      "HM04"))

Exp23_v2 <- Exp23_v2 %>% mutate(TPEN_conc = recode (TPEN_conc,
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
                                                    "12" = "60 uM",
                                                    "14" = "70 uM",
                                                    "16" = "80 uM",
                                                    "18" = "90 uM",
                                                    "20" = "100 uM",
                                                    "22" = "110 uM",
                                                    "24" = "120 uM",
                                                    "26" = "130 uM",
                                                    "28" = "140 uM"))

Exp23_v2$TPEN_conc <- factor(Exp23_v2$TPEN_conc, levels = c("blank",
                                                      "blankZn",
                                                      "0 uM",
                                                      "18 uM",
                                                      "24 uM",
                                                      "30 uM",
                                                      "36 uM",
                                                      "42 uM",
                                                      "48 uM",
                                                      "54 uM",
                                                      "60 uM",
                                                      "70 uM",
                                                      "80 uM",
                                                      "90 uM",
                                                      "100 uM",
                                                      "110 uM",
                                                      "120 uM",
                                                      "130 uM",
                                                      "140 uM"))



# Exp23_v2 %>%
#   ggplot(aes(x= Time_h, y=OD_578, color = Tech_Repl)) +
#   geom_line() +
#   facet_grid (Strain ~ TPEN_conc)


Exp23_v2 %>% filter (Transfer_medium == "mGAM") %>%
  ggplot(aes(x= Time_h, y=OD_578, color = Tech_Repl)) +
  geom_line() +
  facet_grid (Strain ~ TPEN_conc)

Exp23_v2 %>% filter (Transfer_medium == "M8") %>%
  ggplot(aes(x= Time_h, y=OD_578, color = Tech_Repl)) +
  geom_line() +
  facet_grid (Strain ~ TPEN_conc)

Exp23_v3<- Exp23_v2 %>% group_by(Strain, Transfer_medium, Zn_condition, Tech_Repl, Biol_Repl, Well_number) %>%
  mutate(blank = Time_h ==0)

# View(Exp23_v3)

Exp23_v4 <- Exp23_v3 %>%
  group_by(Strain, Transfer_medium, Zn_condition, Tech_Repl, Biol_Repl, Well_number) %>%
  mutate(OD = OD_578 - OD_578[blank]) %>%
  ungroup()


# View(Exp23_v4)

#Obtain the mean and sd
Exp23_v5 <- Exp23_v4 %>% group_by(Time_h, Strain, Transfer_medium, TPEN_conc) %>%
  summarise(ODmean = mean(OD), sdOD = sd(OD), n())


Exp23_v5 <- Exp23_v5 %>% mutate(TPEN_conc = recode (TPEN_conc,
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
                                                    "12" = "60 uM",
                                                    "14" = "70 uM",
                                                    "16" = "80 uM",
                                                    "18" = "90 uM",
                                                    "20" = "100 uM",
                                                    "22" = "110 uM",
                                                    "24" = "120 uM",
                                                    "26" = "130 uM",
                                                    "28" = "140 uM"))


# stocks conc / 0.001(mM --> uM) / d.f (=2)--> conc TPEN (now in the table)


Exp23_v5  %>% filter (Transfer_medium == "mGAM") %>%
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

Exp23_v5  %>% filter (Transfer_medium == "M8") %>%
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
