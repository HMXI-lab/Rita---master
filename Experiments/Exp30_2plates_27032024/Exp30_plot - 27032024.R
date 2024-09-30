setwd("~/ITQB/Plate readings/Exp30_2plates_27032024")

library(tidyverse)

library(readxl)
Exp30 <- Exp30_2plates_27032024 <- read_excel("Exp30_2plates_27032024_final_final.xlsx")
View(Exp30) 

Exp30_v1 <- Exp30 %>%
  pivot_longer(names_to = "Condition", values_to = "OD_578", cols = -Time_h)
Exp30_v1
View(Exp30_v1)


Exp30_v2 <- separate(Exp30_v1, Condition, 
                     c("Community",
                       "Transfer_medium",
                       "Zn_condition", 
                       "Tech_Repl",
                       "Biol_Repl",
                       "TPEN_or_ZnCl2_conc",
                       "Well_number"), 
                     sep = "_", remove = FALSE,
                     convert = FALSE)
View(Exp30_v2)

Exp30_v2$OD_578 <- as.numeric (Exp30_v2$OD_578) #Because OD_578 values are character

Exp30_v2$TPEN_or_ZnCl2_conc <- factor(Exp30_v2$TPEN_or_ZnCl2_conc, levels = c("1.74",
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
                                                            "7.0"))

Exp30_v2 <- Exp30_v2 %>% mutate(TPEN_or_ZnCl2_conc = recode (TPEN_or_ZnCl2_conc,
                                                    "1.74" = "1.74uM",
                                                    "10" = "50uM",
                                                    "0" = "0uM",
                                                    "1.8" = "9uM",
                                                    "2.4" = "12uM",
                                                    "3.0" = "15uM",
                                                    "3.6" = "18uM",
                                                    "4.2" = "21uM",
                                                    "4.8" = "24uM",
                                                    "5.4" = "27uM",
                                                    "6.0" = "30uM",
                                                    "7.0" = "35uM"))

Exp30_v2$Community <- factor(Exp30_v2$Community, levels = c("C1",
                                                            "C2",
                                                            "C3",
                                                            "C4",
                                                            "C5"))
# 
# Exp30_v2 %>%
#   ggplot(aes(x= Time_h, y=OD_578, color = Tech_Repl)) +
#   geom_line() +
#   facet_grid (Community ~ TPEN_or_ZnCl2_conc)


Exp30_v2 %>%
  ggplot(aes(x= Time_h, y=OD_578, color = Tech_Repl)) +
  geom_line() +
  facet_grid (Community ~ TPEN_or_ZnCl2_conc)


Exp30_v3<- Exp30_v2 %>% group_by(Community, Transfer_medium, Zn_condition, Tech_Repl, Biol_Repl, Well_number) %>%
  mutate(blank = Time_h ==0)

View(Exp30_v3)

Exp30_v4 <- Exp30_v3 %>%
  group_by(Community, Transfer_medium, Zn_condition, Tech_Repl, Biol_Repl, Well_number) %>%
  mutate(OD = OD_578 - OD_578[blank]) %>%
  ungroup()


View(Exp30_v4)

#Obtain the mean and sd
Exp30_v5 <- Exp30_v4 %>% group_by(Time_h, Community, TPEN_or_ZnCl2_conc, Zn_condition) %>%
  summarise(ODmean = mean(OD), sdOD = sd(OD), n())

Exp30_v5$Community <- factor(Exp30_v5$Community, levels = c("C1",
                                                            "C2",
                                                            "C3",
                                                            "C4",
                                                            "C5"))


Exp30_v5 <- Exp30_v5 %>% mutate(TPEN_or_ZnCl2_conc =  recode (TPEN_or_ZnCl2_conc,
                                                              ".74uM" = "17.4 uM ZnCl2",
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
                                                              "35uM" = "35 uM TPEN"))





# stocks conc / 0.001(mM --> uM) / d.f (=2)--> conc TPEN (now in the table)

Exp30_v5  %>% ggplot(aes(x = Time_h,
                         y = ODmean,
                         color = TPEN_or_ZnCl2_conc,
                         fill = TPEN_or_ZnCl2_conc
)) +
  geom_point() +
  geom_line() +
  geom_ribbon(aes(ymin=ODmean-sdOD, ymax=ODmean+sdOD),
              alpha = 0.3, colour = NA) +
  facet_wrap (Community ~ .) +
  theme_bw()
