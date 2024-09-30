setwd("~/ITQB/Plate readings/Exp11_3plates_24112023")

library(tidyverse)

library(readxl)
Exp11 <- Exp11_3plates_24112023 <- read_excel("Exp11_3plates_24112023_final_final.xlsx")
View(Exp11) 

Exp11_v1 <- Exp11 %>%
  pivot_longer(names_to = "Condition", values_to = "OD_578", cols = -Time_h)
Exp11_v1
View(Exp11_v1)


Exp11_v2 <- separate(Exp11_v1, Condition, 
                     c("Strain",
                       "Transfer_medium",
                       "Zn_condition", 
                       "Tech_Repl",
                       "Biol_Repl",
                       "TPEN_conc",
                       "Well_number"), 
                     sep = "_", remove = FALSE,
                     convert = FALSE)
View(Exp11_v2)

Exp11_v2$OD_578 <- as.numeric (Exp11_v2$OD_578) #Because OD_578 values are character

Exp11_v2$TPEN_conc <- factor(Exp11_v2$TPEN_conc, levels = c("ctrTPENcells",
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

Exp11_v2 <- Exp11_v2 %>% mutate(TPEN_conc = recode (TPEN_conc,
                                                    "ctrTPENcells" = "blank",
                                                    "ctrTPEN" = "0 uM",
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

Exp11_v2$Strain <- factor(Exp11_v2$Strain, levels = c("5019",
                                                      "5026",
                                                      "5046",
                                                      "5048",
                                                      "5071",
                                                      "5076",
                                                      "5078",
                                                      "HM04"))

Exp11_v2 %>%
  ggplot(aes(x= Time_h, y=OD_578, color = Tech_Repl)) +
  geom_line() +
  facet_grid (Strain ~ TPEN_conc)

Exp11_v3<- Exp11_v2 %>% group_by(Strain, Transfer_medium, Zn_condition, Tech_Repl, Biol_Repl, Well_number) %>%
  mutate(blank = Time_h ==0)

View(Exp11_v3)

Exp11_v4 <- Exp11_v3 %>%
  group_by(Strain, Transfer_medium, Zn_condition, Tech_Repl, Biol_Repl, Well_number) %>%
  mutate(OD = OD_578 - OD_578[blank]) %>%
  ungroup()


View(Exp11_v4)

#Obtain the mean and sd
Exp11_v5 <- Exp11_v4 %>% group_by(Time_h, Strain, Transfer_medium, TPEN_conc) %>%
  summarise(ODmean = mean(OD), sdOD = sd(OD), n())

Exp11_v5$Strain <- factor(Exp11_v5$Strain, levels = c("5019",
                                                     "5026",
                                                     "5046",
                                                     "5048",
                                                     "5071",
                                                     "5076",
                                                     "5078",
                                                     "HM04"               
))


Exp11_v5 <- Exp11_v5 %>% mutate(TPEN_conc = recode (TPEN_conc,
                                                    "ctrTPENcells" = "blank",
                                                    "ctrTPEN" = "0 uM",
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


Exp11_v5  %>% 
  ggplot(aes(x = Time_h,
             y = ODmean,
             color = Transfer_medium,
             fill = "darkturquoise"
  )) +
  geom_point(color = "darkturquoise") + 
  geom_line(color = "darkturquoise") +
  geom_ribbon(aes(ymin=ODmean-sdOD, ymax=ODmean+sdOD, fill = Transfer_medium), 
              alpha = 0.3, colour = NA) + 
  facet_grid (Strain ~ TPEN_conc) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 6)  # Decrease the size of the numbers on the axes
  )
