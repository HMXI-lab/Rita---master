setwd("~/ITQB/Plate readings/Exp8_6plates_16112023")

library(tidyverse)

library(readxl)
Exp8b <- Exp8b_6plates_16112023_final <- read_excel("Exp8b_6plates_16112023_final_final.xlsx")
View(Exp8b) 

Exp8b_v1 <- Exp8b %>%
  pivot_longer(names_to = "Condition", values_to = "OD_578", cols = -Time_h)
Exp8b_v1
View(Exp8b_v1)


Exp8b_v2 <- separate(Exp8b_v1, Condition, 
                     c("Strain",
                       "Transfer_medium",
                       "Zn_condition", 
                       "Tech_Repl",
                       "Biol_Repl",
                       "TPEN_conc",
                       "Well_number"), 
                     sep = "_", remove = FALSE,
                     convert = FALSE)
View(Exp8b_v2)

Exp8b_v2$OD_578 <- as.numeric (Exp8b_v2$OD_578) #Because OD_578 values are character

Exp8b_v2$TPEN_conc <- factor(Exp8b_v2$TPEN_conc, levels = c("ctrTPENcells",
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
                                                            "16"))

Exp8b_v2 <- Exp8b_v2 %>% mutate(TPEN_conc = recode (TPEN_conc,
                                                    "ctrTPENcells" = "blank" ,
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
                                                    "16" = "80 uM"))

Exp8b_v2$Strain <- factor(Exp8b_v2$Strain, levels = c("5026",
                                                      "5037",
                                                      "5038",
                                                      "5071",
                                                      "5076",
                                                      "5078",
                                                      "HM04"
))

Exp8b_v2 %>%
  ggplot(aes(x= Time_h, y=OD_578, color = Tech_Repl)) +
  geom_line() +
  facet_grid (Strain ~ TPEN_conc)

Exp8b_v3<- Exp8b_v2 %>% group_by(Strain, Transfer_medium, Zn_condition, Tech_Repl, Biol_Repl, Well_number) %>%
  mutate(blank = Time_h ==0)

View(Exp8b_v3)

Exp8b_v4 <- Exp8b_v3 %>%
  group_by(Strain, Transfer_medium, Zn_condition, Tech_Repl, Biol_Repl, Well_number) %>%
  mutate(OD = OD_578 - OD_578[blank]) %>%
  ungroup()


View(Exp8b_v4)

#Obtain the mean and sd
Exp8b_v5 <- Exp8b_v4 %>% group_by(Time_h, Strain, Transfer_medium, TPEN_conc) %>%
  summarise(ODmean = mean(OD), sdOD = sd(OD), n())

Exp8b_v5$Strain <- factor(Exp8b_v5$Strain, levels = c("5026",
                                                      "5037",
                                                      "5038",
                                                      "5071",
                                                      "5076",
                                                      "5078",
                                                      "HM04"                     
))


Exp8b_v5 <- Exp8b_v5 %>% mutate(TPEN_conc = recode (TPEN_conc,
                                                    "ctrTPEN" = "blank" ,
                                                    "ctrTPENcells" = "0 uM",
                                                    "3.6"= "18 uM",
                                                    "4.8"= "24 uM",
                                                    "6" = "30 uM",
                                                    "7.2" = "36 uM",
                                                    "8.4" = "42 uM",
                                                    "9.6" = "48 uM",
                                                    "10.8" = "54 uM",
                                                    "12" = " 60 uM",
                                                    "14" = "70 uM",
                                                    "16" = "80 uM"))


# stocks conc / 0.001(mM --> uM) / d.f (=2)--> conc TPEN (now in the table)


Exp8b_v5  %>% 
  ggplot(aes(x = Time_h,
             y = ODmean,
             color = Transfer_medium
  )) +
  geom_point() + 
  geom_line() +
  geom_ribbon(aes(ymin=ODmean-sdOD, ymax=ODmean+sdOD, fill = Transfer_medium), 
              alpha = 0.3, colour = NA) + 
  facet_grid (Strain ~ TPEN_conc) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 6)  # Decrease the size of the numbers on the axes
  )

