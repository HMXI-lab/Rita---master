setwd("~/ITQB/Plate readings/Exp3_6plates_18102023")

library(tidyverse)

library(readxl)
Exp3 <- Exp3_6plates_18102023_final <- read_excel("Exp3_6plates_18102023_final_final.xlsx")
View(Exp3) 

Exp3_v1 <- Exp3 %>%
  pivot_longer(names_to = "Condition", values_to = "OD_578", cols = -Time_h)
Exp3_v1
View(Exp3_v1)


Exp3_v2 <- separate(Exp3_v1, Condition, 
                    c("Strain",
                      "Transfer_medium",
                      "Zn_condition", 
                      "Tech_Repl",
                      "Biol_Repl",
                      "TPEN_conc",
                      "Well_number"), 
                    sep = "_", remove = FALSE,
                    convert = FALSE)
View(Exp3_v2)

Exp3_v2$OD_578 <- as.numeric (Exp3_v2$OD_578) #Because OD_578 values are  character

Exp3_v2 %>%
  ggplot(aes(x= Time_h, y=OD_578, color = Transfer_medium)) +
  geom_line() +
  facet_wrap (Strain ~ TPEN_conc)



Exp3_v3<- Exp3_v2 %>% group_by(Strain, Transfer_medium, Zn_condition, Tech_Repl, Biol_Repl, Well_number) %>%
  mutate(blank = Time_h ==0)

View(Exp3_v3)

Exp3_v4 <- Exp3_v3 %>%
  group_by(Strain, Transfer_medium, Zn_condition, Tech_Repl, Biol_Repl, Well_number) %>%
  mutate(OD = OD_578 - OD_578[blank]) %>%
  ungroup()


View(Exp3_v4)

#Obtain the mean and sd
Exp3_v5 <- Exp3_v4 %>% group_by(Time_h, Strain, Transfer_medium, TPEN_conc) %>%
  summarise(ODmean = mean(OD), sdOD = sd(OD), n())

Exp3_v5$TPEN_conc <- factor(Exp3_v5$TPEN_conc, levels = c("ctrTPEN", 
                                                         "ctrTPENcells",
                                                         "1.2",
                                                         "2.4",
                                                         "3.6",
                                                         "4.8",
                                                         "6",
                                                         "7.2",
                                                         "8.4",
                                                         "9.6",
                                                         "10.8",
                                                         "12"
                                                           ))



Exp3_v5 <- Exp3_v5 %>% mutate(TPEN_conc = recode (TPEN_conc,
                                                  "ctrTPEN" = "blank" ,
                                                  "ctrTPENcells" = "0 uM",
                                                  "1.2" = "6 uM",
                                                  "2.4" = "12 uM",
                                                  "3.6"= "18 uM",
                                                  "4.8"= "24 uM",
                                                  "6" = "30 uM",
                                                  "7.2" = "36 uM",
                                                  "8.4" = "42 uM",
                                                  "9.6" = "48 uM",
                                                  "10.8" = "54 uM",
                                                  "12" = " 60 uM"))
                                 

Exp3_v5  %>% filter (TPEN_conc != "blank") %>%
  ggplot(mapping=aes(x = Time_h,
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

