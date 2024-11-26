setwd("C:/Users/sergi/Downloads/ITQB/Plate readings/Exp2_1plate_12102023")

library(tidyverse)

library(readxl)
Exp2 <- Exp2_1plate_12102023_final <- read_excel("Exp2_1plate_12102023_final_final.xlsx")
View(Exp2) 

Exp2_v1 <- Exp2 %>%
  pivot_longer(names_to = "Condition", values_to = "OD_578", cols = -Time_h)
Exp2_v1
View(Exp2_v1)


Exp2_v2 <- separate(Exp2_v1, Condition, 
                    c("Strain",
                      "Transfer_medium",
                      "Zn_condition", 
                      "Tech_Repl",
                      "Biol_Repl",
                      "Well_number"), 
                    sep = "_", remove = FALSE,
                    convert = FALSE)
View(Exp2_v2)

Exp2_v2$OD_578 <- as.numeric (Exp2_v2$OD_578) #Because OD_578 values are  character

Exp2_v2 %>%
  ggplot(aes(x= Time_h, y=OD_578, color = Transfer_medium)) +
  geom_line() +
  facet_wrap (Strain  ~ Tech_Repl)


Exp2_v3<- Exp2_v2 %>% group_by(Strain, Transfer_medium, Zn_condition, Tech_Repl, Biol_Repl, Well_number) %>%
  mutate(blank = Time_h ==0)

View(Exp2_v3)

Exp2_v4 <- Exp2_v3 %>%
  group_by(Strain, Transfer_medium, Zn_condition, Tech_Repl, Biol_Repl, Well_number) %>%
  mutate(OD = OD_578 - OD_578[blank]) %>%
  ungroup()


View(Exp2_v4)

#Obtain the mean and sd
Exp2_v5 <- Exp2_v4 %>% group_by(Time_h, Strain, Zn_condition, Transfer_medium) %>%
  summarise(ODmean = mean(OD), sdOD = sd(OD), n())

Exp2_v5  %>% 
  ggplot(mapping=aes(x = Time_h,
                     y = ODmean,
                     color = Transfer_medium)) +
  geom_point() + 
  geom_line() +
  geom_ribbon(aes(ymin=ODmean-sdOD, ymax=ODmean+sdOD, fill = Transfer_medium), 
              alpha = 0.3, colour = NA) + 
  facet_grid (Zn_condition ~ Strain) +
  theme_bw()