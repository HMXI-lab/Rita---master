
setwd("C:/Users/sergi/Downloads/ITQB/Plate readings/Exp1_1plate_27092023")

library(tidyverse)

library(readxl)
Exp1 <- Exp1_1plate_27092023_final <- read_excel("Exp1_1plate_27092023_final_final.xlsx")
View(Exp1_1plate_27092023_final_R) 

Exp1_v1 <- Exp1 %>%
  pivot_longer(names_to = "Condition", values_to = "OD_578", cols = -Time_h)
Exp1_v1


Exp1_v2 <- separate(Exp1_v1, Condition, 
                    c("Strain",
                      "Zn_condition", 
                      "Tech_Repl",
                      "Biol_Repl",
                      "Well_number"), 
                    sep = "_", remove = FALSE,
                    convert = FALSE)
Exp1_v2
view (Exp1_v2)

Exp1_v2$OD_578 <- as.numeric (Exp1_v2$OD_578)

Exp1_v2 %>%
  ggplot(aes(x= Time_h, y=OD_578, color = Zn_condition)) +
  geom_line() +
  facet_wrap (Strain  ~ Tech_Repl)


Exp1_v3<- Exp1_v2 %>% group_by(Strain, Zn_condition, Tech_Repl, Biol_Repl, Well_number) %>%
  mutate(blank = Time_h ==0)

View(Exp1_v3)

Exp1_v4 <- Exp1_v3 %>%
  group_by(Strain, Zn_condition, Tech_Repl, Biol_Repl, Well_number) %>%
  mutate(OD = OD_578 - OD_578[blank]) %>%
  ungroup()


View(Exp1_v4)

#Obtain the mean and sd
Exp1_v5 <- Exp1_v4 %>% group_by(Time_h, Strain, Zn_condition) %>%
  summarise(ODmean = mean(OD), sdOD = sd(OD), n())

Exp1_v5  %>% 
  ggplot(mapping=aes(x = Time_h,
                     y = ODmean,
                     color = Zn_condition
  )) +
  geom_point() + 
  geom_line() +
  geom_ribbon(aes(ymin=ODmean-sdOD, ymax=ODmean+sdOD, fill = Zn_condition), 
              alpha = 0.3, colour = NA) + 
  facet_grid (. ~Strain) +
  theme_bw()
