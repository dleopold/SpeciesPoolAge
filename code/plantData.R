# Exploration of the realtionship between fungal species pool characteristics and host plant measurements

library(tidyverse)
library(magrittr)
library(ggtext)

plant.dat <- read.csv("data/PlantData.csv",as.is = T) %>% 
  drop_na %>% filter(poolID!=21) %>% select(-shoot_dry) %>%
  pivot_longer(4:7,values_to = "response", names_to="measurement") %>%
  left_join(read.csv("output/csv/pools.dat.csv",as.is=T) %>% select(1:4) %>% 
              mutate(poolID=as.character(poolID))) %>%
  pivot_longer(6:8, values_to="predictor", names_to="treatment")

plant.dat %<>% mutate(treatment=recode(treatment,
                                       pool_rich="Species pool richness",
                                       age_mean="Mean log_10(age)",
                                       age_var="Variance log_10(age)"))

plant.dat %<>% mutate(measurement=recode(measurement,
                                       pctN="Leaf % N",
                                       pctC="Leaf % C",
                                       shoot_g="Shoot biomass (g)",
                                       root_g="Root biomass (g)",
                                       ))

control.dat <- plant.dat %>% filter(poolID=="C") %>%
  group_by(measurement,treatment) %>%
  summarize(response=median(response))

plant.dat %<>% filter(poolID!="C")

(plant.plot <- ggplot(plant.dat,aes(x=predictor,y=response)) +
  geom_point(alpha=0.5) +
  facet_grid(measurement~treatment,scales="free",switch="both") +
  geom_hline(data=control.dat,aes(yintercept=response),color="red") +
  ggthemes::theme_few() + 
  theme(panel.background = element_rect(colour = 'black'),
        strip.placement = "outside",
        strip.text = element_markdown(),
        axis.title=element_blank()))

ggsave("output/figs/FigS2.jpg",width=18,height=22, units="cm")

