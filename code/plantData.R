
### Plant data exploraion
plant <- dat %>% select(shoot_g,root_g,pctN,pctC,numSpp,age.avg,age.var) %>%
  gather(key=measurement,value=response,shoot_g,root_g,pctN,pctC) %>%
  gather(key=treatment,predictor,numSpp,age.avg,age.var) %>%
  mutate(measurement=factor(measurement),treatment=factor(treatment))

levels(plant$measurement) <- 
  list("Shoot~biomass~(g)"="shoot_g",
       "Root~biomass~(g)"="root_g",
       "Leaf~'%'~N"="pctN",
       "Leaf~'%'~C"="pctC")

levels(plant$treatment) <- 
  list("Mean~log[10](age)"="age.avg",
       "Variance~log[10](age)"="age.var",
       "Species~pool~richness"="numSpp")

controlPlants <- seedlingData %>% filter(poolID=="C")

controlDF <- data.frame(measurement=rep(c("Shoot~biomass~(g)","Root~biomass~(g)","Leaf~'%'~N","Leaf~'%'~C"),each=3),
                        treatment=rep(c("Mean~log[10](age)","Variance~log[10](age)","Species~pool~richness"),4),
                        nofun=rep(c(median(controlPlants$shoot_g),median(controlPlants$root_g),
                                    median(controlPlants$pctN),median(controlPlants$pctC)),each=3),
                        letts=paste("(",letters[1:12],")"))

plant.plot <- ggplot(plant,aes(predictor,response)) +
  geom_point(alpha=0.5) +
  facet_grid(measurement~treatment,scales="free",switch="both", labeller = label_parsed) +
  geom_hline(data=controlDF,aes(yintercept=nofun),color="red") +
  geom_text(data=controlDF,aes(label=letts,x=-Inf,y=Inf),size=3,hjust=-0.25,vjust=2) +
  theme_tufte() + theme(panel.background = element_rect(colour = 'black'),strip.placement = "outside",axis.title=element_blank())

jpeg("figs/plantSupp.jpg",width=6.5,height=8,units="in",res=300);plant.plot;dev.off()
