#supplemental figure showing variation in relative abundance of individual taxa at the end of the experiment.

library(tidyverse)
library(magrittr)
library(phyloseq)
library(Biostrings)

phy <- readRDS("output/rds/phy.rds")
dat <- sample_data(phy) %>% data.frame() %>%
  rownames_to_column("SampID") %>%
  mutate(poolID=as.numeric(poolID))

otu.tab <- phy %>%  rarefy_even_depth %>% 
  otu_table %>% data.frame %>%
  rownames_to_column("SampID") %>%
  pivot_longer(-SampID, names_to="isoID", values_to="Obs") %>%
  filter(Obs>0) %>%
  mutate(Pres=1)

#' keys to isolates used in expt
sites <- read.csv("data/sites.csv", as.is=T)
targets <- readDNAStringSet("data/AllIsolates_ITS2_200bp.fasta")
isoKey <- read.csv("data/cultureData.csv",as.is=T) %>% left_join(sites) %>% mutate(site_logAge=log10(site_age)) %>%
  .[match(names(targets),.$isoID),]

pools.mx <- read.csv("output/csv/pools.mx.csv")
colnames(pools.mx)[1] <- "poolID"

spp.dat <- pools.mx %>%
  pivot_longer(-poolID, names_to="isoID") %>% 
  filter(value==T) %>% select(-value) %>%
  right_join(dat %>% select(
    poolID,
    SampID,
    pool_rich,
    age_mean,
    age_var
  )) %>%
  left_join(otu.tab) %>%
  mutate(Pres=case_when(is.na(Pres) ~ 0,T ~ Pres),
         Obs=case_when(is.na(Obs) ~ 0,T ~ Obs)) %>%
  left_join(isoKey)
  
spp.dat$site_name[spp.dat$site_name=="Ola’a"] <- "Olaa"
spp.dat$site_name %<>% factor(levels=c("Thurston","Olaa","Laupahoehoe","Kohala","Kokee"))

SpSummary <- spp.dat %>% group_by(isoID) %>%
  summarise(n=n(),
            Pres=sum(Pres),
            prop=Pres/n,
            Obs=mean(Obs),
            ObsMod=log(Obs+0.1),
            Site=site_name[1],
            nPool=length(unique(poolID)))
SpSummary$Site[SpSummary$Site=="Ola’a"] <- "Olaa"
SpSummary$Site %<>% factor(levels=c("Thurston","Olaa","Laupahoehoe","Kohala","Kokee"))

Sp.order <- SpSummary$Obs
names(Sp.order) <- SpSummary$isoID
spp.dat$isoID %<>% factor(levels=names(sort(Sp.order)))

ggplot(spp.dat, aes(x=isoID, y=log(Obs+0.5))) +
  geom_boxplot(varwidth = T) +
  facet_grid(.~site_name, space = "free_x", 
             scales = "free_x",
             drop=T) +
  labs(y="log( relative abundance )") +
  ggthemes::theme_few() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle=55, hjust=1, vjust=1, size=8))
ggsave("output/figs/FigS3.jpg", width=7.5, height=3.5 )


fig.relabund <- ggplot(SpSummary, aes(x=Site, y=log(Obs+0.00005))) +
  geom_boxplot(varwidth = T)+
  ylab("log( average abundance )")+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())

kruskal.test(ObsMod~Site, data=SpSummary)

kruskal.test(prop~Site, data=SpSummary)



