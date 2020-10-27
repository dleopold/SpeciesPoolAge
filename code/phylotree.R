# Plot phylogenetic tree

library(tidyverse)
library(magrittr)
library(ggtree)
library(treeio)
library(ape)
library(scico)
library(ggtext)

col.index <- log10(c(300,2100,20000,150000,4100000)) %>% scales::rescale(c(0.001,1)) %>% multiply_by(1000)
siteCols <- scico(1000, palette = 'hawaii', direction = -1, end=0.9)[col.index]

# read in issolate meta data
dat <- read.csv("data/cultureData.csv",as.is=T) %>% 
  dplyr::rename(label=isoID) %>%
  mutate(Taxonomy = ifelse(grepl("ales$",Taxonomy) | grepl("aceae$",Taxonomy) | grepl("mycetes$",Taxonomy),
                           Taxonomy,paste0("*",Taxonomy,"*")),
         label2=paste(label,accession,Taxonomy, sep=" | ")) %>%
  select(label,label2,site) 

# read in tree and merge with meta data
tree <- read.newick("data/LSU_phylo/result.tcs_weighted_phyml/result.tcs_weighted_phyml_tree.txt") %>% 
  root("KY109285.1") %>% 
  drop.tip("KY109285.1") %>%
  full_join(dat)

ggtree(tree) +
  #geom_tiplab(size=3,offset=.015,align=F) +
  geom_richtext(aes(x=x+0.025, label=label2),hjust=0, size=3,
                fill = NA, label.color = NA,label.padding = grid::unit(rep(0, 4), "pt"))+
  geom_tippoint(aes(x=x+0.008, fill=site),shape=21,size=2.5)+
  scale_fill_manual(values=siteCols,name="Chronosequence\nsite age (# of isolates)",na.translate=F,
                    labels=c("3.0&times;10^2 (4)","2.1&times;10^3 (13)","2.0&times;10^4 (13)","1.5&times;10^5 (13)","4.1&times;10^6 (11)")) +
  xlim(0,0.65) +
  geom_treescale(x=0) +
  theme(legend.text = element_markdown(),
        legend.title = element_text(vjust=0),
        legend.position = c(0.82,0.845),
        legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))

ggsave("output/figs/FigS1.jpg",units="cm",width=17,height=20)
