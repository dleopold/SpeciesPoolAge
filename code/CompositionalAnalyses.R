#' ---
#' title: Analyses of functional / phylogenitic diversity
#' author: Devin R Leopold
#' date: May 23, 2019
#' output:
#'    html_document:
#'      toc: true
#'      toc_float: false
#'      self_contained: true
#'      highlight: zenburn
#' ---

#+ warning=F, message=F
#Load packages
library(tidyverse)
library(magrittr)
library(phyloseq)
library(picante)
library(ape)
library(nlme)
library(foreach)
library(doMC)
#use all resources for parallel processing
registerDoMC(cores=parallel::detectCores())

#' # Prepare phylogenetic data
#' Import phylogenetic tree.  
#' Root on outgroup and drop outgroup for analyses
tree <- read.tree("data/LSU_phylo/result.tcs_weighted_phyml/result.tcs_weighted_phyml_tree.txt") %>%
  root("KY109285.1") %>% drop.tip("KY109285.1") 
tree$tip.label %<>% paste0("ISO.",.)

#' Make scaled phylogenetic distance matrix (sqrt transform following Letten & Cornwell 2015)
Pdist <- cophenetic(tree) %>% .[sort(rownames(.)),sort(rownames(.))] %>% 
  sqrt %>% as.dist(upper=T,diag=F) %>% 
  subtract(.,min(.)) %>% divide_by(.,max(.)) %>% as.matrix

#' # Prepare functional data
#' Import phenotypic microarray data
biolog <- read.csv("output/csv/biolog.csv",row.names=1) %>% .[sort(rownames(.)),]

#' Collapse correlated phenotype measurements with PCA
biolog_pca <- biolog %>% prcomp(scale=T,rank=10) %>% .$x

#' Make scaled functional distance matrix
Fdist <- biolog_pca %>% dist(upper=T) %>%
  subtract(.,min(.)) %>% divide_by(.,max(.)) %>% as.matrix

#' # Combined functional/phylogenetic distance
#' Function based on Cadotte et al 2013
FPdist_fun <- function(Pdist,Fdist,a,p=2){
  Pfoo <- a * (Pdist^p) 
  Ffoo <- (1-a) * (Fdist^p)
  (Pfoo + Ffoo)^(1/p)
}
FPdist <- FPdist_fun(Pdist,Fdist,0.5) #equal weighting

#' # Calculate species pool functional / phylogenetic diversity
#' Load matrix of species pool compositions and pool characteristics
pools <- read.csv("output/csv/pools.mx.csv", as.is=T, header=T, row.names = 1) %>% 
  .[,order(colnames(.))]
pools.dat <- read.csv("output/csv/pools.dat.csv", as.is=T, header=T) 

#' Species pool 21 is an extreme outlier in terms of FPdist (Bonferroni p = 5.0199e-05 from car::outlierTest) so we will remove it. Howver, this does not affect significance or magnitude of overall trends.
pools %<>% .[-21, ]
pools.dat %<>% filter(poolID!="21")

#' Add mean pairwise distance metrics
pools.dat$Fdist <- mpd(pools,Fdist)
pools.dat$Pdist <- mpd(pools,Pdist)
pools.dat$FPdist <- mpd(pools,FPdist)

#' ## Test species pool functional / phylogenetic diversty vs. mean age or age variance
lm(FPdist~pool_rich+age_var+age_mean,data = pools.dat) %>% summary
#' Partial R-squared
lm(FPdist~pool_rich+age_var+age_mean,data = pools.dat) %>% rsq::rsq.partial()

#' ## Plot species pool age vs functional / phylogenetic diversity
#+  fig.align="center", fig.asp=8/11, out.width="75%"
pools.dat %>% 
  ggplot(aes(age_mean,FPdist)) +
  stat_smooth(method="lm",colour="black",size=0.5) +
  geom_point(aes(size=pool_rich),shape=19,alpha=0.4) +
  scale_size(range = c(0.5,5.5))+
  labs(x=expression("Mean "*log[10]*"(age)"), 
       y="Mean-pairwise functional /\n phylogenetic distance") +
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_text(size=11),
        axis.text = element_text(size=8),
        axis.line = element_line(size=0.3),
        axis.ticks = element_line(size=0.12),
        axis.ticks.length = unit(.05, "cm"))
ggsave("output/figs/Fig2.pdf",height=8, width=11, units="cm") 

#' ## Plot species pool age or variance vs functional / phylogenetic diversity
pan.lett <- data.frame(metric=c("Mean age","Variance age"),lab=c("(a)","(b)"))
#+ fig.align="center", fig.asp=4/8
pools.dat %>% 
  transmute(richness=pool_rich,"Mean age"=age_mean,"Variance age"=age_var,FPdist=FPdist) %>%
  gather("metric","X",2:3) %>% 
  ggplot(aes(X,FPdist)) +
  geom_point(aes(size=richness),shape=19,alpha=0.4) +
  stat_smooth(method="lm",colour="black") +
  facet_wrap(~metric,scales = "free", strip.position = "bottom") +
  labs(x=NULL, y="Functional / phylogenetic distance") +
  geom_text(data=pan.lett,aes(x=-Inf,y=Inf,label=lab),inherit.aes=F,hjust=-1,vjust=1) +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size=12),
        legend.position = "none")
ggsave("output/figs/Fig2_alt.pdf",height=4, width=8)  


#' # Standardized effect size of mean pairwise distance (ses-mpd)
#' Load community data from saved phyloseq object
phy <- readRDS("output/rds/phy.rds") 
#' extract otu table and apply hellenger transform
otu.tab <- otu_table(phy) %>% data.frame %>% decostand("hellinger")
#' Add a column for missing otus (never observed)
#missing <- colnames(pools)[!(colnames(pools) %in% colnames(otu.tab))]
#for(i in 1:length(missing)){otu.tab[[missing[i]]] <- 0}

#' Filter pools to only include those remaining in the experiment (some lost due to poor sequencing)
pools.remaining <- rownames(pools)[rownames(pools) %in% sample_data(phy)$poolID]
#' Calculate standardized effect size of mpd
ses.func <- function(otu.tab,dist.mx){
  foreach(i=pools.remaining, .combine=rbind) %dopar% {
    samps <- sample_names(phy)[which(sample_data(phy)$poolID==i)]
    spp <- pools[i,] %>% unlist %>% .[.] %>% names
    ses.mpd(otu.tab[samps,spp],dist.mx[spp,spp],null.model="taxa.labels",abundance.weighted=T)
  }
}
ses.dat <- ses.func(otu.tab,FPdist)

#' test variation in ses.mpd
ses.dat %<>% merge(data.frame(sample_data(phy)),by=0) %>% drop_na
(ses.mod <- lme(mpd.obs.z~pool_rich+age_mean+age_var,data=ses.dat, random=~1|poolID,method="ML")) %>% 
  drop1(test="Chisq")

#' summarize means of species pool reps for plotting trend line
ses.dat.means <- ses.dat %>% group_by(poolID,age_mean,pool_rich,age_var) %>% 
  summarize(mpd.obs.z=mean(mpd.obs.z))
#' Plot SES.mpd results
#+  fig.align="center", fig.asp=8/11, out.width="75%"
ggplot(ses.dat, aes(x=age_mean,y=mpd.obs.z))+
  geom_hline(aes(yintercept=0),alpha=0.3,linetype="dashed",size=0.5)+
  xlim(c(min(ses.dat$age_mean),max(ses.dat$age_mean)))+
  stat_smooth(data=ses.dat.means,method="lm",colour="black",size=0.5) +
  geom_point(aes(size=pool_rich),shape=19,alpha=0.3) +
  scale_size(range = c(0.5,3.5))+
  labs(y="Standardized effect size of\nmean pairwise distance",
       x=expression("Mean "*log[10]*"(age)")) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_text(size=11),
        axis.text = element_text(size=8),
        axis.line = element_line(size=0.3),
        axis.ticks = element_line(size=0.12),
        axis.ticks.length = unit(.05, "cm"))
ggsave("output/figs/Fig3.pdf",height=8, width=11, units="cm")  
