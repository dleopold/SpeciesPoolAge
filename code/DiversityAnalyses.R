#' ---
#' title: Analyses of species pool richness / diversity
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
library(pals)
library(ggthemes)
library(vegan)
library(nlme)
library(foreach)
library(doMC)
#use all resources for parallel processing
registerDoMC(cores=parallel::detectCores())

#' ### Load data and calculate richness and diversity metrics.  
#'
#' * Richness = mean richness from 1000 rarefactions.   
#' * Diversity = Inverse Simpson's Diversity of log2 transformed abundances. 
phy <- readRDS("output/rds/phy.rds") #Load from saved phyloseq object
otu.tab <- otu_table(phy) %>% data.frame #extract otu table
#Parallel loop to calculate rarefied richness
richness <- foreach(i=1:1000) %dopar% {
  otu.tab %>% rrarefy(.,min(rowSums(.))) %>% decostand("pa") %>% rowSums()
}
sample_data(phy)$Hill0 <- Reduce("+", richness) / length(richness) 

#' Extract sample data for analyses
dat <- sample_data(phy) %>% data.frame() 
dat$poolID %<>% factor()

#' ### Linear vs. polynomial 

#' Test for curvilinear relationship between pool richness and realized richness
glmm.m0 <- lme(Hill0 ~ pool_rich + I(pool_rich^2), random=~1|poolID, data=dat, weights=varFixed(~pool_rich),method="ML")
glmm.m0lin <- lme(Hill0 ~ pool_rich, random=~1|poolID, data=dat, weights=varFixed(~pool_rich),method="ML")
anova(glmm.m0,glmm.m0lin)

#' Test whether the shape of the curve depends on species pool age and age-variance
#' This shows that the effect of richness depends on species pool age but not age-variance. We will not interrogate this more here because we will focus on the gnls model below.
lme(Hill0 ~ pool_rich + I(pool_rich^2)*age_mean + I(pool_rich^2)*age_var, random=~1|poolID, data=dat, weights=varFixed(~pool_rich),method="ML")  %>% drop1(test="Chisq")

#' ### Fit generalized non-linear least-squares model (Michaelis Menten)
#find starting values
Vm.H0 <- getInitial(Hill0 ~ SSmicmen(pool_rich, max(Hill0), 1), data = dat)[1]
K.H0 <- getInitial(Hill0 ~ SSmicmen(pool_rich, max(Hill0), 1), data = dat)[2]
#fit gnls model
gnls.H0 <- gnls(Hill0 ~ SSmicmen(pool_rich, Vm.H0, K.H0), data = dat, weights=varFixed(~pool_rich))
#Extract residuals from NLS fit of Michaelis-Menten curve
dat$gnls.H0.resid <- residuals(gnls.H0,type="pearson")
#Show model parameters
intervals(gnls.H0)

#' ### Model gnls residuals
#' Full model with richness:age and richness:variance
lme(gnls.H0.resid ~ 0+pool_rich*age_mean+pool_rich*age_var, random=~1|poolID, data=dat,method="ML") %>% drop1(test="Chisq")
#' Drop richness:variance but test direct effect of variance
lme(gnls.H0.resid ~ 0+pool_rich*age_mean+age_var, random=~1|poolID, data=dat,method="ML") %>% drop1(test="Chisq")
#' Refit final model
(gnls.H0.residMod.avg <- lme(gnls.H0.resid ~ 0+pool_rich*age_mean, random=~1|poolID, data=dat,method="REML")) %>% summary

#' ### Make plots
#' 
#' **Plotting function**

gnls.plot.func <- function(gnls.mod, lme.resid.mod,y_frame,y_lab){
  #initiate data frame for model predictions across range of species pool diversity
  predictions <- data.frame(pool_rich=rep(2:30,2))
  #get mean prediction
  predictions$mean <- predict(gnls.mod,predictions)
  #find predictions for high and low mean age
  age_hi <- quantile(dat$age_mean,probs=0.95)
  age_lo <- quantile(dat$age_mean,probs=0.05)
  predictions %<>% mutate(age_mean = rep(c(age_hi,age_lo),each=29),
                          ageLab = rep(c("Old","Young"),each=29),
                          poolID = "foo")
  predictions$resid = predict(lme.resid.mod,predictions,level=0)
  #loess model to convert pearsons to raw residuals for plotting curves
  nls.resid <- residuals(gnls.mod,type="pearson")
  std.function <- loess(attr(nls.resid,"std")~dat$pool_rich)
  predictions$conversion <- predict(std.function,predictions$pool_rich)
  predictions$HiLo = predictions$mean + predictions$resid * predictions$conversion
  
  #prepare raw data
  dat.raw <- data.frame(X=dat$pool_rich, Y=gnls.mod$fitted + gnls.mod$residuals, age_mean=dat$age_mean)
  #make range frame data
  rng.frame <- data.frame(X=c(0,30),Y=c(0,y_frame))
  ###PLOT###
  plot.out <- ggplot(predictions, aes(x=pool_rich, y=mean)) +
    geom_rangeframe(data=rng.frame, aes(x=X,y=Y),size=0.3) +
    geom_point(data=dat.raw, aes(x=X,y=Y,fill=age_mean),
               position=position_jitter(),shape=21,size=1.8,color="black",stroke=0.2) +
    scale_fill_gradientn(colours = coolwarm(1000)[c(rep(1,100),1:1000,rep(1000,100))],name="Species\npool\nage") +
    stat_smooth(method=loess,linetype="dashed",color="black",size=0.5) +
    stat_smooth(aes(color=ageLab,y=predictions$HiLo),method=loess,se=F,size=0.5) +
    scale_color_manual(values=coolwarm(10)[c(9,2)],name="Effect of\nspecies pool\nage") +
    labs(y=y_lab,x="Species pool richness") +
    geom_abline(slope=1,intercept=0,linetype="dashed",color="grey50",size=0.5) +
    guides(colour=F, fill = guide_colorbar(ticks = FALSE)) +
    theme_classic() +
    theme(axis.title = element_text(size=11),
          axis.text = element_text(size=8),
          axis.ticks = element_line(size=0.12),
          axis.ticks.length = unit(.05, "cm"),
          axis.line=element_blank(),
          legend.text = element_text(size=8),
          legend.title = element_text(size=8),
          legend.key.width=unit(0.3,"cm"),
          legend.key.height=unit(0.5,"cm"))
  plot.out
}
#+ fig.align="center", out.width="90%", fig.asp=5.5/11
gnls.plot.func(gnls.H0,gnls.H0.residMod.avg,20,"Observed richness")
ggsave("output/figs/Fig1.pdf",height=5.5, width=11, units="cm") 

#' ### Dependencies
sessionInfo()
