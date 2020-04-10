#Designing mock communities
#uses one of each unique ITS2 97% OTU included in spPool Exp, also including ISO-57 (not in exp).
#All species diluted to 1ng/ul to start
#each spp combined with 1/1000 ISO-57, given unique barcode for Illuimna (maybe include @ 1/2 concentration).
#Each sppecies serially diluted 1:4, 7 times and included in 14 of 16 mock communities, twice at each conc.
  
library(tidyverse)
library(magrittr)

isoIN <- read.csv("~/projects/Dissertation/SpPools/cultures_20160705.csv") %>% filter(InMoc==1)
Spp <- isoIN$ISO_ID %>% as.character() %>% rep(each=7)

set.seed(8080)

Concentrations <- c(10)
for (i in 1:6){
  tmp <- tail(Concentrations,1)/4
  Concentrations %<>% c(tmp)
}
SppConcentrations <- Concentrations %>% rep(times=35)
PoolConcentrations <- Concentrations %>% rep(each=5)


Global <- data.frame(Spp=Spp,Conc=SppConcentrations,index=(1:length(Spp)))
MocComs1 <- data.frame(Spp=NULL,Conc=NULL,index=NULL)
for (i in 1:7){
  print(i)
  repeat{
    tmpOUT <- NULL
    for (j in PoolConcentrations){
      repeat{
        tmp <- Global %>% filter(Conc==j) %>%
          filter(!(Spp %in% tmpOUT$Spp)) 
        if (nrow(tmp)!=0){
          tmp %<>% sample_n(1)} else{
            break}
        if (!(tmp$Spp %in% tmpOUT$Spp)){
          tmpOUT %<>% rbind(tmp);
          break}
      }
    }
    if (is.null(tmpOUT)){next}
    if (nrow(tmpOUT)<35){next}
    if (nrow(tmpOUT)==35){break}
    }
  Global <- Global[!(Global$index %in% tmpOUT$index),]
  tmpOUT$index <- i
  MocComs1 <-  bind_rows(MocComs1,tmpOUT)
}
ComMatrix1 <- MocComs1 %>% spread(index,Conc)

Global <- data.frame(Spp=Spp,Conc=SppConcentrations,index=(1:length(Spp)))
MocComs2 <- data.frame(Spp=NULL,Conc=NULL,index=NULL)
for (i in 1:7){
  print(i)
  repeat{
    tmpOUT <- NULL
    for (j in PoolConcentrations){
      repeat{
        tmp <- Global %>% filter(Conc==j) %>%
          filter(!(Spp %in% tmpOUT$Spp)) 
        if (nrow(tmp)!=0){
          tmp %<>% sample_n(1)} else{
            break}
        if (!(tmp$Spp %in% tmpOUT$Spp)){
          tmpOUT %<>% rbind(tmp);
          break}
      }
    }
    if (is.null(tmpOUT)){next}
    if (nrow(tmpOUT)<35){next}
    if (nrow(tmpOUT)==35){break}
  }
  Global <- Global[!(Global$index %in% tmpOUT$index),]
  tmpOUT$index <- i
  MocComs2 <-  bind_rows(MocComs2,tmpOUT)
}
ComMatrix2 <- MocComs2 %>% spread(index,Conc)

ComMatrix <- bind_cols(ComMatrix1,ComMatrix2[,-1])
colnames(ComMatrix) <- c("Spp",1:14)

#sanity check
colSums(ComMatrix[,-1])
#Total DNA(ng) added if 5 ul of each
colSums(ComMatrix[,-1])[1]*5
#Total volume = 35*5 = 175ul
#Start with 25ul H2O to bring total to 200ul
(colSums(ComMatrix[,-1])[1]*5)/200

#Final Conc = 1.67 ng/ul

#Convert DNA conc to dilution series stage
uniqCon <- unique(PoolConcentrations)
ComMatrix_stage <- ComMatrix
for (i in 1:length(uniqCon)){
  ComMatrix_stage[ComMatrix==uniqCon[i]] <- i
}


#write output to file
write.csv(ComMatrix,"~/projects/Dissertation/SpPools/MocCommunities_conc.csv")
write.csv(ComMatrix_stage,"~/projects/Dissertation/SpPools/MocCommunities_stage.csv")
