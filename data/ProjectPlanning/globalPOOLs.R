library(magrittr)
library(dplyr)
library(tidyr)
library(dissUtils)
library(ggplot2)
library(grid)
library(gridExtra)
library(parallel)
library(gnumeric)
library(gtools)

#read in data and subset only ISOs that will be included in sppool experiment
OTUsims <- read.gnumeric.sheet("/home/ubuntu/data/cultures_20151020.ods", head=T, sheet.name = "ISO.OTUs") %>%
  subset(GlobalPool_trm==1)
rownames(OTUsims) <- NULL

#define integer identifiers
OTUnum <- OTUsims$ITS2_97OTUs %>% substring(4) %>% as.integer()
OTUsims$OTUnum <- OTUnum
IDnum <- seq(1:nrow(OTUsims))
OTUsims$IDnum <- IDnum
uniqueOTUs <- unique(OTUnum)

#LSAG site age
lsag <- list(A=300,B=2100,C=20000,D=150000,F=4100000) %>% lapply(log)
yearsPool <- lsag[OTUsims$site] %>% unlist

#set up pool paramters
maxPool <- 30
minPool <- 2


possReps <- 1000000 #how many pools to select in each chain

chains<-20

#Run parallel chains
parallel.out <- mclapply(1:chains, mc.cores = chains, FUN = function(chain){ 
  out <- matrix(nrow=possReps,ncol=maxPool)
  tmp <- seq(minPool,maxPool) %>% sample(possReps,replace=T)                    #choose random pool sizes
  for (i in 1:possReps){
    out[i,1:tmp[i]] <- sample(uniqueOTUs,tmp[i]) %>% 
      sapply(function(x){ IDnum[OTUnum==x] %>% .[sample.int(length(.),1)] }) %>% 
      sort}
  out %<>% unique
  return(out)
})
ISOids <- do.call(rbind,parallel.out)

ISOids %<>% unique

chains2 <- 5

parallel.out <- mclapply(1:chains2, mc.cores = chains2, FUN = function(chain){
  low <- 1+((chain-1)*ceiling(nrow(ISOids)/chains2))
  if(chain==chains2){high <- nrow(ISOids)} else{high <- chain*ceiling(nrow(ISOids)/chains2)}
  tmp <- ISOids[low:high,]
  out <- apply(tmp,1,FUN=function(x){
    yearsPool[x] %>% var(na.rm=T)
  })
  return(out)})
possVar <- do.call(c,parallel.out)

possVar[possVar %>% is.na] <- 0

parallel.out <- mclapply(1:chains2, mc.cores = chains2, FUN = function(chain){
  low <- 1+((chain-1)*ceiling(nrow(ISOids)/chains2))
  if(chain==chains2){high <- nrow(ISOids)} else{high <- chain*ceiling(nrow(ISOids)/chains2)}
  tmp <- ISOids[low:high,]
  out <- apply(tmp,1,FUN=function(x){
    yearsPool[x] %>% mean(na.rm=T)
  })
  return(out)})
possMean <- do.call(c,parallel.out)

parallel.out <- mclapply(1:chains2, mc.cores = chains2, FUN = function(chain){
  low <- 1+((chain-1)*ceiling(nrow(ISOids)/chains2))
  if(chain==chains2){high <- nrow(ISOids)} else{high <- chain*ceiling(nrow(ISOids)/chains2)}
  tmp <- ISOids[low:high,]
  out <- apply(tmp,1,FUN=function(x){
    maxPool - sum(is.na(x))
  })
  return(out)})
numSPP <- do.call(c,parallel.out)

#calculate bins for 4 reps
poolWidth <- 1.5
varWidth <- 8.9
meanWidth <- 3.9
minNum <- minPool
maxNum <- minPool+poolWidth
Bins <- NULL
for (i in (1:((maxPool-minPool+1)/poolWidth))){
  varTMP <- possVar[numSPP>=minNum & numSPP<maxNum]
  lowerV <- min(varTMP,na.rm=T)
  upperV <- min(varTMP,na.rm=T)+varWidth
  while(max(upperV)<(max(varTMP,na.rm=T)) & max(upperV)<31){
    lowerV <- c(lowerV,tail(upperV,1))
    upperV <- c(upperV,(tail(upperV,1)+varWidth))
  }
  lowerV[1] <- -Inf
  upperV[length(upperV)] <- Inf
  varBins <- cbind(minNum,maxNum,lowerV,upperV)
  for (j in (1:nrow(varBins))){
    meanTMP <- possMean[numSPP>=minNum & numSPP<maxNum & 
                          possVar>=varBins[j,3] & possVar<varBins[j,4]]
    lowerM <- min(meanTMP)
    upperM <- min(meanTMP)+meanWidth
    while(max(upperM)<(max(meanTMP))){
      lowerM <- c(lowerM,tail(upperM,1))
      upperM <- c(upperM,(tail(upperM,1)+meanWidth))
    }
    tmplen <- length(upperM)
    lowerM[1] <- -Inf
    upperM[tmplen] <- Inf
    tmp <- cbind(minNum=rep(minNum,tmplen),maxNum=rep(maxNum,tmplen),
                 minVar=rep(varBins[j,3],tmplen),maxVar=rep(varBins[j,4],tmplen),
                 minMean=lowerM,maxMean=upperM)
    Bins <- do.call(rbind,list(Bins,tmp))  
    
  }
  minNum <- maxNum
  maxNum <- maxNum + poolWidth
}
Bins[Bins[,1]==min(Bins[,1]),1] <- minPool
Bins[Bins[,2]==max(Bins[,2]),2] <- maxPool

#calculate bins for 2 reps
poolWidth <- 1.2
varWidth <- 6.3
meanWidth <- 3.2
minNum <- minPool
maxNum <- minPool+poolWidth
Bins150 <- NULL
for (i in (1:((maxPool-minPool+1)/poolWidth))){
  varTMP <- possVar[numSPP>=minNum & numSPP<maxNum]
  lowerV <- min(varTMP,na.rm=T)
  upperV <- min(varTMP,na.rm=T)+varWidth
  while(max(upperV)<(max(varTMP,na.rm=T)) & max(upperV)<31){
    lowerV <- c(lowerV,tail(upperV,1))
    upperV <- c(upperV,(tail(upperV,1)+varWidth))
  }
  lowerV[1] <- -Inf
  upperV[length(upperV)] <- Inf
  varBins <- cbind(minNum,maxNum,lowerV,upperV)
  for (j in (1:nrow(varBins))){
    meanTMP <- possMean[numSPP>=minNum & numSPP<maxNum & 
                          possVar>=varBins[j,3] & possVar<varBins[j,4]]
    lowerM <- min(meanTMP)
    upperM <- min(meanTMP)+meanWidth
    while(max(upperM)<(max(meanTMP))){
      lowerM <- c(lowerM,tail(upperM,1))
      upperM <- c(upperM,(tail(upperM,1)+meanWidth))
    }
    tmplen <- length(upperM)
    lowerM[1] <- -Inf
    upperM[tmplen] <- Inf
    tmp <- cbind(minNum=rep(minNum,tmplen),maxNum=rep(maxNum,tmplen),
                 minVar=rep(varBins[j,3],tmplen),maxVar=rep(varBins[j,4],tmplen),
                 minMean=lowerM,maxMean=upperM)
    Bins150 <- do.call(rbind,list(Bins150,tmp))  
    
  }
  minNum <- maxNum
  maxNum <- maxNum + poolWidth
}
Bins150[Bins150[,1]==min(Bins150[,1]),1] <- minPool
Bins150[Bins150[,2]==max(Bins150[,2]),2] <- maxPool

#calculate bins for no reps
poolWidth <- 1.3
varWidth <- 3.2
meanWidth <- 2
minNum <- minPool
maxNum <- minPool+poolWidth
Bins300 <- NULL
for (i in (1:((maxPool-minPool+1)/poolWidth))){
  varTMP <- possVar[numSPP>=minNum & numSPP<maxNum]
  lowerV <- min(varTMP,na.rm=T)
  upperV <- min(varTMP,na.rm=T)+varWidth
  while(max(upperV)<(max(varTMP,na.rm=T)) & max(upperV)<31){
    lowerV <- c(lowerV,tail(upperV,1))
    upperV <- c(upperV,(tail(upperV,1)+varWidth))
  }
  lowerV[1] <- -Inf
  upperV[length(upperV)] <- Inf
  varBins <- cbind(minNum,maxNum,lowerV,upperV)
  for (j in (1:nrow(varBins))){
    meanTMP <- possMean[numSPP>=minNum & numSPP<maxNum & 
                          possVar>=varBins[j,3] & possVar<varBins[j,4]]
    lowerM <- min(meanTMP)
    upperM <- min(meanTMP)+meanWidth
    while(max(upperM)<(max(meanTMP))){
      lowerM <- c(lowerM,tail(upperM,1))
      upperM <- c(upperM,(tail(upperM,1)+meanWidth))
    }
    tmplen <- length(upperM)
    lowerM[1] <- -Inf
    upperM[tmplen] <- Inf
    tmp <- cbind(minNum=rep(minNum,tmplen),maxNum=rep(maxNum,tmplen),
                 minVar=rep(varBins[j,3],tmplen),maxVar=rep(varBins[j,4],tmplen),
                 minMean=lowerM,maxMean=upperM)
    Bins300 <- do.call(rbind,list(Bins300,tmp))  
    
  }
  minNum <- maxNum
  maxNum <- maxNum + poolWidth
}
Bins300[Bins300[,1]==min(Bins300[,1]),1] <- minPool
Bins300[Bins300[,2]==max(Bins300[,2]),2] <- maxPool