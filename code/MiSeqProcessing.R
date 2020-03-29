#' ---
#' title: Processing raw MiSeq Data
#' author: Devin R Leopold
#' date: March 22, 2020
#' output:
#'    html_document:
#'      toc: true
#'      toc_float: true
#'      self_contained: true
#'      highlight: zenburn
#' ---

#' Sample types are indicated by the prefix of the sample IDs in the raw reads fastq.gz file names
#'  *Sp - Experimental species pool samples
#'  *N - Negative controls
#'  *M - Mock communities
#"  *F - Single species samples prepared from pure cultures of isolates used in the experiment

#' ## Prepare environment

#' Load packages
#+ warning=F, message=F
library(tidyverse)
library(magrittr)
library(dada2)
library(Biostrings)
library(foreach)
library(knitr)
library(phyloseq)

#' Load fasta of target amplicons for all isolates in experiment
targets <- readDNAStringSet("data/AllIsolates_ITS2_200bp.fasta")
#and the host plant (Vaccinium calycinum)
host <- readDNAStringSet("data/VACCAL_ITS2_200bp.fasta")

#' key to which isolates were included in experimental species pools
pools <- read.csv("data/SpeciesPoolKey.csv", as.is=T, header=T) %>% select(-num,-var,-mean)

#' key linking sample IDs to species pools
sampKey <- read.csv("data/PlantData.csv", as.is=T, header=T)

#' key to mock communities
mock <- read.csv("data/MockCommunities.csv", as.is=T, header=T) 

#' keys to isolates used in expt
sites <- read.csv("data/sites.csv", as.is=T)
isoKey <- read.csv("data/cultureData.csv",as.is=T) %>% left_join(sites) %>% mutate(site_logAge=log10(site_age)) %>%
  .[match(names(targets),.$isoID),]

#' Sample names
samps <- list.files("data/MiSeq/", pattern="R1_001.fastq.gz") %>% strsplit(., "_",fixed=T) %>% sapply(., `[`, 1) %>% sort
#Skip single species sample for ISO.57 which was not included in the experiment
samps %<>% .[samps!="F57"]

#' Paths to fwd and rev raw reads
R1.raw <- list.files("data/MiSeq/", pattern="R1_001.fastq.gz",full.names = T) %>% sort
R2.raw <- list.files("data/MiSeq/", pattern="R2_001.fastq.gz",full.names = T) %>% sort

#' *****

#' ## Trim gene primers

#' make directory for cutadapt output
dir.create("output/MiSeq/cutadapt", recursive = T, showWarnings = F)

#' loop over raw paired read files
for(samp in samps){
  #skip trimming if all files have already been processed
  if(samp==samps[1] & length(list.files("output/MiSeq/cutadapt"))==2*length(samps)){break}
  #path to current raw reads
  FWD <- paste0(samp,"_") %>% grep(R1.raw, value=T)
  REV <- paste0(samp,"_") %>% grep(R2.raw, value=T)
  
  #output paths
  FWD.out <- file.path("output/MiSeq/cutadapt",paste0(samp,".R1.trim.fastq.gz"))
  REV.out <- file.path("output/MiSeq/cutadapt",paste0(samp,".R2.trim.fastq.gz"))
  
  #set cutadapt parameters
  cutadapt.flags <- paste("-g GCATCGATGAAGAACGC", #FWD gene primer
                          "-G TCCGCTTATTGATATGC",  #REV gene primer
                          "--overlap 17",
                          "-e 0.2",
                          "-o ", FWD.out,
                          "-p ", REV.out,
                          "--discard-untrimmed",
                          FWD, REV)
  system2("/home/harpua/miniconda3/bin/cutadapt", args = cutadapt.flags)
}

#' *****

#' ## Filter and crop

#' Initial screening suggests that quality drops off drastically at 150 bp (dramatically for R2) limiting our ability to successfully denoise and merge reads. Since only known sequences are expected and diversity should be low in the microcosms it should be fine to just use R1 and cropping at 200bp, which retains sufficient sequence length to distinguish all isolates. 

# Plot example quality profile
plotQualityProfile(c("output/MiSeq/cutadapt/Sp99.R1.trim.fastq.gz","output/MiSeq/cutadapt/Sp99.R2.trim.fastq.gz"))

#' Identify paths to fwd trimmed reads.
R1.trimmed <- list.files("output/MiSeq/cutadapt", pattern="R1.trim",full.names = T)  %>% sort

#' Create directory and paths for filtering output
dir.create("output/MiSeq/filtered", showWarnings = F)
R1.filtered <- gsub("output/MiSeq/cutadapt","output/MiSeq/filtered", R1.trimmed)  %>% sort
names(R1.filtered) <- samps

#' filter low quality reads and truncate to 200 bp
#+ cache=T
filterAndTrim(R1.trimmed, R1.filtered,
              maxN=0,
              truncLen=200,
              maxEE=2,
              multithread = T)

#' *****

#' ## Denoise with DADA2 

#' create directory for output
dir.create("output/MiSeq/dada2", showWarnings = F)

#' learn errors 
#+ cache=T
errF <- learnErrors(R1.filtered, multithread = T, randomize = T, nbases = 1e9)
plotErrors(errF, nominalQ=TRUE)

#' denoise samples one at a time setting the expected sequences as priors for each sample
#+ cache=T
denoised.list <- lapply(samps,function(x) NULL) %>% setNames(samps)
for(samp in samps){
  
  # Identify expected amplicons (dada priors)
  identifier <- substring(samp,1,1)
  if(identifier=="N"){priors <- host} #Negative controls
  if(identifier=="S"){                #Experimental species pool samples
    pool <- sampKey$poolID[sampKey$sampID==samp]
    expected <- pools %>% filter(poolID==pool) %>% grep("ISO",.,value = T)
    if(pool=="C"){priors <- host} else{priors <- targets[expected] %>% c(host)}
  }
  if(identifier=="M"){                #Mock community sample
    priors <- targets[mock$Spp]%>% c(host)
  }
  if(identifier=="F"){                #Pure culture samples
    expected <- substring(samp,2,4) %>% paste0("ISO.",.)
    priors <- targets[expected] %>% c(host)
  }
  
  #denoise current samp
  samp.path <- paste0(samp,".R1") %>% grep(R1.filtered, value=T)
  # OMEGA_A reduced from default (1e-40) slightly after preliminary check of single species samples to improve denoising
  denoised <- dada(samp.path, err=errF, multithread=TRUE, priors=priors, verbose=F, OMEGA_A=1e-45, MIN_ABUNDANCE=2)
  
  #make OTU table 
  otu.table <- makeSequenceTable(denoised)
  
  #remove chimeras
  #otu.table %<>% removeBimeraDenovo(multithread=T)
  
  #Identify and rename expected sequences
  #Initial trials with only dada2 "exact" sequence variants did not correctly resolve all expected taxa
  #Collapse sequence variants with < 98% similarity (4bp) difference
  for(i in 1:ncol(otu.table)){
    dists <- stringdist::stringdist(colnames(otu.table)[i],as.character(priors))
    if(min(dists)<4){
      colnames(otu.table)[i] <- names(priors)[order(dists)[1]]}
  }
  otu.table %<>% collapseNoMismatch(identicalOnly = T) #collapse otus assigned to the same expected seq
  #add sample ID
  rownames(otu.table) <- samp
   
  #add to list of denoised OTU tables
  denoised.list[[samp]] <- otu.table
}
 
#' merge otu tables
full.otu.table <- mergeSequenceTables(tables=denoised.list)

#' Relabel "unexpected" OTUs. This includes OTUs that were not in any of the inoculations or expected OTUs that showed up somewhere they were not expected (ie in a microcosm where they were not inoculated)
otuCounter <- 1
for(i in 1:ncol(full.otu.table)){
  seq <- colnames(full.otu.table)[i]
  if(nchar(seq)<20){next}
  if(seq %in% as.character(targets)){
    ID <- targets[match(seq, as.character(targets))] %>% names %>% .[1] %>% paste0("x",.)
  } else { 
    ID <- paste0("xOTU.",otuCounter) 
    seq %<>% DNAStringSet
    names(seq) <- ID
    #write unexpected sequences to fasta
    if(otuCounter==1){
      writeXStringSet(seq,"output/MiSeq/dada2/unexpected.fasta")} else{
        writeXStringSet(seq,"output/MiSeq/dada2/unexpected.fasta", append=T)}
    otuCounter %<>% add(1)} #update counter for next unexpected OTU
  colnames(full.otu.table)[i] <- ID
}

#' write the full denoised otu table to file
write.csv(full.otu.table,"output/csv/outTable.denoised.csv")

#' *****

#' ## Single species samples

#' The library included pure-culture samples of each of the 54 isolates used in the experiment. These samples can be used to verify the prior sangar sequencing of the isolates and will also serve as positive controls. Initial processing of these samples identified 3 isolates (F10, F61, F90) whoes ITS2 amplicon inferred from the denoised MiSeq data varied from earlier Sangar sequencing. The library of priors was updated and the MiSeq data processing was rerun.

#' subset single species samples and summarize the % of reads that map to the expected amplicon
singles <- full.otu.table %>% .[grepl("^F",rownames(.)),] %>% .[,colSums(.)>0]
foreach(i = rownames(singles), .combine = bind_rows) %do% { 
  isolate <- paste0("ISO.",substr(i,2,5))
  list(Isolate = isolate,
       NumberOfAmplicons = sum(singles[i,] > 0),
       SeqDepth = sum(singles[i,]),
       PctExpected = round(100*singles[i,isolate]/sum(singles[i,]),2))
} -> singles.summary
#' Almost all reads map to the correct amplicon.
#+ results='asis'
kable(singles.summary, caption = "Summary of single species samples")

#' *****

#' ## Negative controls
#' subset negative control samples
neg <- full.otu.table %>% .[grepl("^N",rownames(.)),] %>% .[,colSums(.)>0]
#' Primarily non-target amplicons present and in very low counts. 
colSums(neg)

#' there are also plant microcosm that were not inoculated 
controlPlants <- full.otu.table %>% .[rownames(.) %in% sampKey$sampID[sampKey$poolID=="C"],] %>% .[,colSums(.)>0]
#' these controls are dominated by host sequences and non-target otus (in low abundance).
colSums(controlPlants)

#' we can also look at how prevalant non-target otus are in the experimental samples
expt <- full.otu.table %>% .[grepl("^Sp",rownames(.)),] %>%
  .[!(rownames(.) %in% sampKey$sampID[sampKey$poolID=="C"]),] %>% .[,colSums(.)>0]
#' These otus make up a very small % of the total reads in any sample
round(100*rowSums(expt[,grepl("xOTU",colnames(expt))])/rowSums(expt),2) %>% sort(decreasing=T) %>% .[.>0]
#' Most are rare. The one that is fairly prevelant (xOTU.3) blasts to Vaccinium, so it is likely host sequences with errors, or some type of chimera.
colSums(expt[,grepl("xOTU",colnames(expt))] > 0) %>% sort(decreasing=T) 

#' Since there are so few unexpected amplicon and they are in such low abundance, we will just remove them (along with the host) for the remainder of the analyses
full.otu.table %<>% .[,grep("^ISO",colnames(.))]

#' *****

#' ## Mock comunity samples
#' The library included mock community samples where genomic DNA from pure culture isolates was mixed in known proportions. For isolates with similar amplicons (which were not included in the same species pools) only one representative was used in the mock communities, resulting in 35 taxa.
#' subset mock community samples and gather to long format
mock.obs <- full.otu.table %>% .[grepl("^M",rownames(.)),] %>% .[,colSums(.)>0] %>%
  apply(1,FUN=function(x){x/sum(x)}) %>% data.frame %>%
  rownames_to_column("Spp") %>%
  pivot_longer(-Spp,names_to = "pool", values_to = "observed")

#' make corresponding table of known concentrations
mock.known <- mock %>% pivot_longer(-Spp,names_to = "pool", values_to = "known")

#' observed vs. expected looks fairly good, with some evidence for possible species specific biases
full_join(mock.obs,mock.known) %>% 
  ggplot(aes(x=known,y=observed,color=Spp))+
  ggbeeswarm::geom_quasirandom()+
  scale_y_log10(name="observed proportional abundance")+
  scale_x_log10(name="known proportional DNA abundance")+
  theme_classic()+
  theme(legend.position = "none")
#' plot for individual species - clearly some taxa are not as easily detected, but overall there is very good evidence that read abundance tracks known abundance for any given taxa used in the experiment. 
full_join(mock.obs,mock.known) %>% 
  ggplot(aes(x=known,y=observed))+
  geom_point()+
  scale_y_log10(name="observed proportional abundance")+
  scale_x_log10(name="known proportional DNA abundance")+
  facet_wrap(~Spp)+
  theme_classic()+
  theme(legend.position = "none")

#' *****

#' ## Experimental samples
#' Convert to phyloseq object
otu_tab <- full.otu.table %>% .[grepl("^Sp",rownames(.)),] %>% otu_table(taxa_are_rows = F)
#' Set of poorly sequenced outlier samples with less than 1500 reads can be removed
data.frame(depth=sort(rowSums(otu_tab))) %>%
  mutate(index=1:nrow(.)) %>%
  ggplot(aes(y=depth,x=index))+
  geom_boxplot()+
  theme_classic()
otu_tab %<>% .[rowSums(otu_tab)>2000,]

#' ### Add metadata
#' convert species pool key into presence absence matrix
pools.mx <- matrix(nrow=nrow(pools),ncol=length(targets), dimnames = list(1:nrow(pools), names(targets)))
for(i in 1:nrow(pools)){
  foo <- pools[i,] %>% c
  pools.mx[i,] <- names(targets) %in% foo
}
#' write pool composition matrix to file
write.csv(pools.mx,"output/csv/pools.mx.csv")

#' summarize pool characteristics
pools.dat <- data.frame <- data.frame(poolID=as.character(pools$poolID),
                                      pool_rich=apply(pools.mx,1,function(x){ length(isoKey$site_logAge[x]) }),
                                      age_mean=apply(pools.mx,1,function(x){ mean(isoKey$site_logAge[x], na.rm=T) }),
                                      age_var=apply(pools.mx,1,function(x){ var(isoKey$site_logAge[x], na.rm=T) }),
                                      stringsAsFactors =F)
#' Add pool data to sample key
sampKey %<>% left_join(pools.dat) %>% column_to_rownames("sampID")
#' Write pool data to file
write.csv(pools.dat,"output/csv/pools.dat.csv",row.names = F)

#' add mothe plant ID (all half/full sibs)
sampKey$plantMother <- "A10"
sampKey$plantMother[(sampKey$plantID > 80 & sampKey$plantID <= 160) | sampKey$plantID == 321] <- "A11"
sampKey$plantMother[(sampKey$plantID > 160 & sampKey$plantID <= 240) | sampKey$plantID == 322] <- "A6"
sampKey$plantMother[(sampKey$plantID > 240 & sampKey$plantID <= 320) | sampKey$plantID == 324] <- "A30"
sampKey$plantMother %<>% factor

#' Combine otu table and sample data into phyloseq object
phy <- phyloseq(
  otu_tab,
  sample_data(sampKey)
)

#' Save phyloseq object
saveRDS(phy,"output/rds/phy.rds")

#' ## Session info
#' cutadapt
system2("/home/harpua/miniconda3/bin/cutadapt", args = "--version", stdout = T) 
#' R 
sessionInfo()
