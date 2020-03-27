# Process raw output from plate reader software for the biolog phenotypic microarray data. 

# Load packages
library(tidyverse)
library(magrittr)
library(foreach)
library(doMC)
# Use all resources for parallel processing
registerDoMC(cores=parallel::detectCores())

# Define function to extract biolog data from messy text output
parse_biolog <- function(x,wave=740){
  if(wave==740){wrng <- 16:27}else{wrng <- 3:15} #column ranges for two wavelengths
  tmp_data <- read.delim(x,header=F,skip=3,row.names=LETTERS[1:11])[-9:-11,wrng]
  colnames(tmp_data) <- paste0("x.",1:12)
  tmp_data <- rownames_to_column(tmp_data,"row")
  tmp_data_long <- gather(tmp_data,"column","abs",2:13)
  tmp_data_long$column <- substring(tmp_data_long$column,3)
  tmp_date <- strsplit(x,"/")[[1]][3] 
  tmp_file <- strsplit(x,"/")[[1]][4] 
  tmp_ISO_ID <- strsplit(tmp_file,".",fixed=T)[[1]][1]
  tmp_rep <- strsplit(tmp_file,".",fixed=T)[[1]][2]
  out <- data.frame(stringsAsFactors = F,
                    ISO_ID = tmp_ISO_ID,
                    rep = tmp_rep,
                    date = tmp_date,
                    tmp_data_long)
  return(out)
}
#paths to all of the biolog data
biolog_paths <- list.files("data/BiologData",full.names=T,recursive=T)
#import and concatenate all biolog data for absorbance at 740nm (not using 490nm)
biolog <- foreach(i=biolog_paths, .combine = rbind) %dopar% parse_biolog(i)

#' Modify columns 
biolog %<>% mutate(plate = paste(ISO_ID,rep,date),
                       date = lubridate::ymd(date),
                       well = paste0(row,column),
                       column = factor(column,levels=1:12),
                       date_num = as.numeric(date),
                       ISO_well=paste0(ISO_ID,"_",well),
                       filt=paste0(ISO_well,"_",rep))

#' Make figures for each isolate (skip if already done)
if(!dir.exists("output/figs/biolog")){
  dir.create("output/figs/biolog")
  foo <- foreach(i=unique(biolog$ISO_ID)[1:10]) %dopar% {
    out.file <- paste0("output/figs/biolog/ISO.",i,".pdf")
    p <- filter(biolog,ISO_ID==i) %>% 
      ggplot(aes(date,abs)) +
      geom_point() +
      stat_smooth(method="loess",se=F) +
      labs(title=paste0("ISO.",i),x="Time",y="Absorbance") +
      facet_grid(row~column) +
      theme_bw() + theme(axis.text = element_blank())
    ggsave(out.file, p, width=8, height = 6)
  }
}

# Remove outlier data for one isolate
biolog %<>% filter(filt!="63_H10_2" & filt!="63_H10_3")

# Calculate area under curve for 47 days of growth for each species / substrate
auc <- function(x){
  tmp_dat <- filter(biolog,ISO_well==x)
  loess_fit <- loess(abs~date_num,tmp_dat)
  loess_fun <- function(xx){predict(loess_fit,xx)}
  start <- min(tmp_dat$date_num)
  stop <- start + 47
  auc <- integrate(loess_fun,start,stop)$value
  ISO_ID <- tmp_dat$ISO_ID[1]
  well <- tmp_dat$well[1]
  out <- data.frame(ISO_ID,well,auc)
  return(out)
}
biolog_auc <- foreach(i=unique(biolog$ISO_well), .combine = rbind) %dopar% auc(i)

#' convert AUC data to wide matrix format
biolog_auc %<>% pivot_wider(names_from = well, values_from = auc) %>% column_to_rownames("ISO_ID")

#' standardize AUC data to blank well (no carbon substrate)
biolog_auc_blank <- biolog_auc[,1]
biolog_auc %<>% dplyr::select(-A1) %>% sweep(1,biolog_auc_blank)
biolog_auc[biolog_auc<=0] <- 0

#' save to file
rownames(biolog_auc) %<>% paste0("ISO.",.)
write.csv(biolog_auc,"output/csv/biolog.csv")


