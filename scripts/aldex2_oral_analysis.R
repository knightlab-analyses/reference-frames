library(tidyverse)
library(ALDEx2)
set.seed(148)
# choose this directory
setwd("/Users/jmorton/Documents/dev/reference-frames/scripts/benchmark_scripts")

# load data ---------------------------------------------------------------

map <- read.delim("../../data/oral_trimmed_metadata.txt", row.names=1)
otu <- read.delim("../../data/oral_trimmed_deblur.txt", row.names=1)
#otu <- read.delim("../../data/oral_trimmed_deblur.txt", row.names=1)

# helper functions --------------------------------------------------------

#' Filter summary for Aldex2
#' @param s output of \code{summary_aldex2}
#' @param pval (pvalue threshold)
sig_aldex2 <- function(s, pval=0.05){
  filter(s, padj < pval)
}

#' Run Aldex2 on simulated dataset
#'
#' @param dat output from either \code{create_true_abundance} or
#'   \code{resample_data}
#' @return output of aldex
run_aldex2 <- function(dat){
  countdata <- t(dat[,-1,drop=F])
  colnames(countdata) <- paste0("n", 1:ncol(countdata))
  aldex.fit <- aldex(countdata, as.character(dat$Condition))
  return(aldex.fit)
}

#' Summarise DE from Aldex2 models
#' @param fit output of run_aldex2
#' @param prob adjusted pvalue threshold 
#' @return data.frame with columns DE, low, and high
summary_aldex2 <- function(fit){
  fit %>% 
    as.data.frame() %>% 
    rownames_to_column("category") %>% 
    select(category, effect, wi.eBH) %>% 
    mutate(padj=wi.eBH) %>% 
    mutate(mean=effect) %>% 
    mutate(low=NA, high=NA)
}


# analysis ----------------------------------------------------------------
d <- data.frame(Condition = map$brushing_event, otu)

fit <- run_aldex2(d)
sfit <- summary_aldex2(fit)
sig_aldex2(sfit)

# Main result is it says nothing is significant

# now just list differential abundance
data.frame(sfit$category, sfit$effect, sfit$padj) %>% 
  write.table(file="../../results/oral-results/aldex2_results.txt")


