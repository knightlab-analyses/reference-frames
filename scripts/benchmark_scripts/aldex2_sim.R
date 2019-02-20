library(tidyverse)
library(phyloseq)
library(ALDEx2)

# this code was originally from Justin Silverman
setwd("/Users/jmorton/Dropbox (Simons Foundation)/oral-saliva/benchmarks")

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


# load data ---------------------------------------------------------------

map <- read.delim("metadata1.txt", row.names=1)
otu <- read.delim("rel_table1.txt", row.names=1)

# analysis ----------------------------------------------------------------

d <- data.frame(Condition = map$labels, otu) 
d
fit <- run_aldex2(d)
sfit <- summary_aldex2(fit)
sig_aldex2(sfit)

# Main result is it says nothing is significant

# now just list differential abundance
sfit %>% write.table(file="aldex2_results1.txt")

# load data ---------------------------------------------------------------

map <- read.delim("metadata2.txt", row.names=1)
otu <- read.delim("rel_table2.txt", row.names=1)

# analysis ----------------------------------------------------------------

d <- data.frame(Condition = map$labels, otu) 
d
fit <- run_aldex2(d)
sfit <- summary_aldex2(fit)
sig_aldex2(sfit)

# Main result is it says nothing is significant

# now just list differential abundance
sfit %>% write.table(file="aldex2_results2.txt")

