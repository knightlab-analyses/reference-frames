#library(tidyverse)
source("ancom.R")

#=======================================================================================================
# Section 0: Meta data and feature table
#=======================================================================================================
# Meta data
meta.dat = read.delim("../../benchmarks/metadata2.txt", row.names=1)
meta.dat = as.data.frame(meta.dat)
meta.dat$Sample.ID = rownames(meta.dat)

# Feature table
feature.table = read.delim("../../benchmarks/rel_table2.txt", row.names=1)
OTU.name=colnames(feature.table)
feature.table=data.frame(Sample.ID=rownames(feature.table), 
                         as.matrix(feature.table), row.names = NULL)
colnames(feature.table)=c("Sample.ID", OTU.name)
feature.table=as.data.frame(feature.table)

# Run ANCOM: Friedman Rank Sum Test
start_time <- Sys.time()
res.W=ANCOM.main(OTUdat=feature.table, Vardat=meta.dat, adjusted=F, repeated=F,
                 main.var="labels", adj.formula=NULL, longitudinal=F,
                 random.formula=NULL, multcorr=2, sig=0.05, prev.cut=0.90)
end_time <- Sys.time()
end_time - start_time

res.ANCOM=res.W$W.taxa

write.csv(res.ANCOM, "../../benchmarks/ancom_rel_results2.csv")
