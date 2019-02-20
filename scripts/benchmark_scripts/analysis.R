library(tidyverse)
source("ancom.R")

#=======================================================================================================
# Section 0: Meta data and feature table
#=======================================================================================================
# Meta data
meta.dat=read_delim("metadata-toothbrush.txt", delim = "\t")
meta.dat=as.data.frame(meta.dat)
meta.dat$Sample.ID=meta.dat$`#SampleID`
meta.dat=meta.dat%>%select(Sample.ID,  brushing_event)%>%
  mutate(subject=rep(paste0("subject", seq(nrow(meta.dat)/2)), each=2))

# Feature table
feature.table.origin=read_delim("trimmed_table.txt", delim = "\t", skip=1)
OTU.name=feature.table.origin$`#OTU ID`
feature.table=data.frame(Sample.ID=colnames(feature.table.origin)[-1], 
                         t(as.matrix(feature.table.origin[, -1])), row.names = NULL)
colnames(feature.table)=c("Sample.ID", OTU.name)
feature.table=as.data.frame(feature.table)

# Run ANCOM: Friedman Rank Sum Test
start_time <- Sys.time()
res.W=ANCOM.main(OTUdat=feature.table, Vardat=meta.dat, adjusted=F, repeated=T,
                 main.var="brushing_event", adj.formula=NULL, repeat.var="subject", longitudinal=F,
                 random.formula=NULL, multcorr=2, sig=0.05, prev.cut=0.90)
end_time <- Sys.time()
end_time - start_time

res.ANCOM=res.W$W.taxa
write_csv(res.ANCOM, "ANCOM results_friedman.csv")

# Run ANCOM: Mixed-Effects Model
start_time <- Sys.time()
res.W=ANCOM.main(OTUdat=feature.table, Vardat=meta.dat, adjusted=F, repeated=T,
                 main.var="brushing_event", adj.formula=NULL, repeat.var="subject", longitudinal=T,
                 random.formula="~1|subject", multcorr=2, sig=0.05, prev.cut=0.90)
end_time <- Sys.time()
end_time - start_time

res.ANCOM=res.W$W.taxa
write_csv(res.ANCOM, "ANCOM results_mixed.csv")
