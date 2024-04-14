#!/usr/bin/env Rscript
library(dplyr)

args<-commandArgs(trailingOnly=TRUE)
nSNPs<-args[1]
outfile<-args[2]
infiles<-args[-c(1,2)]

sd<-read.csv(nSNPs,sep="\t",header=FALSE)
colnames(sd)<-c("nSNPs","Rep","Depth","GlMethod","QsBeta","Platform","GcMethod")

hdr<-c("Sample","GQ","nDiscordant","nHomToHomDiscordant","nHomToHetDiscordant","nHetToHomDiscordant","nHetToHetDiscordant","nConcordant","nHomToHomConcordant","nHetToHetConcordant","nSitesComparedForSample","Rep","Depth","GlMethod","QsBeta","Platform","GcMethod")

d<-NULL
for(i in seq_along(infiles)){
di<-read.csv(infiles[i],sep="\t",header=FALSE)
colnames(di)<-hdr

di$nSitesGQ<-di$nDiscordant+di$nConcordant
di<-merge(di,sd,by=c("Rep","Depth","GlMethod","QsBeta","Platform","GcMethod"))

d<-rbind(d,di)
}

d%>% 
# sum for all samples + reps
group_by(GQ,Depth,GlMethod,QsBeta,Platform,GcMethod) %>%
summarise(
nDiscordant=sum(nDiscordant),
nConcordant=sum(nConcordant),
nSitesGQ=sum(nSitesGQ),
nHomToHomDiscordant=sum(nHomToHomDiscordant),
nHomToHetDiscordant=sum(nHomToHetDiscordant),
nHetToHomDiscordant=sum(nHetToHomDiscordant),
nHetToHetDiscordant=sum(nHetToHetDiscordant),
nHomToHomConcordant=sum(nHomToHomConcordant),
nHetToHetConcordant=sum(nHetToHetConcordant),
nSitesNon0DP=sum(nSitesComparedForSample),
nSNPs=sum(nSNPs)) %>% 
ungroup()-> sumd

save(sumd,file=outfile)

