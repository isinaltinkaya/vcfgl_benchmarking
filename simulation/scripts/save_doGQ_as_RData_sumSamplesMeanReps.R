#!/usr/bin/env Rscript
library(dplyr)

args<-commandArgs(trailingOnly=TRUE)
outfile<-args[1]
infiles<-args[-c(1)]



hdr<-c("Sample","GQ","nDiscordant","nHomToHomDiscordant","nHomToHetDiscordant","nHetToHomDiscordant","nHetToHetDiscordant","nConcordant","nHomToHomConcordant","nHetToHetConcordant","nSitesComparedForSample","Rep","Depth","GlMethod","QsBeta","Platform","GcMethod")

d<-NULL
for(i in seq_along(infiles)){
di<-read.csv(infiles[i],sep="\t",header=FALSE)
colnames(di)<-hdr

di$nSitesGQ<-di$nDiscordant+di$nConcordant

d<-rbind(d,di)
}






d%>% 
# sum for all samples 
group_by(Rep,GQ,Depth,GlMethod,QsBeta,Platform,GcMethod) %>%
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
nSitesNon0DP=sum(nSitesComparedForSample)) %>% 
ungroup()%>%
group_by(GQ,Depth,GlMethod,QsBeta,Platform,GcMethod) %>%
summarise(
nDiscordant=mean(nDiscordant),
nConcordant=mean(nConcordant),
nSitesGQ=mean(nSitesGQ),
nHomToHomDiscordant=mean(nHomToHomDiscordant),
nHomToHetDiscordant=mean(nHomToHetDiscordant),
nHetToHomDiscordant=mean(nHetToHomDiscordant),
nHetToHetDiscordant=mean(nHetToHetDiscordant),
nHomToHomConcordant=mean(nHomToHomConcordant),
nHetToHetConcordant=mean(nHetToHetConcordant),
nSitesNon0DP=mean(nSitesNon0DP)) %>% 
ungroup() -> sumd

save(sumd,file=outfile)



#




