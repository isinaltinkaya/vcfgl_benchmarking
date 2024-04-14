# #!/usr/bin/env Rscript

# args<-commandArgs(trailingOnly=TRUE)
# nSNPs<-args[1]
# outfile<-args[2]
# infiles<-args[-c(1,2)]

# sd<-read.csv(nSNPs,sep="\t",header=FALSE)
# colnames(sd)<-c("nSNPs","Rep","Depth","GlMethod","QsBeta","Platform","GcMethod")


# hdr<-c("Sample","GQ","nDiscordant","nHomToHomDiscordant","nHomToHetDiscordant","nHetToHomDiscordant","nHetToHetDiscordant","nConcordant","nHomToHomConcordant","nHetToHetConcordant","nSitesComparedForSample","Rep","Depth","GlMethod","QsBeta","Platform","GcMethod")

# d<-NULL
# for(i in seq_along(infiles)){
# di<-read.csv(infiles[i],sep="\t",header=FALSE)
# colnames(di)<-hdr

# di$nSitesGQ<-di$nDiscordant+di$nConcordant
# di<-merge(di,sd,by=c("Rep","Depth","GlMethod","QsBeta","Platform","GcMethod"))

# d<-rbind(d,di)
# }

# save(d,file=outfile)



#!/usr/bin/env Rscript

args<-commandArgs(trailingOnly=TRUE)
# nSNPs<-args[1]
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

save(d,file=outfile)


