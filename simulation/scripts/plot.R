# Plot

library(ggplot2)
library(dplyr)



d<-read.csv("~/Mount/vcfgl/vcfgl_paper_analyses/simulation/sim/sim_vcfgl_2310/model_OutOfAfrica_3G09/gc_evaluation/bcftools_gtcheck_e0/sim_vcfgl_2310-OutOfAfrica_3G09.tsv",sep="\t",header=F)


colnames(d)<-c("nDiscordantSites","HWE","nSites","Contig","Rep","Depth","ErrorRate","QsError","nSitesCompared","nSitesNoMatch")

d$Rep <- as.factor(d$Rep)
d$Depth <- as.factor(d$Depth)
d$DiscordanceRate <- d$nDiscordantSites/d$nSites
d$MissingnessRate <- 1-(d$nSites/d$nSitesCompared)
d$BetaVariance <- 10^d$QsError
d$BetaMean <- d$ErrorRate
d$QsError <- as.factor(d$QsError)

d2<-d%>%group_by(Depth,BetaVariance,BetaMean,QsError)%>%summarise_at(vars(DiscordanceRate),list(mean_DiscordanceRate=mean))


# p <- ggplot(d, aes(x = Depth, y = DiscordanceRate, colour = QsError)) +
  # geom_boxplot(notch = FALSE) +
  # labs(x = 'Depth', y = 'Discordance Rate')+
  # labs(colour = 'Beta distribution variance')+
  # scale_colour_brewer(palette = 'Set1') +
  # theme_bw() +
  # theme(
    # legend.position = 'bottom'
  # )
# p
#
#
# p<-ggplot(d, aes(x = MissingnessRate, y = DiscordanceRate, colour = QsError)) +
  # geom_point() +
  # facet_wrap( Depth ~. ,scale="free" ) +
  # theme_bw()
# # >
#

###

#columns:
# HET_SENSITIVITY HET_PPV HET_SPECIFICITY HOMVAR_SENSITIVITY      HOMVAR_PPV      HOMVAR_SPECIFICITY      VAR_SENSITIVITY VAR_PPV VAR_SPECIFICITY GENOTYPE_CONCORDANCE    NON_REF_GENOTYPE_CONCORDANCE

# d<-read.csv("~/Mount/vcfgl/vcfgl_paper_analyses/simulation/delme",sep="\t",header=F)

# colnames(d)<-c("QSERR","HET_SENSITIVITY","HET_PPV","HET_SPECIFICITY","HOMVAR_SENSITIVITY","HOMVAR_PPV","HOMVAR_SPECIFICITY","VAR_SENSITIVITY","VAR_PPV","VAR_SPECIFICITY","GENOTYPE_CONCORDANCE","NON_REF_GENOTYPE_CONCORDANCE")
# d$QSERR<-as.factor(d$QSERR)

