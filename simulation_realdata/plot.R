
## START 241029 ##


rm(list=ls())


library(ggplot2)
library(tidyr)
library(dplyr)

simv<-"sim_v4"
simdir<-"./"

#plotdir<-paste0(simdir,simv,"/results/plots/")
plotdir<-paste0(simdir,simv,"/results/plots241029/")
if (!dir.exists(plotdir)){
  cat("-> Creating plotdir at ",plotdir,"\n")
  dir.create(plotdir)
}

usecommon<-"usecommonv3"




  resdir<-paste0(simdir,simv,"/results/")
  prefix<-paste0(usecommon,"_gt_tgt_")
  indata<-read.csv(paste0(resdir,"HG00096_chr21_",prefix,"discordance.tsv"),sep=',',header = FALSE)

  colnames(indata)<-c("Category","Type","Value","Rep","TargetDepth","ActualDepth","GL")

  indata$Rep<-paste0(indata$Rep,indata$Method)
  indata$Rep<-as.factor(indata$Rep)
  indata$Depth<-indata$TargetDepth
  indata$Depth<-as.factor(indata$Depth)
  indata$GL<-as.factor(indata$GL)

  indata%>%filter(GL==1)->indata

  if(length(levels(indata$GL))>1){
    stop("More than one GL")
  }
  indata$GL<-NULL

  indata$Category<-as.factor(indata$Category)
  indata$Type<-as.factor(indata$Type)


# vcfglv4: no error
# vcfglv2:
indata %>% 
filter(Category=="DiscordanceRateOverCalls")%>%
filter(!Type %in%c("vcfglv2","vcfglv3"))%>%
# rename vcfglv's
mutate(
  Type= case_when(
    Type=="vcfglv1" ~ "Simulation with\nestimated error\n",
    Type=="vcfglv4" ~ "Simulation with\nno error\n",
    Type=="subsamplev1" ~ "Subsampled\nreal data\n",
    TRUE ~ "OTHER"
  )
)%>% 
ggplot(aes(x=Depth, y=Value*100, color=Type))+
  #geom_boxplot(alpha=0.5,outlier.shape = NA)+
  geom_boxplot(alpha=0.5)+
  #geom_point(position=position_jitterdodge(),aes(group=interaction(Type,Rep)))+
  stat_summary(position=position_dodge(width=0.75),fun.y=mean,geom="point",shape=23,size=3,fill="white")+
  stat_summary(position=position_dodge(width=0.75),fun.y=mean,geom="line",aes(group=Type))+
  labs(
    x="Depth",
    y="% Discordance rate",
  )+
  theme_bw()+
  theme(
    legend.position = "top",
    legend.text = element_text(size=8)
  )

ggsave(paste0(plotdir,prefix,
              "x-Depth_y-DiscordanceRateOverCalls_boxplot-replicates_line-boxplotmeans.png"
              ),
              width=6,height=6)


  # (1) REPLICATES
  # (1.1) visualize the distribution of discordance rates for different replicates at same depth and method


  ## "nHetToHetConcordant","nHetToHetDiscordant","nHomToHomConcordant","nHomToHomDiscordant","nHetToHomDiscordant","nHomToHetDiscordant"
  #indata%>%
  #filter(Category=="nHetToHetConcordant")%>%
  #ggplot(aes(x=TargetDepth, y=Value, color=Type,group=interaction(Type,Rep))) +
  #geom_point()+
  #geom_line()+
  #scale_y_continuous(labels = scales::comma)+
  #labs(
  #  x="Depth",
  #  y="Number of heterozygous concordant genotype calls"
  #)+
  #theme_bw()

  #ggsave(paste0(plotdir,prefix,
  #              "x-Depth_y-nHetToHetConcordant.png"
  #              ),
  #              width=10,height=10)


  #indata%>%
  #filter(Category=="nHetToHetDiscordant")%>%
  #ggplot(aes(x=Depth, y=Value, color=Type,group=interaction(Type,Rep))) +
  #geom_point()+
  #geom_line()+
  #scale_y_continuous(labels = scales::comma)+
  #labs(
  #  x="Depth",
  #  y="Number of heterozygous discordant genotype calls"
  #)+
  #theme_bw()


  #ggsave(paste0(plotdir,prefix,
  #              "x-Depth_y-nHetToHetDiscordant.png"
  #              ),
  #              width=10,height=10)

  #indata%>%
  #filter(Category=="nHomToHomConcordant")%>%
  #ggplot(aes(x=Depth, y=Value, color=Type,group=interaction(Type,Rep))) +
  #geom_point()+
  #geom_line()+
  #scale_y_continuous(labels = scales::comma)+
  #labs(
  #  x="Depth",
  #  y="Number of heterozygous concordant genotype calls"
  #)+
  #theme_bw()


  #ggsave(paste0(plotdir,prefix,
  #              "x-Depth_y-nHomToHomConcordant.png"
  #              ),
  #              width=10,height=10)


  indata%>%
  filter(Category=="nHomToHomDiscordant")%>%
  ggplot(aes(x=Depth, y=Value, color=Type,group=interaction(Type,Rep))) +
  geom_point()+
  geom_line()+
  scale_y_continuous(labels = scales::comma)+
  labs(
    x="Depth",
    y="Number of homozygous discordant genotype calls"
  )+
  theme_bw()

  ggsave(paste0(plotdir,prefix,
                "x-Depth_y-nHomToHomDiscordant.png"
                ),
                width=10,height=10)


  indata%>%
  filter(Category=="nHetToHomDiscordant")%>%
  ggplot(aes(x=Depth, y=Value, color=Type,group=interaction(Type,Rep))) +
  geom_point()+
  geom_line()+
  scale_y_continuous(labels = scales::comma)+
  labs(
    x="Depth",
    y="Number of HET→HOM discordant genotype calls"
  )+
  theme_bw()

  ggsave(paste0(plotdir,prefix,
                "x-Depth_y-nHetToHomDiscordant.png"
                ),
                width=10,height=10)


  indata%>%
  filter(Category=="nHomToHetDiscordant")%>%
  ggplot(aes(x=Depth, y=Value, color=Type,group=interaction(Type,Rep))) +
  geom_point()+
  geom_line()+
  scale_y_continuous(labels = scales::comma)+
  labs(
    x="Depth",
    y="Number of HOM→HET discordant genotype calls"
  )+
  theme_bw()

  ggsave(paste0(plotdir,prefix,
                "x-Depth_y-nHomToHetDiscordant.png"
                ),
                width=10,height=10)

  indata%>%
  filter(Category=="nSitesInCall")

  indata%>%
  filter(Category=="nSitesInCall")%>%
  ggplot(aes(x=Depth, y=Value, color=Type,group=interaction(Type,Rep))) +
    geom_point()+
    geom_line()+
    scale_y_continuous(labels = scales::comma)+
    labs(
      x="Depth",
      y="Number of sites with data"
    )+
    theme_bw()

  ggsave(paste0(plotdir,prefix,
                "x-Depth_y-nSitesInCall.png"
                ),
                width=10,height=10)

  indata%>%
  filter(Category=="nSitesInRefNotinCall")%>%
  ggplot(aes(x=Depth, y=Value, color=Type,group=interaction(Type,Rep))) +
    geom_point()+
    geom_line()+
    scale_y_continuous(labels = scales::comma)+
    labs(
      x="Depth",
      y="Number of sites with no data"
    )+
    theme_bw()

  ggsave(paste0(plotdir,prefix,
                "x-Depth_y-nSitesInRefNotinCall.png"
                ),
                width=10,height=10)


  indata%>%
  filter(Category=="nSitesCompared")%>%
  ggplot(aes(x=Depth, y=Value, color=Type,group=interaction(Type,Rep))) +
    geom_point()+
    geom_line()+
    scale_y_continuous(labels = scales::comma)+
    labs(
      x="Depth",
      y="Number of sites compared",
    )+
    theme_bw()

  ggsave(paste0(plotdir,prefix,
                "x-Depth_y-nSitesCompared.png"
                ),
                width=10,height=10)


  indata%>%
  filter(Category=="MissingnessRate")%>%
  ggplot(aes(x = Depth, y = Value, color=Type,group=interaction(Type,Rep))) +
    geom_point()+
    geom_line()+
    labs(
      x="Depth",
      y="Missingness rate",
    )+
    theme_bw()

  ggsave(paste0(plotdir,prefix,
                "x-Depth_y-MissingnessRate.png"
                ),
                width=10,height=10)


  indata%>%
  filter(Category=="DiscordanceRateOverCalls")%>%
  ggplot(aes(x = Depth, y = Value, color=Type,group=interaction(Type,Rep))) +
    geom_point()+
    geom_line()+
    labs(
      x="Depth",
      y="Discordance rate",
    )+
    theme_bw()

  ggsave(paste0(plotdir,prefix,
                "x-Depth_y-DiscordanceRateOverCalls.png"
                ),
                width=10,height=10)


  indata%>%
  filter(Category%in%c("nSitesHetInCall","nSitesHomInCall"))%>%
  group_by(Type,Depth,Category) %>% summarise(Value=mean(Value)) %>%
  ungroup()%>%
  spread(Category,Value)%>%
  mutate(PctHetInCall=nSitesHetInCall/(nSitesHetInCall+nSitesHomInCall))%>%
  ggplot(aes(x = Depth, y = PctHetInCall, color=Type,group=Type)) +
    geom_point()+
    geom_line()+
    # convert to percentage
    scale_y_continuous(labels = scales::percent)+
    labs(
      x="Depth",
      y="% Heterozygous genotype calls"
    )+
    theme_bw()

  ggsave(paste0(plotdir,prefix,
                "x-Depth_y-PctHetGenotypeCalls.png"
                ),
                width=10,height=10)


  indata%>%
  filter(Category%in%c("nSitesHetInCall","nSitesHomInCall"))%>%
  group_by(Type,Depth,Category) %>% summarise(Value=mean(Value)) %>%
  ungroup()%>%
  spread(Category,Value)%>%
  mutate(PctHomInCall=nSitesHomInCall/(nSitesHetInCall+nSitesHomInCall))%>%
  ggplot(aes(x = Depth, y = PctHomInCall, color=Type,group=Type)) +
    geom_point()+
    geom_line()+
    # convert to percentage
    scale_y_continuous(labels = scales::percent)+
    labs(
      x="Depth",
      y="% Homozygous genotype calls"
    )+
    theme_bw()

  ggsave(paste0(plotdir,prefix,
                "x-Depth_y-PctHomGenotypeCalls.png"
                ),
                width=10,height=10)


  # (2) REPLICATE MEANS

  # for each Type named vcfglv* and each Category, calculate the mean of the Value for each Depth


for(vcfglv in c("vcfglv1","vcfglv2","vcfglv3")){
  indata%>% 
  filter(Type %in% vcfglv)%>%
  filter(Category=="DiscordanceRateOverCalls")%>%
  group_by(Type,Depth) %>% summarise(MeanDiscordanceRateOverCalls=mean(Value),ActualDepth=unique(ActualDepth))%>%
  ungroup()%>%
  ggplot(aes(x = ActualDepth, y = MeanDiscordanceRateOverCalls,group = Type, color=Type)) +
    geom_point() +
    geom_line()+
    labs(
      x="Depth",
      y="Discordance Rate",
    )+
    theme_bw()

  ggsave(paste0(plotdir,prefix,
                "x-ActualDepth_y-MeanDiscordanceRateOverCalls_vcfglv",vcfglv,".png"
                ),
                width=10,height=10)
}

#TODO checkme
  indata%>% 
  filter(!Type %in% c("vcfglv1","vcfglv2"))%>%
  filter(Category%in%c("DiscordanceRateOverCalls"))%>%
  group_by(Type,Depth) %>% summarise(MeanDiscordanceRateOverCalls=mean(Value),ActualDepth=unique(ActualDepth))%>%
  ungroup()%>%
  ggplot(aes(x = Depth, y = MeanDiscordanceRateOverCalls,group = Type, color=Type)) +
    geom_point() +
    geom_line()+
    labs(
      x="Depth",
      y="Discordance Rate",
    )+
    theme_bw()

  ggsave(paste0(plotdir,prefix,
                "x-Depth_y-MeanDiscordanceRateOverCalls.png"
                ),
                width=10,height=10)



  indata %>%
  filter(Category=="DiscordanceRateOverCalls")%>%
  ggplot(aes(x = Depth, y = Value,group = interaction(Type,Rep), color=Type)) +
    geom_point() +
    geom_line()+
    labs(
      x="Depth",
      y="Discordance Rate",
    )+
    theme_bw()

  ggsave(paste0(plotdir,prefix,
                "x-Depth_y-DiscordanceRateOverCalls.png"
                ),
                width=10,height=10)
              

  indata %>%
  filter(Category=="DiscordanceRateOverCalls")%>%
  group_by(Type,Depth) %>% summarise(MeanDiscordanceRateOverCalls=mean(Value)) %>% 
  ungroup() %>%
  ggplot(aes(x = Depth, y = MeanDiscordanceRateOverCalls,group = Type, color=Type)) +
    geom_point() +
    geom_line()+
    labs(
      x="Depth",
      y="Discordance Rate",
    )+
    theme_bw()

  ggsave(paste0(plotdir,prefix,
                "x-Depth_y-MeanDiscordanceRateOverCalls.png"
                ),
                width=10,height=10)
