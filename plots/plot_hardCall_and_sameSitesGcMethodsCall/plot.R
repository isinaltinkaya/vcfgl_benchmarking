###############################################################################
# isinaltinkaya
#

rm(list=ls())
gc()

################################################################################
# LOAD LIBRARIES

# require(data.table) 
# require(readr)

# library(gridExtra)
library(ggplot2)
library(reshape2)
# library(ggpubr)
library(dplyr)
# library(tidyr) # pivot_longer

library(ggnewscale)

today<-format(Sys.Date(), "%y%m%d")


# setwd(paste0("plots_",today))
plot_outdir<-paste0(getwd(),"/")

color_betavars<-scale_color_manual(labels=c("0"="0","7"=expression(10^-7), "6"=expression(10^-6),"5"=expression(10^-5)), values=c("0"="black","7"="green","6"="blue","5"="red"))
fill_betavars<-scale_fill_manual(labels=c("0"="0","7"=expression(10^-7), "6"=expression(10^-6),"5"=expression(10^-5)), values=c("0"="black","7"="green","6"="blue","5"="red"))
color_brewer_betavars<-scale_color_brewer(palette="Set1",labels=c("0"="0","5"=expression(10^-5),"6"=expression(10^-6),"7"=expression(10^-7)))
fill_brewer_betavars<-scale_fill_brewer(palette="Set1",labels=c("0"="0","5"=expression(10^-5),"6"=expression(10^-6),"7"=expression(10^-7)))
linetype_gls<-scale_linetype_manual(values=c("gl1"="solid","gl2"="dashed","precise1_gl2"="dotted"),labels=c("gl1"="GL method 1","gl2"="GL method 2","precise1_gl2"="GL method 2 with precise error"))

################################################################################
# FUNCTIONS

assert <- function(x) {
  if (!x) stop("Assertion failed.")
}

save_plt <- function(plt,plot_filename,plot_outdir,overwrite=FALSE,...){
  if(!grepl("/$",plot_outdir)){
    plot_outdir<-paste0(plot_outdir,"/")
  }
  dir.create(plot_outdir, showWarnings = FALSE, recursive = TRUE)
  plot_fullpath<-paste0(plot_outdir,plot_filename)
  if(file.exists(plot_fullpath)){
    if(overwrite){
      print(paste0("Overwriting ",plot_fullpath))
      ggsave(plot_fullpath, plt,...)
    }else{
      stop(paste0("\n\n\n\n\nFile ",plot_fullpath," already exists.\n\n\n\n\n"))
    }
  }else{
    print(paste0("Saving plot to ",plot_fullpath))
    ggsave(plot_fullpath, plt,...)
  }
}


parse_betavar <- function(x){
  for(i in seq_along(levels(x))){
    if(levels(x)[i] != 0){
      levels(x)[i]<-gsub("X","10^-",paste0("X",levels(x)[i]))
    }
  }
  x
}

gcmethod_lut<-c("Genotype calling\nacross populations","Genotype calling\nwithin populations","Basic genotype calling")
names(gcmethod_lut)<-c("genotype_calling","genotype_calling_perpop","hard_genotype_call")

glmethod_lut<-c("GL method 1","GL method 2","GL method 2 \nwith precise error")
names(glmethod_lut)<-c("gl1","gl2","precise1_gl2")

betavar_lut<-c("0","10^-7","10^-6","10^-5")
names(betavar_lut)<-c("0","7","6","5")


labeller1fn <- function(variable,value){
  if (variable=='Betavar') {
      return(label_parsed(ifelse(value==0,"0",paste0("10^-",value))))
  }else if (variable=='Depth'){
    return(as.character(paste0("Depth: ",as.character(value))))
  }else if (variable=='GcMethod'){
    return(as.character(gcmethod_lut[as.character(value)]))
  }else if (variable=='Gl'){
    return(as.character(glmethod_lut[as.character(value)]))
  }else{
    return(as.character(value))
  }
}


glmethod_lut2<-c("GL1","GL2","P_GL2")
names(glmethod_lut2)<-c("gl1","gl2","precise1_gl2")

labeller2fn <- function(variable,value){
  if (variable=='Betavar') {
      return(label_parsed(ifelse(value==0,"0",paste0("10^-",value))))
  }else if (variable=='Depth'){
    return(as.character(paste0("Depth: ",as.character(value))))
  }else if (variable=='GcMethod'){
    return(as.character(gcmethod_lut[as.character(value)]))
  }else if (variable=='Gl'){
    return(as.character(glmethod_lut2[as.character(value)]))
  }else{
    return(as.character(value))
  }
}

################################################################################
# READ DATA




simdir="../../simulation/sim/"

# GlMethod_i=c("gl1","gl2","precise1_gl2")
# hdr<-c("Sample","nSitesTotal","nSitesRetained","nSitesCompared","nSitesCallMis","nDiscordantSites","nSitesInTrueNotCall","nConcordantSites","MissingnessRate","DiscordanceRate","ConcordanceRate","T_Hom","T_Het","F_HomToHom","F_HomToHet","F_HetToHom","F_HetToHet","T_Hom_rate","T_Het_rate","F_HomToHom_rate","F_HomToHet_rate","F_HetToHom_rate","F_HetToHet_rate","Rep","Depth","Betavar")

# d<-NULL
# for(gli in seq_along(GlMethod_i)){

#   di<-read.csv(paste0(simdir, 
#   "sim_vcfgl_2312/model_OutOfAfrica_3G09/gc_evaluation/genotype_discordance_",GlMethod_i[gli],"/sim_vcfgl_2312-OutOfAfrica_3G09-hard_genotype_call.tsv"),sep="\t",header=FALSE)
#   colnames(di)<-hdr
#   di$Gl<-GlMethod_i[gli]
#   d<-rbind(d,di)
# }



# sd<-read.csv(paste0(simdir,"sim_vcfgl_2312/model_OutOfAfrica_3G09/stats/nSites_non0dp/sim_vcfgl_2312-OutOfAfrica_3G09_nSitesNon0DP.csv"),sep=",",header=TRUE)
# d<-merge(d,sd,by=c("Sample","Rep","Depth"))

# save(d,file=paste0("./data_",today,".RData"))

load(file=paste0("./data_",today,".RData"))


################################################################################




rm(d1,d2,d3,plt);gc();


d1<-d 
d1$Sample<-NULL
d1$Rep<-NULL


(d1%>% 
group_by(Gl,Betavar,Depth)%>%
summarize(
  nDiscordant=sum(nDiscordantSites),
  nConcordant=sum(nConcordantSites),
  nSitesRetained=sum(nSitesRetained),
  nSitesCompared=sum(nSitesCompared),
  nSitesCallMis=sum(nSitesCallMis),
  nSitesNon0DP=sum(nSitesNon0DP),
  nSitesInTrueNotCall=sum(nSitesInTrueNotCall),
  meanMissingnessRate=mean(MissingnessRate),
  meanDiscordanceRate=mean(DiscordanceRate),
  meanConcordanceRate=mean(ConcordanceRate)
)%>% 
ungroup()%>%
mutate(DiscordanceRate=nDiscordant/(nDiscordant+nConcordant),
        ConcordanceRate=nConcordant/(nDiscordant+nConcordant),
        MissingnessRate=nSitesCallMis/nSitesCompared,
        DiscordanceRateOverRetained=nDiscordant/nSitesRetained,
        ConcordanceRateOverRetained=nConcordant/nSitesRetained,
        MissingnessRateOverRetained=nSitesCallMis/nSitesRetained,
        DiscordanceRateOverCompared=nDiscordant/nSitesCompared,
        ConcordanceRateOverCompared=nConcordant/nSitesCompared,
        MissingnessRateOverCompared=nSitesCallMis/nSitesCompared,
        DiscordanceRateOverNon0DP=nDiscordant/nSitesNon0DP,
        MissingnessRateOverNon0DP=(nSitesNon0DP-(nDiscordant+nConcordant))/nSitesNon0DP,
        )%>%
mutate(across(c(Betavar,Depth,Gl),as.factor))%>% 
mutate(Betavar=factor(Betavar,levels=c("0","7","6","5"))))->d2



d2%>%
ggplot(aes(x="",y=DiscordanceRate,color=Betavar,group=Betavar))+
  facet_wrap(Gl~Depth,scale="free",labeller=labeller1fn,ncol=6)+
  geom_point()+
  theme_bw()+
  theme(legend.position="top")+
  theme(strip.background = element_rect(colour="black", fill="white"))+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  color_brewer_betavars+
  labs(
    x="",
    y="Discordance rate",
    color="Beta variance",
    fill=""
  )+
  color_brewer_betavars
  
save_plt(plt=last_plot(),plot_filename="plot-point_x-None_y-DiscordanceRate_color-Betavar_group-Betavar_facetwrap-GlDepth.png",plot_outdir=plot_outdir,overwrite=TRUE,width=12,height=10)



bv5<-d[d$Betavar=="5",]
bv0<-d[d$Betavar=="0",]

bv0$diff0Minus5 <- bv0$DiscordanceRate - bv5$DiscordanceRate
sum(bv0$Sample!=bv5$Sample)
sum(bv0$Rep!=bv5$Rep)
sum(bv0$Depth!=bv5$Depth)



bv0 %>%
mutate(across(c(Betavar,Depth,Gl),as.factor))%>%
ggplot(aes(y=diff0Minus5))+
  facet_grid(Gl~Depth, scale="free",labeller=label_both)+
  geom_boxplot()+
  theme_bw()+
  geom_hline(yintercept=0,linetype="dashed")+
  theme(legend.position="top")+
  labs(y="Difference in discordance rate (beta=0 - beta=10^-5)")


ggsave(paste0(plot_outdir,"plot_diff0Minus5.png"),width=10,height=10)


d%>%
group_by(Rep,Depth,Betavar,Gl)%>%
summarize(
  nDiscordant=sum(nDiscordantSites),
  nConcordant=sum(nConcordantSites),
  nSitesRetained=sum(nSitesRetained),
  nSitesCompared=sum(nSitesCompared),
  nSitesCallMis=sum(nSitesCallMis),
  nSitesInTrueNotCall=sum(nSitesInTrueNotCall),
  meanMissingnessRate=mean(MissingnessRate),
  meanDiscordanceRate=mean(DiscordanceRate),
  meanConcordanceRate=mean(ConcordanceRate),
  DiscordanceRate=sum(nDiscordantSites)/(sum(nDiscordantSites)+sum(nConcordantSites)),
)%>% 
ungroup()%>%
mutate(across(c(Betavar,Depth,Gl),as.factor))%>%
mutate(Betavar=factor(Betavar,levels=c("0","7","6","5")))%>%
ggplot(aes(x="",y=DiscordanceRate,color=Betavar,group=Betavar))+
  facet_grid(Depth~Gl,scale="free",labeller=labeller1fn)+
  geom_boxplot()+
  geom_point(position=position_dodge(width=0.75))+
  theme_bw()+
  theme(legend.position="top")+
  labs(
    x="",
    y="Discordance rate",
    color="Beta variance",
    fill="Beta variance",
    subtitle="Each point represents a simulation replicate"
  )+
  theme(strip.background = element_rect(colour="black", fill="white"))+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  color_brewer_betavars


save_plt(plt=last_plot(),plot_filename="plot-boxplot_x-Depth_y-DiscordanceRate_color-Betavar_group-Betavar_facet-GlDepth_point-Rep.png",plot_outdir=plot_outdir,overwrite=TRUE,width=10,height=10)




# each point: one sample from one rep 
d%>%
group_by(Sample,Rep,Depth,Betavar,Gl)%>%
summarize(
  nDiscordant=sum(nDiscordantSites),
  nConcordant=sum(nConcordantSites),
  nSitesRetained=sum(nSitesRetained),
  nSitesCompared=sum(nSitesCompared),
  nSitesCallMis=sum(nSitesCallMis),
  nSitesInTrueNotCall=sum(nSitesInTrueNotCall),
  meanMissingnessRate=mean(MissingnessRate),
  meanDiscordanceRate=mean(DiscordanceRate),
  meanConcordanceRate=mean(ConcordanceRate)
)%>% 
ungroup()%>%
mutate(across(c(Betavar,Depth,Gl),as.factor))%>%
mutate(Betavar=factor(Betavar,levels=c("0","7","6","5")))%>%
ggplot(aes(x=Depth,y=meanDiscordanceRate,color=Betavar,group=Betavar))+
  facet_wrap(Gl~Depth,scale="free",labeller=labeller1fn)+
  geom_boxplot()+
  geom_point(position=position_dodge(width=0.75))+
  theme_bw()+
  theme(legend.position="top")+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  color_brewer_betavars


ggsave(paste0(plot_outdir,"plot-boxplot_x-Depth_y-meanDiscordanceRate_color-Betavar_group-Betavar_facet-GlDepth_point-SampleRep.png"),width=10,height=10)






# ###############################################################################

# READ DATA FROM OTHER GC METHODS
# AND MERGE ALL

GlMethod_i=c("gl1","gl2","precise1_gl2")
hdr<-c("Sample","nSitesTotal","nSitesRetained","nSitesCompared","nSitesCallMis","nDiscordantSites","nSitesInTrueNotCall","nConcordantSites","MissingnessRate","DiscordanceRate","ConcordanceRate","T_Hom","T_Het","F_HomToHom","F_HomToHet","F_HetToHom","F_HetToHet","T_Hom_rate","T_Het_rate","F_HomToHom_rate","F_HomToHet_rate","F_HetToHom_rate","F_HetToHet_rate","Rep","Depth","Betavar")

Method_i=c("genotype_calling", "genotype_calling_perpop", "hard_genotype_call")

d<-NULL
for(mi in seq_along(Method_i)){
for(gli in seq_along(GlMethod_i)){
  di<-read.csv(paste0(simdir, 
  "sim_vcfgl_2312/model_OutOfAfrica_3G09/gc_evaluation/genotype_discordance_",GlMethod_i[gli],"/sim_vcfgl_2312-OutOfAfrica_3G09-",Method_i[mi],".tsv"),sep="\t",header=FALSE)
  colnames(di)<-hdr
  di$Gl<-GlMethod_i[gli]
  di$GcMethod<-Method_i[mi]
  d<-rbind(d,di)
}
}




sd<-read.csv(paste0(simdir,"sim_vcfgl_2312/model_OutOfAfrica_3G09/stats/nSites_non0dp/sim_vcfgl_2312-OutOfAfrica_3G09_nSitesNon0DP.csv"),sep=",",header=TRUE)
d<-merge(d,sd,by=c("Sample","Rep","Depth"))

save(d,file=paste0("./data_",today,".RData"))

load(file=paste0("./data_",today,".RData"))

################################################################################

rm(d1,d2,d3,plt);gc();


d1<-d 
d1$Sample<-NULL
d1$Rep<-NULL


(d1%>% 
group_by(Gl,Betavar,Depth,GcMethod)%>%
summarize(
  nDiscordant=sum(nDiscordantSites),
  nConcordant=sum(nConcordantSites),
  nSitesRetained=sum(nSitesRetained),
  nSitesCompared=sum(nSitesCompared),
  nSitesCallMis=sum(nSitesCallMis),
  nSitesNon0DP=sum(nSitesNon0DP),
  nSitesInTrueNotCall=sum(nSitesInTrueNotCall),
  meanMissingnessRate=mean(MissingnessRate),
  meanDiscordanceRate=mean(DiscordanceRate),
  meanConcordanceRate=mean(ConcordanceRate)
)%>% 
ungroup()%>%
mutate(DiscordanceRate=nDiscordant/(nDiscordant+nConcordant),
        ConcordanceRate=nConcordant/(nDiscordant+nConcordant),
        MissingnessRate=nSitesCallMis/nSitesCompared,
        DiscordanceRateOverRetained=nDiscordant/nSitesRetained,
        ConcordanceRateOverRetained=nConcordant/nSitesRetained,
        MissingnessRateOverRetained=nSitesCallMis/nSitesRetained,
        DiscordanceRateOverCompared=nDiscordant/nSitesCompared,
        ConcordanceRateOverCompared=nConcordant/nSitesCompared,
        MissingnessRateOverCompared=nSitesCallMis/nSitesCompared,
        DiscordanceRateOverNon0DP=nDiscordant/nSitesNon0DP,
        MissingnessRateOverNon0DP=(nSitesNon0DP-(nDiscordant+nConcordant))/nSitesNon0DP,
        )%>%
mutate(across(c(Betavar,Depth,Gl,GcMethod),as.factor))%>% 
mutate(Betavar=factor(Betavar,levels=c("0","7","6","5"))))->d2




d2%>%
ggplot(aes(x=Depth,y=DiscordanceRate,color=Betavar,group=Betavar))+
  facet_wrap(Gl~GcMethod,scale="free",labeller=labeller1fn)+
  geom_point()+
  geom_line()+
  theme_bw()+
  theme(legend.position="top")+
  theme(strip.background = element_rect(colour="black", fill="white"))+
  #rotate x axis labels 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  # theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  color_brewer_betavars+
  labs(
    x="",
    y="Discordance rate",
    color="Beta variance",
    fill=""
  )+
  color_brewer_betavars

save_plt(plt=last_plot(),plot_filename="plot-pointLine_x-Depth_y-DiscordanceRate_color-Betavar_group-Betavar_facetwrap-GlGcMethod.png",plot_outdir=plot_outdir,overwrite=TRUE,width=12,height=10)


d2%>%
filter(Gl%in%c("gl1","gl2"))->d3
d3%>%
ggplot(aes(x=Gl,y=DiscordanceRate,color=Betavar,group=Betavar))+
  facet_wrap(GcMethod~Depth,scale="free",labeller=labeller1fn,ncol=6)+
  geom_point()+
  theme_bw()+
  theme(legend.position="top")+
  theme(strip.background = element_rect(colour="black", fill="white"))+
  #rotate x axis labels 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  # theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  color_brewer_betavars+
  labs(
    x="",
    y="Discordance rate",
    color="Beta variance",
    fill=""
  )+
  color_brewer_betavars+
  scale_x_discrete(labels=c("gl1"="GL 1","gl2"="GL 2","precise1_gl2"="NA"))+
  #make y axis percentage automatically
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.01))

save_plt(plt=last_plot(),plot_filename="plot-point_x-Gl_y-DiscordanceRate_color-Betavar_group-Betavar_facetwrap-GcMethodDepth.png",plot_outdir=plot_outdir,overwrite=TRUE,width=12,height=10)










mypalette <- colorRampPalette(RColorBrewer::brewer.pal(9, "Reds"), space = "Lab")
mypal<-mypalette(9)


#heatmap
d2%>%
filter(Gl%in%c("gl1","gl2"))->d3

myAes1<-aes(y=Betavar,x=Gl,fill=DiscordanceRate)

getTilePanel1<-function(DT,DPVAL,GCMVAL,RMY=FALSE,RMGRIDROWS=TRUE,RMGRIDCOLS=TRUE){
(DT%>%
filter(Depth==DPVAL)%>%
filter(GcMethod==GCMVAL)%>%
ggplot(aes(x=Gl,y=Betavar,fill=DiscordanceRate))+
  geom_tile()+
  # facet_grid(
  #   # rows=vars(GcMethod),
  #   cols=vars(Depth),
  #   # scales="free",
  #   space="free",
  #   labeller=labeller2fn
  # )+
  coord_equal()+
  theme(legend.box = "horizontal",legend.key.size = unit(0.2, "cm"),legend.key.width = unit(0.5, "cm"),legend.key.height = unit(0.5, "cm"),legend.text = element_text(size = 8),legend.title = element_text(size = 8))+
  theme(panel.background = element_rect(fill = "white", colour = "white"))+
  theme_bw()+
  theme(legend.position="top")+
  theme(strip.background = element_rect(colour="black", fill="white"))+
  labs(
    x="",
    y="",
    fill="",
    color=""
  )+
  #make fill labels percentage automatically
  scale_fill_gradientn(
    colours = mypal,
    #round percent by 3 digits
    labels = scales::percent_format(accuracy = 0.01),
    breaks = scales::pretty_breaks(n = 4)
    # breaks = scales::pretty_breaks(n = 4)(DT[DT$Depth==DPVAL,]$DiscordanceRate),
  )+
  theme(legend.position = "bottom")+
  theme(legend.text=element_text(angle=90,hjust=0.5,vjust=0.5,size=8))+
  #make legend smaller
  theme(legend.key.size = unit(0.2, "cm"),legend.key.width = unit(0.5, "cm"),legend.key.height = unit(0.5, "cm"))+
  scale_y_discrete(labels=c("0"="0","7"=expression(10^-7), "6"=expression(10^-6),"5"=expression(10^-5)))+
  {if(RMY)theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())else theme(axis.text.y = element_text(size = 8), axis.ticks.y = element_line(size = 0.2))}+
  {if(RMGRIDROWS)theme(strip.text.y = element_blank())else theme(strip.text.y = element_text(size = 8,angle=0))}+
  {if(RMGRIDCOLS)theme(strip.text.x = element_blank())else theme(strip.text.x = element_text(size = 8,angle=0))}+
  scale_x_discrete(labels=c("gl1"="GL 1","gl2"="GL 2","precise1_gl2"="NA")))
}




gcmvalarr<-c("genotype_calling","genotype_calling_perpop","hard_genotype_call")
allTilePlotsList<-NULL
for(gcmvali in seq_along(gcmvalarr)){
  gcmval<-gcmvalarr[gcmvali]
  if(gcmvali==1)rm_grid_cols=FALSE else rm_grid_cols=TRUE
  p01<-getTilePanel1(d3,0.1,GCMVAL=gcmval,RMY=FALSE,RMGRIDCOLS=rm_grid_cols)
  p05<-getTilePanel1(d3,0.5,GCMVAL=gcmval,RMGRIDCOLS=rm_grid_cols)
  p1<-getTilePanel1(d3,1,GCMVAL=gcmval,RMGRIDCOLS=rm_grid_cols)
  p2<-getTilePanel1(d3,2,GCMVAL=gcmval,RMGRIDCOLS=rm_grid_cols)
  p10<-getTilePanel1(d3,10,GCMVAL=gcmval,RMGRIDCOLS=rm_grid_cols)
  p20<-getTilePanel1(d3,20,GCMVAL=gcmval,RMGRIDROWS=FALSE,RMGRIDCOLS=rm_grid_cols)
  allTilePlotsList<-c(allTilePlotsList,list(p01,NULL,p05,NULL,p1,NULL,p2,NULL,p10,NULL,p20,NULL))
}

ggarrange(plotlist=allTilePlotsList,ncol=12,nrow=3,common.legend=FALSE,legend="bottom",widths=c(100,0),align="hv")

gcmvalarr<-c("hard_genotype_call")
allTilePlotsList<-NULL
for(gcmvali in seq_along(gcmvalarr)){
  gcmval<-gcmvalarr[gcmvali]
  if(gcmvali==1)rm_grid_cols=FALSE else rm_grid_cols=TRUE
  p01<-getTilePanel1(d3,0.1,GCMVAL="hard_genotype_call",RMY=FALSE,RMGRIDCOLS=rm_grid_cols)
  p05<-getTilePanel1(d3,0.5,GCMVAL="hard_genotype_call",RMGRIDCOLS=rm_grid_cols)
  p1<-getTilePanel1(d3,1,GCMVAL="hard_genotype_call",RMGRIDCOLS=rm_grid_cols)
  p2<-getTilePanel1(d3,2,GCMVAL="hard_genotype_call",RMGRIDCOLS=rm_grid_cols)
  p10<-getTilePanel1(d3,10,GCMVAL="hard_genotype_call",RMGRIDCOLS=rm_grid_cols)
  p20<-getTilePanel1(d3,20,GCMVAL="hard_genotype_call",RMGRIDROWS=FALSE,RMGRIDCOLS=rm_grid_cols)
  allTilePlotsList<-c(allTilePlotsList,list(p01,NULL,p05,NULL,p1,NULL,p2,NULL,p10,NULL,p20,NULL))
}

ggarrange(plotlist=allTilePlotsList,ncol=12,nrow=1,common.legend=FALSE,widths=c(100,0),align="hv",legend="bottom", labels = c("0.1","","0.5","","1","","2","","10","","20",""))+bgcolor("white")

save_plt(plt=last_plot(),plot_filename="plot-heatmap_x-Gl_y-Betavar_fill-DiscordanceRate_group-GcMethod_facetwrap-GlDepth.png",plot_outdir=plot_outdir,overwrite=TRUE,width=12,height=3)


p01<-getTilePanel1(d3,0.1,GCMVAL="genotype_calling_perpop",RMY=FALSE)
p05<-getTilePanel1(d3,0.5)
p1<-getTilePanel1(d3,1)
p2<-getTilePanel1(d3,2)
p10<-getTilePanel1(d3,10)
p20<-getTilePanel1(d3,20,RMGRIDROWS=FALSE)


ggarrange(p01,NULL,p05,NULL,p1,NULL,p2,NULL,p10,NULL,p20,ncol=12,common.legend=FALSE,legend="bottom",widths=c(100,0),align="v")




bv5<-d[d$Betavar=="5",]
bv0<-d[d$Betavar=="0",]

bv0$diff0Minus5 <- bv0$DiscordanceRate - bv5$DiscordanceRate
sum(bv0$Sample!=bv5$Sample)
sum(bv0$Rep!=bv5$Rep)
sum(bv0$Depth!=bv5$Depth)



bv0 %>%
mutate(across(c(Betavar,Depth,Gl),as.factor))%>%
ggplot(aes(y=diff0Minus5))+
  facet_grid(Gl~Depth, scale="free",labeller=label_both)+
  geom_boxplot()+
  theme_bw()+
  geom_hline(yintercept=0,linetype="dashed")+
  theme(legend.position="top")+
  labs(y="Difference in discordance rate (beta=0 - beta=10^-5)")


ggsave(paste0(plot_outdir,"plot_diff0Minus5.png"),width=10,height=10)


d%>%
group_by(Rep,Depth,Betavar,Gl)%>%
summarize(
  nDiscordant=sum(nDiscordantSites),
  nConcordant=sum(nConcordantSites),
  nSitesRetained=sum(nSitesRetained),
  nSitesCompared=sum(nSitesCompared),
  nSitesCallMis=sum(nSitesCallMis),
  nSitesInTrueNotCall=sum(nSitesInTrueNotCall),
  meanMissingnessRate=mean(MissingnessRate),
  meanDiscordanceRate=mean(DiscordanceRate),
  meanConcordanceRate=mean(ConcordanceRate),
  DiscordanceRate=sum(nDiscordantSites)/(sum(nDiscordantSites)+sum(nConcordantSites)),
)%>% 
ungroup()%>%
mutate(across(c(Betavar,Depth,Gl),as.factor))%>%
mutate(Betavar=factor(Betavar,levels=c("0","7","6","5")))%>%
ggplot(aes(x="",y=DiscordanceRate,color=Betavar,group=Betavar))+
  facet_grid(Depth~Gl,scale="free",labeller=labeller1fn)+
  geom_boxplot()+
  geom_point(position=position_dodge(width=0.75))+
  theme_bw()+
  theme(legend.position="top")+
  labs(
    x="",
    y="Discordance rate",
    color="Beta variance",
    fill="Beta variance",
    subtitle="Each point represents a simulation replicate"
  )+
  theme(strip.background = element_rect(colour="black", fill="white"))+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  color_brewer_betavars


save_plt(plt=last_plot(),plot_filename="plot-boxplot_x-Depth_y-DiscordanceRate_color-Betavar_group-Betavar_facet-GlDepth_point-Rep.png",plot_outdir=plot_outdir,overwrite=TRUE,width=10,height=10)




# each point: one sample from one rep 
d%>%
group_by(Sample,Rep,Depth,Betavar,Gl)%>%
summarize(
  nDiscordant=sum(nDiscordantSites),
  nConcordant=sum(nConcordantSites),
  nSitesRetained=sum(nSitesRetained),
  nSitesCompared=sum(nSitesCompared),
  nSitesCallMis=sum(nSitesCallMis),
  nSitesInTrueNotCall=sum(nSitesInTrueNotCall),
  meanMissingnessRate=mean(MissingnessRate),
  meanDiscordanceRate=mean(DiscordanceRate),
  meanConcordanceRate=mean(ConcordanceRate)
)%>% 
ungroup()%>%
mutate(across(c(Betavar,Depth,Gl),as.factor))%>%
mutate(Betavar=factor(Betavar,levels=c("0","7","6","5")))%>%
ggplot(aes(x=Depth,y=meanDiscordanceRate,color=Betavar,group=Betavar))+
  facet_wrap(Gl~Depth,scale="free",labeller=labeller1fn)+
  geom_boxplot()+
  geom_point(position=position_dodge(width=0.75))+
  theme_bw()+
  theme(legend.position="top")+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  color_brewer_betavars


ggsave(paste0(plot_outdir,"plot-boxplot_x-Depth_y-meanDiscordanceRate_color-Betavar_group-Betavar_facet-GlDepth_point-SampleRep.png"),width=10,height=10)











################################################################################

